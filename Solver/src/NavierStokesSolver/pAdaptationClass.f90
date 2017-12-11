!
!////////////////////////////////////////////////////////////////////////
!
!      pAdaptationClass.f90
!      Created: December 10, 2017 at 12:56 PM 
!      By: Andrés Rueda
!
!      Class cotaining routines for adapting polynomial orders.
!        -> Currently, the adaptation procedure is only performed with  
!           truncation error estimations (TruncationErrorClass.f90), but other
!           strategies can be added.
!        -> The current implementation is only valid for shared memory 
!           parallelization (OpenMP). TODO: Update to MPI.
!
!////////////////////////////////////////////////////////////////////////
!
module pAdaptationClass
   use SMConstants
   use DGSEMClass
   use PhysicsStorage
   use InterpolationMatrices
   use TruncationErrorClass
   use FTValueDictionaryClass
   
   implicit none
   
   private
   public GetMeshPolynomialOrders, ReadOrderFile, WriteOrderFile
   public pAdaptation_t !, PrintTEmap
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type :: pAdaptation_t
      type(FTValueDictionary)           :: controlVariables ! Own copy of the control variables (needed to construct sem) ! TODO: Change this?
      character(len=LINE_LENGTH)        :: plotFileName     ! Name of file for plotting adaptation information
      real(kind=RP)                     :: reqTE            ! Requested truncation error
      logical                           :: PlotInfo
      logical                           :: Adapt            ! Is the adaptator going to be used??
      logical                           :: increasing       ! Performing an increasing adaptation procedure?
      logical                           :: Constructed      ! 
      integer                           :: NxyzMax(3)       ! Maximum polynomial order in all the directions
      
      type(TruncationError_t), allocatable :: TE(:)         ! Truncation error for every element(:)
      
      contains
         procedure :: construct => ConstructPAdaptator
         procedure :: destruct  => DestructPAdaptator
!~         procedure :: plot      => AdaptationPlotting
         procedure :: pAdaptTE
   end type pAdaptation_t
   
   interface
      character(len=LINE_LENGTH) function getFileName( inputLine )
         use SMConstants
         implicit none
         character(len=*)     :: inputLine
      end function getFileName
   end interface
!
!  ----------------
!  Module variables
!  ----------------
!
   !! Parameters
   integer, parameter    :: NMIN = 1                 ! Minimum polynomial order !! TODO: remove from here and use some criteria!!
   
   
   integer               :: NInc_0 = 4
!!    integer               :: dN_Inc = 3 
   integer               :: fN_Inc = 2
   integer               :: NInc
   integer               :: nelem   ! number of elements in mesh
   
   EXTERNAL                            :: externalStateForBoundaryName
   EXTERNAL                            :: ExternalGradientForBoundaryName
   
!========
 contains
!========
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   subroutine GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
      use ReadMeshFile
      implicit none
      !-------------------------------------------------
      type(FTValueDictionary), intent(in)    :: controlVariables
      integer, allocatable   , intent(inout) :: Nx(:), Ny(:), Nz(:)  
      integer                , intent(out)   :: Nmax
      !-------------------------------------------------
      integer                                :: nelem
      !-------------------------------------------------
      
      if (controlVariables % containsKey("polynomial order file")) then
         call ReadOrderFile( controlVariables % stringValueForKey("polynomial order file", requestedLength = LINE_LENGTH), &
                             Nx, Ny, Nz )
      else
         nelem = NumOfElemsFromMeshFile( controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH) )
         allocate( Nx(nelem), Ny(nelem), Nz(nelem) )
         
         if (controlVariables % containsKey("polynomial order")) then
            Nx = controlVariables % integerValueForKey("polynomial order")
            Ny = Nx
            Nz = Nx
         else
            if (controlVariables % containsKey("polynomial order i") .AND. &
                controlVariables % containsKey("polynomial order j") .AND. &
                controlVariables % containsKey("polynomial order k") ) then
               Nx = controlVariables % integerValueForKey("polynomial order i")
               Ny = controlVariables % integerValueForKey("polynomial order j")
               Nz = controlVariables % integerValueForKey("polynomial order k")
            else
               error stop "The polynomial order(s) must be specified"
            end if
         end if
      end if
      
      Nmax = 0
      if (controlVariables % containsKey("adaptation nmax i")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax i"))
      if (controlVariables % containsKey("adaptation nmax j")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax j"))
      if (controlVariables % containsKey("adaptation nmax k")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax k"))
      
      Nmax = max(Nmax,maxval(Nx),maxval(Ny),maxval(Nz))
      
      end subroutine GetMeshPolynomialOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReadOrderFile(filename, Nx, Ny, Nz)
      implicit none
!
!     ----------------------------------------------------------------------
!     Subroutine that reads input file containing polynomial orders for mesh
!     ----------------------------------------------------------------------
!
      character(len=*), intent(in) :: filename          !<  Name of file containing polynomial orders to initialize
      integer, allocatable         :: Nx(:),Ny(:),Nz(:) !>  Polynomial orders for each element
      !------------------------------------------
      integer                      :: fd       ! File unit
      integer                      :: nelem    ! Number of elements
      integer                      :: i        ! counter
      !------------------------------------------
      
      open(newunit = fd, FILE = filename )   
         READ(fd,*) nelem
         
         allocate(Nx(nelem),Ny(nelem),Nz(nelem))
         
         do i = 1, nelem
            READ(fd,*) Nx(i), Ny(i), Nz(i)
         ENDDO
      close(UNIT=fd)
      
   end subroutine ReadOrderFile
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine WriteOrderFile(sem,filename)
      implicit none
      !------------------------------------------
      type(DGSem)     , intent(in) :: sem      !<  Sem class
      character(len=*), intent(in) :: filename          !<  Name of file containing polynomial orders to initialize
      !------------------------------------------
      integer                      :: fd       ! File unit
      integer                      :: k
      !------------------------------------------
      
      open( newunit = fd , FILE = TRIM(filename), ACTION = 'write')
         
         write(fd,*) size(sem % mesh % elements)
         
         do k=1, size(sem % mesh % elements)
            write(fd,*) sem % mesh % elements(k) % Nxyz
         end do
         
      close (fd)
   end subroutine WriteOrderFile

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ROUTINES FOR ADAPTATION
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for constructing the p-adaptator.
!   -> If increasing (multi-stage) adaptation is selected, the final step is to rewrite the polynomial orders for the sem contruction
!  ----------------------------------------
   subroutine ConstructPAdaptator(this,Nx,Ny,Nz,controlVariables)
      use FTValueDictionaryClass
      implicit none
      !--------------------------------------
      class(pAdaptation_t)                :: this             !>  P-Adaptator
      integer, DIMENSION(:)               :: Nx,Ny,Nz         !<> Input: polynomial orders as read from input files - Output: Polynomial orders to start simulation with (increasing adaptation?)
      type(FTValueDictionary), intent(in) :: controlVariables !<  Input values
      !--------------------------------------
      integer              :: iEl      ! Element counter
      !--------------------------------------
      
!
!     ------------------------------
!     Read variables from input file
!     ------------------------------
!
      this % Adapt        = controlVariables % LogicalValueForKey("padaptation") ! Default false if not present
      
      if (this % Adapt) then
         this % Constructed = .TRUE.
      else
         this % Constructed = .FALSE.
         return
      end if
      
      this % increasing   = controlVariables % LogicalValueForKey("increasing adaptation")
      
      ! The truncation error threshold is required (only adaptation strategy implemented)
      if (.NOT. controlVariables % containsKey("truncation error")) ERROR STOP 'A truncation error must be specified for p-adapt'
      this % reqTE        = controlVariables % RealValueForKey   ("truncation error")
      
      this % controlVariables = controlVariables
      
      this % PlotInfo     = controlVariables % LogicalValueForKey("plot p-adaptation") ! Default false if not present
      if (this % PlotInfo) then
         this % plotFileName = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
         this % plotFileName = trim(this % plotFileName) // '_AdaptInfo.tec'
      end if
!
!     ------------------------------
!     Save maximum polynomial orders
!     ------------------------------
!
      
      if (controlVariables % containsKey("adaptation nmax i")) then
         this % NxyzMax(1) = controlVariables % integerValueForKey("adaptation nmax i")
      else
         this % NxyzMax(1) = MAXVAL(Nx)
      end if
      if (controlVariables % containsKey("adaptation nmax j")) then
         this % NxyzMax(2) = controlVariables % integerValueForKey("adaptation nmax j")
      else
         this % NxyzMax(2) = MAXVAL(Ny)
      end if
      if (controlVariables % containsKey("adaptation nmax k")) then
         this % NxyzMax(3) = controlVariables % integerValueForKey("adaptation nmax k")
      else
         this % NxyzMax(3) = MAXVAL(Nz)
      end if
      
!
!     ----------------------------------------
!     Allocate truncation error array
!     ----------------------------------------
!
      nelem = size(Nx)
      allocate (this % TE(nelem))
      
      do iEl = 1, nelem
         call this % TE(iEl) % construct(NMIN,this % NxyzMax)
      end do
      
!
!     ---------------------------------------------
!     If increasing adaptation is selected, rewrite
!     ---------------------------------------------
!
      if (this % increasing) then
         NInc = NInc_0
!$omp parallel do schedule(runtime)
         do iEl = 1, nelem
            if (Nx(iEl) > NInc) Nx(iEl) = NInc
            if (Ny(iEl) > NInc) Ny(iEl) = NInc
            if (Nz(iEl) > NInc) Nz(iEl) = NInc
         end do
!$omp end parallel do
      end if
      
   end subroutine ConstructPAdaptator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for destructing the p-adaptator
!  ----------------------------------------
   subroutine DestructPAdaptator(this)
      implicit none
      !--------------------------------------
      class(pAdaptation_t) :: this
      !--------------------------------------
      integer              :: iEl
      !--------------------------------------
      
      ! Truncation error
      do iEl = 1, nelem
         call this % TE(iEl) % destruct()
      end do
      
      deallocate (this % TE)
   end subroutine DestructPAdaptator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Main routine for adapting the polynomial order in all elements based on 
!  the truncation error estimation
!  ------------------------------------------------------------------------
   subroutine pAdaptTE(pAdapt,sem,itera,t)
      use AnisFASMultigridClass
      implicit none
      !--------------------------------------
      class(pAdaptation_t)       :: pAdapt            !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      !--------------------------------------
      integer                    :: iEl               !   Element counter
      integer                    :: iEQ               !   Equation counter
      integer                    :: Dir               !   Direction
      integer                    :: NMax              !   Max polynomial order
      integer                    :: p                 !   Polynomial order counter
      integer                    :: NNew(3,nelem)     !   New polynomial orders of mesh (after adaptation!)
      integer                    :: Error(3,nelem)    !   Stores (with ==1) elements where the truncation error has a strange behavior of the truncation error (in every direction)
      integer                    :: Warning(nelem)    !   Stores (with ==1) elements where the specified truncation error was not achieved
                                                      !   New polynomial orders(after adaptation)
      integer                    :: NOld(3)
      type(DGSem)                :: Oldsem
      logical                    :: success
      integer, save              :: Stage = 0         !   Stage of p-adaptation for the increasing method
      !--------------------------------------
      ! For new adaptation algorithm
      integer                    :: Pxyz(3)           !   Polynomial order the estimation was performed with
      real(kind=RP), allocatable :: TEmap(:,:,:)      !   Map of truncation error
      integer                    :: P_1(3)            !   Pxyz-1
      integer                    :: i,j,k             !   Counters
      integer                    :: DOFs, NewDOFs
      logical                    :: notenough(3)
      TYPE(AnisFASMultigrid_t)      :: AnisFASpAdaptSolver
      !--------------------------------------

      write(STD_OUT,*) '****     Performing p-Adaptation      ****'
      
      Error = 0
      Warning= 0
      Stage = Stage + 1
      
!
!     -------------------------------------------------------------
!     Estimate the truncation error using the anisotropic multigrid
!     -------------------------------------------------------------
!
      CALL AnisFASpAdaptSolver % construct(pAdapt % controlVariables,sem)
      CALL AnisFASpAdaptSolver % solve(itera,t,pAdapt % TE)
      CALL AnisFASpAdaptSolver % destruct
!
!     -------------------------------------------------------------
!     Find the polynomial order that fulfills the error requirement
!     -------------------------------------------------------------
!
      ! Loop over all elements
!!!!$OMP PARALLEL do PRIVATE(Pxyz,P_1,TEmap,NewDOFs,i,j,k)  ! TODO: Modify orivate and uncomment
      do iEl = 1, nelem
         
         Pxyz = sem % mesh % elements(iEl) % Nxyz
         P_1  = Pxyz - 1
         
!~         ! TODO: create routine to print this..
!~         !!if (iEl == 12) then
!~         !!   print*, 'Pxyz:',Pxyz
!~         !!   print*, 'P_1: ',P_1
!~         !!   print*, 'TEx: ',TE(iEl) % Dir(1) % maxTE
!~         !!   print*, 'TEy: ',TE(iEl) % Dir(2) % maxTE
!~         !!   print*, 'TEz: ',TE(iEl) % Dir(3) % maxTE
!~         !!end if
         
         where(P_1 < NMIN) P_1 = NMIN ! minimum order
         
         ! Allocate the TEmap with the maximum number of N compinations and initialize it
         allocate(TEmap(NMIN:pAdapt % NxyzMax(1),NMIN:pAdapt % NxyzMax(2),NMIN:pAdapt % NxyzMax(3)))
         TEmap = 0._RP
         
         NNew(:,iEl) = -1 ! Initialized to negative value
         
!
!        -----------------------------------------------------------------------
!        Check if the desired TE can be obtained using the "inner" map (N_i<P_i)
!           Accomplished in 3 merged steps:
!              1. Generate inner TE map
!              2. Check every entry of the map 
!              3. select the one that
!                 fulfills the requirement with the lowest DOFs
!        ----------------------------------------------------------------------
!
         NewDOFs = (P_1(1) + 1) * (P_1(2) + 1) * (P_1(3) + 1) ! Initialized with maximum number of DOFs (for N<P)
         
         do k = NMIN, P_1(3)
            do j = NMIN, P_1(2)
               do i = NMIN, P_1(1) 
                  ! 1. Generate a TEmap entry
                  TEmap(i,j,k) = pAdapt % TE(iEl) % Dir(1) % maxTE(i) + &  !xi   contribution
                                 pAdapt % TE(iEl) % Dir(2) % maxTE(j) + &  !eta  contribution
                                 pAdapt % TE(iEl) % Dir(3) % maxTE(k)      !zeta contribution
                  
                  ! 2. Check if it fulfills the requirement
                  if (TEmap(i,j,k) < pAdapt % reqTE) then
                     DOFs = (i+1) * (j+1) * (k+1)
                     
                     !  3. Select the entry if it minimizes the DOFs
                     if (DOFs <= NewDOFs) then
                        NNew(:,iEl) = [i,j,k]
                        NewDOFs = DOFs
                     end if
                  end if
               end do
            end do
         end do
!
!        -----------------------------------------------------------------------
!        If the desired TE could NOT be achieved using the inner map (N_i<P_i),
!        find it with a higher order
!           Accomplished in 3 steps:
!              1. Perform regression analysis and extrapolate the decoupled TE
!              2. Obtain the extended TE map for higher orders
!              3. Check every entry of the extended map
!              4. Select the N>P that fulfills the requirement with lowest DOFs
!        -----------------------------------------------------------------------
!
         if (any(NNew(:,iEl)<NMIN)) then
            
            ! 1. Regression analysis
            do Dir = 1, 3
               call RegressionIn1Dir(TE    = pAdapt % TE(iEl)     , &
                                     P_1   = P_1(Dir)             , & 
                                     NMax  = pAdapt % NxyzMax(Dir), &
                                     Stage = Stage                , &
                                     iEl   = iEl                  , &
                                     Dir   = Dir                  , &
                                     notenough = notenough(Dir)   , &
                                     error     = Error(Dir,iEl))
            end do
            
            ! If the truncation error behaves as expected, continue, otherwise skip steps 2-3-4 and select maximum N
            if (all(error(:,iEl) < 1)) then
               
               ! 2-3. Obtain extended TE map and search
               
               NewDOFs = (pAdapt % NxyzMax(1) + 1) * (pAdapt % NxyzMax(2) + 1) * (pAdapt % NxyzMax(3) + 1) ! Initialized as maximum DOFs possible
               do k = NMIN, pAdapt % NxyzMax(3)
                  do j = NMIN, pAdapt % NxyzMax(2)
                     do i = NMIN, pAdapt % NxyzMax(1)
                        ! cycle if it is not necessary/possible to compute the TEmap
                        if (k <= P_1(3) .AND. j <= P_1(2) .AND. i <= P_1(1)) cycle ! This is the inner map (already checked)
                        if (notenough(1) .AND. i > Pxyz(1)) cycle ! The regression was not possible in xi   (too few points), hence only checking with known values
                        if (notenough(2) .AND. j > Pxyz(2)) cycle ! The regression was not possible in eta  (too few points), hence only checking with known values
                        if (notenough(3) .AND. k > Pxyz(3)) cycle ! The regression was not possible in zeta (too few points), hence only checking with known values
                        
                        ! 2. Generate a TEmap entry
                        TEmap(i,j,k) = pAdapt % TE(iEl) % Dir(1) % maxTE(i) + &  !x contribution
                                       pAdapt % TE(iEl) % Dir(2) % maxTE(j) + &  !y contribution
                                       pAdapt % TE(iEl) % Dir(3) % maxTE(k)      !z contribution
                        
                        ! 3. Check if TE was achieved
                        if (TEmap(i,j,k) < pAdapt % reqTE) then
                           DOFs = (i+1) * (j+1) * (k+1)
                           ! 4. Select the entry if it minimizes the DOFs
                           if (DOFs <= NewDOFs) then
                              NNew(:,iEl) = [i,j,k]
                              NewDOFs = DOFs
                           end if
                        end if
                        
                     end do
                  end do
               end do
            else
               write(STD_OUT,*) 'p-Adaptation ERROR: Unexpected behavior of truncation error in element',iEl, '. Direction (1,2,3)', error(:,iEl)
               write(STD_OUT,*) '                    Using maximum polynomial order, N=', pAdapt % NxyzMax
               NNew(:,iEl) = pAdapt % NxyzMax
            end if
         end if
         
         ! If requested, print TEmap
!!!!          if (iEl == 3) call PrintTEmap(TEmap,iEl)
         
!
!        ---------------------------------------------------------------------------
!        If the maximum polynomial order was not found, select the maximum available
!        ---------------------------------------------------------------------------
!
         if (any(NNew(:,iEl)<NMIN)) then
            write(STD_OUT,*) 'p-Adaptation WARNING: Desired truncation error not achieved within specified limits in element', iEl
            write(STD_OUT,*) '                      Using max polynomial order instead, N=', pAdapt % NxyzMax
            Warning(iEl) = 1
            NNew(:,iEl) = pAdapt % NxyzMax
         end if
         
         deallocate(TEmap)
         
      end do
      
!!!!$OMP END PARALLEL DO
!
!     ----------------------------------------------------------------------------
!     In case of increasing adaptator, modify polynomial orders according to stage
!      And decide if it is necessary to continue adapting
!     ----------------------------------------------------------------------------
!
      if (pAdapt % increasing) then
!!          Stage = Stage + 1
!!          NInc = NInc + dN_Inc
         NInc = NInc * fN_Inc
         
         if (MAXVAL(NNew) > NInc) then
            where(NNew > NInc) NNew = NInc
!!             where(NNew > 1 .AND. NNew < NInc - dN_Inc) NNew = NInc - dN_Inc   ! This is new only for testing
         else
!!             Stage = 0
            pAdapt % Adapt = .FALSE.
         end if
         
      else ! Only adapt once
         pAdapt % Adapt = .FALSE.
      end if
      
      call ReorganizePolOrders(sem % mesh % faces,NNew)
      
!
!     -----------------------
!     Plot files if requested
!     -----------------------
!
!~      if (pAdapt % PlotInfo) call pAdapt % plot(sem,TE,Stage,NNew,Error,Warning)
!
!     ----------------------------------
!     Adapt sem to new polynomial orders
!     ----------------------------------
!
      Oldsem = sem
      call sem % destruct
      call sem % construct (  controlVariables  = pAdapt % controlVariables                              ,  &
                              externalState     = externalStateForBoundaryName,                             &
                              externalGradients = ExternalGradientForBoundaryName,                          &
                              Nx_ = NNew(1,:) ,     Ny_ = NNew(2,:),     Nz_ = NNew(3,:),                   &
                              success           = success)
      IF(.NOT. success)   ERROR STOP "Error constructing adapted DGSEM"
      
!
!     ------------------------------------
!     Save the solution in the adapted sem 
!     ------------------------------------
!
      ! Loop over all elements
      do iEl = 1, size(sem % mesh % elements)
         
         NOld = Oldsem % mesh % elements(iEl) % Nxyz
         
         ! Copy the solution if the polynomial orders are the same, if not, interpolate
         if (ALL(NOld == NNew(:,iEl))) then
            sem % mesh % elements(iEl) % storage % Q = Oldsem % mesh % elements(iEl) % storage % Q
         else
            
            !------------------------------------------------------------------
            ! Construct the interpolation matrices in every direction if needed
            !------------------------------------------------------------------
            call ConstructInterpolationMatrices( NOld(1),NNew(1,iEl) )  ! Xi
            call ConstructInterpolationMatrices( NOld(2),NNew(2,iEl) )  ! Eta
            call ConstructInterpolationMatrices( NOld(3),NNew(3,iEl) )  ! Zeta
            
            !---------------------------------------------
            ! Interpolate solution to new solution storage
            !---------------------------------------------
            
            call Interp3DArrays  (Nvars      = N_EQN                                       , &
                                  Nin        = NOld                                        , &
                                  inArray    = Oldsem % mesh % elements(iEl) % storage % Q , &
                                  Nout       = NNew(:,iEl)                                 , &
                                  outArray   = sem % mesh % elements(iEl) % storage % Q)
            
         end if
         
      end do
      
      call Oldsem % destruct
      
      ! Update residuals
      call ComputeTimeDerivative(sem, t)
      
      write(STD_OUT,*) 'p-Adaptation done, DOFs=', SUM((NNew(1,:)+1)*(NNew(2,:)+1)*(NNew(3,:)+1))
   end subroutine pAdaptTE
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for reorganizing the polynomial orders with neighbor contraints
!  -> So far, this only works with shared memory ...and computes everything in serial! (TODO: Update to MPI!)
!  -----------------------------------------------------------------------
   subroutine ReorganizePolOrders(faces,NNew)
      implicit none
      !------------------------------------------------------------
      type(Face), intent(in)    :: faces(:)
      integer   , intent(inout) :: NNew (:,:)
      !------------------------------------------------------------
      integer :: sweep
      integer :: fID
      integer :: eIDL            ! Element ID on the left
      integer :: eIDR            ! Element ID on the right
      integer :: indL(2)         ! Index of face polorders in 3D polorders (left element)
      integer :: indR(2)         ! Index of face polorders in 3D polorders (right element)
      integer :: NL,NR           ! Polynomial order on the left/right
      
      logical :: finalsweep
      !------------------------------------------------------------
!~      ! Vars for order perpendicular to face.
!~      integer :: indZL
!~      integer :: indZR
!~      integer :: NzL
!~      integer :: NzR
!~      !------------------------------------------------------------
      
      
      sweep = 0
      
      ! perform succesive (serial) sweeps until no further elements have to be modified 
      do 
         sweep = sweep + 1
         finalsweep = .TRUE. ! let's first assume this is the final sweep
         
         ! Modify the elements on both sides of each face according to their polynomial order (so far, this is performed in serial)
         do fID = 1, size(faces)
            
            !Cycle if this is a boundary face!!
            if (faces(fID) % elementIDs(2) == 0) cycle
            
            eIDL = faces(fID) % elementIDs(1)
            eIDR = faces(fID) % elementIDs(2)
            
            indL = axisMap(:, faces(fID) % elementSide(1))
            
            select case ( faces(fID) % rotation )
               case ( 0, 2, 5, 7 ) ! Local x and y axis are parallel or antiparallel
                  indR(1) = axisMap(1, faces(fID) % elementSide(2))
                  indR(2) = axisMap(2, faces(fID) % elementSide(2))
               case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
                  indR(2) = axisMap(1, faces(fID) % elementSide(2))
                  indR(1) = axisMap(2, faces(fID) % elementSide(2))
            end select
            
            !! Compare the polynomial order in the x-direction
            NL = NNew(indL(1),eIDL)
            NR = NNew(indR(1),eIDR)
            
            if (MIN(NL,NR) < GetOrderAcrossFace(MAX(NL,NR))) then
               finalsweep = .FALSE.
               if (NL<NR) then
                  NNew(indL(1),eIDL) = GetOrderAcrossFace(NR)
               else
                  NNew(indR(1),eIDR) = GetOrderAcrossFace(NL)
               end if
            end if
            
            !! Compare the polynomial order in the y-direction
            NL = NNew(indL(2),eIDL)
            NR = NNew(indR(2),eIDR)
            
            if (MIN(NL,NR) < GetOrderAcrossFace(MAX(NL,NR))) then
               finalsweep = .FALSE.
               if (NL<NR) then
                  NNew(indL(2),eIDL) = GetOrderAcrossFace(NR)
               else
                  NNew(indR(2),eIDR) = GetOrderAcrossFace(NL)
               end if
            end if
            
            !! Compare the polynomial order in the z-direction (needed???)
!             indZL = RemainingIndex(indL)
!             indZR = RemainingIndex(indR)
!             NzL = NNew(indZL,eIDL)
!             NzR = NNew(indZR,eIDR)
            
!             if (MIN(NzL,NzR) < GetOrderAcrossFace(MAX(NzL,NzR))) then
!                finalsweep = .FALSE.
!                if (NzL<NzR) then
!                   NNew(indZL,eIDL) = GetOrderAcrossFace(NzR) !(NyR*2)/3
!                else
!                   NNew(indZR,eIDR) = GetOrderAcrossFace(NzL) !(NyL*2)/3
!                end if
!             end if
            
         end do
         
         write(STD_OUT,*) 'Finishing "ReorganizePolOrders" sweep', sweep, '. Final = ', finalsweep
         if (finalsweep) exit
      end do
      
   end subroutine ReorganizePolOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
!  Subroutine to specify the minimum polynomial order of an element, when the element
!  accross the face is of order a. I.e., we specify the allowed jump of polynomial order
!  across a face.
!  -------------------------------------------------------------------------------------
   integer function GetOrderAcrossFace(a)
      integer :: a
      
!       GetOrderAcrossFace = (a*2)/3
       GetOrderAcrossFace = a-1
   end function GetOrderAcrossFace
!~!
!~!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!~!
!~   function RemainingIndex(a) RESULT(b)
!~      integer :: a(2)
!~      integer :: b
      
!~      if (any(a==1)) then
!~         if (any(a==2)) b = 3
!~         if (any(a==3)) b = 2
!~      else
!~         b = 1
!~      end if
      
!~   end function
!~!
!~!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!~!
!~!  -----------------------------------------------------------------------
!~!  Subroutine for printing the TE map(s) of one element
!~!  -----------------------------------------------------------------------
!~   subroutine PrintTEmap(TEmap,iEl)
!~      implicit none
!~      !-------------------------------------------
!~      real(kind=RP)  :: TEmap(NMIN:,NMIN:,NMIN:)
!~      integer        :: iEl
!~      !-------------------------------------------
!~      integer                :: k, i, l
!~      integer                :: fd
!~      character(LINE_LENGTH) :: TEmapfile
!~      !-------------------------------------------
      
!~      do k = NMIN, size(TEmap,3)
!~         write(TEmapfile,'(A,I4.4,A,I2.2,A)') 'RegressionFiles/TEmapXY-Elem_',iEl,'-Nz_',k,'.dat'
   
!~         open(newunit = fd, file=TRIM(TEmapfile), action='write')
!~            do i = NMIN, size(TEmap, 1)
!~               write(fd,*) (TEmap(i,l,k),l=NMIN,size(TEmap,2))
!~            end do
!~         close(fd)
!~      end do
!~   end subroutine PrintTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine that extrapolates the behavior of the directional components
!  of the truncation error.
!  -----------------------------------------------------------------------
   subroutine RegressionIn1Dir(TE,P_1,NMax,Dir,notenough,error,Stage,iEl)
      implicit none
      !---------------------------------------
      type(TruncationError_t)    :: TE                !<> Decaupled truncation error for one element
      integer                    :: P_1               !<  P-1 (max polynomial order with tau estimation for regression)
      integer                    :: NMax
      integer                    :: Dir
      logical                    :: notenough         !>  .TRUE. if there are not enough points in every direction for regression 
      integer                    :: error             !>  error=1 if line behavior is not as expected
      ! Additional arguments for IO
      integer                    :: Stage
      integer                    :: iEl
      !---------------------------------------
      
      real(kind=RP)              :: x   (P_1-NMIN+1)
      real(kind=RP)              :: y   (P_1-NMIN+1)
      integer                    :: N
      integer                    :: i
      real(kind=RP)              :: C,eta             ! Regression variables
      character(len=LINE_LENGTH) :: RegfileName
      integer                    :: fd
      !---------------------------------------
      
      ! Initializations
      error     = 0
      notenough = .FALSE.
      
      ! Check if there are enough points for regression
      if (P_1 < NMIN + 1) then
         notenough = .TRUE.
         return
      end if
      
      ! Perform regression analysis   
      N = P_1 - NMIN + 1
      y = LOG10(TE % Dir(Dir) % maxTE (NMIN:P_1))
      x = (/ (real(i,RP), i=NMIN,P_1) /)
      call C_and_eta_estimation(N,x,y,C,eta)
      
!~      ! Save files
!~      write(RegfileName,'(A,I2.2,A,I4.4,A,I1,A)') 'RegressionFiles/Stage_',Stage,'/Regression_Elem_',iEl,'Dir_',Dir,'.dat'
      
!~      open(newunit = fd, file=RegfileName, action='write')
!~         write(fd,*) x
!~         write(fd,*) y
!~         write(fd,*) C, eta
!~         write(fd,*) P_1+1, NMax
!~         write(fd,*) iEl, Dir
!~      close(fd)
      
      ! Check if there is an unexpected behavior
      if (C >= 0) then
         Error = 1 
         return
      end if
      
      ! Extrapolate the TE
      do i = P_1+1, NMax
         TE % Dir(Dir) % maxTE(i) = 10**(C + eta*i)
      end do
      
      
   end subroutine RegressionIn1Dir
!~!
!~!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!~!
!~!  -----------------------------------------------------------------------
!~!  Subroutine for performing least square regression and giving up the coefficients
!~!  -----------------------------------------------------------------------
!~   function ExtrapolateNTau(TE,NMIN,NMax,NMaxF,reqTE,iEl,Dir,Error,Warn,Stage) RESULT(NNew)
!~      implicit none
!~      !---------------------------------
!~      real(kind=RP)  :: TE(NMax-NMIN)
!~      integer        :: NMIN
!~      integer        :: NMax
!~      integer        :: NMaxF          ! Maximum allowed by user
!~      real(kind=RP)  :: reqTE
!~      integer        :: iEl            ! Number of element (only for info purposes)
!~      integer        :: Dir            ! Direction         (only for info purposes)
!~      integer        :: NNew           !>  New polynomial ordr
!~      Integer        :: Error
!~      integer        :: Warn
!~      integer        :: Stage
!~      !---------------------------------
!~      real(kind=RP)  :: x   (NMax-NMIN)
!~      real(kind=RP)  :: y   (NMax-NMIN)
!~      integer        :: N
!~      integer        :: i
!~      real(kind=RP)  :: C,eta             ! Regression variables
!~      integer        :: fd
!~      character(len=LINE_LENGTH) :: RegfileName
!~      !---------------------------------
      
!~      N = NMax - NMIN
!~      y = LOG10(TE)
!~      x = (/ (real(i,RP), i=NMIN,NMax-1) /)
      
!~      call C_and_eta_estimation(N,x,y,C,eta)
      
!~      !!!!!!! New here **************************
      
!~      if (C > 0) then
!~         print*, 'ERROR: Unexpected behavior of truncation error... This is very bad. El =', iEl,'Dir =', Dir
!~         NNew = NMaxF
!~         print*, '        ** Leaving with', NNew
!~         Error = 1
!~      else
!~         NNew = INT((LOG10(reqTE) - C)/eta) + 1
!~         if (NNew > NMaxF) then
!~            print*, 'Warning: Requested truncation error not found within specified limits. El =', iEl,'Dir =', Dir
!~            NNew = NMaxF
!~            print*, '        ** Leaving with', NNew
!~            Warn = 1
!~         else if (NNew < NMIN) then
!~            NNew = NMIN
!~         end if
!~      end if
      
!~      ! Save files
      
!~      write(RegfileName,'(A,I2.2,A,I4.4,A,I1,A)') 'RegressionFiles/Stage_',Stage,'/Regression_Elem_',iEl,'Dir_',Dir,'.dat'
      
!~      open(newunit = fd, file=RegfileName, action='write')
!~         write(fd,*) x
!~         write(fd,*) y
!~         write(fd,*) C, eta
!~         write(fd,*) NMax, NMaxF
!~         write(fd,*) iEl, Dir, NNew
!~      close(fd)
      
!~   end function
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for performing least square regression and giving up the coefficients
!  -----------------------------------------------------------------------
   pure subroutine C_and_eta_estimation(N,x,y,C,eta)
      implicit none
      !--------------------------------------
      integer      , intent(in)  :: N        !< Number of points
      real(kind=RP), intent(in)  :: x(N) 
      real(kind=RP), intent(in)  :: y(N)
      real(kind=RP), intent(out) :: C
      real(kind=RP), intent(out) :: eta
      !--------------------------------------
      real(kind=RP)              :: sumx,sumy,sumsqx,sumxy,deno
      !--------------------------------------
      
      sumx = sum(x)
      sumy = sum(y)
      sumsqx = sum(x*x)
      sumxy  = sum(x*y)
      
      deno=n*sumsqx-sumx*sumx
      
      eta = (n*sumxy-sumx*sumy)/deno
      C   = (sumsqx*sumy-sumx*sumxy)/deno
      
   end subroutine C_and_eta_estimation

!~!
!~!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!~!
!~!  -----------------------------------------------------------------------
!~!  Subroutine for plotting the adaptation information:
!~!     1. Plots the polynomial orders after adaptation
!~!     2. Plots the truncation error in every direction for the highest 
!~!        approximated order.
!~!  -----------------------------------------------------------------------
!~   subroutine AdaptationPlotting(this,sem,TE,Stage,NNew,Error,Warning)
!~      implicit none
!~      !--------------------------------------
!~      class(pAdaptation_t)   , intent(in) :: this      !<  p-Adaptator
!~      type(DGSem)            , intent(in) :: sem       !<  sem
!~      type(TruncationError_t), intent(in) :: TE(:)     !<  Type containing the element-wise truncation error estimation
!~      integer                , intent(in) :: Stage     !<  Stage of the adaptation process
!~      integer                , intent(in) :: NNew(:,:) !<  New polynomial orders
!~      integer                , intent(in) :: Error(:,:)
!~      integer                , intent(in) :: Warning(:,:)
!~      !--------------------------------------
!~      integer                            :: nelem
!~      integer                            :: fd
!~      integer                            :: id
!~      integer                            :: N(3)      ! Polynomial order in element
!~      integer                            :: i,j,k
!~      real(kind=RP)                      :: Tau(3)    ! Truncation error for plotting (L_inf per element in the finest approximation level)
!~      character(len=LINE_LENGTH)         :: plotFileName
!~      !--------------------------------------
      
!~      write(plotFileName,'(A,A,I3.3,A)') trim(this % plotFileName),'_AdaptInfo_Stage_', Stage, '.tec'
      
!~      open(newunit = fd, file=plotFileName, action='write')
      
!~         Nelem = size(sem % mesh % elements)
         
!~         write(fd,*) 'TITLE = "p-Adaptation information (HORSES3D)"'
!~      write(fd,*) 'VARIABLES = "X","Y","Z","TauX","TauY","TauZ","ErrorX","ErrorY","ErrorZ","WarnX","WarnY","WarnZ","Nx","Ny","Nz" '
         
!~         do id = 1, Nelem
!~            N = sem % mesh % elements(id) % Nxyz
!~            write(fd,*) 'ZONE I=', N(1)+1, ",J=",N(2)+1, ",K=",N(3)+1,", F=POINT"
            
!~            do i = 1, 3
!~               if (N(i) <= NMIN) then
!~                  Tau(i) = TE(id) % Dir(i) % maxTE(N(i))
!~               else
!~                  Tau(i) = TE(id) % Dir(i) % maxTE(N(i)-1)
!~               end if
!~            end do
            
!~            !-------
!~            ! Plot!
!~            !-------
!~            do k = 0, N(3)
!~               do j= 0, N(2)
!~                  do i = 0, N(1)
!~                     write(fd,'(3E13.5,x,3E13.5,x,3I10,x,3I10,x,3I10)') &
!~                                sem % mesh % elements(id) % geom % x(1,i,j,k), &
!~                                sem % mesh % elements(id) % geom % x(2,i,j,k), &
!~                                sem % mesh % elements(id) % geom % x(3,i,j,k), &
!~                                Tau,Error(:,id),Warning(:,id), NNew(:,id)
!~                  end do
!~               end do
!~            enddo
!~         end do
      
!~      close(fd)
      
!~   end subroutine AdaptationPlotting
   
end module pAdaptationClass
