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
   use InterpolationMatrices, only: Tset, Interp3DArrays, ConstructInterpolationMatrices
   use PhysicsStorage       , only: NCONS, N_GRAD_EQN
   use FaceClass            , only: Face
   use ElementClass
   use DGSEMClass
   use TruncationErrorClass
   use FTValueDictionaryClass
   use StorageClass
   use SharedBCModule
   use FileReadingUtilities , only: RemovePath, getFileName, getArrayFromString
   use FileReaders          , only: ReadOrderFile
   use ParamfileRegions     , only: readValueInRegion, getSquashedLine
   use HexMeshClass         , only: HexMesh
#if defined(CAHNHILLIARD)
   use BoundaryConditionFunctions, only: C_BC, MU_BC
#endif
   
   implicit none
   
#include "Includes.h"
   private
   public GetMeshPolynomialOrders, ReadOrderFile
   public pAdaptation_t
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type :: overenriching_t
      integer           :: ID
      integer           :: order
      real(kind=RP)     :: x_span(2)
      real(kind=RP)     :: y_span(2)
      real(kind=RP)     :: z_span(2)
      
      contains
         procedure      :: initialize => OverEnriching_Initialize
   end type overenriching_t
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type :: pAdaptation_t
      character(len=LINE_LENGTH)        :: solutionFileName ! Name of file for plotting adaptation information
      real(kind=RP)                     :: reqTE            ! Requested truncation error
      logical                           :: RegressionFiles  ! Write regression files?
      logical                           :: saveGradients    ! Save gradients in pre-adapt and p-adapted solution files?
      logical                           :: PlotInfo
      logical                           :: Adapt            ! Is the adaptator going to be used??
      logical                           :: increasing       ! Performing an increasing adaptation procedure?
      logical                           :: Constructed      ! 
      integer                           :: NxyzMax(3)       ! Maximum polynomial order in all the directions
      integer                           :: TruncErrorType   ! Truncation error type (either ISOLATED_TE or NON_ISOLATED_TE)
      type(TruncationError_t), allocatable :: TE(:)         ! Truncation error for every element(:)
      
      type(overenriching_t)  , allocatable :: overenriching(:)
      
      contains
         procedure :: construct => ConstructPAdaptator
         procedure :: destruct  => DestructPAdaptator
!~         procedure :: plot      => AdaptationPlotting
         procedure :: pAdaptTE
   end type pAdaptation_t
!
!  ----------------
!  Module variables
!  ----------------
!
   !! Variables
   integer    :: NMIN = 1
   integer    :: NInc_0 = 4
!!    integer               :: dN_Inc = 3 
   integer    :: fN_Inc = 2
   integer    :: NInc
   integer    :: nelem           ! number of elements in mesh
   
#if defined(NAVIERSTOKES)
   procedure(BCState_FCN)   :: externalStateForBoundaryName
   procedure(BCGradients_FCN)   :: ExternalGradientForBoundaryName
#elif defined(CAHNHILLIARD)
   procedure(BCState_FCN)   :: externalStateForBoundaryName
   procedure(BCGradients_FCN)   :: ExternalChemicalPotentialGradientForBoundaryName
   procedure(BCGradients_FCN)   :: ExternalConcentrationGradientForBoundaryName
#endif
   
   ! Here we define the input variables that can be changed after p-adaptation
   character(len=18), parameter :: ReplacedInputVars(4) = (/'mg sweeps         ', &
                                                            'mg sweeps pre     ', &
                                                            'mg sweeps post    ', &
                                                            'mg sweeps coarsest'/)
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
      integer              :: nelem
      integer, allocatable :: Nx_r(:), Ny_r(:), Nz_r(:)  
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
      
!
!     ********************************************************
!     Set maximum polynomial order for NodalStorage allocation
!     ********************************************************
!
      Nmax = 0
      
      ! Adaptation
      if (controlVariables % containsKey("adaptation nmax i")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax i"))
      if (controlVariables % containsKey("adaptation nmax j")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax j"))
      if (controlVariables % containsKey("adaptation nmax k")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax k"))
      
      ! Restart polynomial order
      if (controlVariables % containsKey("restart polorder" )) Nmax = max(Nmax,controlVariables % integerValueForKey("restart polorder" ))
      
      if (controlVariables % containsKey("restart polorder file" )) then
         call ReadOrderFile( controlVariables % stringValueForKey("restart polorder file", requestedLength = LINE_LENGTH), &
                             Nx_r, Ny_r, Nz_r )
         
         Nmax = max(Nmax,maxval(Nx_r),maxval(Ny_r),maxval(Nz_r))
      end if
      
      Nmax = max(Nmax,maxval(Nx),maxval(Ny),maxval(Nz))
      
      end subroutine GetMeshPolynomialOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ROUTINES FOR ADAPTATION
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine OverEnriching_Initialize(this,oID)
         implicit none
         !----------------------------------
         class(overenriching_t) :: this
         integer, intent(in)    :: oID
         !----------------------------------
         character(LINE_LENGTH) :: paramFile
         character(LINE_LENGTH) :: in_label
         character(LINE_LENGTH) :: x_span
         character(LINE_LENGTH) :: y_span
         character(LINE_LENGTH) :: z_span
         integer, allocatable   :: order
         !----------------------------------
         
         call get_command_argument(1, paramFile)
!
!        Get overenriching region ID
!        ---------------------------
         this % ID = oID
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define overenriching box " , this % ID
         
         call get_command_argument(1, paramFile) !
         call readValueInRegion ( trim ( paramFile )  , "order"  , order       , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "x span" , x_span      , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "y span" , y_span      , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "z span" , z_span      , in_label , "# end" ) 
         
         if (allocated(order)) then
            this % order = order
         else
            this % order = 1
         end if
         
         this % x_span = getArrayFromString(x_span)
         this % y_span = getArrayFromString(y_span)
         this % z_span = getArrayFromString(z_span)
         
      end subroutine OverEnriching_Initialize
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine getNoOfOverEnrichingBoxes(no_of_regions)
      implicit none
      integer, intent(out)    :: no_of_regions
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
!
!     Initialize
!     ----------
      no_of_regions = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find monitors
!     ------------------------------------
readloop:do 
         read ( fID , '(A)' , iostat = io ) line

         if ( io .lt. 0 ) then
!
!           End of file
!           -----------
            line = ""
            exit readloop

         elseif ( io .gt. 0 ) then
!
!           Error
!           -----
            errorMessage(STD_OUT)
            stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineoverenrichingbox' ) .gt. 0 ) then
               no_of_regions = no_of_regions + 1
            end if
            
         end if
         
      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

   end subroutine getNoOfOverEnrichingBoxes
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine OverEnrichRegions(overenriching,mesh,NNew,Nmax)
      implicit none
      !---------------------------------------
      type(overenriching_t), allocatable :: overenriching(:)
      type(HexMesh), intent(in)          :: mesh
      integer                            :: NNew(:,:)
      integer      , intent(in)          :: Nmax(3)
      !---------------------------------------
      integer :: oID       ! Overenriching region counter
      integer :: eID       ! Element counter
      integer :: cornerID  ! Corner counter
      logical :: enriched(mesh % no_of_elements)
      !---------------------------------------
      
      if (.not. allocated(overenriching) ) return
      
      enriched = .FALSE.
      
      do oID = 1, size(overenriching)
         associate (region => overenriching(oID) )
         
         element_loop: do eID=1, mesh % no_of_elements
            
            if ( enriched(eID) ) cycle element_loop
            
            associate ( corners => mesh % elements(eID) % hexMap % corners )
            
!
!           Enrich element if any of the corners is inside the region
!           ---------------------------------------------------------
            corner_loop: do cornerID=1, 8
               if( (corners(1,cornerID) > region % x_span(1) .and. corners(1,cornerID) < region % x_span(2)) .and. &
                   (corners(2,cornerID) > region % y_span(1) .and. corners(2,cornerID) < region % y_span(2)) .and. &
                   (corners(3,cornerID) > region % z_span(1) .and. corners(3,cornerID) < region % z_span(2)) ) then
                   
                  if ( NNew(1,eID) + region % order >= Nmax(1) ) then
                     NNew(1,eID) = Nmax(1)
                  else
                     NNew(1,eID) = NNew(1,eID) + region % order
                  end if
                  if ( NNew(2,eID) + region % order >= Nmax(2) ) then
                     NNew(2,eID) = Nmax(2)
                  else
                     NNew(2,eID) = NNew(2,eID) + region % order
                  end if
                  if ( NNew(3,eID) + region % order >= Nmax(3) ) then
                     NNew(3,eID) = Nmax(3)
                  else
                     NNew(3,eID) = NNew(3,eID) + region % order
                  end if
                  
                  enriched(eID) = .TRUE.
                  exit corner_loop
               end if 
            end do corner_loop
            
            end associate
         end do element_loop
         
         end associate
      end do
      
   end subroutine OverEnrichRegions
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
      implicit none
      !--------------------------------------
      class(pAdaptation_t)                :: this             !>  P-Adaptator
      integer, DIMENSION(:)               :: Nx,Ny,Nz         !<> Input: polynomial orders as read from input files - Output: Polynomial orders to start simulation with (increasing adaptation?)
      type(FTValueDictionary), intent(in) :: controlVariables !<  Input values
      !--------------------------------------
      integer              :: i      ! Element counter
      integer              :: no_of_overen_boxes
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
      this % reqTE        = controlVariables % doubleprecisionValueForKey   ("truncation error")
      
      this % PlotInfo        = controlVariables % LogicalValueForKey("plot p-adaptation") ! Default false if not present
      this % RegressionFiles = controlVariables % LogicalValueForKey("write regression files") ! Default false if not present
      
      this % solutionFileName = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
      this % saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
      
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
      if ( controlVariables % containsKey("truncation error type") ) then
         select case ( trim(controlVariables % stringValueForKey("truncation error type", requestedLength = LINE_LENGTH)) )
            case ('isolated')
               this % TruncErrorType = ISOLATED_TE
            case ('non-isolated')
               this % TruncErrorType = NON_ISOLATED_TE
            case default
               write(STD_OUT,*) "Not recognized 'truncation error type'. Defaulting to isolated."
               this % TruncErrorType = ISOLATED_TE
         end select
      else
         this % TruncErrorType = ISOLATED_TE
      end if
      
      nelem = size(Nx)
      allocate (this % TE(nelem))
      
!
!     If this is a p-anisotropic 3D case, the minimum polynomial order is 2
!     ---------------------------------------------------------------------
      
      if ( controlVariables % containsKey("adaptation nmin") ) then
         NMIN = controlVariables % integerValueForKey("adaptation nmin")
      else
         NMIN = 1
      end if
      
      do i = 1, nelem
         call this % TE(i) % construct(NMIN,this % NxyzMax)
      end do
      
      
!     Adaptation overenriching
!     *************************
      
      call getNoOfOverEnrichingBoxes(no_of_overen_boxes)
         
      if (no_of_overen_boxes > 0) then
         allocate ( this % overenriching(no_of_overen_boxes) )
         
         do i = 1, no_of_overen_boxes
            call this % overenriching(i) % initialize (i)
         end do
      end if
!
!     ---------------------------------------------
!     If increasing adaptation is selected, rewrite
!     ---------------------------------------------
!
      if (this % increasing) then
         NInc = NInc_0
!$omp parallel do schedule(runtime)
         do i = 1, nelem
            if (Nx(i) > NInc) Nx(i) = NInc
            if (Ny(i) > NInc) Ny(i) = NInc
            if (Nz(i) > NInc) Nz(i) = NInc
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
   subroutine pAdaptTE(pAdapt,sem,itera,t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
      use AnisFASMultigridClass
      use StopwatchClass
      implicit none
      !--------------------------------------
      class(pAdaptation_t)       :: pAdapt            !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      procedure(ComputeQDot_FCN) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables  !<> Input vaiables (that can be modified depending on the user input)
      !--------------------------------------
      integer                    :: iEl               !   Element counter
      integer                    :: iEQ               !   Equation counter
      integer                    :: Dir               !   Direction
      integer                    :: NMax              !   Max polynomial order
      integer                    :: p                 !   Polynomial order counter
      integer                    :: NNew(3,nelem)     !   New polynomial orders of mesh (after adaptation!)
      integer                    :: Error(3,nelem)    !   Stores (with ==1) elements where the truncation error has a strange behavior of the truncation error (in every direction)
      integer                    :: Warning(nelem)    !   Stores (with ==1) elements where the specified truncation error was not achieved
      integer                    :: NOld(3,nelem)     !   Old polynomial orders of mesh
      type(ElementStorage_t), allocatable :: TempStorage(:) ! Temporary variable to store the solution before the adaptation procedure 
      logical                    :: success
      integer, save              :: Stage = 0         !   Stage of p-adaptation for the increasing method
      CHARACTER(LEN=LINE_LENGTH) :: newInput          !   Variable used to change the input in controlVariables after p-adaptation 
      !--------------------------------------
      ! For new adaptation algorithm
      integer                    :: Pxyz(3)           !   Polynomial order the estimation was performed with
      real(kind=RP), allocatable :: TEmap(:,:,:)      !   Map of truncation error
      integer                    :: P_1(3)            !   Pxyz-1
      integer                    :: i,j,k             !   Counters
      integer                    :: DOFs, NewDOFs
      logical                    :: notenough(3)
      TYPE(AnisFASMultigrid_t)   :: AnisFASpAdaptSolver
      character(len=LINE_LENGTH) :: AdaptedMeshFile
      type(BCFunctions_t)        :: BCFunctions(no_of_BCsets)
      logical                    :: last
      interface
         character(len=LINE_LENGTH) function getFileName( inputLine )
            use SMConstants
            implicit none
            character(len=*)     :: inputLine
         end function getFileName
      end interface
      !--------------------------------------
      
      write(STD_OUT,*)
      write(STD_OUT,*)
      select case (pAdapt % TruncErrorType)
         case (ISOLATED_TE)
            write(STD_OUT,*) '****     Performing p-Adaptation with isolated truncation error estimates      ****'
         case (NON_ISOLATED_TE)
            write(STD_OUT,*) '****     Performing p-Adaptation with non-isolated truncation error estimates      ****'
         case default
            error stop ':: Non-defined truncation error type.'
      end select
      write(STD_OUT,*)
      
      Error = 0
      Warning= 0
      Stage = Stage + 1
      
!
!     -------------------------------------------
!     Get exact truncation error map if requested
!     -------------------------------------------
!
      if ( controlVariables % containsKey("get exact temap elem") ) then
         iEl = controlVariables % integerValueForKey("get exact temap elem")
         call GenerateExactTEmap(sem, NMIN, pAdapt % NxyzMax, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, iEl, pAdapt % TruncErrorType)
      end if
      
!
!     --------------------------------------
!     Write pre-adaptation mesh and solution
!     --------------------------------------
!
      write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( pAdapt % solutionFileName ), '_pre-Adapt_Stage_', Stage, '.hsol'
      call sem % mesh % Export(AdaptedMeshFile)
      
      call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),pAdapt % saveGradients)
!
!     -------------------------------------------------------------
!     Estimate the truncation error using the anisotropic multigrid
!     -------------------------------------------------------------
!
      CALL AnisFASpAdaptSolver % construct(controlVariables,sem,estimator=.TRUE.,NMINestim = NMIN)
      CALL AnisFASpAdaptSolver % solve(itera,t,computeTimeDerivative,ComputeTimeDerivativeIsolated,pAdapt % TE, pAdapt % TruncErrorType)
      CALL AnisFASpAdaptSolver % destruct
!
!     -------------------------------------------------------------
!     Find the polynomial order that fulfills the error requirement
!     -------------------------------------------------------------
!
      ! Allocate the TEmap with the maximum number of N compinations and initialize it
      allocate(TEmap(NMIN:pAdapt % NxyzMax(1),NMIN:pAdapt % NxyzMax(2),NMIN:pAdapt % NxyzMax(3)))
      
      ! Loop over all elements
!!!!$OMP PARALLEL do PRIVATE(Pxyz,P_1,TEmap,NewDOFs,i,j,k)  ! TODO: Modify private and uncomment
      do iEl = 1, nelem
         
         Pxyz = sem % mesh % elements(iEl) % Nxyz
         P_1  = Pxyz - 1
         
         where(P_1 < NMIN) P_1 = NMIN ! minimum order
         
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
               call RegressionIn1Dir(Adaptator = pAdapt               , &
                                     P_1   = P_1(Dir)             , & 
                                     NMax  = pAdapt % NxyzMax(Dir), &
                                     Stage = Stage                , &
                                     iEl   = iEl                  , &
                                     Dir   = Dir                  , &
                                     notenough = notenough(Dir)   , &
                                     error     = Error(Dir,iEl)   )
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
         
!~         if (iEl==1) call PrintTEmap(TEmap,1,NMIN,pAdapt % solutionFileName)
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
         
      end do
!!!!$OMP END PARALLEL DO
      
      deallocate(TEmap)
      
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
         else
            pAdapt % Adapt = .FALSE.
         end if
         
      else ! Only adapt once
         pAdapt % Adapt = .FALSE.
      end if
      
      call OverEnrichRegions(pAdapt % overenriching,sem % mesh,NNew, pAdapt % NxyzMax)
      
      last = .FALSE.
      do while (.not. last)
         last = .TRUE.
         call makeBoundariesPConforming(sem % mesh,NNew,last)
         call ReorganizePolOrders(sem % mesh % faces,NNew,last)
      end do
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
      call Stopwatch % Pause("Solver")
      call Stopwatch % Start("Preprocessing")
      
#if defined(NAVIERSTOKES)
      BCFunctions(1) % externalState => externalStateForBoundaryName
      BCFunctions(1) % externalGradients => externalGradientForBoundaryName
#elif defined(CAHNHILLIARD)
      BCFunctions(C_BC) % externalState      => externalStateForBoundaryName
      BCFunctions(C_BC) % externalGradients  => externalConcentrationGradientForBoundaryName

      BCFunctions(MU_BC) % externalState     => externalStateForBoundaryName
      BCFunctions(MU_BC) % externalGradients => externalChemicalPotentialGradientForBoundaryName
#endif
!
!     ---------------------------
!     Store the previous solution
!     ---------------------------
!
      allocate (TempStorage(nelem))
      do iEl = 1, nelem
         NOld (:,iEl) = sem % mesh % elements(iEl) % Nxyz
         call TempStorage(iEl) % Construct(NOld (1,iEl), NOld (2,iEl), NOld (3,iEl), NCONS, N_GRAD_EQN, .FALSE.)
         TempStorage(iEl) % Q = sem % mesh % elements(iEl) % storage % Q
      end do
!
!     -----------------
!     Construct new sem
!     -----------------
!
      call sem % destruct
      call sem % construct (  controlVariables  = controlVariables                       ,   &
                              BCFunctions       = BCFunctions, &
                              Nx_ = NNew(1,:) ,     Ny_ = NNew(2,:),     Nz_ = NNew(3,:),    &
                              success           = success)
      IF(.NOT. success)   ERROR STOP "Error constructing adapted DGSEM"
      call Stopwatch % Pause("Preprocessing")
      call Stopwatch % Start("Solver")
      
!
!     ------------------------------------
!     Save the solution in the adapted sem 
!     ------------------------------------
!
      ! Loop over all elements
      do iEl = 1, size(sem % mesh % elements)
         
         ! Copy the solution if the polynomial orders are the same, if not, interpolate
         if (ALL(NOld(:,iEl) == NNew(:,iEl))) then
            sem % mesh % elements(iEl) % storage % Q = TempStorage(iEl) % Q
         else
            
            !------------------------------------------------------------------
            ! Construct the interpolation matrices in every direction if needed
            !------------------------------------------------------------------
            call ConstructInterpolationMatrices( NOld(1,iEl),NNew(1,iEl) )  ! Xi
            call ConstructInterpolationMatrices( NOld(2,iEl),NNew(2,iEl) )  ! Eta
            call ConstructInterpolationMatrices( NOld(3,iEl),NNew(3,iEl) )  ! Zeta
            
            !---------------------------------------------
            ! Interpolate solution to new solution storage
            !---------------------------------------------
            
            call Interp3DArrays  (Nvars      = NCONS                                       , &
                                  Nin        = NOld(:,iEl)                                 , &
                                  inArray    = TempStorage(iEl) % Q , &
                                  Nout       = NNew(:,iEl)                                 , &
                                  outArray   = sem % mesh % elements(iEl) % storage % Q)
            
         end if
         
      end do

!
!     ------------------------------
!     Destruct the temporary storage
!     ------------------------------
!
      do iEl = 1, nelem
         call TempStorage(iEl) % destruct
      end do
      deallocate (TempStorage)
!
!     ---------------------------------------------------
!     Write post-adaptation mesh, solution and order file
!     ---------------------------------------------------
!
      write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( pAdapt % solutionFileName ), '_p-Adapted_Stage_', Stage, '.hsol'
      call sem % mesh % Export(AdaptedMeshFile)
      call sem % mesh % ExportOrders(AdaptedMeshFile)
      
      call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),pAdapt % saveGradients)
      
!
!     ------------------------------------------------
!     Rewrite controlVariables according to user input
!     ------------------------------------------------
!
      do i = 1, size(ReplacedInputVars)
         if ( controlVariables % containsKey( "padapted " // trim(ReplacedInputVars(i)) ) ) then
            newInput = controlVariables % stringValueForKey("padapted " // trim(ReplacedInputVars(i)), requestedLength = LINE_LENGTH)
            call controlVariables % removeObjectForKey( trim(ReplacedInputVars(i)) )
            call controlVariables % addValueForKey( TRIM(newInput) , trim(ReplacedInputVars(i)) )
         end if
      end do
!
!
!
!     ----------------
!     Update residuals
!     ----------------
!
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, sem % BCFunctions)
      
      write(STD_OUT,*) '****    p-Adaptation done, DOFs=', SUM((NNew(1,:)+1)*(NNew(2,:)+1)*(NNew(3,:)+1)), '****'
   end subroutine pAdaptTE
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to make the p-representation on certain boundaries conforming
!  -----------------------------------------------------------------------
   subroutine makeBoundariesPConforming(mesh,NNew,last)
      use ElementConnectivityDefinitions
      implicit none
      !------------------------------------------------------------
      type(HexMesh), intent(in)    :: mesh
      integer      , intent(inout) :: NNew (:,:)
      logical :: last
      !------------------------------------------------------------
      integer :: zoneID    ! Zone counters
      integer :: fID       ! Index of face on boundary (in partition)
      integer :: fIdx      ! Index of face on boundary (in zone)
      integer :: eID       ! Index of element on bountary
      integer :: eSide     ! Side of element on boundary
      integer :: f_conf
      integer :: sweep     ! Sweep counter
      logical :: finalsweep
      character(len=LINE_LENGTH), allocatable :: boundaryNames(:)
      !------------------------------------------------------------
      
      allocate ( boundaryNames( conformingBoundariesDic % COUNT() )  ) 
      boundaryNames = conformingBoundariesDic % allKeys()
      
      if ( size(boundaryNames) < 1 ) return
      
      write(STD_OUT,*) '## Forcing p-conforming boundaries ##'
      
!     ************************
!     Loop over the boundaries
!     ************************
      do zoneID = 1, size(mesh % zones)
         
         if ( all ( boundaryNames /= trim(mesh % zones(zoneID) % Name) ) ) cycle
         
         write(STD_OUT,*) '## Boundary:', mesh % zones(zoneID) % Name
         write(STD_OUT,*) '   Sweep   |   Last'
         sweep = 0
         
         ! Perform a number of sweeps until the representation is conforming
         do
            sweep = sweep + 1
            finalsweep = .TRUE. ! let's first assume this is the final sweep
            
            ! loop over the faces on every boundary
            do fIdx = 1, mesh % zones(zoneID) % no_of_faces
               
               fID   = mesh % zones(zoneID) % faces(fIdx)
               eID   = mesh % faces(fID) % elementIDs(1)
               eSide = mesh % faces(fID) % elementSide(1)
               
               ! loop over the faces that are shares between boundary elements
               do f_conf = 1, 4
                  
                  associate (f => mesh % faces ( mesh % elements(eID) % faceIDs (neighborFaces(f_conf,eSide) ) ) )
                  
                  call AdjustPAcrossFace( f, NNew, SameNumber, finalsweep)
                  
                  end associate
               end do
               
            end do
            
            write(STD_OUT,'(I10,X,A,X,L)') sweep ,'|', finalsweep 
            
            if (finalsweep) then
               if (sweep > 1) last = .FALSE.
               exit
            end if
         end do
      end do
      
      write(STD_OUT,*) '#####################################'
      
   end subroutine makeBoundariesPConforming
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for reorganizing the polynomial orders with neighbor contraints
!  -> So far, this only works with shared memory ...and computes everything in serial! (TODO: Update to MPI!)
!  -----------------------------------------------------------------------
   subroutine ReorganizePolOrders(faces,NNew,last)
      implicit none
      !------------------------------------------------------------
      type(Face), intent(in)    :: faces(:)
      integer   , intent(inout) :: NNew (:,:)
      logical :: last
      !------------------------------------------------------------
      integer :: sweep
      integer :: fID
      logical :: finalsweep
      !------------------------------------------------------------
      
      sweep = 0
      
      ! perform succesive (serial) sweeps until no further elements have to be modified 
      do 
         sweep = sweep + 1
         finalsweep = .TRUE. ! let's first assume this is the final sweep
         
         ! Modify the elements on both sides of each face according to their polynomial order (so far, this is performed in serial)
         do fID = 1, size(faces)
            
            !Cycle if this is a boundary face!!
            if (faces(fID) % elementIDs(2) == 0) cycle
            
            call AdjustPAcrossFace(faces(fID),NNew,GetOrderAcrossFace,finalsweep)
            
         end do
         
         write(STD_OUT,*) 'Finishing "ReorganizePolOrders" sweep', sweep, '. Final = ', finalsweep
         if (finalsweep) then
            if (sweep > 1) last = .FALSE.
            exit
         end if
      end do
      
   end subroutine ReorganizePolOrders
   
   subroutine AdjustPAcrossFace(f,NNew,OrderAcrossFace,same)
      implicit none
      !-arguments--------------------------------------------
      type(Face), intent(in)    :: f
      integer   , intent(inout) :: NNew(:,:)
      logical                   :: same   !<> Returns false if the function changes a polynomial order
      !------------------------------------------------------
      interface
         integer function OrderAcrossFace(a)
            integer :: a
         end function
      end interface
      !-local-variables--------------------------------------
      integer :: eIDL            ! Element ID on the left
      integer :: eIDR            ! Element ID on the right
      integer :: indL(2)         ! Index of face polorders in 3D polorders (left element)
      integer :: indR(2)         ! Index of face polorders in 3D polorders (right element)
      integer :: NL,NR           ! Polynomial order on the left/right
      !------------------------------------------------------
      
      eIDL = f % elementIDs(1)
      eIDR = f % elementIDs(2)
      
      if (eIDL < 1 .or. eIDR < 1) return
      
      indL = axisMap(:, f % elementSide(1))
      
      select case ( f % rotation )
         case ( 0, 2, 5, 7 ) ! Local x and y axis are parallel or antiparallel
            indR(1) = axisMap(1, f % elementSide(2))
            indR(2) = axisMap(2, f % elementSide(2))
         case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
            indR(2) = axisMap(1, f % elementSide(2))
            indR(1) = axisMap(2, f % elementSide(2))
      end select
      
      !! Compare the polynomial order in the x-direction
      NL = NNew(indL(1),eIDL)
      NR = NNew(indR(1),eIDR)
      
      if (MIN(NL,NR) < OrderAcrossFace(MAX(NL,NR))) then
         same = .FALSE.
         if (NL<NR) then
            NNew(indL(1),eIDL) = OrderAcrossFace(NR)
         else
            NNew(indR(1),eIDR) = OrderAcrossFace(NL)
         end if
      end if
      
      !! Compare the polynomial order in the y-direction
      NL = NNew(indL(2),eIDL)
      NR = NNew(indR(2),eIDR)
      
      if (MIN(NL,NR) < OrderAcrossFace(MAX(NL,NR))) then
         same = .FALSE.
         if (NL<NR) then
            NNew(indL(2),eIDL) = OrderAcrossFace(NR)
         else
            NNew(indR(2),eIDR) = OrderAcrossFace(NL)
         end if
      end if
      
   end subroutine AdjustPAcrossFace
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
      
      !GetOrderAcrossFace = (a*2)/3
      GetOrderAcrossFace = a-1
   end function GetOrderAcrossFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
   integer function SameNumber(a)
      implicit none
      integer :: a
      
      SameNumber = a
   end function SameNumber
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine that extrapolates the behavior of the directional components
!  of the truncation error.
!  -----------------------------------------------------------------------
   subroutine RegressionIn1Dir(Adaptator,P_1,NMax,Dir,notenough,error,Stage,iEl)
      implicit none
      !---------------------------------------
      type(pAdaptation_t)        :: Adaptator
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
      y = LOG10(Adaptator % TE(iEl) % Dir(Dir) % maxTE (NMIN:P_1))
      x = (/ (real(i,RP), i=NMIN,P_1) /)
      call C_and_eta_estimation(N,x,y,C,eta)
      
      ! Extrapolate the TE
      do i = P_1+1, NMax
         Adaptator % TE(iEl) % Dir(Dir) % maxTE(i) = 10**(C + eta*i)
      end do
      
      ! Write regression files
      if (Adaptator % RegressionFiles) then
         write(RegfileName,'(A,I2.2,3A,I7.7,A,I1,A)') 'RegressionFiles/Stage_',Stage,'/',trim(removePath(Adaptator % solutionFileName)),'_Elem_',iEl,'_Dir_',Dir,'.dat'
      
         open(newunit = fd, file=RegfileName, action='WRITE')
            WRITE(fd,*) x
            WRITE(fd,*) y
            WRITE(fd,*) C, eta
            WRITE(fd,*) NMax
            WRITE(fd,*) iEl, Dir
         close(fd)
      end if
      
      ! Check if there is an unexpected behavior
      if (eta >= 0) then
         Error = 1 
         return
      end if
      
   end subroutine RegressionIn1Dir
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
   
end module pAdaptationClass
