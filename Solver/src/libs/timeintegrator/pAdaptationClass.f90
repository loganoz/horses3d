!
!//////////////////////////////////////////////////////
!
!      Base class containing routines for adapting polynomial orders.
!        -> The current implementation is only valid for shared memory 
!           parallelization (OpenMP). TODO: Update to MPI.
!
!////////////////////////////////////////////////////////////////////////
!
module pAdaptationClass
   use SMConstants
   use FaceClass                       , only: Face
   use ElementClass
   use ElementConnectivityDefinitions  , only: axisMap
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f
   use FTValueDictionaryClass          , only: FTValueDictionary
   use StorageClass
   use FileReadingUtilities            , only: getRealArrayFromString, getIntArrayFromString
   use FileReaders                     , only: ReadOrderFile
   use ParamfileRegions                , only: readValueInRegion, getSquashedLine
   use HexMeshClass                    , only: HexMesh
   use ElementConnectivityDefinitions  , only: neighborFaces
   use ReadMeshFile                    , only: NumOfElemsFromMeshFile
   implicit none
   
#include "Includes.h"
   ! private
   ! public GetMeshPolynomialOrders, ReadOrderFile
   ! public pAdaptation_t, ADAPT_DYNAMIC_TIME, ADAPT_STATIC, ADAPT_DYNAMIC_TIME, NO_ADAPTATION
   ! public NInc_0, NInc, reorganize_Nz, GetOrderAcrossFace
   public
!
!  ----------------
!  Module parameters
!  ----------------
!
   integer, parameter :: ADAPT_STATIC = 0
   integer, parameter :: ADAPT_DYNAMIC_ITER = 1
   integer, parameter :: ADAPT_DYNAMIC_TIME = 2
   integer, parameter :: NO_ADAPTATION = 3
!
!  ----------------
!  Module variables
!  ----------------
!
   integer    :: NInc_0     = 4
   integer    :: NInc
   logical    :: reorganize_Nz = .FALSE.
   
   procedure(OrderAcrossFace_f), pointer :: GetOrderAcrossFace
   
   !-----------------------------------------
   ! Type for storing the overenriching areas
   !-----------------------------------------
   type :: overenriching_t
      integer           :: ID
      integer           :: order
      real(kind=RP)     :: x_span(2)
      real(kind=RP)     :: y_span(2)
      real(kind=RP)     :: z_span(2)
      
      contains
         procedure      :: initialize => OverEnriching_Initialize
   end type overenriching_t
   
   !-------------------------------------------------------
   ! Main base type for performing a p-adaptation procedure
   !-------------------------------------------------------
   type, abstract :: pAdaptation_t
      real(kind=RP)                     :: reqTE = 0.01_RP                 ! Requested truncation error
      character(len=LINE_LENGTH)        :: solutionFileName                ! Name of file for plotting adaptation information
      logical                           :: saveGradients                   ! Save gradients in pre-adapt and p-adapted solution files?
      logical                           :: saveSensor                      ! Save sensor in pre-adapt and p-adapted solution files?
      logical                           :: Adapt                           ! Is the adaptator going to be used??
      logical                           :: Constructed      ! 
      logical                           :: restartFiles    = .FALSE.
      logical                           :: UnSteady
      integer                           :: NxyzMax(3)                      ! Maximum polynomial order in all the directions
      integer                           :: adaptation_mode = NO_ADAPTATION ! Adaptation mode 
      real(kind=RP)                     :: time_interval
      integer                           :: iter_interval
      logical                           :: performPAdaptationT
      real(kind=RP)                     :: nextAdaptationTime = huge(1._RP)
      character(len=BC_STRING_LENGTH), allocatable :: conformingBoundaries(:) ! Stores the conforming boundaries (the length depends on FTDictionaryClass)
      type(overenriching_t)  , allocatable :: overenriching(:)
      
      contains
         procedure(constructInterface), deferred :: pAdaptation_Construct
         generic :: construct  => pAdaptation_Construct
         procedure(destructInterface), deferred :: pAdaptation_Destruct
         generic :: destruct   => pAdaptation_Destruct
         procedure, non_overridable :: makeBoundariesPConforming
         procedure(pAdaptInterface), deferred :: pAdaptation_pAdapt
         generic :: pAdapt   => pAdaptation_pAdapt
         procedure, non_overridable :: hasToAdapt => pAdaptation_hasToAdapt
   end type pAdaptation_t

!  ------------------------------------------
!  Interface for constructing the p-adaptator
!  ------------------------------------------
   interface
   subroutine constructInterface(this, controlVariables, t0)
      import pAdaptation_t
      import FTValueDictionary
      import RP
      !-------------------------------------------------
      implicit none
      !-------------------------------------------------
      class(pAdaptation_t)   , intent(inout) :: this             !>  P-Adaptator
      type(FTValueDictionary), intent(in)    :: controlVariables !<  Input values
      real(kind=RP)          , intent(in)    :: t0
   end subroutine constructInterface
   end interface

!  -----------------------------------------
!  Interface for destructing the p-adaptator
!  -----------------------------------------
   interface
   subroutine destructInterface(this)
      import pAdaptation_t
      !--------------------------------------
      implicit none
      !--------------------------------------
      class(pAdaptation_t) :: this
   end subroutine destructInterface
   end interface

!  -----------------------------------------------------------------
!  Main interface for adapting the polynomial order in all elements 
!  -----------------------------------------------------------------
   interface
   subroutine pAdaptInterface(this,sem,itera,t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
      import pAdaptation_t
      import DGSem
      import FTValueDictionary
      import RP
      import ComputeTimeDerivative_f
      !--------------------------------------
      implicit none
      !--------------------------------------
      class(pAdaptation_t)       :: this              !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables  !<> Input variables (that can be modified depending on the user input)
   end subroutine pAdaptInterface
   end interface 

!  -------------------
!  Abstract Interfaces
!  -------------------
   abstract interface
      pure integer function OrderAcrossFace_f(a)
         integer, intent(in) :: a
      end function
   end interface

!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     GENERAL ROUTINES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
      implicit none
      !-------------------------------------------------
      type(FTValueDictionary), intent(in)    :: controlVariables
      integer, allocatable   , intent(inout) :: Nx(:), Ny(:), Nz(:)  
      integer                , intent(out)   :: Nmax
      !-------------------------------------------------
      integer              :: nelem, i
      integer, allocatable :: Nx_r(:), Ny_r(:), Nz_r(:)  
      character(LINE_LENGTH)         :: paramFile
      character(LINE_LENGTH)         :: in_label
      character(LINE_LENGTH)         :: R_Nmax
      logical, allocatable           :: R_increasing
      !-------------------------------------------------
      
!
!     *************************************
!     Read the simulation polynomial orders
!     *************************************
!
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
      
!     If it is specified by the p-adaptation, read block
!     --------------------------------------------------
      
      write(in_label , '(A)') "#define p-adaptation"
      call get_command_argument(1, paramFile) !
      call readValueInRegion ( trim ( paramFile )  , "nmax"       , R_Nmax      , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "increasing" , R_increasing, in_label , "# end" ) 
      if ( R_Nmax /= "" ) then
         Nmax = maxval (getIntArrayFromString(R_Nmax))
      end if
      
!     Restart polynomial order
!     ------------------------
      
      if (controlVariables % containsKey("restart polorder" )) Nmax = max(Nmax,controlVariables % integerValueForKey("restart polorder" ))
      
      if (controlVariables % containsKey("restart polorder file" )) then
         call ReadOrderFile( controlVariables % stringValueForKey("restart polorder file", requestedLength = LINE_LENGTH), &
                             Nx_r, Ny_r, Nz_r )
         
         Nmax = max(Nmax,maxval(Nx_r),maxval(Ny_r),maxval(Nz_r))
         deallocate (Nx_r , Ny_r , Nz_r )
      end if
      
      Nmax = max(Nmax,maxval(Nx),maxval(Ny),maxval(Nz))
      
!
!     *****************************************************************************
!     Correct the simulation polynomial orders if increasing adaptation is selected
!     *****************************************************************************
!
      if ( allocated(R_increasing) ) then
         if (R_increasing) then
            NInc = NInc_0
!$omp parallel do schedule(runtime)
            do i = 1, size(Nx)
               if (Nx(i) > NInc) Nx(i) = NInc
               if (Ny(i) > NInc) Ny(i) = NInc
               if (Nz(i) > NInc) Nz(i) = NInc
            end do
!$omp end parallel do
         end if
      end if
      
   end subroutine GetMeshPolynomialOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ROUTINES FOR OVERENRICHING AREAS
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
         
         this % x_span = getRealArrayFromString(x_span)
         this % y_span = getRealArrayFromString(y_span)
         this % z_span = getRealArrayFromString(z_span)
         
         deallocate(order)
         
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
!     Read the whole file to find overenriching regions
!     -------------------------------------------------
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
            error stop "error stopped."

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
   pure function pAdaptation_hasToAdapt(this,iter) result(hasToAdapt)
      implicit none
      class(pAdaptation_t), intent(in) :: this
      integer             , intent(in) :: iter
      logical                          :: hasToAdapt
      
      select case (this % adaptation_mode)
         
         case (ADAPT_STATIC)
            hasToAdapt = .FALSE.
            
         case (ADAPT_DYNAMIC_ITER)
            if ( (mod(iter, this % iter_interval) == 0) .or. (iter == 1) ) then
               hasToAdapt = .TRUE.
            else
               hasToAdapt = .FALSE.
            end if
            
         case (ADAPT_DYNAMIC_TIME)
            hasToAdapt = this % performPAdaptationT
            
         case default
            hasToAdapt = .FALSE.
      end select
      
   end function pAdaptation_hasToAdapt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   logical function pAdaptationIsDefined()
      implicit none
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
      pAdaptationIsDefined = .FALSE.
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find overenriching regions
!     -------------------------------------------------
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
            error stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#definep-adaptation' ) .gt. 0 ) then
               pAdaptationIsDefined = .TRUE.
               close(fID)  
               return
            end if
            
         end if
         
      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

   end function pAdaptationIsDefined 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to make the p-representation on certain boundaries conforming
!  -----------------------------------------------------------------------
   subroutine makeBoundariesPConforming(this,mesh,NNew,last)
      implicit none
      !-arguments--------------------------------------------------
      class(pAdaptation_t), intent(inout) :: this            !<> Adaptation class
      type(HexMesh)       , intent(in)    :: mesh
      integer             , intent(inout) :: NNew (:,:)
      logical             , intent(inout) :: last
      !-local-variables--------------------------------------------
      integer :: zoneID    ! Zone counters
      integer :: fID       ! Index of face on boundary (in partition)
      integer :: fIdx      ! Index of face on boundary (in zone)
      integer :: eID       ! Index of element on bountary
      integer :: eSide     ! Side of element on boundary
      integer :: f_conf
      integer :: sweep     ! Sweep counter
      logical :: finalsweep
      !------------------------------------------------------------
      ! New definition of neighborFaces in order to consider the face
      ! that is not in contact with the boundary. In order to disable, 
      ! just comment the following
      integer, parameter :: neighborFaces(5,6) = reshape (  (/ 2, 3, 4, 5, 6, &
                                                               3, 4, 5, 6, 1, &
                                                               1, 2, 4, 5, 6, &
                                                               1, 2, 3, 5, 6, &
                                                               1, 2, 3, 4, 6, &
                                                               1, 2, 3, 4, 5 /) , (/5,6/) )
      !------------------------------------------------------------
      
      if ( .not. allocated(this % conformingBoundaries) ) return
      
      ! write(STD_OUT,*) '## Forcing p-conforming boundaries ##'
      
!     ************************
!     Loop over the boundaries
!     ************************
      do zoneID = 1, size(mesh % zones)
         
         if ( all ( this % conformingBoundaries /= trim(mesh % zones(zoneID) % Name) ) ) cycle
         
         ! write(STD_OUT,*) '## Boundary:', mesh % zones(zoneID) % Name
         ! write(STD_OUT,*) '   Sweep   |   Last'
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
               do f_conf = 1, size(neighborFaces,1)
                  
                  associate (f => mesh % faces ( mesh % elements(eID) % faceIDs (neighborFaces(f_conf,eSide) ) ) )
                  
                  call AdjustPAcrossFace( f, NNew, SameNumber, finalsweep)
                  
                  end associate
               end do
               
            end do
            
            ! write(STD_OUT,'(I10,X,A,X,L)') sweep ,'|', finalsweep 
            
            if (finalsweep) then
               if (sweep > 1) last = .FALSE.
               exit
            end if
         end do
      end do
      
      ! write(STD_OUT,*) '#####################################'
      
   end subroutine makeBoundariesPConforming
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for reorganizing the polynomial orders with neighbor constraints
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
      
      ! perform successive (serial) sweeps until no further elements have to be modified 
      do 
         sweep = sweep + 1
         finalsweep = .TRUE. ! let's first assume this is the final sweep
         
         ! Modify the elements on both sides of each face according to their polynomial order (so far, this is performed in serial)
         do fID = 1, size(faces)
            
            !Cycle if this is a boundary face!!
            if (faces(fID) % elementIDs(2) == 0) cycle
            
            call AdjustPAcrossFace(faces(fID),NNew,GetOrderAcrossFace,finalsweep,reorganize_Nz)
            
         end do
         
         write(STD_OUT,*) 'Finishing "ReorganizePolOrders" sweep', sweep, '. Final = ', finalsweep
         if (finalsweep) then
            if (sweep > 1) last = .FALSE.
            exit
         end if
      end do
      
   end subroutine ReorganizePolOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine AdjustPAcrossFace(f,NNew,OrderAcrossFace,same,z_dir)
      implicit none
      !-arguments--------------------------------------------
      type(Face), intent(in)    :: f
      integer   , intent(inout) :: NNew(:,:)
      logical                   :: same   !<> Returns false if the function changes a polynomial order
      logical   , optional      :: z_dir
      !------------------------------------------------------
      procedure(OrderAcrossFace_f) :: OrderAcrossFace
      !-local-variables--------------------------------------
      integer :: eIDL            ! Element ID on the left
      integer :: eIDR            ! Element ID on the right
      integer :: indL(2)         ! Index of face polorders in 3D polorders (left element)
      integer :: indR(2)         ! Index of face polorders in 3D polorders (right element)
      integer :: NL,NR           ! Polynomial order on the left/right
      integer :: indZL, indZR
      logical :: direction_z
      !------------------------------------------------------
      
      if (present(z_dir) ) then
         direction_z = z_dir
      else
         direction_z = .FALSE.
      end if
      
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
      
      if ( direction_z ) then
         !! Compare the polynomial order in the z-direction (needed???) 
         indZL = RemainingIndex(indL) 
         indZR = RemainingIndex(indR) 
         NL = NNew(indZL,eIDL) 
         NR = NNew(indZR,eIDR) 
                
         if (MIN(NL,NR) < OrderAcrossFace(MAX(NL,NR))) then 
            same = .FALSE. 
            if (NL<NR) then 
               NNew(indZL,eIDL) = OrderAcrossFace(NR) 
            else 
               NNew(indZR,eIDR) = OrderAcrossFace(NL) 
            end if 
         end if 
      end if
   end subroutine AdjustPAcrossFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to make the p-representation on certain boundaries conforming
!  -----------------------------------------------------------------------
   function getAdaptationType() RESULT(adaptationType)
      !--------------------------------------
      implicit none
      !--------------------------------------
      integer                        :: adaptationType  ! 0 for Truncation Error and 1 for Reinforcement Learning (VI algorithm)
      character(LINE_LENGTH)         :: paramFile
      character(LINE_LENGTH)         :: in_label
      character(LINE_LENGTH)         :: R_Type
      !--------------------------------------
      write(in_label , '(A)') "#define p-adaptation"
      
      call get_command_argument(1, paramFile) !
      
      call readValueInRegion ( trim ( paramFile )  , "adaptation type"  , R_Type    , in_label , "# end" )

      if ( R_Type /= "" ) then
         select case ( trim (R_Type) )
            case ("te")
               adaptationType = 0
            case ("rl")
               adaptationType = 1
            case default
               adaptationType = 0
               write(STD_OUT,*) "Unknown type. Default type used for p-adaptation: Truncation Error"
         end select
      else
         adaptationType = 0
      end if

   end function
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Subroutines to specify the maximum polynomial order jump across a face
!  The element across the face is of order a
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   pure integer function NumberN_1(a)
      implicit none
      integer, intent(in) :: a
      
      NumberN_1 = a-1
   end function NumberN_1
! 
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
! 
   pure integer function NumberN23(a)
      implicit none
      integer, intent(in) :: a
      
      NumberN23 = (a*2)/3
   end function NumberN23
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
   pure integer function SameNumber(a)
      implicit none
      integer, intent(in) :: a
      
      SameNumber = a
   end function SameNumber
   
   
! 
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
! 
   function RemainingIndex(a) RESULT(b) 
      integer :: a(2) 
      integer :: b 
       
      if (any(a==1)) then 
         if (any(a==2)) b = 3 
         if (any(a==3)) b = 2 
      else 
         b = 1 
      end if 
       
   end function 

   
end module pAdaptationClass