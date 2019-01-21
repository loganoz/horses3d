!
!//////////////////////////////////////////////////////
!
!   @File:    pAdaptationClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Sun Dec 10 12:57:00 2017
!   @Last revision date: Mon Jan 21 01:04:17 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 713d5b12bf01d53752c5c65264e93945cf78090e
!
!//////////////////////////////////////////////////////
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
#ifdef NAVIERSTOKES
   use PhysicsStorage                  , only: NTOTALVARS, CTD_IGNORE_MODE, flowIsNavierStokes
#else
   use PhysicsStorage                  , only: NTOTALVARS, CTD_IGNORE_MODE
#endif
   use FaceClass                       , only: Face
   use ElementClass
   use ElementConnectivityDefinitions  , only: axisMap
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f, MaxTimeStep, ComputeMaxResiduals
   use TruncationErrorClass
   use FTValueDictionaryClass          , only: FTValueDictionary
   use StorageClass
   use FileReadingUtilities            , only: RemovePath, getFileName, getIntArrayFromString, getRealArrayFromString, getCharArrayFromString, GetRealValue, GetIntValue
   use FileReaders                     , only: ReadOrderFile
   use ParamfileRegions                , only: readValueInRegion, getSquashedLine
   use HexMeshClass                    , only: HexMesh
   use ElementConnectivityDefinitions  , only: neighborFaces
   use Utilities                       , only: toLower
   use ReadMeshFile                    , only: NumOfElemsFromMeshFile
   use ExplicitMethods                 , only: TakeRK3Step
   use InterpolationMatrices           , only: Interp3DArrays
   implicit none
   
#include "Includes.h"
   private
   public GetMeshPolynomialOrders, ReadOrderFile
   public pAdaptation_t, ADAPT_UNSTEADY_TIME
   
   
   integer, parameter :: ADAPT_STEADY = 0
   integer, parameter :: ADAPT_UNSTEADY_ITER = 1
   integer, parameter :: ADAPT_UNSTEADY_TIME = 2
   integer, parameter :: NO_ADAPTATION = 3
   integer, parameter :: MAX_STEPS_SMOOTHING = 1000
   
   
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
   
   !---------------------------------
   ! Type for storing the source term
   !---------------------------------
   type :: SourceStorage_t
      integer                    :: Nold(NDIM)
      real(kind=RP), allocatable :: S(:,:,:,:)
   end type SourceStorage_t
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type :: pAdaptation_t
      character(len=LINE_LENGTH)        :: solutionFileName                ! Name of file for plotting adaptation information
      real(kind=RP)                     :: reqTE                           ! Requested truncation error
      logical                           :: RegressionFiles = .FALSE.       ! Write regression files?
      logical                           :: saveGradients                   ! Save gradients in pre-adapt and p-adapted solution files?
      logical                           :: Adapt                           ! Is the adaptator going to be used??
      logical                           :: increasing      = .FALSE.       ! Performing an increasing adaptation procedure?
      logical                           :: Constructed      ! 
      logical                           :: restartFiles    = .FALSE.
      logical                           :: UnSteady
      integer                           :: NxyzMax(3)                      ! Maximum polynomial order in all the directions
      integer                           :: TruncErrorType                  ! Truncation error type (either ISOLATED_TE or NON_ISOLATED_TE)
      integer                           :: adaptation_mode = NO_ADAPTATION ! Adaptation mode 
      integer, allocatable              :: maxNdecrease
      real(kind=RP)                     :: time_interval
      integer                           :: iter_interval
      logical                           :: performPAdaptationT
      real(kind=RP)                     :: nextAdaptationTime = huge(1._RP)
      ! Variables for post-smoothing
      logical                           :: postSmoothing   = .FALSE.       ! Signals if a smoothing operation must be performed post p-adaptation
      logical                           :: Compute_dt  = .FALSE.
      real(kind=RP)                     :: postSmoothRes                   ! Post smoothing residual
      real(kind=RP)                     :: cfl, dcfl
      real(kind=RP)                     :: dt
      type(SourceStorage_t),allocatable :: Source(:)
      !
      character(len=BC_STRING_LENGTH), allocatable :: conformingBoundaries(:) ! Stores the conforming boundaries (the length depends on FTDictionaryClass)
      type(TruncationError_t), allocatable :: TE(:)         ! Truncation error for every element(:)
      
      type(overenriching_t)  , allocatable :: overenriching(:)
      
      contains
         procedure :: construct  => pAdaptation_Construct
         procedure :: destruct   => pAdaptation_Destruct
         procedure :: makeBoundariesPConforming
         procedure :: pAdaptTE   => pAdaptation_pAdaptTE
         procedure :: hasToAdapt => pAdaptation_hasToAdapt
         procedure :: StoreQdot  => pAdaptation_StoreQdot
         procedure :: postSmooth => pAdaptation_postSmooth
   end type pAdaptation_t
!
!  ----------
!  Interfaces
!  ----------
!   
   abstract interface
      pure integer function OrderAcrossFace_f(a)
         integer, intent(in) :: a
      end function
   end interface
!
!  ----------------
!  Module variables
!  ----------------
!
   integer    :: NMIN(NDIM) = 1
   integer    :: NMINest    = 1     ! Minimum polynomil order used for estimation 
   integer    :: NInc_0     = 4
!!    integer               :: dN_Inc = 3 
   integer    :: fN_Inc = 2
   integer    :: NInc
   integer    :: nelem           ! number of elements in mesh
   logical    :: reorganize_Nz = .FALSE.
   
   procedure(OrderAcrossFace_f), pointer :: GetOrderAcrossFace
   
   ! Here we define the input variables that can be changed after p-adaptation
   character(len=18), parameter :: ReplaceableVars(4) = (/'mg sweeps         ', &
                                                          'mg sweeps pre     ', &
                                                          'mg sweeps post    ', &
                                                          'mg sweeps coarsest'/)
   type(FTValueDictionary) :: ReplacedVariables
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
!  Routine for constructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Construct(this,controlVariables,t0)
      implicit none
      !--------------------------------------
      class(pAdaptation_t)   , intent(inout) :: this             !>  P-Adaptator
      type(FTValueDictionary), intent(in)    :: controlVariables !<  Input values
      real(kind=RP)          , intent(in)    :: t0
      !--------------------------------------
      ! For block reading
      character(LINE_LENGTH)         :: paramFile
      character(LINE_LENGTH)         :: in_label
      character(20*BC_STRING_LENGTH) :: confBoundaries
      character(LINE_LENGTH)         :: R_Nmax, R_Nmin, R_TruncErrorType, R_OrderAcrossFaces, replacedValue, R_mode, R_interval
      logical      , allocatable     :: R_increasing, regressionFiles, reorganize_z, R_restart
      real(kind=RP), allocatable     :: TruncError, R_pSmoothing
      ! Extra vars
      integer                        :: i      ! Element counter
      integer                        :: no_of_overen_boxes
      !--------------------------------------
      
      this % Adapt  = pAdaptationIsDefined()
      
      if (this % Adapt) then
         this % Constructed = .TRUE.
      else
         this % Constructed = .FALSE.
         return
      end if
      
      nelem = NumOfElemsFromMeshFile( controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH) )
!
!     **************************************************
!     * p-adaptation is defined - Proceed to construct *     
!     **************************************************
!

!     Read block
!     **********
      write(in_label , '(A)') "#define p-adaptation"
      
      call get_command_argument(1, paramFile) !
      
      call readValueInRegion ( trim ( paramFile )  , "conforming boundaries"  , confBoundaries     , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "increasing"             , R_increasing       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "truncation error"       , TruncError         , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "regression files"       , regressionFiles    , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "adjust nz"              , reorganize_z       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmax"                   , R_Nmax             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmin"                   , R_Nmin             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "truncation error type"  , R_TruncErrorType   , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "order across faces"     , R_OrderAcrossFaces , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "max n decrease"         , this % maxNdecrease, in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "mode"                   , R_mode             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "interval"               , R_interval         , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "restart files"          , R_restart          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "post smoothing residual", R_pSmoothing       , in_label , "# end" )
      
!     Conforming boundaries?
!     ----------------------
      if ( confBoundaries /= "" ) then
         call getCharArrayFromString (confBoundaries,BC_STRING_LENGTH,this % conformingBoundaries)
         do i=1, size(this % conformingBoundaries)
            call toLower(this % conformingBoundaries(i))
         end do
      end if
      
!     Increasing adaptation?
!     ----------------------
      if ( allocated(R_increasing) ) then
         this % increasing = R_increasing
      end if
      
!     Truncation error threshold
!     --------------------------
      if ( allocated(TruncError) ) then
         this % reqTE = TruncError
      else
         ERROR STOP 'A truncation error must be specified for p-adapt'
      end if
      
!     Regression files?
!     -----------------
      if ( allocated(regressionFiles) ) then
         this % RegressionFiles = regressionFiles
      end if
      
!     Adjust Nz?
!     ----------
      if ( allocated(reorganize_z) ) then
         reorganize_Nz = reorganize_z
      end if
      
!     Nmax
!     ----
      if ( R_Nmax /= "" ) then
         this % NxyzMax = getIntArrayFromString(R_Nmax)
      else
         ERROR STOP 'Keyword Nmax is mandatory for p-adaptation'
      end if
      
!     Nmin -> If this is a p-nonconforming 3D case, it should be 2
!     ----
      if ( R_Nmin /= "" ) then
         NMIN = getIntArrayFromString(R_Nmin)
      end if
      
!     Truncation error type
!     ---------------------
      if ( R_TruncErrorType /= "" ) then
         call toLower (R_TruncErrorType)
         select case ( trim(R_TruncErrorType) )
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
      
!     Polynomial order jump
!     ---------------------
      if ( R_OrderAcrossFaces /= "" ) then
         select case ( trim (R_OrderAcrossFaces) )
            case ("n*2/3")
               GetOrderAcrossFace => NumberN23
            case default
               GetOrderAcrossFace => NumberN_1
         end select
      else
         GetOrderAcrossFace => NumberN_1
      end if
      
!     Adaptation mode: Steady(default) or unsteady
!     --------------------------------------------
      
      if ( R_mode == "" ) R_mode = "steady"
      select case ( trim(R_mode) )
         case ("time")
            this % adaptation_mode = ADAPT_UNSTEADY_TIME
            this % time_interval   = GetRealValue(R_interval)
            this % iter_interval   = huge(this % iter_interval)
            this % nextAdaptationTime = t0   ! + this % time_interval
         case ("iteration")
            this % adaptation_mode = ADAPT_UNSTEADY_ITER
            this % time_interval   = huge(this % time_interval)
            this % iter_interval   = GetIntValue(R_interval)
         case ("steady")
            this % adaptation_mode = ADAPT_STEADY
            this % time_interval   = huge(this % time_interval)
            this % iter_interval   = huge(this % iter_interval)
         case default
            WRITE(STD_OUT,*) 'Not recognized adaptation mode. Options are:'
            WRITE(STD_OUT,*) '   * time'
            WRITE(STD_OUT,*) '   * iteration'
            WRITE(STD_OUT,*) '   * steady'
      end select
      
!     Restart files
!     -------------
      if ( allocated(R_restart) ) then
         this % restartFiles = R_restart
      end if
      
!     Post smoothing residual and related variables
!     ---------------------------------------------
      if ( allocated(R_pSmoothing) ) then
         this % postSmoothing = .TRUE.
         this % postSmoothRes = R_pSmoothing
         allocate ( this % Source(nelem) )
      
         if (controlVariables % containsKey("cfl")) then
#if defined(NAVIERSTOKES) || defined(INCNS)
            this % Compute_dt = .TRUE.
            this % cfl        = controlVariables % doublePrecisionValueForKey("cfl")
#if defined(NAVIERSTOKES)
            if (flowIsNavierStokes) then
               if (controlVariables % containsKey("dcfl")) then
                  this % dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
               else
                  ERROR STOP '"cfl" and "dcfl", or "dt" keyword must be specified for the time integrator'
               end if
            end if
#endif
#elif defined(CAHNHILLIARD)
            print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
            errorMessage(STD_OUT)
            stop
#endif
         elseif (controlVariables % containsKey("dt")) then
            this % Compute_dt = .FALSE.
            this % dt         = controlVariables % doublePrecisionValueForKey("dt")
         else
            ERROR stop '"cfl" or "dt" keyword must be specified for the time integrator'
         end if
      end if
      
!     Check replaceable control variables
!     -----------------------------------
      call ReplacedVariables % initWithSize(16)
      
      do i=1, size(ReplaceableVars)
         call readValueInRegion ( trim ( paramFile )  , "padapted " // trim(ReplaceableVars(i)), replacedValue, in_label , "# end" )
         if (replacedValue == "") cycle
         call ReplacedVariables % addValueForKey(trim(replacedValue),TRIM(ReplaceableVars(i)))
      end do
!
!     Some things are read from the control file
!     ******************************************
!      
      this % solutionFileName = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
      this % saveGradients    = controlVariables % logicalValueForKey("save gradients with solution")
      if ( trim( controlVariables % StringValueForKey("simulation type",LINE_LENGTH) ) == 'time-accurate' ) this % UnSteady = .TRUE.
      
      
!     Construct the truncation error
!     ******************************    
      
      allocate (this % TE(nelem))
      do i = 1, nelem
         call this % TE(i) % construct(NMINest,this % NxyzMax)
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
      
   end subroutine pAdaptation_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   pure function pAdaptation_hasToAdapt(this,iter) result(hasToAdapt)
      implicit none
      class(pAdaptation_t), intent(in) :: this
      integer             , intent(in) :: iter
      logical                          :: hasToAdapt
      
      select case (this % adaptation_mode)
         
         case (ADAPT_STEADY)
            hasToAdapt = .FALSE.
            
         case (ADAPT_UNSTEADY_ITER)
            if ( (mod(iter, this % iter_interval) == 0) .or. (iter == 1) ) then
               hasToAdapt = .TRUE.
            else
               hasToAdapt = .FALSE.
            end if
            
         case (ADAPT_UNSTEADY_TIME)
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
            stop "Stopped."

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
!  ----------------------------------------
!  Routine for destructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Destruct(this)
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
      
      safedeallocate  (this % conformingBoundaries)
      safedeallocate  (this % overenriching)
   end subroutine pAdaptation_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Main routine for adapting the polynomial order in all elements based on 
!  the truncation error estimation
!  ------------------------------------------------------------------------
   subroutine pAdaptation_pAdaptTE(pAdapt,sem,itera,t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
      use AnisFASMultigridClass
      use StopwatchClass
      implicit none
      !--------------------------------------
      class(pAdaptation_t)       :: pAdapt            !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
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
      logical                    :: success
      integer, save              :: Stage = 0         !   Stage of p-adaptation for the increasing method
      CHARACTER(LEN=LINE_LENGTH) :: newInput          !   Variable used to change the input in controlVariables after p-adaptation 
      !--------------------------------------
      ! For new adaptation algorithm
      integer                    :: Pxyz(3)           !   Polynomial order the estimation was performed with
      real(kind=RP), allocatable :: TEmap(:,:,:)      !   Map of truncation error
      real(kind=RP)              :: TE_L
      integer                    :: L(3)
      integer                    :: P_1(3)            !   Pxyz-1
      integer                    :: i,j,k             !   Counters
      integer                    :: DOFs, NewDOFs
      logical                    :: notenough(3)
      TYPE(AnisFASMultigrid_t)   :: AnisFASpAdaptSolver
      character(len=LINE_LENGTH) :: AdaptedMeshFile
      logical                    :: last
      !--------------------------------------
#if defined(NAVIERSTOKES)
      
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
      if (pAdapt % restartFiles) then
         write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( pAdapt % solutionFileName ), '_pre-Adapt_Stage_', Stage, '.hsol'
         call sem % mesh % Export(AdaptedMeshFile)         
         call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),pAdapt % saveGradients)
      end if
!
!     -------------------------------------------------------------
!     Estimate the truncation error using the anisotropic multigrid
!     -------------------------------------------------------------
!
      CALL AnisFASpAdaptSolver % construct(controlVariables,sem,estimator=.TRUE.,NMINestim = NMINest)
      CALL AnisFASpAdaptSolver % solve(itera,t,computeTimeDerivative,ComputeTimeDerivativeIsolated,pAdapt % TE, pAdapt % TruncErrorType)
      CALL AnisFASpAdaptSolver % destruct
!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!
      if ( pAdapt % postSmoothing) then
         call pAdapt % StoreQdot(sem, t, ComputeTimeDerivative)
      end if
!
!     -------------------------------------------------------------
!     Find the polynomial order that fulfills the error requirement
!     -------------------------------------------------------------
!
      ! Allocate the TEmap with the maximum number of N compinations and initialize it
      allocate(TEmap(NMIN(1):pAdapt % NxyzMax(1),NMIN(2):pAdapt % NxyzMax(2),NMIN(3):pAdapt % NxyzMax(3)))
      
      ! Loop over all elements
!!!!$OMP PARALLEL do PRIVATE(Pxyz,P_1,TEmap,NewDOFs,i,j,k)  ! TODO: Modify private and uncomment
      do iEl = 1, nelem
         
         Pxyz = sem % mesh % elements(iEl) % Nxyz
         P_1  = Pxyz - 1
         
         do dir=1, NDIM
            if ( (P_1(dir) < NMIN(dir)) .and. (P_1(dir) < NMINest+1) ) P_1(dir) = NMIN(dir)
         end do
         
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
         
         do k = NMIN(3), P_1(3)
            do j = NMIN(2), P_1(2)
               do i = NMIN(1), P_1(1) 
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
!        Extra check:
!           When NMIN > P_1, it is not always necessary to extrapolate.
!           If Tau(P_1) < threshold, then NMIN should fulfill it too
!        -----------------------------------------------------------------------
!
         if ( any(NMIN > P_1) .and. all( P_1 >= NMINest ) ) then
            do dir=1, NDIM
               if (P_1(dir) <= NMIN(dir)) then
                  L(dir) = P_1(dir)
               else
                  L(dir) = NMIN(dir)
               end if
            end do
            
            TE_L= pAdapt % TE(iEl) % Dir(1) % maxTE(L(1)) + &  !xi   contribution
                  pAdapt % TE(iEl) % Dir(2) % maxTE(L(2)) + &  !eta  contribution
                  pAdapt % TE(iEl) % Dir(3) % maxTE(L(3))      !zeta contribution
            
            if (TE_L < pAdapt % reqTE) NNew(:,iEl) = NMIN
         end if
               
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
               do k = NMIN(3), pAdapt % NxyzMax(3)
                  do j = NMIN(2), pAdapt % NxyzMax(2)
                     do i = NMIN(1), pAdapt % NxyzMax(1)
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
         
!~         if (iEl==1) call PrintTEmap(NMIN,TEmap,1,pAdapt % solutionFileName)
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
      
!
!     ----------------------------
!     Overenrich specified regions
!     ----------------------------
!
      call OverEnrichRegions(pAdapt % overenriching,sem % mesh,NNew, pAdapt % NxyzMax)
!
!     ---------------------------------------------------
!     Restrict polynomial order decrease in every element
!     ---------------------------------------------------
!
      if ( allocated(pAdapt % maxNdecrease) ) then
         do iEl=1, nelem
            associate (e => sem % mesh % elements(iEl) )
            do dir=1, NDIM
               
               if ( NNew(dir,iEl) < (e % Nxyz(dir) - pAdapt % maxNdecrease) ) NNew(dir,iEl) = e % Nxyz(dir) - pAdapt % maxNdecrease
               
            end do
            end associate
         end do
      end if
!
!     ------------------------------
!     Restrict polynomial order jump
!     ------------------------------
!
      last = .FALSE.
      do while (.not. last)
         last = .TRUE.
         call pAdapt % makeBoundariesPConforming(sem % mesh,NNew,last)
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
      
      call sem % mesh % pAdapt (NNew, controlVariables)
      
      ! Reconstruct probes
      do i=1, sem % monitors % no_of_probes
         call sem % monitors % probes(i) % Initialization (sem % mesh, i, trim(sem % monitors % solution_file), .FALSE.)
      end do
      
      call Stopwatch % Pause("Preprocessing")
      call Stopwatch % Start("Solver")
      
!
!     ------------------------------------------------------
!     Perform the post p-adaptation smoothing (if requested)
!     ------------------------------------------------------
!
      if ( pAdapt % postSmoothing ) then
         call pAdapt % postSmooth(sem, t, ComputeTimeDerivative)
      end if
!
!     ---------------------------------------------------
!     Write post-adaptation mesh, solution and order file
!     ---------------------------------------------------
!
      if ( pAdapt % UnSteady) then
         write(AdaptedMeshFile,'(A,A,I10.10,A)')  trim( pAdapt % solutionFileName ), '_', itera+1, '.hsol'
      else
         write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( pAdapt % solutionFileName ), '_p-Adapted_Stage_', Stage, '.hsol'
      end if
      
      call sem % mesh % Export(AdaptedMeshFile)
      call sem % mesh % ExportOrders(AdaptedMeshFile)
      
      if (pAdapt % restartFiles) call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),pAdapt % saveGradients)
      
!
!     ------------------------------------------------
!     Rewrite controlVariables according to user input
!     ------------------------------------------------
!
      do i = 1, size(ReplaceableVars)
         if ( ReplacedVariables % containsKey( trim(ReplaceableVars(i)) ) ) then
            newInput = ReplacedVariables % stringValueForKey(trim(ReplaceableVars(i)), requestedLength = LINE_LENGTH)
            call controlVariables % removeObjectForKey( trim(ReplaceableVars(i)) )
            call controlVariables % addValueForKey( TRIM(newInput) , trim(ReplaceableVars(i)) )
         end if
      end do
!
!     ----------------
!     Update residuals
!     ----------------
!
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      
      write(STD_OUT,*) '****    p-Adaptation done, DOFs=', SUM((NNew(1,:)+1)*(NNew(2,:)+1)*(NNew(3,:)+1)), '****'

#endif
   end subroutine pAdaptation_pAdaptTE
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------
!  Subroutine to post-smooth the solution after p-adaptation
!  ---------------------------------------------------------
   subroutine pAdaptation_storeQdot(this, sem, t, ComputeTimeDerivative)
      implicit none
      !-arguments--------------------------------------------------
      class(pAdaptation_t), intent(inout) :: this            !<> Adaptation class
      type(DGSem)         , intent(inout) :: sem
      real(kind=RP)       , intent(in)    :: t
      procedure(ComputeTimeDerivative_f)  :: ComputeTimeDerivative
      !-local-variables--------------------------------------------
      integer :: eID
      integer :: N(NDIM)
      integer :: nEqn
      !------------------------------------------------------------
      
      nEqn = size(sem % mesh % storage % elements(1) % Qdot,1)
      
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      
      do eID = 1, sem % mesh % no_of_elements
         N = sem % mesh % elements(eID) % Nxyz
         
         this % Source(eID) % Nold = N
         allocate (this % Source(eID) % S ( nEqn, N(1), N(2), N(3) ) )
         
         this % Source(eID) % S = - sem % mesh % storage % elements(eID) % Qdot
      end do
      
   end subroutine pAdaptation_storeQdot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------
!  Subroutine to post-smooth the solution after p-adaptation
!  ---------------------------------------------------------
   subroutine pAdaptation_postSmooth(this, sem, t, ComputeTimeDerivative)
      implicit none
      !-arguments--------------------------------------------------
      class(pAdaptation_t), intent(inout) :: this            !<> Adaptation class
      type(DGSem)         , intent(inout) :: sem
      real(kind=RP)       , intent(in)    :: t
      procedure(ComputeTimeDerivative_f)  :: ComputeTimeDerivative
      !-local-variables--------------------------------------------
      integer       :: k
      integer       :: nEqn
      real(kind=RP),allocatable :: maxResidual(:)
      !------------------------------------------------------------
      
      nEqn = size(sem % mesh % storage % elements(1) % Qdot,1)
      
      allocate ( maxResidual(nEqn) )
      
!     Assign source term
!     ------------------
#if defined(NAVIERSTOKES) || defined(INCNS)
      do k=1, sem % mesh % no_of_elements
         
         call Interp3DArrays  (  nEqn, &
                                 this % Source(k) % Nold        , this % Source(k) % S, &
                                 sem % mesh % elements(k) % Nxyz, sem % mesh % storage % elements(k) % S_NS )
         
      end do
#endif
      
!     Smooth solution
!     ---------------
      do k=1, MAX_STEPS_SMOOTHING
         if ( this % Compute_dt ) this % dt = MaxTimeStep( sem, this % cfl, this % dcfl )
         
         call TakeRK3Step( sem % mesh, sem % particles, t, this % dt, ComputeTimeDerivative )
         
         maxResidual = ComputeMaxResiduals(sem % mesh)
         if (maxval(maxResidual) <= this % postSmoothRes )  then
            print*, 'post-smoothed in', k, 'iterations'
            exit
         end if
      end do
      
!     Remove source term
!     ------------------
#if defined(NAVIERSTOKES) || defined(INCNS)
      do k=1, sem % mesh % no_of_elements
         sem % mesh % storage % elements(k) % S_NS = 0._RP
         deallocate ( this % Source(k) % S )
      end do
#endif
   end subroutine pAdaptation_postSmooth   
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
      
      write(STD_OUT,*) '## Forcing p-conforming boundaries ##'
      
!     ************************
!     Loop over the boundaries
!     ************************
      do zoneID = 1, size(mesh % zones)
         
         if ( all ( this % conformingBoundaries /= trim(mesh % zones(zoneID) % Name) ) ) cycle
         
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
               do f_conf = 1, size(neighborFaces,1)
                  
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
!  Subroutines to specify the maximum polynomial order jump across a face
!  The element accross the face is of order a
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
      
      real(kind=RP)              :: x   (P_1-NMINest+1)
      real(kind=RP)              :: y   (P_1-NMINest+1)
      integer                    :: N
      integer                    :: i
      real(kind=RP)              :: C,eta,r             ! Regression variables
      character(len=LINE_LENGTH) :: RegfileName
      integer                    :: fd
      !---------------------------------------
      
      ! Initializations
      error     = 0
      notenough = .FALSE.
      
      ! Check if there are enough points for regression
      if (P_1 < NMINest + 1) then
         notenough = .TRUE.
         return
      end if
      
      ! Perform regression analysis   
      N = P_1 - NMINest + 1
      y = LOG10(Adaptator % TE(iEl) % Dir(Dir) % maxTE (NMINest:P_1))
      x = (/ (real(i,RP), i=NMINest,P_1) /)
      call C_and_eta_estimation(N,x,y,C,eta,r)
      
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
   pure subroutine C_and_eta_estimation(N,x,y,C,eta,r)
      implicit none
      !--------------------------------------
      integer      , intent(in)  :: N        !< Number of points
      real(kind=RP), intent(in)  :: x(N) 
      real(kind=RP), intent(in)  :: y(N)
      real(kind=RP), intent(out) :: C
      real(kind=RP), intent(out) :: eta
      real(kind=RP), intent(out) :: r
      !--------------------------------------
      real(kind=RP)              :: sumx,sumy,sumsqx,sumxy,deno,sumsqy
      !--------------------------------------
      
      sumx = sum(x)
      sumy = sum(y)
      sumsqx = sum(x*x)
      sumsqy = sum(y*y)
      sumxy  = sum(x*y)
      
      deno=n*sumsqx-sumx*sumx
      
      eta = (n*sumxy-sumx*sumy)/deno
      C   = (sumsqx*sumy-sumx*sumxy)/deno
      r   = (  sumxy-sumx*sumy/n) / sqrt((sumsqx - sumx*sumx/n) * (sumsqy - sumy*sumy/n))
      
   end subroutine C_and_eta_estimation
   
end module pAdaptationClass
