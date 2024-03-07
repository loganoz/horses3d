!
!//////////////////////////////////////////////////////
!
!      Class cotaining routines for adapting polynomial orders based on Reinforcement Learning.
!        -> The adaptation procedure is performed with a RL agent trined with a Value Iteration algorithm.
!        -> The current implementation is compatible with OpenMP and MPI.
!
!////////////////////////////////////////////////////////////////////////
!
module pAdaptationClassRL
   use SMConstants
#ifdef NAVIERSTOKES
   use PhysicsStorage                  , only: CTD_IGNORE_MODE, flowIsNavierStokes
#else
   use PhysicsStorage                  , only: CTD_IGNORE_MODE
#endif
   use ElementClass
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f, MaxTimeStep, ComputeMaxResiduals
   use FTValueDictionaryClass          , only: FTValueDictionary
   use StorageClass
   use MPI_Process_Info
   use FileReadingUtilities            , only: getFileName, getIntArrayFromString, getCharArrayFromString, GetRealValue, GetIntValue
   use ParamfileRegions                , only: readValueInRegion
   use Utilities                       , only: toLower
   use ReadMeshFile                    , only: NumOfElemsFromMeshFile
   use ExplicitMethods                 , only: TakeRK3Step
   use InterpolationMatrices           , only: Interp3DArrays
   use MultiTauEstimationClass         , only: MultiTauEstim_t
   use StopwatchClass                  , only: Stopwatch
   use pAdaptationClass      
   use ReinforcementLearning           , only: pAdaptationAgent_t   
   
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
#include "Includes.h"
   private
   public pAdaptationRL_t
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type, extends(pAdaptation_t) :: pAdaptationRL_t
      type(pAdaptationAgent_t) :: agent
      real(kind=RP)            :: threshold
      real(kind=RP)            :: tol = 1e-2_RP
      
      contains
         ! Base class procedures
         procedure :: pAdaptation_Construct
         procedure :: pAdaptation_Destruct
         procedure :: pAdaptation_pAdapt
   end type pAdaptationRL_t
!
!  ----------------
!  Module variables
!  ----------------
!
   integer    :: NMIN(NDIM) = 1

!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for constructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Construct(this,controlVariables,t0)
      implicit none
      !--------------------------------------
      class(pAdaptationRL_t) , intent(inout) :: this             !>  P-Adaptator
      type(FTValueDictionary), intent(in)    :: controlVariables !<  Input values
      real(kind=RP)          , intent(in)    :: t0
      !--------------------------------------
      ! For block reading
      character(LINE_LENGTH)         :: paramFile
      character(LINE_LENGTH)         :: in_label
      character(LINE_LENGTH)         :: agentFile
      character(20*BC_STRING_LENGTH) :: confBoundaries
      character(LINE_LENGTH)         :: R_Nmax, R_Nmin, R_OrderAcrossFaces, replacedValue, R_mode, R_interval, cwd
      logical      , allocatable     :: R_increasing, reorganize_z, R_restart
      real(kind=RP), allocatable     :: R_tolerance, R_threshold
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
      call readValueInRegion ( trim ( paramFile )  , "tolerance"              , R_tolerance        , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "adjust nz"              , reorganize_z       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmax"                   , R_Nmax             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmin"                   , R_Nmin             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "order across faces"     , R_OrderAcrossFaces , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "mode"                   , R_mode             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "interval"               , R_interval         , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "restart files"          , R_restart          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "agent file"             , agentFile          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "threshold"              , R_threshold        , in_label , "# end" )
      
!     Conforming boundaries
!     ----------------------
      if ( confBoundaries /= "" ) then
         call getCharArrayFromString (confBoundaries,BC_STRING_LENGTH,this % conformingBoundaries)
         do i=1, size(this % conformingBoundaries)
            call toLower(this % conformingBoundaries(i))
         end do
      end if
      
!     Truncation error tolerance
!     --------------------------
      if ( allocated(R_tolerance) ) then
         this % reqTE = R_tolerance !Required for static adaptation mode
         this % tol = R_tolerance
      end if

!     Adaptation threshold
!     ---------------------
      if ( allocated(R_threshold) ) then
         this % threshold = R_threshold
      else
         this % threshold = 0.0_RP
      end if
      
!     Adjust Nz
!     ----------
      if ( allocated(reorganize_z) ) then
         reorganize_Nz = reorganize_z
      end if
      
!     Nmax
!     ----
      if ( R_Nmax /= "" ) then
         this % NxyzMax = getIntArrayFromString(R_Nmax)
      else
         error stop 'Keyword Nmax is mandatory for p-adaptation'
      end if
      
!     Nmin -> If this is a p-nonconforming 3D case, it should be 2
!     ----
      if ( R_Nmin /= "" ) then
         NMIN = getIntArrayFromString(R_Nmin)
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
      
      if ( R_mode == "" ) R_mode = "static"
      select case ( trim(R_mode) )
         case ("time")
            this % adaptation_mode = ADAPT_DYNAMIC_TIME
            if ( R_interval /= "" ) then
               this % time_interval = GetRealValue(R_interval)
            else
               error stop 'Keyword interval is mandatory for p-adaptation if adaptation mode is time'
            end if
            this % iter_interval   = huge(this % iter_interval)
            this % nextAdaptationTime = t0   ! + this % time_interval
         case ("iteration")
            this % adaptation_mode = ADAPT_DYNAMIC_ITER
            this % time_interval   = huge(this % time_interval)
            if ( R_interval /= "" ) then
               this % iter_interval   = GetIntValue(R_interval)
            else
               error stop 'Keyword interval is mandatory for p-adaptation if adaptation mode is iteration'
            end if
         case ("static")
            this % adaptation_mode = ADAPT_STATIC
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
      
!
!     Stopwatch events
!     ****************
!  
      call Stopwatch % CreateNewEvent("pAdapt: PolOrder selection")
      call Stopwatch % CreateNewEvent("pAdapt: Adaptation")
      
      
!
!     Some things are read from the control file
!     ******************************************
!      
      this % solutionFileName = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
      this % saveGradients    = controlVariables % logicalValueForKey("save gradients with solution")
      this % saveSensor       = controlVariables % logicalValueForKey("save sensor with solution")
      if ( trim( controlVariables % StringValueForKey("simulation type",LINE_LENGTH) ) == 'time-accurate' ) this % UnSteady = .TRUE.
      
      
!     Adaptation overenriching
!     *************************
      
      call getNoOfOverEnrichingBoxes(no_of_overen_boxes)
         
      if (no_of_overen_boxes > 0) then
         allocate ( this % overenriching(no_of_overen_boxes) )
         
         do i = 1, no_of_overen_boxes
            call this % overenriching(i) % initialize (i)
         end do
      end if

!     Policy definition for the RL agent
!     **********************************
      if ( agentFile /= "" ) then
         call this % agent % construct(trim(agentFile))
         NMIN(1) = max(NMIN(1), this % agent % pmin)
         NMIN(2) = max(NMIN(2), this % agent % pmin)
         NMIN(3) = max(NMIN(3), this % agent % pmin)
         this % NxyzMax(1) = min(this % NxyzMax(1), this % agent % pmax)
         this % NxyzMax(2) = min(this % NxyzMax(2), this % agent % pmax)
         this % NxyzMax(3) = min(this % NxyzMax(3), this % agent % pmax)
      else
         error stop 'Keyword agentFile is mandatory for p-adaptation with RL'
      end if
      
   end subroutine pAdaptation_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for destructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Destruct(this)
      implicit none
      !--------------------------------------
      class(pAdaptationRL_t) :: this
      !--------------------------------------
      
      call this % agent % destruct()
      
      safedeallocate  (this % conformingBoundaries)
      safedeallocate  (this % overenriching)
      
   end subroutine pAdaptation_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Main routine for adapting the polynomial order in all elements based on 
!  the Value Iteration RL agent
!  ------------------------------------------------------------------------
   subroutine pAdaptation_pAdapt(this, sem, itera, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
      use AnisFASMultigridClass
      implicit none
      !-arguments----------------------------
      class(pAdaptationRL_t)     :: this              !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables  !<> Input variables (that can be modified depending on the user input)
      !-local-variables----------------------
      integer                    :: eID                                     !   Element counter
      integer                    :: Dir                                     !   Direction
      integer                    :: NNew(3,sem % mesh % no_of_elements)     !   New polynomial orders of mesh (after adaptation!)
      integer, save              :: Stage = 0                               !   Stage of p-adaptation for the increasing method
      CHARACTER(LEN=LINE_LENGTH) :: newInput                                !   Variable used to change the input in controlVariables after p-adaptation 
      character(len=LINE_LENGTH) :: RegfileName
      integer                    :: i                                       !   Counters
      TYPE(AnisFASMultigrid_t)   :: AnisFASpAdaptSolver
      character(len=LINE_LENGTH) :: AdaptedMeshFile
      logical                    :: last
      !-mpi-variables-------------------------
      integer                    :: ierr
      integer                    :: local_DOFs, global_DOFs
      integer                    :: adaptedElements, allAdaptedElements
      real(kind=RP)              :: adaptationPercentage
      !--------------------------------------
#if defined(NAVIERSTOKES)
      
      Stage = Stage + 1    
!
!     -------------------------------------------------------------
!     Find the polynomial order that fulfills the error requirement
!     -------------------------------------------------------------
!
      
      call Stopwatch % Start("pAdapt: PolOrder selection")
      adaptedElements = 0
!$omp parallel do schedule(runtime)
      do eID = 1, sem % mesh % no_of_elements
         call pAdaptation_pAdaptRL_SelectElemPolorders (this, sem % mesh % elements(eID) , NNew(:,eID))
         if ( .not. all( sem % mesh % elements(eID)  % Nxyz == NNew(:,eID)) ) then
!$omp critical
            adaptedElements = adaptedElements + 1
!$omp end critical
         end if
      end do
!$omp end parallel do
      call Stopwatch % Pause("pAdapt: PolOrder selection")

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         call mpi_allreduce ( adaptedElements, allAdaptedElements , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
         adaptationPercentage = 100.0_RP * allAdaptedElements / sem % mesh % no_of_allElements
#endif
      else
         adaptationPercentage = 100.0_RP * adaptedElements / sem % mesh % no_of_elements
      end if  

      ! Only adapt once
      this % Adapt = .FALSE.
      
!
!     -----------------------------------------------------------------------------
!     Adapt only if the percentage of elements to be adapted is above the threshold
!     -----------------------------------------------------------------------------
!
      if (adaptationPercentage > this % threshold) then

         if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,*)
               write(STD_OUT,*)
               write(STD_OUT,'(A)') '****     Performing p-Adaptation with a VI Reinforcement Learning agent    ****'
               write(STD_OUT,*)
            end if
#endif
         else
            write(STD_OUT,*)
            write(STD_OUT,*)
            write(STD_OUT,'(A)') '****     Performing p-Adaptation with a VI Reinforcement Learning agent    ****'
            write(STD_OUT,*)
         end if

         if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,*) '****    Percentage of elements to be adapted: ', adaptationPercentage, '% ****'
               write(STD_OUT,*)
            end if
#endif
         else
            write(STD_OUT,*) '****    Percentage of elements to be adapted: ', adaptationPercentage, '% ****'
            write(STD_OUT,*)
         end if  


!
!     --------------------------------------
!     Write pre-adaptation mesh and solution
!     --------------------------------------
!
      if (this % restartFiles) then
         write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( this % solutionFileName ), '_pre-Adapt_Stage_', Stage, '.hsol'
         call sem % mesh % Export(AdaptedMeshFile)         
         call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),this % saveGradients,this % saveSensor)
      end if
!
!     ----------------------------
!     Overenrich specified regions
!     ----------------------------
!
      call OverEnrichRegions(this % overenriching, sem % mesh, NNew, this % NxyzMax)

!
!     ---------------------------------------------------------------
!     Restrict polynomial order jump and make boundaries p-conforming
!     ---------------------------------------------------------------
!
      last = .FALSE.
      do while (.not. last)
         last = .TRUE.
         call this % makeBoundariesPConforming(sem % mesh, NNew, last)
         ! call ReorganizePolOrders(sem % mesh % faces, NNew, last)  #Not implemented for MPI, but not required if pmax<=6 and pmin>=2
      end do

!
!     ----------------------------------
!     Adapt sem to new polynomial orders
!     ----------------------------------
!
      call Stopwatch % Start("pAdapt: Adaptation")
      call sem % mesh % pAdapt_MPI (NNew, controlVariables)
      call Stopwatch % Pause("pAdapt: Adaptation")
      
      ! Reconstruct probes
      do i=1, sem % monitors % no_of_probes
         call sem % monitors % probes(i) % Initialization (sem % mesh, i, trim(sem % monitors % solution_file), .FALSE.)
      end do
          
!
!     ---------------------------------------------------
!     Write post-adaptation mesh, solution and order file
!     ---------------------------------------------------
!
      if ( this % UnSteady) then
         write(AdaptedMeshFile,'(A,A,I10.10,A)')  trim( this % solutionFileName ), '_', itera+1, '.hsol'
      else
         write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( this % solutionFileName ), '_p-Adapted_Stage_', Stage, '.hsol'
      end if
      
      call sem % mesh % Export(AdaptedMeshFile)
      call sem % mesh % ExportOrders(AdaptedMeshFile)
      
      if (this % restartFiles) call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),this % saveGradients,this % saveSensor)
      
!
!     ----------------
!     Update residuals
!     ----------------
!
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)

!
!     --------------------------------------------------------------------------
!     Perform a reduction to know how many DOFs are in each process
!     --------------------------------------------------------------------------
      local_DOFs = SUM((NNew(1,:)+1)*(NNew(2,:)+1)*(NNew(3,:)+1))

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
      call mpi_reduce ( local_DOFs, global_DOFs , 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
      if ( MPI_Process % isRoot ) then
         write(STD_OUT,*) '****    p-Adaptation done, DOFs=', global_DOFs, '****'
         write(STD_OUT,*)
         write(STD_OUT,*)
      end if
#endif
      else
         write(STD_OUT,*) '****    p-Adaptation done, DOFs=', local_DOFs, '****'
         write(STD_OUT,*)
         write(STD_OUT,*)
      end if

   ! End of adapting
   end if

!End NAVIERSTOKES
#endif 
   end subroutine pAdaptation_pAdapt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------
!  pAdaptation_pAdaptRL_SelectElemPolorders:
!  Select the polynomial orders for one element based on the VI Reinforcement Learning agent
!  -----------------------------------------------------------------------------------------
   subroutine pAdaptation_pAdaptRL_SelectElemPolorders (this, e, NNew)
      implicit none
      !-arguments----------------------------------
      class(pAdaptationRL_t), intent(in)    :: this        !<> Adaptation class
      type(Element)         , intent(in)    :: e
      integer               , intent(out)   :: NNew(NDIM)
      !-local-variables----------------------------
      integer                    :: Pxyz(NDIM)     ! Initial polynomial order
      integer                    :: i, j, k     ! Coordinate counters
      integer                    :: dir         ! Coordinate direction
      integer(kind=1)            :: action      ! Action
      integer                    :: indices_dir1(e % Nxyz(1)+1), indices_dir2(e % Nxyz(2)+1), indices_dir3(e % Nxyz(3)+1)
      real(kind=RP)              :: Q_dir1(NDIM, 0:e % Nxyz(1)), Q_dir2(NDIM, 0:e % Nxyz(2)), Q_dir3(NDIM, 0:e % Nxyz(3))
      real(kind=RP)              :: minQ, maxQ
      !--------------------------------------------
      
!     Initialization of P
!     ---------------------------
      Pxyz = e % Nxyz 
      
!     Initialization of NNew
!     -------------------------------
      NNew = -1 ! Initialized to negative value

!     --------------------------------
!     Select the polynomial order in Direction 1
!     --------------------------------
      action = -1 ! Initialized to minimum value
      do i = 0, Pxyz(3)
         do j = 0, Pxyz(2)
            !Compute the state variables for each Gauss-node in axis 1
            do k = 0, Pxyz(1)
               do dir = 1, 3
                  Q_dir1(dir, k) = e % storage % Q(dir+1, k, j, i) !RHOU, RHOV, RHOW
               enddo
            enddo
            !Select the most restrictive action among RHOU, RHOV, RHOW
            do dir = 1, 3
               ! Find minimum and maximum values in axis 1
               minQ = minval(Q_dir1(dir, :))
               maxQ = maxval(Q_dir1(dir, :))
               if (maxQ - minQ < this % tol) then
                  indices_dir1(:) = this % agent % smax + 1
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 1
                  indices_dir1(:) = nint(2 * this % agent % smax * (Q_dir1(dir, :) - minQ) / (maxQ - minQ) - this % agent % smax) + this % agent % smax + 1
               end if
               ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
               action = max(action, this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getData(indices_dir1))
            enddo
         enddo
      enddo
      ! Update NNew and make sure it is within the boundaries
      NNew(1) = min(max(Pxyz(1) + action, NMIN(1)), this % NxyzMax(1))

!     --------------------------------
!     Select the polynomial order in Direction 2
!     --------------------------------
      action = -1 ! Initialized to minimum value
      do i = 0, Pxyz(3)
         do j = 0, Pxyz(1)
            !Compute the state variables for each Gauss-node in axis 2
            do k = 0, Pxyz(2)
               do dir = 1, 3
                  Q_dir2(dir, k) = e % storage % Q(dir+1, j, k, i) !RHOU, RHOV, RHOW
               enddo
            enddo
            !Select the most restrictive action among RHOU, RHOV, RHOW
            do dir = 1, 3
               ! Find minimum and maximum values in axis 2
               minQ = minval(Q_dir2(dir, :))
               maxQ = maxval(Q_dir2(dir, :))
               if (maxQ - minQ < this % tol) then
                  indices_dir2(:) = this % agent % smax + 1
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 2
                  indices_dir2(:) = nint(2 * this % agent % smax * (Q_dir2(dir, :) - minQ) / (maxQ - minQ) - this % agent % smax) + this % agent % smax + 1
               end if
               ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
               action = max(action, this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getData(indices_dir2))
            enddo
         enddo
      enddo
      ! Update NNew and make sure it is within the boundaries
      NNew(2) = min(max(Pxyz(2) + action, NMIN(2)), this % NxyzMax(2))

!     --------------------------------
!     Select the polynomial order in Direction 3
!     --------------------------------
      action = -1 ! Initialized to minimum value
      do i = 0, Pxyz(2)
         do j = 0, Pxyz(1)
            !Compute the state variables for each Gauss-node in axis 3
            do k = 0, Pxyz(3)
               do dir = 1, 3
                  Q_dir3(dir, k) = e % storage % Q(dir+1, j, i, k) !RHOU, RHOV, RHOW
               enddo
            enddo
            !Select the most restrictive action among RHOU, RHOV, RHOW
            do dir = 1, 3
               ! Find minimum and maximum values in axis 3
               minQ = minval(Q_dir3(dir, :))
               maxQ = maxval(Q_dir3(dir, :))
               if (maxQ - minQ < this % tol) then
                  indices_dir3(:) = this % agent % smax + 1
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 3
                  indices_dir3(:) = nint(2 * this % agent % smax * (Q_dir3(dir, :) - minQ) / (maxQ - minQ) - this % agent % smax) + this % agent % smax + 1
               end if
               ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
               action = max(action, this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getData(indices_dir3))
            enddo
         enddo
      enddo
      ! Update NNew and make sure it is within the boundaries
      NNew(3) = min(max(Pxyz(3) + action, NMIN(3)), this % NxyzMax(3))   
      
   end subroutine pAdaptation_pAdaptRL_SelectElemPolorders  
   
end module pAdaptationClassRL