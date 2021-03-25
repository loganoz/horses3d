!
!//////////////////////////////////////////////////////
!
!   @File:    FASMultigridClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Sun Apr 27 12:57:00 2017
!   @Last revision date: Tue Feb  9 22:12:58 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: cd254b7fc25795caf3a285825faa6b46b79c6e79
!
!//////////////////////////////////////////////////////
!
!
!      FAS Multigrid Class
!        Provides the routines for solving a time step with nonlinear multigrid procedures.
!        Available smoothers are:
!           -> RK3:         Explicit smoother based on the Williamson's low-storage 3rd order Runge-Kutta. Only valid for STEADY_STATE cases since the time-stepping is advanced in every level.
!           -> BlockJacobi: Implicit smoother based on the matrix-free version of the classical Block-Jacobi method. STEADY_STATE or TIME_ACCURATE.
!           -> GMRES:       Implicit smoother based on the matrix-free GMRES method. STEADY_STATE or TIME_ACCURATE.
!
!        **NOTE: The implicit smoothers are currently using BDF time-integration. This can be changed...           
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "Includes.h"
module FASMultigridClass
   use SMConstants
   use ExplicitMethods
   use DGSEMClass
   use PhysicsStorage
   use Physics
   use InterpolationMatrices
   use MultigridTypes
   use TimeIntegratorDefinitions
   use LinearSolverClass
   use BDFTimeIntegrator
   use FileReadingUtilities      , only: getFileName
   use MPI_Process_Info          , only: MPI_Process
   use FileReadingUtilities      , only: getIntArrayFromString 
#if defined(NAVIERSTOKES)
   use ManufacturedSolutions
#endif
   
   implicit none
   
   private
   public FASMultigrid_t
   
!
!  Multigrid class
!  ---------------
   type :: FASMultigrid_t
      type(DGSem)              , pointer      :: p_sem                 ! Pointer to DGSem class variable of current system
      type(FASMultigrid_t)     , pointer      :: Child                 ! Next coarser multigrid solver
      type(FASMultigrid_t)     , pointer      :: Parent                ! Next finer multigrid solver
      type(MGSolStorage_t)     , allocatable  :: MGStorage(:)          ! Storage
      class(GenericLinSolver_t), allocatable  :: linsolver             ! Linear solver for implicit smoothing
      integer                                 :: MGlevel               ! Current Multigrid level
      logical                                 :: computeA              !< Compute A in this level?
      real(kind=RP),             allocatable  :: lts_dt(:)             ! dt array for LTS
      contains
         procedure :: construct
         procedure :: solve
         procedure :: Smooth
         procedure :: destruct
         procedure :: SetPreviousSolution => FAS_SetPreviousSolution    ! For implicit smoothing, it's necessary to store the previous solution(s) in all levels
         procedure :: TakePseudoStep ! solve for Dual Time Stepping 

   end type FASMultigrid_t
!
!  ----------------
!  Module variables
!  ----------------
!
   procedure(SmoothIt_t), pointer :: SmoothIt
   
   !! Parameters
   integer, parameter :: MAX_SWEEPS_DEFAULT = 10000
   
   ! Other variables
   integer        :: MaxN           ! Maximum polynomial order of the mesh
   integer        :: NMIN           ! Minimum polynomial order allowed
   integer        :: MGlevels       ! Total number of multigrid levels
   integer        :: deltaN         ! 
   integer        :: nelem          ! Number of elements
   integer        :: num_of_allElems
   integer        :: Smoother       ! Current smoother being used
   integer        :: MaxSweeps      ! Maximum number of sweeps in a smoothing process
   logical        :: MGOutput       ! Display output?
   logical        :: FMG = .FALSE.  ! Use Full Multigrid algorithm?
   logical        :: PostFCycle,PostSmooth ! Post smoothing options
   logical        :: SmoothFine     !      
   logical        :: ManSol         ! Does this case have manufactured solutions?
   logical        :: Compute_dt
   logical        :: SaveFMGFile =.FALSE.
   character(len=LINE_LENGTH)              :: FMGSolutionFile
   logical                                 :: saveGradients
   real(kind=RP)  :: SmoothFineFrac ! Fraction that must be smoothed in fine before going to coarser level
   real(kind=RP), target  :: cfl            ! Advective cfl number
   real(kind=RP), target  :: dcfl           ! Diffusive cfl number
   real(kind=RP), target  :: own_dt             ! dt
   integer, allocatable :: MGSweeps(:) ! Number of pre- and post- smoothings operations on each level
   integer        :: Preconditioner       ! Current smoother being used
   integer        :: CurrentMGCycle 
!-----CFL-ramping-variables-----------------------------------------------------------
   logical        :: CFLboost  = .false.
   logical        :: DCFLboost = .false.
   real(kind=RP)  :: cfl_ini            ! Inital advective cfl number (for CFL boost)
   real(kind=RP)  :: dcfl_ini           ! Initial diffusive cfl number (for CFL boost)
   integer        :: erk_order = 5      ! Steady-state OptERK type 
!-----DTS-variables-------------------------------------------------------------------
   logical        :: DualTimeStepping = .false.
   real(kind=RP), target  :: p_cfl            ! Pseudo advective cfl number
   real(kind=RP), target  :: p_dcfl           ! Pseudo diffusive cfl number
   real(kind=RP), target  :: p_dt             ! Pseudo dt
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,controlVariables,sem)
      use FTValueDictionaryClass
      use StopwatchClass
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t)  , intent(inout), target :: this               !<> Anisotropic FAS multigrid solver to be constructed
      type(FTValueDictionary), intent(in)            :: controlVariables   !<  Input variables
      type(DGSem)            , intent(in)   , target :: sem                !<  Fine sem class
      !-----------------------------------------------------------
      character(len=LINE_LENGTH)                     :: PostSmoothOptions
      integer                                        :: zoneID                ! Zone counter
      logical                                        :: conformingBoundaries  ! Is the mesh conforming on all boundaries?
      character(len=LINE_LENGTH)                     :: tmpc
      integer                                        :: i
      !-----------------------------------------------------------
      
      call Stopwatch % Pause("Solver")
      call Stopwatch % Start("Preprocessing")
      
!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
!
!     Read cfl and dcfl numbers
!     -------------------------
      if (controlVariables % containsKey("cfl")) then
#if defined(NAVIERSTOKES)      
         Compute_dt = .TRUE.
         cfl = controlVariables % doublePrecisionValueForKey("cfl")
         ! if there is no dcfl, use the same number as for cfl 
         if (flowIsNavierStokes) then
            if (controlVariables % containsKey("dcfl")) then
               dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
            else
               dcfl = cfl
            end if
         end if
         ! CFL boosting
         if (controlVariables % containsKey("cfl boost")) then
            CFLboost = controlVariables % logicalValueForKey("cfl boost")
         end if
         ! DCFL boosting
         if (controlVariables % containsKey("dcfl boost")) then
            DCFLboost = controlVariables % logicalValueForKey("dcfl boost")
         end if
#elif defined(CAHNHILLIARD)
         print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
         errorMessage(STD_OUT)
         stop
#endif
      elseif (controlVariables % containsKey("dt")) then
         Compute_dt = .FALSE.
         own_dt = controlVariables % doublePrecisionValueForKey("dt")
      else
         ERROR STOP '"cfl" (and "dcfl" if Navier-Stokes) or "dt" keywords must be specified for the FAS integrator'
      end if

!
!     Pseudo time stepping variables
!     -------------------------
      if ( trim(controlVariables % StringValueForKey("simulation type",LINE_LENGTH)) == "time-accurate" ) then
         ! is the simulation type is time-accurate and we use FAS
         DualTimeStepping = .true.

         if (controlVariables % containsKey("pseudo cfl")) then
#if defined(NAVIERSTOKES)      
            p_cfl = controlVariables % doublePrecisionValueForKey("pseudo cfl")
            if (flowIsNavierStokes) then
               if (controlVariables % containsKey("pseudo dcfl")) then
                  p_dcfl       = controlVariables % doublePrecisionValueForKey("pseudo dcfl")
               else
                  p_dcfl = p_cfl
               end if
            end if
#elif defined(CAHNHILLIARD)
            print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
            errorMessage(STD_OUT)
            stop
#endif
         elseif (controlVariables % containsKey("pseudo dt")) then
            p_dt = controlVariables % doublePrecisionValueForKey("pseudo dt")
         else
            ERROR STOP '"pseudo cfl" or "pseudo dt" keywords must be specified for the time-accurate FAS integrator'
         end if
      end if ! time-accurate

      
!
!     Read multigrid variables
!     -------------------------
      if (.NOT. controlVariables % containsKey("multigrid levels")) then
         print*, 'Fatal error: "multigrid levels" keyword is needed by the FASMultigrid solver'
         STOP
      end if
      
      MGlevels  = controlVariables % IntegerValueForKey("multigrid levels")
      
      if (controlVariables % containsKey("delta n")) then
         deltaN = controlVariables % IntegerValueForKey("delta n")
      else
         deltaN = 1
      end if

      allocate(MGSweeps(MGlevels)) 
      if ( controlVariables % containsKey("mg sweeps") ) then
         tmpc = controlVariables % StringValueForKey("mg sweeps",LINE_LENGTH)
         MGSweeps = getIntArrayFromString(tmpc)
      else
         MGSweeps= 1
         if (MGlevels .gt. 1) then
            do i=2,MGlevels
               MGSweeps(i) = 2*MGSweeps(i-1)
            end do
         end if
      end if
!
!     Select the smoother
!     -------------------
      select case (controlVariables % StringValueForKey("mg smoother",LINE_LENGTH))
         case('Euler')
            Smoother = Euler_SMOOTHER
         case('RK3')
            Smoother = RK3_SMOOTHER
         case('RK5')
            Smoother = RK5_SMOOTHER
         case('RKOpt')
            Smoother = RKOpt_SMOOTHER
            ! Select order of ERK scheme
            if ( controlVariables % containsKey("rk order") ) then
               erk_order = controlVariables % IntegerValueForKey("rk order")
               if (erk_order .gt. 7) then
                  erk_order = 7
                  print *, "FASMultigrid :: ERK Order too high, switching to 7."
               else if (erk_order .lt. 2) then
                  erk_order = 2
                  print *, "FASMultigrid :: ERK Order too low, switching to 2."
               end if
            end if
         case('BlockJacobi')
            Smoother = BJ_SMOOTHER
            call BDF_SetOrder( controlVariables % integerValueForKey("bdf order") )
         case('GMRES')
            Smoother = JFGMRES_SMOOTHER
            call BDF_SetOrder( controlVariables % integerValueForKey("bdf order") )
         case('SIRK')
            !! SmoothIt => TakeSIRKStep
            error stop ':: SIRK smoother not implemented yet'
         case default 
            if (MPI_Process % isRoot) write(STD_OUT,*) '"mg smoother" not recognized. Defaulting to RK3.'
            Smoother = RK3_SMOOTHER
      end select

!
!     Select the preconditioner for smoothing
!     ---------------------------------------
      if (controlVariables % containsKey("mg preconditioner")) then
         select case (controlVariables % StringValueForKey("mg preconditioner",LINE_LENGTH))
         case('LTS')
            Preconditioner = PRECONDIIONER_LTS
         case default
            Preconditioner = PRECONDIIONER_NONE
         end select
      else
         Preconditioner = PRECONDIIONER_NONE
      end if
      
!     Check that the BDF order is consistent
!        (only valid for implicit smoothers)
!     --------------------------------------
      select case (Smoother)
         case(BJ_SMOOTHER)
            if (bdf_order > 1) then
               if (.not. controlVariables % containsKey("dt") ) then
                  ERROR stop ':: "bdf order">1 is only valid with fixed time-step sizes'
               end if
            end if
      end select
      
!
!     Control parameters for mg cycle
!     -------------------------------
      PostSmoothOptions = controlVariables % StringValueForKey("postsmooth option",LINE_LENGTH)
      if (trim(PostSmoothOptions) == 'f-cycle') then
         PostFCycle = .true.
      elseif (trim(PostSmoothOptions) == 'smooth') then
         PostSmooth = .true.
      end if
      
      if (controlVariables % containsKey("smooth fine")) then
         SmoothFine = .TRUE.
         SmoothFineFrac = controlVariables % doublePrecisionValueForKey("smooth fine")
      else
         SmoothFine = .FALSE.
      end if
      
      if (controlVariables % containsKey("max mg sweeps")) then
         MaxSweeps = controlVariables % IntegerValueForKey("max mg sweeps")
      else
         MaxSweeps = MAX_SWEEPS_DEFAULT
      end if
      
      if (controlVariables % logicalValueForKey("fasfmg save solutions") ) then
         SaveFMGFile = .TRUE.
         saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
         FMGSolutionFile = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
      end if
!
!     ------------------------------------------
!     Get the minimum multigrid polynomial order
!     ------------------------------------------
!
      NMIN = 1     
!
!     -----------------------
!     Update module variables
!     -----------------------
!
      MGOutput       = controlVariables % logicalValueForKey("multigrid output")
      plotInterval   = controlVariables % integerValueForKey("output interval")
      ManSol         = sem % ManufacturedSol
      MaxN           = MAX(MAXVAL(sem % mesh % Nx),MAXVAL(sem % mesh % Ny),MAXVAL(sem % mesh % Nz))
      MGlevels       = MIN (MGlevels,MaxN - NMIN + 1)
      
      if (MPI_Process % isRoot) then
         write(STD_OUT,*) 'Constructing FAS Multigrid'
         write(STD_OUT,*) 'Number of levels:', MGlevels
      end if
      
      this % p_sem => sem
      
      nelem = SIZE(sem % mesh % elements)
      num_of_allElems = sem % mesh % no_of_allElements
!
!     --------------------------
!     Create linked solvers list
!     --------------------------
!
      call RecursiveConstructor(this, sem % mesh % Nx, sem % mesh % Ny, sem % mesh % Nz, MGlevels, controlVariables)
      
      call Stopwatch % Pause("Preprocessing")
      call Stopwatch % Start("Solver")
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine RecursiveConstructor(Solver, N1x, N1y, N1z, lvl, controlVariables)
#if defined(NAVIERSTOKES)
      use ManufacturedSolutions
#endif
      use FTValueDictionaryClass
      implicit none
      type(FASMultigrid_t), target  :: Solver
      integer, dimension(:)         :: N1x,N1y,N1z      !<  Order of approximation for every element in current solver
      integer                       :: lvl              !<  Current multigrid level
      type(FTValueDictionary)       :: controlVariables !< Control variables (for the construction of coarse sems
      !----------------------------------------------
      integer, dimension(nelem)           :: N2x,N2y,N2z            !   Order of approximation for every element in child solver
      integer, dimension(num_of_allElems) :: N2xAll,N2yAll,N2zAll   !   Order of approximation for every element in child solver
      integer                             :: i,j,k, iEl             !   Counter
      logical                             :: success                ! Did the creation of sem succeed?
      type(FASMultigrid_t), pointer       :: Child_p                ! Pointer to Child
      integer                             :: Q1,Q2,Q3,Q4            ! Sizes of vector Q (conserved solution) used for allocation. In this version the argument MOLD of ALLOCATE is not used since several versions of gfortran don't support it yet...
      !----------------------------------------------
      !
      integer :: Nxyz(3), fd, l
      
      Solver % MGlevel = lvl
!
!     --------------------------
!     Allocate Multigrid storage
!     --------------------------
!
      ALLOCATE (Solver % MGStorage(nelem))
!$omp parallel do private(Q1,Q2,Q3,Q4) schedule(runtime)
      DO k = 1, nelem
         Q1 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,1)
         Q2 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,2) - 1
         Q3 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,3) - 1
         Q4 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,4) - 1
         ALLOCATE(Solver % MGStorage(k) % Q    (Q1,0:Q2,0:Q3,0:Q4))
         ALLOCATE(Solver % MGStorage(k) % E    (Q1,0:Q2,0:Q3,0:Q4))
         ALLOCATE(Solver % MGStorage(k) % S    (Q1,0:Q2,0:Q3,0:Q4))
         ALLOCATE(Solver % MGStorage(k) % Scase(Q1,0:Q2,0:Q3,0:Q4))

         if (DualTimeStepping) then
            ALLOCATE(Solver % MGStorage(k) % R    (Q1,0:Q2,0:Q3,0:Q4))
            ALLOCATE(Solver % MGStorage(k) % Q0   (Q1,0:Q2,0:Q3,0:Q4))
         end if
         
         if (Smoother >= IMPLICIT_SMOOTHER_IDX) then
            ALLOCATE(Solver % MGStorage(k) % dQ   (Q1,0:Q2,0:Q3,0:Q4))
         end if
         Solver % MGStorage(k) % Scase = 0._RP
      end DO   
!$omp end parallel do

      if (Smoother >= IMPLICIT_SMOOTHER_IDX) then
         if ( .not. allocated(Solver % p_sem % mesh % storage % QdotNS  ) ) then 
            allocate(Solver % p_sem % mesh % storage % QdotNS(Solver % p_sem % mesh % storage % NDOF * NCONS))
         end if 
         if ( .not. allocated(Solver % p_sem % mesh % storage % QNS  ) ) then 
            allocate(Solver % p_sem % mesh % storage % QNS(Solver % p_sem % mesh % storage % NDOF * NCONS))
         end if 
         if ( .not. allocated(Solver % p_sem % mesh % storage % PrevQNS  ) ) then 
            allocate(Solver % p_sem % mesh % storage % PrevQNS(Solver % p_sem % mesh % storage % NDOF * NCONS,1))
         end if 
      end if

      ! allocate array for LTS
      if (Preconditioner .eq. PRECONDIIONER_LTS) allocate( Solver % lts_dt(nelem))
      ! allocate( Solver % lts_dt(nelem))
         
!
!     --------------------------------------------------------------
!     Fill MGStorage(iEl) % Scase if required (manufactured solutions)
!        (only for lower meshes)
!     --------------------------------------------------------------
!
#if defined(NAVIERSTOKES)
      if (ManSol) then
         DO iEl = 1, nelem
            
            DO k=0, Solver % p_sem % mesh % Nz(iEl)
               DO j=0, Solver % p_sem % mesh % Ny(iEl)
                  DO i=0, Solver % p_sem % mesh % Nx(iEl)
                     if (flowIsNavierStokes) then
                        call ManufacturedSolutionSourceNS(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                     else
                        call ManufacturedSolutionSourceEuler(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                     end if
                  end DO
               end DO
            end DO
         end DO
      end if
#endif
!
!     -------------------------------------------
!     Create linear solver for implicit smoothing
!                                 (if needed)
!     -------------------------------------------
!
      if ( Smoother >= IMPLICIT_SMOOTHER_IDX ) then
!
!        Solver initialization
!        ---------------------
         select case ( Smoother )
            case (BJ_SMOOTHER)
               allocate ( MatFreeSmooth_t :: Solver % linsolver )
            case (JFGMRES_SMOOTHER)
               allocate ( MatFreeGMRES_t  :: Solver % linsolver )
         end select
!
!        Solver construction
!        ---------------------
         call Solver % linsolver % construct(Solver % p_sem % NDOF*NCONS, Solver % p_sem % totalNDOF*NCONS, NCONS, controlVariables, Solver % p_sem, BDF_MatrixShift)
         Solver % computeA = .TRUE.
      end if
      
      if (lvl > 1) then
         ALLOCATE  (Solver % Child)
         Child_p => Solver % Child
         Solver % Child % Parent => Solver
         
!
!        ---------------------------------------------
!        Create restriction and prolongation operators
!        ---------------------------------------------
!
         DO k=1, nelem
            call CreateInterpolationOperators(N1x(k), N2x(k), MaxN, MGlevels, lvl-1, DeltaN, Solver % p_sem % nodes)
            call CreateInterpolationOperators(N1y(k), N2y(k), MaxN, MGlevels, lvl-1, DeltaN, Solver % p_sem % nodes)
            call CreateInterpolationOperators(N1z(k), N2z(k), MaxN, MGlevels, lvl-1, DeltaN, Solver % p_sem % nodes)
         end DO
         
         ! Create DGSEM class for child
         ALLOCATE (Child_p % p_sem)
         
         !<old
         ! Get the polynomial orders for constructing the child sem
         N2xAll = 0
         N2yAll = 0
         N2zAll = 0
         do k=1, nelem
            N2xAll( Solver % p_sem % mesh % elements(k) % globID ) = N2x(k)
            N2yAll( Solver % p_sem % mesh % elements(k) % globID ) = N2y(k)
            N2zAll( Solver % p_sem % mesh % elements(k) % globID ) = N2z(k)
         end do
         call Child_p % p_sem % construct (controlVariables = controlVariables,                                          &
                                           Nx_ = N2xAll,    Ny_ = N2yAll,    Nz_ = N2zAll,                               &
                                           success = success,                                                            &
                                           ChildSem = .TRUE. )
         if (.NOT. success) ERROR STOP "Multigrid: Problem creating coarse solver."
         !old>
         
!~!<New
!~         N2(:,1) = N2x
!~         N2(:,2) = N2y
!~         N2(:,3) = N2z
         
!~         ! Copy the sem
!~         Child_p % p_sem = Solver % p_sem
         
!~         ! Mark the mesh as a child mesh
!~         Child_p % p_sem % mesh % child = .TRUE. 
         
!~         ! Adapt the mesh to the new polynomial orders
!~         N2trans = transpose(N2)
!~         call Child_p % p_sem % mesh % pAdapt (N2trans, controlVariables)
!~         call Child_p % p_sem % mesh % storage % PointStorage
!New>
         call RecursiveConstructor(Solver % Child, N2x, N2y, N2z, lvl - 1, controlVariables)
      end if
      
   end subroutine RecursiveConstructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the FAS multigrid solving procedure
!  ---------------------------------------------
   subroutine solve(this, timestep, t, dt, ComputeTimeDerivative, FullMG, tol)
      implicit none
      !-------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      integer                              :: timestep
      real(kind=RP)        , intent(in)    :: t
      real(kind=RP)        , intent(in)    :: dt
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      logical           , OPTIONAL         :: FullMG
      real(kind=RP)     , OPTIONAL         :: tol        !<  Tolerance for full multigrid
      !-------------------------------------------------
      integer :: maxVcycles = 40, i
      real(kind=RP) :: rnorm, xnorm
      
      ThisTimeStep = timestep
      
      if (PRESENT(FullMG) .AND. FullMG) then
         if (.NOT. PRESENT(tol)) ERROR STOP 'FASFMG needs tolerance'
         FMG = .TRUE.
      else
         FMG = .FALSE.
      end if

      
      ! if (Smoother >= IMPLICIT_SMOOTHER_IDX) call this % linsolver % SetOperatorDt(dt)
!
!     -----------------------
!     Perform multigrid cycle
!     -----------------------
!
      if (Smoother >= IMPLICIT_SMOOTHER_IDX) call FAS_SetPreviousSolution(this,MGlevels) ! FINDME - why are we doing this?
      
      if (FMG) then
         call FASFMGCycle(this,t,tol,MGlevels, ComputeTimeDerivative)
      else
         do i = 1, maxVcycles
            CurrentMGCycle = timestep
            call FASVCycle(this,t,dt,MGlevels,MGlevels, ComputetimeDerivative)
            select case(Smoother)
               case ( : (IMPLICIT_SMOOTHER_IDX-1)) ! Only one iteration per pseudo time-step for RK smoothers
                  exit 
               case (IMPLICIT_SMOOTHER_IDX : )  ! Check if the nonlinear problem was solved to a given tolerance
                  rnorm = this % linsolver % Getrnorm()
                  xnorm = this % linsolver % Getxnorm('infinity')
                  print*, 'V-Cycle', i, 'rnorm=', rnorm, 'xnorm', xnorm 
                  if (xnorm<1.e-6_RP) exit
            end select
            if (rnorm > 1e-2_RP) call computeA_AllLevels(this,MGlevels)
         end do ! i
         if (i > 10) call computeA_AllLevels(this,MGlevels) ! Hard-coded: 4
      end if
      
   end subroutine solve  
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the Dual time-stepping procedure:
!  Q_{m+1} = Q_m + d\tau <(> (Q_m - Q_n)/dt + R(Q_m) <)>
!
!  Q_n    => this % MGStorage(id) % Q0
!  Q_m    => this % MGStorage(id) % Qm
!  R(Q_m) => this % MGStorage(id) % R
!  ---------------------------------------------
   subroutine TakePseudoStep(this, timestep, t, dt, ComputeTimeDerivative, FullMG, tol)
      implicit none
      !-------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      integer                              :: timestep
      real(kind=RP)        , intent(in)    :: t
      real(kind=RP)        , intent(in)    :: dt
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      logical           , OPTIONAL         :: FullMG
      real(kind=RP)     , OPTIONAL         :: tol        !<  Tolerance for full multigrid
      !-------------------------------------------------
      integer :: i, id
      integer :: tau_maxit = 10000
      real(kind=RP) :: dQdtau_norm, Qdot_norm
      real(kind=RP) :: tk
      real(kind=RP) :: tau_tol

!
!     -----------------------
!     Solve local steady-state problem
!     -----------------------
!
      write(STD_OUT,*) "Starting pseudo step."

      tk = t

      if (Compute_dt) call MaxTimeStep( self=this % p_sem, cfl=cfl, dcfl=dcfl , MaxDt=own_dt)

      call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)
      Qdot_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
      call ComputePseudoTimeDerivative(this % p_sem % mesh, tk, own_dt)
      dQdtau_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
      tau_tol = min(1e-5_RP,Qdot_norm * 1e-1_RP)

!$omp parallel do schedule(runtime)
      do id = 1, SIZE(this % p_sem % mesh % elements )
         this % p_sem % mesh % elements(id) % storage % prevQ(1) % Q = this % p_sem % mesh % elements(id) % storage % Q
      end do ! id
!$omp end parallel do

      do i = 1, tau_maxit ! if i > 1/dts_factor, then we don't have speed up

         

         ! if (Preconditioner .eq. PRECONDIIONER_LTS) then
         !    if (Compute_dt) call MaxTimeStep( self=this % p_sem, cfl=p_cfl, dcfl=p_dcfl , MaxDt=p_dt, MaxDtVec=this % lts_dt)
         !    call TakeRK3Step( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=tk, deltaT=p_dt, &
         !       ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=.true., global_dt=own_dt )
         ! else
         !    if (Compute_dt) call MaxTimeStep( self=this % p_sem, cfl=p_cfl, dcfl=p_dcfl , MaxDt=p_dt)
         !    call TakeRK3Step( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=tk, deltaT=p_dt, &
         !       ComputeTimeDerivative=ComputeTimeDerivative, dts=.true., global_dt=own_dt )
         ! end if

         
         call this % solve(i, tk, p_dt, ComputeTimeDerivative)

         

!$omp parallel do schedule(runtime)
      do id = 1, SIZE(this % p_sem % mesh % elements )
         if ( any(isnan(this % p_sem % mesh % elements(id) % storage % Q))) error stop "TBC"
      end do ! id
!$omp end parallel do

         ! call Smooth(this,1,tk,p_dt, ComputeTimeDerivative)

         ! tk = tk + p_dt

         call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)
         Qdot_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
         call ComputePseudoTimeDerivative(this % p_sem % mesh, tk, own_dt)
         dQdtau_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))

         
         write(STD_OUT,*) "---------------------------"
         write(STD_OUT,*) "Iteration: ", i
         write(STD_OUT,*) "Time: ", tk
         write(STD_OUT,*) "dQdt res: ", Qdot_norm
         write(STD_OUT,*) "dQdtau res: ", dQdtau_norm
         if (dQdtau_norm .le. tau_tol) exit

      end do

      tk = tk + own_dt

      write(STD_OUT,*) "Finishing pseudo step."
      
   end subroutine TakePseudoStep  
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Signals computeA=.TRUE. in all MG levels
!  -----------------------------------------
   recursive subroutine computeA_AllLevels(this,lvl)
      implicit none
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      integer :: lvl
      
      this % computeA = .TRUE.
      
      if (lvl>1) call computeA_AllLevels(this % Child,lvl-1)
      
   end subroutine computeA_AllLevels
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Recursive subroutine to perform a v-cycle
!  -----------------------------------------
   recursive subroutine FASVCycle(this,t,dt,lvl,MGlevels, ComputeTimeDerivative)
      implicit none
      !----------------------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      real(kind=RP)        , intent(in)    :: t        !<  Simulation time
      real(kind=RP)        , intent(in)    :: dt       !<  Time-step
      integer              , intent(in)    :: lvl      !<  Current multigrid level
      integer              , intent(in)    :: MGlevels !<  Number of finest multigrid level
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      !----------------------------------------------------------------------------
      integer                       :: iEl,iEQ              !Element/equation counter
      type(FASMultigrid_t), pointer :: Child_p              !Pointer to child
      integer                       :: N1(3), N2(3)
      real(kind=RP)                 :: maxResidual(NCONS)
      real(kind=RP)                 :: PrevRes
      real(kind=RP)                 :: NewRes
      integer                       :: sweepcount           ! Number of sweeps done in a point in time
      !----------------------------------------------------------------------------
#if defined(NAVIERSTOKES)      
!
!     -----------------------
!     Pre-smoothing procedure
!     -----------------------
!
!~      this % computeA = .TRUE.
      sweepcount = 0
      DO
         call this % Smooth(MGSweeps(lvl),t,dt, ComputeTimeDerivative)
         sweepcount = sweepcount + MGSweeps(lvl)
         
         if (MGOutput) call PlotResiduals( lvl , sweepcount,this % p_sem % mesh)
         
         if (SmoothFine .AND. lvl > 1) then ! .AND. .not. FMG
            if (FMG .and. MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) < 0.1_RP) exit
            call MGRestrictToChild(this,lvl-1,t, ComputeTimeDerivative)
            call ComputeTimeDerivative(this % Child % p_sem % mesh,this % Child % p_sem % particles, t, CTD_IGNORE_MODE)
            if (DualTimeStepping) call ComputePseudoTimeDerivative(this % Child % p_sem % mesh, t, own_dt)
            
            if (MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) < SmoothFineFrac * MAXVAL(ComputeMaxResiduals(this % Child % p_sem % mesh))) exit
         else
            exit
         end if
         
         if (sweepcount .ge. MaxSweeps) exit
      end DO
      
      PrevRes = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
      
      
      if (lvl > 1) then
         if (.not. SmoothFine) call MGRestrictToChild(this,lvl-1,t, ComputeTimeDerivative)
!
!        --------------------
!        Perform V-Cycle here
!        --------------------
!
         call FASVCycle(this % Child, t, dt, lvl-1,MGlevels, ComputeTimeDerivative)
         
         Child_p => this % Child
!
!        -------------------------------------------
!        Interpolate coarse-grid error to this level
!        -------------------------------------------
!
!$omp parallel 
!$omp do private(N1,N2) schedule(runtime)
         DO iEl = 1, nelem
            N1 = Child_p % p_sem % mesh % elements (iEl) % Nxyz
            N2 = this    % p_sem % mesh % elements (iEl) % Nxyz
            call Interp3DArrays(NCONS, N1, Child_p % MGStorage(iEl) % E, N2, this % MGStorage(iEl) % E )
         end DO
!$omp end do
!
!        -----------------------------------------------
!        Correct solution with coarse-grid approximation
!        -----------------------------------------------
!
!$omp barrier
!$omp do schedule(runtime)
         DO iEl = 1, nelem
            this % p_sem % mesh % elements(iEl) % storage % Q = &
                              this % p_sem % mesh % elements(iEl) % storage % Q + this % MGStorage(iEl) % E
         end DO
!$omp end do
!$omp end parallel
      
      end if
!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!
      sweepcount = 0
      DO
         call this % Smooth(MGSweeps(lvl), t, dt, ComputeTimeDerivative)
         NewRes = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))

         call CFLRamp(cfl_ini,cfl,CurrentMGCycle,PrevRes,NewRes,CFLboost)
         call CFLRamp(dcfl_ini,dcfl,CurrentMGCycle,PrevRes,NewRes,CFLboost)

         sweepcount = sweepcount + MGSweeps(lvl)
         if (MGOutput) call PlotResiduals( lvl, sweepcount , this % p_sem % mesh)
         
         if (sweepcount .ge. MaxSweeps) exit
         
         if (lvl > 1 .and. PostFCycle) then
            if (NewRes > PrevRes) then
               call MGRestrictToChild(this,lvl-1,t, ComputeTimeDerivative)
               call FASVCycle(this,t,dt,lvl-1,lvl, ComputeTimeDerivative)
            else
               exit
            end if
         elseif (PostSmooth .or. PostFCycle) then
            !if (FMG .and. MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) < 0.1_RP) exit
            if (NewRes < PrevRes) exit
         else
            exit
         end if
         
      end DO
      
!
!     -------------------------
!     Compute coarse-grid error
!     -------------------------
!
      if (lvl < MGlevels) then
!$omp parallel do schedule(runtime)
         DO iEl = 1, nelem
            this % MGStorage(iEl) % E = this % p_sem % mesh % elements(iEl) % storage % Q - this % MGStorage(iEl) % Q
         end DO
!$omp end parallel do
      end if
      
#endif 
   end subroutine FASVCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------
!  Recursive subroutine to perform a full multigrid cycle
!  ------------------------------------------------------
   recursive subroutine FASFMGCycle(this,t,tol,lvl, ComputeTimeDerivative)
      implicit none
      !----------------------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this    !<> Current level solver
      real(kind=RP)        , intent(in)    :: t       !<  Simulation time
      real(kind=RP)        , intent(in)    :: tol     !<  Convergence tolerance
      integer              , intent(in)    :: lvl     !<  Current multigrid level
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      !----------------------------------------------------------------------------
      integer        :: iEl, iEQ             ! Element and equation counters
      integer        :: N1(3), N2(3)
      real(kind=RP)  :: maxResidual(NCONS)   ! Maximum residual in each equation
      integer        :: counter              ! Iteration counter
      character(len=LINE_LENGTH) :: FMGFile
      !----------------------------------------------------------------------------
#if defined(NAVIERSTOKES)
!
!     ------------------------------------------
!     At the beginning, go to the coarsest level
!        (the initial condition must be passed)
!     ------------------------------------------
!
      if (lvl > 1) then
!$omp parallel do private(N1,N2) schedule(runtime)
         DO iEl = 1, nelem
            N1 = this         % p_sem % mesh % elements (iEl) % Nxyz
            N2 = this % Child % p_sem % mesh % elements (iEl) % Nxyz
            call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                       N2, this % Child % p_sem % mesh % elements(iEl) % storage % Q )
         end DO
!$omp end parallel do

         call FASFMGCycle(this % Child,t,tol,lvl-1, ComputeTimeDerivative)
      end if
!
!     ------------------------------
!     Save FMG solution if requested
!     ------------------------------
!
      
      if (SaveFMGFile) then
         write(FMGFile,'(A,A,I2.2,A)')  trim( FMGSolutionFile ), '_FMG_', lvl, '.hsol'
         call this % p_sem % mesh % SaveSolution(0,0._RP,trim(FMGFile),saveGradients)
      end if

!
!     ----------------------
!     Perform a V-Cycle here
!     ----------------------
!
      counter = 0
      if (lvl > 1 ) then
         DO
            counter = counter + 1
            call FASVCycle(this,t, own_dt,lvl,lvl, ComputeTimeDerivative) ! FMG is still for STEADY_STATE. TODO: Change that
            maxResidual = ComputeMaxResiduals(this % p_sem % mesh)
            if (maxval(maxResidual) <= tol) exit
         end DO
      else
         DO
            counter = counter + 1
            call this % Smooth(1,t,own_dt,ComputeTimeDerivative)

            maxResidual = ComputeMaxResiduals(this % p_sem % mesh)
            
            if (MOD(counter,100)==0) call PlotResiduals( lvl ,counter, this % p_sem % mesh)
            if (maxval(maxResidual) <= tol) exit
         end DO
      end if
      call PlotResiduals( lvl ,counter, this % p_sem % mesh,.TRUE.)
!
!     --------------------------------------------------
!     If not on finest, Interpolate to next (finer) grid
!     --------------------------------------------------
! 
      if (lvl < MGlevels) then
!$omp parallel do private(N1,N2) schedule(runtime)
         DO iEl = 1, nelem
            N1 = this          % p_sem % mesh % elements (iEl) % Nxyz
            N2 = this % Parent % p_sem % mesh % elements (iEl) % Nxyz
            call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                       N2, this % Parent % p_sem % mesh % elements(iEl) % storage % Q )
         end DO
!$omp end parallel do
      end if
#endif
   end subroutine FASFMGCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine for restricting the solution and residual to the child solver
!  ------------------------------------------
   subroutine MGRestrictToChild(this,lvl,t, ComputeTimeDerivative)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      integer              , intent(IN)    :: lvl
      real(kind=RP)        , intent(IN)    :: t
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      !-------------------------------------------------------------
      class(FASMultigrid_t), pointer       :: Child_p  ! The child
      integer  :: iEl
      integer  :: iEQ
      integer  :: N1(3), N2(3)
      !-------------------------------------------------------------
#if defined(NAVIERSTOKES)      

      Child_p => this % Child

!$omp parallel
!$omp do private(N1,N2) schedule(runtime)
      DO iEl = 1, nelem
         N1 = this    % p_sem % mesh % elements (iEl) % Nxyz
         N2 = Child_p % p_sem % mesh % elements (iEl) % Nxyz
         
!           Restrict solution
!           -----------------
         call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                    N2, Child_p % p_sem % mesh % elements(iEl) % storage % Q )
                                    
!           Restrict residual
!           -----------------
         call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Qdot, &
                                    N2, Child_p % MGStorage(iEl) % S)
      end DO
!$omp end do

!
!     **********************************************************************
!     **********************************************************************
!              Now arrange all the storage in the child solver
!     **********************************************************************
!     **********************************************************************
!

!     ------------------------------------
!     Copy solution from fine grid to MGStorage
!        ... and clear source term
!     ------------------------------------
!
!$omp barrier
!$omp do schedule(runtime)
      DO iEl = 1, nelem
         Child_p % MGStorage(iEl) % Q = Child_p % p_sem % mesh % elements(iEl) % storage % Q
         Child_p % p_sem % mesh % elements(iEl) % storage % S_NS = 0._RP
      end DO
!$omp end do
!$omp end parallel
!
!     -------------------------------------------
!     If not on finest level, correct source term
!     -------------------------------------------
!      
      call ComputeTimeDerivative(Child_p % p_sem % mesh,Child_p % p_sem % particles, t, CTD_IGNORE_MODE) 
      if (DualTimeStepping) call ComputePseudoTimeDerivative(Child_p % p_sem % mesh, t, own_dt)
      
!$omp parallel do schedule(runtime)
      DO iEl = 1, nelem
         Child_p % p_sem % mesh % elements(iEl) % storage % S_NS = Child_p % MGStorage(iEl) % S - &
                                                                Child_p % p_sem % mesh % elements(iEl) % storage % Qdot
      end DO
!$omp end parallel do
      
#endif      
   end subroutine MGRestrictToChild
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that destructs a FAS integrator
!  ------------------------------------------
   subroutine destruct(this)       
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      !-----------------------------------------------------------
      
      call RecursiveDestructor(this,MGlevels)
      
   end subroutine destruct
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine RecursiveDestructor(Solver,lvl)
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: Solver
      integer                              :: lvl
      !-----------------------------------------------------------
      
      ! First go to coarsest level
      if (lvl > 1) call RecursiveDestructor(Solver % Child,lvl-1)

      ! Deallocate LTS
      if (allocated(Solver % lts_dt)) deallocate(Solver % lts_dt)
      
      !Destruct Multigrid storage
      deallocate (Solver % MGStorage) ! allocatable components are automatically deallocated
      
      ! Destruct linear solver (if present)
      if (Smoother >= IMPLICIT_SMOOTHER_IDX) then
         call Solver % linsolver % destroy
         deallocate ( Solver % linsolver )
      end if
      
      if (lvl < MGlevels) then
         call Solver % p_sem % destruct()
         deallocate (Solver % p_sem)
         nullify    (Solver % Parent)
      else
         nullify    (Solver % p_sem)
      end if
      
      if (lvl > 1) then
         deallocate (Solver % Child)
      end if
      
   end subroutine RecursiveDestructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Smooth(this,SmoothSweeps,t,dt, ComputeTimeDerivative)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t)  , intent(inout), target :: this     !<> Anisotropic FAS multigrid 
      integer                , intent(in)            :: SmoothSweeps
      real(kind=RP)          , intent(in)            :: t
      real(kind=RP)          , intent(in)            :: dt
      procedure(ComputeTimeDerivative_f)                     :: ComputeTimeDerivative
      !-------------------------------------------------------------
      real(kind=RP), pointer :: smoother_dt, smoother_cfl, smoother_dcfl
      integer :: sweep
      !-------------------------------------------------------------
      
      if (DualTimeStepping) then
         smoother_dt   => p_dt
         smoother_cfl  => p_cfl
         smoother_dcfl => p_dcfl
      else 
         smoother_dt   => own_dt
         smoother_cfl  => cfl
         smoother_dcfl => dcfl
      end if

      select case (Preconditioner)
      case (PRECONDIIONER_LTS)
         ! compute LTS
         if (Compute_dt) then
            call MaxTimeStep( self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl , MaxDt=smoother_dt, MaxDtVec=this % lts_dt)
         else
            error stop "FASMultigrid :: LTS needs cfd & dcfl."
         end if

         select case (Smoother)
            ! Euler Smoother
            case (Euler_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  ! call TakeExplicitEulerStep (this % p_sem % mesh, this % p_sem % particles, t, &
                  !    smoother_dt, ComputeTimeDerivative, this % lts_dt )
                  call TakeExplicitEulerStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=own_dt )
               end do
!
!           3rd order Runge-Kutta smoother
!           -> Has its own dt, since it's for steady-state simulations
!           ----------------------------------------------------------
            case (RK3_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  call TakeRK3Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=own_dt )
               end do
            ! RK5 smoother
            case (RK5_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  call TakeRK5Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=own_dt )
               end do
            ! RK5 smoother opt for Steady State
            case (RKOpt_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  call TakeRKOptStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, N_STAGES=erk_order, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=own_dt )
               end do
            case default
               error stop "FASMultigrid :: No smoother defined for the multigrid."
         end select ! Smoother
      case (PRECONDIIONER_NONE)
         select case (Smoother)
            ! Euler Smoother
            case (Euler_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  call TakeExplicitEulerStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=own_dt )
               end do
!
!           3rd order Runge-Kutta smoother
!           -> Has its own dt, since it's for steady-state simulations
!           ----------------------------------------------------------
            case (RK3_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )
                  call TakeRK3Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=own_dt )
               end do
            ! RK5 smoother
            case (RK5_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )
                  call TakeRK5Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=own_dt )
               end do
            ! RK5 opt smoother
            case (RKOpt_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )
                  call TakeRKOptStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, N_STAGES=erk_order, dts=DualTimeStepping, global_dt=own_dt )
               end do
!
!           Implicit smoothers
!           ------------------
            ! case (666)
            case (BJ_SMOOTHER)
               call this % p_sem % mesh % storage % local2GlobalQdot(this % p_sem % NDOF)
               call this % p_sem % mesh % storage % local2GlobalQ(this % p_sem % NDOF)

               ! print *, "QNS: ", this % p_sem % mesh % storage % QNS(1)
               ! print *, "Q  : ", this % p_sem % mesh % elements(1) % storage % Q(1,1,1,1)

               call ComputeRHS(this % p_sem, t, dt, this % linsolver, ComputeTimeDerivative )

               do sweep = 1, SmoothSweeps
                  call TakeBJSweep (this , t, &
                                ComputeTimeDerivative )
               end do
               call this % p_sem % mesh % storage % global2LocalQ

               ! print *, "QNS: ", this % p_sem % mesh % storage % QNS(1)
               ! print *, "Q  : ", this % p_sem % mesh % elements(1) % storage % Q(1,1,1,1)
               ! error stop "TBC"
            case default
               call ComputeRHS(this % p_sem, t, dt, this % linsolver, ComputeTimeDerivative )               ! Computes b (RHS) and stores it into linsolver

!~               this % computeA = .TRUE.
               call this % linsolver % solve(NCONS, NGRAD, maxiter=SmoothSweeps, time= t, dt = dt, &
                                                ComputeTimeDerivative = ComputeTimeDerivative, computeA = this % computeA) ! 
               call UpdateNewtonSol(this % p_sem, this % linsolver)
         end select ! Smoother
      end select ! Preconditioner
      
   end subroutine Smooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine FAS_SetPreviousSolution(this,lvl)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), target, intent(inout) :: this     !<> Anisotropic FAS multigrid 
      integer                      , intent(in)    :: lvl      !<  Current level
      !-------------------------------------------------------------
      integer :: N1(3), N2(3), eID
      !-------------------------------------------------------------
      
!
!     Set the previous solution in this level
!     ---------------------------------------
      
      if (lvl == MGlevels) then
         call BDF_SetPreviousSolution(this % p_sem % mesh % storage)
      else
         call BDF_SetPreviousSolution(this % p_sem % mesh % storage, NotANewStep = .TRUE. )
      end if
!
!     Send the solution to the next (coarser) level
!     ---------------------------------------------

      if (lvl > 1) then
!$omp parallel do private(N1,N2) schedule(runtime)
         do eID = 1, nelem
            N1 = this % p_sem % mesh % elements (eID) % Nxyz
            N2 = this % Child % p_sem % mesh % elements (eID) % Nxyz
         
!           Restrict solution
!           -----------------
            call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(eID) % storage % Q, &
                                       N2, this % Child % p_sem % mesh % elements(eID) % storage % Q )
         end do
!$omp end parallel do
         
         call FAS_SetPreviousSolution(this % Child,lvl-1)
      end if
   end subroutine FAS_SetPreviousSolution

   SUBROUTINE TakeBJSweep( this, t, ComputeTimeDerivative )
         implicit none
         class(FASMultigrid_t), target, intent(inout) :: this     !<> Anisotropic FAS multigrid 
         REAL(KIND=RP)                   :: t, tk
         procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
         !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer                    :: id, k
         real(kind=RP) :: dalpha=1e-3_RP

         associate ( mesh => this % p_sem % mesh)
   
         mesh % storage % PrevQNS(:,1) = mesh % storage % QNS

         select type (GMRES => this % linsolver)
            type is (MatFreeGMRES_t)
               ! call GMRES % BlockPreco % SolveBlocks_LU(mesh % storage % QNS,-mesh % storage % QdotNS)
               call GMRES % BlockPreco % SolveBlocks_LU(mesh % storage % QNS, GMRES % RHS)
         end select 

         DO id = 1, SIZE( mesh % storage % QNS )
            dalpha = min(0.1_RP * abs( mesh % storage % PrevQNS(id,1) / mesh % storage % QNS(id) ),1.0)
         END DO

         ! print *, mesh % storage % PrevQNS(1,1)
         ! print *, mesh % storage % QNS(1)

         mesh % storage % QNS = mesh % storage % PrevQNS(:,1) + dalpha * mesh % storage % QNS

         ! print *, mesh % storage % QNS(1)

         end associate
   
         end subroutine TakeBJSweep

end module FASMultigridClass