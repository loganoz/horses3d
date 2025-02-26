!
!   Module for integrating in time using the Backward Differentiation Formulas (B  D   F)
!     using the split steps form
!
!////////////////////////////////////////////////////////////////////////
MODULE SplitStepsTimeIntegrator
   use SMConstants
   USE PhysicsStorage
   use HexMeshClass
   USE LinearSolverClass
   USE CSRMatrixClass
   USE FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use DGSEMClass
   use StorageClass              , only: SolutionStorage_t
   use StopwatchClass            , only: Stopwatch
   use Utilities                 , only: AlmostEqual
   use MPI_Process_Info          , only: MPI_Process
   implicit none
   
   PRIVATE                          
   PUBLIC SplitStepsIntegrator_t, TakeSplitStepsStep, ComputeRHSSteps2, UpdateNewtonSol, SplitSteps_SetPreviousSolution, bdf_order, SplitSteps_MatrixShift, SplitSteps_SetOrder
   public SplitStepsCoeff, SplitStepsInitialiseQ
   public ComputeAdvectStep1
   public SplitSteps_MatrixShift_step2, SplitSteps_MatrixShift_step3
   
!
!  ********************
!  SplitSteps integrator class
!  ********************
   type SplitStepsIntegrator_t
      
      class(GenericLinSolver_t), allocatable :: linsolver     !  Linear solver
      class(GenericLinSolver_t), allocatable :: linsolverStep2     !  Linear solver
      class(GenericLinSolver_t), allocatable :: linsolverStep3     !  Linear solver
      integer                                :: maxLinSolIter
      real(kind=RP)                          :: firstNorm
      real(kind=RP)                          :: NewtonFactor
      real(kind=RP)                          :: LinSolTolFactor
      integer                                :: StepsForJac   !· Maximum number of steps that should be taken for computing a new Jacobian matrix
      integer                                :: StepsSinceJac !  
      integer                                :: MaxNewtonIter
      logical                                :: JacByConv     !· .TRUE. if the Jacobian must be computed only when the convergence is bad
      logical                                :: TimeAccurate  !· .TRUE. if this is a time-accurate simulation
      logical                                :: UserNewtonTol !· .TRUE. if the newton tolerance is specified by the user
      real(kind=RP)                          :: NewtonTol     !  Specified Newton tolerance
      real(kind=RP)                          :: inner_dt   
      
      contains
         procedure :: construct => ConstructSplitStepsIntegrator
         procedure :: destruct  => DestructSplitStepsIntegrator
         procedure :: TakeStep  => TakeSplitStepsStep
         procedure :: TakeSplit_3_Steps
   end type SplitStepsIntegrator_t
   
!
!  Module variables
!  ----------------
   
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
   logical       :: computeA_step2 = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
   logical       :: computeA_step3 = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
   logical       :: Adaptive_dt = .TRUE.
   ! integer       :: bdf_order       ! SplitSteps order specified by user
   integer       :: order           ! SplitSteps order to be used
   integer       :: StepsTaken = 0

!
!  SplitSteps coefficients for constant time-step
!  ---------------------------------------
   real(kind=RP), parameter :: SplitStepsCoeff(6,5) = &
!                    a_1             a_2     a_3           a_4             a_5          a_6
         reshape( (/ 1.0_RP        , -1._RP, 0._RP       , 0._RP         , 0._RP      , 0._RP        ,  &   ! SplitSteps1
                     1.5_RP        , -2._RP, 0.5_RP      , 0._RP         , 0._RP      , 0._RP        ,  &   ! SplitSteps2
                     11._RP/6_RP   , -3._RP, 3._RP/2._RP , -1._RP/3._RP  , 0._RP      , 0._RP        ,  &   ! SplitSteps3
                     25._RP/12_RP  , -4._RP, 3._RP       , -4._RP/3._RP  , 1._RP/4._RP, 0._RP        ,  &   ! SplitSteps4
                     137._RP/60_RP , -5._RP, 5._RP       , -10._RP/3._RP , 5._RP/4._RP, -1._RP/5._RP /) &   ! SplitSteps5
                                                                                                      , (/6,5/) )
   integer, parameter :: MAX_ORDER_CONS_DT = 5
   
!
!  Default parameters for Newton iterative procedure
!  -------------------------------------------------
   real(kind=RP), parameter   :: NEWTON_MIN_CONVRATE = 0.1_RP     ! Minimum convergence rate for Newton method... If newton loop convergence rate passes this value, inner_dt is decreased
   real(kind=RP), parameter   :: NEWTON_MAX_CONVRATE = 1.7_RP     ! Maximum convergence rate for Newton method... If newton loop convergence rate passes this value, inner_dt is increased
   real(kind=RP), parameter   :: NEWTON_TOL_DEFAULT  = 1.e-6_RP   ! Default convergence tolerance
   integer      , parameter   :: MAX_NEWTON_ITER = 30          ! If newton iter reaches this limit, this iteration is marked as  not converged 
   integer      , parameter   :: LIM_NEWTON_ITER = 12          ! If Newton converges but this limit is reached, jacobian matrix will be recomputed
   logical                    :: PRINT_NEWTON_INFO
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine ConstructSplitStepsIntegrator(this,controlVariables,sem)
      implicit none
      !--------------------------------------------------------
      class(SplitStepsIntegrator_t) , intent(inout) :: this
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(DGSem)            , intent(in)    :: sem
      !--------------------------------------------------------
      integer :: DimPrb, globalDimPrb
      !--------------------------------------------------------
      
!
!     Get general definitions
!     -----------------------

      Adaptive_dt = controlVariables % logicalValueForKey("implicit adaptive dt")
      this % JacByConv = controlVariables % LogicalValueForKey("jacobian by convergence")
      if (controlVariables % StringValueForKey("simulation type",LINE_LENGTH) == 'time-accurate') then
         this % TimeAccurate = .TRUE.
      else
         this % TimeAccurate = .FALSE.
      end if
      
      if (controlVariables % containsKey("newton max iter")) then
         this % MaxNewtonIter = controlVariables % integerValueForKey("newton max iter")
      else
         this % MaxNewtonIter = MAX_NEWTON_ITER
      end if
      
      if (controlVariables % containsKey("linsolver max iter")) then
         this % maxLinSolIter = controlVariables % integerValueForKey("linsolver max iter")
      else
         this % maxLinSolIter = 500
      end if

      if (controlVariables % containsKey("newton first norm")) then
         this % firstNorm = controlVariables % doublePrecisionValueForKey("newton first norm")
      else
         this % firstNorm = 2.d-1
      end if

      if (controlVariables % containsKey("newton factor")) then
         this % NewtonFactor = controlVariables % doublePrecisionValueForKey("newton factor")
      else
         this % NewtonFactor = 1e-3_RP
      end if

      if (controlVariables % containsKey("linsolver tol factor")) then
         this % LinSolTolFactor = controlVariables % doublePrecisionValueForKey("linsolver tol factor")
      else
         this % LinSolTolFactor = 0.5_RP
      end if

      PRINT_NEWTON_INFO = controlVariables % logicalValueForKey("print newton info") .and. MPI_Process % isRoot

      if (controlVariables % containsKey("newton tolerance")) then
         this % UserNewtonTol = .TRUE.
         this % NewtonTol = controlVariables % doublePrecisionValueForKey("newton tolerance")
      else
         this % UserNewtonTol = .FALSE.
         this % NewtonTol = NEWTON_TOL_DEFAULT
      end if
      
      this % StepsForJac = controlVariables % integerValueForKey("compute jacobian every")
!
!     Setup linear solver
!     -------------------
      DimPrb = sem % NDOF * NCONS
      globalDimPrb = sem % totalNDOF * NCONS
      
      select case ( trim(controlVariables % StringValueForKey("linear solver",LINE_LENGTH)) )
         ! case('petsc')
         !    allocate (PetscKspLinearSolver_t :: this % linsolver)
         case('pardiso')
            allocate (MKLPardisoSolver_t     :: this % linsolver)
            allocate (MKLPardisoSolver_t     :: this % linsolverStep2)
            allocate (MKLPardisoSolver_t     :: this % linsolverStep3)
            call this % linsolverStep2 % construct (sem % NDOF,         sem % totalNDOF,              1, controlVariables, sem, SplitSteps_MatrixShift_step2)             !Constructs linear solver 
            call this % linsolverStep3 % construct (sem % NDOF * N_INS, sem % totalNDOF * N_INS,  N_INS, controlVariables, sem, SplitSteps_MatrixShift_step3)             !Constructs linear solver 

         ! case('matrix-free smooth')
         !    allocate (MatFreeSmooth_t        :: this % linsolver)
         ! case('matrix-free gmres')
         !    allocate (MatFreeGMRES_t         :: this % linsolver)
         ! case('smooth')
         !    allocate (IterativeSolver_t      :: this % linsolver)
         ! case('multigrid')
         !    allocate (LinearMultigridSolver_t      :: this % linsolver)
         ! case('static-condensation')
         !    allocate (StaticCondSolver_t     :: this % linsolver)
         case default
            write(STD_OUT,*) "Keyword 'linear solver' missing... Using PETSc as default"
            ! allocate (PetscKspLinearSolver_t :: this % linsolver)
      end select
      
      ! call this % linsolver % construct (DimPrb,globalDimPrb,NCONS,controlVariables,sem,SplitSteps_MatrixShift_step3)             !Constructs linear solver 
      call this % linsolver % construct (DimPrb,globalDimPrb,NCONS,controlVariables,sem,SplitSteps_MatrixShift)             !Constructs linear solver 
      
!
!     Setup SplitSteps methods
!     -----------------
      call SplitSteps_SetOrder( controlVariables % integerValueForKey("SplitSteps order") )
      
      ! Check that the SplitSteps order is consistent
      if (bdf_order > 1) then
         if ( (.not. controlVariables % containsKey("dt") ) .or. Adaptive_dt) then
            error stop ':: "SplitSteps order">1 is only valid with fixed time-step sizes'
         end if
      end if
      
!     Create Stopwatch for solving
!     ----------------------------
      call Stopwatch % CreateNewEvent("SplitSteps Newton-Solve")
      
   end subroutine ConstructSplitStepsIntegrator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   subroutine DestructSplitStepsIntegrator(this)
      implicit none
      !--------------------------------------------------------
      class(SplitStepsIntegrator_t), intent(inout) :: this
      !--------------------------------------------------------
      
      call this % linsolver % destroy
      call this % linsolverStep2 % destroy
      call this % linsolverStep3 % destroy

      deallocate (this % linsolver)
      deallocate (this % linsolverStep2)
      deallocate (this % linsolverStep3)
      
   end subroutine DestructSplitStepsIntegrator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   subroutine TakeSplitStepsStep (this, sem, t , dt, ComputeTimeDerivative)
      implicit none
      !----------------------------------------------------------------------
      class(SplitStepsIntegrator_t), intent(inout) :: this
      type(DGSem),                  INTENT(inout)           :: sem                  !<>DGSem class with solution storage 
      real(kind=RP),                INTENT(in)              :: t                    !< Time at the beginning of time step
      real(kind=RP),                INTENT(in)              :: dt                   !< Initial (outer) time step (the subroutine can use a smaller one depending on convergence)
      procedure(ComputeNonlinearStep1_f)                            :: ComputeTimeDerivative
      !----------------------------------------------------------------------
      
      real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
      integer                                               :: k, newtonit
      
      real(kind=RP)                                         :: ConvRate
      logical                                               :: CONVERGED
      !----------------------------------------------------------------------
      
      if ((.not. this % TimeAccurate) .and. (.not. this % UserNewtonTol)) then
         this % NewtonTol = sem % MaxResidual* this % NewtonFactor
      end if
      
      !**************************
      ! If the Jacobian must only be computed sometimes
      if (this % JacByConv) then
         if ( (.not. computeA) .and. (.not. AlmostEqual (this % inner_dt,dt) ) ) then
            call this % linsolver % ReSetOperatorDt(this % inner_dt)
         end if
      end if
      ! 
      !**************************
      
      this % inner_dt = dt            ! first inner_dt is the outer step dt 
      time = t
      ! write (*,*) 1_RP/(t-t)
      
      call sem % mesh % storage % local2GlobalQ(sem % NDOF)
!
!     ********************
!     Sub-time-step solver
!     ********************
      do
         
!        Set previous solution for inner time-step
!        -----------------------------------------
         
         call SplitSteps_SetPreviousSolution(sem % mesh % storage)

         ! write (*,*) "--------sem % mesh % storage % prevSol_num-------------", sem % mesh % storage % prevSol_num


         ! ComputeAdvectStep1(sem, t, dt, this % linsolver, ComputeNonlinearStep1)
  
!        Perform Newton iterative procedure
!        ----------------------------------
         
         if (computeA) then
            this % StepsSinceJac = 0
         else
            this % StepsSinceJac = this % StepsSinceJac + 1
            if (this % StepsSinceJac == this % StepsForJac) then
               computeA = .TRUE.
               this % StepsSinceJac = 0
            end if
         end if

         write (*,*) "DirectSolve(sem, time+this % inner_dt, this % inner_dt, this % linsolver, this % NewtonTol, ***"

         call DirectSolve(sem, time+this % inner_dt, this % inner_dt, this % linsolver, this % NewtonTol, &
                           this % MaxNewtonIter, this % maxLinSolIter,                                      &
                           this % firstNorm, this % LinSolTolFactor,                                        &
                          this % JacByConv,ConvRate, newtonit,CONVERGED, ComputeTimeDerivative)

         time = time + this % inner_dt

         write (*,*) "=== time = time + this % inner_dt ===========",time
                                     
!          Check if the sub time-stepping is done
!          --------------------------------------
         if (ABS((time)-(t+dt)) < 10 * EPSILON(1._RP)) then       ! If outer t+dt is reached, the time integration is done
            exit                                            
         end if
       
      end do
 
      if (PRINT_NEWTON_INFO) write(*,'(A10,f5.2)') "ConvRate: ", ConvRate
      
      !**************************
      ! for computing sometimes
      if (this % JacByConv .AND. ConvRate <0.65_RP ) then
         computeA = .TRUE.
      end if
      ! for computing sometimes
      !**************************
      
      
      call sem % mesh % storage % global2LocalQ
      call sem % mesh % storage % global2LocalQdot
      
   end subroutine TakeSplitStepsStep
!   
! =======================================================================***********************************************
! =======================================================================***********************************************
! =======================================================================***********************************************

   subroutine TakeSplit_3_Steps (this, sem, t , dt, ComputeTimeDerivative, ComputeTimeDerivative_Second, &
                                 ComputeTimeDerivative_Third, ComputeNonlinearStep1 &
                                 ,kNum, FirstStep, SecondStep                                      & 
      )
      implicit none
      !----------------------------------------------------------------------
      class(SplitStepsIntegrator_t), intent(inout) :: this
      type(DGSem),                  INTENT(inout)           :: sem                  !<>DGSem class with solution storage 
      real(kind=RP),                INTENT(in)              :: t                    !< Time at the beginning of time step
      real(kind=RP),                INTENT(in)              :: dt                   !< Initial (outer) time step (the subroutine can use a smaller one depending on convergence)
      procedure(ComputeTimeDerivative_f)     :: ComputeTimeDerivative
      procedure(ComputeNonlinearStep1_f)     :: ComputeTimeDerivative_Second
      procedure(ComputeNonlinearStep1_f)     :: ComputeTimeDerivative_Third
      procedure(ComputeNonlinearStep1_f)     :: ComputeNonlinearStep1
      integer       , optional      , INTENT(in)       :: kNum
      logical      , optional      , intent(inout)          :: FirstStep
      logical      , optional      , intent(inout)          :: SecondStep
      !----------------------------------------------------------------------
      
      real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
      integer                                               :: k, newtonit
      integer     :: iterNum
      real(kind=RP) :: gamma, alpha0, alpha1, beta0, beta1         

      
      real(kind=RP)                                         :: ConvRate
      logical                                               :: CONVERGED
      !----------------------------------------------------------------------
      
      if ((.not. this % TimeAccurate) .and. (.not. this % UserNewtonTol)) then
         this % NewtonTol = sem % MaxResidual* this % NewtonFactor
      end if
      
      iterNum = 0 
      ! !**************************
      ! ! If the Jacobian must only be computed sometimes
      ! if (this % JacByConv) then
      !    if ( (.not. computeA) .and. (.not. AlmostEqual (this % inner_dt,dt) ) ) then
      !       call this % linsolver % ReSetOperatorDt(this % inner_dt)
      !    end if
      ! end if

      ! computeA_step2
      ! 
      !**************************
      
      this % inner_dt = dt            ! first inner_dt is the outer step dt 
      time = t
      ! write (*,*) 1_RP/(t-t)
      
      call sem % mesh % storage % local2GlobalQ(sem % NDOF)

      ConvRate = 1.0_RP
      ! call sem % mesh % storage % local2GlobalQ(sem % NDOF)


      write (*,*), "iterNum ==0========= ", iterNum, kNum

      if (FirstStep .eqv. .True.) then
         gamma  =  1.0_RP 
         alpha0 =  1.0_RP
         alpha1 =  0.0_RP 
         beta0  =  1.0_RP
         beta1  =  0.0_RP
      else if (.not. FirstStep) then
         ! gamma  =  1.5_RP 
         ! alpha0 =  2.0_RP
         ! alpha1 = -0.5_RP 
         ! beta0  =  2.0_RP
         ! beta1  = -1.0_RP

         gamma  =  1.0_RP 
         alpha0 =  1.0_RP
         alpha1 =  0.0_RP 
         beta0  =  1.0_RP
         beta1  =  0.0_RP
      endif

      if (kNum == 1) then 
         computeA_step3 = .true.
      endif

!
!     ********************
!     Sub-time-step solver
!     ********************
      do
         write (*,*), "iterNum ==1========= ", iterNum+1, kNum, sem% numberOfTimeSteps, FirstStep
         
!        Set previous solution for inner time-step
!        -----------------------------------------
         
         call SplitSteps_SetPreviousSolution(sem % mesh % storage)

         ! write (*,*) "--------sem % mesh % storage % prevSol_num-------------", sem % mesh % storage % prevSol_num


         call ComputeAdvectStep1(sem, t, dt, this % linsolver, ComputeNonlinearStep1, gamma, alpha0, alpha1, beta0, beta1)

         call sem % mesh % storage % local2GlobalQ(sem % NDOF)

         ! =====================================================================================================================
         ! =====================================================================================================================

         ! write (*,*) "=========computeA 1==============", computeA_step2
!          ! =====================================================================================================================

         ! computeA_step2 = .True.
         if (computeA_step2) then
            this % StepsSinceJac = 0
         else
            this % StepsSinceJac = this % StepsSinceJac + 1
            if (this % StepsSinceJac == this % StepsForJac) then
               computeA_step2 = .TRUE.
               this % StepsSinceJac = 0
            end if
         end if

         ! write (*,*) "DirectSolveSecond(sem, time+this % inner_dt, this % inner_dt, this % linsolver, this % NewtonTol, ***"
         ! write (*,*) "=========computeA 2==============", computeA_step2

         call this % linsolverStep2 % MatrixShift_T(0.0_RP)
         
         call DirectSolveSecond(sem, time+this % inner_dt, this % inner_dt, this % linsolverStep2, this % NewtonTol, &
                           this % MaxNewtonIter, this % maxLinSolIter,                                      &
                           this % firstNorm, this % LinSolTolFactor,                                        &
                          this % JacByConv,ConvRate, newtonit,CONVERGED, ComputeTimeDerivative_Second,  &
                          gamma, alpha0, alpha1, beta0, beta1         &
                          , computeA_step2   &
                          )
         ! =====================================================================================================================
         ! =====================================================================================================================
         ! write (*,*) "=========computeA 3==============", computeA
         ! write (*,*) "=========computeA_step2 3==============", computeA_step2
         ! write (*,*) "=========computeA_step3 3==============", computeA_step3


         ! computeA_step3 = .True.
         if (computeA_step3) then
            this % StepsSinceJac = 0
         else
            this % StepsSinceJac = this % StepsSinceJac + 1
            if (this % StepsSinceJac == this % StepsForJac) then
               computeA_step3 = .TRUE.
               this % StepsSinceJac = 0
            end if
         end if

         ! write (*,*) "=========computeA 4==============", computeA_step3

         ! write (*,*) "DirectSolveThird(sem, time+this % inner_dt, this % inner_dt, this % linsolver, this % NewtonTol, ***"

         ! call this % linsolverStep3 % MatrixShift_T(0.0_RP)
         ! call this % linsolverStep3 % MatrixShift_T(1.0_RP/dt)
         call this % linsolverStep3 % MatrixShift_T(1.0_RP*gamma/dt)

         call DirectSolveThird(sem, time+this % inner_dt, this % inner_dt, this % linsolverStep3, this % NewtonTol, &
                           this % MaxNewtonIter, this % maxLinSolIter,                                      &
                           this % firstNorm, this % LinSolTolFactor,                                        &
                          this % JacByConv,ConvRate, newtonit,CONVERGED, ComputeTimeDerivative_Third, &
                          gamma, alpha0, alpha1, beta0, beta1         &
                          , computeA_step3                               &
                          )
!          ! =====================================================================================================================
         ! write (*,*) "=========computeA 5==============", computeA_step3


!          ! =====================================================================================================================
!          ! =====================================================================================================================
!          ! =====================================================================================================================

         iterNum = iterNum + 1
         time = time + this % inner_dt

         write (*,*) "=== time = time + this % inner_dt, iteration numbers ===========",time, iterNum
                                     
!          Check if the sub time-stepping is done
!          --------------------------------------
         if (ABS((time)-(t+dt)) < 10 * EPSILON(1._RP)) then       ! If outer t+dt is reached, the time integration is done
            exit                                            
         end if
      
      end do
 
      if (PRINT_NEWTON_INFO) write(*,'(A10,f5.2)') "ConvRate: ", ConvRate
      
      !**************************
      ! for computing sometimes
      if (this % JacByConv .AND. ConvRate <0.65_RP ) then
         computeA = .TRUE.
      end if
      ! for computing sometimes
      !**************************
      

     
      
      call sem % mesh % storage % global2LocalQ

      call sem % mesh % storage % SolutionStorage_global2LocalPress
      call sem % mesh % storage % SolutionStorage_global2LocalVelINS


      call sem % mesh % storage % global2LocalQdot
      
   end subroutine TakeSplit_3_Steps

   !
!  -------------------------------------------------------------
!  Routine for performing a nonlinear Newton iterative procedure
!  -> This can be taken out of the SplitStepsTimeIntegrator if needed
!        (but careful with Adaptive_dt and the Newton vars)
!  -------------------------------------------------------------
   subroutine DirectSolveSecond(sem, t, dt, linsolver, NEWTON_TOLERANCE, MaxNewtonIter, maxLinSolIter, firstNorm,  &
                                 LinSolTolFactor, JacByConv,ConvRate, niter,CONVERGED, ComputeTimeDerivative_Second, &
                                 gamma, alpha0, alpha1, beta0, beta1,  &
                                 computeA)
      implicit none
      !----------------------------------------------------------------------
      type(DGSem),                  intent(inout)           :: sem
      real(kind=RP),                intent(in)              :: t
      real(kind=RP),                intent(in)              :: dt              !< Inner dt
      class(GenericLinSolver_t),    intent(inout)           :: linsolver       !Linear operator is calculate outside this subroutine
      real(kind=RP),                intent(in)              :: NEWTON_TOLERANCE
      integer,                      intent(in)              :: MaxNewtonIter
      integer,                      intent(in)              :: maxLinSolIter
      real(kind=RP),                intent(in)              :: firstNorm
      real(kind=RP),                intent(in)              :: LinSolTolFactor
      logical,                      intent(in)              :: JacByConv         !< Must the Jacobian be computed for bad convergence? if .false., the Jacobian is computed at the beginning of every newton it
      real(kind=RP),                intent(out)             :: ConvRate
      integer,                      intent(out)             :: niter
      logical,                      intent(out)             :: CONVERGED   
      procedure(ComputeNonlinearStep1_f)                    :: ComputeTimeDerivative_Second
      logical      , optional      , intent(inout)          :: ComputeA
      real(kind=RP),              intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         
      !----------------------------------------------------------------------           
      integer              :: newtonit
      real(kind=RP)        :: norm, norm_old, rel_tol
      real(kind=RP)        :: linsolver_tol
      logical      , save  :: isfirst = .TRUE.
      real(kind=RP), save  :: norm1
      !----------------------------------------------------------------------
      
!
!     Initializations
!     ---------------
      
      if (isfirst) then
         norm = firstNorm  ! A value to define the initial linsolver_tol
         isfirst = .FALSE.
      else
         norm = norm1
      end if

      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      ConvRate = 1.0_RP
   
!     -----------
         call ComputeRHSSteps2Second(sem, t, dt, linsolver, ComputeTimeDerivative_Second,gamma, alpha0, alpha1, beta0, beta1 )               ! Computes b (RHS) and stores it into linsolver
         
         
         call linsolver%solve( nEqn=1, nGradEqn=1, tol = linsolver_tol, maxiter=maxLinSolIter, time= t, dt=dt, &
                              ComputeTimeDerivative = ComputeTimeDerivative_Second, computeA = computeA, startNum=4)     


         sem % mesh % storage % pressINS = linsolver % GetX()

         norm = linsolver%Getxnorm('infinity')
         ! write (*,*) "norm ============== ** ==============",norm
         
   
   end subroutine DirectSolveSecond


    !
!  -------------------------------------------------------------
!  Routine for performing a nonlinear Newton iterative procedure
!  -> This can be taken out of the SplitStepsTimeIntegrator if needed
!        (but careful with Adaptive_dt and the Newton vars)
!  -------------------------------------------------------------
   subroutine DirectSolveThird(sem, t, dt, linsolver, NEWTON_TOLERANCE, MaxNewtonIter, maxLinSolIter, &
      firstNorm, LinSolTolFactor, JacByConv,ConvRate, niter,CONVERGED, ComputeTimeDerivative_Third,&
      gamma, alpha0, alpha1, beta0, beta1         ,  &
      computeA)
      implicit none
      !----------------------------------------------------------------------
      type(DGSem),                  intent(inout)           :: sem
      real(kind=RP),                intent(in)              :: t
      real(kind=RP),                intent(in)              :: dt              !< Inner dt
      class(GenericLinSolver_t),    intent(inout)           :: linsolver       !Linear operator is calculate outside this subroutine
      real(kind=RP),                intent(in)              :: NEWTON_TOLERANCE
      integer,                      intent(in)              :: MaxNewtonIter
      integer,                      intent(in)              :: maxLinSolIter
      real(kind=RP),                intent(in)              :: firstNorm
      real(kind=RP),                intent(in)              :: LinSolTolFactor
      logical,                      intent(in)              :: JacByConv         !< Must the Jacobian be computed for bad convergence? if .false., the Jacobian is computed at the beginning of every newton it
      real(kind=RP),                intent(out)             :: ConvRate
      integer,                      intent(out)             :: niter
      logical,                      intent(out)             :: CONVERGED   
      procedure(ComputeNonlinearStep1_f)                    :: ComputeTimeDerivative_Third
      logical      , optional      , intent(inout)          :: ComputeA
      real(kind=RP),              intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         
      !----------------------------------------------------------------------           
      integer              :: newtonit
      real(kind=RP)        :: norm, norm_old, rel_tol
      real(kind=RP)        :: linsolver_tol
      logical      , save  :: isfirst = .TRUE.
      real(kind=RP), save  :: norm1
      !----------------------------------------------------------------------
      
!
!     Initializations
!     ---------------
      
      if (isfirst) then
         norm = firstNorm  ! A value to define the initial linsolver_tol
         isfirst = .FALSE.
      else
         norm = norm1
      end if

      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      ConvRate = 1.0_RP
   
!     -----------
         call ComputeRHSSteps3Third(sem, t, dt, linsolver, ComputeTimeDerivative_Third,gamma, alpha0, alpha1, beta0, beta1 )               ! Computes b (RHS) and stores it into linsolver
         
         ! write(*,*) " hello  DirectSolveThirdr SplitStepsTimeIntegrator.f90, after ComputeRHSSteps2Third ============"
         
         call linsolver%solve( nEqn=N_INS, nGradEqn=N_INS, tol = linsolver_tol, maxiter=maxLinSolIter, time= t, dt=dt, &
                              ComputeTimeDerivative = ComputeTimeDerivative_Third, computeA = computeA, startNum=5)     


         sem % mesh % storage % vel__INS = linsolver % GetX()

         norm = linsolver%Getxnorm('infinity')
         ! write (*,*) "norm ============== * Third * ==============",norm
         
   
   end subroutine DirectSolveThird
!/////////////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------
!  Routine for performing a nonlinear Newton iterative procedure
!  -> This can be taken out of the SplitStepsTimeIntegrator if needed
!        (but careful with Adaptive_dt and the Newton vars)
!  -------------------------------------------------------------
   subroutine DirectSolve(sem, t, dt, linsolver, NEWTON_TOLERANCE, MaxNewtonIter, maxLinSolIter, firstNorm, LinSolTolFactor, JacByConv,ConvRate, niter,CONVERGED, ComputeTimeDerivative)
      implicit none
      !----------------------------------------------------------------------
      type(DGSem),                  intent(inout)           :: sem
      real(kind=RP),                intent(in)              :: t
      real(kind=RP),                intent(in)              :: dt              !< Inner dt
      class(GenericLinSolver_t),    intent(inout)           :: linsolver       !Linear operator is calculate outside this subroutine
      real(kind=RP),                intent(in)              :: NEWTON_TOLERANCE
      integer,                      intent(in)              :: MaxNewtonIter
      integer,                      intent(in)              :: maxLinSolIter
      real(kind=RP),                intent(in)              :: firstNorm
      real(kind=RP),                intent(in)              :: LinSolTolFactor
      logical,                      intent(in)              :: JacByConv         !< Must the Jacobian be computed for bad convergence? if .false., the Jacobian is computed at the beginning of every newton it
      real(kind=RP),                intent(out)             :: ConvRate
      integer,                      intent(out)             :: niter
      logical,                      intent(out)             :: CONVERGED   
      procedure(ComputeNonlinearStep1_f)                    :: ComputeTimeDerivative
      ! procedure(ComputeTimeDerivative_f)                    :: ComputeTimeDerivative
      !----------------------------------------------------------------------           
      integer              :: newtonit
      real(kind=RP)        :: norm, norm_old, rel_tol
      real(kind=RP)        :: linsolver_tol
      logical      , save  :: isfirst = .TRUE.
      real(kind=RP), save  :: norm1
      !----------------------------------------------------------------------
      
!
!     Initializations
!     ---------------
      
      if (isfirst) then
         norm = firstNorm  ! A value to define the initial linsolver_tol
         isfirst = .FALSE.
      else
         norm = norm1
      end if

      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      ConvRate = 1.0_RP
   
! !     -----------
!          call ComputeRHSSteps2(sem, t, dt, linsolver, ComputeTimeDerivative )               ! Computes b (RHS) and stores it into linsolver
         
!          write(*,*) " hello  DirectSolver SplitStepsTimeIntegrator.f90, after ComputeRHSSteps2 ============"
         
!          call linsolver%solve( nEqn=NCONS, nGradEqn=NGRAD, tol = linsolver_tol, maxiter=maxLinSolIter, time= t, dt=dt, &
!                               ComputeTimeDerivative = ComputeTimeDerivative, computeA = computeA)     


!          sem % mesh % storage % Q = linsolver % GetX()

!          ! call UpdateNewtonSol(sem, linsolver)                    ! Q_r+1 = Q_r + x
         
!          norm = linsolver%Getxnorm('infinity')
         
   
   end subroutine DirectSolve
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeRHSSteps2(sem, t, dt, linsolver, ComputeTimeDerivative )
      implicit none
      !----------------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      real(kind=RP),              intent(in)       :: t
      real(kind=RP),              intent(in)       :: dt
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      procedure(ComputeTimeDerivative_f)                   :: ComputeTimeDerivative
      !----------------------------------------------------------------
      real(kind=RP)                                :: value
      real(kind=RP)  :: RHS(NCONS*sem % NDOF)
      !----------------------------------------------------------------
      
      ! =========================== ZY modify ===============================
      call sem % mesh % storage % global2LocalQ
      call ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      ! call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      ! =========================== ZY modify ===============================
      
      RHS = SplitSteps_GetRHS(sem, linsolver, dt)
      
      call linsolver % SetRHS(RHS)
      
      call linsolver % AssemblyRHS     ! b must be assembled before using
   end subroutine ComputeRHSSteps2


   ! ////////////////////////////////////////////////////////////////////////////////////////////
   subroutine ComputeRHSSteps3Third(sem, t, dt, linsolver, ComputeTimeDerivative_Third,gamma, alpha0, alpha1, beta0, beta1 )
      implicit none
      !----------------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      real(kind=RP),              intent(in)       :: t
      real(kind=RP),              intent(in)       :: dt
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      procedure(ComputeNonlinearStep1_f)                   :: ComputeTimeDerivative_Third
      real(kind=RP),              intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         
      !----------------------------------------------------------------
      real(kind=RP)                                :: value
      real(kind=RP)  :: RHS(N_INS * sem % NDOF)
      !----------------------------------------------------------------
      
      ! =========================== ZY modify ===============================
      ! call sem % mesh % storage % global2LocalQ
      ! call ComputeTimeDerivative_Third( sem % mesh, sem % particles, dt, CTD_IGNORE_MODE)
      call ComputeTimeDerivative_Third( sem % mesh, sem % particles,t, dt, CTD_IGNORE_MODE, gamma, alpha0, alpha1, beta0, beta1)
      ! call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      ! =========================== ZY modify ===============================
      
      RHS = SplitSteps_GetRHS_Third(sem, linsolver, dt)
      
      call linsolver % SetRHS(RHS)
      
      call linsolver % AssemblyRHS     ! b must be assembled before using
   end subroutine ComputeRHSSteps3Third

   ! ////////////////////////////////////////////////////////////////////////////////////////////

   subroutine ComputeRHSSteps2Second(sem, t, dt, linsolver, ComputeTimeDerivative_Second, gamma, alpha0, alpha1, beta0, beta1)
      implicit none
      !----------------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      real(kind=RP),              intent(in)       :: t
      real(kind=RP),              intent(in)       :: dt
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      procedure(ComputeNonlinearStep1_f)                   :: ComputeTimeDerivative_Second
      real(kind=RP),              intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         
      !----------------------------------------------------------------
      real(kind=RP)                                :: value
      real(kind=RP)  :: RHS(1*sem % NDOF)
      !----------------------------------------------------------------
      
      ! =========================== ZY modify ===============================
      ! call sem % mesh % storage % global2LocalQ
      call ComputeTimeDerivative_Second( sem % mesh, sem % particles,t, dt, CTD_IGNORE_MODE, gamma, alpha0, alpha1, beta0, beta1)
      ! call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      ! =========================== ZY modify ===============================
      
      RHS = SplitSteps_GetRHS_second(sem, linsolver, dt)

      RHS(1) = 100000_RP
      
      call linsolver % SetRHS(RHS)
      
      call linsolver % AssemblyRHS     ! b must be assembled before using
   end subroutine ComputeRHSSteps2Second

   subroutine ComputeAdvectStep1(sem, t, dt, linsolver, ComputeNonlinearStep1, gamma, alpha0, alpha1, beta0, beta1)
      implicit none
      !----------------------------------------------------------------
         ! type(HexMesh), target    , intent(inout) :: mesh
      ! type(DGSem)                                  :: sem
      type(DGSem),             intent(inout)    :: sem
      real(kind=RP),              intent(in)       :: t
      real(kind=RP),              intent(in)       :: dt
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      procedure(ComputeNonlinearStep1_f)     :: ComputeNonlinearStep1

      real(kind=RP),              intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         

      !----------------------------------------------------------------
      real(kind=RP)                                :: value
      real(kind=RP)  :: RHS(NCONS*sem % NDOF)
      !----------------------------------------------------------------

      !     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, k
      

      call ComputeNonlinearStep1( sem % mesh, sem % particles, t,dt, mode = CTD_IGNORE_MODE,&
                                 gamma= gamma, alpha0= alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1)
      
   
   end subroutine ComputeAdvectStep1
   !  
!/////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine UpdateNewtonSol(sem, linsolver)

      type(DGSem),                     intent(inout)    :: sem
      class(GenericLinSolver_t),       intent(inout)    :: linsolver
      
      sem % mesh % storage % Q = sem % mesh % storage % Q  + linsolver % GetX()
      
   end subroutine UpdateNewtonSol
!
!////////////////////////////////////////////////////////////////////////////////////////////
!  TODO: Move from here....
   subroutine WriteEigenFiles(Mat,sem,FileName)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Writes files for performing eigenvalue analysis using TAUev
!        This only works for isotropic order meshes.........................TODO: Change that
!     -----------------------------------------------------------
!
      type(csrMat_t)    :: Mat      !< Jacobian matrix
      type(DGSem)       :: sem      !< DGSem class containing mesh
      character(len=*)  :: FileName !< ...
!     -----------------------------------------------------------
      integer           :: fd
!     -----------------------------------------------------------
      
      ! .frm file
      OPEN(newunit=fd, file=TRIM(FileName)//'.frm', action='write')
         write(fd,*)
         write(fd,*) SIZE(Mat % Values), SIZE(Mat % Rows)-1, 1, NCONS, 1
         write(fd,*) sem % mesh % elements(1) % Nxyz(1), SIZE(sem % mesh % elements)
      CLOSE (fd)
      
      ! .amg file
      call Mat % Visualize(TRIM(FileName)//'.amg',FirstRow=.FALSE.)
      
      ! .coo file
      call sem % mesh % WriteCoordFile(NCONS, TRIM(FileName)//'.coo')
      
      
   end subroutine WriteEigenFiles
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SplitSteps_SetOrder(order)
      implicit none
      !------------------------------------------------------
      integer, intent(in) :: order
      !------------------------------------------------------
      
      if (order > MAX_ORDER_CONS_DT) then
         write(STD_OUT,*) 'WARNING :: Maximum SplitSteps order for constant time-step is 5. Using 1 by default.'
         bdf_order = 1
      else
         bdf_order = order
      end if
      
   end subroutine SplitSteps_SetOrder
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SplitSteps_SetPreviousSolution(storage,NotANewStep)
      implicit none
      !------------------------------------------------------
      type(SolutionStorage_t), intent(inout) :: storage
      logical, optional :: NotANewStep
      !------------------------------------------------------
      integer :: i      ! Counter
      !------------------------------------------------------
      
      if (present(NotANewStep)) then
         if (.not. NotANewStep) StepsTaken = StepsTaken + 1
      else
         StepsTaken = StepsTaken + 1
      end if
      
      order = min(StepsTaken, bdf_order)
      
      call storage % SetGlobalPrevQ(storage % Q)
   end subroutine SplitSteps_SetPreviousSolution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function SplitSteps_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = 0
      ! Ashift = -SplitStepsCoeff(1,order)/dt
      write (*,*) "SplitSteps_MatrixShift_step1", Ashift , dt 
      
   end function SplitSteps_MatrixShift


   function SplitSteps_MatrixShift_step2(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      write (*,*) "SplitSteps_MatrixShift_step2", Ashift , dt 
      Ashift = 0
      ! Ashift = -SplitStepsCoeff(1,order)/dt
      
   end function SplitSteps_MatrixShift_step2

   function SplitSteps_MatrixShift_step3(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = 0.0_RP 
      write (*,*) "SplitSteps_MatrixShift_step3", Ashift , dt 
      ! write (*,*) "SplitSteps_MatrixShift_step3 bug", Ashift / (Ashift-Ashift)  
      ! Ashift = -SplitStepsCoeff(1,order)/dt
      
   end function SplitSteps_MatrixShift_step3
!

   !
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function SplitSteps_GetRHS(sem,linsolver, dt) result(RHS)
      implicit none
      !------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      ! type(HexMesh), target    , intent(inout) :: mesh
      ! type(SolutionStorage_t), intent(in) :: storage
      real(kind=RP)          , intent(in) :: dt
      real(kind=RP)                       :: RHS(size(sem % mesh % storage % Q))
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      !------------------------------------------------------
      integer :: k
      real(kind=RP) :: invdt
      !------------------------------------------------------
      
      invdt = 1._RP/dt


      write (*,*) "size(storage % Q) =========== ", size(sem % mesh % storage % Q)


      

      call linsolver % Jacobian % TimeDerivative_ComputeRHS( sem %  mesh, N_INS, dt, RHS, sem % NDOF )

      call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      call sem % mesh % storage % SolutionStorage_local2Global_Qdot1(sem % NDOF)

      write (*,*) "  sem %  mesh % storage % Qdot =======", sem %  mesh % storage % Qdot
      write (*,*) "  sem %  mesh % storage % slrDot1=======", sem %  mesh % storage % slrDot1


      
      RHS =  sem %  mesh % storage % Qdot
      ! call sem % mesh % storage % local2GlobalQdot(sem % NDOF)

      !  (this % p_sem, nEqn, time, this % PETScA, ComputeTimeDerivative, ComputeTimeDerivative)

      ! RHS =   storage % Qdot
      ! RHS = - storage % Qdot
      ! RHS = storage % Q * SplitStepsCoeff(1,order)*invdt - storage % Qdot
      
      ! do k=1, order
      !    RHS = RHS + SplitStepsCoeff(k+1,order) * storage % PrevQ(:,storage % prevSol_index(k)) * invdt
      ! end do
      
   end function SplitSteps_GetRHS
!
   function SplitSteps_GetRHS_Third(sem,linsolver, dt) result(RHS)
      implicit none
      !------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      ! type(HexMesh), target    , intent(inout) :: mesh
      ! type(SolutionStorage_t), intent(in) :: storage
      real(kind=RP)          , intent(in) :: dt
      real(kind=RP)                       :: RHS(size(sem % mesh % storage % vel_source))
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      !------------------------------------------------------
      integer :: k
      real(kind=RP) :: invdt
      !------------------------------------------------------
      
      invdt = 1._RP/dt


      ! write (*,*) "size(storage % Q) =========== ", size(sem % mesh % storage % Q)
      ! write (*,*) "  sem %  mesh % storage % vel_source =======3", sem %  mesh % storage % vel_source

      call linsolver % Jacobian % TimeDerivative_ComputeRHS_velocity_specific( sem %  mesh, N_INS, dt, RHS, sem % NDOF, &
                                                                      5                                    &
                                                                     )

      call sem % mesh % storage % SolutionStorage_local2Global_vel_source(sem % NDOF)
      ! call sem % mesh % storage % local2GlobalQdot(sem % NDOF)

      ! write (*,*) "  sem %  mesh % storage % vel_source =======3", sem %  mesh % storage % vel_source

      RHS =  sem %  mesh % storage % vel_source * 1.0_RP
      ! write (*,*) "RHS (=***********) =======", RHS
 
      
   end function SplitSteps_GetRHS_Third

   function SplitSteps_GetRHS_second(sem,linsolver, dt) result(RHS)
      implicit none
      !------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      ! type(HexMesh), target    , intent(inout) :: mesh
      ! type(SolutionStorage_t), intent(in) :: storage
      real(kind=RP)          , intent(in) :: dt
      real(kind=RP)                       :: RHS(size(sem % mesh % storage % pre_source))
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      !------------------------------------------------------
      integer :: k
      real(kind=RP) :: invdt
      !------------------------------------------------------
      
      invdt = 1._RP/dt


      ! write (*,*) "size(storage % Q) =========== ", size(sem % mesh % storage % Q)
      ! write (*,*) "  sem %  mesh % storage % pre_source =======1", sem %  mesh % storage % pre_source

      call linsolver % Jacobian % TimeDerivative_ComputeRHS_pressure_specific( sem %  mesh, 1, dt, RHS, sem % NDOF, &
                                                                      4                                    &
                                                                     )

      call sem % mesh % storage % SolutionStorage_local2Global_pre_source(sem % NDOF)

      call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      ! call sem % mesh % storage % SolutionStorage_local2Global_Qdot1(sem % NDOF)

      ! write (*,*) "  sem %  mesh % storage % pre_source =======2", sem %  mesh % storage % pre_source
      ! write (*,*) "  sem %  mesh % storage % slrDot1=======", sem %  mesh % storage % slrDot1


      RHS =  sem %  mesh % storage % pre_source
      ! RHS(1) = 0.0_RP
      ! write (*,*) "RHS (=***********) =======", RHS
 
      
   end function SplitSteps_GetRHS_second





   !
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SplitStepsInitialiseQ(mesh)
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                                        :: mesh
!-----Local-Variables-----------------------------------------------------
      integer :: i, id
!-------------------------------------------------------------------------
      
      do i = 1, bdf_order
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % prevQ(i) % Q = mesh % elements(id) % storage % Q
         end do ! id
!$omp end parallel do
      end do
      
   end subroutine SplitStepsInitialiseQ
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE SplitStepsTimeIntegrator