!
!////////////////////////////////////////////////////////////////////////
!
!      Module for integrating in time using Rosenbrock type implicit Runge-Kutta schemes
!
!////////////////////////////////////////////////////////////////////////
#include "Includes.h"
module RosenbrockTimeIntegrator
   use DGSEMClass
   use SMConstants
   use LinearSolverClass
   use PhysicsStorage
   implicit none
   
   private
   public RosenbrockIntegrator_t
   
!
!  ***************************
!  Rosenbrock integrator class
!  ***************************
   type RosenbrockIntegrator_t
      
      class(GenericLinSolver_t), allocatable :: linsolver      ! Linear solver
      integer                                :: NumStages      ! Number of stages of the Rosenbrock scheme
      
      real(kind=RP)            , allocatable :: Y(:,:)         ! Intermediate solutions
      logical                                :: JacByConv     ! .TRUE. if the Jacobian must be computed only when the convergence is bad
      logical                                :: TimeAccurate  ! .TRUE. if this is a time-accurate simulation
      
      contains
         procedure :: construct
         procedure :: destruct
         procedure :: TakeStep
         procedure :: SetupCoefficients
         procedure :: ComputeRHS
   
   end type RosenbrockIntegrator_t
   
!
!  ***********************
!  Rosenbrock coefficients
!  ***********************
   real(kind=RP)              :: Ros_gamma
   real(kind=RP), allocatable :: Ros_a(:,:)
   real(kind=RP), allocatable :: Ros_c(:,:)
   real(kind=RP), allocatable :: Ros_m(:)
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   subroutine construct(this,controlVariables,sem)
      use FTValueDictionaryClass
      implicit none
      !---------------------------------------------------------------
      class(RosenbrockIntegrator_t)       :: this
      type(FTValueDictionary), intent(in) :: controlVariables
      type(DGSem)            , intent(in) :: sem
      !---------------------------------------------------------------
      integer :: DimPrb, globalDimPrb
      !---------------------------------------------------------------
      
!
!     Setup Rosenbrock coefficients
!     -----------------------------
      call this % SetupCoefficients ( trim(controlVariables % StringValueForKey("rosenbrock scheme",LINE_LENGTH)) )
      
      allocate ( this % Y (sem % NDOF * NCONS,this % NumStages) )
!
!     Setup linear solver
!     -------------------
      DimPrb = sem % NDOF * NCONS
      globalDimPrb = sem % totalNDOF * NCONS
      
      select case ( trim(controlVariables % StringValueForKey("linear solver",LINE_LENGTH)) )
         case('petsc')
            allocate (PetscKspLinearSolver_t :: this % linsolver)
         case('pardiso')
            allocate (MKLPardisoSolver_t     :: this % linsolver)
         case('matrix-free smooth')
            allocate (MatFreeSmooth_t        :: this % linsolver)
         case('matrix-free gmres')
            allocate (MatFreeGMRES_t         :: this % linsolver)
         case('smooth')
            allocate (IterativeSolver_t      :: this % linsolver)
         case('multigrid')
            allocate (LinearMultigridSolver_t      :: this % linsolver)
         case default
            print*, "Keyword 'linear solver' missing... Using PETSc as default"
            allocate (PetscKspLinearSolver_t :: this % linsolver)
      end select
      
      call this % linsolver % construct (DimPrb,globalDimPrb,NCONS,controlVariables,sem,Rosenbrock_MatrixShift)
      
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !--------------------------------------------------------
      class(RosenbrockIntegrator_t) :: this
      !--------------------------------------------------------
      
      call this % linsolver % destroy
      deallocate (this % linsolver)
      
      safedeallocate( Ros_a )
      safedeallocate( Ros_c )
      safedeallocate( Ros_m )
      safedeallocate( this % Y )
      
   end subroutine destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine TakeStep(this,sem,t,dt,ComputeTimeDerivative)
      implicit none
      !--------------------------------------------------------
      class(RosenbrockIntegrator_t) :: this
      type(DGSem),    intent(inout) :: sem                  !<>DGSem class with solution storage 
      real(kind=RP),  intent(in)    :: t                    !< Time at the beginning of time step
      real(kind=RP),  intent(in)    :: dt                   !< Initial (outer) time step (the subroutine can use a smaller one depending on convergence)
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      !--------------------------------------------------------
      integer :: stage     ! Current stage
      logical :: computeA  ! Must the linear solver compute the Jacobian?
      !--------------------------------------------------------
      
      computeA = .TRUE. ! Check this!
      this % Y = 0._RP
      
      write(STD_OUT, "(A8,x,A18,x,A10,x,A9)") "Stage", "LinSolverErr","iterations","CONVERGED"
      
      do stage = 1, this % NumStages
         
         call this % ComputeRHS(sem, t, dt, this % linsolver, ComputeTimeDerivative, stage)
         
         CALL this % linsolver % solve ( nEqn=NCONS, nGradEqn=NGRAD, tol = 1e-6_RP, maxiter=500, time= t, dt=dt, &
                                          ComputeTimeDerivative = ComputeTimeDerivative, computeA = computeA)        ! Solve (J-I/dt)·x = (Q_r- U_n)/dt - Qdot_r
         
         this % Y (:,stage) = this % linsolver % GetX()
         
         write(STD_OUT, "(I8,x,E18.3,x,I10,x,L9)") stage, this % linsolver % Getrnorm(), this % linsolver % niter, this % linsolver % CONVERGED
      end do
      
!     Dump solution to sem
!     --------------------
      associate (Q => sem % mesh % storage % Q)
      do stage = 1, this % NumStages
         Q = Q + Ros_m(stage) * this % Y(:,stage)
      end do
      end associate
      
   end subroutine TakeStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeRHS(this, sem, t, dt, linsolver, ComputeTimeDerivative, stage)
      implicit none
      !--------------------------------------------------------
      class(RosenbrockIntegrator_t) :: this
      type(DGSem),    intent(inout) :: sem
      real(kind=RP),  intent(in)    :: t                       !<  Time
      real(kind=RP),  intent(in)    :: dt                      !<  Time-step
      class(GenericLinSolver_t)     :: linsolver               !<> Linear solver to load solution to
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      integer,        intent(in)    :: stage                   !<  Current stage
      !--------------------------------------------------------
      real(kind=RP) :: Qn (sem % NDOF * NCONS) ! Buffer to store the previous solution
      real(kind=RP) :: RHS(sem % NDOF * NCONS) ! Right-hand side for this stage
      integer       :: j               ! Counter
      !--------------------------------------------------------
      
      associate (Q    => sem % mesh % storage % Q, &
                 Qdot => sem % mesh % storage % Qdot )
      
!     Initializations
!     ---------------
      Qn    = Q
      RHS   = 0._RP
      
!     Compute RHS
!     -----------
      
      do j = 1, stage - 1
         Q   = Q   + Ros_a(j,stage) * this % Y(:,j)
         RHS = RHS - Ros_c(j,stage) * this % Y(:,j)
      end do
      
      call sem % mesh % storage % global2LocalQ
      call ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      
      RHS = RHS/dt - Qdot
      
!     Restore previous solution
!     -------------------------      
      Q = Qn
      end associate
      
!     Load RHS into solver
!     --------------------
      CALL linsolver % SetRHS(RHS)
      
   end subroutine ComputeRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetUpCoefficients(this,RosType)
      implicit none
      !--------------------------------------------------------
      class(RosenbrockIntegrator_t) :: this
      character(len=*), intent(in)  :: RosType
      !--------------------------------------------------------
      
      select case ( trim(RosType) )
!        *************************************
!        Six-stage sixth order accurate scheme     !> Bassi et al. Linearly implicit Rosenbrock-type Runge–Kutta schemes applied to the Discontinuous Galerkin solution of compressible and incompressible unsteady flows
!        *************************************
         case('RO6-6')
            this % NumStages = 6
            
            Ros_gamma = 3.341423670680504e-1_RP
            
            allocate (Ros_a(1:5,2:6))
            Ros_a = 0._RP
            
            Ros_a(1,2) = 2.0_RP
            Ros_a(1:2,3) = (/  1.751493065942685_RP    , -1.454290536332865e-01_RP /)
            Ros_a(1:3,4) = (/ -1.847093912231436_RP    , -2.513756792158473_RP    ,  1.874707432337999_RP /)
            Ros_a(1:4,5) = (/  1.059634783677141e1_RP  ,  1.974951525952609_RP    , -1.905211286263863_RP    , -3.575118228830491_RP /)
            Ros_a(1:5,6) = (/  2.417642067883312_RP    ,  3.050984437044573e-01_RP, -2.346208879122501e-01_RP, -1.327038464607418e-01_RP , 3.912922779645768e-02_RP /)
            
            allocate (Ros_c(1:5,2:6))
            Ros_c = 0._RP
            
            Ros_c(1,2) = -1.745029492512995e1_RP 
            Ros_c(1:2,3) = (/  -1.202359936227844e1_RP, 1.315910110742745_RP /)
            Ros_c(1:3,4) = (/  2.311230597159272e1_RP,  1.297893129565445e1_RP, -8.445374594562038_RP /)
            Ros_c(1:4,5) = (/  -3.147228891330713_RP,  -1.761332622909965_RP  ,  6.115295934038585_RP, 1.499319950457112e1_RP /)
            Ros_c(1:5,6) = (/  -2.015840911262880e1_RP,-1.603923799800133_RP  ,  1.155870096920252_RP, 6.304639815292044e-01_RP, -1.602510215637174e-01_RP /)
            
            allocate ( Ros_m(6) )
            Ros_m = (/ 3.399347452674165e1_RP, -2.091829882847333e1_RP, -1.375688477471081e1_RP, -1.113925929930077e1_RP, 2.873406527609468_RP, 3.876609945620840e1_RP /)
            
!        ***********************
!        Not implemented scheme!
!        ***********************
         case default
            error stop ':: Requested Rosenbrock scheme is not yet implemented'
      end select
   
   end subroutine SetUpCoefficients
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Rosenbrock_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = -1._RP / (Ros_gamma * dt)
      
   end function Rosenbrock_MatrixShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module RosenbrockTimeIntegrator