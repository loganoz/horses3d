!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!      Class for solving a linear system obtained from a DGSEM discretization using p-Multigrid
!
!        This class is not finished yet!!!!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE MultigridSolverClass
   USE GenericLinSolverClass
   USE CSRMatrixClass
   USE SMConstants
   USE PhysicsStorage
   
   USE PolynomialInterpAndDerivsModule
   USE GaussQuadrature
   use DGSEMClass
   use TimeIntegratorDefinitions
   use MatrixClass
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC MultigridSolver_t, GenericLinSolver_t
   
   TYPE, EXTENDS(GenericLinSolver_t) :: MultigridSolver_t
      TYPE(csrMat_t)                             :: A                                  ! Jacobian matrix
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: F_Ur                               ! Qdot at the beginning of solving procedure
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: Ur                                 ! Q at the beginning of solving procedure
      REAL(KIND=RP)                              :: rnorm                              ! L2 norm of residual
      REAL(KIND=RP)                              :: Ashift                             ! Shift that the Jacobian matrix currently(!) has
      LOGICAL                                    :: AIsPrealloc                        ! Has A already been preallocated? (to be deprecated)
      
      ! Variables that are specially needed for Multigrid
      TYPE(MultigridSolver_t), POINTER           :: Child                 ! Next coarser multigrid solver
      TYPE(MultigridSolver_t), POINTER           :: Parent                ! Next finer multigrid solver
      INTEGER                                    :: MGlevel               ! Current Multigrid level
!~       INTEGER                    , ALLOCATABLE   :: OrderList(:,:,:)      ! List containing the polynomial orders of all elements in current solver (not needed now.. Maybe when using p-adaptation)
      TYPE(Interpolator_t)       , ALLOCATABLE   :: Restriction(:,:,:)    ! Restriction operators (element level)
      TYPE(Interpolator_t)       , ALLOCATABLE   :: Prolongation(:,:,:)   ! Prolongation operators (element level)
      
   CONTAINS
      !Subroutines:
      PROCEDURE                                  :: construct
      PROCEDURE                                  :: SetRHSValue
      PROCEDURE                                  :: SetRHSValues
      procedure                                  :: SetRHS           => MG_SetRHS
      PROCEDURE                                  :: solve
      PROCEDURE                                  :: GetXValue
      PROCEDURE                                  :: destroy
      PROCEDURE                                  :: SetOperatorDt
      PROCEDURE                                  :: ReSetOperatorDt
      !Functions:
      PROCEDURE                                  :: Getxnorm    !Get solution norm
      PROCEDURE                                  :: Getrnorm    !Get residual norm
      
      !! Internal procedures
      PROCEDURE                                  :: AxMult                  ! arueda: This is needed here for lower multigrid levels
      PROCEDURE                                  :: WeightedJacobiSmoother
   END TYPE MultigridSolver_t
   
!
!  ----------------
!  Module variables
!  ----------------
!
   REAL(KIND=RP)  :: timesolve ! Time
   REAL(KIND=RP)  :: dtsolve   ! dt   
   REAL(KIND=RP)  :: eps       ! Size of perturbation for matrix-free vector product
   
   ! Multigrid
   INTEGER        :: MGlevels  ! Total number of multigrid levels
   INTEGER        :: deltaN    ! 
   INTEGER        :: nelem     ! Number of elements (this is a p-multigrid implementation)
   
CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE construct(this,DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t) , INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      INTEGER                  , INTENT(IN)            :: globalDimPrb
      integer                  , intent(in)            :: nEqn
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-----------------------------------------------------------
      INTEGER                            :: lvl              ! Multigrid level
      INTEGER                            :: k
      INTEGER, DIMENSION(:), ALLOCATABLE :: Nx, Ny, Nz       ! Dimensions of fine mesh
      !-----------------------------------------------------------
      !Module variables: MGlevels, deltaN
      
      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      
      MatrixShift => MatrixShiftFunc
      
      IF (.NOT. PRESENT(sem)) error stop 'Fatal error: MultigridSolver needs sem.'
      IF (.NOT. PRESENT(controlVariables)) error stop 'Fatal error: MultigridSolver needs controlVariables.'
      
!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
      IF (.NOT. controlVariables % containsKey("multigrid levels")) THEN
         print*, 'Fatal error: "multigrid levels" keyword is needed by the multigrid solver'
         error stop
      END IF
      MGlevels  = controlVariables % IntegerValueForKey("multigrid levels")
      
      IF (controlVariables % containsKey("delta n")) THEN
         deltaN = controlVariables % IntegerValueForKey("delta n")
      ELSE
         deltaN = 1
      END IF
      
      this % p_sem => sem
      this % DimPrb = DimPrb
      
      nelem = SIZE(sem % mesh % elements)
      ALLOCATE(Nx(nelem),Ny(nelem),Nz(nelem))
      DO k=1, nelem
         Nx(k) = sem % mesh % elements(k) % Nxyz(1)
         Ny(k) = sem % mesh % elements(k) % Nxyz(2)
         Nz(k) = sem % mesh % elements(k) % Nxyz(3)
      END DO
!
!     --------------------------
!     Create linked solvers list
!     --------------------------
!
      CALL RecursiveConstructor(this, Nx, Ny, Nz, MGlevels, controlVariables)
      
      CALL this % A % construct (num_of_Rows = DimPrb)
      
      DEALLOCATE (Nx,Ny,Nz)
      
   END SUBROUTINE construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   RECURSIVE SUBROUTINE RecursiveConstructor(Solver, N1x, N1y, N1z, lvl, controlVariables)
      IMPLICIT NONE
      TYPE(MultigridSolver_t), TARGET  :: Solver
      INTEGER, DIMENSION(:)            :: N1x,N1y,N1z      !<  Order of approximation for every element in current solver
      INTEGER                          :: lvl              !<  Current multigrid level
      TYPE(FTValueDictionary)          :: controlVariables !< Control variables (for the construction of coarse sems
      !----------------------------------------------
      INTEGER                   :: DimPrb                 !   Dimension of problem for child solver
      INTEGER, DIMENSION(nelem) :: N2x,N2y,N2z            !   Order of approximation for every element in child solver
      INTEGER                   :: N1xMAX,N1yMAX,N1zMAX   !   Maximum polynomial orders for current (fine) grid
      INTEGER                   :: N2xMAX,N2yMAX,N2zMAX   !   Maximum polynomial orders for child (coarse) grid
      INTEGER                   :: k                      !   Counter
      LOGICAL                   :: success                ! Did the creation of sem succeed?
      TYPE(MultigridSolver_t) , POINTER :: Child_p          ! Pointer to Child
      !----------------------------------------------
      !
#if defined(NAVIERSTOKES)      
      
      Solver % MGlevel = lvl
      
      IF (lvl > 1) THEN
         ALLOCATE  (Solver % Child)
         Child_p => Solver % Child
         Solver % Child % Parent => Solver
!
!        -----------------------------------------------
!        Allocate restriction and prolongation operators
!        -----------------------------------------------
!         
         N1xMAX = MAXVAL(N1x)
         N1yMAX = MAXVAL(N1y)
         N1zMAX = MAXVAL(N1z)
         
         N2xMAX = N1xMAX - deltaN
         N2yMAX = N1yMAX - deltaN
         N2zMAX = N1zMAX - deltaN
         IF (N2xMAX < 0) N2xMAX = 0             ! TODO: Complete for Lobatto quadrature (max order = 1)
         IF (N2yMAX < 0) N2yMAX = 0
         IF (N2zMAX < 0) N2zMAX = 0
         
         ALLOCATE (Solver  % Restriction (0:N1xMAX,0:N1yMAX,0:N1zMAX))
         ALLOCATE (Child_p % Prolongation(0:N2xMAX,0:N2yMAX,0:N2zMAX))
         
!
!        ---------------------------------------------
!        Create restriction and prolongation operators
!        ---------------------------------------------
!
         DimPrb = 0
         DO k=1, nelem
            CALL CreateInterpolationOperators(Solver % Restriction, Child_p % Prolongation, &
                                              N1x(k),N1y(k),N1z(k),                         &
                                              N2x(k),N2y(k),N2z(k), DeltaN)    ! TODO: Add lobatto flag if required
            
            DimPrb = DimPrb + NCONS * (N2x(k) + 1) * (N2y(k) + 1) * (N2z(k) + 1)
         END DO
         Solver % Child % DimPrb = DimPrb
         
         ! Create DGSEM class for child
         ALLOCATE (Child_p % p_sem)
         CALL Child_p % p_sem % construct (controlVariables = controlVariables,              &
                                           Nx_ = N2x,    Ny_ = N2y,    Nz_ = N2z,            &
                                           success = success )
         IF (.NOT. success) error stop "Multigrid: Problem creating coarse solver."
         
         
         CALL RecursiveConstructor(Solver % Child, N2x, N2y, N2z, lvl - 1, controlVariables)
      END IF
      
      ALLOCATE(Solver % x   (Solver % DimPrb))
      ALLOCATE(Solver % b   (Solver % DimPrb))
      ALLOCATE(Solver % F_Ur(Solver % DimPrb))
      ALLOCATE(Solver % Ur  (Solver % DimPrb))
      
#endif
   END SUBROUTINE RecursiveConstructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetRHSValue(this, irow, value)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(IN)    :: value
      !-----------------------------------------------------------
      
      this % b (irow) = value
      
   END SUBROUTINE SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetRHSValues(this, nvalues, irow, values)
      IMPLICIT NONE
      CLASS(MultigridSolver_t)   , INTENT(INOUT)     :: this
      INTEGER                     , INTENT(IN)        :: nvalues
      INTEGER      , DIMENSION(1:), INTENT(IN)        :: irow
      REAL(KIND=RP), DIMENSION(1:), INTENT(IN)        :: values
      !------------------------------------------------------
      INTEGER                                        :: i
      
      DO i=1, nvalues
         IF (irow(i)<0) CYCLE
         this % b(irow(i)) = values(i)
      END DO
      
   END SUBROUTINE SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(MultigridSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % b = RHS
      
   end subroutine MG_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE solve(this, nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt, ComputeA)
      IMPLICIT NONE
      CLASS(MultigridSolver_t), target, INTENT(INOUT) :: this
      integer,       intent(in)               :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      REAL(KIND=RP), OPTIONAL                 :: tol
      INTEGER      , OPTIONAL                 :: maxiter
      REAL(KIND=RP), OPTIONAL                 :: time
      REAL(KIND=RP), OPTIONAL                 :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-------------------------------------------------
      INTEGER                                 :: niter
      INTEGER                                 :: i
!~      LOGICAL, SAVE :: isfirst = .TRUE.
      !-------------------------------------------------
      
      IF (.NOT. PRESENT(time) .OR. .NOT. PRESENT(dt)) STOP 'time and dt needed for Multigrid solver'
      
!
!     Compute Jacobian matrix if needed
!     -----------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
            call this % SetOperatorDt(dt) 
            ComputeA = .FALSE.
         end if
      else 
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
         call this % SetOperatorDt(dt) 
      end if
      
      timesolve= time
      dtsolve  = dt
      
!~      IF (isfirst) THEN
         this % F_Ur = this % p_sem % mesh % storage % Qdot
         this % Ur   = this % p_sem % mesh % storage % Q
!~         isfirst = .FALSE.
!~      END IF
      
      ! Initialize x
      DO i=1, this % DimPrb
         this % x(i) = this % b(i) / this % A % Values(this%A%Diag(i))
      END DO
      
      error stop 'incomplete solver'
      
      CALL this % WeightedJacobiSmoother( this%A%Values(this%A%Diag), maxiter, tol, this % niter)
      
      this % p_sem % mesh % storage % Q = this % Ur
      
!~      IF (this % niter < maxiter) THEN
         this % CONVERGED = .TRUE.
!~      ELSE
!~         this % CONVERGED = .FALSE.
!~      END IF
      
   END SUBROUTINE solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetXValue(this,irow,x_i)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------
      
      x_i = this % x(irow)
      
   END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destroy(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      CALL this % A % destruct()
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      
      this % AIsPrealloc = .FALSE.
      
    END SUBROUTINE destroy
    
    
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = MatrixShift(dt)
      call this % A % Shift(this % Ashift)
      
    END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ReSetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      
      call this % A % Shift(-this % Ashift)
      call this % A % Shift(shift)
      
      this % Ashift = shift
      
    END SUBROUTINE ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      SELECT CASE (TypeOfNorm)
         CASE ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         CASE ('l2')
            xnorm = NORM2(this % x)
         CASE DEFAULT
            error stop 'MKLPardisoSolverClass ERROR: Norm not implemented yet'
      END SELECT
   END FUNCTION Getxnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getrnorm(this) RESULT(rnorm)
      IMPLICIT NONE
!
!     ----------------------------------------
!     Currently implemented with infinity norm
!     ----------------------------------------
!
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------
      
      rnorm = this % rnorm
      
      
   END FUNCTION Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION ComputeANextStep(this) RESULT(ComputeA)
      IMPLICIT NONE
      CLASS(MultigridSolver_t), INTENT(IN) :: this
      LOGICAL                              :: ComputeA
      
      ComputeA = .FALSE.
   END FUNCTION ComputeANextStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  Internal procedures
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   FUNCTION AxMult(this,nEqn,x, ComputeTimeDerivative) RESULT(Ax)
      IMPLICIT NONE
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
      integer, intent(in)                     :: nEqn
      REAL(KIND=RP)                           :: x (:)
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      REAL(KIND=RP)                           :: Ax(size(x))
      !--------------------------------------------------
!~      REAL(KIND=RP)                           :: eps
      REAL(KIND=RP)                           :: F (size(x)), shift
!~      REAL(KIND=RP)                           :: buffer (size(x))
      !--------------------------------------------------
      shift = MatrixShift(dtsolve)
!~      eps = 1e-8_RP * (1._RP + NORM2(x))                           ! ~2e-5 2e-4
!~      eps = 1e-8_RP * (1._RP + NORM2(this % Ur))                   ! better: ~6e-7
      eps = SQRT(EPSILON(eps)) * (1._RP + NORM2(this % Ur))        !slightly better: ~4e-7 
!~      eps = SQRT(EPSILON(eps)) * (1._RP + MAXVAL(ABS(this % Ur)))  !slightly worse: ~1e-5 9e-6
!~      eps = SQRT(EPSILON(eps))                                     !worse:        : ~1e-4
      
!~      buffer = this % p_sem % mesh % storage % Q
      this % p_sem % mesh % storage % Q = this % Ur + x*eps
      call this % p_sem % mesh % storage % global2LocalQ
      CALL ComputeTimeDerivative(this % p_sem % mesh, this % p_sem % particles, timesolve, CTD_IGNORE_MODE)
      call this % p_sem % mesh % storage % local2GlobalQdot (this % p_sem % NDOF)
      
      F = this % p_sem % mesh % storage % Qdot
!~      this % p_sem % mesh % storage % Q = buffer
      Ax = ( F - this % F_Ur) / eps + shift * x                          !First order
      
   END FUNCTION AxMult
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  Smoothers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   SUBROUTINE WeightedJacobiSmoother( this, Diag, SmoothIters, tol, niter)
      IMPLICIT NONE
      !--------------------------------------------
      CLASS(MultigridSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
      INTEGER                                         :: SmoothIters     !<  Number of smoothing operations
      REAL(KIND=RP)                                   :: Diag(:)         !   Matrix diagonal
      REAL(KIND=RP), OPTIONAL                         :: tol             !   Relative AND absolute tolerance of the method
      INTEGER                         , INTENT(OUT)   :: niter           !>  Number of iterations needed
      !--------------------------------------------
      INTEGER                                 :: n                ! System size
      REAL(KIND=RP)                           :: r(this % DimPrb) ! Residual
      REAL(KIND=RP), POINTER                  :: x(:)             ! Solution
      REAL(KIND=RP), POINTER                  :: b(:)             ! Right-hand-side
      REAL(KIND=RP), PARAMETER                :: w = 2._RP/3._RP  ! Weight (optimal for FD laplacian... but DGSEM?)
      INTEGER                                 :: i,j              ! Counters
      
      REAL(KIND=RP)                           :: bnorm, rnorm     ! Norm of b and r vectors
      REAL(KIND=RP)                           :: endtol           ! Final tolerance that will be used to evaluate convergence 
!~       REAL(KIND=RP) :: res(size(x,1),1), LinChange
      !--------------------------------------------
      
      n =  this % DimPrb
      x => this % x
      b => this % b
      
      IF(PRESENT(tol)) THEN
         bnorm = NORM2(b)
         endtol = MAX(bnorm*tol,tol)  ! rtol and atol are taken as the same value
      END IF
      
!~      print*, 'bnorm = ', bnorm
!~      print*, '    iter      residual'
      
      DO i=1,SmoothIters
!~         r = this % AxMult(x)        ! Matrix free mult
         r = CSR_MatVecMul(this%A,x) ! CSR matrix product
         
!$omp parallel do
         DO j=1,n
            r(j) = b(j) - r(j)
            x(j) = x(j) + w * r(j) / Diag(j)
         END DO
!$omp end parallel do
         
         IF (PRESENT(tol)) THEN
            rnorm = NORM2(r)       ! Saves relative tolerance (one iteration behind)
            print*, i, rnorm
            read(*,*)
            IF (rnorm < endtol) THEN
               this % rnorm = rnorm
               EXIT
            END IF
         END IF
        
!~         IF (i==1) call xyplot(Sol(:,1))
!~         print*, x
!~         read(*,*)
      END DO
      
      this % rnorm = NORM2(r)
      niter=i
      
   END SUBROUTINE WeightedJacobiSmoother
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  Routines for interpolation procedures (will probably be moved to another module)
!  
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE CreateInterpolationOperators(Restriction,Prolongation,N1x,N1y,N1z,N2x,N2y,N2z,DeltaN,Lobatto)
      IMPLICIT NONE
!
!     ------------------------------------------------------------------------
!     Creates the restriction and prolongation operators for a certain element
!     for multigrid. Takes into account order anisotropy, but the coarse grid 
!     is constructed by reducing the polynomial order uniformly.
!     ------------------------------------------------------------------------
!
      !-----------------------------------------------------
      TYPE(Interpolator_t), TARGET  :: Restriction (0:,0:,0:)  !>  Restriction operator
      TYPE(Interpolator_t), TARGET  :: Prolongation(0:,0:,0:)  !>  Prolongation operator
      INTEGER                       :: N1x, N1y, N1z           !<  Fine grid order(anisotropic) of the element
      INTEGER                       :: N2x, N2y, N2z           !>  Coarse grid order(anisotropic) of the element
      INTEGER                       :: DeltaN                  !<  Interval of reduction of polynomial order for coarser level
      LOGICAL, OPTIONAL             :: Lobatto                 !<  Is the quadrature a Legendre-Gauss-Lobatto representation?
      !-----------------------------------------------------
      TYPE(Interpolator_t), POINTER :: rest                    ! Pointer to constructed restriction interpolator
      TYPE(Interpolator_t), POINTER :: prol                    ! Pointer to constructed prolongation interpolator
      LOGICAL                       :: LGL = .FALSE.           ! Is the quadrature a Legendre-Gauss-Lobatto representation? (false is default)
      INTEGER                       :: i,j,k,l,m,n             ! Index counters
      INTEGER                       :: s,r                     ! Row/column counters for operators
      REAL(KIND=RP), ALLOCATABLE    :: x1 (:), y1 (:), z1 (:)  ! Position of quadrature points on mesh 1
      REAL(KIND=RP), ALLOCATABLE    :: w1x(:), w1y(:), w1z(:)  ! Weights for quadrature points on mesh 1
      REAL(KIND=RP), ALLOCATABLE    :: x2 (:), y2 (:), z2 (:)  ! Position of quadrature points on mesh 2
      REAL(KIND=RP), ALLOCATABLE    :: w2x(:), w2y(:), w2z(:)  ! Weights for quadrature points on mesh 2
      !-----------------------------------------------------
      
      IF (PRESENT(Lobatto)) THEN
         if (Lobatto) LGL = .TRUE.
      END IF
!
!     --------------------------------------
!     Compute order of coarse representation
!     --------------------------------------
!
      N2x = N1x - DeltaN
      N2y = N1y - DeltaN
      N2z = N1z - DeltaN
      
      ! The order must be greater or equal to 0 (Legendre-Gauss quadrature) or 1 (Legendre-Gauss-Lobatto)
      IF (LGL) THEN
         IF (N2x < 1) N2x = 1
         IF (N2y < 1) N2y = 1
         IF (N2z < 1) N2z = 1
      ELSE
         IF (N2x < 0) N2x = 0
         IF (N2y < 0) N2y = 0
         IF (N2z < 0) N2z = 0
      END IF
      
      ! Return if the operators were already created
      IF (Restriction(N1x,N1y,N1z) % Created) RETURN
      
      rest => Restriction (N1x,N1y,N1z)
      prol => Prolongation(N2x,N2y,N2z)
!
!     ----------------------------
!     Allocate important variables
!     ----------------------------
!
      !Nodes and weights
      ALLOCATE(x1 (0:N1x), y1 (0:N1y), z1 (0:N1z), &
               w1x(0:N1x), w1y(0:N1y), w1z(0:N1z), &
               x2 (0:N2x), y2 (0:N2y), z2 (0:N2z), &
               w2x(0:N2x), w2y(0:N2y), w2z(0:N2z))
!
!     ------------------------------------------
!     Obtain the quadrature nodes on (1) and (2)
!     ------------------------------------------
!
      IF (LGL) THEN
         CALL LegendreLobattoNodesAndWeights(N1x, x1, w1x)
         CALL LegendreLobattoNodesAndWeights(N1y, y1, w1y)
         CALL LegendreLobattoNodesAndWeights(N1z, z1, w1z)
         CALL LegendreLobattoNodesAndWeights(N2x, x2, w2x)
         CALL LegendreLobattoNodesAndWeights(N2y, y2, w2y)
         CALL LegendreLobattoNodesAndWeights(N2z, z2, w2z)
      ELSE
         CALL GaussLegendreNodesAndWeights(N1x, x1, w1x)
         CALL GaussLegendreNodesAndWeights(N1y, y1, w1y)
         CALL GaussLegendreNodesAndWeights(N1z, z1, w1z)
         CALL GaussLegendreNodesAndWeights(N2x, x2, w2x)
         CALL GaussLegendreNodesAndWeights(N2y, y2, w2y)
         CALL GaussLegendreNodesAndWeights(N2z, z2, w2z)
      END IF
!
!     -----------------------------
!     Fill the restriction operator
!     -----------------------------
!
      CALL Create3DInterpolationMatrix(rest % Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2)
!
!     ------------------------------
!     Fill the prolongation operator
!     ------------------------------
!
      CALL Create3DInterpolationMatrix(prol % Mat,N2x,N2y,N2z,N1x,N1y,N1z,x2,y2,z2,x1,y1,z1)
      
      ! All done
      rest % Created = .TRUE.
      prol % Created = .TRUE.
      
   END SUBROUTINE CreateInterpolationOperators
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!!
!!    This function interpolates a vector U1 of the solution to the lower level (can be used for restriction or prolongation)
!!
!!    CALL: U2 = this % InterpolateSol(Prolongation(Nx(iEL),Ny(iEL),Nz(iEL)),U1)
!!          
   FUNCTION InterpolateSol(this,Interp,U1) RESULT(U2)
      IMPLICIT NONE
      !--------------------------------------------------------------
      CLASS(MultigridSolver_t), TARGET, INTENT(INOUT) :: this  !<  Iterative solver class
      REAL(KIND=RP)                                   :: Interp(:,:)
      REAL(KIND=RP)           , TARGET, INTENT(IN)    :: U1(:) ! Vector to be projected
      REAL(KIND=RP)           , TARGET                :: U2(SIZE(Interp,1)) ! Projected vector
      !--------------------------------------------------------------
      REAL(KIND=RP), POINTER :: U1_p(:)               ! U1 pointer
      REAL(KIND=RP), POINTER :: U2_p(:)               ! U2 pointer
      INTEGER      , POINTER :: N1x(:),N1y(:),N1z(:)  ! Pointers to element orders of (1)
      INTEGER      , POINTER :: N2x(:),N2y(:),N2z(:)  ! Pointers to element orders of (2)
      INTEGER                :: iEQ       ! Equation counter
      INTEGER                :: Idx1      ! First index of the solution of this element ( - 1)
      !--------------------------------------------------------------
#if defined(NAVIERSTOKES)      
      Idx1 = 0
      
      DO iEQ = 1, NCONS
!~             U1_p => U1(Idx1+iEq::NCONS)
!~             U2_p => U2()
         
         U2_p = MATMUL(Interp,U1_p)
         
      END DO
!
!     ------------
!     Clean memory
!     ------------
!
      NULLIFY(U1_p,U2_p,N1x,N1y,N1z,N2x,N2y,N2z)
#endif      
   END FUNCTION InterpolateSol
   
   !SUBROUTINE InterpolateJac
END MODULE MultigridSolverClass