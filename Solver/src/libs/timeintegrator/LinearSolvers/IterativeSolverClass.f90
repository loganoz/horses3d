!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!      Class for solving a system with simple smoothers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE IterativeSolverClass
   use MatrixClass
   USE GenericLinSolverClass
   USE CSRMatrixClass
   USE SMConstants
   use DGSEMClass
   use TimeIntegratorDefinitions
   use NumericalJacobian
   use AnalyticalJacobian
   use PhysicsStorage
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC IterativeSolver_t, GenericLinSolver_t
   
   TYPE :: BlockPreco_t
      real(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: PLU        ! LU factorization of elemental preconditioner matrix
      integer      , dimension(:)  , allocatable :: LUpivots   ! LU pivots
   END TYPE BlockPreco_t
   
   TYPE, EXTENDS(GenericLinSolver_t) :: IterativeSolver_t
      TYPE(csrMat_t)                             :: A                                  ! Jacobian matrix
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      
      REAL(KIND=RP)                              :: rnorm                              ! L2 norm of residual
      REAL(KIND=RP)                              :: Ashift                             ! Shift that the Jacobian matrix currently(!) has
      
      CHARACTER(LEN=LINE_LENGTH)                 :: Smoother
      TYPE(BlockPreco_t), ALLOCATABLE            :: BlockPreco(:)
   CONTAINS
      !Subroutines:
      procedure                                  :: construct
      procedure                                  :: SetRHSValue
      procedure                                  :: SetRHSValues
      procedure                                  :: SetRHS => Iter_SetRHS
      procedure                                  :: solve
      procedure                                  :: GetCSRMatrix
      procedure                                  :: GetXValue
      procedure                                  :: destroy
      procedure                                  :: SetOperatorDt
      procedure                                  :: ReSetOperatorDt
      !Functions:
      procedure                                  :: GetX
      procedure                                  :: Getxnorm    !Get solution norm
      procedure                                  :: Getrnorm    !Get residual norm
      
      !! Internal procedures
      procedure                                  :: WeightedJacobiSmoother
      procedure                                  :: ComputeBlockPreco
   END TYPE IterativeSolver_t
   
!
!  ----------------
!  Module variables
!  ----------------
!
   REAL(KIND=RP)  :: timesolve ! Time
   REAL(KIND=RP)  :: dtsolve   ! dt   
   REAL(KIND=RP)  :: eps       ! Size of perturbation for matrix-free vector product
   
CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE construct(this,DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t) , INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      INTEGER                  , INTENT(IN)            :: globalDimPrb
      integer,       intent(in)                :: nEqn
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-----------------------------------------------------------
      INTEGER :: nelem      ! Number of elements
      INTEGER :: Nx,Ny,Nz   ! Polynomial orders for element
      INTEGER :: ndofelm    ! Number of degrees of freedom of element
      INTEGER :: k          ! Counter  
      INTEGER :: NCONS                                       
      !-----------------------------------------------------------
      
      call this % GenericLinSolver_t % construct(DimPrb,globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      
      IF (.NOT. PRESENT(sem)) error stop 'Fatal error: IterativeSolver needs sem.'
      
      MatrixShift => MatrixShiftFunc
      
      this % DimPrb = DimPrb
      this % Smoother = controlVariables % StringValueForKey("smoother",LINE_LENGTH)
      
      ALLOCATE(this % x   (DimPrb))
      ALLOCATE(this % b   (DimPrb))
      
      CALL this % A % construct (num_of_Rows = DimPrb, withMPI = .FALSE.)
      
      this % p_sem => sem
!
!     ------------------------------------------------
!     Allocate important variables for preconditioners
!     ------------------------------------------------
!
      SELECT CASE (this % Smoother)
         CASE('Block-Jacobi')
            nelem = SIZE(sem % mesh % elements)
            NCONS = SIZE(sem % mesh % elements(1) % storage % Q,4)
            ALLOCATE (this % BlockPreco(nelem))
            DO k = 1, nelem
               Nx = sem % mesh % elements(k) % Nxyz(1)
               Ny = sem % mesh % elements(k) % Nxyz(2)
               Nz = sem % mesh % elements(k) % Nxyz(3)
               ndofelm = NCONS*(Nx+1)*(Ny+1)*(Nz+1)
               allocate (this % BlockPreco(k) % PLU(ndofelm,ndofelm) )
               allocate (this % BlockPreco(k) % LUpivots   (ndofelm) )
            END DO
      END SELECT
   END SUBROUTINE construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetRHSValue(this, irow, value)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
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
      CLASS(IterativeSolver_t)   , INTENT(INOUT)     :: this
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
   subroutine Iter_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(IterativeSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % b = RHS
      
   end subroutine Iter_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE solve(this, nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt, ComputeA)
      IMPLICIT NONE
      CLASS(IterativeSolver_t), target, INTENT(INOUT) :: this
      integer,       intent(in)               :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      REAL(KIND=RP), OPTIONAL                 :: tol
      INTEGER      , OPTIONAL                 :: maxiter
      REAL(KIND=RP), OPTIONAL                 :: time
      REAL(KIND=RP), OPTIONAL                 :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-------------------------------------------------
      INTEGER                                 :: i
      !-------------------------------------------------
      
      IF (.NOT. PRESENT(time) .OR. .NOT. PRESENT(dt)) ERROR STOP 'time and dt needed for iterative solver'
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
            call this % SetOperatorDt(dt) ! Block preco computed inside
            ComputeA = .FALSE.
         end if
      else 
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
         call this % SetOperatorDt(dt) ! Block preco computed inside
      end if
      
      timesolve= time
      dtsolve  = dt
      
      ! Initialize x
      DO i=1, this % DimPrb
         this % x(i) = this % b(i) / this % A % Values(this%A%Diag(i))
      END DO
      
      
      SELECT CASE (this % Smoother)
         CASE('WeightedJacobi')
            CALL this % WeightedJacobiSmoother( this%A%Values(this%A%Diag), maxiter, tol, this % niter, ComputeTimeDerivative)
         CASE('Block-Jacobi')
            CALL BlockJacobiSmoother(this, maxiter, tol, this % niter, ComputeTimeDerivative)
      END SELECT
      
      
      IF (this % niter < maxiter) THEN
         this % CONVERGED = .TRUE.
      ELSE
         this % CONVERGED = .FALSE.
      END IF
      
   END SUBROUTINE solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetCSRMatrix(this,Acsr)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(IN)  :: this
      TYPE(csrMat_t)          , INTENT(OUT) :: Acsr 
      !-----------------------------------------------------------
      
      Acsr = this % A
   END SUBROUTINE GetCSRMatrix
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetXValue(this,irow,x_i)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------
      
      x_i = this % x(irow)
      
   END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function GetX(this) result(x)
      implicit none
      !-----------------------------------------------------------
      class(IterativeSolver_t), intent(inout)   :: this
      real(kind=RP)                             :: x(this % DimPrb)
      !-----------------------------------------------------------
      
      x = this % x
      
   end function GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destroy(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      CALL this % A % destruct()
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      
      IF(this % Smoother == 'Block-Jacobi') deallocate (this % BlockPreco)
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
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = MatrixShift(dt)
      CALL this % A % Shift(this % Ashift)
      
      IF(this % Smoother == 'Block-Jacobi') CALL this % ComputeBlockPreco
      
    END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ReSetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      CALL this % A % Shift(-this % Ashift)
      CALL this % A % Shift(shift)
      
      this % Ashift = shift
      
      IF(this % Smoother == 'Block-Jacobi') CALL this % ComputeBlockPreco
      
    END SUBROUTINE ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
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
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
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
      CLASS(IterativeSolver_t), INTENT(IN) :: this
      LOGICAL                              :: ComputeA
      
      ComputeA = .FALSE.
   END FUNCTION ComputeANextStep
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  Smoothers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   SUBROUTINE WeightedJacobiSmoother( this, Diag, SmoothIters, tol, niter, ComputeTimeDerivative)
      IMPLICIT NONE
      !--------------------------------------------
      CLASS(IterativeSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
      INTEGER                                         :: SmoothIters     !<  Number of smoothing operations
      REAL(KIND=RP)                                   :: Diag(:)         !   Matrix diagonal
      REAL(KIND=RP), OPTIONAL                         :: tol             !   Relative AND absolute tolerance of the method
      INTEGER                         , INTENT(OUT)   :: niter           !>  Number of iterations needed
      procedure(ComputeTimeDerivative_f)                      :: ComputeTimeDerivative
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
          r = CSR_MatVecMul(this%A,x) ! CSR matrix product
         
!$omp parallel do
         DO j=1,n
            r(j) = b(j) - r(j)
            x(j) = x(j) + w * r(j) / Diag(j)
         END DO
!$omp end parallel do
         
         IF (PRESENT(tol)) THEN
            rnorm = NORM2(r)       ! Saves relative tolerance (one iteration behind)
!~             print*, i, rnorm
!~             read(*,*)
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
!
   SUBROUTINE BlockJacobiSmoother(this, SmoothIters, tol, niter, ComputeTimeDerivative)
      USE DenseMatUtilities
      IMPLICIT NONE
      !--------------------------------------------
      CLASS(IterativeSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
      INTEGER                                         :: SmoothIters     !<  Number of smoothing operations
      REAL(KIND=RP), OPTIONAL                         :: tol             !   Relative AND absolute tolerance of the method
      INTEGER                         , INTENT(OUT)   :: niter           !>  Number of iterations needed
      procedure(ComputeTimeDerivative_f)                      :: ComputeTimeDerivative
      !--------------------------------------------
       INTEGER                                 :: n                ! System size
      REAL(KIND=RP)                           :: r(this % DimPrb) ! Residual
      REAL(KIND=RP)                           :: P_1r(this % DimPrb) ! Residual
      REAL(KIND=RP), POINTER                  :: x(:)             ! Solution
      REAL(KIND=RP), POINTER                  :: b(:)             ! Right-hand-side
      INTEGER                                 :: i,j              ! Counters
      INTEGER                                 :: idx1, idx2       ! Indexes of block
      
      REAL(KIND=RP)                           :: bnorm, rnorm, oldrnorm, ConvRate     ! Norm of b and r vectors
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
      
      oldrnorm = -1._RP
      ConvRate = 1._RP
      do i=1,SmoothIters
         
         r = CSR_MatVecMul(this%A,x) ! CSR matrix product
!$omp parallel do private(idx1,idx2)
         DO j=1,size (this % p_sem % mesh % elements)
            idx1 = this % A % BlockIdx(j)
            idx2 = this % A % BlockIdx(j+1) -1
            r(idx1:idx2) = b(idx1:idx2) - r(idx1:idx2)
            call SolveLU(ALU      = this%BlockPreco(j) % PLU, &
                         LUpivots = this%BlockPreco(j) % LUpivots, &
                         x = P_1r(idx1:idx2), &
                         b = r   (idx1:idx2))
            
            x(idx1:idx2) = x(idx1:idx2) + P_1r(idx1:idx2)
         END DO
!$omp end parallel do
         
         IF (PRESENT(tol)) THEN
            rnorm = NORM2(r)       ! Saves relative tolerance (one iteration behind)
!~             print*, '\x1b[1;34m', i, rnorm, rnorm/oldrnorm ,'\x1b[0m'
!~             read(*,*)
            IF (oldrnorm .NE. -1.0_RP) THEN
               ConvRate = ConvRate + (LOG10(oldrnorm/rnorm)-ConvRate)/i 
            ENDIF
            IF (rnorm < endtol .OR. ABS(rnorm/oldrnorm-1._RP) < 0.01_RP) THEN
               this % rnorm = rnorm
               oldrnorm     = rnorm
               EXIT
            END IF
            oldrnorm     = rnorm
         END IF
        
!~         IF (i==1) call xyplot(Sol(:,1))
!~         print*, x
!~         read(*,*)
      END DO
!~      print*, '\x1b[1;34mSmoother ConvRate:', ConvRate ,'\x1b[0m'
      this % rnorm = NORM2(r)
      niter=i
      
   END SUBROUTINE BlockJacobiSmoother
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ComputeBlockPreco(this)
      USE DenseMatUtilities
      IMPLICIT NONE
      !-------------------------------------------------------------
      CLASS(IterativeSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
      !-------------------------------------------------------------
      INTEGER :: k, N     ! Counter
      !-------------------------------------------------------------
      
!$omp parallel do schedule(runtime)
      DO k=1, size (this % p_sem % mesh % elements)
         N = this % A % BlockSizes(k)
         call ComputeLU (A        = this % A % GetBlock(k,N), &
                         ALU      = this % BlockPreco(k) % PLU, &
                         LUpivots = this % BlockPreco(k) % LUpivots)
      END DO
!$omp end parallel do
      
   END SUBROUTINE ComputeBlockPreco
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE IterativeSolverClass