!
!      Class for solving a system with a simple BlockJacobi smoother and matrix free operations
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE MatrixFreeSmootherClass
   use MatrixClass
   USE GenericLinSolverClass
   USE CSRMatrixClass
   USE SMConstants
   use DGSEMClass
   use TimeIntegratorDefinitions
   use PhysicsStorage
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC MatFreeSmooth_t, GenericLinSolver_t
   
   TYPE :: BlockPreco_t
      real(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: PLU        ! LU factorization of elemental preconditioner matrix
      integer      , dimension(:)  , allocatable :: LUpivots   ! LU pivots
   END TYPE BlockPreco_t
   
   TYPE, EXTENDS(GenericLinSolver_t) :: MatFreeSmooth_t
      TYPE(DenseBlockDiagMatrix_t)               :: A                                  ! Jacobian matrix
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: F_Ur                               ! Qdot at the beginning of solving procedure
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: Ur                                 ! Q at the beginning of solving procedure
      REAL(KIND=RP)                              :: rnorm                              ! L2 norm of residual
      REAL(KIND=RP)                              :: Ashift                             ! Shift that the Jacobian matrix currently(!) has
      
      CHARACTER(LEN=LINE_LENGTH)                 :: Smoother
      TYPE(BlockPreco_t), ALLOCATABLE            :: BlockPreco(:)
   CONTAINS
      !Subroutines:
      PROCEDURE                                  :: construct
      PROCEDURE                                  :: SetRHSValue
      PROCEDURE                                  :: SetRHSValues
      procedure                                  :: SetRHS => Smooth_SetRHS
      PROCEDURE                                  :: solve
      PROCEDURE                                  :: GetXValue
      PROCEDURE                                  :: destroy
      PROCEDURE                                  :: SetOperatorDt
      PROCEDURE                                  :: ReSetOperatorDt
      !Functions:
      procedure                                  :: GetX
      PROCEDURE                                  :: Getxnorm    !Get solution norm
      PROCEDURE                                  :: Getrnorm    !Get residual norm
      
      !! Internal procedures
      PROCEDURE                                  :: AxMult
      PROCEDURE                                  :: ComputeBlockPreco
      procedure                                  :: SetInitialGuess
      
      PROCEDURE                                  :: p_F
   END TYPE MatFreeSmooth_t
   
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
      CLASS(MatFreeSmooth_t) , INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      INTEGER                  , INTENT(IN)            :: globalDimPrb
      integer                , intent(in)            :: nEqn
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-----------------------------------------------------------
      INTEGER :: nelem      ! Number of elements
      INTEGER :: Nx,Ny,Nz   ! Polynomial orders for element
      INTEGER :: ndofelm    ! Number of degrees of freedom of element
      INTEGER :: k          ! Counter                     
      !-----------------------------------------------------------
      
      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      
      IF (.NOT. PRESENT(sem)) error stop 'Fatal error: IterativeSolver needs sem.'
      
      MatrixShift => MatrixShiftFunc
      
      this % DimPrb = DimPrb
      this % Smoother = controlVariables % StringValueForKey("smoother",LINE_LENGTH)
      
      ALLOCATE(this % x   (DimPrb))
      ALLOCATE(this % b   (DimPrb))
      ALLOCATE(this % F_Ur(DimPrb))
      ALLOCATE(this % Ur  (DimPrb))
      
      this % p_sem => sem
      nelem = SIZE(sem % mesh % elements)
      
      call this % A % construct(num_of_Blocks = nelem)
!
!     ------------------------------------------------
!     Allocate important variables for preconditioners
!     ------------------------------------------------
!
      SELECT CASE (this % Smoother)
         CASE('Block-Jacobi')
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
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
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
      CLASS(MatFreeSmooth_t)   , INTENT(INOUT)     :: this
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
   subroutine Smooth_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(MatFreeSmooth_t), intent(inout) :: this
      real(kind=RP)         , intent(in)    :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % b = RHS
      
   end subroutine Smooth_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE solve(this, nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt, ComputeA)
      use DenseMatUtilities
      IMPLICIT NONE
      CLASS(MatFreeSmooth_t), target, INTENT(INOUT) :: this
      integer, intent(in)                     :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      REAL(KIND=RP), OPTIONAL                 :: tol
      INTEGER      , OPTIONAL                 :: maxiter
      REAL(KIND=RP), OPTIONAL                 :: time
      REAL(KIND=RP), OPTIONAL                 :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-------------------------------------------------
      INTEGER                                 :: i, k
      logical :: TolPresent
      !-------------------------------------------------
      
      IF (.NOT. PRESENT(time) .OR. .NOT. PRESENT(dt)) error stop 'time and dt needed for iterative solver'
      TolPresent = present(tol)
!
!     Compute Jacobian matrix if needed
!     -----------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
            call this % SetOperatorDt(dt) ! the matrix is factorized inside
            
            ComputeA = .FALSE.
         end if
      else 
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
         call this % SetOperatorDt(dt) ! the matrix is factorized inside
      end if
      
      timesolve= time
      dtsolve  = dt
      
!~      IF (isfirst) THEN
         this % F_Ur = this % p_sem % mesh % storage % Qdot
         this % Ur   = this % p_sem % mesh % storage % Q
!~         isfirst = .FALSE.
!~      END IF
      
      ! Initialize x
      call this % SetInitialGuess
      
      SELECT CASE (this % Smoother)
         CASE('Block-Jacobi')
            CALL BlockJacobiSmoother(this, maxiter, this % niter, ComputeTimeDerivative, TolPresent, tol)
      END SELECT
      
      this % p_sem % mesh % storage % Q = this % Ur
      
      IF (this % niter <= maxiter) THEN
         this % CONVERGED = .TRUE.
      ELSE
         this % CONVERGED = .FALSE.
      END IF
      
   END SUBROUTINE solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetInitialGuess(this)
      implicit none
      !-----------------------------------------------------------
      CLASS(MatFreeSmooth_t),target, INTENT(INOUT) :: this
      !-----------------------------------------------------------
      integer :: k,i, firstIdx, lastIdx
      real(kind=RP), pointer :: x_p(:), b_p(:), Mat_p(:,:)
      !-----------------------------------------------------------
      
!$omp parallel do private(i,x_p,b_p,Mat_p,firstIdx,lastIdx) schedule(runtime)
      do k=1, this % A % num_of_Blocks
         firstIdx = this % A % BlockIdx(k)
         lastIdx  = this % A % BlockIdx(k+1) - 1
         Mat_p => this % A % Blocks(k) % Matrix
         x_p => this % x(firstIdx:lastIdx)
         b_p => this % b(firstIdx:lastIdx)
         do i=1, size(this % A % Blocks(k) % Matrix,1)
            x_p(i) = b_p(i) / Mat_p(i,i)
         end do
      end do
!$omp end parallel do
   end subroutine SetInitialGuess
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetXValue(this,irow,x_i)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
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
      class(MatFreeSmooth_t), intent(inout)  :: this
      real(kind=RP)                          :: x(this % DimPrb)
      !-----------------------------------------------------------
      
      x = this % x
      
   end function GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destroy(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      CALL this % A % destruct()
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      
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
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = MatrixShift(dt)
      CALL this % A % Shift( this % Ashift )
      
      IF(this % Smoother == 'Block-Jacobi') CALL this % ComputeBlockPreco
      
    END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ReSetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      
      CALL this % A % Shift( -this % Ashift )
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
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      SELECT CASE (TypeOfNorm)
         CASE ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         CASE ('l2')
            xnorm = NORM2(this % x)
         CASE DEFAULT
            error stop 'MatFreeSmoothClass ERROR: Norm not implemented yet'
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
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
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
      CLASS(MatFreeSmooth_t), INTENT(IN) :: this
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

   FUNCTION AxMult(this,x, computeTimeDerivative) RESULT(Ax)
      IMPLICIT NONE
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                           :: x (:)
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      REAL(KIND=RP)                           :: Ax(size(x))
      !--------------------------------------------------
      real(kind=RP)                           :: shift
!~      REAL(KIND=RP)                           :: eps
      REAL(KIND=RP)                           :: F (size(x))
      
!~       REAL(KIND=RP)                           :: xxx (size(x)) !x vector... But normalized!!
!~      REAL(KIND=RP)                           :: buffer (size(x))
      !--------------------------------------------------
      
      shift = MatrixShift(dtsolve)
      
!~      eps = 1e-8_RP * (1._RP + NORM2(x))                           ! ~2e-5 2e-4
!~      eps = 1e-8_RP * (1._RP + NORM2(this % Ur))                   ! better: ~6e-7
!~       eps = SQRT(EPSILON(eps)) * (1._RP + NORM2(this % Ur))        !slightly better: ~4e-7 
!~       eps = SQRT(EPSILON(eps)) * (NORM2(this % Ur))
!~      eps = SQRT(EPSILON(eps)) * (1._RP + MAXVAL(ABS(this % Ur)))  !slightly worse: ~1e-5 9e-6
!~      eps = SQRT(EPSILON(eps))                                     !worse:        : ~1e-4
      
!~       eps = SQRT(EPSILON(eps)) * NORM2(this % Ur) / NORM2(x) ! hillewaert2013 ! Best performance... but eps too big for small x?
      
!~      eps = SQRT(EPSILON(eps)) * NORM2(this % Ur) / NORM2(x)**2 ! This doesn't work at all
!~      eps = SQRT(EPSILON(eps)) * SIGN(MAX(DOT_PRODUCT(this % Ur,x),MAXVAL(ABS(x))),DOT_PRODUCT(this % Ur,x)) / (NORM2(x)) ! Saad with typical value u~1
!~      eps = SQRT(EPSILON(eps)) * SIGN(NORM2(this % Ur),DOT_PRODUCT(this % Ur,x)) / NORM2(x) ! hillewaert2003 using different sign (same behavior)
!~       eps = SQRT(EPSILON(eps)) * (1._RP + NORM2(this % Ur)) / NORM2(x) !Combining hillawaert with Sipp
!~       eps = SQRT(EPSILON(eps)) * 1._RP + NORM2(this % Ur) / (NORM2(x))
!~       eps = SQRT(EPSILON(eps) * (1._RP + NORM2(this % Ur))) / (NORM2(x)) ! NISTOL Package
      
!~      eps = SQRT(EPSILON(eps) * (1._RP + NORM2(this % Ur))) * DOT_PRODUCT(this % Ur,x) / (NORM2(x)) ! NISTOL Package modified    !! This works down to 10⁻¹⁰, but slowly (linear residual does not always go as low as desired).. At the end short tie-steps are needed
      
       eps = SQRT(EPSILON(eps)) * (NORM2(this % Ur)*DOT_PRODUCT(this % Ur,x)) / NORM2(x) ! My recipe.. goes lower but slower
      
!~      buffer = this % p_sem % mesh % storage % Q

!~       xxx = x / NORM2(x)

!~       this % p_sem % mesh % storage % Q = this % Ur + x*eps
!~       CALL ComputeTimeDerivative(this % p_sem,timesolve)
!~       F =this % p_sem % mesh % storage % Qdot
!~      this % p_sem % mesh % storage % Q = buffer
      Ax = ( this % p_F(this % Ur + x * eps, computeTimeDerivative) - this % F_Ur) / eps + shift * x
!~       Ax = ( this % p_F(this % Ur + x * eps) - this % p_F(this % Ur - x * eps))  /(2._RP * eps)  - x / dtsolve   !Second order
      
      ! *NORM2(x)
   END FUNCTION AxMult
   
   FUNCTION p_F(this,u, computeTimeDerivative) RESULT(F)
      IMPLICIT NONE
      CLASS(MatFreeSmooth_t), INTENT(INOUT) :: this
      REAL(KIND = RP), INTENT(IN)             :: u(:)
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      REAL(KIND = RP)                         :: F(size(u))
      
      this % p_sem % mesh % storage % Q = u
      call this % p_sem % mesh % storage % global2LocalQ
      CALL ComputeTimeDerivative(this % p_sem % mesh, this % p_sem % particles, timesolve, CTD_IGNORE_MODE)
      call this % p_sem % mesh % storage % local2GlobalQdot(this % p_sem % NDOF)
      F = this % p_sem % mesh % storage % Qdot
      
   END FUNCTION p_F
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  Smoothers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE BlockJacobiSmoother(this, SmoothIters, niter, ComputeTimeDerivative, TolPresent, tol)
      USE DenseMatUtilities
      IMPLICIT NONE
      !--------------------------------------------
      CLASS(MatFreeSmooth_t), TARGET, intent(INOUT) :: this            !<  Iterative solver class
      INTEGER                                       :: SmoothIters     !<  Number of smoothing operations
      INTEGER                       , intent(OUT)   :: niter           !>  Number of iterations needed
      procedure(ComputeTimeDerivative_f)                    :: ComputeTimeDerivative
      logical                       , intent(in)    :: TolPresent      !   
      REAL(KIND=RP), OPTIONAL                       :: tol             !   Relative AND absolute tolerance of the method
      !--------------------------------------------
       INTEGER                                 :: n                ! System size
      REAL(KIND=RP)                           :: r   (this % DimPrb) ! Residual
      REAL(KIND=RP)                           :: P_1r(this % DimPrb) ! Residual
      REAL(KIND=RP), POINTER                  :: x(:)             ! Solution
      REAL(KIND=RP), POINTER                  :: b(:)             ! Right-hand-side
      INTEGER                                 :: i,j              ! Counters
      INTEGER                                 :: idx1, idx2       ! Indexes of block
      
      REAL(KIND=RP)                           :: bnorm, rnorm, oldrnorm, ConvRate     ! Norm of b and r vectors
      REAL(KIND=RP)                           :: endtol           ! Final tolerance that will be used to evaluate convergence 
      !--------------------------------------------
      
      n =  this % DimPrb
      x => this % x
      b => this % b
      
      IF(TolPresent) THEN
         bnorm = NORM2(b)
         endtol = MAX(bnorm*tol,tol)  ! rtol and atol are taken as the same value
      END IF
      
!~      print*, 'bnorm = ', bnorm
!~      print*, '    iter      residual'
      
      oldrnorm = -1._RP
      ConvRate = 1._RP
      DO i=1,SmoothIters
         r = this % AxMult(x, computeTimeDerivative)        ! Matrix free mult
         
!$omp parallel do private(idx1,idx2) schedule(runtime)
         DO j=1, this % A % num_of_Blocks
            idx1 = this % A % BlockIdx(j)
            idx2 = this % A % BlockIdx(j+1)-1

            r(idx1:idx2) = b(idx1:idx2) - r(idx1:idx2)
            call SolveLU(ALU      = this%BlockPreco(j) % PLU, &
                         LUpivots = this%BlockPreco(j) % LUpivots, &
                         x = P_1r(idx1:idx2), &
                         b = r   (idx1:idx2))
            
            x(idx1:idx2) = x(idx1:idx2) + P_1r(idx1:idx2)
         END DO
!$omp end parallel do
         
         IF (TolPresent) THEN
            rnorm = NORM2(r)       ! Saves relative tolerance (one iteration behind)
!~             print*, '\x1b[1;34m', i, rnorm, rnorm/oldrnorm ,'\x1b[0m'
!~             read(*,*)
            IF (oldrnorm > 0._RP) THEN
               ConvRate = ConvRate + (LOG10(oldrnorm/rnorm)-ConvRate)/i 
            ENDIF
            IF (rnorm < endtol .or. ConvRate <= 8e-3_RP) then ! .OR. ABS(rnorm/oldrnorm-1._RP) < 0.01_RP) THEN
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
      CLASS(MatFreeSmooth_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
      !-------------------------------------------------------------
      INTEGER :: k      ! Counter
      !-------------------------------------------------------------
!$omp parallel do schedule(runtime)
      DO k=1, this % A % num_of_Blocks
         call ComputeLU (A        = this % A % Blocks(k) % Matrix, &
                         ALU      = this % BlockPreco(k) % PLU, &
                         LUpivots = this % BlockPreco(k) % LUpivots)
      END DO
!$omp end parallel do
      
   END SUBROUTINE ComputeBlockPreco
END MODULE MatrixFreeSmootherClass