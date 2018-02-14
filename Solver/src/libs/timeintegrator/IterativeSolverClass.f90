!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      IterativeClass.f90
!      Created: 2017-04-XX 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for solving a system with simple smoothers and matrix free operations
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE IterativeSolverClass
   USE GenericLinSolverClass
   USE CSR_Matrices
   USE SMConstants
   USE PetscSolverClass   ! For allocating Jacobian matrix
   use DGSEMClass
   use TimeIntegratorDefinitions
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC IterativeSolver_t, GenericLinSolver_t
   
   TYPE :: BlockPreco_t
      REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: Pinv     ! Inverse of elemental preconditioner matrix
   END TYPE BlockPreco_t
   
   TYPE, EXTENDS(GenericLinSolver_t) :: IterativeSolver_t
      TYPE(csrMat_t)                             :: A                                  ! Jacobian matrix
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: F_Ur                               ! Qdot at the beginning of solving procedure
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: Ur                                 ! Q at the beginning of solving procedure
      REAL(KIND=RP)                              :: rnorm                              ! L2 norm of residual
      REAL(KIND=RP)                              :: Ashift                             ! Shift that the Jacobian matrix currently(!) has
      LOGICAL                                    :: AIsPrealloc                        ! Has A already been preallocated? (to be deprecated)
      
      !Variables for creating Jacobian in PETSc context:
      TYPE(PetscKspLinearSolver_t)               :: PetscSolver                        ! PETSc solver (created only for allocating the matrix -to be deprecated?)
      LOGICAL                                    :: AIsPetsc = .TRUE.
      
      TYPE(DGSem), POINTER                       :: p_sem                              ! Pointer to DGSem class variable of current system
      CHARACTER(LEN=LINE_LENGTH)                 :: Smoother
      TYPE(BlockPreco_t), ALLOCATABLE            :: BlockPreco(:)
   CONTAINS
      !Subroutines:
      PROCEDURE                                  :: construct
      PROCEDURE                                  :: PreallocateA
      PROCEDURE                                  :: ResetA
      PROCEDURE                                  :: SetAColumn 
      PROCEDURE                                  :: AssemblyA
      PROCEDURE                                  :: SetBValue
      PROCEDURE                                  :: SetBValues
      PROCEDURE                                  :: solve
      PROCEDURE                                  :: GetCSRMatrix
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
      PROCEDURE                                  :: ComputeBlockPreco
      
      PROCEDURE                                  :: p_F
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
   SUBROUTINE construct(this,DimPrb,controlVariables,sem)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t) , INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      !-----------------------------------------------------------
      INTEGER :: nelem      ! Number of elements
      INTEGER :: Nx,Ny,Nz   ! Polynomial orders for element
      INTEGER :: ndofelm    ! Number of degrees of freedom of element
      INTEGER :: k          ! Counter  
      INTEGER :: N_EQN                                       
      !-----------------------------------------------------------
      
      IF (.NOT. PRESENT(sem)) stop 'Fatal error: IterativeSolver needs sem.'
      
      
      this % DimPrb = DimPrb
      this % Smoother = controlVariables % StringValueForKey("smoother",LINE_LENGTH)
      
      ALLOCATE(this % x   (DimPrb))
      ALLOCATE(this % b   (DimPrb))
      ALLOCATE(this % F_Ur(DimPrb))
      ALLOCATE(this % Ur  (DimPrb))
      
      IF(this % AIsPetsc) CALL this % PetscSolver % construct (DimPrb)
      
      this % p_sem => sem
!
!     ------------------------------------------------
!     Allocate important variables for preconditioners
!     ------------------------------------------------
!
      SELECT CASE (this % Smoother)
         CASE('BlockJacobi')
            nelem = SIZE(sem % mesh % elements)
            N_EQN = SIZE(sem % mesh % elements(1) % storage % Q,4)
            ALLOCATE (this % BlockPreco(nelem))
            DO k = 1, nelem
               Nx = sem % mesh % elements(k) % Nxyz(1)
               Ny = sem % mesh % elements(k) % Nxyz(2)
               Nz = sem % mesh % elements(k) % Nxyz(3)
               ndofelm = N_EQN*(Nx+1)*(Ny+1)*(Nz+1)
               ALLOCATE (this % BlockPreco(k) % Pinv(ndofelm,ndofelm))
            END DO
      END SELECT
   END SUBROUTINE construct
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE PreallocateA(this,nnz)
      IMPLICIT NONE
!
!  --------------------------------------------
!  Subroutine for preallocating Jacobian matrix
!     IMPORTANT:
!        - To date, it's using the PETSc library for preallocating due to its efficiency
!        - A proposed alternative to emulate PETSc allocation is to use 
!          a linked list matrix and then transform it to CSR format
!  --------------------------------------------
!
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT)  :: this
      INTEGER                                  :: nnz
      !-----------------------------------------------------------
      LOGICAL, SAVE                            :: isfirst = .TRUE.
      !-----------------------------------------------------------
      
      this % AIsPetsc = .TRUE.      ! This is always the case when creating a new Jacobian
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % PreallocateA(nnz)
      ELSE
         IF (.NOT. this % AIsPrealloc) THEN
            stop 'routines for preallocating matrix not implemented'
         ELSE
            print*, 'Matrix already preallocated'
         END IF
      END IF
      this % AIsPrealloc = .TRUE.
      
   END SUBROUTINE PreallocateA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ResetA(this)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % ResetA
      ELSE
         this % A % Values = 0._RP
      END IF
      
   END SUBROUTINE ResetA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetAColumn(this,nvalues,irow,icol,values)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT)  :: this
      INTEGER       , INTENT(IN)                :: nvalues
      INTEGER       , INTENT(IN), DIMENSION(:)  :: irow
      INTEGER       , INTENT(IN)                :: icol
      REAL(KIND=RP) , INTENT(IN), DIMENSION(:)  :: values
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % SetAColumn(nvalues,irow,icol,values)
      ELSE
         CALL this % A % SetColumn(irow+1,icol+1,values)
      END IF
      
   END SUBROUTINE SetAColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE AssemblyA(this,BlockIdx,BlockSize)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t) , INTENT(INOUT) :: this
      INTEGER, TARGET, OPTIONAL, INTENT(IN)    :: BlockIdx(:)
      INTEGER, TARGET, OPTIONAL, INTENT(IN)    :: BlockSize(:)
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % AssemblyA
         CALL this % PetscSolver % GetCSRMatrix(this % A)
!~         CALL this % PetscSolver % destroy
         this % AIsPetsc = .FALSE.
      ELSE
         print*, 'A is assembled'
      END IF
      
      IF (this % Smoother == 'BlockJacobi') THEN
         IF (.NOT. PRESENT(BlockIdx) .OR. .NOT. PRESENT(BlockSize)) THEN
            STOP 'Iterative Class :AssemblyA :: Block Jacobi smoother needs block information'
         ELSE
            this % A % BlockIdx  => BlockIdx
            this % A % BlockSize => BlockSize
         END IF
      END IF
   end SUBROUTINE AssemblyA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetBValue(this, irow, value)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(IN)    :: value
      !-----------------------------------------------------------
      
      this % b (irow+1) = value
      
   END SUBROUTINE SetBValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SetBValues(this, nvalues, irow, values)
      IMPLICIT NONE
      CLASS(IterativeSolver_t)   , INTENT(INOUT)     :: this
      INTEGER                     , INTENT(IN)        :: nvalues
      INTEGER      , DIMENSION(1:), INTENT(IN)        :: irow
      REAL(KIND=RP), DIMENSION(1:), INTENT(IN)        :: values
      !------------------------------------------------------
      INTEGER                                        :: i
      
      DO i=1, nvalues
         IF (irow(i)<0) CYCLE
         this % b(irow(i)+1) = values(i)
      END DO
      
   END SUBROUTINE SetBValues
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE solve(this,tol,maxiter,time,dt, ComputeTimeDerivative)
      IMPLICIT NONE
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP), OPTIONAL                 :: tol
      INTEGER      , OPTIONAL                 :: maxiter
      REAL(KIND=RP), OPTIONAL                 :: time
      REAL(KIND=RP), OPTIONAL                 :: dt
      procedure(ComputeQDot_FCN)              :: ComputeTimeDerivative
      !-------------------------------------------------
      INTEGER                                 :: i
!~      LOGICAL, SAVE :: isfirst = .TRUE.
      !-------------------------------------------------
      
      IF (.NOT. PRESENT(time) .OR. .NOT. PRESENT(dt)) STOP 'time and dt needed for iterative solver'
      
      timesolve= time
      dtsolve  = dt
      
!~      IF (isfirst) THEN
         CALL this % p_sem % GetQdot(this % F_Ur)
         CALL this % p_sem % GetQ   (this % Ur)
!~         isfirst = .FALSE.
!~      END IF
      
      ! Initialize x
      DO i=1, this % DimPrb
         this % x(i) = this % b(i) / this % A % Values(this%A%Diag(i))
      END DO
      
      SELECT CASE (this % Smoother)
         CASE('WeightedJacobi')
            CALL this % WeightedJacobiSmoother( this%A%Values(this%A%Diag), maxiter, tol, this % niter, ComputeTimeDerivative)
         CASE('BlockJacobi')
            CALL BlockJacobiSmoother(this, maxiter, tol, this % niter, ComputeTimeDerivative)
      END SELECT
      
      CALL this % p_sem % SetQ   (this % Ur)
      
!~      IF (this % niter < maxiter) THEN
         this % CONVERGED = .TRUE.
!~      ELSE
!~         this % CONVERGED = .FALSE.
!~      END IF
      
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
      
      x_i = this % x(irow+1)
      
   END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destroy(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      CALL this % A % destroy()
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      
      IF (this % AIsPetsc) CALL this % PetscSolver % destroy()
      this % AIsPetsc    = .TRUE.
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
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = -1._RP/dt
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % SetOperatorDt(dt)
      ELSE
         CALL this % A % SetMatShift(this % Ashift)
      END IF
      
      IF(this % Smoother == 'BlockJacobi') CALL this % ComputeBlockPreco
      
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
      
      shift = -1._RP/dt
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % ReSetOperatorDt(dt)
      ELSE
         CALL this % A % SetMatShift(-this % Ashift)
         CALL this % A % SetMatShift(shift)
      END IF
      this % Ashift = shift
      
      IF(this % Smoother == 'BlockJacobi') CALL this % ComputeBlockPreco
      
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
            STOP 'MKLPardisoSolverClass ERROR: Norm not implemented yet'
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
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                           :: x (:)
      procedure(ComputeQDot_FCN)              :: ComputeTimeDerivative
      REAL(KIND=RP)                           :: Ax(size(x))
      
      !--------------------------------------------------
!~      REAL(KIND=RP)                           :: eps
      REAL(KIND=RP)                           :: F (size(x))
      
!~       REAL(KIND=RP)                           :: xxx (size(x)) !x vector... But normalized!!
!~      REAL(KIND=RP)                           :: buffer (size(x))
      !--------------------------------------------------
      
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
      
      eps = SQRT(EPSILON(eps) * (1._RP + NORM2(this % Ur))) * DOT_PRODUCT(this % Ur,x) / (NORM2(x)) ! NISTOL Package modified    !! This works down to 10⁻¹⁰, but slowly (linear residual does not always go as low as desired).. At the end short tie-steps are needed
      
!~       eps = SQRT(EPSILON(eps)) * (NORM2(this % Ur)*DOT_PRODUCT(this % Ur,x)) / NORM2(x) ! My recipe.. goes lower but slower
      
!~      CALL this % p_sem % GetQ(buffer)

!~       xxx = x / NORM2(x)

!~       CALL this % p_sem % SetQ(this % Ur + x*eps)
!~       CALL ComputeTimeDerivative(this % p_sem,timesolve)
!~       CALL this % p_sem % GetQdot(F)
!~      CALL this % p_sem % SetQ(buffer)
!~       Ax = ( F - this % F_Ur) / eps - x / (dtsolve)                          !First order   ! arueda: this is defined only for BDF1
      Ax = ( this % p_F(this % Ur + x * eps, computeTimeDerivative) - this % F_Ur) / eps - x / dtsolve
!~       Ax = ( this % p_F(this % Ur + x * eps) - this % p_F(this % Ur - x * eps))  /(2._RP * eps)  - x / dtsolve   !Second order
      
      ! *NORM2(x)
   END FUNCTION AxMult
   
   FUNCTION p_F(this,u, computeTimeDerivative) RESULT(F)
      IMPLICIT NONE
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND = RP), INTENT(IN)             :: u(:)
      procedure(ComputeQDot_FCN)              :: ComputeTimeDerivative
      REAL(KIND = RP)                         :: F(size(u))
      
      CALL this % p_sem % SetQ(u)
      CALL ComputeTimeDerivative(this % p_sem % mesh, timesolve, this % p_sem % externalState, this % p_sem % externalGradients)
      CALL this % p_sem % GetQdot(F)
      
   END FUNCTION p_F
   
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
      procedure(ComputeQDot_FCN)                      :: ComputeTimeDerivative
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
        r = this % AxMult(x, ComputeTimeDerivative)        ! Matrix free mult
!~          r = CSR_MatVecMul(this%A,x) ! CSR matrix product
         
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
      IMPLICIT NONE
      !--------------------------------------------
      CLASS(IterativeSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
      INTEGER                                         :: SmoothIters     !<  Number of smoothing operations
      REAL(KIND=RP), OPTIONAL                         :: tol             !   Relative AND absolute tolerance of the method
      INTEGER                         , INTENT(OUT)   :: niter           !>  Number of iterations needed
      procedure(ComputeQDot_FCN)                      :: ComputeTimeDerivative
      !--------------------------------------------
       INTEGER                                 :: n                ! System size
      REAL(KIND=RP)                           :: r(this % DimPrb) ! Residual
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
      DO i=1,SmoothIters
         r = this % AxMult(x, computeTimeDerivative)        ! Matrix free mult
!~          r = CSR_MatVecMul(this%A,x) ! CSR matrix product
         
!$omp parallel do private(idx1,idx2)
         DO j=1,SIZE(this % A % BlockSize)
            idx1 = this % A % BlockIdx(j) + 1   ! zero-based indexing
            idx2 = this % A % BlockIdx(j+1)
            
            r(idx1:idx2) = b(idx1:idx2) - r(idx1:idx2)
            x(idx1:idx2) = x(idx1:idx2) + MATMUL(this%BlockPreco(j)%Pinv,r(idx1:idx2))
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
      print*, '\x1b[1;34mSmoother ConvRate:', ConvRate ,'\x1b[0m'
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
      INTEGER :: nelem
      INTEGER :: k, N      ! Counter / size of block
      !-------------------------------------------------------------
      
      nelem = SIZE (this % A % BlockSize)
      DO k=1, nelem
         N = this % A % BlockSize(k)
         this % BlockPreco(k) % Pinv = inverse(this % A % GetBlock(k,N))
      END DO
      
   END SUBROUTINE ComputeBlockPreco
END MODULE IterativeSolverClass
