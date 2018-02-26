!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      IterativeClass.f90
!      Created: 2017-04-XX 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for solving a system with simple smoothers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE IterativeSolverClass
   use MatrixClass
   USE GenericLinSolverClass
   USE CSRMatrixClass
   USE SMConstants
   USE PetscSolverClass   ! For allocating Jacobian matrix
   use DGSEMClass
   use TimeIntegratorDefinitions
   use NumericalJacobian
   use AnalyticalJacobian
   use PhysicsStorage
   IMPLICIT NONE
#ifdef HAS_PETSC
#include <petsc.h>
#endif
   PRIVATE
   PUBLIC IterativeSolver_t, GenericLinSolver_t
   
   TYPE :: BlockPreco_t
      real(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: PLU        ! LU factorization of elemental preconditioner matrix
      integer      , dimension(:)  , allocatable :: LUpivots   ! LU pivots
   END TYPE BlockPreco_t
   
   TYPE, EXTENDS(GenericLinSolver_t) :: IterativeSolver_t
      TYPE(csrMat_t)                             :: A                                  ! Jacobian matrix
      type(PETSCMatrix_t)                        :: PETScA
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
      
      TYPE(DGSem), POINTER                       :: p_sem                          ! Pointer to DGSem class variable of current system
      CHARACTER(LEN=LINE_LENGTH)                 :: Smoother
      TYPE(BlockPreco_t), ALLOCATABLE            :: BlockPreco(:)
   CONTAINS
      !Subroutines:
      PROCEDURE                                  :: construct
      PROCEDURE                                  :: SetBValue
      PROCEDURE                                  :: SetBValues
      PROCEDURE                                  :: solve
      PROCEDURE                                  :: GetCSRMatrix
      PROCEDURE                                  :: GetXValue
      PROCEDURE                                  :: destroy
      PROCEDURE                                  :: SetOperatorDt
      PROCEDURE                                  :: ReSetOperatorDt
      procedure                                  :: FillAInfo
      !Functions:
      PROCEDURE                                  :: Getxnorm    !Get solution norm
      PROCEDURE                                  :: Getrnorm    !Get residual norm
      
      !! Internal procedures
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
#ifdef HAS_PETSC
      PetscErrorCode :: ierr
      CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#endif
      IF (.NOT. PRESENT(sem)) stop 'Fatal error: IterativeSolver needs sem.'
      
      
      this % DimPrb = DimPrb
      this % Smoother = controlVariables % StringValueForKey("smoother",LINE_LENGTH)
      
      ALLOCATE(this % x   (DimPrb))
      ALLOCATE(this % b   (DimPrb))
      ALLOCATE(this % F_Ur(DimPrb))
      ALLOCATE(this % Ur  (DimPrb))

      IF(this % AIsPetsc) CALL this % PETScA % construct (DimPrb,.FALSE.)
      
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
               ndofelm = NCONS*(Nx+1)*(Ny+1)*(Nz+1)
               allocate (this % BlockPreco(k) % PLU(ndofelm,ndofelm) )
               allocate (this % BlockPreco(k) % LUpivots   (ndofelm) )
            END DO
      END SELECT
   END SUBROUTINE construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FillAInfo(this)
      implicit none
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      integer :: i
      integer, target, allocatable, dimension(:) :: Nx, Ny, Nz, ndofelm
      integer, target, allocatable :: firstIdx(:)
      integer :: nelem
      !-----------------------------------------------------------
      nelem= size(this % p_sem % mesh % elements)
      
      allocate ( Nx(nelem), Ny(nelem), Nz(nelem), ndofelm(nelem) )
      allocate ( firstIdx(nelem+1) )
      
      firstIdx(1) = 1
      
      nelem = size(this % p_sem % mesh % elements)
      do i=1, nelem
            Nx(i) = this % p_sem % mesh % elements(i) % Nxyz(1)
            Ny(i) = this % p_sem % mesh % elements(i) % Nxyz(2)
            Nz(i) = this % p_sem % mesh % elements(i) % Nxyz(3)
            
!           Get block sizes and position in matrix
!              TODO: change to store the permutation indexes in the element
!           --------------------------------------
! 
            ndofelm(i)  = N_EQN * (Nx(i)+1) * (Ny(i)+1) * (Nz(i)+1)
            IF (i>1) firstIdx(i) = firstIdx(i-1) + ndofelm(i-1)
      end do
      firstIdx(nelem+1) = firstIdx(nelem) + ndofelm(nelem)
      
      allocate (this % A % BlockIdx(nelem+1))
      allocate (this % A % BlockSize(nelem))
      this % A % BlockIdx  = firstIdx
      this % A % BlockSize = ndofelm
   end subroutine
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
   SUBROUTINE solve(this, ComputeTimeDerivative,tol,maxiter,time,dt, ComputeA)
   use DenseMatUtilities
      IMPLICIT NONE
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      procedure(ComputeQDot_FCN)              :: ComputeTimeDerivative
      REAL(KIND=RP), OPTIONAL                 :: tol
      INTEGER      , OPTIONAL                 :: maxiter
      REAL(KIND=RP), OPTIONAL                 :: time
      REAL(KIND=RP), OPTIONAL                 :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-------------------------------------------------
      INTEGER                                 :: i
!~      LOGICAL, SAVE :: isfirst = .TRUE.
      !-------------------------------------------------
      !<temp
!~      real(kind=RP), allocatable :: A(:,:)
!~      integer :: blockn
      !temp>
      IF (.NOT. PRESENT(time) .OR. .NOT. PRESENT(dt)) STOP 'time and dt needed for iterative solver'
      
!
!     Compute Jacobian matrix if needed
!        (done in petsc format and then transformed to CSR since the CSR cannot be filled by the Jacobian calculators)
!     -----------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
!~            call AnalyticalJacobian_Compute(this % p_sem,time,this % PETScA,.TRUE.)
            call NumericalJacobian_Compute(this % p_sem, time, this % PETScA, ComputeTimeDerivative, .TRUE. )
            call this % PETScA % shift(-1._RP/dt)
            call this % PETScA % GetCSRMatrix(this % A)
            if (this % Smoother == 'BlockJacobi') call this % FillAInfo
            IF(this % Smoother == 'BlockJacobi') CALL this % ComputeBlockPreco
            this % AIsPetsc = .FALSE.
            ComputeA = .FALSE.
         end if
      else 
         call NumericalJacobian_Compute(this % p_sem, time, this % A, ComputeTimeDerivative, .TRUE. )
         call this % PETScA % shift(-1._RP/dt)
         call this % PETScA % GetCSRMatrix(this % A)
         if (this % Smoother == 'BlockJacobi') call this % FillAInfo
      end if
      
!~      blockn = 1
!~      allocate ( A(this % A % BlockSize(blockn),this % A % BlockSize(blockn)) )
!~      A = this % A % GetBlock(blockn,this % A % BlockSize(blockn))
!~      call Mat2File(A,'NewNum.dat')
!~      stop
      
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
      
      CALL this % A % destruct()
      
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
         CALL this % A % Shift(this % Ashift)
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
         CALL this % A % Shift(-this % Ashift)
         CALL this % A % Shift(shift)
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
      procedure(ComputeQDot_FCN)                      :: ComputeTimeDerivative
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
      DO i=1,SmoothIters
            print*, 'MatMult'
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
!~   SUBROUTINE ComputeBlockPreco(this)
!~      USE DenseMatUtilities
!~      IMPLICIT NONE
!~      !-------------------------------------------------------------
!~      CLASS(IterativeSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
!~      !-------------------------------------------------------------
!~      INTEGER :: nelem
!~      INTEGER :: k, N      ! Counter / size of block
!~      !-------------------------------------------------------------
      
!~      nelem = SIZE (this % A % BlockSize)
!~      DO k=1, nelem
!~         N = this % A % BlockSize(k)
!~         this % BlockPreco(k) % Pinv = inverse(this % A % GetBlock(k,N))
!~      END DO
      
!~   END SUBROUTINE ComputeBlockPreco
   
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
         N = this % A % BlockSize(k)
         call ComputeLU (A        = this % A % GetBlock(k,N), &
                         ALU      = this % BlockPreco(k) % PLU, &
                         LUpivots = this % BlockPreco(k) % LUpivots)
      END DO
!$omp end parallel do
      
   END SUBROUTINE ComputeBlockPreco
END MODULE IterativeSolverClass
