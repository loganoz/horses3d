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
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC IterativeSolver_t, GenericLinSolver_t
   
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
      
      !!PROCEDURE(matmultfun), POINTER, NOPASS     :: MatMult => NULL()                  ! arueda: This is needed here for lower multigrid levels
      !TYPE(IterativeSolver_t), POINTER :: Child ! This will be needed for Multigrid (remove from here)
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
   SUBROUTINE construct(this,DimPrb,sem)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t) , INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: DimPrb
      TYPE(DGSem), TARGET      , OPTIONAL      :: sem
      !-----------------------------------------------------------
      
      IF (.NOT. PRESENT(sem)) stop 'Fatal error: IterativeSolver needs sem.'
      
      
      this % DimPrb = DimPrb
      
      ALLOCATE(this % x   (DimPrb))
      ALLOCATE(this % b   (DimPrb))
      ALLOCATE(this % F_Ur(DimPrb))
      ALLOCATE(this % Ur  (DimPrb))
      
      IF(this % AIsPetsc) CALL this % PetscSolver % construct (DimPrb)
      
      this % p_sem => sem
      
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
   SUBROUTINE AssemblyA(this)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % AssemblyA
         CALL this % PetscSolver % GetCSRMatrix(this % A)
         CALL this % PetscSolver % destroy
         this % AIsPetsc = .FALSE.
      ELSE
         print*, 'A is assembled'
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
   SUBROUTINE solve(this,tol,maxiter,time,dt)
      IMPLICIT NONE
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP), OPTIONAL                 :: tol
      INTEGER      , OPTIONAL                 :: maxiter
      REAL(KIND=RP), OPTIONAL                 :: time
      REAL(KIND=RP), OPTIONAL                 :: dt
      !-------------------------------------------------
      INTEGER                                 :: niter
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
      
      CALL this % WeightedJacobiSmoother( this%A%Values(this%A%Diag), maxiter, tol, this % niter)
      
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

   FUNCTION AxMult(this,x) RESULT(Ax)
      IMPLICIT NONE
      CLASS(IterativeSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                           :: x (:)
      REAL(KIND=RP)                           :: Ax(size(x))
      !--------------------------------------------------
!~      REAL(KIND=RP)                           :: eps
      REAL(KIND=RP)                           :: F (size(x))
!~      REAL(KIND=RP)                           :: buffer (size(x))
      !--------------------------------------------------
      
!~      eps = 1e-8_RP * (1._RP + NORM2(x))                           ! ~2e-5 2e-4
!~      eps = 1e-8_RP * (1._RP + NORM2(this % Ur))                   ! better: ~6e-7
      eps = SQRT(EPSILON(eps)) * (1._RP + NORM2(this % Ur))        !slightly better: ~4e-7 
!~      eps = SQRT(EPSILON(eps)) * (1._RP + MAXVAL(ABS(this % Ur)))  !slightly worse: ~1e-5 9e-6
!~      eps = SQRT(EPSILON(eps))                                     !worse:        : ~1e-4
      
!~      CALL this % p_sem % GetQ(buffer)
      CALL this % p_sem % SetQ(this % Ur + x*eps)
      CALL ComputeTimeDerivative(this % p_sem,timesolve)
      CALL this % p_sem % GetQdot(F)
!~      CALL this % p_sem % SetQ(buffer)
      Ax = ( F - this % F_Ur) / eps - x / (dtsolve)                          !First order   ! arueda: this is defined only for BDF1
      
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
      CLASS(IterativeSolver_t), TARGET, INTENT(INOUT) :: this            !<  Iterative solver class
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
!~            print*, i, rnorm
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
   
END MODULE IterativeSolverClass
