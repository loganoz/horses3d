!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      IterativeClass.f90
!      Created: 2017-04-XX 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for solving a system with simple smoothers and matrix free operations
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE MultigridSolverClass
   USE GenericLinSolverClass
   USE CSR_Matrices
   USE SMConstants
   USE PetscSolverClass   ! For allocating Jacobian matrix
   
   ! To be moved to other module?
   USE PolynomialInterpAndDerivsModule
   USE GaussQuadrature
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
      
      !Variables for creating Jacobian in PETSc context:
      TYPE(PetscKspLinearSolver_t)               :: PetscSolver                        ! PETSc solver (created only for allocating the matrix -to be deprecated?)
      LOGICAL                                    :: AIsPetsc = .TRUE.
      
      TYPE(DGSem)            , POINTER           :: p_sem                              ! Pointer to DGSem class variable of current system
      
      ! Variables that are specially needed for Multigrid
      TYPE(MultigridSolver_t), POINTER           :: Child                 ! Next (coarser) multigrid solver
      INTEGER                                    :: MGlevel               ! Current Multigrid level
      TYPE(Interpolator_t)       , ALLOCATABLE   :: Restriction(:,:,:)    ! Restriction operators (element level)
      TYPE(Interpolator_t)       , ALLOCATABLE   :: Prolongation(:,:,:)   ! Prolongation operators (element level)
      
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
   
CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE construct(this,DimPrb,controlVariables,sem)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MultigridSolver_t) , INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      !-----------------------------------------------------------
      INTEGER                           :: lvl              ! Multigrid level
      TYPE(MultigridSolver_t) , POINTER :: Solver_p         ! Pointer to the current solver (being constructed)
      !-----------------------------------------------------------
      !Module variables: MGlevels, deltaN
      
      IF (.NOT. PRESENT(sem)) stop 'Fatal error: MultigridSolver needs sem.'
      IF (.NOT. PRESENT(controlVariables)) stop 'Fatal error: MultigridSolver needs controlVariables.'
      
!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
      IF (.NOT. controlVariables % containsKey("multigrid levels")) THEN
         print*, 'Fatal error: "multigrid levels" keyword is needed by the multigrid solver'
         STOP
      END IF
      MGlevels  = controlVariables % IntegerValueForKey("multigrid levels")
      
      IF (controlVariables % containsKey("delta n")) THEN
         deltaN = controlVariables % IntegerValueForKey("delta n")
      ELSE
         deltaN = 1
      END IF
!
!     --------------------------
!     Create linked solvers list
!     --------------------------
!
      Solver_p => this
      
      DO lvl = MGlevels, 1, -1
         
         IF (lvl /= 1) THEN
            ALLOCATE (Solver_p % Child)
            Solver_p => Solver_p % Child
         END IF
      END DO
      
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
      CLASS(MultigridSolver_t), INTENT(INOUT)  :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT)  :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t)   , INTENT(INOUT)     :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t), INTENT(IN)  :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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

   FUNCTION AxMult(this,x) RESULT(Ax)
      IMPLICIT NONE
      CLASS(MultigridSolver_t), INTENT(INOUT) :: this
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
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  Routines for interpolation procedures (will probably be moved to another module)
!  
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE CreateInterpolationOperators(Restriction,Prolongation,N1x,N1y,N1z,DeltaN,Lobatto)
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
      INTEGER                       :: DeltaN                  !<  Interval of reduction of polynomial order for coarser level
      LOGICAL, OPTIONAL             :: Lobatto                 !<  Is the quadrature a Legendre-Gauss-Lobatto representation?
      !-----------------------------------------------------
      TYPE(Interpolator_t), POINTER :: rest                    ! Pointer to constructed restriction interpolator
      TYPE(Interpolator_t), POINTER :: prol                    ! Pointer to constructed prolongation interpolator
      INTEGER                       :: N2x, N2y, N2z           ! Coarse grid order(anisotropic) of the element
      LOGICAL                       :: LGL = .FALSE.           ! Is the quadrature a Legendre-Gauss-Lobatto representation? (false is default)
      INTEGER                       :: i,j,k,l,m,n             ! Index counters
      INTEGER                       :: s,r                     ! Row/column counters for operators
      REAL(KIND=RP), ALLOCATABLE    :: x1 (:), y1 (:), z1 (:)  ! Position of quadrature points on mesh 1
      REAL(KIND=RP), ALLOCATABLE    :: w1x(:), w1y(:), w1z(:)  ! Weights for quadrature points on mesh 1
      REAL(KIND=RP), ALLOCATABLE    :: x2 (:), y2 (:), z2 (:)  ! Position of quadrature points on mesh 2
      REAL(KIND=RP), ALLOCATABLE    :: w2x(:), w2y(:), w2z(:)  ! Weights for quadrature points on mesh 2
      !-----------------------------------------------------
      
      IF (PRESENT(Lobatto) .AND. Lobatto) LGL = .TRUE.
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
! TODO: Remove this from here... The implementation is in Interpolation module
   SUBROUTINE Interpolate1Eqn(Q1, Q2, Interp, N1x, N1y, N1z, N2x, N2y, N2z)
      IMPLICIT NONE
      INTEGER        :: N1x, N1y, N1z
      INTEGER        :: N2x, N2y, N2z                      !<  Polynomial orders
      REAL(KIND=RP)  :: Q1((N1x+1)*(N1y+1)*(N1z+1))        !<  Solution to be interpolated (grid (1))
      REAL(KIND=RP)  :: Q2((N2x+1)*(N2y+1)*(N2z+1))        !>  Interpolated solution       (grid (2))
      REAL(KIND=RP)  :: Interp((N2x+1)*(N2y+1)*(N2z+1), & 
                               (N1x+1)*(N1y+1)*(N1z+1))    !<  Interpolation matrix
                               
      Q2 = MATMUL(Interp,Q1)
   END SUBROUTINE Interpolate1Eqn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   
END MODULE MultigridSolverClass
