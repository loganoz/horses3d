!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      PetscSolverClass.f90
!      Created: 2017-04-10 10:006:00 +0100 
!      By: Carlos Redondo
!          Andrés Rueda (small modifications)
!
!      Class for solving linear systems using the Krylov Subspace Methods of PETSc library
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE PetscSolverClass
   USE GenericLinSolverClass
   use MatrixClass
   USE CSRMatrixClass
   USE SMConstants
   use DGSEMClass
   use TimeIntegratorDefinitions
   use NumericalJacobian
   IMPLICIT NONE
#ifdef HAS_PETSC
#include <petsc.h>
#endif
   TYPE, EXTENDS(GenericLinSolver_t) :: PetscKspLinearSolver_t
      type(PETSCMatrix_t), allocatable              :: A
      TYPE(DGSem), POINTER                          :: p_sem   
#ifdef HAS_PETSC
      Mat                                           :: AA                                 ! Jacobian matrix
      Vec                                           :: x                                  ! Solution vector
      Vec                                           :: b                                  ! Right hand side
      KSP                                           :: ksp                                ! 
      PC                                            :: pc
      PetscReal                                     :: tol = 1d-15
      PetscInt                                      :: maxiter = 50
      PetscInt                                      :: nz = 0
      PetscScalar                                   :: Ashift                              ! Stores the shift to the Jacobian due to time integration
      PetscBool                                     :: init_context = PETSC_FALSE
      PetscBool                                     :: setpreco = PETSC_TRUE
      PetscBool                                     :: withMPI = PETSC_FALSE
#endif
      CONTAINS
         !Subroutines
         PROCEDURE                                  :: construct => ConstructPetscContext
         PROCEDURE                                  :: SetBValues
         PROCEDURE                                  :: SetBValue
         PROCEDURE                                  :: GetXValues
         PROCEDURE                                  :: GetXValue
         PROCEDURE                                  :: SetOperatorDt
         PROCEDURE                                  :: ReSetOperatorDt
         PROCEDURE                                  :: AssemblyB
         PROCEDURE                                  :: SaveMat
         PROCEDURE                                  :: solve   =>   SolveLinPrb
         PROCEDURE                                  :: destroy =>   DestroyPetscObjects
         !Functions
         PROCEDURE                                  :: Getxnorm
         PROCEDURE                                  :: Getrnorm
   END TYPE PetscKspLinearSolver_t
   
  
   PRIVATE                                          
   PUBLIC                                           :: PetscKspLinearSolver_t, GenericLinSolver_t
   
   

!========
 CONTAINS
!========

!///////////////////////////////////////////////////////////////////////////////
   SUBROUTINE CheckPetscErr(ierr,msg)
      IMPLICIT NONE
      CHARACTER(LEN=*), OPTIONAL                   :: msg
#ifdef HAS_PETSC
      PetscErrorCode, INTENT(IN)                   :: ierr
      
      IF(ierr .EQ. 0) THEN
         RETURN
      ELSE
         IF (.NOT. PRESENT(msg)) msg = 'error in petsc'
         WRITE(*,*) msg,' **** Petsc call returned an error. Code: ' ,ierr
         STOP
      ENDIF
#else
      INTEGER                                      :: ierr
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE CheckPetscErr
!
!/////////////////////////////////////////////////////////////////////////////// 
!
   SUBROUTINE ConstructPetscContext(this, DimPrb,controlVariables,sem)
      IMPLICIT NONE
      !--------------------------------------------------------------
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT), TARGET :: this
      TYPE(FTValueDictionary)      , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                      , OPTIONAL  :: sem
      !--------------------------------------------------------------
#ifdef HAS_PETSC
      PetscInt, INTENT(IN)                                 :: DimPrb
      PetscErrorCode                                       :: ierr
      
      this % p_sem => sem
      
      !Initialisation of the PETSc variables
      CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)

!     PETSc matrix A 
      allocate (PETSCMatrix_t :: this % A)
      call this % A % construct(DimPrb,this % WithMPI)

!     Petsc vectors x and b (of A x = b)
      CALL VecCreate(PETSC_COMM_WORLD,this%x,ierr)          ; CALL CheckPetscErr(ierr,'error creating Petsc vector')
      CALL VecSetSizes(this%x,PETSC_DECIDE,dimPrb,ierr)     ; CALL CheckPetscErr(ierr,'error setting Petsc vector options')
      CALL VecSetFromOptions(this%x,ierr)                   ; CALL CheckPetscErr(ierr,'error setting Petsc vector options')
      CALL VecDuplicate(this%x,this%b,ierr)                 ; CALL CheckPetscErr(ierr,'error creating Petsc vector')

!     Petsc ksp solver context      
      CALL KSPCreate(PETSC_COMM_WORLD,this%ksp,ierr)                    ; CALL CheckPetscErr(ierr,'error in KSPCreate')

!     Petsc preconditioner 
      CALL KSPGetPC(this%ksp,this%pc,ierr)                              ; CALL CheckPetscErr(ierr,'error in KSPGetPC')

!     Preset ksp solver
      CALL KSPSetOperators(this%ksp, this%A%A, this%A%A, ierr)              ; CALL CheckPetscErr(ierr,'error in KSPSetOperators')
!~      CALL KSPSetTolerances(this%ksp,PETSC_DEFAULT_REAL,this%tol,PETSC_DEFAULT_REAL,this%maxiter,ierr)
!~      CALL CheckPetscErr(ierr,'error in KSPSetTolerances')
      this%init_context = PETSC_TRUE
      this%dimprb = DimPrb
#else
      INTEGER, INTENT(IN)                       :: DimPrb
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE SetPetscPreco(this)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT)      :: this
#ifdef HAS_PETSC
      PetscErrorCode                                  :: ierr
      
      CALL PCSetType(this%pc,PCILU,ierr)       ;CALL CheckPetscErr(ierr, 'error in PCSetType')
      this%setpreco = PETSC_FALSE
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SetPetscPreco 
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SolveLinPrb(this, ComputeTimeDerivative, tol, maxiter, time,dt, ComputeA)
      IMPLICIT NONE
      !-------------------------------------------------------------
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT) :: this
      procedure(ComputeQDot_FCN)                   :: ComputeTimeDerivative
      REAL(KIND=RP), OPTIONAL                      :: time
      REAL(KIND=RP), OPTIONAL                      :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-------------------------------------------------------------
#ifdef HAS_PETSC
      PetscReal, OPTIONAL                             :: tol
      PetscInt,  OPTIONAL                             :: maxiter
      PetscErrorCode                                  :: ierr
      !-------------------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call NumericalJacobian_Compute(this % p_sem, time, this % A, ComputeTimeDerivative, .TRUE. )
            call this % A % shift(-1._RP/dt)
            ComputeA = .FALSE.
         end if
      else 
         call NumericalJacobian_Compute(this % p_sem, time, this % A, ComputeTimeDerivative, .TRUE. )
         call this % A % shift(-1._RP/dt)
      end if
      
      ! Set , if given, solver tolerance and max number of iterations
      IF (PRESENT(tol)) THEN
         this%tol = tol
      ELSE
         this%tol = PETSC_DEFAULT_REAL
      ENDIF
      
      IF (PRESENT(maxiter)) THEN
         this%maxiter = maxiter
      ELSE
         this%maxiter = PETSC_DEFAULT_INTEGER
      ENDIF
      
      CALL KSPSetTolerances(this%ksp,PETSC_DEFAULT_REAL,this%tol,PETSC_DEFAULT_REAL,this%maxiter,ierr)
      CALL CheckPetscErr(ierr, 'error in KSPSetTolerances')
      
      IF (this%setpreco) THEN    ! if setpreco = true sets the default preconditioner PCJACOBI. If PC already set, setpreco = false
         CALL PCSetType(this%pc,PCILU,ierr)                 ; CALL CheckPetscErr(ierr, 'error in PCSetType')
         this%setpreco = PETSC_FALSE
         !Summary (arueda):
         !  - PCILU: Bad results
         !  - PCJACOBI: Slow convergence
         !  - PCPBJACOBI: (point Jacobi) Almost same convergence as PCJACOBI
         !  - PCBJACOBI: (block Jacobi) Bad results 
      ENDIF

      CALL KSPSetOperators(this%ksp, this%A%A, this%A%A, ierr)     ; CALL CheckPetscErr(ierr, 'error in KSPSetOperators')
      CALL KSPSolve(this%ksp,this%b,this%x,ierr)               ; CALL CheckPetscErr(ierr, 'error in KSPSolve')
      CALL KSPGetIterationNumber(this%ksp,this%niter,ierr)     ; CALL CheckPetscErr(ierr,'error in KSPGetIterationNumber')
!~       CALL KSPGetResidualNorm(this%ksp, this%residual, ierr)   ; CALL CheckPetscErr(ierr,'error in KSPGetResidualNorm')
!~       CALL VecNorm(this%x,NORM_INFINITY,this%xnorm,ierr)       ; CALL CheckPetscErr(ierr,'error in VecNorm')
      IF (this%niter < maxiter) THEN
         this%converged = .TRUE.
      ELSE
         this%converged = .FALSE.
      END IF
#else
      REAL*8 , OPTIONAL                     :: tol
      INTEGER, OPTIONAL                     :: maxiter
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SolveLinPrb
!
!/////////////////////////////////////////////////////////////////////////////////////////////////   
!
   SUBROUTINE SetOperatorDt(this, dt)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      PetscScalar,                     INTENT(IN)        :: dt
      PetscScalar                                        :: shift
      PetscScalar                                        :: eps = 1e-10

      shift = -1._RP/dt !
      IF (ABS(shift) .GT. eps) THEN                  
         call this % A % shift(shift) ! A = A + shift * I
         this % Ashift = shift
      ENDIF
#else
      REAL*8,                     INTENT(IN)        :: dt
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////// 
!   
   SUBROUTINE ReSetOperatorDt(this, dt)
      IMPLICIT NONE
!
!     --------------------------------------------------
!     Removes previous shift in order to insert new one 
!                 (important when Jacobian is reused)
!     --------------------------------------------------
!
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      PetscScalar,                     INTENT(IN)        :: dt
      PetscScalar                                        :: shift
      PetscScalar                                        :: eps = 1e-10

      shift = -1._RP/dt !
      IF (ABS(shift) .GT. eps) THEN
         call this % A % Reshift (shift) ! A = A + shift * I
         this % Ashift = shift
      ENDIF
#else
      REAL*8,                     INTENT(IN)        :: dt
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE ReSetOperatorDt
!
!/////////////////////////////////////////////////////////////////////////////////////////////////   
!
   SUBROUTINE SetBValues(this, nvalues, irow, values)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      PetscInt,                        INTENT(IN)        :: nvalues
      PetscInt, DIMENSION(:),          INTENT(IN)        :: irow
      PetscScalar, DIMENSION(:),       INTENT(IN)        :: values
      PetscErrorCode                                     :: ierr
       
      CALL VecSetValues(this%b,nvalues, irow,values,INSERT_VALUES, ierr)
      CALL CheckPetscErr(ierr, 'error in VecSetValues')
#else
      INTEGER,                        INTENT(IN)        :: nvalues
      INTEGER, DIMENSION(:),          INTENT(IN)        :: irow
      REAL*8    , DIMENSION(:),       INTENT(IN)        :: values
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SetBValues
!/////////////////////////////////////////////////////////////////////////////////////////////////   
   SUBROUTINE SetBValue(this, irow, value)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      PetscInt,                        INTENT(IN)        :: irow
      PetscScalar,                     INTENT(IN)        :: value
      PetscErrorCode                                     :: ierr
        
      CALL VecSetValue(this%b, irow,value,INSERT_VALUES, ierr)
      CALL CheckPetscErr(ierr, 'error in VecSetValues')
#else
      INTEGER,           INTENT(IN)        :: irow
      REAL*8    ,        INTENT(IN)        :: value
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SetBValue
!
!/////////////////////////////////////////////////////////////////////////////////////////////////     
!
   SUBROUTINE AssemblyB(this)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)   :: this
#ifdef HAS_PETSC
      PetscErrorCode                                     :: ierr
      
      CALL VecAssemblyBegin(this%b, ierr);  CALL CheckPetscErr(ierr," Assembly B in PETSc Begin")      
      CALL VecAssemblyEnd(this%b, ierr)  ;  CALL CheckPetscErr(ierr," Assembly B in PETSc End")  
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE GetXValues(this, nvalues, irow, values)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)      :: this
#ifdef HAS_PETSC
      PetscInt,                        INTENT(IN)         :: nvalues
      PetscInt, DIMENSION(:),          INTENT(IN)         :: irow
      PetscScalar, DIMENSION(:),       INTENT(OUT)        :: values
      PetscErrorCode                                      :: ierr
      
      CALL VecGetValues(this%x,nvalues,irow,values, ierr)
      CALL CheckPetscErr(ierr, 'error in VecGetValues')
#else
      INTEGER,                        INTENT(IN)        :: nvalues
      INTEGER, DIMENSION(:),          INTENT(IN)        :: irow
      REAL*8    , DIMENSION(:),       INTENT(IN)        :: values
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE GetXValues
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE GetXValue(this, irow, x_i)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t),     INTENT(INOUT)      :: this
#ifdef HAS_PETSC
      PetscInt,                        INTENT(IN)         :: irow
      PetscScalar,                     INTENT(OUT)        :: x_i
      PetscErrorCode                                      :: ierr
      
      CALL VecGetValues(this%x,1 ,irow,x_i, ierr)
      CALL CheckPetscErr(ierr, 'error in VecGetValue')
#else
      INTEGER,           INTENT(IN)        :: irow
      REAL*8    ,        INTENT(OUT)       :: x_i
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE GetXValue
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE SaveMat(this,filename)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT)         :: this
      CHARACTER(LEN=*), OPTIONAL                         :: filename
#ifdef HAS_PETSC
      PetscViewer                                        :: viewer
      PetscErrorCode                                     :: ierr
      
      
      CALL MatView(this % A % A,PETSC_VIEWER_DRAW_SELF)
      read(*,*)
!~       IF (.NOT. PRESENT(filename)) filename = &
!~                             '/home/andresrueda/Dropbox/PhD/03_Initial_Codes/3D/Implicit/nslite3d/Tests/Euler/NumJac/MatMatlab.dat'
!~       CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename , viewer, ierr)    ; CALL CheckPetscErr(ierr)
!~       CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB , ierr)      ; CALL CheckPetscErr(ierr)
!~       CALL MatView(this%A, viewer, ierr)                                      ; CALL CheckPetscErr(ierr)
!~       CALL PetscViewerDestroy(viewer, ierr)                                   ; CALL CheckPetscErr(ierr)
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SaveMat
!////////////////////////////////////////////////////////////////////////////////////////////////// 
   SUBROUTINE DestroyPetscObjects(this)
      IMPLICIT NONE
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT)       :: this
#ifdef HAS_PETSC
      PetscErrorCode                                   :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6
      CALL VecDestroy(this%x,ierr1)
      CALL VecDestroy(this%b,ierr2)
      call this % A % destruct
!~       CALL PCDestroy (this%pc, ierr5) !! Outcommented: PC destroyed in KSP environment (KSPDestroy)
      CALL KSPDestroy(this%ksp,ierr4)
      
      CALL PetscFinalize(ierr6)
      CALL CheckPetscErr(ierr1,'error in VecDestroy x')
      CALL CheckPetscErr(ierr2,'error in VecDestroy b')
      CALL CheckPetscErr(ierr3,'error in MatDestroy_A')
!~       CALL CheckPetscErr(ierr5,'error in PCDestroy')
      CALL CheckPetscErr(ierr4,'error in KSPDestroy')
      CALL CheckPetscErr(ierr6,'error in PetscFinalize')
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE DestroyPetscObjects
!
!////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      !--------------------------------------------------------------
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                             :: TypeOfNorm
      REAL(KIND=RP)                                :: xnorm
      !--------------------------------------------------------------
#ifdef HAS_PETSC
      PetscErrorCode                               :: ierr
      !--------------------------------------------------------------
      
      SELECT CASE(TypeOfNorm)
         CASE('infinity')
            CALL VecNorm(this%x,NORM_INFINITY,xnorm,ierr)       ; CALL CheckPetscErr(ierr,'error in VecNorm')
         CASE DEFAULT
            STOP 'PetscSolverClass error: Type of Norm not defined'
      END SELECT
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END FUNCTION Getxnorm
!
!////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   FUNCTION Getrnorm(this) RESULT(rnorm)
      IMPLICIT NONE
      !--------------------------------------------------------------
      CLASS(PetscKspLinearSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                                :: rnorm
      !--------------------------------------------------------------
#ifdef HAS_PETSC
      PetscErrorCode                               :: ierr
      !--------------------------------------------------------------
      
      ! I don't know which type of norm PETSc computes!
      CALL KSPGetResidualNorm(this%ksp, rnorm, ierr)   ; CALL CheckPetscErr(ierr,'error in KSPGetResidualNorm')
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END FUNCTION Getrnorm
END MODULE PetscSolverClass
