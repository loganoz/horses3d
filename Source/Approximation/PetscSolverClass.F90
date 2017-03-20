MODULE PetscSolverClass

   IMPLICIT NONE

#include <petsc.h>

   TYPE PetscKspLinearSolver
      Mat                                           :: A
      Vec                                           :: x
      Vec                                           :: b
      KSP                                           :: ksp
      PC                                            :: pc
      PetscInt                                      :: dimprb
      PetscReal                                     :: tol = 1d-15
      PetscInt                                      :: maxiter = 50
      PetscInt                                      :: nz = 0
      PetscScalar                                   :: Ashift = 0.0d0
      PetscReal                                     :: residual
      PetscReal                                     :: xnorm
      PetscInt                                      :: niter
      PetscBool                                     :: converged = PETSC_FALSE
      PetscBool                                     :: init_context = PETSC_FALSE
      PetscBool                                     :: setpreco = PETSC_TRUE
      PetscBool                                     :: withMPI = PETSC_FALSE
      CONTAINS
         PROCEDURE                                  :: construct => ConstructPetscContext
         PROCEDURE                                  :: SetAColumn
         PROCEDURE                                  :: SetBValues
         PROCEDURE                                  :: SetBValue
         PROCEDURE                                  :: GetXValues
         PROCEDURE                                  :: GetXValue
         PROCEDURE                                  :: PreallocateA
         PROCEDURE                                  :: SetOperatorDt
         PROCEDURE                                  :: AssemblyA
         PROCEDURE                                  :: AssemblyB
         PROCEDURE                                  :: ResetA
         PROCEDURE                                  :: SaveMat
         PROCEDURE                                  :: solve   =>   SolveLinPrb
         PROCEDURE                                  :: clean   =>   DestroyPetscObjects
   END TYPE PetscKspLinearSolver
   
  
   PRIVATE                                          
   PUBLIC                                           :: PetscKspLinearSolver
   PUBLIC                                           :: ConstructPetscContext, SolveLinPrb, PreallocateA, SetAColumn, SetBValue
   PUBLIC                                           :: SetBValues, DestroyPetscObjects, ShiftA, AssemblyA, AssemblyB, GetXValues
   PUBLIC                                           :: GetXValue, ResetA, SaveMat, SetOperatorDt
   
   CONTAINS  
!///////////////////////////////////////////////////////////////////////////////
   SUBROUTINE CheckPetscErr(ierr,msg)
      PetscErrorCode, INTENT(IN)                   :: ierr
      CHARACTER(LEN=*), OPTIONAL                   :: msg
      IF(ierr .EQ. 0) THEN
         RETURN
      ELSE
         IF (.NOT. PRESENT(msg)) msg = 'error in petsc'
         WRITE(*,*) msg,' **** Petsc call returned an error. Code: ' ,ierr
         STOP
      ENDIF
   END SUBROUTINE CheckPetscErr
!///////////////////////////////////////////////////////////////////////////////     
   SUBROUTINE ConstructPetscContext(this, DimPrb)
      CLASS(PetscKspLinearSolver), INTENT(INOUT) :: this
      PetscInt, INTENT(IN)                       :: DimPrb
      PetscErrorCode                             :: ierr
      
      !Initialisation of the PETSc variables
      CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)

!     PETSc matrix A 
      CALL MatCreate(PETSC_COMM_WORLD,this%A,ierr)                           ; CALL CheckPetscErr(ierr,'error creating A matrix')
      IF (this%withMPI) THEN
         CALL MatSetSizes(this%A,PETSC_DECIDE,PETSC_DECIDE,dimPrb,dimPrb,ierr)
         CALL CheckPetscErr(ierr,'error setting mat size')
         CALL MatSetType(this%A,MATMPIAIJ, ierr)                              
         CALL CheckPetscErr(ierr,'error in MatSetType')
         CALL MatSetFromOptions(this%A,ierr)                                  
         CALL CheckPetscErr(ierr,'error in MatSetFromOptions')
      ELSE
         CALL MatSetSizes(this%A,PETSC_DECIDE,PETSC_DECIDE,dimPrb,dimPrb,ierr)
         CALL CheckPetscErr(ierr,'error setting mat size')
         CALL CheckPetscErr(ierr,'error in MatSetType')
         CALL MatSetFromOptions(this%A,ierr)                                  
         CALL CheckPetscErr(ierr,'error in MatSetFromOptions')
!         CALL MatSetOption(this%A,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)  
         CALL MatSetOption(this%A,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE,ierr)                 
         CALL CheckPetscErr(ierr,'error in MatSetOption')
       ENDIF

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
      CALL KSPSetOperators(this%ksp, this%A, this%A, ierr)              ; CALL CheckPetscErr(ierr,'error in KSPSetOperators')
!~      CALL KSPSetTolerances(this%ksp,PETSC_DEFAULT_REAL,this%tol,PETSC_DEFAULT_REAL,this%maxiter,ierr)
!~      CALL CheckPetscErr(ierr,'error in KSPSetTolerances')
      this%init_context = PETSC_TRUE
      this%dimprb = DimPrb
      
   END SUBROUTINE
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE SetPetscPreco(this)
      CLASS(PetscKspLinearSolver), INTENT(INOUT)      :: this
      PetscErrorCode                                  :: ierr
      
      CALL PCSetType(this%pc,PCILU,ierr)       ;CALL CheckPetscErr(ierr, 'error in PCSetType')
      this%setpreco = PETSC_FALSE

   END SUBROUTINE SetPetscPreco 
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE SolveLinPrb(this, tol, maxiter)
      CLASS(PetscKspLinearSolver), INTENT(INOUT)      :: this
      PetscReal, OPTIONAL                             :: tol
      PetscInt,  OPTIONAL                             :: maxiter
      PetscErrorCode                                  :: ierr
      
      ! Set , if given, solver tolerance and max number of iterations
      IF (PRESENT(tol)) THEN
         this%tol = tol
      ENDIF
      IF (PRESENT(maxiter)) THEN
         this%maxiter = maxiter
      ENDIF
      
      CALL KSPSetTolerances(this%ksp,this%tol,this%tol,PETSC_DEFAULT_REAL,this%maxiter,ierr)
      CALL CheckPetscErr(ierr, 'error in KSPSetTolerances')
      
      IF (this%setpreco) THEN    ! if setpreco = true sets the default preconditioner PCJACOBI. If PC already set, setpreco = false
         CALL PCSetType(this%pc,PCJACOBI,ierr)                 ; CALL CheckPetscErr(ierr, 'error in PCSetType')
         this%setpreco = PETSC_FALSE
         !Summary (arueda):
         !  - PCILU: Bad results
         !  - PCJACOBI: Slow convergence
         !  - PCPBJACOBI: (point Jacobi) Almost same convergence as PCJACOBI
         !  - PCBJACOBI: (block Jacobi) Bad results 
      ENDIF

      CALL KSPSetOperators(this%ksp, this%A, this%A, ierr)     ; CALL CheckPetscErr(ierr, 'error in KSPSetOperators')
      CALL KSPSolve(this%ksp,this%b,this%x,ierr)               ; CALL CheckPetscErr(ierr, 'error in KSPSolve')
      CALL KSPGetIterationNumber(this%ksp,this%niter,ierr)     ; CALL CheckPetscErr(ierr,'error in KSPGetIterationNumber')
      CALL KSPGetResidualNorm(this%ksp, this%residual, ierr)   ; CALL CheckPetscErr(ierr,'error in KSPGetResidualNorm')
      CALL VecNorm(this%x,NORM_INFINITY,this%xnorm,ierr)       ; CALL CheckPetscErr(ierr,'error in VecNorm')
      IF (this%niter < maxiter) THEN
         this%converged = PETSC_TRUE
      ELSE
         this%converged = PETSC_FALSE
      END IF
      
   END SUBROUTINE SolveLinPrb
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE PreallocateA(this, nnz)
      CLASS(PetscKspLinearSolver), INTENT(INOUT)       :: this
      PetscInt                                         :: nnz
      PetscErrorCode                                   :: ierr      
      
      this%nz = nnz
      IF(this%withMPI) THEN
         CALL MatMPIAIJSetPreallocation(this%A, this%nz/5, PETSC_NULL_INTEGER,4 * this%nz/5, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatMPIAIJSetPreallocation')
      ELSE
         CALL MatSeqAIJSetPreallocation(this%A, this%nz, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatSeqAIJSetPreallocation')
      ENDIF
         
   END SUBROUTINE
!/////////////////////////////////////////////////////////////////////////////////////////////////   
   SUBROUTINE SetAColumn(this,nvalues, irow, icol, values )
      CLASS(PetscKspLinearSolver), INTENT(INOUT)         :: this
      PetscInt, INTENT(IN)                               :: nvalues
      PetscInt, DIMENSION(:), INTENT(IN)                 :: irow
      PetscInt, INTENT(IN)                               :: icol
      PetscScalar, DIMENSION(:), INTENT(IN)              :: values
      PetscErrorCode                                     :: ierr
   
      CALL MatSetValues(this%A,nvalues,irow,1,icol,values,INSERT_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
   
   END SUBROUTINE SetAColumn
!/////////////////////////////////////////////////////////////////////////////////////////////////   
   SUBROUTINE ShiftA(this, shift)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscScalar                                        :: shift
      PetscErrorCode                                     :: ierr
      
      CALL MatShift(this%A,shift, ierr)                       ! A = A + shift * I
      this%AShift = this%AShift + shift                       ! Updates AShift value in linsolver structure
      CALL CheckPetscErr(ierr)

   END SUBROUTINE ShiftA
!/////////////////////////////////////////////////////////////////////////////////////////////////      
   SUBROUTINE ResetA(this)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscErrorCode                                     :: ierr

      CALL MatZeroEntries(this%A, ierr)
      CALL CheckPetscErr(ierr,'error in MatZeroEntries')
      
   END SUBROUTINE
!/////////////////////////////////////////////////////////////////////////////////////////////////   
   SUBROUTINE SetOperatorDt(this, dt, coeff)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscScalar,                     INTENT(IN)        :: dt
      PetscScalar,                     OPTIONAL          :: coeff
      PetscScalar                                        :: shift
      PetscScalar                                        :: eps = 1e-10                                       
      PetscErrorCode                                     :: ierr

      shift = -1/dt - this%AShift
      IF (ABS(shift) .GT. eps) THEN
         CALL MatShift(this%A,shift, ierr)                  ! A = A + shift * I
         CALL CheckPetscErr(ierr)
         this%Ashift = shift
      ENDIF

   END SUBROUTINE SetOperatorDt
!/////////////////////////////////////////////////////////////////////////////////////////////////   
   SUBROUTINE SetBValues(this, nvalues, irow, values)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscInt,                        INTENT(IN)        :: nvalues
      PetscInt, DIMENSION(:),          INTENT(IN)        :: irow
      PetscScalar, DIMENSION(:),       INTENT(IN)        :: values
      PetscErrorCode                                     :: ierr
       
      CALL VecSetValues(this%b,nvalues, irow,values,INSERT_VALUES, ierr)
      CALL CheckPetscErr(ierr, 'error in VecSetValues')
    
   END SUBROUTINE SetBValues
!/////////////////////////////////////////////////////////////////////////////////////////////////   
   SUBROUTINE SetBValue(this, irow, value)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscInt,                        INTENT(IN)        :: irow
      PetscScalar,                     INTENT(IN)        :: value
      PetscErrorCode                                     :: ierr
        
      CALL VecSetValue(this%b, irow,value,INSERT_VALUES, ierr)
      CALL CheckPetscErr(ierr, 'error in VecSetValues')
    
   END SUBROUTINE SetBValue
!/////////////////////////////////////////////////////////////////////////////////////////////////     
   SUBROUTINE AssemblyA(this)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscErrorCode                                     :: ierr
 
      CALL MatAssemblyBegin(this%A,MAT_FINAL_ASSEMBLY,ierr);  CALL CheckPetscErr(ierr," Assembly A in PETSc Begin")      
      CALL MatAssemblyEnd(this%A,MAT_FINAL_ASSEMBLY,ierr)  ;  CALL CheckPetscErr(ierr," Assembly A in PETSc End")                  
  
   END SUBROUTINE
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE AssemblyB(this)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)     :: this
      PetscErrorCode                                     :: ierr
      
      CALL VecAssemblyBegin(this%b, ierr);  CALL CheckPetscErr(ierr," Assembly B in PETSc Begin")      
      CALL VecAssemblyEnd(this%b, ierr)  ;  CALL CheckPetscErr(ierr," Assembly B in PETSc End")                  
   END SUBROUTINE
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE GetXValues(this, nvalues, irow, values)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)      :: this
      PetscInt,                        INTENT(IN)         :: nvalues
      PetscInt, DIMENSION(:),          INTENT(IN)         :: irow
      PetscScalar, DIMENSION(:),       INTENT(OUT)        :: values
      PetscErrorCode                                      :: ierr
      
      CALL VecGetValues(this%x,nvalues,irow,values, ierr)
      CALL CheckPetscErr(ierr, 'error in VecGetValues')
      
   END SUBROUTINE GetXValues
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE GetXValue(this, irow, value)
      CLASS(PetscKspLinearSolver),     INTENT(INOUT)      :: this
      PetscInt,                        INTENT(IN)         :: irow
      PetscScalar,                     INTENT(OUT)        :: value
      PetscErrorCode                                      :: ierr
      
      CALL VecGetValues(this%x,1 ,irow,value, ierr)
      CALL CheckPetscErr(ierr, 'error in VecGetValue')
      
   END SUBROUTINE GetXValue
!//////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE SaveMat(this,filename)
      CLASS(PetscKspLinearSolver), INTENT(INOUT)         :: this
      CHARACTER(LEN=*), OPTIONAL                         :: filename
      
      PetscViewer                                        :: viewer
      PetscErrorCode                                     :: ierr
      
      
!~       CALL MatView(this%A,PETSC_VIEWER_DRAW_SELF)
!~       read(*,*)
      IF (.NOT. PRESENT(filename)) filename = &
                            '/home/andresrueda/Dropbox/PhD/03_Initial_Codes/3D/Implicit/nslite3d/Tests/Euler/NumJac/MatMatlab.dat'
      CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename , viewer, ierr)    ; CALL CheckPetscErr(ierr)
      CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB , ierr)      ; CALL CheckPetscErr(ierr)
      CALL MatView(this%A, viewer, ierr)                                      ; CALL CheckPetscErr(ierr)
      CALL PetscViewerDestroy(viewer, ierr)                                   ; CALL CheckPetscErr(ierr)
   END SUBROUTINE SaveMat
!////////////////////////////////////////////////////////////////////////////////////////////////// 
   SUBROUTINE DestroyPetscObjects(this)
      CLASS(PetscKspLinearSolver), INTENT(INOUT)       :: this
      PetscErrorCode                                   :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6

      CALL VecDestroy(this%x,ierr1)
      CALL VecDestroy(this%b,ierr2)
      CALL MatDestroy(this%A,ierr3)
      CALL KSPDestroy(this%ksp,ierr4)
      CALL PCDestroy(this%pc, ierr5)
      CALL PetscFinalize(ierr6)
      CALL CheckPetscErr(ierr1,'error in VecDestroy x')
      CALL CheckPetscErr(ierr2,'error in VecDestroy b')
      CALL CheckPetscErr(ierr3,'error in MatDestroy_A')
      CALL CheckPetscErr(ierr4,'error in KSPDestroy')
      CALL CheckPetscErr(ierr5,'error in PCDestroy')
      CALL CheckPetscErr(ierr6,'error in PetscFinalize')
   END SUBROUTINE DestroyPetscObjects
!////////////////////////////////////////////////////////////////////////////////////////////////// 
END MODULE PetscSolverClass
