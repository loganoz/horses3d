!////////////////////////////////////////////////////////////////////////
!
!      PETScMatrixClass.f90
!      Created: 2018-02-18 17:07:00 +0100 
!      By: Andrés Rueda
!
!      Class for sparse csr matrices in PETSc context
!
!////////////////////////////////////////////////////////////////////////
module PETScMatrixClass
   use SMConstants
   use GenericMatrixClass
   use CSRMatrixClass
   implicit none
#ifdef HAS_PETSC
#include <petsc.h>
#endif
   private
   public PETSCMatrix_t, Matrix_t
   
   type, extends(Matrix_t) :: PETSCMatrix_t
#ifdef HAS_PETSC
      Mat         :: A        ! Matrix in PETSc context 
      PetscScalar :: Ashift   ! Stores the current shift to the matrix
      PetscBool   :: withMPI
#endif
      contains
         procedure :: construct
         procedure :: destruct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetColumn
         procedure :: AddToColumn
         procedure :: Shift
         procedure :: ReShift
         procedure :: PreAssembly
         procedure :: Assembly
         procedure :: GetCSRMatrix
   end type PETSCMatrix_t
   
!
!  Module variables
!  ----------------
#ifdef HAS_PETSC
   PetscErrorCode :: ierr
#endif
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,dimPrb,withMPI)
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t)  :: this
#ifdef HAS_PETSC
      PetscInt, INTENT(IN)  :: DimPrb
      PetscBool, optional, INTENT(IN) :: withMPI
      !---------------------------------------------
      
      this % NumRows = dimPrb
      this % withMPI = withMPI
      
      !     PETSc matrix A 
      CALL MatCreate(PETSC_COMM_WORLD,this%A,ierr)                           ; CALL CheckPetscErr(ierr,'error creating A matrix')
      
      IF (withMPI) THEN
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
#else
      integer, INTENT(IN)              :: dimPrb
      logical, optional, INTENT(IN)    :: WithMPI
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Preallocate(this, nnz, nnzs)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t), INTENT(INOUT)       :: this
#ifdef HAS_PETSC
      PetscInt, optional, intent(in)            :: nnz
      INTEGER, optional, intent(in)  :: nnzs(:)
      !---------------------------------------------
      
      if (.not. present(nnz)) ERROR stop ':: PETSc matrix needs nnz'
      
      IF(this%withMPI) THEN
         CALL MatMPIAIJSetPreallocation(this%A, nnz/25, PETSC_NULL_INTEGER,24 * nnz/25, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatMPIAIJSetPreallocation')
      ELSE
         CALL MatSeqAIJSetPreallocation(this%A, nnz, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatSeqAIJSetPreallocation')
      ENDIF
#else
      INTEGER, optional, intent(in)  :: nnz
      INTEGER, optional, intent(in)  :: nnzs(:)
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Reset(this)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      !---------------------------------------------
      integer :: i
      
      CALL MatZeroEntries(this%A, ierr)
      CALL CheckPetscErr(ierr,'error in MatZeroEntries')
      
      ! secure diagonal entries
      do i=0, this % NumRows-1
         CALL MatSetValues(this%A,1,i,1,i,0._RP ,INSERT_VALUES,ierr)
      end do
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetColumn(this,nvalues, irow, icol, values )
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t), INTENT(INOUT)         :: this
#ifdef HAS_PETSC
      PetscInt, INTENT(IN)                               :: nvalues
      PetscInt, DIMENSION(:), INTENT(IN)                 :: irow
      PetscInt, INTENT(IN)                               :: icol
      PetscScalar, DIMENSION(:), INTENT(IN)              :: values
      !---------------------------------------------
   
      CALL MatSetValues(this%A,nvalues,irow-1,1,icol-1,values,INSERT_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
#else
      INTEGER, INTENT(IN)                               :: nvalues
      INTEGER, DIMENSION(:), INTENT(IN)                 :: irow
      INTEGER, INTENT(IN)                               :: icol
      REAL*8 , DIMENSION(:), INTENT(IN)                 :: values
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   SUBROUTINE AddToColumn(this,nvalues, irow, icol, values )
      IMPLICIT NONE
      CLASS(PETSCMatrix_t), INTENT(INOUT)         :: this
#ifdef HAS_PETSC
      PetscInt, INTENT(IN)                               :: nvalues
      PetscInt, DIMENSION(:), INTENT(IN)                 :: irow
      PetscInt, INTENT(IN)                               :: icol
      PetscScalar, DIMENSION(:), INTENT(IN)              :: values
   
      CALL MatSetValues(this%A,nvalues,irow-1,1,icol-1,values,ADD_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
#else
      INTEGER, INTENT(IN)                               :: nvalues
      INTEGER, DIMENSION(:), INTENT(IN)                 :: irow
      INTEGER, INTENT(IN)                               :: icol
      REAL*8 , DIMENSION(:), INTENT(IN)                 :: values
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE AddToColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Shift(this, shiftval)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t), INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      PetscScalar,                   INTENT(IN)        :: shiftval
      !---------------------------------------------
      
      CALL MatShift(this%A,shiftval, ierr)                  ! A = A + shift * I
      CALL CheckPetscErr(ierr)
      this % Ashift = shiftval
      
#else
      REAL*8,                     INTENT(IN)        :: shiftval
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE Shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Removes previous shift in order to insert new one 
!              (important when Jacobian is reused)
!  --------------------------------------------------
   SUBROUTINE ReShift(this, shiftval)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t), INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      PetscScalar,                   INTENT(IN)        :: shiftval
      !---------------------------------------------
      
      CALL MatShift(this%A,shiftval - this % Ashift, ierr)                  ! A = A + shift * I
      CALL CheckPetscErr(ierr)
      this % Ashift = shiftval
#else
      REAL*8,                     INTENT(IN)        :: shiftval
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE ReShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE PreAssembly(this)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     INTENT(INOUT)   :: this
#ifdef HAS_PETSC
      !---------------------------------------------
      
      CALL MatAssemblyBegin(this%A,MAT_FLUSH_ASSEMBLY ,ierr);  CALL CheckPetscErr(ierr," PreAssembly A in PETSc Begin")      
      CALL MatAssemblyEnd(this%A,MAT_FLUSH_ASSEMBLY,ierr)  ;  CALL CheckPetscErr(ierr," PreAssembly A in PETSc End")                  
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE PreAssembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Assembly(this,BlockIdx,BlockSize)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     INTENT(INOUT)   :: this
      INTEGER, TARGET, OPTIONAL    ,     INTENT(IN)      :: BlockIdx(:)
      INTEGER, TARGET, OPTIONAL, INTENT(IN)    :: BlockSize(:)
      !---------------------------------------------
#ifdef HAS_PETSC
 
      CALL MatAssemblyBegin(this%A,MAT_FINAL_ASSEMBLY,ierr);  CALL CheckPetscErr(ierr," Assembly A in PETSc Begin")      
      CALL MatAssemblyEnd(this%A,MAT_FINAL_ASSEMBLY,ierr)  ;  CALL CheckPetscErr(ierr," Assembly A in PETSc End")                  
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE Assembly
!
!/////////////////////////////////////////////////////////////////////////////////////////////////     
!
!
!////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   SUBROUTINE GetCSRMatrix(this,Acsr)
      IMPLICIT NONE
      !---------------------------------------------------------------------------------
      CLASS(PETSCMatrix_t), INTENT(IN)  :: this        !<    PETSc matrix
      TYPE(csrMat_t)      , INTENT(OUT) :: Acsr        !>    CSR matrix
      !---------------------------------------------------------------------------------
#ifdef HAS_PETSC
      INTEGER                                  :: i, ncols
      INTEGER                                  :: nnz_row(this % NumRows)
      INTEGER      , ALLOCATABLE, DIMENSION(:) :: ACols 
      REAL(KIND=RP), ALLOCATABLE, DIMENSION(:) :: AVals
      PetscErrorCode                           :: ierr
      !---------------------------------------------------------------------------------
      
      !We first get the number of nonzero entries in each row
      DO i = 1, this % NumRows
         CALL MatGetRow(this % A,i-1,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
         CALL CheckPetscErr(ierr,'error in Petsc MatGetRow')
         
         nnz_row(i) = ncols
         
         CALL MatRestoreRow(this % A,i-1,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
         CALL CheckPetscErr(ierr,'error in Petsc MatRestoreRow')
      END DO
      
      CALL Acsr % construct(this % NumRows)
      call Acsr % Preallocate(nnzs= nnz_row)
      call Acsr % Reset
      
      DO i = 1, this % NumRows
         
         ALLOCATE(AVals(nnz_row(i)))
         ALLOCATE(ACols(nnz_row(i)))
         
         CALL MatGetRow(this % A,i-1,ncols,ACols,AVals,ierr)      ;  CALL CheckPetscErr(ierr,'error in Petsc MatGetRow')
         
         Acsr % Values  (Acsr % Rows (i) : Acsr % Rows (i+1) -1) = AVals
         Acsr % Cols    (Acsr % Rows (i) : Acsr % Rows (i+1) -1) = ACols + 1
         
         CALL MatRestoreRow(this % A,i-1,ncols,ACols,AVals,ierr)  ;  CALL CheckPetscErr(ierr,'error in Petsc MatRestoreRow')
         
         DEALLOCATE(AVals)
         DEALLOCATE(ACols)
      END DO
      
      CALL Acsr % assigndiag
      
#else
      STOP ':: PETSc is not linked correctly'
#endif
   END SUBROUTINE GetCSRMatrix
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destruct(this)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     INTENT(INOUT)     :: this
#ifdef HAS_PETSC
      CALL MatDestroy(this%A,ierr)
      CALL CheckPetscErr(ierr," A destruction")  
#endif
   END SUBROUTINE destruct
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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
   
end module PETScMatrixClass
