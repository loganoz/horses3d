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
#ifdef HAS_PETSC
   use petsc
#endif
   implicit none
   private
   public PETSCMatrix_t, Matrix_t

#ifdef HAS_PETSC
!#include <petsc.h>
#include "petsc/finclude/petsc.h"
#endif
   
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
         procedure :: SetEntry
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
!
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
      PetscInt, intent(in)  :: DimPrb
      PetscBool, optional, intent(in) :: withMPI
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
      integer, intent(in)              :: dimPrb
      logical, optional, intent(in)    :: WithMPI
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs, ForceDiagonal)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t), intent(inout)       :: this
      integer, optional, intent(in)  :: nnzs(:)
      logical, optional, intent(in)  :: ForceDiagonal
#ifdef HAS_PETSC
      PetscInt, optional, intent(in)            :: nnz
      !---------------------------------------------
      logical :: mustForceDiagonal
      !---------------------------------------------
      
      if (.not. present(nnz)) ERROR stop ':: PETSc matrix needs nnz'
      
      if ( present(ForceDiagonal) ) then
         mustForceDiagonal = ForceDiagonal
      else
         mustForceDiagonal = .FALSE.
      end if
      
      IF(this%withMPI) THEN
         CALL MatMPIAIJSetPreallocation(this%A, nnz/25, PETSC_NULL_INTEGER,24 * nnz/25, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatMPIAIJSetPreallocation')
      ELSE
         CALL MatSeqAIJSetPreallocation(this%A, nnz, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatSeqAIJSetPreallocation')
      ENDIF
      
      if (mustForceDiagonal) then
         do i = 1, this % NumRows
            call this % SetEntry(i,i,0._RP)
         end do
      end if
#else
      INTEGER, optional, intent(in)  :: nnz
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     intent(inout)     :: this
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
   end subroutine Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetEntry(this, row, col, value )
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t), intent(inout) :: this
#ifdef HAS_PETSC
      PetscInt            , intent(in)    :: row
      PetscInt            , intent(in)    :: col
      PetscScalar         , intent(in)    :: value
      !---------------------------------------------
   
      CALL MatSetValues(this%A, 1 ,row-1,1,col-1,value,INSERT_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
#else
      INTEGER        , intent(in) :: row
      INTEGER        , intent(in) :: col
      real(kind=RP)  , intent(in) :: value
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine SetEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetColumn(this,nvalues, irow, icol, values )
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t), intent(inout)         :: this
#ifdef HAS_PETSC
      PetscInt, intent(in)                               :: nvalues
      PetscInt, DIMENSION(:), intent(in)                 :: irow
      PetscInt, intent(in)                               :: icol
      PetscScalar, DIMENSION(:), intent(in)              :: values
      !---------------------------------------------
   
      CALL MatSetValues(this%A,nvalues,irow-1,1,icol-1,values,INSERT_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
#else
      INTEGER, intent(in)                               :: nvalues
      INTEGER, DIMENSION(:), intent(in)                 :: irow
      INTEGER, intent(in)                               :: icol
      REAL*8 , DIMENSION(:), intent(in)                 :: values
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   subroutine AddToColumn(this,nvalues, irow, icol, values )
      implicit none
      CLASS(PETSCMatrix_t), intent(inout)         :: this
#ifdef HAS_PETSC
      PetscInt, intent(in)                               :: nvalues
      PetscInt, DIMENSION(:), intent(in)                 :: irow
      PetscInt, intent(in)                               :: icol
      PetscScalar, DIMENSION(:), intent(in)              :: values
   
      CALL MatSetValues(this%A,nvalues,irow-1,1,icol-1,values,ADD_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
#else
      INTEGER, intent(in)                               :: nvalues
      INTEGER, DIMENSION(:), intent(in)                 :: irow
      INTEGER, intent(in)                               :: icol
      REAL*8 , DIMENSION(:), intent(in)                 :: values
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine AddToColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Shift(this, shiftval)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t), intent(inout)     :: this
#ifdef HAS_PETSC
      PetscScalar,                   intent(in)        :: shiftval
      !---------------------------------------------
      
      CALL MatShift(this%A,shiftval, ierr)                  ! A = A + shift * I
      CALL CheckPetscErr(ierr)
      this % Ashift = shiftval
      
#else
      REAL*8,                     intent(in)        :: shiftval
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine Shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Removes previous shift in order to insert new one 
!              (important when Jacobian is reused)
!  --------------------------------------------------
   subroutine ReShift(this, shiftval)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t), intent(inout)     :: this
#ifdef HAS_PETSC
      PetscScalar,                   intent(in)        :: shiftval
      !---------------------------------------------
      
      CALL MatShift(this%A,shiftval - this % Ashift, ierr)                  ! A = A + shift * I
      CALL CheckPetscErr(ierr)
      this % Ashift = shiftval
#else
      REAL*8,                     intent(in)        :: shiftval
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine ReShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PreAssembly(this)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     intent(inout)   :: this
#ifdef HAS_PETSC
      !---------------------------------------------
      
      CALL MatAssemblyBegin(this%A,MAT_FLUSH_ASSEMBLY ,ierr);  CALL CheckPetscErr(ierr," PreAssembly A in PETSc Begin")      
      CALL MatAssemblyEnd(this%A,MAT_FLUSH_ASSEMBLY,ierr)  ;  CALL CheckPetscErr(ierr," PreAssembly A in PETSc End")                  
#else
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine PreAssembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Assembly(this)
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t), intent(inout)   :: this
      !---------------------------------------------
#ifdef HAS_PETSC
 
      CALL MatAssemblyBegin(this%A,MAT_FINAL_ASSEMBLY,ierr);  CALL CheckPetscErr(ierr," Assembly A in PETSc Begin")      
      CALL MatAssemblyEnd(this%A,MAT_FINAL_ASSEMBLY,ierr)  ;  CALL CheckPetscErr(ierr," Assembly A in PETSc End")                  
#else
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine Assembly
!
!/////////////////////////////////////////////////////////////////////////////////////////////////     
!
!
!////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   subroutine GetCSRMatrix(this,Acsr)
      implicit none
      !---------------------------------------------------------------------------------
      CLASS(PETSCMatrix_t), intent(in)  :: this        !<    PETSc matrix
      TYPE(csrMat_t)      , intent(OUT) :: Acsr        !>    CSR matrix
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
   end subroutine GetCSRMatrix
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     intent(inout)     :: this
#ifdef HAS_PETSC
      CALL MatDestroy(this%A,ierr)
      CALL CheckPetscErr(ierr," A destruction")  
#endif
   end subroutine destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CheckPetscErr(ierr,msg)
      implicit none
      CHARACTER(LEN=*), OPTIONAL                   :: msg
#ifdef HAS_PETSC
      PetscErrorCode, intent(in)                   :: ierr
      
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
   end subroutine CheckPetscErr
   
end module PETScMatrixClass
