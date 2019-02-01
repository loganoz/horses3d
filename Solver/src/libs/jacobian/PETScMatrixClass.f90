!
!//////////////////////////////////////////////////////
!
!   @File:    PETScMatrixClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Sun Feb 18 14:00:00 2018
!   @Last revision date: Fri Feb  1 17:25:00 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 0bf6bde04abec1f8f9eb04f644c9cac0cc0df9e9
!
!//////////////////////////////////////////////////////
!
!      Class for sparse csr matrices in PETSc context
!
!////////////////////////////////////////////////////////////////////////
#include "Includes.h"
#ifdef HAS_PETSC
#include "petsc/finclude/petsc.h"
#endif
module PETScMatrixClass
   use SMConstants
   use GenericMatrixClass
   use CSRMatrixClass
   use Jacobian            , only: JACEPS
#ifdef HAS_PETSC
   use petsc
#endif
   implicit none
   private
   public PETSCMatrix_t, Matrix_t
   
   type, extends(Matrix_t) :: PETSCMatrix_t
#ifdef HAS_PETSC
      Mat         :: A        ! Matrix in PETSc context 
      PetscScalar :: Ashift   ! Stores the current shift to the matrix
      PetscBool   :: withMPI
#endif
      ! Variables for matrices with blocks
      integer      ,  allocatable :: BlockIdx(:)  ! Index of first element of block
      integer      ,  allocatable :: BlockSize(:) ! Size of each block
      
      contains
         procedure :: construct
         procedure :: destruct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetEntry
         procedure :: AddToEntry       => PETScMat_AddToEntry
         procedure :: SetColumn
         procedure :: AddToColumn
         procedure :: Shift
         procedure :: ReShift
         procedure :: PreAssembly
         procedure :: Assembly
         procedure :: GetCSRMatrix
         procedure :: SpecifyBlockInfo => PETSCMat_SpecifyBlockInfo
         procedure :: AddToBlockEntry  => PETScMat_AddToBlockEntry
         procedure :: SetBlockEntry    => PETScMat_SetBlockEntry
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
   subroutine construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,withMPI)
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t)  :: this
      integer, optional, intent(in) :: num_of_Cols
      integer, optional, intent(in) :: num_of_Blocks
#ifdef HAS_PETSC
      PetscInt, optional, intent(in) :: num_of_Rows
      PetscBool, optional, intent(in) :: withMPI
      !---------------------------------------------
      
      if ( .not. present(num_of_Rows) ) then
         ERROR stop 'PETSCMatrix_t needs num_of_Rows'
      end if
      
      this % num_of_Rows = num_of_Rows
      this % withMPI = withMPI
      
      !     PETSc matrix A 
      CALL MatCreate(PETSC_COMM_WORLD,this%A,ierr)                           ; CALL CheckPetscErr(ierr,'error creating A matrix')
      
      IF (withMPI) THEN
         CALL MatSetSizes(this%A,PETSC_DECIDE,PETSC_DECIDE,num_of_Rows,num_of_Rows,ierr)
         CALL CheckPetscErr(ierr,'error setting mat size')
         CALL MatSetType(this%A,MATMPIAIJ, ierr)                              
         CALL CheckPetscErr(ierr,'error in MatSetType')
         CALL MatSetFromOptions(this%A,ierr)                                  
         CALL CheckPetscErr(ierr,'error in MatSetFromOptions')
      ELSE
         CALL MatSetSizes(this%A,PETSC_DECIDE,PETSC_DECIDE,num_of_Rows,num_of_Rows,ierr)
         CALL CheckPetscErr(ierr,'error setting mat size')
         CALL CheckPetscErr(ierr,'error in MatSetType')
         CALL MatSetFromOptions(this%A,ierr)                                  
         CALL CheckPetscErr(ierr,'error in MatSetFromOptions')
!         CALL MatSetOption(this%A,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)  
         CALL MatSetOption(this%A,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE,ierr)                 
         CALL CheckPetscErr(ierr,'error in MatSetOption')
       ENDIF
#else
      integer, optional, intent(in) :: num_of_Rows
      logical, optional, intent(in)    :: WithMPI
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs)
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t) , intent(inout) :: this
      integer, optional    , intent(in)    :: nnzs(:)
#ifdef HAS_PETSC
      PetscInt, optional, intent(in)            :: nnz
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
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this, ForceDiagonal)
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t), intent(inout) :: this
      logical, optional   , intent(in)    :: ForceDiagonal
#ifdef HAS_PETSC
      !---------------------------------------------
      integer :: i
      logical :: mustForceDiagonal
      !---------------------------------------------
      
      if ( present(ForceDiagonal) ) then
         mustForceDiagonal = ForceDiagonal
      else
         mustForceDiagonal = .FALSE.
      end if
      
      call MatZeroEntries(this%A, ierr)
      call CheckPetscErr(ierr,'error in MatZeroEntries')
      
      if (mustForceDiagonal) then
!$omp critical
         do i = 1, this % num_of_Rows
            CALL MatSetValues(this%A, 1 ,i-1,1,i-1,0._RP,ADD_VALUES,ierr)
         end do
!$omp end critical
      end if
      
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
      
      if (abs(value) < JACEPS) return
      
!$omp critical
      CALL MatSetValues(this%A, 1 ,row-1,1,col-1,value,INSERT_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
!$omp end critical
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
   subroutine PETScMat_AddToEntry(this, row, col, value )
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t), intent(inout) :: this
#ifdef HAS_PETSC
      PetscInt            , intent(in)    :: row
      PetscInt            , intent(in)    :: col
      PetscScalar         , intent(in)    :: value
      !---------------------------------------------
      
      if (abs(value) < JACEPS) return
      
!$omp critical
      CALL MatSetValues(this%A, 1 ,row-1,1,col-1,value,ADD_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
!$omp end critical
#else
      INTEGER        , intent(in) :: row
      INTEGER        , intent(in) :: col
      real(kind=RP)  , intent(in) :: value
      STOP ':: PETSc is not linked correctly'
#endif
   end subroutine PETScMat_AddToEntry
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
!$omp critical
      CALL MatSetValues(this%A,nvalues,irow-1,1,icol-1,values,INSERT_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
!$omp end critical
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
      
!$omp critical
      CALL MatSetValues(this%A,nvalues,irow-1,1,icol-1,values,ADD_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
!$omp end critical
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
      INTEGER                                  :: nnz_row(this % num_of_Rows)
      INTEGER      , ALLOCATABLE, DIMENSION(:) :: ACols 
      REAL(KIND=RP), ALLOCATABLE, DIMENSION(:) :: AVals
      PetscErrorCode                           :: ierr
      !---------------------------------------------------------------------------------
      
      !We first get the number of nonzero entries in each row
      DO i = 1, this % num_of_Rows
         CALL MatGetRow(this % A,i-1,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
         CALL CheckPetscErr(ierr,'error in Petsc MatGetRow')
         
         nnz_row(i) = ncols
         
         CALL MatRestoreRow(this % A,i-1,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
         CALL CheckPetscErr(ierr,'error in Petsc MatRestoreRow')
      END DO
      
      CALL Acsr % construct(num_of_Rows = this % num_of_Rows)
      call Acsr % Preallocate(nnzs= nnz_row)
      call Acsr % Reset
      
      DO i = 1, this % num_of_Rows
         
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
   subroutine PETScMat_SpecifyBlockInfo(this,BlockIdx,BlockSize)
      implicit none
      !-arguments-----------------------------------
      CLASS(PETSCMatrix_t), intent(inout) :: this        !<    PETSc matrix
      integer             , intent(in)    :: BlockIdx(:)
      integer             , intent(in)    :: BlockSize(:)
      !---------------------------------------------
      
      safedeallocate(this % BlockIdx)  ; allocate (this % BlockIdx (size(BlockIdx )) )
      safedeallocate(this % BlockSize) ; allocate (this % BlockSize(size(BlockSize)) )
      this % BlockIdx  = BlockIdx
      this % BlockSize = BlockSize
      
   end subroutine PETScMat_SpecifyBlockInfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  --------------------------------------------------------------
   subroutine PETScMat_SetBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      CLASS(PETSCMatrix_t), intent(inout) :: this        !<    PETSc matrix
      integer             , intent(in)    :: iBlock, jBlock
      integer             , intent(in)    :: i, j
      real(kind=RP)       , intent(in)    :: value
      !-local-variables-----------------------------
      integer :: row, col
      !---------------------------------------------
      
      if (.not. allocated(this % BlockIdx)) then
         write(STD_OUT,*) 'PETSCMatrix :: Error '
         write(STD_OUT,*) '            :: PETScMat_SetBlockEntry only available after CSR_SpecifyBlockInfo has been called'
         stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % SetEntry(row, col, value)
      
   end subroutine PETScMat_SetBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to add a value to the entries of a block with relative index
!  -----------------------------------------------------------------------
   subroutine PETScMat_AddToBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      CLASS(PETSCMatrix_t), intent(inout) :: this        !<    PETSc matrix
      integer             , intent(in)    :: iBlock, jBlock
      integer             , intent(in)    :: i, j
      real(kind=RP)       , intent(in)    :: value
      !-local-variables-----------------------------
      integer :: row, col
      !---------------------------------------------
      
      if (.not. allocated(this % BlockIdx)) then
         write(STD_OUT,*) 'PETSCMatrix :: Error '
         write(STD_OUT,*) '            :: PETScMat_AddToBlockEntry only available after CSR_SpecifyBlockInfo has been called'
         stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % AddToEntry(row, col, value)
      
   end subroutine PETScMat_AddToBlockEntry
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
