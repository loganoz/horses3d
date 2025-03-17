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
   use JacobianDefinitions , only: JACEPS
   use MPI_Process_Info    , only: MPI_Process
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
      Vec         :: rowvec   ! Auxiliary vector with size = num_of_Rows
      Vec         :: colvec   ! Auxiliary vector with size = num_of_Cols
      PetscInt    :: num_of_totalRows
#endif
      
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
         procedure :: SpecifyBlockInfo       => PETSCMat_SpecifyBlockInfo
         procedure :: AddToBlockEntry        => PETScMat_AddToBlockEntry
         procedure :: SetBlockEntry          => PETScMat_SetBlockEntry
         procedure :: ForceAddToEntry        => PETScMat_ForceAddToEntry
         procedure :: ForceAddToBlockEntry   => PETScMat_ForceAddToBlockEntry
         procedure :: MatMatMul        => PETScMat_MatMatMul
         procedure :: MatVecMul        => PETScMat_MatVecMul
         procedure :: MatAdd           => PETScMat_MatAdd
         procedure :: ConstructFromDiagBlocks   => PETScMat_ConstructFromDiagBlocks
         procedure :: constructWithCSRArrays => PETScMat_constructWithCSRArrays
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
   subroutine construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,num_of_TotalRows,withMPI)
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t)  :: this
      integer, optional, intent(in) :: num_of_Blocks
#ifdef HAS_PETSC
      PetscInt, optional, intent(in) :: num_of_Rows
      PetscInt, optional, intent(in) :: num_of_Cols
      PetscInt, optional, intent(in) :: num_of_TotalRows
      PetscBool, optional, intent(in) :: withMPI
      !---------------------------------------------
      PetscBool :: hasMPI
      PetscInt  :: num_of_totalCols
      !---------------------------------------------
!
!     Initialize PETSc environment... If it was already done by the solver, it's alright
!     ----------------------------------------------------------------------------------
      call PetscInitialize(PETSC_NULL_character,ierr)
      
      if ( .not. present(num_of_Rows) ) then
         error stop 'PETSCMatrix_t needs num_of_Rows'
      end if
      
      if ( present(num_of_Cols) ) then
         this % num_of_Cols = num_of_Cols
      else
         this % num_of_Cols = num_of_Rows
      end if
!~      if ( present(withMPI) ) then
!~         hasMPI = withMPI
!~      else
!~         hasMPI = PETSC_FALSE
!~      end if
      
      this % num_of_Rows = num_of_Rows
      
      if ( MPI_Process % doMPIAction .and. present(num_of_totalRows) ) then
         this % withMPI = PETSC_TRUE
         this % num_of_totalRows = num_of_totalRows
      else
         this % withMPI = PETSC_FALSE
         this % num_of_totalRows = num_of_Rows
      end if
      
      !     PETSc matrix A 
      CALL MatCreate(PETSC_COMM_WORLD,this%A,ierr)                           ; CALL CheckPetscErr(ierr,'error creating A matrix')
      
      if (this % withMPI) then
         num_of_totalCols = num_of_totalRows
         CALL MatSetSizes(this%A,num_of_Rows,PETSC_DECIDE,num_of_totalRows,num_of_totalRows,ierr) ! At the moment only for square matrices...
         CALL CheckPetscErr(ierr,'error setting mat size')
         CALL MatSetType(this%A,MATMPIAIJ, ierr)                       
         CALL CheckPetscErr(ierr,'error in MatSetType')
!~         CALL MatSetFromOptions(this%A,ierr)                                  
!~         CALL CheckPetscErr(ierr,'error in MatSetFromOptions')
      else
         num_of_totalCols = this % num_of_Cols
         CALL MatSetSizes(this%A,num_of_Rows,this % num_of_Cols,num_of_Rows,this % num_of_Cols,ierr)
         CALL CheckPetscErr(ierr,'error setting mat size')
         CALL MatSetType(this%A,MATSEQAIJ, ierr)
         CALL CheckPetscErr(ierr,'error in MatSetType')
!~         CALL MatSetFromOptions(this%A,ierr)                                  
!~         CALL CheckPetscErr(ierr,'error in MatSetFromOptions')
!~         CALL MatSetOption(this%A,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)  
!~         CALL MatSetOption(this%A,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE,ierr)                 
!~         CALL CheckPetscErr(ierr,'error in MatSetOption')
      end if
      
!     Construct auxiliary vectors
!     --------------------------
      call VecCreate(PETSC_COMM_WORLD,this % rowvec,ierr)                     ; call CheckPetscErr(ierr,'error creating Petsc vector')
      call VecSetSizes(this % rowvec,this % num_of_Rows,this % num_of_totalRows,ierr)    ; call CheckPetscErr(ierr,'error setting Petsc vector options')
      call VecSetFromOptions(this % rowvec,ierr)                              ; call CheckPetscErr(ierr,'error setting Petsc vector options')
      
      call VecCreate(PETSC_COMM_WORLD,this % colvec,ierr)                     ; call CheckPetscErr(ierr,'error creating Petsc vector')
      call VecSetSizes(this % colvec,PETSC_DECIDE,num_of_totalCols,ierr)    ; call CheckPetscErr(ierr,'error setting Petsc vector options')
      call VecSetFromOptions(this % colvec,ierr)                              ; call CheckPetscErr(ierr,'error setting Petsc vector options')
      
#else
      integer, optional, intent(in) :: num_of_Rows
      integer, optional, intent(in) :: num_of_Cols
      integer, optional, intent(in) :: num_of_TotalRows
      logical, optional, intent(in)    :: WithMPI
      error stop ':: PETSc is not linked correctly'
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
      
      if (.not. present(nnz)) error stop ':: PETSc matrix needs nnz'
      
      IF(this%withMPI) THEN
         CALL MatMPIAIJSetPreallocation(this%A, nnz/7, PETSC_NULL_INTEGER, nnz*6/7, PETSC_NULL_INTEGER,ierr) ! hard-coded: 6 neighbors... Changhe!
         CALL CheckPetscErr(ierr, 'error in MatMPIAIJSetPreallocation')
      ELSE
         CALL MatSeqAIJSetPreallocation(this%A, nnz, PETSC_NULL_INTEGER,ierr)
         CALL CheckPetscErr(ierr, 'error in MatSeqAIJSetPreallocation')
      ENDIF
      
#else
      INTEGER, optional, intent(in)  :: nnz
      error stop ':: PETSc is not linked correctly'
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
         if (this % withMPI) then
            call VecZeroEntries(this % rowvec, ierr)
            
            call VecAssemblyBegin(this % rowvec, ierr);  call CheckPetscErr(ierr," Assembly B in PETSc Begin")      
            call VecAssemblyEnd  (this % rowvec, ierr)  ;  call CheckPetscErr(ierr," Assembly B in PETSc End")  
            
            call MatDiagonalSet(this % A, this % rowvec, ADD_VALUES, ierr)
            ! This is failing: SOLVE!!
         else
            do i = 1, this % num_of_Rows
               CALL MatSetValues(this%A, 1 ,i-1,1,i-1,0._RP,ADD_VALUES,ierr)
            end do
         end if
!$omp end critical
      end if
      
#else
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETScMat_AddToEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETScMat_ForceAddToEntry(this, row, col, value )
      implicit none
      !---------------------------------------------
      class(PETSCMatrix_t), intent(inout) :: this
#ifdef HAS_PETSC
      PetscInt            , intent(in)    :: row
      PetscInt            , intent(in)    :: col
      PetscScalar         , intent(in)    :: value
      !---------------------------------------------
      
!$omp critical
      CALL MatSetValues(this%A, 1 ,row-1,1,col-1,value,ADD_VALUES,ierr)
      CALL CheckPetscErr(ierr, 'error in MatSetValues')
!$omp end critical
#else
      INTEGER        , intent(in) :: row
      INTEGER        , intent(in) :: col
      real(kind=RP)  , intent(in) :: value
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETScMat_ForceAddToEntry
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      error stop ':: PETSc is not linked correctly'
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
      safedeallocate(this % BlockSizes) ; allocate (this % BlockSizes(size(BlockSize)) )
      this % BlockIdx  = BlockIdx
      this % BlockSizes = BlockSize
      
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
         error stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % SetEntry(row, col, value)
      
   end subroutine PETScMat_SetBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  --------------------------------------------------------------
   subroutine PETScMat_MatMatMul(A,B,Cmat,trans)
      implicit none
      !-arguments--------------------------------------------------------------------
      class(PETSCMatrix_t) , intent(in)      :: A       !< Structure holding matrix
      class(Matrix_t)      , intent(in)      :: B       !< Structure holding matrix
      class(Matrix_t)      , intent(inout)   :: Cmat    !< Structure holding matrix
      logical, optional    , intent(in)      :: trans   !< A matrix is transposed?
      !-local-variables-----------------------------
      integer :: row, col
      !---------------------------------------------
#ifdef HAS_PETSC
!
!     Since the arguments must be class(Matrix_t), an extra check is needed
!     ---------------------------------------------------------------------
      select type(B) ; class is(PETSCMatrix_t) ; select type (Cmat) ; class is (PETSCMatrix_t)
!
!     Perform MatMatMul
!     -----------------
      if ( present(trans) ) then
         if (trans) error stop 'PETScMat_MatMatMul :: ERROR: trans not implemented'
      end if 
      call MatMatMult( A % A, B % A,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL, Cmat % A, ierr)
      call CheckPetscErr(ierr,"PETScMat_MatMatMul: Problem doing MatMatMult")
!
!     Finish extra check
!     ------------------
      class default
         error stop ':: Wrong type of arguments in CSR_MatMatMul'
      end select ; end select
#endif
   end subroutine PETScMat_MatMatMul
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------
!  MatVecMul:
!  Matrix vector product (v = Au) being A a CSR matrix
!  -> v needs to be allocated beforehand
!  ----------------------------------------------------
   function PETScMat_MatVecMul( A,u, trans) result(v)
      !-arguments--------------------------------------------------------------------
      class(PETSCMatrix_t) , intent(inout)  :: A  !< Structure holding matrix
      logical, optional    , intent(in)     :: trans   !< A matrix is transposed?
#ifdef HAS_PETSC
      PetscScalar          , intent(in)     :: u(A % num_of_Cols)  !< Vector to be multiplied
      PetscScalar                           :: v(A % num_of_Rows)  !> Result vector 
      !------------------------------------------------------------------------------
      PetscInt :: i
      !------------------------------------------------------------------------------
      
      call VecSetValues    (A % colvec, A % num_of_Cols, [(i, i=0, A % num_of_Cols-1)] , u, INSERT_VALUES, ierr)
      call CheckPetscErr(ierr,"PETScMat_MatVecMul: VecSetValues A % colvec in PETSc Begin")      
      
      call VecAssemblyBegin(A % colvec, ierr)   ;  call CheckPetscErr(ierr,"PETScMat_MatVecMul: Assembly this % colvec in PETSc Begin")      
      call VecAssemblyEnd  (A % colvec, ierr)   ;  call CheckPetscErr(ierr,"PETScMat_MatVecMul: Assembly this % colvec in PETSc End")  
      
      call MatMult(A % A, A % colvec, A % rowvec, ierr)
      call CheckPetscErr(ierr,"PETScMat_MatVecMul: Doing MatMult")
      
      call VecGetValues(A % rowvec, A % num_of_Rows ,[(i, i=0, A % num_of_Rows-1)],v, ierr)
      call CheckPetscErr(ierr, 'PETScMat_MatVecMul: error in VecGetValue v')
      
#else
      real(kind=RP)    , intent(in)  :: u(A % num_of_Cols)  !< Vector to be multiplied
      real(kind=RP)                  :: v(A % num_of_Rows)  !> Result vector 
#endif
   end function PETScMat_MatVecMul
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------
!  MatAdd:
!  Matrix addition: Cmat = A + Factor*B
!  ----------------------
   subroutine PETScMat_MatAdd(A,B,Cmat,Factor)
      implicit none
      !-arguments--------------------------------------------------------------------
      class(PETSCMatrix_t) , intent(in)      :: A       !< Structure holding matrix
      class(Matrix_t)      , intent(in)      :: B       !< Structure holding matrix
      class(Matrix_t)      , intent(inout)   :: Cmat    !< Structure holding matrix
#ifdef HAS_PETSC
      PetscScalar          , intent(in)      :: Factor  !< Factor for addition
      !------------------------------------------------------------------------------
!
!     Since the arguments must be class(Matrix_t), an extra check is needed
!     ---------------------------------------------------------------------
      select type(B) ; class is(PETSCMatrix_t) ; select type (Cmat) ; class is (PETSCMatrix_t)
!
!     Perform MatAdd
!     --------------
      
      ! Copy matrix A in Cmat
      call MatDuplicate(A % A,MAT_COPY_VALUES,Cmat % A, ierr)
      call CheckPetscErr(ierr,"PETScMat_MatAdd: Problem copying A in Cmat")
      
      ! Perform operation
      call MatAXPY(Cmat % A, Factor, B % A, DIFFERENT_NONZERO_PATTERN, ierr)
      call CheckPetscErr(ierr,"PETScMat_MatAdd: Problem Adding matrices")
!
!     Finish extra check
!     ------------------
      class default
         error stop ':: Wrong type of arguments in CSR_MatMatMul'
      end select ; end select
#else
      real(kind=RP)  , intent(in)  :: Factor  !< Factor for addition
#endif
   end subroutine PETScMat_MatAdd
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
         error stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % AddToEntry(row, col, value)
      
   end subroutine PETScMat_AddToBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to add a value to the entries of a block with relative index (Force)
!  -----------------------------------------------------------------------
   subroutine PETScMat_ForceAddToBlockEntry(this, iBlock, jBlock, i, j, value )
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
         error stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % ForceAddToEntry(row, col, value)
      
   end subroutine PETScMat_ForceAddToBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------
!  Construct a general sparse matrix from diagonal blocks
!  ------------------------------------------------------
   subroutine PETScMat_ConstructFromDiagBlocks(this, num_of_Blocks, Blocks, BlockIdx, BlockSizes)
      implicit none
      !-arguments-----------------------------------
      class(PETSCMatrix_t) , intent(inout) :: this
      integer              , intent(in)    :: num_of_Blocks
      type(DenseBlock_t)   , intent(in)    :: Blocks(num_of_Blocks)
      integer              , intent(in)    :: BlockIdx(num_of_Blocks+1)
      integer              , intent(in)    :: BlockSizes(num_of_Blocks)
      !-local-variables-----------------------------
      integer :: bID, j, i, j_offset
      !---------------------------------------------
      
!     Construct matrix
!     ----------------
      call this % construct (num_of_Rows = sum(BlockSizes))
      call this % PreAllocate( nnz = maxval(BlockSizes) )
      call this % reset
      call this % SpecifyBlockInfo (BlockIdx, BlockSizes)
      
!     Fill it with the right entries
!     ------------------------------
      
      j_offset = 0
      do bID=1, num_of_Blocks
         do j=1, BlockSizes(bID)
            
            do i=1, BlockSizes(bID)
               call this % SetBlockEntry (bID, bID, i, j, Blocks(bID) % Matrix(i,j) )
            end do
            
         end do
      end do
      
      call this % assembly
      
   end subroutine PETScMat_ConstructFromDiagBlocks
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Construct a PETSc matrix from CSR arrays
!  ----------------------------------------
   subroutine PETScMat_constructWithCSRArrays(this,Rows,Cols,Values,num_of_Cols)
      !-arguments-----------------------------------
      class(PETSCMatrix_t)          :: this       !<> Matrix to be Created
      integer, optional, intent(in) :: num_of_Cols
#ifdef HAS_PETSC
      PetscInt         , intent(in) :: Rows(:)    ! Row indices (index of first value of each row)
      PetscInt         , intent(in) :: Cols(:)    ! Column indices that correspond to each value
      PetscScalar      , intent(in) :: Values(:)  ! Values of nonzero entries of matrix
      !---------------------------------------------
      
      if (this % withMPI) then
         call MatMPIAIJSetPreallocationCSR(this % A, Rows-1, Cols-1, Values, ierr)
         call CheckPetscErr(ierr," MatMPIAIJSetPreallocationCSR")
      else
         call MatSeqAIJSetPreallocationCSR(this % A, Rows-1, Cols-1, Values, ierr)
         call CheckPetscErr(ierr," MatMPIAIJSetPreallocationCSR")
      end if
#else
      integer          , intent(in) :: Rows(:)    ! Row indices (index of first value of each row)
      integer          , intent(in) :: Cols(:)    ! Column indices that correspond to each value
      real(kind=RP)    , intent(in) :: Values(:)  ! Values of nonzero entries of matrix
#endif
   end subroutine PETScMat_constructWithCSRArrays
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !---------------------------------------------
      CLASS(PETSCMatrix_t),     intent(inout)     :: this
      !---------------------------------------------
      
      safedeallocate(this % BlockIdx)
      safedeallocate(this % BlockSizes)
#ifdef HAS_PETSC
      CALL MatDestroy(this%A,ierr)
      CALL CheckPetscErr(ierr," A destruction")  
      
      call VecDestroy(this % rowvec,ierr) ; call CheckPetscErr(ierr,"Problem destructung vector")  
      call VecDestroy(this % colvec,ierr) ; call CheckPetscErr(ierr,"Problem destructung vector")  
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
         error stop
      ENDIF
#else
      INTEGER                                      :: ierr
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine CheckPetscErr
   
end module PETScMatrixClass