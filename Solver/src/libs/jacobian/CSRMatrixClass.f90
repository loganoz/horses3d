!////////////////////////////////////////////////////////////////////////
!
!      CSRMatrixClass.f90
!      Created: 2017-03-21 17:07:00 +0100 
!      By: Andrés Rueda
!
!      Class for sparse Compressed Sparse Row (CSR) matrices
! TODO: change contiguous pointers by allocatables... Memory leaking here!!!!
!////////////////////////////////////////////////////////////////////////
MODULE CSRMatrixClass
   USE SMConstants,                 ONLY: RP    
   use GenericMatrixClass              
   IMPLICIT NONE
   
   !-----------------------------------------------------------------------------
   TYPE, extends(Matrix_t) :: csrMat_t
      REAL(KIND=RP),  POINTER, CONTIGUOUS :: Values(:)   =>NULL()  ! Values of nonzero entries of matrix
      INTEGER,        POINTER, CONTIGUOUS :: Cols(:)     =>NULL()  ! Column indices that correspond to each value
      INTEGER,        POINTER, CONTIGUOUS :: Rows(:)     =>NULL()  ! Row indices (index of first value of each row)
      INTEGER,        POINTER, CONTIGUOUS :: Diag(:)     =>NULL()  ! Array containing position of the diagonal entry (handy for some calculations)
      
      INTEGER                             :: NumCols               ! Number of colunms of matrix
      
      ! Variables for matrices with blocks
      INTEGER,        POINTER, CONTIGUOUS :: BlockIdx(:) =>NULL()  ! Index of first element of block (this is used by the routine CSR_GetBlock).. Note that in the DGSEM, the Jacobian matrices have a block diagonal with the Jacobian information of each element      
      INTEGER,        POINTER, CONTIGUOUS :: BlockSize(:)=>NULL()  ! Size of each block
      integer                             :: n_max_elements
      integer,        allocatable         :: firstIdx(:,:)         ! For each row, specifies the position of the beginning of each element column
   CONTAINS
   
      procedure                           :: construct   => CSR_CreateMat
      procedure                           :: PreAllocate => CSR_PreAllocate
      procedure                           :: Reset       => CSR_Reset
      PROCEDURE                           :: assigndiag  => CSR_AssignDiag
      PROCEDURE                           :: Visualize   => CSR2Visualize
      PROCEDURE                           :: destruct
      PROCEDURE                           :: Shift       => SetMatShift
      PROCEDURE                           :: SetColumn
      PROCEDURE                           :: SetElem
      PROCEDURE                           :: GetDense => CSR2Dense
      PROCEDURE                           :: GetBlock => CSR_GetBlock
      procedure                           :: Assembly
!      procedure                           :: SetFirstIdx => CSR_SetFirstIdx
      procedure                           :: PreAllocateWithStructure => CSR_PreAllocateWithStructure
   END TYPE
   !-----------------------------------------------------------------------------   
   
   PRIVATE
   PUBLIC                                 :: csrMat_t, CSR_MatVecMul, Matrix_t
!
!========
 CONTAINS
!========
!
   !------------------------------------------------------------------------------
   SUBROUTINE CSR_CreateMat(this,DimPrb,WithMPI)
   !  Creates a matrix in the CSR format using the same number of nonzero entries
   !  for all rows (nnz)
   !------------------------------------------------------------------------------
      CLASS(csrMat_t)     :: this             !<> Matrix to be Created
      INTEGER, INTENT(IN) :: DimPrb          !<  Num of rows in the matrix
      logical, optional, intent(in) :: WithMPI
      !------------------------------------------------------------------------------
      INTEGER             :: istat
      !------------------------------------------------------------------------------
      
      ALLOCATE( this % Rows(DimPrb+1),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
      
      ALLOCATE( this % Diag(DimPrb),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
      
      this % NumRows = DimPrb
      this % NumCols = DimPrb !The matrix can be a rectangular matrix
      
   !------------------------------------------------------------------------------
   END SUBROUTINE CSR_CreateMat
   !------------------------------------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CSR_PreAllocate(this,nnz,nnzs)
      implicit none
      !-----------------------------------
      CLASS(csrMat_t), INTENT(INOUT):: this             !<> Matrix to be preallocated
      INTEGER, optional, intent(in) :: nnz        !< Num of nonzero entries all rows
      INTEGER, optional, intent(in) :: nnzs(:)     !< Num of nonzero entries in every row 
      !-----------------------------------
      integer :: i,k,istat
      integer :: total_nnz
      !-----------------------------------
      
      if (present(nnz)) then
         total_nnz = this % NumRows * nnz
      elseif (present(nnzs)) then
         if (size(nnzs) /= this % NumRows) ERROR stop ':: CSRMatrix: Not consistent nnzs'
         total_nnz = sum(nnzs)
      else
         ERROR stop ':: CSR matrix needs nnz or nnzs'
      end if
      
      IF(total_nnz < 1) STOP ':: Invalid nnz' 
      
      k = this % NumRows * total_nnz
       
      ALLOCATE( this % Cols(total_nnz),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      ALLOCATE( this % Values(total_nnz), STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
      
      this % Rows(1) = 1
      if (present(nnz)) then
         DO i=2, this % NumRows + 1
            this % Rows(i) = this % Rows(i-1) + nnz
         END DO
      else
         DO i=2, this % NumRows + 1
            this % Rows(i) = this % Rows(i-1) + nnzs(i-1)
         END DO
      end if
      
   end subroutine CSR_PreAllocate

   subroutine CSR_PreAllocateWithStructure(self, nnz, rows, cols, diag)
      class(CSRMat_t),  intent(inout)  :: self
      integer,          intent(in)     :: nnz
      integer,          intent(in)     :: rows(self % NumRows+1)
      integer,          intent(in)     :: cols(nnz)
      integer,          intent(in)     :: diag(self % NumRows)
!
!     ---------------
!     Local variables
!     ---------------
!
      allocate(self % Cols(nnz)) 
      allocate(self % Values(nnz))

      self % Rows = rows
      self % Cols = cols
      self % Values = 0.0_RP
      self % diag = diag

   end subroutine CSR_PreAllocateWithStructure
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CSR_Reset(this)
      implicit none
      !-----------------------------------
      CLASS(csrMat_t), INTENT(INOUT) :: this           !<> Matrix
      !-----------------------------------
      
      this % Cols = 0
      this % Diag = 0
      this % Values = 0._RP
   end subroutine CSR_Reset
      
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !----------------------------------------------------------------------------------
   SUBROUTINE CSR_AssignDiag(A)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !-----------------------------------
      CLASS(csrMat_t)       :: A           !<> Matrix
      !-----------------------------------
      INTEGER               :: i, j       !   Counters
      INTEGER, POINTER      :: row_ptr(:) !   Row pointer
      !-----------------------------------
      
      DO i=1, A % NumRows
         row_ptr => A % Cols (A % Rows (i) : A % Rows (i+1) -1)
         DO j=A % Rows (i), A % Rows (i+1) -1
         IF (A % Cols(j) == i) THEN
            A % Diag(i) = j
            EXIT
         END IF
         IF (j == A % Rows (i+1)) THEN
            WRITE(*,*) 'CSR_AssignDiag: ERROR? - No diagonal entry found in matrix'
         END IF
         END DO
      END DO
   
   !----------------------------------------------------------------------------------
   END SUBROUTINE CSR_AssignDiag
   !----------------------------------------------------------------------------------
      
   !------------------------------------------------------------------------------
   SUBROUTINE SetColumn(this, nvalues, irow, icol, values )
   !    Adds values to (part of) a column of a csr matrix
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      CLASS(csrMat_t), INTENT(INOUT) :: this              !<> Global matrix
      INTEGER        , INTENT(IN)    :: nvalues
      INTEGER        , INTENT(IN)    :: irow(1:)    !< Different positions of Column
      INTEGER        , INTENT(IN)    :: icol            !< Number of Row/Column
      REAL(KIND=RP)  , INTENT(IN)    :: values(:)         !< Values to be added to global matrivx¿x
      !------------------------------------------------------------------------------ 
      INTEGER :: i !,k,l,c,
      !------------------------------------------------------------------------------
      
      IF (nvalues .NE. SIZE(Values)) THEN
         WRITE (*,*) 'CSR_AddToCol: Dimension error (Values-RowIndexes)'
         STOP
      END IF
      
      IF ( icol <= 0 ) THEN
         WRITE (*,*) 'CSR_AddToCol: icol error'
         STOP
      END IF
      
      DO i=1,nvalues
         IF ( irow(i) <= 0 ) CYCLE
         
         CALL this % SetElem(irow(i),icol,values(i))
      END DO
      
   !------------------------------------------------------------------------------
   END SUBROUTINE SetColumn
   !------------------------------------------------------------------------------
  
   !------------------------------------------------------------------------------
   SUBROUTINE SetElem( A,i,j,Aij)
   !   Set given value to an element of a CSR matrix
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      CLASS(csrMat_t), INTENT(INOUT)  :: A     !<>Matrix to be changed
      INTEGER      , INTENT(IN)       :: i     !< row of the matrix element
      INTEGER      , INTENT(IN)       :: j     !< column number of the matrix element
      REAL(KIND=RP), INTENT(IN)       :: Aij   !< new value of the matrix element
      !------------------------------------------------------------------------------ 
      INTEGER                         :: k
      !------------------------------------------------------------------------------
   
      IF(i .GT. A % NumRows .OR. j .GT. A % NumCols ) THEN
         WRITE (*,*) 'CSR_SetElement: Dimension error'
         STOP
      END IF
      
      k = CSR_Search(A,i,j)
      A % Values(k) = Aij
   !------------------------------------------------------------------------------
   END SUBROUTINE SetElem
   !------------------------------------------------------------------------------
   
   !------------------------------------------------------------------------------
   FUNCTION CSR_Search (A,i,j) RESULT(k)
   !    Obtains the position k for the information of A(i,j) --> A % Values (k)
   !    If the position is not contained in the sparse matrix A, k = 0 is returned
   !       (supposing that A is ordered)
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      TYPE(csrMat_t), INTENT(IN) :: A         !< Matrix to be searched
      INTEGER       , INTENT(IN) :: i, j      !< Position to bi searched
      INTEGER                    :: k         !> Position in A%Values of A(i,j)
      !------------------------------------------------------------
    
      IF(i .GT. A % NumRows .OR. j .GT. A % NumCols ) THEN
         WRITE (*,*) 'CSR_Search: Dimension error'
         STOP
      END IF

      DO k = A % Rows(i), A % Rows(i+1)-1
         IF(A % Cols(k) == j) EXIT
      END DO

      IF (k == A % Rows(i+1)) THEN
         k = 0
      END IF
   !------------------------------------------------------------------------------
   END FUNCTION CSR_Search
   !------------------------------------------------------------------------------
   
   !----------------------------------------------------------------------------------
   SUBROUTINE CSR2Visualize(this,FileName,FirstRow)
   !  Creates a file <<FileName>> with the entries of the csr matrix
   !  (for visualizing using python spy)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat_t)             :: this
      CHARACTER(len=*)            :: FileName
      LOGICAL, OPTIONAL           :: FirstRow   !< Write First row?
      !------------------------------------------
      INTEGER                     :: n, nnz, i, fd
      LOGICAL                     :: First
      !------------------------------------------
      
      IF (.NOT. PRESENT(FirstRow)) First = .TRUE.
      
      n = this % NumRows
      nnz = SIZE(this % Values)
      
      OPEN(newunit=fd,file=FileName)
         IF (First) WRITE(fd,'(I0,X,I0)') n, nnz
         
         DO i=1, n + 1
            WRITE(fd,'(I0)') this % Rows(i)
         END DO
         
         DO i=1, nnz
            WRITE(fd,'(I0)') this % Cols(i)
         END DO
         
         DO i=1, nnz
            WRITE(fd,*) this % Values(i)
         END DO
         
      ClOSE(fd)
      
   !----------------------------------------------------------------------------------
   END SUBROUTINE CSR2Visualize
   !----------------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------------
   SUBROUTINE destruct(this)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat_t), INTENT(INOUT) :: this
      !------------------------------------------
      
      NULLIFY(this % Rows)
      NULLIFY(this % Cols)
      NULLIFY(this % Values)
      NULLIFY(this % Diag)
   !----------------------------------------------------------------------------------
   END SUBROUTINE destruct
   !----------------------------------------------------------------------------------
   
   !----------------------------------------------------------------------------------
   SUBROUTINE SetMatShift(this,shiftval)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat_t), INTENT(INOUT)     :: this
      real(kind=RP)  , INTENT(IN)        :: shiftval
      !------------------------------------------
      INTEGER                        :: i
      !------------------------------------------
      
      DO i=1, this % NumRows
         this % Values(this % Diag(i)) = this % Values(this % Diag(i)) + shiftval
      END DO 
      
   !----------------------------------------------------------------------------------
   END SUBROUTINE SetMatShift
   !----------------------------------------------------------------------------------
   
   !------------------------------------------------------------------------------
   FUNCTION CSR_MatVecMul( A,u) RESULT(v)
   ! Matrix vector product (v = Au) for a CSR matrix .... v needs to be allocated beforehand
   !   Assuming there's MKL
   !------------------------------------------------------------------------------
      REAL(KIND=RP), DIMENSION(:)  , INTENT(IN)  :: u  !< Vector to be multiplied
      TYPE(csrMat_t)               , INTENT(IN)  :: A  !< Structure holding matrix
      REAL(KIND=RP), DIMENSION(A % NumRows)      :: v  !> Result vector 
      !------------------------------------------------------------------------------
      INTEGER       :: i,j
      REAL(KIND=RP) :: rsum
      !------------------------------------------------------------------------------
#ifdef HAS_MKL
      !The interface is not really necessary, but it's better to have it
      INTERFACE
         SUBROUTINE mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
            USE SMConstants
            
            CHARACTER      :: transa
            INTEGER        :: m
            REAL(KIND=RP)  :: a(*)
            INTEGER        :: ia(*), ja(*)
            REAL(KIND=RP)  :: x(*), y(*)
         END SUBROUTINE mkl_dcsrgemv
      END INTERFACE
#endif
      !------------------------------------------------------------------------------
    
      IF (A % NumCols .NE. SIZE(U) .OR. A % NumRows .NE. SIZE(v)) THEN
         STOP 'CSR_MatVecMul: Error - u dimensions mismatch'
      END IF
    
#ifdef HAS_MKL
      CALL mkl_dcsrgemv('N', A % NumRows, A % Values, A % Rows, A % Cols, u, v)
#else
      !$omp parallel do private(j,rsum)
      DO i=1,A % NumRows
         rsum = 0.0d0
         DO j=A % Rows(i),A % Rows(i+1)-1
            rsum = rsum + u(A % Cols(j)) * A % Values(j)
         END DO
         v(i) = rsum
      END DO
      !$omp end parallel do
#endif
   !------------------------------------------------------------------------------
   END FUNCTION
   !------------------------------------------------------------------------------
   
   !----------------------------------------------------------------------------------
   SUBROUTINE CSR2Dense(A,Mat) 
   !     Transforms a matrix of type Matrix_t into a common real(dp) matrix
   !----------------------------------------------------------------------------------
      CLASS(csrMat_t)           , INTENT(IN)  :: A          !< CSR matric
      REAL(KIND=RP), ALLOCATABLE, INTENT(OUT) :: Mat(:,:)   !> Dense matrix
      !------------------------------------------------
      INTEGER                                 :: i, k       !  Counters
      !------------------------------------------------
      
      IF (ALLOCATED(Mat)) DEALLOCATE(Mat)
      ALLOCATE(Mat(A%NumRows,A%NumCols))
      
      Mat = 0.d0
      DO i=1, A % NumRows
         DO k=A % Rows(i), A % Rows(i+1)-1
            Mat(i,A%Cols(k)) = A % Values (k)
         END DO
      END DO
   !----------------------------------------------------------------------------------
   END SUBROUTINE CSR2Dense
   !----------------------------------------------------------------------------------
   
   !------------------------------------------------------------------------------
   FUNCTION CSR_GetBlock( A,Num,N) RESULT(Mat)
   !   Get a block of the diagonal of matrix A and store it as a dense matrix
   !  IMPORTANT: The CSR
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      CLASS(csrMat_t)           :: A          !< Matrix to be read
      INTEGER, INTENT(IN)       :: Num        !< Number of diagonal element
      INTEGER, INTENT(IN)       :: N          !< Size of block
      REAL(KIND=RP)             :: Mat(N,N)   !> Value of the matrix element
      !------------------------------------------------------------------------------ 
      INTEGER                   :: RC0      ! Index of first row/column of block in global sparse matrix
      INTEGER                   :: RCf      ! Index of last row/column of block in global sparse matrix
      INTEGER                   :: i, j     ! Row/Column index of output matrix
      INTEGER                   :: ii, jj   ! Row/Column index in global sparse matrix
      INTEGER                   :: k, l     ! Variables containing positions in sparse matrix arrays
      !------------------------------------------------------------------------------
      
      RC0 = A % BlockIdx(Num) 
      RCf  = RC0 + N - 1
      
      i=0
      
      Mat = 0._RP
      
      ! Loop over the rows where the desired block is contained
      DO ii = RC0, RCf
         i = i+1
         
         ! Search first nonzero of this row in block
         k = 0
         DO jj = RC0, RCf
            k = CSR_Search(A,ii,jj)
            IF (k /= 0) EXIT
         END DO
         IF (k == 0) CYCLE
         
         ! Loop over the nonzeros following the first of the row to get all the columns
         DO l = k, A % Rows(ii+1) - 1
            IF (A % Cols(l) > RCf) EXIT
            
            j = A % Cols(l) - RC0 + 1
            Mat(i,j) = A % Values(l)
         END DO
      END DO
       
   !------------------------------------------------------------------------------
   END FUNCTION CSR_GetBlock
   !------------------------------------------------------------------------------
   
   subroutine Assembly(this,BlockIdx,BlockSize)
      implicit none
      !---------------------------------------------
      class(csrMat_t),     intent(inout)   :: this
      integer, target, optional    ,     intent(in)      :: BlockIdx(:)
      integer, target, optional, intent(in)    :: BlockSize(:)
      !---------------------------------------------
      
      allocate (this % BlockIdx (size(BlockIdx)) )
      allocate (this % BlockSize(size(BlockSize)) )
      this % BlockIdx = BlockIdx
      this % BlockSize = BlockSize
   end subroutine Assembly

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   
END MODULE CSRMatrixClass
