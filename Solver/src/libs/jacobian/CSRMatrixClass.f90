!////////////////////////////////////////////////////////////////////////
!
!      CSRMatrixClass.f90
!      Created: 2017-03-21 17:07:00 +0100 
!      By: Andrés Rueda
!
!      Class for sparse Compressed Sparse Row (CSR) matrices
!////////////////////////////////////////////////////////////////////////
MODULE CSRMatrixClass
   USE SMConstants          , only: RP    
   use GenericMatrixClass              
   use LinkedListMatrixClass, only: LinkedListMatrix_t
#include "Includes.h"
   IMPLICIT NONE
   
   !-----------------------------------------------------------------------------
   TYPE, extends(Matrix_t) :: csrMat_t
      REAL(KIND=RP),  allocatable :: Values(:)  ! Values of nonzero entries of matrix
      INTEGER      ,  allocatable :: Cols(:)    ! Column indices that correspond to each value
      INTEGER      ,  allocatable :: Rows(:)    ! Row indices (index of first value of each row)
      INTEGER      ,  allocatable :: Diag(:)    ! Array containing position of the diagonal entry (handy for some calculations)
      
      INTEGER                     :: NumCols               ! Number of colunms of matrix
      
      ! Variables for matrices with blocks
      INTEGER      ,  allocatable :: BlockIdx(:)  ! Index of first element of block (this is used by the routine CSR_GetBlock).. Note that in the DGSEM, the Jacobian matrices have a block diagonal with the Jacobian information of each element      
      INTEGER      ,  allocatable :: BlockSize(:) ! Size of each block
      integer                     :: n_max_elements
      integer,        allocatable :: firstIdx(:,:)         ! For each row, specifies the position of the beginning of each element column
      type(LinkedListMatrix_t)    :: ListMatrix
      logical                     :: usingListMat
   CONTAINS
   
      procedure :: construct   => CSR_CreateMat
      procedure :: PreAllocate => CSR_PreAllocate
      procedure :: Reset       => CSR_Reset
      procedure :: assigndiag  => CSR_AssignDiag
      procedure :: Visualize   => CSR2Visualize
      procedure :: destruct
      procedure :: Shift       => SetMatShift
      procedure :: SetColumn
      procedure :: SetEntry
      procedure :: GetDense => CSR2Dense
      procedure :: GetBlock => CSR_GetBlock
      procedure :: Assembly
!      procedure                           :: SetFirstIdx => CSR_SetFirstIdx
      procedure :: PreAllocateWithStructure => CSR_PreAllocateWithStructure
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
!
!        Constant number of non-zeros per row
!        ------------------------------------
         
         total_nnz = this % NumRows * nnz
      elseif (present(nnzs)) then
!
!        Number of non-zeros different in every row
!        ------------------------------------------
         
         if (size(nnzs) /= this % NumRows) ERROR stop ':: CSRMatrix: Not consistent nnzs'
         total_nnz = sum(nnzs)
      else
!
!        Unknown number of nonzeros per row
!        -> Preallocating with LinkedListMatrix
!        --------------------------------------
         
         this % usingListMat = .TRUE.
         call this % ListMatrix % construct(this % NumRows)
         return
      end if
      
      IF(total_nnz < 1) STOP ':: Invalid nnz' 
      
      k = this % NumRows * total_nnz
       
      safedeallocate(this % Cols)   ; ALLOCATE( this % Cols(total_nnz),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      safedeallocate(this % Values) ; ALLOCATE( this % Values(total_nnz), STAT=istat )
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
      safedeallocate(self % Cols)   ; allocate(self % Cols(nnz)) 
      safedeallocate(self % Values) ; allocate(self % Values(nnz))

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
!  ----------------------------------------------------------------------------------
   SUBROUTINE CSR_AssignDiag(A)
      IMPLICIT NONE
      !-----------------------------------
      CLASS(csrMat_t) :: A           !<> Matrix
      !-----------------------------------
      INTEGER         :: i, j       !   Counters
      !-----------------------------------
      do i=1, A % NumRows
         
         do j=A % Rows (i), A % Rows (i+1) -1
            if (A % Cols(j) == i) then
               A % Diag(i) = j
               exit
            end if
            if (j == A % Rows (i+1)) then
               write(*,*) 'CSR_AssignDiag: ERROR? - No diagonal entry found in matrix'
            end if
         end do
      end do
   
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
      
      if (this % usingListMat) then
         call this % ListMatrix % SetColumn ( nvalues, irow, icol, values )
      else
         DO i=1,nvalues
            IF ( irow(i) <= 0 ) CYCLE
            
            CALL this % SetEntry(irow(i),icol,values(i))
         END DO
      end if
   !------------------------------------------------------------------------------
   END SUBROUTINE SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Set given value to an element of a CSR matrix
!  ---------------------------------------------
   subroutine SetEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !-local-variables-----------------------------
      integer                         :: k
      !------------------------------------------------------------------------------
      
      IF(row .GT. this % NumRows .OR. col .GT. this % NumCols ) THEN
         WRITE (*,*) 'CSR_SetEntry: Dimension error'
         STOP
      END IF
      
      if (this % usingListMat) then
         call this % ListMatrix % SetEntry(row,col,value)
      else
         k = CSR_Search(this,row,col)
         this % Values(k) = value
      end if
   !------------------------------------------------------------------------------
   end subroutine SetEntry
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
      
      deallocate(this % Rows)
      deallocate(this % Cols)
      deallocate(this % Values)
      deallocate(this % Diag)
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
      
      if ( this % usingListMat ) then
         call this % ListMatrix % shift(shiftval)
      else
         DO i=1, this % NumRows
            this % Values(this % Diag(i)) = this % Values(this % Diag(i)) + shiftval
         END DO 
      end if
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
      
      safedeallocate(this % BlockIdx)  ; allocate (this % BlockIdx (size(BlockIdx)) )
      safedeallocate(this % BlockSize) ; allocate (this % BlockSize(size(BlockSize)) )
      this % BlockIdx = BlockIdx
      this % BlockSize = BlockSize
      
      if ( this % usingListMat ) then
         call this % ListMatrix % getCSRarrays(this % Values, this % Cols, this % Rows)
         call this % AssignDiag()
         call this % ListMatrix % destruct()
         this % usingListMat = .FALSE.
      end if
   end subroutine Assembly

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   
END MODULE CSRMatrixClass
