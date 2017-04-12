!////////////////////////////////////////////////////////////////////////
!
!      CSR_Matrices.f90
!      Created: 2017-03-21 17:07:00 +0100 
!      By: Andrés Rueda
!
!      Class for sparse csr matrices
!
!////////////////////////////////////////////////////////////////////////
MODULE CSR_Matrices
   USE SMConstants,                 ONLY: RP                  
   IMPLICIT NONE
   
   !-----------------------------------------------------------------------------
   TYPE csrMat
      REAL(KIND=RP),  POINTER, CONTIGUOUS :: Values(:)=>NULL()  ! Values of nonzero entries of matrix
      INTEGER,        POINTER, CONTIGUOUS :: Cols(:)  =>NULL()  ! Column indices that correspond to each value
      INTEGER,        POINTER, CONTIGUOUS :: Rows(:)  =>NULL()  ! Row indices (index of first value of each row)
      INTEGER,        POINTER, CONTIGUOUS :: Diag(:)  =>NULL()  ! Array containing position of the diagonal entry (handy for some calculations)
      
      INTEGER                             :: NumRows            ! Number of rows of matrix
      INTEGER                             :: NumCols            ! Number of colunms of matrix
      
   CONTAINS
   
      GENERIC, PUBLIC                     :: construct  => CSR_CreateMatnnz, CSR_CreateMatnnzs
      PROCEDURE                           :: assigndiag => CSR_AssignDiag
      PROCEDURE                           :: Visualize  => CSR2Visualize
      PROCEDURE                           :: destroy
      PROCEDURE                           :: SetMatShift
      PROCEDURE                           :: SetColumn
      PROCEDURE                           :: SetElem
      
      PROCEDURE, PRIVATE                  :: CSR_CreateMatnnz, CSR_CreateMatnnzs
   END TYPE
   !-----------------------------------------------------------------------------   
   
   PRIVATE
   PUBLIC                                 :: csrMat, CSR_MatVecMul
!
!========
 CONTAINS
!========
!
   !------------------------------------------------------------------------------
   SUBROUTINE CSR_CreateMatnnz(this,NumRows,NumCols,nnz)
   !  Creates a matrix in the CSR format using the same number of nonzero entries
   !  for all rows (nnz)
   !------------------------------------------------------------------------------
      CLASS(csrMat)       :: this             !<> Matrix to be Created
      INTEGER, INTENT(IN) :: NumRows          !<  Num of rows in the matrix
      INTEGER, INTENT(IN) :: NumCols          !<  Num of cols in the matrix
      INTEGER, INTENT(IN) :: nnz              !<  Num of nonzero entries in every row (extend to non-constant value with an interface)
   
    
      !------------------------------------------------------------------------------
      INTEGER             :: i,k,istat
      !------------------------------------------------------------------------------
      
      IF(nnz < 1) STOP 'CSR_Create: Invalid nnz'
      
      ALLOCATE( this % Rows(NumRows+1),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
      
      ALLOCATE( this % Diag(NumRows),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      k = NumRows*nnz
       
      ALLOCATE( this % Cols(k),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      ALLOCATE( this % Values(k), STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      this % NumRows = NumRows
      this % NumCols = NumCols !The matrix can be a rectangular matrix
      this % Rows(1) = 1
       
      DO i=2, NumRows + 1
         this % Rows(i) = this % Rows(i-1) + nnz
      END DO
      
      this % Cols = 0
      this % Diag = 0
       
      this % Values = 0._RP
   !------------------------------------------------------------------------------
   END SUBROUTINE CSR_CreateMatnnz
   !------------------------------------------------------------------------------
   
   !------------------------------------------------------------------------------
   SUBROUTINE CSR_CreateMatnnzs(this,NumRows,NumCols,nnz_row)
   !  Creates a matrix in the CSR format using a different number of nonzero entries
   !  for each row (nnz_row)
   !------------------------------------------------------------------------------
      CLASS(csrMat)       :: this             !<> Matrix to be Created
      INTEGER, INTENT(IN) :: NumRows          !<  Num of rows in the matrix
      INTEGER, INTENT(IN) :: NumCols          !<  Num of cols in the matrix
      INTEGER, INTENT(IN) :: nnz_row(:)       !<  Num of nonzero entries in every row (extend to non-constant value with an interface)
   
    
      !------------------------------------------------------------------------------
      INTEGER             :: i,k,istat
      INTEGER             :: nnz
      !------------------------------------------------------------------------------
      
      nnz = SUM(nnz_row)
      IF(nnz < 1) STOP 'CSR_Create: Invalid nnz'
      
      ALLOCATE( this % Rows(NumRows+1),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
       
      ALLOCATE( this % Diag(NumRows),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      ALLOCATE( this % Cols(nnz),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      ALLOCATE( this % Values(nnz), STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      this % NumRows = NumRows
      this % NumCols = NumCols !The matrix can be a rectangular matrix
      this % Rows(1) = 1
       
      DO i=2, NumRows + 1
         this % Rows(i) = this % Rows(i-1) + nnz_row(i-1)
      END DO
      
      this % Cols = 0
      this % Diag = 0
       
      this % Values = 0._RP
   !------------------------------------------------------------------------------
   END SUBROUTINE CSR_CreateMatnnzs
   !------------------------------------------------------------------------------
  
   !----------------------------------------------------------------------------------
   SUBROUTINE CSR_AssignDiag(A)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !-----------------------------------
      CLASS(csrMat)        :: A           !<> Matrix
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
   SUBROUTINE SetColumn(A, RowIndexes, ColNum, Values )
   !    Adds values to (part of) a column of a csr matrix
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      CLASS(csrMat)                 :: A                 !<> Global matrix
      INTEGER       , INTENT(IN)    :: RowIndexes(1:)    !< Different positions of Column
      INTEGER       , INTENT(IN)    :: ColNum            !< Number of Row/Column
      REAL(KIND=RP) , INTENT(IN)    :: Values(:)         !< Values to be added to global matrivx¿x
      !------------------------------------------------------------------------------ 
      INTEGER :: NumRows, i,Row !,k,l,c,
      !------------------------------------------------------------------------------
      
      NumRows = SIZE(RowIndexes)
      
      IF (NumRows .NE. SIZE(Values)) THEN
         WRITE (*,*) 'CSR_AddToCol: Dimension error (Values-RowIndexes)'
         STOP
      END IF
      
      IF ( ColNum <= 0 ) THEN
         WRITE (*,*) 'CSR_AddToCol: ColNum error'
         STOP
      END IF
      
      DO i=1,NumRows
         Row = RowIndexes(i)
         IF ( Row <=0 ) CYCLE
         
         CALL A % SetElem(Row,ColNum,Values(i))
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
      CLASS(csrMat), INTENT(INOUT)    :: A     !<>Matrix to be changed
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
   !       (supposing that A is ordered)
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      TYPE(csrMat)  , INTENT(IN) :: A         !< Matrix to be searched
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
         WRITE(*,*) 'CSR_Search: Error. Searched position not defined in matrix.'
         STOP
      END IF
   !------------------------------------------------------------------------------
   END FUNCTION CSR_Search
   !------------------------------------------------------------------------------
   
   !----------------------------------------------------------------------------------
   SUBROUTINE CSR2Visualize(this,FileName,FirstRow)
   !  Creates a file <<FileName>> with the entries of the csr matrix
   !  (foor visualizing using python spy)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat)               :: this
      CHARACTER(len=*)            :: FileName
      LOGICAL, OPTIONAL           :: FirstRow   !< Write First row?
      !------------------------------------------
      INTEGER                     :: n, nnz, i, fd
      !------------------------------------------
      
      IF (.NOT. PRESENT(FirstRow)) FirstRow = .TRUE.
      
      n = this % NumRows
      nnz = SIZE(this % Values)
      
      OPEN(newunit=fd,file=FileName)
         IF (FirstRow) WRITE(fd,'(I0,X,I0)') n, nnz
         
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
   SUBROUTINE destroy(this)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat), INTENT(INOUT) :: this
      !------------------------------------------
      
      NULLIFY(this % Rows)
      NULLIFY(this % Cols)
      NULLIFY(this % Values)
      NULLIFY(this % Diag)
   !----------------------------------------------------------------------------------
   END SUBROUTINE destroy
   !----------------------------------------------------------------------------------
   
   !----------------------------------------------------------------------------------
   SUBROUTINE SetMatShift(this,shift)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat), INTENT(INOUT) :: this
      REAL(KIND=RP), INTENT(IN)    :: shift
      !------------------------------------------
      INTEGER                      :: i
      !------------------------------------------
      
      DO i=1, this % NumRows
         this % Values(this % Diag(i)) = this % Values(this % Diag(i)) + shift
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
      TYPE(csrMat)                 , INTENT(IN)  :: A  !< Structure holding matrix
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
END MODULE
