!
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
      INTEGER,        POINTER, CONTIGUOUS :: Cols(:)=>NULL()    ! Column indices that correspond to each value
      INTEGER,        POINTER, CONTIGUOUS :: Rows(:)=>NULL()    ! Row indices (index of first value of each row)
      INTEGER,        POINTER, CONTIGUOUS :: Diag(:)=>NULL()    ! Array containing position of the diagonal entry (handy for some calculations)
      
      INTEGER                             :: NumRows            ! Number of rows of matrix
      INTEGER                             :: NumCols            ! Number of colunms of matrix
      
   CONTAINS
   
      GENERIC, PUBLIC                     :: construct  => CSR_CreateMatnnz, CSR_CreateMatnnzs
      PROCEDURE                           :: assigndiag => CSR_AssignDiag
      PROCEDURE                           :: Visualize  => CSR2Visualize
      
      PROCEDURE, PRIVATE                  :: CSR_CreateMatnnz, CSR_CreateMatnnzs
   END TYPE
   !-----------------------------------------------------------------------------   

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
   
   !----------------------------------------------------------------------------------
   SUBROUTINE CSR2Visualize(this,FileName)
   !  Creates a file <<FileName>> with the entries of the csr matrix
   !  (foor visualizing using python spy)
   !----------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      CLASS(csrMat)               :: this
      CHARACTER(len=*)            :: FileName
      !------------------------------------------
      INTEGER                     :: n, nnz, i
      !------------------------------------------
      
      n = this % NumRows
      nnz = SIZE(this % Values)
      
      OPEN(unit=30,file=FileName)
         WRITE(30,'(I0,X,I0)') n, nnz
         
         DO i=1, n + 1
            WRITE(30,'(I0)') this % Rows(i)
         END DO
         
         DO i=1, nnz
            WRITE(30,'(I0)') this % Cols(i)
         END DO
         
         DO i=1, nnz
            WRITE(30,'(F20.5)') this % Values(i)
         END DO
         
      ClOSE(30)
      
  !----------------------------------------------------------------------------------
  END SUBROUTINE CSR2Visualize
  !----------------------------------------------------------------------------------
  
END MODULE
