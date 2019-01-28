!
!//////////////////////////////////////////////////////
!
!   @File:    CSRMatrixClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: 
!   @Last revision date: Mon Jan 28 16:45:16 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: d7287eb57ed10f49c4457713557b614ea83bee6f
!
!//////////////////////////////////////////////////////
!
MODULE CSRMatrixClass
   USE SMConstants          , only: RP, STD_OUT   
   use GenericMatrixClass              
   use LinkedListMatrixClass, only: LinkedListMatrix_t
#include "Includes.h"
   IMPLICIT NONE
   
   !-----------------------------------------------------------------------------
   TYPE, extends(Matrix_t) :: csrMat_t
      real(kind=RP),  allocatable :: Values(:)  ! Values of nonzero entries of matrix
      integer      ,  allocatable :: Cols(:)    ! Column indices that correspond to each value
      integer      ,  allocatable :: Rows(:)    ! Row indices (index of first value of each row)
      integer      ,  allocatable :: Diag(:)    ! Array containing position of the diagonal entry (handy for some calculations)
      
      integer                     :: num_of_Cols               ! Number of colunms of matrix
      
      ! Variables for matrices with blocks
      integer      ,  allocatable :: BlockIdx(:)  ! Index of first element of block (this is used by the routine CSR_GetBlock).. Note that in the DGSEM, the Jacobian matrices have a block diagonal with the Jacobian information of each element      
      integer      ,  allocatable :: BlockSize(:) ! Size of each block
      
      integer,        allocatable :: firstIdx(:,:)         ! For each row, specifies the position of the beginning of each element column
      type(LinkedListMatrix_t)    :: ListMatrix
      logical                     :: usingListMat
   contains
   
      procedure :: construct           => CSR_Construct
      procedure :: constructWithArrays => CSR_constructWithArrays
      procedure :: PreAllocate         => CSR_PreAllocate
      procedure :: Reset               => CSR_Reset
      procedure :: ResetBlock          => CSR_ResetBlock
      procedure :: assigndiag          => CSR_AssignDiag
      procedure :: Visualize           => CSR2Visualize
      procedure :: destruct
      procedure :: Shift               => SetMatShift
      procedure :: SetColumn
      procedure :: SetEntry            => CSR_SetEntry
      procedure :: AddToEntry          => CSR_AddToEntry
      procedure :: GetDense            => CSR2Dense
      procedure :: GetBlock            => CSR_GetBlock
      procedure :: Assembly            => CSR_Assembly
      procedure :: SpecifyBlockInfo    => CSR_SpecifyBlockInfo
      procedure :: AddToBlockEntry     => CSR_AddToBlockEntry
      procedure :: SetBlockEntry       => CSR_SetBlockEntry
!      procedure                           :: SetFirstIdx => CSR_SetFirstIdx
      procedure :: PreAllocateWithStructure => CSR_PreAllocateWithStructure
   END TYPE
   !-----------------------------------------------------------------------------   
   
   PRIVATE
   public :: csrMat_t, CSR_MatVecMul, CSR_MatMatMul, CSR_MatAdd, Matrix_t
   
   interface
      subroutine mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
         use SMConstants
         character      :: transa
         integer        :: m
         real(kind=RP)  :: a(*)
         integer        :: ia(*), ja(*)
         real(kind=RP)  :: x(*), y(*)
      end subroutine mkl_dcsrgemv
   end interface
   
!
!========
 CONTAINS
!========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------
!  Common constructor for CSR matrices
!  -----------------------------------
   subroutine CSR_Construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,num_of_rows_reduced,WithMPI)
      !-arguments-----------------------------------
      class(csrMat_t)               :: this             !<> Matrix to be Created
      integer, optional, intent(in) :: num_of_Rows
      integer, optional, intent(in) :: num_of_Cols
      integer, optional, intent(in) :: num_of_Blocks
      integer, optional, intent(in) :: num_of_rows_reduced
      logical, optional, intent(in) :: WithMPI
      !-local-variables-----------------------------
      integer             :: istat
      integer             :: NumCols
      !---------------------------------------------
      
      if ( .not. present(num_of_Rows) ) then
         ERROR stop 'csrMat_t needs num_of_Rows'
      end if
      if ( present(num_of_Cols) ) then
         NumCols = num_of_Cols
      else
         NumCols = num_of_Rows
      end if
      
      allocate( this % Rows(num_of_Rows+1),stat=istat )
      if ( istat .NE. 0 ) write(*,*) 'CSR_construct: Memory allocation error'
      
      allocate( this % Diag(num_of_Rows), stat=istat )
      if ( istat .NE. 0 ) write(*,*) 'CSR_construct: Memory allocation error'
      
      this % num_of_Rows = num_of_Rows
      this % num_of_Cols = NumCols     !The matrix can be a rectangular matrix
      
   end subroutine CSR_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------
!  Constructor that is fed with the arrays directly
!  ------------------------------------------------
   subroutine CSR_constructWithArrays(this,Rows,Cols,Values,num_of_Cols)
      !-arguments-----------------------------------
      class(csrMat_t)               :: this       !<> Matrix to be Created
      integer          , intent(in) :: Rows(:)    ! Row indices (index of first value of each row)
      integer          , intent(in) :: Cols(:)    ! Column indices that correspond to each value
      real(kind=RP)    , intent(in) :: Values(:)  ! Values of nonzero entries of matrix
      integer, optional, intent(in) :: num_of_Cols
      !-local-variables-----------------------------
      integer             :: istat
      !---------------------------------------------
      
      this % num_of_Rows = size(Rows) - 1
      
      if ( present(num_of_Cols) ) then
         this % num_of_Cols = num_of_Cols
      else
         this % num_of_Cols = this % num_of_Rows
      end if
      
      if ( maxval(Cols) > this % num_of_Cols) then
         print*, 'CSR_constructWithArrays :: WARNING: Increasing num_of_Cols'
         this % num_of_Cols = maxval(Cols)
      end if
      
!     Memory allocation
!     -----------------
      
      safedeallocate(this % Rows) ; allocate( this % Rows  ( size(Rows)   ), stat=istat )
      if ( istat .NE. 0 ) write(*,*) 'CSR_construct: Memory allocation error'
      
      safedeallocate(this % Diag) ; allocate( this % Diag  (this % num_of_Rows), stat=istat )
      if ( istat .NE. 0 ) write(*,*) 'CSR_construct: Memory allocation error'
      
      safedeallocate(this % Cols) ; allocate( this % Cols  ( size(Cols)   ), stat=istat )
      if ( istat .NE. 0 ) write(*,*) 'CSR_construct: Memory allocation error'
      
      safedeallocate(this % Values) ; allocate( this % Values( size(Values) ), stat=istat )
      if ( istat .NE. 0 ) write(*,*) 'CSR_construct: Memory allocation error'
      
      this % Rows   = Rows
      this % Cols   = Cols
      this % Values = Values
      
      call this % assigndiag()
      
   end subroutine CSR_constructWithArrays
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CSR_PreAllocate(this,nnz,nnzs, ForceDiagonal)
      implicit none
      !-----------------------------------
      CLASS(csrMat_t), INTENT(INOUT) :: this             !<> Matrix to be preallocated
      integer, optional, intent(in)  :: nnz        !< Num of nonzero entries all rows
      integer, optional, intent(in)  :: nnzs(:)     !< Num of nonzero entries in every row 
      logical, optional, intent(in)  :: ForceDiagonal
      !-----------------------------------
      integer :: i,k,istat
      integer :: total_nnz
      logical :: mustForceDiagonal
      !-----------------------------------
      
      if ( present(ForceDiagonal) ) then
         mustForceDiagonal = ForceDiagonal
      else
         mustForceDiagonal = .FALSE.
      end if
      
      if (present(nnz)) then
!
!        Constant number of non-zeros per row
!        ------------------------------------
         
         total_nnz = this % num_of_Rows * nnz
      elseif (present(nnzs)) then
!
!        Number of non-zeros different in every row
!        ------------------------------------------
         
         if (size(nnzs) /= this % num_of_Rows) ERROR stop ':: CSRMatrix: Not consistent nnzs'
         total_nnz = sum(nnzs)
      else
!
!        Unknown number of nonzeros per row
!        -> Preallocating with LinkedListMatrix
!        --------------------------------------
         
         this % usingListMat = .TRUE.
         call this % ListMatrix % construct(num_of_Rows = this % num_of_Rows)
         
         if (mustForceDiagonal) then
            do i = 1, this % num_of_Rows
               call this % ListMatrix % ForceAddToEntry(i,i,0._RP)
            end do
         end if
         return
      end if
      
      IF(total_nnz < 1) STOP ':: Invalid nnz' 
      
      k = this % num_of_Rows * total_nnz
       
      safedeallocate(this % Cols)   ; ALLOCATE( this % Cols(total_nnz),STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
       
      safedeallocate(this % Values) ; ALLOCATE( this % Values(total_nnz), STAT=istat )
      IF ( istat .NE. 0 ) WRITE(*,*) 'CSR_construct: Memory allocation error'
      
      this % Rows(1) = 1
      if (present(nnz)) then
         DO i=2, this % num_of_Rows + 1
            this % Rows(i) = this % Rows(i-1) + nnz
         END DO
      else
         DO i=2, this % num_of_Rows + 1
            this % Rows(i) = this % Rows(i-1) + nnzs(i-1)
         END DO
      end if
      
   end subroutine CSR_PreAllocate

   subroutine CSR_PreAllocateWithStructure(self, nnz, rows, cols, diag)
      class(CSRMat_t),  intent(inout)  :: self
      integer,          intent(in)     :: nnz
      integer,          intent(in)     :: rows(self % num_of_Rows+1)
      integer,          intent(in)     :: cols(nnz)
      integer,          intent(in)     :: diag(self % num_of_Rows)
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
      
      if (this % usingListMat) then
         call this % ListMatrix % Reset()
      else
         this % Cols = 0
         this % Diag = 0
         this % Values = 0._RP
      end if
   end subroutine CSR_Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  Subroutine to reset (=0._RP) all the entries of a block
!  -------------------------------------------------------
   subroutine CSR_ResetBlock(this, iBlock, jBlock )
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: iBlock, jBlock
      !---------------------------------------------
      integer :: row, col, pos
      !---------------------------------------------
      
      if ( (.not. allocated(this % BlockIdx)) .or. (.not. allocated(this % BlockSize)) ) then
         write(STD_OUT,*) 'CSRMatrixClass :: Error '
         write(STD_OUT,*) '               :: CSR_SetBlockEntry only available after CSR_SpecifyBlockInfo has been called'
         stop 99
      end if
      
      if (this % usingListMat) then
         do row = this % BlockIdx(iBlock), this % BlockIdx(iBlock) + this % BlockSize(iBlock) - 1
            do col = this % BlockIdx(jBlock), this % BlockIdx(jBlock) + this % BlockSize(jBlock) - 1
               call this % ListMatrix % ResetEntry(row,col)
            end do
         end do
      else
         ! TODO: This can be improved...
         do row = this % BlockIdx(iBlock), this % BlockIdx(iBlock) + this % BlockSize(iBlock) - 1
            do col = this % BlockIdx(jBlock), this % BlockIdx(jBlock) + this % BlockSize(jBlock) - 1
               
               pos = CSR_Search (this,row,col)
               if (pos == 0) cycle
               
               this % Values(pos) = 0._RP
               
            end do
         end do
      end if
      
      ! TODO: Do thisssss!!!
!~      if (this % usingListMat) then
!~         call this % ListMatrix % Reset()
!~      else
!~         this % Cols = 0
!~         this % Diag = 0
!~         this % Values = 0._RP
!~      end if
      
   end subroutine CSR_ResetBlock 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------------------
   SUBROUTINE CSR_AssignDiag(A)
      IMPLICIT NONE
      !-----------------------------------
      CLASS(csrMat_t) :: A           !<> Matrix
      !-----------------------------------
      integer         :: i, j       !   Counters
      !-----------------------------------
      do i=1, A % num_of_Rows
         
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
      integer        , INTENT(IN)    :: nvalues
      integer        , INTENT(IN)    :: irow(1:)    !< Different positions of Column
      integer        , INTENT(IN)    :: icol            !< Number of Row/Column
      REAL(KIND=RP)  , INTENT(IN)    :: values(:)         !< Values to be added to global matrivx¿x
      !------------------------------------------------------------------------------ 
      integer :: i !,k,l,c,
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
   subroutine CSR_SetEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !-local-variables-----------------------------
      integer                         :: k
      !---------------------------------------------
      
      if ( (row > this % num_of_Rows) .or. (col > this % num_of_Cols) ) then
         write (*,*) 'CSR_SetEntry: Dimension error. [row,col]=', row, col
         stop
      end if
      
      if (this % usingListMat) then
         call this % ListMatrix % SetEntry(row,col,value)
      else
         k = CSR_Search(this,row,col)
         this % Values(k) = value
      end if
      
   end subroutine CSR_SetEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Add a given value to an element of a CSR matrix
!  ---------------------------------------------
   subroutine CSR_AddToEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !-local-variables-----------------------------
      integer                         :: k
      !---------------------------------------------
      
      if ( (row > this % num_of_Rows) .or. (col > this % num_of_Cols) ) then
         write (*,*) 'CSR_SetEntry: Dimension error'
         stop
      end if
      
      if (this % usingListMat) then
         call this % ListMatrix % AddToEntry(row,col,value)
      else
         k = CSR_Search(this,row,col)
         this % Values(k) = this % Values(k) + value
      end if
      
   end subroutine CSR_AddToEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION CSR_Search (A,i,j) RESULT(k)
   !    Obtains the position k for the information of A(i,j) --> A % Values (k)
   !    If the position is not contained in the sparse matrix A, k = 0 is returned
   !       (supposing that A is ordered)
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------------------------ 
      TYPE(csrMat_t), INTENT(IN) :: A         !< Matrix to be searched
      integer       , INTENT(IN) :: i, j      !< Position to bi searched
      integer                    :: k         !> Position in A%Values of A(i,j)
      !------------------------------------------------------------
    
      IF(i .GT. A % num_of_Rows .OR. j .GT. A % num_of_Cols ) THEN
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
      integer                     :: n, nnz, i, fd
      LOGICAL                     :: First
      !------------------------------------------
      
      IF (.NOT. PRESENT(FirstRow)) First = .TRUE.
      
      n = this % num_of_Rows
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
      
      safedeallocate(this % Rows)
      safedeallocate(this % Cols)
      safedeallocate(this % Values)
      safedeallocate(this % Diag)
      safedeallocate(this % BlockIdx)
      safedeallocate(this % BlockIdx)
      safedeallocate(this % BlockSize)
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
      integer                        :: i
      !------------------------------------------
      
      if ( this % usingListMat ) then
         call this % ListMatrix % shift(shiftval)
      else
         DO i=1, this % num_of_Rows
            this % Values(this % Diag(i)) = this % Values(this % Diag(i)) + shiftval
         END DO 
      end if
   !----------------------------------------------------------------------------------
   END SUBROUTINE SetMatShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------
!  CSR Matrix addition:
!     Cmat = A + Factor*B
!  ----------------------
   function CSR_MatAdd(A,B,Factor) result(Cmat)
      implicit none
      !-arguments--------------------------------------------------------------------
      type(csrMat_t), intent(in)  :: A    !< Structure holding matrix
      type(csrMat_t), intent(in)  :: B    !< Structure holding matrix
      real(kind=RP) , intent(in)  :: Factor  !< 
      type(csrMat_t)              :: Cmat !< Structure holding matrix
      !-local-variables--------------------------------------------------------------
      real(kind=RP), allocatable :: c(:)
      integer,       allocatable :: ic(:), jc(:)
      integer                    :: info
      !------------------------------------------------------------------------------
#ifdef HAS_MKL
      
      allocate( ic(A % num_of_Rows+1), jc(A % num_of_Rows+1), c(A % num_of_Rows+1) ) ! For some reason, mkl_dcsradd needs jc and c to be allocated, even with request 1 - 2
      
      call mkl_dcsradd  ('n', 1, 4, A % num_of_Rows, A % num_of_Cols, &
                         A % Values, A % Cols, A % Rows, Factor, &
                         B % Values, B % Cols, B % Rows, c, jc, ic, 0, info)
      
      deallocate(jc,c)
      allocate( c(ic(A % num_of_Rows+1) - 1), jc(ic(A % num_of_Rows+1) - 1) )
      
      call mkl_dcsradd  ('n', 2, 4, A % num_of_Rows, A % num_of_Cols, &
                         A % Values, A % Cols, A % Rows, Factor, &
                         B % Values, B % Cols, B % Rows, c, jc, ic, 0, info)
      
      if (info /= 0) then
         print*, 'ERROR :: CSR_MatAdd stopped at line', info
         stop
      end if
                              
      call Cmat % constructWithArrays (ic,jc,c, A % num_of_Cols)
#else
      stop ':: CSR_MatAdd needs MKL'
#endif
   end function CSR_MatAdd
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Matrix-matrix multiplication for CSR matrices
!     Cmat = A * B
!  ---------------------------------------------
   function CSR_MatMatMul(A,B,trans) result(Cmat)
      implicit none
      !-arguments--------------------------------------------------------------------
      type(csrMat_t)   , intent(in) :: A       !< Structure holding matrix
      type(csrMat_t)   , intent(in) :: B       !< Structure holding matrix
      logical, optional, intent(in) :: trans   !< A matrix is transposed?
      type(csrMat_t)                :: Cmat    !< Structure holding matrix
      !-local-variables--------------------------------------------------------------
      real(kind=RP), allocatable :: c(:)
      integer,       allocatable :: ic(:), jc(:)
      integer                    :: info
      character(len=1)           :: transInfo
      !------------------------------------------------------------------------------
      
      if (A % num_of_Cols /= B % num_Of_Rows) then
         write(STD_OUT,'(A,I0,A,I0,A,I0,A,I0,A)') 'CSR_MatMatMul :: ERROR: Matrix dimensions mismatch: A(', A % num_of_Rows,',', A % num_Of_Cols, ') ; B(', B % num_Of_Rows, ',', B % num_Of_Cols,')'
         stop
      end if
      
      if ( present(trans) ) then
         if (trans) then
            transInfo = 't'
         else
            transInfo = 'n'
         end if
      else
         transInfo = 'n'
      end if
      
#ifdef HAS_MKL
      
      allocate( ic(A % num_of_Rows+1), jc(A % num_of_Rows+1), c(A % num_of_Rows+1) ) ! For some reason, mkl_dcsrmultcsr needs jc and c to be allocated, even with request 1 - 2
      
      call mkl_dcsrmultcsr  (transInfo, 1, 8, A % num_of_Rows, A % num_of_Cols, B % num_of_Cols, &
                              A % Values, A % Cols, A % Rows, &
                              B % Values, B % Cols, B % Rows, &
                              c, jc, ic, 0, info)
      
      deallocate(jc,c)
      allocate( c(ic(A % num_of_Rows+1) - 1), jc(ic(A % num_of_Rows+1) - 1) )
      
      call mkl_dcsrmultcsr  (transInfo, 2, 8, A % num_of_Rows, A % num_of_Cols, B % num_of_Cols, &
                              A % Values, A % Cols, A % Rows, &
                              B % Values, B % Cols, B % Rows, &
                              c, jc, ic, 0, info)
      
      if (info /= 0) then
         print*, 'ERROR :: mkl_dcsrmultcsr stopped at line', info
         stop
      end if
      
      call Cmat % constructWithArrays (ic,jc,c, B % num_of_Cols)
#else
      stop ':: CSR_MatMatMul needs MKL'
#endif
   end function CSR_MatMatMul
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION CSR_MatVecMul( A,u, trans) RESULT(v)
   ! Matrix vector product (v = Au) for a CSR matrix .... v needs to be allocated beforehand
   !------------------------------------------------------------------------------
      type(csrMat_t)   , intent(in)  :: A  !< Structure holding matrix
      real(kind=RP)    , intent(in)  :: u(A % num_of_Cols)  !< Vector to be multiplied
      logical, optional, intent(in) :: trans   !< A matrix is transposed?
      real(kind=RP)               :: v(A % num_of_Rows)  !> Result vector 
      !------------------------------------------------------------------------------
      integer           :: i,j
      REAL(KIND=RP)     :: rsum
      character(len=1)  :: transInfo
      !------------------------------------------------------------------------------
    
      IF (A % num_of_Cols .NE. SIZE(U) .OR. A % num_of_Rows .NE. SIZE(v)) THEN
         STOP 'CSR_MatVecMul: Error - u dimensions mismatch'
      END IF
      
      if ( present(trans) ) then
         if (trans) then
            transInfo = 't'
         else
            transInfo = 'n'
         end if
      else
         transInfo = 'n'
      end if
    
#ifdef HAS_MKL
      CALL mkl_dcsrgemv(transInfo, A % num_of_Rows, A % Values, A % Rows, A % Cols, u, v)
#else
      if (transInfo = 't') ERROR stop "CSR_MatVecMul with 't' only with MKL"
!$omp parallel do private(j,rsum)
      DO i=1,A % num_of_Rows
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
      integer                                 :: i, k       !  Counters
      !------------------------------------------------
      
      IF (ALLOCATED(Mat)) DEALLOCATE(Mat)
      ALLOCATE(Mat(A%num_of_Rows,A%num_of_Cols))
      
      Mat = 0.d0
      DO i=1, A % num_of_Rows
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
      integer, INTENT(IN)       :: Num        !< Number of diagonal element
      integer, INTENT(IN)       :: N          !< Size of block
      REAL(KIND=RP)             :: Mat(N,N)   !> Value of the matrix element
      !------------------------------------------------------------------------------ 
      integer                   :: RC0      ! Index of first row/column of block in global sparse matrix
      integer                   :: RCf      ! Index of last row/column of block in global sparse matrix
      integer                   :: i, j     ! Row/Column index of output matrix
      integer                   :: ii, jj   ! Row/Column index in global sparse matrix
      integer                   :: k, l     ! Variables containing positions in sparse matrix arrays
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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CSR_Assembly(this)
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout)    :: this
      !---------------------------------------------
      
      if ( this % usingListMat ) then
         call this % ListMatrix % getCSRarrays(this % Values, this % Cols, this % Rows)
         call this % AssignDiag()
         call this % ListMatrix % destruct()
         this % usingListMat = .FALSE.
      end if
   end subroutine CSR_Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CSR_SpecifyBlockInfo(this,BlockIdx,BlockSize)
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: BlockIdx(:)
      integer        , intent(in)    :: BlockSize(:)
      !---------------------------------------------
      
      safedeallocate(this % BlockIdx)  ; allocate (this % BlockIdx (size(BlockIdx )) )
      safedeallocate(this % BlockSize) ; allocate (this % BlockSize(size(BlockSize)) )
      this % BlockIdx  = BlockIdx
      this % BlockSize = BlockSize
   end subroutine CSR_SpecifyBlockInfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  --------------------------------------------------------------
   subroutine CSR_SetBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: iBlock, jBlock
      integer        , intent(in)    :: i, j
      real(kind=RP)  , intent(in)    :: value
      !-local-variables-----------------------------
      integer :: row, col
      !---------------------------------------------
      
      if (.not. allocated(this % BlockIdx)) then
         write(STD_OUT,*) 'CSRMatrixClass :: Error '
         write(STD_OUT,*) '               :: CSR_SetBlockEntry only available after CSR_SpecifyBlockInfo has been called'
         stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % SetEntry(row, col, value)
      
   end subroutine CSR_SetBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to add a value to the entries of a block with relative index
!  -----------------------------------------------------------------------
   subroutine CSR_AddToBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(csrMat_t), intent(inout) :: this
      integer        , intent(in)    :: iBlock, jBlock
      integer        , intent(in)    :: i, j
      real(kind=RP)  , intent(in)    :: value
      !-local-variables-----------------------------
      integer :: row, col
      !---------------------------------------------
      
      if (.not. allocated(this % BlockIdx)) then
         write(STD_OUT,*) 'CSRMatrixClass :: Error '
         write(STD_OUT,*) '               :: CSR_AddToBlockEntry only available after CSR_SpecifyBlockInfo has been called'
         stop 99
      end if
      
      row = this % BlockIdx(iBlock) + i - 1
      col = this % BlockIdx(jBlock) + j - 1
      
      call this % AddToEntry(row, col, value)
      
   end subroutine CSR_AddToBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE CSRMatrixClass
