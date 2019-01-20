!
!////////////////////////////////////////////////////////////////////////
!
!      DenseBlockDiagonalMatrixClass.f90
!      Created: 2018-02-19 17:07:00 +0100 
!      By: Andrés Rueda
!
!      Class for sparse block diagonal matrices
!      -> The matrix is not dense, the block is!
!
!////////////////////////////////////////////////////////////////////////
module DenseBlockDiagonalMatrixClass
   use SMConstants
   use GenericMatrixClass
   use CSRMatrixClass, only: csrMat_t
#include "Includes.h"
   implicit none
   
   private
   public DenseBlockDiagMatrix_t, Matrix_t
   
   type Block_t
      real(kind=RP), allocatable :: Matrix(:,:)
      integer      , allocatable :: Indexes(:)
   end type Block_t
   
   type, extends(Matrix_t) :: DenseBlockDiagMatrix_t
      type(Block_t), allocatable :: Blocks(:)   ! Array containing each block in a dense matrix
      integer                    :: NumOfBlocks ! Number of blocks in matrix
      integer      , allocatable :: BlockSizes(:)
      integer      , allocatable :: BlockIdx(:)
      contains
         procedure :: construct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetColumn
         procedure :: SetDiagonalBlock
         procedure :: SetBlockEntry
         procedure :: AddToBlockEntry
         procedure :: ResetBlock
         procedure :: Assembly
         procedure :: SpecifyBlockInfo
         procedure :: shift
         procedure :: destruct
         procedure :: FactorizeBlocks_LU
         procedure :: SolveBlocks_LU
         procedure :: InvertBlocks_LU
         procedure :: getCSR
   end type DenseBlockDiagMatrix_t
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,num_of_rows_reduced,withMPI)
      implicit none
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t) :: this     !<> This matrix
      integer, optional, intent(in) :: num_of_Rows
      integer, optional, intent(in) :: num_of_Cols
      integer, optional, intent(in) :: num_of_Blocks
      integer, optional, intent(in) :: num_of_rows_reduced
      logical, optional, intent(in) :: WithMPI
      !---------------------------------------------
      
      if ( .not. present(num_of_Blocks) ) then
         ERROR stop 'DenseBlockDiagMatrix_t needs num_of_Blocks'
      end if
      
      allocate ( this % Blocks(num_of_Blocks) )
      this % NumOfBlocks = num_of_Blocks
      allocate ( this % BlockSizes(num_of_Blocks) )
      allocate ( this % BlockIdx(num_of_Blocks+1) )
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs, ForceDiagonal)
      IMPLICIT NONE
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this    !<> This matrix
      integer, optional            , intent(in)    :: nnz     !<  Not needed here
      integer, optional            , intent(in)    :: nnzs(:) !<  nnzs contains the block sizes!
      logical, optional, intent(in)  :: ForceDiagonal
      !---------------------------------------------
      integer :: i, k ! counters
      !---------------------------------------------
      
      if (.not. present(nnzs) ) ERROR stop ':: DenseBlockDiagMatrix needs the block sizes'
      if ( size(nnzs) /= this % NumOfBlocks) ERROR stop ':: DenseBlockDiagMatrix: wrong dimension for the block sizes'
      
      this % BlockSizes = nnzs
      this % num_of_Rows = sum(nnzs)
      
      this % BlockIdx(1) = 1
      do i=2, this % NumOfBlocks + 1
         this % BlockIdx(i) = this % BlockIdx(i-1) + nnzs(i-1)
      end do
      
!$omp parallel do private(k) schedule(runtime)
      do i=1, this % NumOfBlocks
         safedeallocate (this % Blocks(i) % Matrix ) ; allocate ( this % Blocks(i) % Matrix(nnzs(i),nnzs(i)) )
         safedeallocate (this % Blocks(i) % Indexes) ; allocate ( this % Blocks(i) % Indexes(nnzs(i)) )
         
         this % Blocks(i) % Indexes = (/ (k, k=this % BlockIdx(i),this % BlockIdx(i+1) - 1 ) /)
         
      end do
!$omp end parallel do      
   end subroutine Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this)
      IMPLICIT NONE
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this     !<> This matrix
      !---------------------------------------------
      integer :: i
      !---------------------------------------------
      
      do i=1, this % NumOfBlocks
         this % Blocks(i) % Matrix = 0._RP
      end do
      
   end subroutine Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetColumn(this,nvalues, irow, icol, values )
      implicit none
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: nvalues
      integer, dimension(:)        , intent(in)    :: irow
      integer                      , intent(in)    :: icol
      real(kind=RP), dimension(:)  , intent(in)    :: values
      !---------------------------------------------
      integer :: thisblock, thiscol, thisrow, firstIdx, lastIdx
      integer :: i
!~      integer, pointer :: indexes(:)
      !---------------------------------------------
      
      if ( (icol > this % num_of_Rows) .or. (icol < 1) ) ERROR stop ':: DenseBlockDiagMatrix: icol value is out of bounds'
      
      ! Search the corresponding block (they are ordered)
      do thisblock=1, this % NumOfBlocks
         if (icol <= this % BlockIdx(thisblock+1) -1) exit
      end do
      
      associate (indexes => this % Blocks(thisblock) % Indexes)
      firstIdx = this % BlockIdx(thisblock)
      lastIdx  = this % BlockIdx(thisblock+1) - 1
      
      ! Get relative position of column
      do thiscol=1, this % BlockSizes(thisblock)
         if (icol == indexes(thiscol)) exit
      end do
      
      ! Fill the column info
      do i=1, nvalues
         if ( irow(i) < firstIdx .or. irow(i) > lastIdx ) cycle
         ! Get relative row
         do thisrow=1, this % BlockSizes(thisblock)
            if (irow(i) == indexes(thisrow)) exit
         end do
         this % Blocks(thisblock) % Matrix(thisrow,thiscol) = values(i)
      
      end do
      
      end associate
      
   end subroutine SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetDiagonalBlock(this,BlockNum, values )
      implicit none
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: BlockNum
      real(kind=RP), dimension(:,:), intent(in)    :: values
      !---------------------------------------------
      
      if ( size(values,1) /= this % BlockSizes(BlockNum) .or. &
           size(values,2) /= this % BlockSizes(BlockNum) ) ERROR stop ':: DenseBlockDiagMatrix_t % DenseBlockDiagMatrix_t. Block size is not consistent'
      
      this % Blocks(BlockNum) % Matrix = values
      
   end subroutine SetDiagonalBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  --------------------------------------------------------------
   subroutine SetBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: iBlock, jBlock
      integer                      , intent(in)    :: i, j
      real(kind=RP)                , intent(in)    :: value
      !---------------------------------------------
      
      if (iBlock /= jBlock) return ! Only diagonal blocks here!
      
      this % Blocks(iBlock) % Matrix(i,j) = value
      
   end subroutine SetBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to add a value to the entries of a block with relative index
!  -----------------------------------------------------------------------
   subroutine AddToBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: iBlock, jBlock
      integer                      , intent(in)    :: i, j
      real(kind=RP)                , intent(in)    :: value
      !---------------------------------------------
      
      if (iBlock /= jBlock) return ! Only diagonal blocks here!
      
      this % Blocks(iBlock) % Matrix(i,j) = this % Blocks(iBlock) % Matrix(i,j) + value
      
   end subroutine AddToBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  Subroutine to reset (=0._RP) all the entries of a block
!  -------------------------------------------------------
   subroutine ResetBlock(this, iBlock, jBlock )
      implicit none
      !-arguments-----------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: iBlock, jBlock
      !---------------------------------------------
      
      if (iBlock /= jBlock) return ! Only diagonal blocks here!
      
      this % Blocks(iBlock) % Matrix = 0._RP
      
   end subroutine ResetBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Assembly(this)
      implicit none
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      !---------------------------------------------
      
      ! Do nothing
      
   end subroutine Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SpecifyBlockInfo(this,BlockIdx,BlockSize)
      implicit none
      !-arguments-----------------------------------
      class(DenseBlockDiagMatrix_t) , intent(inout) :: this
      integer                       , intent(in)    :: BlockIdx(:)
      integer                       , intent(in)    :: BlockSize(:)
      !---------------------------------------------
      ! Do nothing
      ! Currently not needed since this info is provided at construction time
   end subroutine SpecifyBlockInfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine shift(this,shiftval)
      implicit none
      !------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this
      real(kind=RP)                , intent(in)    :: shiftval
      !------------------------------------------
      INTEGER                :: i, iBL
      !------------------------------------------
      
!$omp parallel do private(i) schedule(runtime)
      do iBL=1, this % NumOfBlocks
         do i=1, size(this % Blocks(iBL) % Matrix,1)
            this % Blocks(iBL) % Matrix(i,i) = this % Blocks(iBL) % Matrix(i,i) + shiftval
         end do
      end do
!$omp end parallel do
      
   end subroutine shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t),     intent(inout)     :: this
      !---------------------------------------------
      integer :: i
      !---------------------------------------------
      
      do i = 1, this % NumOfBlocks
         deallocate (this % Blocks(i) % Matrix )
         deallocate (this % Blocks(i) % Indexes)
      end do
      
      deallocate (this % Blocks)
      deallocate (this % BlockSizes)
      deallocate (this % BlockIdx)
      
   end subroutine destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  Factorizes the blocks of the matrix as A=LU
!  -> Stores L+U in Factorized % Blocks(:) % Matrix
!  -> Stores the pivots in Factorized % Blocks(:) % Indexes
!  --------------------------------------------------------
   subroutine FactorizeBlocks_LU(this,Factorized)
      use DenseMatUtilities
      implicit none
      !-------------------------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(in)    :: this            !<  This matrix
      class(Matrix_t), intent(inout) :: Factorized      !<  Facorized matrix
      !-------------------------------------------------------------
      integer :: k      ! Counter
      !-------------------------------------------------------------
      
      select type (Factorized)
         class is(DenseBlockDiagMatrix_t)
!$omp parallel do schedule(runtime)
            do k=1, this % NumOfBlocks
               call ComputeLU (A        = this       % Blocks(k) % Matrix, &
                               ALU      = Factorized % Blocks(k) % Matrix, &
                               LUpivots = Factorized % Blocks(k) % Indexes)
            end do
!$omp end parallel do
         class default
            write(STD_OUT,*) 'DenseBlockDiagonalMatrixClass :: Wrong type For factorized matrix in FactorizeBlocks_LU'
            stop
      end select
   end subroutine FactorizeBlocks_LU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  Solves a system of the for Ax=b using LU factorization in each block
!  -> To be used with a matrix factorized using FactorizeBlocks_LU (Factorized)
!  -> 
!  --------------------------------------------------------
   subroutine SolveBlocks_LU(this,x,b)
      use DenseMatUtilities
      implicit none
      !-------------------------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(in)    :: this                 !<  FACTORIZED matrix for solving the problem
      real(kind=RP)                , intent(in)    :: b(this % num_of_Rows)    !<  RHS
      real(kind=RP)                , intent(inout) :: x(this % num_of_Rows)    !<  Solution
      !-------------------------------------------------------------
      integer                    :: k        ! Counter
      real(kind=RP), allocatable :: x_loc(:) ! Local x
      !-------------------------------------------------------------
      
!$omp parallel do private(x_loc) schedule(runtime)
      do k = 1, this % NumOfBlocks
         allocate( x_loc(this % BlockSizes(k)) )
         call SolveLU(ALU      = this % Blocks(k) % Matrix, &
                      LUpivots = this % Blocks(k) % Indexes, &
                      x = x_loc, &
                      b = b ( this % BlockIdx(k):this % BlockIdx(k+1)-1 ) )
!$omp critical
         x( this % BlockIdx(k):this % BlockIdx(k+1)-1 ) = x_loc
!$omp end critical
         deallocate(x_loc)
      end do
!$omp end parallel do
      
   end subroutine SolveBlocks_LU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  Inverts the blocks of a matrix using LU factorization
!  --------------------------------------------------------
   subroutine InvertBlocks_LU(this,Inverted)
      use DenseMatUtilities
      implicit none
      !-arguments---------------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(in)    :: this            !<  This matrix
      class(Matrix_t)              , intent(inout) :: Inverted      !<  Facorized matrix
      !-local-variables---------------------------------------------
      integer :: k      ! Counter
      !-------------------------------------------------------------
      
      select type (Inverted)
         class is(DenseBlockDiagMatrix_t)
!$omp parallel do schedule(runtime)
            do k=1, this % NumOfBlocks
               Inverted % Blocks(k) % Matrix = inverse (this % Blocks(k) % Matrix)
            end do
!$omp end parallel do
         class default
            write(STD_OUT,*) 'DenseBlockDiagonalMatrixClass :: Wrong type For factorized matrix in FactorizeBlocks_LU'
            stop
      end select
   end subroutine InvertBlocks_LU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Convert the DBD matrix to a CSR matrix
!  --------------------------------------
   subroutine getCSR(this,Acsr)
      implicit none
      !-arguments---------------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(in)    :: this          !<  This matrix
      class(csrMat_t)              , intent(inout) :: Acsr      !<  Facorized matrix
      !-local-variables---------------------------------------------
      integer :: ii, jj
      integer :: bID
      !-------------------------------------------------------------
      
      if (this % num_of_Rows /= Acsr % num_of_Rows) then
         print*, 'DBD_getCSR :: ERROR: Matrix dimensions mismatch:', this % num_of_Rows, Acsr % num_of_Rows
         stop
      end if
      
      call Acsr % PreAllocate()
      call Acsr % Reset
      
      call Acsr % SpecifyBlockInfo(this % BlockIdx, this % BlockSizes)
      
      
!     Fill the Matrix
!     ---------------
      
!$omp parallel do private(ii,jj)
      do bID=1, this % NumOfBlocks
               
         do jj=1, this % BlockSizes(bID)
            do ii=1, this % BlockSizes(bID)
                  call Acsr % SetBlockEntry(bID,bID,ii,jj, this % Blocks (bID) % Matrix(ii,jj))
            end do
         end do
         
      end do
!$omp end parallel do
      
      call Acsr % assembly()
      
   end subroutine getCSR   
end module DenseBlockDiagonalMatrixClass
