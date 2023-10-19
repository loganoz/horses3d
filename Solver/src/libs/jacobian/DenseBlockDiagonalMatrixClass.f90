!
!////////////////////////////////////////////////////////////////////////
!
!      Class for sparse block diagonal matrices
!      -> The matrix is not dense, the block is!
!
!////////////////////////////////////////////////////////////////////////
module DenseBlockDiagonalMatrixClass
   use SMConstants
   use GenericMatrixClass  , only: Matrix_t, DenseBlock_t
   use CSRMatrixClass      , only: csrMat_t
   use JacobianDefinitions , only: JACEPS
   use PartitionedMeshClass, only: mpi_partition ! for MPI
   use MPI_Process_Info    , only: MPI_Process
#include "Includes.h"
   implicit none
   
   private
   public DenseBlockDiagMatrix_t, Matrix_t
   
   type, extends(Matrix_t) :: DenseBlockDiagMatrix_t
      type(DenseBlock_t), allocatable :: Blocks(:)   ! Array containing each block in a dense matrix
      contains
         procedure :: construct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetColumn           => DBD_SetColumn
         procedure :: AddToColumn         => DBD_AddToColumn
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
         procedure :: getTransCSR
         procedure :: Visualize
   end type DenseBlockDiagMatrix_t
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,num_of_TotalRows,withMPI)
      implicit none
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t) :: this     !<> This matrix
      integer, optional, intent(in) :: num_of_Rows
      integer, optional, intent(in) :: num_of_Cols
      integer, optional, intent(in) :: num_of_Blocks
      integer, optional, intent(in) :: num_of_TotalRows
      logical, optional, intent(in) :: WithMPI
      !---------------------------------------------
      
      if ( .not. present(num_of_Blocks) ) then
         error stop 'DenseBlockDiagMatrix_t needs num_of_Blocks'
      end if
      
      allocate ( this % Blocks(num_of_Blocks) )
      this % num_of_Blocks = num_of_Blocks
      allocate ( this % BlockSizes(num_of_Blocks) )
      allocate ( this % BlockIdx(num_of_Blocks+1) )
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs)
      IMPLICIT NONE
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this    !<> This matrix
      integer, optional            , intent(in)    :: nnz     !<  Not needed here
      integer, optional            , intent(in)    :: nnzs(:) !<  nnzs contains the block sizes!
      !---------------------------------------------
      integer :: i, k ! counters
      !---------------------------------------------
      
      if (.not. present(nnzs) ) error stop ':: DenseBlockDiagMatrix needs the block sizes'
      if ( size(nnzs) /= this % num_of_Blocks) error stop ':: DenseBlockDiagMatrix: wrong dimension for the block sizes'
      
      this % BlockSizes = nnzs
      this % num_of_Rows = sum(nnzs)
      
      this % BlockIdx(1) = 1
      do i=2, this % num_of_Blocks + 1
         this % BlockIdx(i) = this % BlockIdx(i-1) + nnzs(i-1)
      end do
      
!$omp parallel do private(k) schedule(runtime)
      do i=1, this % num_of_Blocks
         safedeallocate (this % Blocks(i) % Matrix ) ; allocate ( this % Blocks(i) % Matrix(nnzs(i),nnzs(i)) )
         safedeallocate (this % Blocks(i) % Indexes) ; allocate ( this % Blocks(i) % Indexes(nnzs(i)) )
         
         this % Blocks(i) % Indexes = (/ (k, k=this % BlockIdx(i),this % BlockIdx(i+1) - 1 ) /)
         
      end do
!$omp end parallel do      
   end subroutine Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this, ForceDiagonal)
      IMPLICIT NONE
      !---------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(inout) :: this     !<> This matrix
      logical, optional, intent(in)  :: ForceDiagonal
      !---------------------------------------------
      integer :: i
      !---------------------------------------------
      
      do i=1, this % num_of_Blocks
         this % Blocks(i) % Matrix = 0._RP
      end do
      
   end subroutine Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine DBD_SetColumn(this,nvalues, irow, icol, values )
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
      
      if (MPI_Process % doMPIAction) error stop 'DBD_SetColumn not ready for MPI'
      
      if ( (icol > this % num_of_Rows) .or. (icol < 1) ) error stop ':: DenseBlockDiagMatrix: icol value is out of bounds'
      
      ! Search the corresponding block (they are ordered)
      do thisblock=1, this % num_of_Blocks
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
      
   end subroutine DBD_SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine DBD_AddToColumn(this,nvalues, irow, icol, values )
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
      
      if (MPI_Process % doMPIAction) error stop 'DBD_AddToColumn not ready for MPI'
      
      if ( (icol > this % num_of_Rows) .or. (icol < 1) ) error stop ':: DenseBlockDiagMatrix: icol value is out of bounds'
      
      ! Search the corresponding block (they are ordered)
      do thisblock=1, this % num_of_Blocks
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
         this % Blocks(thisblock) % Matrix(thisrow,thiscol) = this % Blocks(thisblock) % Matrix(thisrow,thiscol) + values(i)
      
      end do
      
      end associate
      
   end subroutine DBD_AddToColumn
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
           size(values,2) /= this % BlockSizes(BlockNum) ) error stop ':: DenseBlockDiagMatrix_t % DenseBlockDiagMatrix_t. Block size is not consistent'
      
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
      
      this % Blocks( mpi_partition % global2localeID(iBlock) ) % Matrix(i,j) = value
      
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
      
      this % Blocks( mpi_partition % global2localeID(iBlock) ) % Matrix(i,j) = this % Blocks( mpi_partition % global2localeID(iBlock) ) % Matrix(i,j) + value
      
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
      
      this % Blocks( mpi_partition % global2localeID(iBlock) ) % Matrix = 0._RP
      
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
      do iBL=1, this % num_of_Blocks
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
      
      do i = 1, this % num_of_Blocks
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
            do k=1, this % num_of_Blocks
               call ComputeLU (A        = this       % Blocks(k) % Matrix, &
                               ALU      = Factorized % Blocks(k) % Matrix, &
                               LUpivots = Factorized % Blocks(k) % Indexes)
            end do
!$omp end parallel do
         class default
            write(STD_OUT,*) 'DenseBlockDiagonalMatrixClass :: Wrong type For factorized matrix in FactorizeBlocks_LU'
            error stop
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
      do k = 1, this % num_of_Blocks
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
            do k=1, this % num_of_Blocks
               Inverted % Blocks(k) % Matrix = inverse (this % Blocks(k) % Matrix)
            end do
!$omp end parallel do
         class default
            write(STD_OUT,*) 'DenseBlockDiagonalMatrixClass :: Wrong type For factorized matrix in FactorizeBlocks_LU'
            error stop
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
         error stop
      end if
      
      call Acsr % PreAllocate()
      call Acsr % Reset
      
      call Acsr % SpecifyBlockInfo(this % BlockIdx, this % BlockSizes)
      
      
!     Fill the Matrix
!     ---------------
      
!$omp parallel do private(ii,jj)
      do bID=1, this % num_of_Blocks
               
         do jj=1, this % BlockSizes(bID)
            do ii=1, this % BlockSizes(bID)
                  call Acsr % SetBlockEntry(bID,bID,ii,jj, this % Blocks (bID) % Matrix(ii,jj))
            end do
         end do
         
      end do
!$omp end parallel do
      
      call Acsr % assembly()
      
   end subroutine getCSR
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Convert the DBD matrix to a CSR matrix
!  --------------------------------------
   subroutine getTransCSR(this,Acsr)
      implicit none
      !-arguments---------------------------------------------------
      class(DenseBlockDiagMatrix_t), intent(in)    :: this          !<  This matrix
      class(csrMat_t)              , intent(inout) :: Acsr      !<  Facorized matrix
      !-local-variables---------------------------------------------
      integer :: k, ii, jj, j_offset, i
      integer :: bID, rowsize, nnz_0, nnz
      real(kind=RP), allocatable :: Vals(:)
      integer      , allocatable :: Cols(:), Rows(:)
      !-------------------------------------------------------------
      
      if (this % num_of_Rows /= Acsr % num_of_Rows) then
         print*, 'DBD_getCSR :: ERROR: Matrix dimensions mismatch:', this % num_of_Rows, Acsr % num_of_Rows
         error stop
      end if
      
      nnz_0 = sum(this % BlockSizes**2)
      
      allocate ( Rows(this % num_of_Rows+1) )
      allocate ( Cols(nnz_0), Vals(nnz_0) )
      
      
      
      Rows(1) = 1
      
      i=1
      k=1
      j_offset = 0
      nnz = 0
      do bID = 1, this % num_of_Blocks
         do jj = 1, this % BlockSizes(bID)
            rowsize = 0
            do ii = 1, this % BlockSizes(bID) 
               if (abs(this % Blocks(bID) % Matrix(ii,jj)) < JACEPS) cycle
               
               Vals(k) = this % Blocks(bID) % Matrix(ii,jj)
               Cols(k) = ii + j_offset
               k = k + 1
               rowsize = rowsize + 1
               nnz = nnz + 1
            end do
            Rows(i+1) = Rows(i) + rowsize
            i = i + 1
         end do
         
         j_offset = j_offset + this % BlockSizes(bID)
      end do
      
      call Acsr % constructWithCSRArrays( Rows, Cols(1:nnz), Vals(1:nnz), this % num_of_Rows )
      
   end subroutine getTransCSR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Visualize(this, filename, FirstRow)
!  ---------------------------------------------------------
!  Save Jacobian matrix in the ASCII file. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(DenseBlockDiagMatrix_t)  , intent(in) :: this
      character(len=*)               , intent(in) :: filename
      logical, optional              , intent(in) :: FirstRow   !< Write First row?
!-----Local-Variables-----------------------------------------------------
      class(csrMat_t) , allocatable :: Acsr
!-------------------------------------------------------------------------

      allocate(Acsr)
      call Acsr % Construct(num_of_Rows = this % num_of_Rows)
      call this % getCSR(Acsr)
      if (present(FirstRow)) then
         call Acsr % Visualize(trim(filename),FirstRow)
      else
         call Acsr % Visualize(trim(filename))
      end if
      call Acsr % destruct
      write(STD_OUT,'(20X,A,A)') "Jacobian saved to: ",filename

   end subroutine Visualize
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module DenseBlockDiagonalMatrixClass