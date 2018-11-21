module GenericMatrixClass
   use SMConstants
   implicit none
   
   private
   public Matrix_t
   
   type Matrix_t
      integer :: NumRows
      contains
         procedure :: construct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetColumn
         procedure :: AddToColumn
         procedure :: SetDiagonalBlock
         procedure :: SetEntry
         procedure :: AddToEntry
         procedure :: SetBlockEntry
         procedure :: AddToBlockEntry
         procedure :: ResetBlock
         procedure :: Shift
         procedure :: ReShift
         procedure :: PreAssembly
         procedure :: Assembly
         procedure :: SpecifyBlockInfo
         procedure :: destruct
         procedure :: SolveBlocks_LU
         procedure :: FactorizeBlocks_LU
   end type Matrix_t
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,dimPrb,withMPI)
      implicit none
      !---------------------------------------------
      class(Matrix_t) :: this
      integer, intent(in)           :: dimPrb
      logical, optional, intent(in) :: WithMPI
      !---------------------------------------------
      ERROR stop ' :: construct not implemented for current matrix type'
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs)
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout) :: this
      integer, optional, intent(in)  :: nnz
      integer, optional, intent(in)  :: nnzs(:)
      !---------------------------------------------
      ERROR stop ' :: Preallocate not implemented for current matrix type'
   end subroutine Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)     :: this
      !---------------------------------------------
      ERROR stop ' :: Reset not implemented for current matrix type'
   end subroutine Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetColumn(this,nvalues, irow, icol, values )
      implicit none
      !---------------------------------------------
      class(Matrix_t)            , intent(inout) :: this
      integer                    , intent(in)    :: nvalues
      integer, dimension(:)      , intent(in)    :: irow
      integer                    , intent(in)    :: icol
      real(kind=RP), dimension(:), intent(in)    :: values
      !---------------------------------------------
      ERROR stop ' :: SetColumn not implemented for current matrix type'
   end subroutine SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   subroutine AddToColumn(this,nvalues, irow, icol, values )
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout)         :: this
      integer, intent(in)                               :: nvalues
      integer, dimension(:), intent(in)                 :: irow
      integer, intent(in)                               :: icol
      real(kind=RP) , dimension(:), intent(in)                 :: values
      !---------------------------------------------
      ERROR stop ' :: AddToColumn not implemented for current matrix type'
   end subroutine AddToColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetDiagonalBlock(this,BlockNum, values )
      implicit none
      !---------------------------------------------
      class(Matrix_t)              , intent(inout) :: this
      integer                      , intent(in)    :: BlockNum
      real(kind=RP), dimension(:,:), intent(in)    :: values
      !---------------------------------------------
      ERROR stop ' :: SetDiagonalBlock not implemented for current matrix type'
   end subroutine SetDiagonalBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   subroutine SetEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(Matrix_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !---------------------------------------------
      ERROR stop ' :: SetEntry not implemented for current matrix type'
   end subroutine SetEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   subroutine AddToEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(Matrix_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !---------------------------------------------
      ERROR stop ' :: AddToEntry not implemented for current matrix type'
   end subroutine AddToEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  ------------------------------------------------------------
   subroutine SetBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(Matrix_t), intent(inout) :: this
      integer        , intent(in)    :: iBlock, jBlock
      integer        , intent(in)    :: i, j
      real(kind=RP)  , intent(in)    :: value
      !---------------------------------------------
      ERROR stop ' :: SetBlockEntry not implemented for current matrix type'
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
      class(Matrix_t), intent(inout) :: this
      integer        , intent(in)    :: iBlock, jBlock
      integer        , intent(in)    :: i, j
      real(kind=RP)  , intent(in)    :: value
      !---------------------------------------------
      ERROR stop ' :: AddToBlockEntry not implemented for current matrix type'
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
      class(Matrix_t), intent(inout) :: this
      integer        , intent(in)    :: iBlock, jBlock
      !---------------------------------------------
      ERROR stop ' :: ResetBlock not implemented for current matrix type'
   end subroutine ResetBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Shift(this, shiftval)
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout)     :: this
      real(kind=RP),                     intent(in)        :: shiftval
      !---------------------------------------------
      ERROR stop ' :: Shift not implemented for current matrix type'
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
      class(Matrix_t), intent(inout)     :: this
      real(kind=RP),               intent(in)        :: shiftval
      !---------------------------------------------
      ERROR stop ' :: ReShift not implemented for current matrix type'
   end subroutine ReShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PreAssembly(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)   :: this
      !---------------------------------------------
      ERROR stop ' :: PreAssembly not implemented for current matrix type'
   end subroutine PreAssembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Assembly(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout) :: this
      !---------------------------------------------
      ERROR stop ' :: Assembly not implemented for current matrix type'
   end subroutine Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SpecifyBlockInfo(this,BlockIdx,BlockSize)
      implicit none
      !-arguments-----------------------------------
      class(Matrix_t), intent(inout) :: this
      integer        , intent(in)    :: BlockIdx(:)
      integer        , intent(in)    :: BlockSize(:)
      !---------------------------------------------
      ERROR stop ' :: SpecifyBlockInfo not implemented for current matrix type'
   end subroutine SpecifyBlockInfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)     :: this
      !---------------------------------------------
      ERROR stop ' :: destruct not implemented for current matrix type'
   end subroutine destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------
!  Only for block-diaginal matrices: Factorizes the blocks of the matrix as A=LU
!  -----------------------------------------------------------------------------
   subroutine FactorizeBlocks_LU(this,Factorized)
      implicit none
      !-------------------------------------------------------------
      class(Matrix_t), intent(in)    :: this            !<  This matrix
      class(Matrix_t), intent(inout) :: Factorized      !<  Facorized matrix
      !-------------------------------------------------------------
      integer :: k      ! Counter
      !-------------------------------------------------------------
      ERROR stop ' :: FactorizeBlocks_LU not implemented for current matrix type'
   end subroutine FactorizeBlocks_LU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------------------
!  Only for block-diaginal matrices: Solves linear systems with the factorized blocks
!  ----------------------------------------------------------------------------------
   subroutine SolveBlocks_LU(this,x,b)
      implicit none
      !-------------------------------------------------------------
      class(Matrix_t), intent(in)    :: this                 !<  FACTORIZED matrix for solving the problem
      real(kind=RP)  , intent(in)    :: b(this % NumRows)    !<  RHS
      real(kind=RP)  , intent(inout) :: x(this % NumRows)    !<  Solution
      !-------------------------------------------------------------
      ERROR stop ' :: SolveBlocks_LU not implemented for current matrix type'
   end subroutine SolveBlocks_LU
end module GenericMatrixClass
