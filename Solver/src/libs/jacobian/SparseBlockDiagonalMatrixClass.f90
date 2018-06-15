!////////////////////////////////////////////////////////////////////////
!
!      SparseBlockDiagonalMatrixClass.f90
!      Created: 2018-02-19 17:07:00 +0100 
!      By: Andrés Rueda
!
!      Class for sparse block diagonal matrices with sparse CSR blocks
!
!     -> If block factorization and system solving is needed, MKL pardiso is needed
!
!////////////////////////////////////////////////////////////////////////
module SparseBlockDiagonalMatrixClass
   use SMConstants
   use GenericMatrixClass
   use CSRMatrixClass, only: csrMat_t
#include "Includes.h"
   implicit none
   
   private
   public SparseBlockDiagMatrix_t, Matrix_t, Block_t
   
   type Block_t
      type(csrMat_t)             :: Matrix
      
      ! Pardiso variables
      integer                        :: mtype                              ! Matrix type. See construct
      integer              , pointer :: Pardiso_iparm(:) => NULL()         ! Parameters for mkl version of pardiso
      integer(kind=AddrInt), pointer :: Pardiso_pt(:)    => NULL()  
      integer      , allocatable     :: Indexes(:)
   end type Block_t
   
   type, extends(Matrix_t) :: SparseBlockDiagMatrix_t
      type(Block_t), allocatable :: Blocks(:)   ! Array containing each block in a dense matrix
      integer                    :: NumOfBlocks ! Number of blocks in matrix
      integer      , allocatable :: BlockSizes(:)
      integer      , allocatable :: BlockIdx(:)
      
      contains
         procedure :: construct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: ResetBlock
         procedure :: SetBlockEntry
         procedure :: AddToBlockEntry
         procedure :: shift
         procedure :: Assembly
         procedure :: destruct
         procedure :: FactorizeBlocks_LU
         procedure :: SolveBlocks_LU
   end type SparseBlockDiagMatrix_t
   
!
!  Interfaces for pardiso routines
!  *******************************
   
   interface
      subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, values, rows, cols, perm, nrhs, iparm, msglvl, b, x, ierror)
         use SMConstants
         implicit none
         real(kind=RP) :: values(*), b(*), x(*)
         integer(kind=AddrInt) :: pt(*)
         integer :: perm(*), nrhs, iparm(*), msglvl, ierror
         integer :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
      end subroutine pardiso
      subroutine pardisoinit(pt, mtype, iparm)
         USE SMConstants
         implicit none
         integer(kind=AddrInt) :: pt(*)
         integer :: mtype
         integer :: iparm(*)
      end subroutine pardisoinit
   end interface
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,dimPrb,withMPI)
      implicit none
      !-arguments-----------------------------------
      class(SparseBlockDiagMatrix_t) :: this     !<> This matrix
      integer          , intent(in) :: dimPrb   !<  Number of blocks of the matrix!
      logical, optional, intent(in) :: WithMPI
      !-local-variables-----------------------------
      integer :: iBL
      !---------------------------------------------
      
      allocate ( this % Blocks(dimPrb) )
      this % NumOfBlocks = dimPrb
      allocate ( this % BlockSizes(dimPrb) )
      allocate ( this % BlockIdx(dimPrb+1) )
      
#ifdef HAS_MKL
      do iBL = 1, dimPrb
         this % Blocks(iBL) % mtype = 11 !Set matrix type to real unsymmetric (change?)
         allocate(this % Blocks(iBL) % Pardiso_pt   (64))
         allocate(this % Blocks(iBL) % Pardiso_iparm(64))
      
         call pardisoinit(this % Blocks(iBL) % Pardiso_pt, this % Blocks(iBL) % mtype, this % Blocks(iBL) % Pardiso_iparm)
      end do
#endif
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs)
      IMPLICIT NONE
      !---------------------------------------------
      class(SparseBlockDiagMatrix_t), intent(inout) :: this    !<> This matrix
      integer, optional            , intent(in)    :: nnz     !<  Not needed here
      integer, optional            , intent(in)    :: nnzs(:) !<  nnzs contains the block sizes!
      !---------------------------------------------
      integer :: i, k ! counters
      !---------------------------------------------
      
      if (.not. present(nnzs) ) ERROR stop ':: SparseBlockDiagMatrix needs the block sizes'
      if ( size(nnzs) /= this % NumOfBlocks) ERROR stop ':: SparseBlockDiagMatrix: wrong dimension for the block sizes'
      
      this % BlockSizes = nnzs
      this % NumRows = sum(nnzs)
      
      this % BlockIdx(1) = 1
      do i=2, this % NumOfBlocks + 1
         this % BlockIdx(i) = this % BlockIdx(i-1) + nnzs(i-1)
      end do
      
!$omp parallel do private(k) schedule(runtime)
      do i=1, this % NumOfBlocks
         
         call this % Blocks(i) % Matrix % construct ( nnzs(i) )
         call this % Blocks(i) % Matrix % PreAllocate()
         
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
      class(SparseBlockDiagMatrix_t), intent(inout) :: this     !<> This matrix
      !---------------------------------------------
      integer :: i
      !---------------------------------------------
      
      do i=1, this % NumOfBlocks
         call this % Blocks(i) % Matrix % Reset
      end do
      
   end subroutine Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  --------------------------------------------------------------
   subroutine SetBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(SparseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: iBlock, jBlock
      integer                      , intent(in)    :: i, j
      real(kind=RP)                , intent(in)    :: value
      !---------------------------------------------
      
      if (iBlock /= jBlock) return ! Only diagonal blocks here!
      
      call this % Blocks(iBlock) % Matrix % SetEntry (i, j, value)
      
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
      class(SparseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: iBlock, jBlock
      integer                      , intent(in)    :: i, j
      real(kind=RP)                , intent(in)    :: value
      !---------------------------------------------
      
      if (iBlock /= jBlock) return ! Only diagonal blocks here!
      
      call this % Blocks(iBlock) % Matrix % AddToEntry (i, j, value)
      
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
      class(SparseBlockDiagMatrix_t), intent(inout) :: this
      integer                      , intent(in)    :: iBlock, jBlock
      !---------------------------------------------
      
      if (iBlock /= jBlock) return ! Only diagonal blocks here!
      
      call this % Blocks(iBlock) % Matrix % Reset
      
   end subroutine ResetBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine shift(this,shiftval)
      implicit none
      !------------------------------------------
      class(SparseBlockDiagMatrix_t), intent(inout) :: this
      real(kind=RP)                , intent(in)    :: shiftval
      !------------------------------------------
      INTEGER                :: iBL
      !------------------------------------------
      
!$omp parallel do schedule(runtime)
      do iBL=1, this % NumOfBlocks
         call this % Blocks(iBL) % Matrix % shift (shiftval)
      end do
!$omp end parallel do
      
   end subroutine shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Assembly(this,BlockIdx,BlockSize)
      implicit none
      !---------------------------------------------
      class(SparseBlockDiagMatrix_t), intent(inout) :: this
      integer, target, optional     , intent(in)    :: BlockIdx(:)
      integer, target, optional     , intent(in)    :: BlockSize(:)
      !-local-variables-----------------------------
      integer :: iBL
      !---------------------------------------------
      
      !$omp parallel do schedule(runtime)
      do iBL=1, this % NumOfBlocks
         call this % Blocks(iBL) % Matrix % Assembly
      end do
!$omp end parallel do
      
   end subroutine Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !---------------------------------------------
      class(SparseBlockDiagMatrix_t),     intent(inout)     :: this
      !---------------------------------------------
      integer :: i
      !---------------------------------------------
      
      do i = 1, this % NumOfBlocks
         call this % Blocks(i) % Matrix % destruct
         deallocate (this % Blocks(i) % Indexes)
      end do
      
      deallocate (this % Blocks)
      deallocate (this % BlockSizes)
      deallocate (this % BlockIdx)
      
#ifdef HAS_MKL
      do i = 1, this % NumOfBlocks
         deallocate(this % Blocks(i) % Pardiso_pt)
         deallocate(this % Blocks(i) % Pardiso_iparm)
      end do
#endif
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
      class(SparseBlockDiagMatrix_t), intent(in)    :: this            !<  This matrix
      class(SparseBlockDiagMatrix_t), intent(inout) :: Factorized      !<  Facorized matrix
      !-------------------------------------------------------------
      integer                    :: iBL, error      ! Counter
      real(kind=RP), allocatable :: x_loc(:)
      real(kind=RP), allocatable :: b_loc(:)
      !-------------------------------------------------------------
      
#ifdef HAS_MKL
      
!$omp parallel do private(x_loc,b_loc) schedule(runtime)
      do iBL=1, this % NumOfBlocks
         associate (A      => this       % Blocks(iBL) % Matrix, &
                    CBlock => Factorized % Blocks(iBL)   )
         
         call CBlock % Matrix % destruct
         call CBlock % Matrix % constructWithArrays(A % Rows, &
                                                    A % Cols, &
                                                    A % Values)
         
         allocate( x_loc(CBlock % Matrix % NumRows) )
         allocate( b_loc(CBlock % Matrix % NumRows) )
         
         CALL pardiso(  pt      = CBlock % Pardiso_pt      ,     &
                        maxfct  = 1                        ,     &     ! Set up space for 1 matrix at most
                        mnum    = 1                        ,     &     ! Matrix to use in the solution phase (1st and only one)
                        mtype   = CBlock % mtype           ,     &
                        phase   = 12                       ,     &     !  12 for factorization
                        n       = CBlock % Matrix % NumRows,     &     ! Number of equations
                        values  = CBlock % Matrix % Values ,     & 
                        rows    = CBlock % Matrix % Rows   ,     &
                        cols    = CBlock % Matrix % Cols   ,     &
                        perm    = CBlock % Indexes         ,     &     ! ...
                        nrhs    = 1                        ,     &     ! Only one right hand side 
                        iparm   = CBlock % Pardiso_iparm   ,     &
                        msglvl  = 0                        ,     &     ! 1: verbose... Too much printing
                        b       = b_loc                    ,     &
                        x       = x_loc                    ,     &
                        ierror  = error                     )
         
         if (error .NE. 0) then
            write(*,*) 'MKL Pardiso ERROR:', error
            stop
         end if
         
         deallocate (x_loc, b_loc)
         end associate
      end do
!$omp end parallel do
#else
      stop 'SparseBlockDiagMat :: MKL is not linked correctly'
#endif
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
      class(SparseBlockDiagMatrix_t), intent(in)    :: this                 !<  FACTORIZED matrix for solving the problem
      real(kind=RP)                , intent(in)    :: b(this % NumRows)    !<  RHS
      real(kind=RP)                , intent(inout) :: x(this % NumRows)    !<  Solution
      !-------------------------------------------------------------
      integer                    :: iBL, error ! Counter
      real(kind=RP), allocatable :: x_loc(:) ! Local x
      !-------------------------------------------------------------
      
#ifdef HAS_MKL
!$omp parallel do private(x_loc) schedule(runtime)
      do iBL = 1, this % NumOfBlocks
         associate (CBlock => this % Blocks(iBL) )
         allocate( x_loc(this % BlockSizes(iBL)) )
         
         CALL pardiso(  pt      = CBlock % Pardiso_pt      ,     &
                        maxfct  = 1                        ,     &     ! Set up space for 1 matrix at most
                        mnum    = 1                        ,     &     ! Matrix to use in the solution phase (1st and only one)
                        mtype   = CBlock % mtype           ,     &
                        phase   = 33                       ,     &     !  12 for factorization
                        n       = CBlock % Matrix % NumRows,     &     ! Number of equations
                        values  = CBlock % Matrix % Values ,     & 
                        rows    = CBlock % Matrix % Rows   ,     &
                        cols    = CBlock % Matrix % Cols   ,     &
                        perm    = CBlock % Indexes         ,     &     ! ...
                        nrhs    = 1                        ,     &     ! Only one right hand side 
                        iparm   = CBlock % Pardiso_iparm   ,     &
                        msglvl  = 0                        ,     &     ! 1: verbose... Too much printing
                        b       = b (this % BlockIdx(iBL):this % BlockIdx(iBL+1)-1 ),     &
                        x       = x_loc                    ,     &
                        ierror  = error                     )
         
         if (error .NE. 0) then
            write(*,*) 'MKL Pardiso ERROR:', error
            stop
         end if
         
!$omp critical
         x( this % BlockIdx(iBL):this % BlockIdx(iBL+1)-1 ) = x_loc
!$omp end critical
         deallocate(x_loc)
         end associate
      end do
!$omp end parallel do
#else
      stop 'SparseBlockDiagMat :: MKL is not linked correctly'
#endif
      
   end subroutine SolveBlocks_LU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module SparseBlockDiagonalMatrixClass
