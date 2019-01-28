!
!//////////////////////////////////////////////////////
!
!   @File:    StaticCondensedMatrixClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Dec  4 16:26:02 2018
!   @Last revision date: Mon Jan 28 18:24:36 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: dee8ab8d5e78ed3d128a4350162a11d8455b300d
!
!//////////////////////////////////////////////////////
!
!  StaticCondensedMatrixClass:
!  -> This type of matrix reorganizes the Jacobian matrix as 
!
!      +-----+-------+
!      | Mbb | Mib   |
!      |     |       |
!      +-----+-------+
!      |     |       |
!      | Mbi | Mii   |
!      |     |       |
!      +-----+-------+
!
#include "Includes.h"
module StaticCondensedMatrixClass
   use SMConstants
   use GenericMatrixClass
   use CSRMatrixClass               , only: csrMat_t
   use DenseBlockDiagonalMatrixClass, only: DenseBlockDiagMatrix_t
   use HexMeshClass                 , only: HexMesh
   use ElementClass                 , only: Element
   use MeshTypes                    , only: EFRONT, EBACK, EBOTTOM, ERIGHT, ETOP, ELEFT
   
   implicit none
   
   private
   public StaticCondensedMatrix_t, INNER_DOF, BOUNDARY_DOF
   
   type ElemInfo_t
      integer, allocatable :: dof_association(:)   ! Whether it is an INNER_DOF or a BOUNDARY_DOF
      integer, allocatable :: perm_Indexes(:)      ! Permutation indexes for the map: element DOF -> Condensed system DOF (Mii- or Mbb-based according to dof_association )
      integer, allocatable :: perm_Indexes_i(:)    ! Permutation indexes for the map: element DOF -> relative index of the corresponding block in Mii
   end type ElemInfo_t
   
   type, extends(Matrix_t) :: StaticCondensedMatrix_t
      type(csrMat_t)                :: Mbb         ! Boundary to boundary matrix
      type(csrMat_t)                :: Mib         ! Interior to boundary matrix
      type(csrMat_t)                :: Mbi         ! Boundary to interior matrix
      type(DenseBlockDiagMatrix_t)  :: Mii         ! Interior to interior matirx
      
      type(ElemInfo_t), allocatable :: ElemInfo(:)
      integer, allocatable          :: inner_blockSizes(:)
      integer, allocatable          :: BlockSizes(:)
      
      integer                       :: size_b      ! Size of condensed system (size of Mbb)
      integer                       :: size_i      ! Size of inner system (size of Mii)
      integer                       :: num_of_Blocks
      integer                       :: num_of_Cols
      
      contains
         procedure :: construct                    => Static_construct
         procedure :: destruct                     => Static_destruct
         procedure :: Preallocate                  => Static_Preallocate
         procedure :: reset                        => Static_reset
         procedure :: SpecifyBlockInfo             => Static_SpecifyBlockInfo
         procedure :: SetBlockEntry                => Static_SetBlockEntry
         procedure :: AddToBlockEntry              => Static_AddToBlockEntry
         procedure :: Assembly                     => Static_Assembly
         procedure :: shift                        => Static_Shift
         
         procedure :: constructPermutationArrays   => Static_constructPermutationArrays
         procedure :: getCSR                       => Static_getCSR
   end type StaticCondensedMatrix_t
   
   integer, parameter :: INNER_DOF = 1
   integer, parameter :: BOUNDARY_DOF = 2
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,withMPI)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t) :: this     !<> This matrix
      integer, optional, intent(in)  :: num_of_Rows
      integer, optional, intent(in)  :: num_of_Cols
      integer, optional, intent(in)  :: num_of_Blocks
      logical, optional, intent(in)  :: WithMPI
      !---------------------------------------------
      
      if ( .not. present(num_of_Rows) ) then
         ERROR stop 'StaticCondensedMatrix_t needs num_of_Rows'
      end if
      if ( .not. present(num_of_Blocks) ) then
         ERROR stop 'StaticCondensedMatrix_t needs num_of_Blocks'
      end if
      
      this % num_of_Rows   = num_of_Rows
      this % num_of_Cols   = num_of_Rows     ! Only for square matrices
      this % num_of_Blocks = num_of_Blocks
      
   end subroutine Static_construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_destruct(this)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this    !<> This matrix
      !---------------------------------------------
      
      call this % Mii % destruct
      call this % Mbb % destruct
      call this % Mib % destruct
      call this % Mbi % destruct
      
      deallocate (this % ElemInfo)
      deallocate (this % inner_blockSizes)
      deallocate (this % BlockSizes)
      
      this % size_b        = 0
      this % size_i        = 0
      this % num_of_Rows   = 0
      this % num_of_Cols   = 0
      this % num_of_Blocks = 0
   end subroutine Static_destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_Preallocate(this, nnz, nnzs, ForceDiagonal)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this    !<> This matrix
      integer, optional             , intent(in)    :: nnz     !<  Not needed here
      integer, optional             , intent(in)    :: nnzs(:) !<  nnzs contains the block sizes!
      logical, optional             , intent(in)    :: ForceDiagonal
      !-local-variables-----------------------------
      logical :: mustForceDiagonal
      !---------------------------------------------
      
      if ( present(ForceDiagonal) ) then
         mustForceDiagonal = ForceDiagonal
      else
         mustForceDiagonal = .FALSE.
      end if
      
      call this % Mbb % Preallocate(ForceDiagonal = mustForceDiagonal)
      call this % Mib % Preallocate()
      call this % Mbi % Preallocate()
      call this % Mii % Preallocate(nnzs = this % inner_blockSizes, ForceDiagonal = mustForceDiagonal)
      
   end subroutine Static_Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_Reset(this)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      !---------------------------------------------
      
      call this % Mii % reset
      call this % Mbb % reset
      call this % Mbi % reset
      call this % Mib % reset
      
   end subroutine Static_Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_SpecifyBlockInfo(this,BlockIdx,BlockSize)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      integer                       , intent(in)    :: BlockIdx(:)
      integer                       , intent(in)    :: BlockSize(:)
      !---------------------------------------------
      ! do nothing
   end subroutine Static_SpecifyBlockInfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------
!  Subroutine to set the entries of a block with relative index
!  ------------------------------------------------------------
   subroutine Static_SetBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      integer                       , intent(in)    :: iBlock, jBlock
      integer                       , intent(in)    :: i, j
      real(kind=RP)                 , intent(in)    :: value
      !-local-variables-----------------------------
      integer :: dof_association(2) ! local DOF associtions
      integer :: perm_Indexes(2)    ! local permutation indexes
      !---------------------------------------------
      
      dof_association = [ this % ElemInfo(iBlock) % dof_association(i) , this % ElemInfo(jBlock) % dof_association(j) ]
      
      select case ( dof_association(1) )
         case (INNER_DOF)
            select case ( dof_association(2) )
               case (INNER_DOF)
!
!                 Contribution to Mii matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes_i (i) , this % ElemInfo(jBlock) % perm_Indexes_i (j) ]
                  call this % Mii % SetBlockEntry (iBlock, jBlock, perm_Indexes(1), perm_Indexes(2), value)
                  
               case (BOUNDARY_DOF)
!
!                 Contribution to Mbi matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes   (i) , this % ElemInfo(jBlock) % perm_Indexes   (j) ]
                  call this % Mbi % SetEntry(perm_Indexes(1), perm_Indexes(2), value )
                  
               case default
                  ERROR stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
            end select
         case (BOUNDARY_DOF)
            select case ( dof_association(2) )
               case (INNER_DOF)
!
!                 Contribution to Mib matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes   (i) , this % ElemInfo(jBlock) % perm_Indexes   (j) ]
                  call this % Mib % SetEntry(perm_Indexes(1), perm_Indexes(2), value )
                  
               case (BOUNDARY_DOF)
!
!                 Contribution to Mbb matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes   (i) , this % ElemInfo(jBlock) % perm_Indexes   (j) ]
                  call this % Mbb % SetEntry(perm_Indexes(1), perm_Indexes(2), value )
                  
               case default
                  ERROR stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
            end select
         
         case default
            ERROR stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
            
      end select
   end subroutine Static_SetBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine to add a value to the entries of a block with relative index
!  -----------------------------------------------------------------------
   subroutine Static_AddToBlockEntry(this, iBlock, jBlock, i, j, value )
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      integer                       , intent(in)    :: iBlock, jBlock
      integer                       , intent(in)    :: i, j
      real(kind=RP)                 , intent(in)    :: value
      !-local-variables-----------------------------
      integer :: dof_association(2) ! local DOF associtions
      integer :: perm_Indexes(2)    ! local permutation indexes
      !---------------------------------------------
      
      dof_association = [ this % ElemInfo(iBlock) % dof_association(i) , this % ElemInfo(jBlock) % dof_association(j) ]
      
      select case ( dof_association(1) )
         case (INNER_DOF)
            select case ( dof_association(2) )
               case (INNER_DOF)
!
!                 Contribution to Mii matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes_i (i) , this % ElemInfo(jBlock) % perm_Indexes_i (j) ]
                  call this % Mii % AddToBlockEntry (iBlock, jBlock, perm_Indexes(1), perm_Indexes(2), value)
                  
               case (BOUNDARY_DOF)
!
!                 Contribution to Mbi matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes   (i) , this % ElemInfo(jBlock) % perm_Indexes   (j) ]
                  call this % Mbi % AddToEntry(perm_Indexes(1), perm_Indexes(2), value )
                  
               case default
                  ERROR stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
            end select
         case (BOUNDARY_DOF)
            select case ( dof_association(2) )
               case (INNER_DOF)
!
!                 Contribution to Mib matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes   (i) , this % ElemInfo(jBlock) % perm_Indexes   (j) ]
                  call this % Mib % AddToEntry(perm_Indexes(1), perm_Indexes(2), value )
                  
               case (BOUNDARY_DOF)
!
!                 Contribution to Mbb matrix
!                 --------------------------
                  perm_Indexes = [ this % ElemInfo(iBlock) % perm_Indexes   (i) , this % ElemInfo(jBlock) % perm_Indexes   (j) ]
                  call this % Mbb % AddToEntry(perm_Indexes(1), perm_Indexes(2), value )
                  
               case default
                  ERROR stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
            end select
         
         case default
            ERROR stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
            
      end select
      
   end subroutine Static_AddToBlockEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_Assembly(this)
      implicit none
      !---------------------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      !---------------------------------------------
      
      call this % Mii % Assembly
      call this % Mbb % Assembly
      call this % Mib % Assembly
      call this % Mbi % Assembly
      
   end subroutine Static_Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_Shift(this, shiftval)
      implicit none
      !---------------------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      real(kind=RP)                 , intent(in)    :: shiftval
      !---------------------------------------------
      
      call this % Mbb % shift(shiftval)
      call this % Mii % shift(shiftval)
      
   end subroutine Static_Shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  constructPermutationArrays:
!     Routine to construct the "ElemInfo"
!  --------------------------------------
   subroutine Static_constructPermutationArrays(this,mesh,nEqn)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this    !<> This matrix
      type(HexMesh), target         , intent(in)    :: mesh
      integer                       , intent(in)    :: nEqn
      !-local-variables-----------------------------
      integer :: eID
      integer :: i,j,k,eq
      integer :: count_el  ! counter for the element's DOFs
      integer :: count_i   ! counter for the inner    (elemental) DOFs
      integer :: count_b   ! counter for the boundary (elemental) DOFs
      integer :: count_ii  ! counter for the inner    (elemental) DOFs, relative to the corresponding block of the Mii matrix
      integer :: nelem
      integer :: NDOF
      integer        , pointer :: Nx(:)
      integer        , pointer :: Ny(:)
      integer        , pointer :: Nz(:)
      type(Element)  , pointer :: e
      !---------------------------------------------
      
      count_i = 0
      count_b = 0
      
      nelem = mesh % no_of_elements
      Nx => mesh % Nx
      Ny => mesh % Ny
      Nz => mesh % Nz
      
      safedeallocate(this % ElemInfo)
      safedeallocate(this % inner_blockSizes)
      safedeallocate(this % BlockSizes)
      
      allocate ( this % ElemInfo(nelem) )
      allocate ( this % inner_blockSizes(nelem) )
      allocate ( this % BlockSizes(nelem) )
      
      do eID = 1, nelem
         e => mesh % elements(eID)
         NDOF = nEqn * (Nx(eID) + 1) * (Ny(eID) + 1) * (Nz(eID) + 1)
         this % inner_blockSizes(eID) = 0 !nEqn * (Nx(eID) - 1) * (Ny(eID) - 1) * (Nz(eID) - 1)
         this %       BlockSizes(eID) = nEqn * (Nx(eID) + 1) * (Ny(eID) + 1) * (Nz(eID) + 1)
         
         allocate ( this % ElemInfo(eID) % perm_Indexes   (NDOF) )
         allocate ( this % ElemInfo(eID) % perm_Indexes_i (NDOF) )
         allocate ( this % ElemInfo(eID) % dof_association(NDOF) )
         this % ElemInfo(eID) % perm_Indexes_i = -1
         
         count_ii = 0
         count_el = 0
         do k=0, Nz(eID) ; do j=0, Ny(eID) ; do i=0, Nx(eID) ; do eq=1, nEqn
            count_el = count_el + 1
            
            if (     (j==0       .and. e % NumberOfConnections(EFRONT ) > 0 ) &
                .or. (j==Ny(eID) .and. e % NumberOfConnections(EBACK  ) > 0 ) &
                .or. (k==0       .and. e % NumberOfConnections(EBOTTOM) > 0 ) &
                .or. (i==Nx(eID) .and. e % NumberOfConnections(ERIGHT ) > 0 ) &
                .or. (k==Nz(eID) .and. e % NumberOfConnections(ETOP   ) > 0 ) &
                .or. (i==0       .and. e % NumberOfConnections(ELEFT  ) > 0 )    ) then
               this % ElemInfo(eID) % dof_association(count_el) = BOUNDARY_DOF
               count_b = count_b + 1
               this % ElemInfo(eID) % perm_Indexes(count_el) = count_b
            else
               this % ElemInfo(eID) % dof_association(count_el) = INNER_DOF
               count_i = count_i + 1
               this % ElemInfo(eID) % perm_Indexes(count_el) = count_i
               count_ii = count_ii + 1
               this % ElemInfo(eID) % perm_Indexes_i(count_el) = count_ii
               this % inner_blockSizes(eID) = this % inner_blockSizes(eID) + 1
            end if
            
         end do          ; end do          ; end do          ; end do
         
      end do
      
      if ( (count_i+count_b) /= this % num_of_Rows ) then
         ERROR stop 'StaticCondensedMatrixClass :: Ivalid arguments in constructPermutationArrays'
      end if
      
      this % size_i        = count_i
      this % size_b        = count_b
      
      call this % Mii % construct (num_of_Blocks = this % num_of_Blocks)
      call this % Mbb % construct (num_of_Rows   = this % size_b)
      call this % Mib % construct (num_of_Rows   = this % size_b , num_of_Cols = this % size_i)
      call this % Mbi % construct (num_of_Rows   = this % size_i , num_of_Cols = this % size_b)
      
   end subroutine Static_constructPermutationArrays
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Convert static-condensed matrix to CSR matrix
!  ---------------------------------------------
   subroutine Static_getCSR(this,Acsr)
      implicit none
      !-arguments------------------------------------------------------
      class(StaticCondensedMatrix_t), intent(inout)  :: this
      type(csrMat_t)                , intent(inout)  :: Acsr
      !-local-variables------------------------------------------------
      integer :: i, j         ! Indexes of global matrix
      integer :: ii           ! Indexes of submatrices
      integer :: k            ! Current value number of global matrix
      integer :: kk           ! Current value number of submatrix
      integer :: bi, bj       ! Row and column index in a block
      integer :: bID, lastbID
      integer :: num_of_entries
      real(kind=RP), allocatable  :: Values(:)
      integer      , allocatable  :: Cols(:), Rows(:)
      !----------------------------------------------------------------
      
      num_of_entries = size(this % Mbb % Values) + size(this % Mib % Values) + size(this % Mbi % Values) + sum(this % Mii % BlockSizes**2)
      
      allocate ( Rows  (this % num_of_Rows + 1) )
      allocate ( Cols  (num_of_entries) )
      allocate ( Values(num_of_entries) )
      
      k = 1
      
!     Fill the boundary degrees of freedom
!     ------------------------------------
      
      do i=1, this % size_b
         Rows(i) = k
         
!        Mbb contribution
!        ----------------
         do kk=this % Mbb % Rows(i), this % Mbb % Rows(i+1)-1
            Cols(k)   = this % Mbb % Cols  (kk)
            Values(k) = this % Mbb % Values(kk)
            k = k + 1
         end do
         
!        Mib contribution
!        ----------------
         do kk=this % Mib % Rows(i), this % Mib % Rows(i+1)-1
            Cols(k)   = this % Mib % Cols  (kk) + this % size_b
            Values(k) = this % Mib % Values(kk)
            k = k + 1
         end do
      end do
      
!     Fill the inner degrees of freedom
!     ---------------------------------
      lastbID = 1
      
      do ii=1, this % size_i
         i = ii + this % size_b
         Rows(i) = k
         
!        Mbi contribution
!        ----------------
         do kk=this % Mbi % Rows(ii), this % Mbi % Rows(ii+1)-1
            Cols(k)   = this % Mbi % Cols  (kk)
            Values(k) = this % Mbi % Values(kk)
            k = k + 1
         end do
         
!        Mii contribution
!        ----------------
         do bID=lastbID, this % Mii % NumOfBlocks   ! Search the block where this row is contained
            if (this % Mii % BlockIdx(bID+1) > ii) then
               
               do bi=1, this % Mii % BlockSizes(bID)  ! Search the block row that corresponds to this row
                  
                  if ( this % Mii % BlockIdx(bID) + bi -1 == ii) then
                     
                     do bj=1, this % Mii % BlockSizes(bID)
                        Cols(k)   = (this % Mii % BlockIdx(bID) + bj - 1) + this % size_b
                        Values(k) = this % Mii % Blocks(bID) % Matrix(bi,bj)
                        k = k + 1
                     end do
                     
                     exit
                  end if
               end do
               
               lastbID = bID
               exit
            end if
         end do
      end do
      Rows (this % num_of_Rows + 1) = k
      
!     Finish constructing the matrix
!     ------------------------------
      call Acsr % constructWithArrays(Rows, Cols, Values, this % num_of_Cols)
      
   end subroutine Static_getCSR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module StaticCondensedMatrixClass
