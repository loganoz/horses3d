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
!  -> Mii is DenseBlockDiagMatrix_t
!  -> Mbb, Mib and Mbi can be either csrMat_t or PETSCMatrix_t
!
#include "Includes.h"
module StaticCondensedMatrixClass
   use SMConstants
   use GenericMatrixClass
   use CSRMatrixClass               , only: csrMat_t
   use PETScMatrixClass             , only: PETSCMatrix_t
   use DenseBlockDiagonalMatrixClass, only: DenseBlockDiagMatrix_t
   use HexMeshClass                 , only: HexMesh
   use ElementClass                 , only: Element
   use MeshTypes                    , only: EFRONT, EBACK, EBOTTOM, ERIGHT, ETOP, ELEFT
   
   implicit none
   
   private
   public StaticCondensedMatrix_t, INNER_DOF, BOUNDARY_DOF, SC_MATRIX_CSR, SC_MATRIX_PETSC
   
   type ElemInfo_t
      integer, allocatable :: dof_association(:)   ! Whether it is an INNER_DOF or a BOUNDARY_DOF
      integer, allocatable :: perm_Indexes(:)      ! Permutation indexes for the map: element DOF -> Condensed system DOF (Mii- or Mbb-based according to dof_association )
      integer, allocatable :: perm_Indexes_i(:)    ! Permutation indexes for the map: element DOF -> relative index of the corresponding block in Mii
   end type ElemInfo_t
   
   type, extends(Matrix_t) :: StaticCondensedMatrix_t
      class(Matrix_t), allocatable  :: Mbb         ! Boundary to boundary matrix
      class(Matrix_t), allocatable  :: Mib         ! Interior to boundary matrix
      class(Matrix_t), allocatable  :: Mbi         ! Boundary to interior matrix
      type(DenseBlockDiagMatrix_t)  :: Mii         ! Interior to interior matirx
      
      type(ElemInfo_t), allocatable :: ElemInfo(:)
      integer, allocatable          :: inner_blockSizes(:)
      
      integer                       :: MatrixType
      integer                       :: size_b      ! Size of condensed system (size of Mbb)
      integer                       :: size_i      ! Size of inner system (size of Mii)
      integer                       :: maxnumCon   ! Maximum number of connections of any element
      
      logical                       :: ignore_boundaries = .FALSE. ! When .TRUE., Mii, does not contain the DOFs on the physical boundaries
      
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
         procedure :: getSchurComplement           => Static_getSchurComplement
         
         procedure :: constructPermutationArrays   => Static_constructPermutationArrays
         procedure :: getCSR                       => Static_getCSR
   end type StaticCondensedMatrix_t
   
!
!  **********
!  Parameters
!  **********
!
!  DOF types
!  ---------
   integer, parameter :: INNER_DOF = 1
   integer, parameter :: BOUNDARY_DOF = 2
!
!  Matrix types
!  ------------
   integer, parameter :: SC_MATRIX_CSR   = 1
   integer, parameter :: SC_MATRIX_PETSC = 2
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,num_of_TotalRows,withMPI)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t) :: this     !<> This matrix
      integer, optional, intent(in)  :: num_of_Rows
      integer, optional, intent(in)  :: num_of_Cols
      integer, optional, intent(in)  :: num_of_Blocks
      integer, optional, intent(in)  :: num_of_TotalRows
      logical, optional, intent(in)  :: WithMPI
      !---------------------------------------------
      
      if ( .not. present(num_of_Rows) ) then
         error stop 'StaticCondensedMatrix_t needs num_of_Rows'
      end if
      if ( .not. present(num_of_Blocks) ) then
         error stop 'StaticCondensedMatrix_t needs num_of_Blocks'
      end if
      if ( present(num_of_Cols) ) then
         if (num_of_Cols /= num_of_Rows) then
            error stop 'StaticCondensedMatrix_t must be a square matrix'
         end if
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
      deallocate (this % Mbb)
      call this % Mib % destruct
      deallocate (this % Mib)
      call this % Mbi % destruct
      deallocate (this % Mbi)
      
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
   subroutine Static_Preallocate(this, nnz, nnzs)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this    !<> This matrix
      integer, optional             , intent(in)    :: nnz     !<  Not needed here
      integer, optional             , intent(in)    :: nnzs(:) !<  nnzs contains the block sizes!
      !---------------------------------------------
      
      select case (this % MatrixType)
         case (SC_MATRIX_CSR)
!
!           CSR: Preallocate with standard linked-list
!           ------------------------------------------
            call this % Mbb % Preallocate()
            call this % Mib % Preallocate()
            call this % Mbi % Preallocate()
         case (SC_MATRIX_PETSC)
!
!           PETSc: Preallocate with a (very conservative) estimate of the maximum nnz
!           -------------------------------------------------------------------------
            call this % Mbb % Preallocate( nnz = this % maxnumCon * maxval(this % BlockSizes - this % inner_blockSizes) )
            call this % Mib % Preallocate( nnz = this % maxnumCon * maxval(this % inner_blockSizes) )
            call this % Mbi % Preallocate( nnz = this % maxnumCon * maxval(this % BlockSizes - this % inner_blockSizes) )
      end select
      call this % Mii % Preallocate(nnzs = this % inner_blockSizes)
      
   end subroutine Static_Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Static_Reset(this, ForceDiagonal)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this
      logical, optional             , intent(in)    :: ForceDiagonal
      !-local-variables-----------------------------
      logical :: mustForceDiagonal
      !---------------------------------------------
      
      if ( present(ForceDiagonal) ) then
         mustForceDiagonal = ForceDiagonal
      else
         mustForceDiagonal = .FALSE.
      end if
      
      call this % Mii % reset(ForceDiagonal = mustForceDiagonal)
      call this % Mbb % reset(ForceDiagonal = mustForceDiagonal)
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
      integer :: dof_association(2) ! local DOF associations
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
                  error stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
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
                  error stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
            end select
         
         case default
            error stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
            
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
      integer :: dof_association(2) ! local DOF associations
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
                  error stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
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
                  error stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
                  
            end select
         
         case default
            error stop 'StaticCondensedMatrix_t :: wrong permutation indexes'
            
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
   subroutine Static_constructPermutationArrays(this,mesh,nEqn,MatrixType,ignore_boundaries)
      implicit none
      !-arguments-----------------------------------
      class(StaticCondensedMatrix_t), intent(inout) :: this    !<> This matrix
      type(HexMesh), target         , intent(in)    :: mesh
      integer                       , intent(in)    :: nEqn
      integer                       , intent(in)    :: MatrixType
      logical      , optional       , intent(in)    :: ignore_boundaries
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
      
      if ( present(ignore_boundaries) ) then
         this % ignore_boundaries = ignore_boundaries
      end if
      
      this % MatrixType = MatrixType
      
      count_i = 0
      count_b = 0
      
      nelem = mesh % no_of_elements
      this % num_of_Blocks = nelem
      Nx => mesh % Nx
      Ny => mesh % Ny
      Nz => mesh % Nz
      
      safedeallocate(this % ElemInfo)
      safedeallocate(this % inner_blockSizes)
      safedeallocate(this % BlockSizes)
      
      allocate ( this % ElemInfo(nelem) )
      allocate ( this % inner_blockSizes(nelem) )
      allocate ( this % BlockSizes(nelem) )
      
      this % maxnumCon = 0
      
      if (this % ignore_boundaries) then
!
!        Do not assign (physical) boundary DOFs to Mii
!        ---------------------------------------------
         do eID = 1, nelem
            e => mesh % elements(eID)
            NDOF = nEqn * (Nx(eID) + 1) * (Ny(eID) + 1) * (Nz(eID) + 1)
            this % inner_blockSizes(eID) = nEqn * (Nx(eID) - 1) * (Ny(eID) - 1) * (Nz(eID) - 1)
            this %       BlockSizes(eID) = nEqn * (Nx(eID) + 1) * (Ny(eID) + 1) * (Nz(eID) + 1)
            
            allocate ( this % ElemInfo(eID) % perm_Indexes   (NDOF) )
            allocate ( this % ElemInfo(eID) % perm_Indexes_i (NDOF) )
            allocate ( this % ElemInfo(eID) % dof_association(NDOF) )
            this % ElemInfo(eID) % perm_Indexes_i = -1
            
            this % maxnumCon = max(this % maxnumCon,sum(e % NumberOfConnections))
            
            count_ii = 0
            count_el = 0
            do k=0, Nz(eID) ; do j=0, Ny(eID) ; do i=0, Nx(eID) ; do eq=1, nEqn
               count_el = count_el + 1
               
               if ( (i==0) .or. (i==Nx(eID)) .or. (j==0) .or. (j==Ny(eID)) .or. (k==0) .or. (k==Nz(eID)) ) then
                  this % ElemInfo(eID) % dof_association(count_el) = BOUNDARY_DOF
                  count_b = count_b + 1
                  this % ElemInfo(eID) % perm_Indexes(count_el) = count_b
               else
                  this % ElemInfo(eID) % dof_association(count_el) = INNER_DOF
                  count_i = count_i + 1
                  this % ElemInfo(eID) % perm_Indexes(count_el) = count_i
                  count_ii = count_ii + 1
                  this % ElemInfo(eID) % perm_Indexes_i(count_el) = count_ii
               end if
               
            end do          ; end do          ; end do          ; end do
            
         end do
      else
!
!        Assign (physical) boundary DOFs to Mii (DEFAULT)
!        ------------------------------------------------
         do eID = 1, nelem
            e => mesh % elements(eID)
            NDOF = nEqn * (Nx(eID) + 1) * (Ny(eID) + 1) * (Nz(eID) + 1)
            this % inner_blockSizes(eID) = 0
            this %       BlockSizes(eID) = nEqn * (Nx(eID) + 1) * (Ny(eID) + 1) * (Nz(eID) + 1)
            
            allocate ( this % ElemInfo(eID) % perm_Indexes   (NDOF) )
            allocate ( this % ElemInfo(eID) % perm_Indexes_i (NDOF) )
            allocate ( this % ElemInfo(eID) % dof_association(NDOF) )
            this % ElemInfo(eID) % perm_Indexes_i = -1
            
            this % maxnumCon = max(this % maxnumCon,sum(e % NumberOfConnections))
            
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
      end if
      
      if ( (count_i+count_b) /= this % num_of_Rows ) then
         error stop 'StaticCondensedMatrixClass :: Invalid arguments in constructPermutationArrays'
      end if
      
      this % size_i        = count_i
      this % size_b        = count_b
!
!     Assign matrix types
!     -------------------
!
      select case (this % MatrixType)
         case (SC_MATRIX_CSR)
            allocate(csrMat_t :: this % Mbb)
            allocate(csrMat_t :: this % Mbi)
            allocate(csrMat_t :: this % Mib)
         case (SC_MATRIX_PETSC)
            allocate(PETSCMatrix_t :: this % Mbb)
            allocate(PETSCMatrix_t :: this % Mbi)
            allocate(PETSCMatrix_t :: this % Mib)
      end select
!
!     Construct matrices
!     ------------------
!
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
      
      select type(Mbb => this % Mbb) ; class is(csrMat_t) ; select type(Mbi => this % Mbi) ; class is(csrMat_t) ; select type(Mib => this % Mib) ; class is(csrMat_t)
      
      num_of_entries = size(Mbb % Values) + size(Mib % Values) + size(Mbi % Values) + sum(this % Mii % BlockSizes**2)
      
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
         do kk=Mbb % Rows(i), Mbb % Rows(i+1)-1
            Cols(k)   = Mbb % Cols  (kk)
            Values(k) = Mbb % Values(kk)
            k = k + 1
         end do
         
!        Mib contribution
!        ----------------
         do kk=Mib % Rows(i), Mib % Rows(i+1)-1
            Cols(k)   = Mib % Cols  (kk) + this % size_b
            Values(k) = Mib % Values(kk)
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
         do kk=Mbi % Rows(ii), Mbi % Rows(ii+1)-1
            Cols(k)   = Mbi % Cols  (kk)
            Values(k) = Mbi % Values(kk)
            k = k + 1
         end do
         
!        Mii contribution
!        ----------------
         do bID=lastbID, this % Mii % num_of_Blocks   ! Search the block where this row is contained
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
      call Acsr % constructWithCSRArrays(Rows, Cols, Values, this % num_of_Cols)
      
      class default
         error stop 'ERROR :: Static_getCSR not defined for submatrices /= cstMat_t'
      end select ; end select ; end select
      
   end subroutine Static_getCSR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Static_getSchurComplement:
!  Compute the Schur complement with the inverse of Mii in correct format
!  ----------------------------------------------------------------------
   subroutine Static_getSchurComplement(this, Mii_inv, SchurComp)
      implicit none
      !-arguments------------------------------------------------------
      class(StaticCondensedMatrix_t), intent(in)  :: this
      class(Matrix_t)               , intent(in)  :: Mii_inv
      class(Matrix_t)               , intent(inout) :: SchurComp
      !-local-variables------------------------------------------------
      class(Matrix_t), allocatable :: Mii_inv_Mbi 
      class(Matrix_t), allocatable :: Mib_Mii_inv_Mbi
      integer :: BlockIdx(this % num_of_Blocks+1), BlockSizes(this % num_of_Blocks), bID
      !----------------------------------------------------------------
      
      
!     Construct auxiliary matrices
!     ---------------------------
      select case (this % MatrixType)
         case (SC_MATRIX_CSR)
            allocate (csrMat_t :: Mii_inv_Mbi)
            allocate (csrMat_t :: Mib_Mii_inv_Mbi)
         case (SC_MATRIX_PETSC)
            allocate (PETSCMatrix_t :: Mii_inv_Mbi)
            allocate (PETSCMatrix_t :: Mib_Mii_inv_Mbi)
      end select
      call Mii_inv_Mbi % construct (num_of_Rows = this % size_i, num_of_Cols = this % size_b)
      call Mib_Mii_inv_Mbi % construct (num_of_Rows = this % size_b, num_of_Cols = this % size_b)
      
!     Perform operations
!     ------------------
      call Mii_inv % MatMatMul( this % Mbi, Mii_inv_Mbi)
      call this % Mib % MatMatMul( Mii_inv_Mbi, Mib_Mii_inv_Mbi )
      call Mii_inv_Mbi % destruct ; deallocate (Mii_inv_Mbi)
      
      call this % Mbb % MatAdd (Mib_Mii_inv_Mbi, SchurComp, -1._RP)
      call Mib_Mii_inv_Mbi % destruct ; deallocate(Mib_Mii_inv_Mbi)
      
!     Specify info of matrix
!     ----------------------
      BlockSizes = this % BlockSizes - this % inner_blockSizes
      BlockIdx(1) = 1
      do bID=1, this % num_of_Blocks
         BlockIdx(bID+1) = BlockIdx(bID) + BlockSizes(bID)
      end do
      call SchurComp % SpecifyBlockInfo(BlockIdx, BlockSizes)
      
   end subroutine Static_getSchurComplement
end module StaticCondensedMatrixClass