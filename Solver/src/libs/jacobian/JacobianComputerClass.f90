!
!//////////////////////////////////////////////////////
!
! Base class for the Jacobian computers
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module JacobianComputerClass

   use SMConstants
   use HexMeshClass
   use PhysicsStorage 
   use MPI_Process_Info                , only: MPI_Process
   use Utilities                       , only: QsortWithFriend
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f
   use GenericMatrixClass              , only: Matrix_t
   use ParamfileRegions                , only: readValueInRegion
   use DenseBlockDiagonalMatrixClass   , only: DenseBlockDiagMatrix_t
   use SparseBlockDiagonalMatrixClass  , only: SparseBlockDiagMatrix_t
   use FTValueDictionaryClass
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none 

   private
   public Look_for_neighbour, ijk2local, local2ijk
   public JacobianComputer_t, GetJacobianFlag
   
!  ------------------
!  Main Jacobian type
!  ------------------
   type :: JacobianComputer_t
      logical              :: preAllocate = .FALSE.
      logical              :: verbose
      integer, allocatable :: ndofelm_l(:)   ! Array containing the number of degrees of freedom for all the LOCAL elements
      integer, allocatable :: ndofelm (:)    ! Array containing the number of degrees of freedom for all the elements
      integer, allocatable :: firstIdx(:)    ! Array containing the first index of the matrix for every element (supposing ordered numbering)
      integer, allocatable :: globIDs_l(:)
      integer              :: nEqn
      contains
         procedure :: Construct           => Jacobian_Construct
         procedure :: Configure           => Jacobian_Configure
         procedure :: GetCSRAllocVectors  => Jacobian_GetCSRVectorsForAllocation
         procedure :: Destruct            => Jacobian_Destruct
         procedure :: Compute             => Jacobian_Compute
   end type JacobianComputer_t
   
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------
!     GetJacobianFlag:
!     Returns the Jacobian flag that is defined in the control file.
!     If not found, returns 0
!     ------------------------------------
      integer function GetJacobianFlag()
         implicit none
         !-local-variables----------------------
         integer, allocatable       :: flag
         character(len=LINE_LENGTH) :: in_label
         character(len=LINE_LENGTH) :: paramFile
         !--------------------------------------
         
         write(in_label , '(A)') "#define jacobian"
         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "type"  , flag     , in_label , "# end" )
         
         if ( allocated(flag) ) then
            GetJacobianFlag = flag
            deallocate(flag)
         else
            GetJacobianFlag = 0
         end if
         
      end function GetJacobianFlag
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------
!     Construct the JacobianInfo variables
!     ------------------------------------
      subroutine Jacobian_Construct(this, mesh, nEqn, controlVariables)
         implicit none
         !-arguments-----------------------------------------
         class(JacobianComputer_t), intent(inout) :: this
         type(HexMesh)            , intent(inout) :: mesh
         integer                  , intent(in)    :: nEqn
         type(FTValueDictionary)  , intent(in)    :: controlVariables
         !-local-variables-----------------------------------
         integer :: eID, ierr
         integer :: all_globIDs(mesh % no_of_allElements)
         integer :: all_ndofelm(mesh % no_of_allElements)
         integer :: no_of_elems(MPI_Process % nProcs)
         integer :: displs     (MPI_Process % nProcs)
         logical, allocatable       :: verbose_in, prealloc_in
         character(len=LINE_LENGTH) :: in_label
         character(len=LINE_LENGTH) :: paramFile
         !---------------------------------------------------
         
         this % nEqn = nEqn
         
!
!        ***************
!        Allocate memory
!        ***************
!
         allocate ( this % ndofelm_l(mesh % no_of_elements       ) )
         allocate ( this % ndofelm  (mesh % no_of_allElements    ) )
         allocate ( this % firstIdx (mesh % no_of_allElements + 1) )
         allocate ( this % globIDs_l(mesh % no_of_elements       ) )
!
!        *********************
!        Get local block sizes
!        *********************
!
         do eID=1, mesh % no_of_elements
            this % ndofelm_l(eID)  = nEqn * (mesh % Nx(eID)+1) * (mesh % Ny(eID)+1) * (mesh % Nz(eID)+1)
         end do
!
!        ***********
!        MPI sharing
!        ***********
!
         if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
!
!           Get global element IDs in all partitions
!           ----------------------------------------
            do eID=1, mesh % no_of_elements
               this % globIDs_l(eID) = mesh % elements(eID) % globID
            end do
!
!           Share info with other processes
!           -------------------------------

            call mpi_allgather(mesh % no_of_elements, 1, MPI_INT, no_of_elems, 1, MPI_INT, MPI_COMM_WORLD, ierr)
            
            displs(1) = 0
            do eID = 1, MPI_Process % nProcs-1
               displs(eID+1) = displs(eID) + no_of_elems(eID)
            end do
            
            call mpi_allgatherv(this % globIDs_l, mesh % no_of_elements, MPI_INT, all_globIDs, no_of_elems , displs, MPI_INT, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(this % ndofelm_l, mesh % no_of_elements, MPI_INT, all_ndofelm, no_of_elems , displs, MPI_INT, MPI_COMM_WORLD, ierr)
!
!           Reorganize the ndofelm
!           ----------------------
            call QsortWithFriend(all_globIDs, all_ndofelm)
            
            this % ndofelm = all_ndofelm
#endif
         else
            this % ndofelm = this % ndofelm_l
         end if
         
!
!        ************
!        Get firstIdx
!        ************
!
         this % firstIdx(1) = 1
         do eID=1, mesh % no_of_allElements
            this % firstIdx(eID+1) = this % firstIdx(eID) + this % ndofelm(eID)
         end do
         
!        Read block
!        **********
         write(in_label , '(A)') "#define jacobian"
         call get_command_argument(1, paramFile) !
         call readValueInRegion ( trim ( paramFile )  , "print info"  , verbose_in     , in_label , "# end" )
         call readValueInRegion ( trim ( paramFile )  , "preallocate" , prealloc_in    , in_label , "# end" )
         
         if ( allocated(verbose_in) ) then
            this % verbose = MPI_Process % isRoot .and. verbose_in
         else
            this % verbose = MPI_Process % isRoot
         end if
         
         if ( allocated(prealloc_in) ) then
            this % preAllocate = prealloc_in
         end if
         
      end subroutine Jacobian_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------
!     Jacobian_Configure:
!     Final setup of the JacobianComputer:
!     -> Preallocation of system matrix (if necessary)
!     ------------------------------------
      subroutine Jacobian_Configure(this, mesh, nEqn, Matrix)
         implicit none
         !-arguments-----------------------------------------
         class(JacobianComputer_t), intent(inout) :: this
         type(HexMesh)            , intent(in)    :: mesh
         integer                  , intent(in)    :: nEqn
         class(Matrix_t)          , intent(inout) :: Matrix
         !-local-variables-----------------------------------
         integer      , allocatable :: rows(:), cols(:)
         real(kind=RP), allocatable :: vals(:)
         !---------------------------------------------------
         
!        Preallocate matrix if requested
!        -------------------------------
         if (this % preAllocate) then
            
            select type(Matrix_p => Matrix)
!
!              If block-diagonal matrix, construct with size of blocks
!              -------------------------------------------------------
               type is(DenseBlockDiagMatrix_t)
                  call Matrix_p % Preallocate(nnzs=this % ndofelm_l)
               type is(SparseBlockDiagMatrix_t)
                  call Matrix_p % Preallocate(nnzs=this % ndofelm_l)
               class default
                  call this % GetCSRAllocVectors (mesh, nEqn, rows, cols) 
                  allocate ( vals(size(cols)) )
                  vals = 0._RP
                  call Matrix_p % constructWithCSRArrays(rows,cols,vals)
                  
                  deallocate (rows,cols,vals)
            end select
            
            call Matrix % SpecifyBlockInfo(this % firstIdx,this % ndofelm)
         end if
         
      end subroutine Jacobian_Configure
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------
!     Destruct the JacobianInfo variables
!     ------------------------------------
      subroutine Jacobian_Destruct(this)
         implicit none
         !-arguments-----------------------------------------
         class(JacobianComputer_t), intent(inout) :: this
         !---------------------------------------------------
         
         deallocate (this % ndofelm)
         deallocate (this % firstIdx)
         
      end subroutine Jacobian_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------
!     Generic subroutine for computing the Jacobian matrix
!     ----------------------------------------------------
      subroutine Jacobian_Compute(this, sem, nEqn, time, Matrix, TimeDerivative, TimeDerivativeIsolated, eps_in, BlockDiagonalized, mode)
         implicit none
         !-arguments----------------------------------
         class(JacobianComputer_t)        , intent(inout)     :: this
         type(DGSem)              , intent(inout)     :: sem
         integer,                   intent(in)        :: nEqn               ! TODO:  Deprecate
         real(kind=RP)            , intent(in)        :: time
         class(Matrix_t)          , intent(inout)     :: Matrix
         procedure(ComputeTimeDerivative_f), optional :: TimeDerivative
         procedure(ComputeTimeDerivative_f), optional :: TimeDerivativeIsolated 
         real(kind=RP)  , optional, intent(in)        :: eps_in
         logical        , optional, intent(in)        :: BlockDiagonalized  !<? Construct only the block diagonal? (Only for AnJacobian_t)
         integer        , optional, intent(in)        :: mode
         !--------------------------------------------
         
         error stop 'JacobianComputer_t must be cast to an extended type (e.g. AnJacobian_r or NumJacobian_t) for computation'
         
      end subroutine Jacobian_Compute
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------------------
!     Routine to get the CSR vectors to preallocate a Jacobian matrix
!     -> Should this be moved to discretization?
!     -> Only for NAVIERSTOKES (temporarily: this can be made general for any advection/diffusion physics)
!     -> Only for compact viscous numerical fluxes
!     ---------------------------------------------------------------
      subroutine Jacobian_GetCSRVectorsForAllocation(this, mesh, nEqn, rows, cols)
         use IntegerDataLinkedList           , only: IntegerDataLinkedList_t
         use ElementConnectivityDefinitions  , only: axisMap, normalAxis
         use MeshTypes                       , only: indexesOnOtherFace
         use NodalStorageClass               , only: GAUSSLOBATTO
         implicit none
         !-arguments----------------------------------
         class(JacobianComputer_t)        , intent(in)     :: this
         type(HexMesh)            , intent(in)     :: mesh
         integer                  , intent(in)     :: nEqn
         integer, allocatable     , intent(inout)  :: rows(:)
         integer, allocatable     , intent(inout)  :: cols(:)
         !-local-variables----------------------------
         type(IntegerDataLinkedList_t), allocatable :: colList(:) ! A linked list containing the cols of the nonzero values for every DOF
         integer, allocatable :: colArray(:)
         integer :: eID
         integer :: i1, j1, k1   ! row degrees of freedom
         integer :: i2, j2, k2   ! col degrees of freedom
         integer :: eq           ! Equation counter
         integer :: thisdof      ! Number of SPATIAL DOF
         integer :: counter      ! Multi-use counter
         integer :: nXi_d, nEta_d         ! Number of nodes in every direction (diagonal blocks)
         integer :: EtaSpa_d, ZetaSpa_d   ! Spacing for these two coordinate directions (diagonal blocks)
         integer :: nXi_od, nEta_od       ! Number of nodes in every direction (diagonal blocks)
         integer :: EtaSpa_od, ZetaSpa_od ! Spacing for these two coordinate directions (diagonal blocks)
         integer :: deltas          ! Number of deltas for a specific dof combination
         integer :: baseCol         ! Base column for a specific DOF
         integer :: nnz             ! Total number of nonzeros
         integer :: nnz_row         ! Number of nonzeros for specific row
         integer :: elSide          ! Side of element
         integer :: elInd_plus(3)            ! Element indexes on e⁺
         integer :: elInd_minus(3)           ! Element indexes on e⁻
         integer :: faceInd_plus (2)         ! Face indexes on e⁺
         integer :: faceInd_minus(2)         ! Face indexes on e⁻
         integer :: tanAx_plus (2)           ! Tangent axes to f on e⁺
         integer :: tanAx_minus(2)           ! Tangent axes to f on e⁻
         integer :: faceInd_minus2plus(2)    ! Face indexes on e⁻ passed to the reference frame of e⁺ (only needed for viscous fluxes)
         integer :: NxyFace_minus(2)         ! Polynomial orders of face on element e⁻ (only needed for viscous fluxes)
         integer :: side                     ! face side
         integer :: normAx_minus
         integer :: normAxInd_minus
         integer :: normAx_plus
         integer :: normAxInd_plus
         integer, parameter :: other(2) = [2, 1]
         !--------------------------------------------
#ifdef NAVIERSTOKES
!        Initializations
!        ***************
         ! Allocations
         allocate (rows (mesh % NDOF * nEqn + 1) )
         
         ! Construct linked lists
         allocate( colList(mesh % NDOF) )
         do thisdof=1, mesh % NDOF
            colList(thisdof) = IntegerDataLinkedList_t(.FALSE.)
         end do
         
         ! Other initializations
         thisdof = 0
         
!        Get matrix sparsity
!        *******************
         do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
            nXi_d     = e % Nxyz(1) + 1
            nEta_d    = e % Nxyz(2) + 1
            EtaSpa_d  = nEqn * nXi_d
            ZetaSpa_d = nEqn * nXi_d * nEta_d
            
!           Diagonal blocks
!           ---------------
            
            do k1 = 0, e % Nxyz(3) ; do j1 = 0, e % Nxyz(2) ; do i1 = 0, e % Nxyz(1)      ! Loop over dofs (rows)
               thisdof = thisdof + 1
               
               do k2 = 0, e % Nxyz(3) ; do j2 = 0, e % Nxyz(2) ; do i2 = 0, e % Nxyz(1)
                  ! Get number of deltas
                  deltas = 0
                  if (i1 == i2) deltas = deltas + 1
                  if (j1 == j2) deltas = deltas + 1
                  if (k1 == k2) deltas = deltas + 1
                  
                  ! Physics specific cycling
                  if (flowIsNavierStokes) then  ! Elliptic terms
                     if (deltas < 1) cycle
                  else                          ! Hyperbolic terms
                     if (deltas < 2) cycle
                  end if
                  
                  ! Add entries to linked list
                  baseCol = (this % firstIdx(e % globID) - 1) + i2*nEqn + j2*EtaSpa_d + k2*ZetaSpa_d
                  do eq = 1, nEqn
                     call colList(thisdof) % Add (baseCol + eq)
                  end do
                  
               end do                 ; end do                 ; end do
               
!              Off-diagonal blocks (it is not super efficient to do this here... but it's easy)
!              --------------------------------------------------------------------------------
            
               do elSide=1, 6
                  if (e % NumberOfConnections(elSide) == 0) cycle
                  
                  associate ( e_minus => e % Connection(elSide), &
                              f       => mesh % faces(e % faceIDs(elSide) ) )
                  
                  nXi_od     = e_minus % Nxyz(1) + 1
                  nEta_od    = e_minus % Nxyz(2) + 1
                  EtaSpa_od  = nEqn * nXi_od
                  ZetaSpa_od = nEqn * nXi_od * nEta_od
                  side = e % faceSide(elSide)
                  
                  tanAx_plus   = axisMap(:,f % elementSide(    side     ) )
                  tanAx_minus  = axisMap(:,f % elementSide( other(side) ) )
                  NxyFace_minus = e_minus % Nxyz ( tanAx_minus )
                  
                  normAx_minus = normalAxis ( f % elementSide( other(side) ) )
                  if (normAx_minus < 0) then
                     normAxInd_minus = 0
                  else
                     normAxInd_minus = e_minus % Nxyz(abs(normAx_minus))
                  end if
                  normAx_minus = abs(normAx_minus)
                  
                  normAx_plus = normalAxis ( elSide )
                  if (normAx_plus < 0) then
                     normAxInd_plus = 0
                  else
                     normAxInd_plus = e % Nxyz(abs(normAx_plus))
                  end if
                  normAx_plus = abs(normAx_plus)
                  
                  do k2 = 0, e_minus % Nxyz(3) ; do j2 = 0, e_minus % Nxyz(2) ; do i2 = 0, e_minus % Nxyz(1)
                     
                     elInd_plus  = [i1, j1, k1]
                     elInd_minus = [i2, j2, k2]
                     
                     ! Gauss-Lobatto cycling
                     
                     if (flowIsNavierStokes) then  ! Elliptic terms
                        if (mesh % nodeType == GAUSSLOBATTO) then
                           if ( elInd_minus(normAx_minus) /= normAxInd_minus .and. &
                                elInd_plus (normAx_plus)  /= normAxInd_plus       ) cycle  ! TODO: This gives a sparsity pattern that works but they may be too dense!!
                        end if
                     else                          ! Hyperbolic terms
                        if (mesh % nodeType == GAUSSLOBATTO) then
                           if ( elInd_minus(normAx_minus) /= normAxInd_minus .or. &
                                elInd_plus (normAx_plus)  /= normAxInd_plus       ) cycle
                        end if
                     end if
                     
                     faceInd_plus  = elInd_plus ( tanAx_plus  )
                     faceInd_minus = elInd_minus( tanAx_minus )
                     
                     call indexesOnOtherFace(faceInd_minus(1),faceInd_minus(2), NxyFace_minus(1),NxyFace_minus(2), f % rotation, other(side), faceInd_minus2plus(1),faceInd_minus2plus(2))
                     
                     ! Get number of deltas
                     Deltas = 0
                     if ( faceInd_plus(1) == faceInd_minus2plus(1) ) then
                        deltas = deltas + 1
                     end if
                     if ( faceInd_plus(2) == faceInd_minus2plus(2) ) then
                        deltas = deltas + 1
                     end if
                     
                     ! Physics specific cycling
                     if (flowIsNavierStokes) then  ! Elliptic terms
                        if (deltas < 1) cycle
                     else                          ! Hyperbolic terms
                        if (deltas < 2) cycle
                     end if
                     
                     ! Add entries to linked list
                     baseCol = (this % firstIdx(e_minus % globID) - 1) + i2*nEqn + j2*EtaSpa_od + k2*ZetaSpa_od
                     do eq = 1, nEqn
                        call colList(thisdof) % Add (baseCol + eq)
                     end do
                  end do                       ; end do                       ; end do
                  end associate
               end do
            end do                 ; end do                 ; end do
            

            
            end associate
         end do
         
!        Gather info and create arrays
!        *****************************
         
         ! Get number of nonzeros and construct rows array
         nnz = 0
         rows(1) = 1
         counter = 1 ! row number
         do thisdof=1, mesh % NDOF
            nnz_row = colList(thisdof) % no_of_entries
            nnz = nnz + nnz_row
            ! One row for each equation of this DOF
            do eq=1, nEqn
               counter = counter + 1
               rows(counter) = rows(counter-1) + nnz_row
            end do
         end do
         nnz = nnz * nEqn
         
         ! Construct cols array
         allocate ( cols(nnz) )
         counter = 1 ! first nonzero of row
         do thisdof=1, mesh % NDOF
            nnz_row = colList(thisdof) % no_of_entries
            call colList(thisdof) % ExportToArray (colArray,.TRUE.)
            
            ! One row for each equation of this DOF
            do eq=1, nEqn
               cols(counter:counter+nnz_row-1) = colArray
               counter = counter + nnz_row
            end do
            
            deallocate (colArray)
         end do
         
!        Finish up
!        *********
         call colList % destruct
         deallocate (colList)
#else
         error stop 'Jacobian_GetCSRVectorsForAllocation only for NS'
#endif
      end subroutine Jacobian_GetCSRVectorsForAllocation
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Look_for_neighbour(this, mesh)
      IMPLICIT NONE 
!
!     -----------------------------------------------------------------------------
!     This subroutine finds the neighbours of all elements (conforming meshes only)
!        current implementation only finds neighbour 3D elements, but not faces
!     -----------------------------------------------------------------------------
!
!
!     ----------------
!     Input parameters
!     ----------------
!
      TYPE(Neighbor_t)                  :: this(:)
      TYPE(HexMesh)                     :: mesh
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                         :: i,j,iEl

         DO iEl = 1, SIZE(mesh%elements)
            this(iEl)%elmnt(7) = mesh % elements(iEl) % globID  ! The last one is itself
            DO j = 1, 6
               IF (mesh % elements(iEl) % NumberOfConnections(j) == 0) THEN
                  this(iEl)%elmnt(j) = 0
               ELSE
                  this(iEl)%elmnt(j) = mesh % elements(iEl) % Connection(j) % globID 
               ENDIF
            ENDDO
         ENDDO

   END SUBROUTINE 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   FUNCTION ijk2local(i,j,k,l,N_EQN,Nx,Ny,Nz) RESULT(idx)
!  -----------------------------------------------------------------------------------------------------------------------------
!  Returns the local index relative to an element from the local coordinates: i(lagrange node x), j(lagrange node y), 
!  k(lagrange node z), l(equation number)
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
!  -----------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)   :: i, j, k, l, Nx, Ny, Nz, N_EQN
      INTEGER               :: idx
      
      IF (l < 1 .OR. l > N_EQN)  error stop 'error in ijk2local, l has a wrong value'
      IF (i < 0 .OR. i > Nx)     error stop 'error in ijk2local, i has a wrong value'
      IF (j < 0 .OR. j > Ny)     error stop 'error in ijk2local, j has a wrong value'
      IF (k < 0 .OR. k > Nz)     error stop 'error in ijk2local, k has a wrong value'
      
      idx = k*(Nx+1)*(Ny+1)*N_EQN + j*(Nx+1)*N_EQN + i*N_EQN + l
   END FUNCTION
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION local2ijk(idx,N_EQN,Nx,Ny,Nz) RESULT (indices)
!  -----------------------------------------------------------------------------------------------------------------------------
!  Returns the coordinates relative to an element: l(equation number), i(lagrange node x), j(lagrange node y), k(lagrange node z)
!  from the local index  
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
!  -----------------------------------------------------------------------------------------------------------------------------
      INTEGER, INTENT(IN)   :: idx, Nx, Ny, Nz, N_EQN
      INTEGER               :: indices(4)
      INTEGER               :: tmp1, tmp2
      
      IF (idx < 1 .OR. idx > (Nx+1)*(Ny+1)*(Nz+1)*N_EQN) error stop 'error in local2ijk, idx has wrong value'
      
      indices(4) = (idx-1) / ((Nx+1)*(Ny+1) * N_EQN)
      tmp1       = MOD((idx-1),((Nx+1)*(Ny+1) * N_EQN) )
      indices(3) = tmp1 / ((Nx+1)*N_EQN)
      tmp2       = MOD(tmp1,((Nx+1)*N_EQN))
      indices(2) = tmp2 / (N_EQN)
      indices(1) = MOD(tmp2, N_EQN) + 1
   END FUNCTION
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module JacobianComputerClass

