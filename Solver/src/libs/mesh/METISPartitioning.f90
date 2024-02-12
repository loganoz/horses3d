   subroutine GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain, useWeights)
!
!     *********************************************************************
!        This subroutine performs the mesh partitioning using METIS
!     *********************************************************************
!        
      use HexMeshClass
      use SMConstants
      implicit none
      type(HexMesh), intent(in)              :: mesh
      integer,       intent(in)              :: no_of_domains
      integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
      integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
      logical,       intent(in)              :: useWeights
!
!     ---------------
!     Local Variables
!     ---------------
!
      integer                :: i
      integer                :: ielem
      integer                :: ne                    ! # elements
      integer                :: nn                    ! # nodes
      integer                :: nvertex               ! # vertices per element
      integer, allocatable   :: eptr(:)               ! index in eind --> eptr(ne+1)
      integer, allocatable   :: eind(:)               ! vertices of each element   --> eind(nvertex*ne)
      integer, pointer       :: vwgt(:) => null()     ! vertices weights
      integer, pointer       :: vsize(:) => null()    !
      integer, parameter     :: ncommon = 4           ! common faces for dual nesh (4 for hexahedrals)
      real(kind=RP), pointer :: tpwgt(:) => null()    ! partitions' weights --> tpwgt(no_of_domains)
      integer, pointer       :: opts(:) => null()     ! metis options
      integer, parameter     :: metis_noptions = 39   ! number of metis options
      integer                :: objval                ! objective calculated value
      integer, parameter     :: MIN_EDGE_CUT = -1     ! option to minimize edge-cut
      integer, parameter     :: MIN_COMM_VOL = 1      ! option to minimize the communication volume
      integer, allocatable, target   :: weights(:)
#ifdef _HAS_METIS_
!
!     **************
!     Initialization    
!     **************
!
      ne = mesh % no_of_elements
      nn = size(mesh % nodes)
      nvertex = 8

      allocate(eptr(ne+1))
      allocate(eind(nvertex*ne))
!
!     Gather each element nodes: they are stored using C indexes (starting on 0)
!     -------------------------
      i = 1
      do ielem=1,ne
!
!        Save each element nodes
!        -----------------------
         eind(i:i+nvertex-1) = mesh % elements(ielem) % nodeIDs - 1
!
!        Save each element ID
!        --------------------
         eptr(ielem) = i - 1
         i = i + nvertex
      end do
!
!     Termination: set the last element position in the N+1 entry
!     -----------
      eptr(ne+1) = i - 1
!
!     *****************
!     Set METIS options
!     *****************
!
      allocate(opts(0:metis_noptions))
      call METIS_SetDefaultOptions(opts)
!
!     First option chooses the method: edge-cut / communication volume      
!     -------------------------------
      opts(1) = MIN_EDGE_CUT
!
!     Disable verbosity
!     -----------------
      opts(5) = 0
!
!     *******************************
!     Calculate weights based on NDOF
!     *******************************
!
      if (useWeights) then
          allocate(weights(ne))
          do ielem=1,ne
              weights(ielem) = product(mesh % elements(ielem) % Nxyz + 1)
          end do
          ! weights(ne+1) = product(mesh % elements(ielem) % Nxyz + 1)
          ! eptr(ne+1) = i - 1
          if (maxval(weights) .ne. minval(weights)) then
              vwgt => weights
          endif
      end if 
!     **********************
!     Perform the partitions
!     **********************
!
      call METIS_PartMeshDual(ne, nn,  eptr, eind,  vwgt, vsize, ncommon, no_of_domains, tpwgt,  opts,  objval, elementsDomain, nodesDomain)
!
!     Recover FORTRAN displacements by adding 1
!     -----------------------------------------
      elementsDomain = elementsDomain + 1
      nodesDomain = nodesDomain + 1
!
!     ****
!     Free
!     ****
!
      deallocate(eptr)
      deallocate(eind)
#endif
   end subroutine GetMETISElementsPartition
