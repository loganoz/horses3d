!
!//////////////////////////////////////////////////////
!
!   @File:    METISPartitioning.f90
!   @Author:  Mariola Gomez and Marta Cordero (marta.cordero@upm.es / mariola.gomez@upm.es)
!   @Created: Wed Dec 20 19:57:13 2017
!   @Last revision date: Wed Dec 20 20:22:04 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 6b682d113059e1a767187182f6324a94681ce47a
!
!//////////////////////////////////////////////////////
!
   subroutine GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain)
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
      integer,       intent(out)             :: elementsDomain(mesh % no_of_allElements)
      integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
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
#ifdef _HAS_MPI_
!
!     **************
!     Initialization    
!     **************
!
      ne = mesh % no_of_allElements
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
