!
!//////////////////////////////////////////////////////
!
!   @File:    MeshPartitioning.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:08 2017
!   @Last revision date: Tue Dec 19 16:59:00 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 388a9acc3eac24b19e7273b113b07639cafbd3cf
!
!//////////////////////////////////////////////////////
!
module MeshPartitioning
   use SMConstants
   use HexMeshClass
   use PartitionedMeshClass

   private
   public   PerformMeshPartitioning

   contains
      subroutine PerformMeshPartitioning(mesh, no_of_domains, partitions)
         implicit none
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_domains
         type(PartitionedMesh_t)    :: partitions(no_of_domains)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer               :: fID, domain
         integer               :: elementsDomain(mesh % no_of_elements)
!
!        ************************************************
!        Now, they will just be ordered in a hard-coded
!        way, and just valid for the TaylorGreen geometry
!        It is required to consider using METIS or a 
!        universal partitioner later.
!        ************************************************
!
!        Initialize partitions
!        ---------------------
         do domain = 1, no_of_domains
            partitions(domain) = PartitionedMesh_t(domain)
         end do
!
!        Get each domain elements and nodes
!        ----------------------------------
         call GetElementsDomain(mesh, no_of_domains, elementsDomain, partitions)
!
!        Get the partition boundary faces
!        --------------------------------
         call GetPartitionBoundaryFaces(mesh, no_of_domains, elementsDomain, partitions)

      end subroutine PerformMeshPartitioning

      subroutine GetElementsDomain(mesh, no_of_domains, elementsDomain, partitions)
!
!        *******************************************
!        Here is where the METIS partitioner would
!        divide the mesh.
!        *******************************************
!
         use IntegerDataLinkedList
         use MPI_Process_Info
         implicit none
         type(HexMesh), intent(in)              :: mesh
         integer,       intent(in)              :: no_of_domains
         integer,       intent(out)             :: elementsDomain(mesh % no_of_allElements)
         type(PartitionedMesh_t), intent(inout) :: partitions(no_of_domains)      
!
!        ---------------
!        Local variables
!        ---------------
!
       integer, allocatable :: nodesDomain(:)

!--- End of header ------------------------------------------------------

       allocate(nodesDomain(size(mesh % nodes)))

       ! call METIS to set elementsDomains
       call GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain)

       ! set partitions
       call GetNodesPartition(mesh, no_of_domains, elementsDomain, nodesDomain, partitions)   

       deallocate(nodesDomain)   
         
      end subroutine GetElementsDomain

!
!////////////////////////////////////////////////////////////////////////
!
     subroutine GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain)
!
!      ---------
!      Arguments
!      ---------
!
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
         integer                :: ne               ! # elements
       integer                :: nn               ! # nodes
       integer                :: nvertex            ! # vertices per element
       integer, allocatable   :: eptr(:)               ! index in eind --> eptr(ne+1)
       integer, allocatable   :: eind(:)            ! vertices of each element   --> eind(nvertex*ne)
       integer, pointer       :: vwgt(:) => null()   ! vertices weights
       integer, pointer       :: vsize(:) => null()   !
       integer                :: ncommon            ! common faces for dual nesh
       real(kind=RP), pointer :: tpwgt(:) => null()   ! partitions' weights --> tpwgt(no_of_domains)
       integer, pointer       :: opts(:) => null()   ! metis options
       integer, parameter     :: metis_noptions = 39

       ! output METIS variables
       integer              :: objval               ! objective calculated value
     
!--- End of header ------------------------------------------------------

       ne = mesh % no_of_allElements
        nn = size(mesh % nodes)
       nvertex = 8

       allocate(eptr(ne+1))
       allocate(eind(nvertex*ne))

        ! C-index   
       i = 1
       do ielem=1,ne
           eind(i:i+nvertex-1) = mesh % elements(ielem) % nodeIDs - 1
          eptr(ielem) = i - 1
          i = i + nvertex
       end do
       eptr(ne+1) = i - 1

       allocate(opts(0:metis_noptions))
       call METIS_SetDefaultOptions(opts)
       opts(1) = 1                         ! OBJTYPE: -1  minimizing edge-cut | 1 minimizing communication volume
                                          ! TODO fichero parametros

       ! opts(5) = 1     to enable verbosity 

       ncommon = 4                     ! for hexaeder elements
       call METIS_PartMeshDual(ne, nn,  eptr, eind,  vwgt, vsize, ncommon, no_of_domains, tpwgt,  opts,  objval, elementsDomain, nodesDomain)

       ! rectify idomain
       elementsDomain = elementsDomain + 1
       nodesDomain = nodesDomain + 1

       deallocate(eptr)
       deallocate(eind)

     end subroutine GetMETISElementsPartition
!
!////////////////////////////////////////////////////////////////////////
!
     subroutine GetNodesPartition(mesh, no_of_domains, elementsDomain, nodesDomain, partitions)
!
!      ---------
!      Arguments
!      ---------
!
         implicit none
         type(HexMesh), intent(in)              :: mesh
       integer,       intent(in)              :: no_of_domains
         integer,       intent(in)              :: elementsDomain(mesh % no_of_allElements)
         integer,       intent(in)              :: nodesDomain(size(mesh % nodes))
         type(PartitionedMesh_t), intent(inout) :: partitions(no_of_domains)   
!
!     ---------------
!     Local Variables
!     ---------------
!
        integer :: nvertex
      integer :: i
      integer :: j
      integer :: k
      integer :: ipoint
      integer :: jpoint
      integer :: idomain
      integer :: npoints
      integer :: ielem
      logical :: isnewpoint
      integer, allocatable :: points(:)

!--- End of header ------------------------------------------------------
   
       nvertex = 8

      do idomain=1,no_of_domains

         ! carga nelements
         partitions(idomain)%no_of_elements = count(elementsDomain == idomain)

         ! alocata elementIDs
         allocate(partitions(idomain)%elementIDs(partitions(idomain)%no_of_elements))

         ! dummy variable
         allocate(points(nvertex*partitions(idomain)%no_of_elements))
         points = 0
      
         ! set the elements and the nodes of the domain
         k = 0
         npoints = 0
         do ielem=1,mesh % no_of_allElements

            if (elementsDomain(ielem) == idomain) then
               k = k + 1
               partitions(idomain)%elementIDs(k) = ielem                  ! set the element
            
               ! recorre los nodos de ese elemento para ver si ya esta introducido o no
               do j=1,nvertex

                  jpoint = mesh % elements(ielem) % nodeIDs(j)

                  isnewpoint = .true.
                  do i=1,npoints
                     ipoint = points(i)
                     if (jpoint == ipoint) then
                        isnewpoint = .false.
                        exit
                     end if
                  end do
      
                  if (isnewpoint) then
                     npoints = npoints + 1
                     points(npoints) = jpoint
                  end if                  

               end do

            end if      
         end do
      
         allocate(partitions(idomain)%nodeIDs(npoints))

         ! Fill the nodes vector
         partitions(idomain)%nodeIDs(:) = points(1:npoints)

         ! sort the vector
         call sort(partitions(idomain)%nodeIDs)

         partitions(idomain)%no_of_nodes = npoints

         deallocate(points)

      end do

     end subroutine GetNodesPartition
!
!////////////////////////////////////////////////////////////////////////
!
     subroutine sort(U)
!
!      ---------
!      Arguments
!      ---------
!
       integer, intent(inout) :: U(:)
!
!        ---------------
!        Local variables
!        ---------------
!
       integer :: i
       integer :: n
       integer :: temp
       integer :: v(1)

!--- End of header ------------------------------------------------------

       n = size(U)
       do i=1,n
         v = minloc(U(i:n)) + i - 1
         temp = U(i)
         U(i) = U(v(1))
         U(v(1)) = temp
       end do

     end subroutine sort
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine GetPartitionBoundaryFaces(mesh, no_of_domains, elementsDomain, partitions)
 
         implicit none
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_domains
         integer,       intent(in)  :: elementsDomain(mesh % no_of_elements)
         type(PartitionedMesh_t)    :: partitions(no_of_domains)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fID, eL, eR, dL, dR, domain
         integer  :: bfaceID(no_of_domains)
!
!        *******************************************
!        Get the number of boundary faces per domain 
!        *******************************************
!
         do fID = 1, size(mesh % faces)
            associate(f => mesh % faces(fID))
!
!           Cycle non-interior faces
!           ------------------------
            if ( f % faceType .ne. HMESH_INTERIOR ) cycle
!
!           Create references to left and right elements
!           --------------------------------------------
            associate( eL => mesh % elements(f % elementIDs(1)), &
                       eR => mesh % elements(f % elementIDs(2)))
!
!           Get each elements domain
!           ------------------------
            dL = elementsDomain(eL % eID)
            dR = elementsDomain(eR % eID)
!
!           Cycle if both elements belong to the same domain
!           ------------------------------------------------
            if ( dL .eq. dR ) cycle 
!
!           Otherwise, the face is a domain boundary face for domains dL and dR
!           -------------------------------------------------------------------
            partitions(dL) % no_of_mpifaces = partitions(dL) % no_of_mpifaces + 1
            partitions(dR) % no_of_mpifaces = partitions(dR) % no_of_mpifaces + 1

            end associate
            end associate
         end do
!
!        **************************************
!        Allocate boundary faces-related memory
!        **************************************
!
         do domain = 1, no_of_domains
            associate(nFaces => partitions(domain) % no_of_mpifaces)
            allocate(partitions(domain) % mpiface_elements(nFaces))
            allocate(partitions(domain) % element_mpifaceSide(nFaces))
            allocate(partitions(domain) % mpiface_rotation(nFaces))
            allocate(partitions(domain) % mpiface_elementSide(nFaces))
            allocate(partitions(domain) % mpiface_sharedDomain(nFaces))
            end associate
         end do
!
!        ***************************
!        Get each boundary face data
!        ***************************
!
         bfaceID = 0
         do fID = 1, size(mesh % faces)
            associate(f => mesh % faces(fID))
!
!           Cycle non-interior faces
!           ------------------------
            if ( f % faceType .ne. HMESH_INTERIOR ) cycle
!
!           Create references to left and right elements
!           --------------------------------------------
            associate( eL => mesh % elements(f % elementIDs(1)), &
                       eR => mesh % elements(f % elementIDs(2)))
!
!           Get each elements domain
!           ------------------------
            dL = elementsDomain(eL % eID)
            dR = elementsDomain(eR % eID)
!
!           Cycle if both elements belong to the same domain
!           ------------------------------------------------
            if ( dL .eq. dR ) cycle 
!
!           Otherwise, the face is a domain boundary face for domains dL and dR
!           -------------------------------------------------------------------
            bfaceID(dL) = bfaceID(dL) + 1 
            bfaceID(dR) = bfaceID(dR) + 1 
!
!           Get the elements
!           ----------------
            partitions(dL) % mpiface_elements(bfaceID(dL)) = eL % eID
            partitions(dR) % mpiface_elements(bfaceID(dR)) = eR % eID
!
!           Get the face sides in the elements
!           ----------------------------------
            partitions(dL) % element_mpifaceSide(bfaceID(dL)) = f % elementSide(1)
            partitions(dR) % element_mpifaceSide(bfaceID(dR)) = f % elementSide(2)
!
!           Get the face rotation
!           ---------------------
            partitions(dL) % mpiface_rotation(bfaceID(dL)) = f % rotation 
            partitions(dR) % mpiface_rotation(bfaceID(dR)) = f % rotation 
!
!           Get the element face side
!           -------------------------
            partitions(dL) %  mpiface_elementSide(bfaceID(dL)) = 1
            partitions(dR) %  mpiface_elementSide(bfaceID(dR)) = 2
!
!           Get the shared domain
!           ---------------------
            partitions(dL) % mpiface_sharedDomain(bfaceID(dL)) = dR
            partitions(dR) % mpiface_sharedDomain(bfaceID(dR)) = dL

            end associate
            end associate
         end do

      end subroutine GetPartitionBoundaryFaces

end module MeshPartitioning
