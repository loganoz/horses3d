!
!//////////////////////////////////////////////////////
!
!   @File:    MeshPartitioning.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:08 2017
!   @Last revision date: Sat Nov 25 14:01:19 2017
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 3b34c7a95eb684f3c89837d445c6430f5449f298
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
         implicit none
         type(HexMesh), intent(in)          :: mesh
         integer,       intent(in)          :: no_of_domains
         integer,       intent(out)         :: elementsDomain(mesh % no_of_elements)
         type(PartitionedMesh_t), intent(out) :: partitions(no_of_domains)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: domain, eID, nID 
         type(IntegerDataLinkedList_t) :: elements(no_of_domains)
         type(IntegerDataLinkedList_t) :: nodes(no_of_domains)
         real(kind=RP)  :: xc(3)
!
!        Initialize the linked lists         
!        ---------------------------
         do domain = 1, no_of_domains
            elements(domain) = IntegerDataLinkedList_t(allowRepetitions = .false.)
            nodes(domain) = IntegerDataLinkedList_t(allowRepetitions = .false.)
         end do
!
!        Loop in elements and fill the domains
!        -------------------------------------         
         do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
!
!           Get the centroid
!           ----------------
            xc = 0.0_RP
            do nID = 1, 8
               xc = xc + mesh % nodes(e % nodeIDs(nID)) % x
            end do
            xc = xc / 8.0_RP
!
!           Inquire to which domain it belongs (only valid for the TGV with 8 domains)
!           ----------------------------------
            xc = xc - PI
          
            if ( xc(2) .lt. 0 ) then
               domain = 1
            else
               domain = 2
            end if

            if ( xc(1) .lt. 0 ) domain = domain + 2
            if ( xc(3) .gt. 0 ) domain = domain + 4
!
!           Add elements and nodes to the domain linked list
!           ------------------------------------------------
            elementsDomain(eID) = domain
            call elements(domain) % Add(eID)
            do nID = 1, 8
               call nodes(domain) % Add(e % nodeIDs(nID))
            end do
            end associate           
         end do
!
!        ----------------------------------
!        Build the MeshPartition structures      
!        ----------------------------------
!
         do domain = 1, no_of_domains
            partitions(domain) % no_of_nodes = nodes(domain) % no_of_entries      
            partitions(domain) % no_of_elements = elements(domain) % no_of_entries
            
            allocate(partitions(domain) % nodeIDs( partitions(domain) % no_of_nodes ))
            allocate(partitions(domain) % elementIDs( partitions(domain) % no_of_elements ))

            call nodes(domain)    % ExportToArray(partitions(domain) % nodeIDs)
            call elements(domain) % ExportToArray(partitions(domain) % elementIDs)
         end do
         
      end subroutine GetElementsDomain

      subroutine GetPartitionBoundaryFaces(mesh, no_of_domains, elementsDomain, partitions)
         implicit none
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_domains
         integer,       intent(in)  :: elementsDomain(no_of_domains)
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
            partitions(dL) % no_of_bdryfaces = partitions(dL) % no_of_bdryfaces + 1
            partitions(dR) % no_of_bdryfaces = partitions(dR) % no_of_bdryfaces + 1

            end associate
            end associate
         end do
!
!        **************************************
!        Allocate boundary faces-related memory
!        **************************************
!
         do domain = 1, no_of_domains
            associate(nFaces => partitions(domain) % no_of_bdryfaces)
            allocate(partitions(domain) % bdryface_elements(nFaces))
            allocate(partitions(domain) % element_bdryfaceSide(nFaces))
            allocate(partitions(domain) % bdryface_rotation(nFaces))
            allocate(partitions(domain) % bdryface_elementSide(nFaces))
            allocate(partitions(domain) % bdryface_sharedDomain(nFaces))
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
            partitions(dL) % bdryface_elements(bfaceID(dL)) = eL % eID
            partitions(dR) % bdryface_elements(bfaceID(dR)) = eR % eID
!
!           Get the face sides in the elements
!           ----------------------------------
            partitions(dL) % element_bdryfaceSide(bfaceID(dL)) = f % elementSide(1)
            partitions(dR) % element_bdryfaceSide(bfaceID(dR)) = f % elementSide(2)
!
!           Get the face rotation
!           ---------------------
            partitions(dL) % bdryface_rotation(bfaceID(dL)) = f % rotation 
            partitions(dR) % bdryface_rotation(bfaceID(dR)) = f % rotation 
!
!           Get the element face side
!           -------------------------
            partitions(dL) %  bdryface_elementSide(bfaceID(dL)) = 1
            partitions(dR) %  bdryface_elementSide(bfaceID(dR)) = 2
!
!           Get the shared domain
!           ---------------------
            partitions(dL) % bdryface_sharedDomain(bfaceID(dL)) = dR
            partitions(dR) % bdryface_sharedDomain(bfaceID(dR)) = dL

            end associate
            end associate
         end do

      end subroutine GetPartitionBoundaryFaces

end module MeshPartitioning
