!
!//////////////////////////////////////////////////////
!
!   @File:    MeshPartitioning.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:08 2017
!   @Last revision date: Wed Sep 12 13:12:42 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 06ee9cb1c578e7283cc508c9aacf873836377025
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module MeshPartitioning
   use SMConstants
   use MeshTypes
   use HexMeshClass
   use PartitionedMeshClass
   use FileReadingUtilities            , only: RemovePath, getFileName

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

!
!        Export partitions file
!        ----------------------
         call WritePartitionsFile(mesh, elementsDomain)
      end subroutine PerformMeshPartitioning

      subroutine GetElementsDomain(mesh, no_of_domains, elementsDomain, partitions)
         use IntegerDataLinkedList
         use MPI_Process_Info
         implicit none
         type(HexMesh), intent(in)              :: mesh
         integer,       intent(in)              :: no_of_domains
         integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
         type(PartitionedMesh_t), intent(inout) :: partitions(no_of_domains)      
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, allocatable :: nodesDomain(:)


         allocate(nodesDomain(size(mesh % nodes)))
!
!        **********************************************
!        Set which elements belong to each domain using
!        **********************************************
!
         select case (MPI_Partitioning)
!     
!           Space-filling curve partitioning
!           --------------------------------
            case (SFC_PARTITIONING)
               call GetSFCElementsPartition(no_of_domains, mesh % no_of_elements, elementsDomain)
!     
!           METIS partitioning
!           ------------------
            case (METIS_PARTITIONING)
               call GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain)
         end select
!
!        ****************************************
!        Get which nodes belong to each partition
!        ****************************************
!
         call GetNodesPartition(mesh, no_of_domains, elementsDomain, partitions)   

         deallocate(nodesDomain)   
         
      end subroutine GetElementsDomain

!
!////////////////////////////////////////////////////////////////////////
!
     subroutine GetNodesPartition(mesh, no_of_domains, elementsDomain, partitions)
        use Utilities, only: Qsort
        implicit none
        type(HexMesh), intent(in)              :: mesh
        integer,       intent(in)              :: no_of_domains
        integer,       intent(in)              :: elementsDomain(mesh % no_of_elements)
        type(PartitionedMesh_t), intent(inout) :: partitions(no_of_domains)   
!
!       ---------------
!       Local Variables
!       ---------------
!
        integer              :: nvertex
        integer              :: i
        integer              :: j
        integer              :: k
        integer              :: ipoint
        integer              :: jpoint
        integer              :: idomain
        integer              :: npoints
        integer              :: ielem
        logical              :: isnewpoint
        logical              :: meshIsHOPR
        integer, allocatable :: points(:)
        integer, allocatable :: HOPRpoints(:)

        nvertex = 8
        
        meshIsHOPR = allocated (mesh % HOPRnodeIDs)
        
        do idomain=1,no_of_domains
!
!       Get the number of elements for the partition
!       --------------------------------------------
        partitions(idomain)%no_of_elements = count(elementsDomain == idomain)
        allocate(partitions(idomain)%elementIDs(partitions(idomain)%no_of_elements))
!
!       This will store the partition nodes (allocated as 8 * no_of_elements)
!       ---------------------------------------------------------------------
        allocate(points(nvertex*partitions(idomain)%no_of_elements))
        points = 0
        if (meshIsHOPR) then
            allocate(HOPRpoints(nvertex*partitions(idomain)%no_of_elements))
            HOPRpoints = 0
        end if
!
!       ****************************************
!       Gather each partition nodes and elements      
!       ****************************************
!
        k = 0
        npoints = 0
        do ielem=1,mesh % no_of_elements
           if (elementsDomain(ielem) == idomain) then
!
!             Append a new element
!             --------------------
              k = k + 1
              partitions(idomain)%elementIDs(k) = ielem
!
!             Append its nodes
!             ----------------           
              do j=1,nvertex
!
!                Get the node ID
!                ---------------
                 jpoint = mesh % elements(ielem) % nodeIDs(j)
!
!                Check if it is already stored
!                -----------------------------
                 isnewpoint = .true.
                 do i=1,npoints
                    ipoint = points(i)
                    if (jpoint == ipoint) then
                       isnewpoint = .false.
                       exit
                    end if
                 end do
!
!                Store the node
!                --------------      
                 if (isnewpoint) then
                    npoints = npoints + 1
                    points(npoints) = jpoint
                    if (meshIsHOPR) HOPRpoints(npoints) = mesh % HOPRnodeIDs(jpoint)
                 end if                  
              end do
            end if      
         end do
!
!        Put the nodeIDs into the partitions structure
!        ---------------------------------------------      
         allocate(partitions(idomain)%nodeIDs(npoints))
         partitions(idomain)%nodeIDs(:) = points(1:npoints)
         
         if (meshIsHOPR) then
            allocate(partitions(idomain)%HOPRnodeIDs(npoints))
            partitions(idomain)%HOPRnodeIDs(:) = HOPRpoints(1:npoints)
         end if
!
!        Sort the nodeIDs to read the mesh file accordingly (only needed for SpecMesh)
!        --------------------------------------------------
         if (.not. meshIsHOPR) call Qsort(partitions(idomain)%nodeIDs)

         partitions(idomain)%no_of_nodes = npoints
!
!        ****
!        Free
!        ****
!
         deallocate(points)
         safedeallocate(HOPRpoints)
      end do

     end subroutine GetNodesPartition
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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------
!     Space-filling curve partitioning
!     --------------------------------
      subroutine GetSFCElementsPartition(no_of_domains, nelem, elementsDomain)
         implicit none
         !-arguments--------------------------------------------------
         integer, intent(in)    :: no_of_domains
         integer, intent(in)    :: nelem
         integer, intent(inout) :: elementsDomain(nelem)
         !-local-variables--------------------------------------------
         integer :: elems_per_domain(no_of_domains)
         integer :: biggerdomains
         integer :: first, last, domain
         !------------------------------------------------------------
         
         elems_per_domain = nelem / no_of_domains
         
         biggerdomains = mod(nelem,no_of_domains)
         elems_per_domain(1:biggerdomains) = elems_per_domain(1:biggerdomains) + 1
         
         first = 1
         do domain = 1, no_of_domains
            last = first + elems_per_domain(domain) - 1
            elementsDomain(first:last) = domain
            first = last + 1
         end do
         
      end subroutine GetSFCElementsPartition
      
      subroutine WritePartitionsFile(mesh,elementsDomain)
         implicit none
         !-arguments--------------------------------------------------
         type(HexMesh), intent(in)  :: mesh
         integer      , intent(in)  :: elementsDomain(mesh % no_of_elements)
         !-local-variables--------------------------------------------
         character(LINE_LENGTH)     :: pmeshName
         integer                    :: fID, eID
         !------------------------------------------------------------
         
         pmeshName = "./MESH/" // trim(removePath(getFileName(mesh % meshFileName))) // ".pmesh"
         
         open(newunit = fID, file=trim(pmeshName),action='write')
         
         write(fID,*) mesh % no_of_elements
         do eID = 1, mesh % no_of_elements
            write(fID,*) elementsDomain(eID)
         end do
            
         close(fID)
         
      end subroutine WritePartitionsFile
         
end module MeshPartitioning
