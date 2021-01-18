!
!////////////////////////////////////////////////////////////////////////
!
!   @File:    Read_GMSH.f90
!   @Author:  Wojciech Laskowski (wj.laskowski@upm.es)
!   @Created: Thu Dec 10 10:05:22 2020
!   @Last revision date: Mon Jan 18 14:52:59 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: d6a6ff4a9a1ef873f3a649179fd701c9720fbe91
!
!//////////////////////////////////////////////////////
!
!  Module for reading hexahedral conforming meshes in GMSH (https://gmsh.info/) mesh format.
!  Current supported versions: GMSH Mesh Format 4.1.
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE Read_GMSH
      use SMConstants
      use MeshTypes
      use ElementConnectivityDefinitions
      USE TransfiniteMapClass
      use MappedGeometryClass
      use NodeClass
      use ElementClass
      use HexMeshClass
      use sharedBCModule
      use PhysicsStorage
      use FileReadingUtilities      , only: getFileName, getRealArrayFromStringNoCommas
      use Utilities, only: UnusedUnit, toLower
      implicit none
      
      private
      public ConstructMesh_FromGMSHFile_, NumOfElems_GMSH

!
!  ------------------------------------------------
!  Local temporary element storage.
!  ------------------------------------------------
!
   type :: MSH_BCinfo_t
!-----Variables-----------------------------------------------------------
      integer                             :: dim
      integer                             :: tag
      character(len=127)                  :: name
      integer                             :: no_of_elements
      integer                             :: no_of_nodes
      integer                             :: no_of_surfaces
      integer                             :: no_of_curves
      integer                             :: no_of_points
      integer, dimension(:), allocatable  :: element_tags
      integer, dimension(:), allocatable  :: wall_side
      integer, dimension(:), allocatable  :: node_tags
      integer, dimension(:), allocatable  :: volume_tags
      integer, dimension(:), allocatable  :: surface_tags
      integer, dimension(:), allocatable  :: curve_tags
      integer, dimension(:), allocatable  :: point_tags
      contains
!--   ---Subroutines-----------------------------------------------------------
      procedure                                  :: Destruct  => MSH_DestructBCStorage
   end type MSH_BCinfo_t

   type :: MSH_point_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: x(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
   end type MSH_point_t

   type :: MSH_curve_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: minX(3)
      real(kind=RP)       :: maxX(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
      integer             :: no_bps
      integer             :: bps(2)
   end type MSH_curve_t

   type :: MSH_surf_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: minX(3)
      real(kind=RP)       :: maxX(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
      integer             :: no_bps
      integer             :: bps(4)
   end type MSH_surf_t

   type :: MSH_vol_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: minX(3)
      real(kind=RP)       :: maxX(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
      integer             :: no_bps
      integer             :: bps(6)
   end type MSH_vol_t

   type :: MSH_node_block_t
!-----Variables-----------------------------------------------------------
      integer             :: entity_dim
      integer             :: entity_tag
      logical             :: parametric
      integer             :: no_nodes
      integer, dimension(:), allocatable         :: tags
      real(kind=RP), dimension(:,:), allocatable :: cords
      contains
!--   ---Subroutines-----------------------------------------------------------
      procedure                                  :: Construct => MSH_ConstructNodeBlock
      procedure                                  :: Destruct  => MSH_DestructNodeBlock
   end type MSH_node_block_t

   type :: MSH_element_block_t
!-----Variables-----------------------------------------------------------
      integer             :: el_type
      integer             :: no_els
      integer, dimension(:,:), allocatable   :: BCs
      integer, dimension(:), allocatable     :: tags
      integer, dimension(:,:), allocatable   :: nodes
      contains
!--   ---Subroutines-----------------------------------------------------------
      procedure                                  :: Construct => MSH_ConstructElementBlock
      procedure                                  :: Destruct  => MSH_DestructElementBlock
   end type MSH_element_block_t
!
!
!     ========
      CONTAINS
!     ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructMesh_FromGMSHFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success )
!  ---------------------------------------------------------
!  Build mesh from GMSH file. 
!  ---------------------------------------------------------
      USE Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                   :: self
      integer                         :: nodes
      character(len=*)                :: fileName
      integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      integer                         :: dir2D
      logical           , intent(out) :: success
!-----Local-Variables-----------------------------------------------------
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      logical                         :: tmpb
      integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2
      character(len=127),    allocatable :: tmps_vec(:)
      type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
      integer                         :: msh_no_BCs
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      type(MSH_point_t), dimension(:), allocatable  :: msh_points
      type(MSH_curve_t), dimension(:), allocatable  :: msh_curves
      type(MSH_surf_t) , dimension(:), allocatable  :: msh_surfaces
      type(MSH_vol_t)  , dimension(:), allocatable  :: msh_volumes
      integer                         :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el

      type(MSH_node_block_t)  , dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t)  , dimension(:), allocatable  :: msh_element_blocks

      character(len=1024) :: msh_entity
      real(kind=RP), allocatable :: msh_entity_vec(:)
      integer, dimension(4)      :: check_eltype

      integer                         :: numberOfElements
      integer                         :: numberOfNodes
      integer                         :: numberOfBoundaryFaces
      integer                         :: numberOfFaces
       
      integer                         :: bFaceOrder, numBFacePoints
      integer                         :: i, j, k, l
      integer                         :: fUnit, fileStat
      integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
      real(kind=RP)                   :: x(NDIM)
      integer                         :: faceFlags(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
      real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
!-----Curved-patches------------------------------------------------------
      real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
      real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!-----Flat-patches--------------------------------------------------------
      real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
      real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
      real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
!  -----------------------------------------------------------------------

!-----Check-if-a-mesh-partition-exists-----------------------------------
      if ( MPI_Process % doMPIAction ) then
         if ( mpi_partition % Constructed ) then
            ! ERROR stop "Read_GMSH :: No support for MPI yet."
            call ConstructMeshPartition_FromGMSHFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success ) 
         else         
            ! ERROR Stop "Read_GMSH :: No MPI partitions constructed." 
            call ConstructSimplestMesh_FromGMSHFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success ) 
         end if
         return
      end if
       
      numberOfBoundaryFaces = 0
      success               = .TRUE.

      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         PRINT *, "Error opening file: ", fileName
         success = .FALSE.
         RETURN 
      END IF
!------------------------------------------------------------------------

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file." 
      read(fUnit,*) tmpd, tmpi, tmpi 
      read(fUnit,*) tmps
!------------------------------------------------------------------------

!-----Read-BC-info-------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined." 
      read(fUnit,*) tmpi

      allocate(tmpi_vec1(tmpi))
      allocate(tmpi_vec2(tmpi))
      allocate(tmps_vec(tmpi))

      msh_no_BCs = 0 ! only surfaces count!
      do i=1, tmpi
         read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
         if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
      end do ! tmpi
      if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

      allocate(msh_bcs(msh_no_BCs))

      j = 0
      do i=1, msh_no_BCs
         if(tmpi_vec1(i) .eq. 2) then
            j = j + 1
            msh_bcs(j)%dim  = tmpi_vec1(i)
            msh_bcs(j)%tag  = tmpi_vec2(i)
            msh_bcs(j)%name = tmps_vec(i)
         end if
      end do ! msh_no_BCs
 
      deallocate(tmpi_vec1)
      deallocate(tmpi_vec2)
      deallocate(tmps_vec)

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected." 
!------------------------------------------------------------------------

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found." 
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      allocate(msh_points(msh_no_points))
      allocate(msh_curves(msh_no_curves))
      allocate(msh_surfaces(msh_no_surfaces))
      allocate(msh_volumes(msh_no_volumes))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:6+msh_points(i)%no_ptags)) 
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_curves(i)%tag      = int(msh_entity_vec(1))
         msh_curves(i)%minX     = msh_entity_vec(2:4)
         msh_curves(i)%maxX     = msh_entity_vec(5:7)
         msh_curves(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_curves(i)%no_ptags .gt. 0) msh_curves(i)%ptags(1:msh_curves(i)%no_ptags) = int(msh_entity_vec(9:9+msh_curves(i)%no_ptags)) 
         msh_curves(i)%no_bps = msh_entity_vec(9+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = msh_entity_vec(10+msh_curves(i)%no_ptags:11+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = abs(msh_curves(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_curves
!------------------------------------------------------------------------

!-----Read-surfaces------------------------------------------------------
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_surfaces(i)%tag      = int(msh_entity_vec(1))
         msh_surfaces(i)%minX     = msh_entity_vec(2:4)
         msh_surfaces(i)%maxX     = msh_entity_vec(5:7)
         msh_surfaces(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_surfaces(i)%no_ptags .gt. 0) msh_surfaces(i)%ptags(1:msh_surfaces(i)%no_ptags) = int(msh_entity_vec(9:9+msh_surfaces(i)%no_ptags)) 
         msh_surfaces(i)%no_bps = msh_entity_vec(9+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = msh_entity_vec(10+msh_surfaces(i)%no_ptags:11+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = abs(msh_surfaces(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_surfaces
!------------------------------------------------------------------------

!-----Read-volumes-------------------------------------------------------
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_volumes(i)%tag      = int(msh_entity_vec(1))
         msh_volumes(i)%minX     = msh_entity_vec(2:4)
         msh_volumes(i)%maxX     = msh_entity_vec(5:7)
         msh_volumes(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_volumes(i)%no_ptags .gt. 0) msh_volumes(i)%ptags(1:msh_volumes(i)%no_ptags) = int(msh_entity_vec(9:9+msh_volumes(i)%no_ptags)) 
         msh_volumes(i)%no_bps = msh_entity_vec(9+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps = msh_entity_vec(10+msh_volumes(i)%no_ptags:11+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps = abs(msh_volumes(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_volumes
!------------------------------------------------------------------------

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected." 
!------------------------------------------------------------------------

!-----Read-nodes---------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found." 

      read(fUnit,*) msh_no_nodeblocks, numberOfNodes, tmpi1, tmpi2
      if (numberOfNodes .ne. tmpi2) error stop "READ_gmsh :: Incoherent node numbering."

      ! allocate nodes storage
      allocate(msh_node_blocks(msh_no_nodeblocks))

      msh_global_node = 0
      do msh_nodeblock=1, msh_no_nodeblocks
         read(fUnit,*) tmpi1, tmpi2, tmpi3, msh_nodes_per_block
         tmpb=tmpi3
         if (tmpb) error stop "READ_gmsh :: Parametric nodes not supported."

         call msh_node_blocks(msh_nodeblock) % Construct(tmpi1, tmpi2, tmpb, msh_nodes_per_block)

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % tags(msh_node)
         end do ! msh_nodes_per_block /for tags

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % cords(msh_node,:)
         end do ! msh_nodes_per_block /for coordinates

         msh_global_node = msh_global_node + msh_nodes_per_block

      end do ! msh_no_nodeblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected." 
!------------------------------------------------------------------------

!-----Read-elements------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found." 

      read(fUnit,*) msh_no_elblocks, numberOfElements, tmpi1, tmpi2
      if (numberOfElements .ne. tmpi2) error stop "READ_gmsh :: Incoherent element numbering."

      allocate(msh_element_blocks(msh_no_elblocks))

      msh_global_el = 0
      do msh_elblock=1, msh_no_elblocks

         read(fUnit,*) tmpi1, tmpi2, element_type, msh_els_per_block

         call msh_element_blocks(msh_elblock) % Construct(element_type,msh_els_per_block)

         do msh_el=1, msh_els_per_block
            read(fUnit,*) msh_element_blocks(msh_elblock) % tags(msh_el), msh_element_blocks(msh_elblock) % nodes(msh_el,:)
         end do ! msh_els_per_block /for tags

         msh_global_el = msh_global_el + msh_els_per_block

      end do ! msh_no_elblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected." 
      close( fUnit )
!------------------------------------------------------------------------

!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype(1) = count(msh_element_blocks(:) % el_type .eq.  5)
      check_eltype(2) = count(msh_element_blocks(:) % el_type .eq. 12)
      check_eltype(3) = count(msh_element_blocks(:) % el_type .eq. 92)
      check_eltype(4) = count(msh_element_blocks(:) % el_type .eq. 93)
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      check_eltype = (/5,12,92,93/)
      org_element_type = check_eltype(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do
!------------------------------------------------------------------------

!-----Reorder-nodes-in-elements------------------------------------------
      select case(org_element_type)
      case (5)
         allocate(tmpi_vec1(8))
      case (12)
         error stop "READ_GMSH :: TBC"
      case (92)
         error stop "READ_GMSH :: TBC"
      case (93)
         error stop "READ_GMSH :: TBC"
      case default 
         error stop "READ_GMSH :: No 3D elements detected in the mesh."
      end select
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els 
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               msh_element_blocks(msh_elblock) % nodes(j,1) = tmpi_vec1(4)
               msh_element_blocks(msh_elblock) % nodes(j,2:4) = tmpi_vec1(1:3)
               msh_element_blocks(msh_elblock) % nodes(j,5) = tmpi_vec1(8)
               msh_element_blocks(msh_elblock) % nodes(j,6:8) = tmpi_vec1(5:7)
            end do
         end if
      end do
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-BCs---------------------------------------------------------
      do i=1,msh_no_BCs
         ! find no surfaces
         msh_bcs(i) % no_of_surfaces = 0
         do j=1, msh_no_surfaces
            tmpi = findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
         end do ! msh_no_surfaces

         ! assign surfaces and curves
         allocate(tmpi_vec1(4 * msh_bcs(i) % no_of_surfaces))
         allocate(msh_bcs(i) % surface_tags(msh_bcs(i) % no_of_surfaces))
         k = 0
         do j=1, msh_no_surfaces
            tmpi = findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               msh_bcs(i) % surface_tags(k) = msh_surfaces(j) % tag
               tmpi_vec1(1 + 4*(k-1) : 4*k) = msh_surfaces(j) % bps 
            end if ! tmpi
         end do ! msh_no_surfaces
         call unique(tmpi_vec1,msh_bcs(i) % curve_tags)
         msh_bcs(i) % no_of_curves = size(msh_bcs(i) % curve_tags,1)
         deallocate(tmpi_vec1)

         ! assign points
         k = 0
         allocate(tmpi_vec1(2 * msh_bcs(i) % no_of_curves))
         do j=1, msh_no_curves
            tmpi = findloc(msh_bcs(i) % curve_tags, msh_curves(j)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               tmpi_vec1(1 + (k-1)*2: k*2) = msh_curves(j) % bps
            end if ! tmpi
         end do ! msh_no_points
         call unique(tmpi_vec1,msh_bcs(i) % point_tags)
         msh_bcs(i) % no_of_points = size(msh_bcs(i) % point_tags,1)
         deallocate(tmpi_vec1)

         ! assign nodes
         k = 0
         allocate(tmpi_vec1(numberOfNodes)) ! approximate upper bound
         tmpi_vec1 = 0
         do j=1, msh_no_nodeblocks

            select case (msh_node_blocks(j) % entity_dim)
            case (0) ! points 
               tmpi = findloc(msh_bcs(i) % point_tags, msh_node_blocks(j) % entity_tag, 1)
            case (1) ! curves 
               tmpi = findloc(msh_bcs(i) % curve_tags, msh_node_blocks(j) % entity_tag, 1)
            case (2) ! surfaces
               tmpi = findloc(msh_bcs(i) % surface_tags, msh_node_blocks(j) % entity_tag, 1)
            case (3) ! volume
               tmpi = 0
            end select ! msh_node_blocks(j) % entity_dim

            if (tmpi .gt. 0) then 
               tmpi_vec1(k+1 : k+msh_node_blocks(j) % no_nodes) = msh_node_blocks(j) % tags
               k = k + msh_node_blocks(j) % no_nodes
            end if ! tmpi

         end do ! msh_no_nodeblocks
         call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
         msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
         deallocate(tmpi_vec1)

      end do ! msh_no_BCs
!------------------------------------------------------------------------

!-----Assign-BC-to-nodes-------------------------------------------------
      allocate(tmpi_vec1(8)) ! help array for the walls
      do i=1, msh_no_BCs
         do msh_elblock=1, msh_no_elblocks
            if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
               do j=1, msh_element_blocks(msh_elblock) % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = findloc(msh_bcs(i) % node_tags,msh_element_blocks(msh_elblock) % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes 
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,4) = msh_bcs(i) % tag ! 4:front 
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_element_blocks(msh_elblock) % no_els 
            end if
         end do ! msh_no_elblocks
      end do ! msh_no_BCs
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-new-tags----------------------------------------------------
      tmp_eltag = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els 
               tmp_eltag = tmp_eltag + 1
               msh_element_blocks(msh_elblock) % tags(j) = tmp_eltag
            end do
         end if
      end do
      if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!------------------------------------------------------------------------

!-----Build-nodes--------------------------------------------------------
      self % nodeType = nodes
      self % no_of_elements = numberOfElements
      self % no_of_allElements = numberOfElements
!------------------------------------------------------------------------

!-----Set-up-face-patches------------------------------------------------
      numBFacePoints = bFaceOrder + 1
      allocate(uNodes(numBFacePoints))
      allocate(vNodes(numBFacePoints))
      allocate(values(3,numBFacePoints,numBFacePoints))
      do i = 1, numBFacePoints
         uNodes(i) = -cos((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         vNodes(i) = -cos((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
      end do
!------------------------------------------------------------------------

!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(numberOfelements) )
      allocate( self % nodes(numberOfNodes) )
      allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
!------------------------------------------------------------------------

!----Set-nodes-----------------------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref
            call ConstructNode( self % nodes(msh_node_blocks(msh_nodeblock) % tags(msh_node)), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-elements-----------------------------------------------------------
      j = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do msh_el = 1, msh_element_blocks(msh_elblock) % no_els 
               j = j + 1
               ! setting l'th element
               l =  msh_element_blocks(msh_elblock) % tags(msh_el)
               nodeIDs = msh_element_blocks(msh_elblock) % nodes(msh_el,1:8) ! only non-curved nodes

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(l) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(l) % SurfInfo % corners = corners
               else ! curved mesh

               end if

               call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = findloc( msh_bcs % tag,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(l) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(l) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(l) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(l) % boundaryName(k) )
                     numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(l) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(l) % boundaryName(k)), trim(self % elements(l) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k
            end do ! msh_element_blocks(msh_elblock) % no_els
         end if ! if el_type .eq. org_element_type
      end do ! msh_no_elblocks
      if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."
!------------------------------------------------------------------------

!-----Deallocate-msh-vairables-------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         call msh_node_blocks(msh_nodeblock) % Destruct()
      end do
      deallocate(msh_node_blocks)

      do msh_elblock=1, msh_no_elblocks
         call msh_element_blocks(msh_elblock) % Destruct()
      end do
      deallocate(msh_element_blocks)

      do i=1,msh_no_BCs
         call msh_bcs(i) % Destruct()
      end do
      deallocate(msh_bcs)

      deallocate(msh_points)
      deallocate(msh_curves)
      deallocate(msh_surfaces)
      deallocate(msh_volumes)
      deallocate(msh_entity_vec)

      deallocate(uNodes) ! Check if we can do it! FIXME
      deallocate(vNodes) ! Check if we can do it! FIXME
      deallocate(values) ! Check if we can do it! FIXME
!------------------------------------------------------------------------
!
!     ---------------------------
!     Construct the element faces
!     ---------------------------
!
      numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
      self % numberOfFaces = numberOfFaces
      allocate( self % faces(self % numberOfFaces) )
      CALL ConstructFaces( self, success )
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % ConstructZones()
!
!     ---------------------------
!     Construct periodic faces
!     ---------------------------
!
      CALL ConstructPeriodicFaces( self ) 
!
!     ---------------------------
!     Delete periodic- faces
!     ---------------------------
!
      CALL DeletePeriodicMinusFaces( self ) 
!
!     ---------------------------
!     Assign faces ID to elements
!     ---------------------------
!
      CALL getElementsFaceIDs(self)
!
!     ---------------------
!     Define boundary faces
!     ---------------------
!
      call self % DefineAsBoundaryFaces()
!
!     -----------------------------------
!     Check if this is a 2D extruded mesh
!     -----------------------------------
!
      call self % CheckIfMeshIs2D()
!
!     -------------------------------
!     Set the mesh as 2D if requested
!     -------------------------------
!
      self % dir2D_ctrl = dir2D
      if ( dir2D .ne. 0 ) then
         call SetMappingsToCrossProduct
         call self % CorrectOrderFor2DMesh(dir2D,1)
      end if
!
!     ------------------------------
!     Set the element connectivities
!     ------------------------------
!
      call self % SetConnectivitiesAndLinkFaces(nodes)
!
!     ---------------------------------------
!     Construct elements' and faces' geometry
!     ---------------------------------------
!
      call self % ConstructGeometry()
!
!     -------------------------------
!     Set the mesh as 2D if requested
!     -------------------------------
!
      if ( dir2D .ne. 0 ) then
         call self % CorrectOrderFor2DMesh(dir2D,0)
      end if
!
!     ---------------------------------
!     Describe mesh and prepare for I/O
!     ---------------------------------
!
      if (.not. self % child) CALL self % Describe( trim(fileName), bFaceOrder )
      call self % PrepareForIO
      
      call self % ExportBoundaryMesh (trim(fileName))
         
   end subroutine ConstructMesh_FromGMSHFile_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------
!     Constructor of mesh partitions
!     ------------------------------
   SUBROUTINE ConstructMeshPartition_FromGMSHFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success )
!  ---------------------------------------------------------
!  Build mesh from GMSH file. 
!  ---------------------------------------------------------
      USE Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      use MPI_Face_Class
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                   :: self
      integer                         :: nodes
      character(len=*)                :: fileName
      integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      integer                         :: dir2D
      logical           , intent(out) :: success
!-----Local-Variables-----------------------------------------------------
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      logical                         :: tmpb
      integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2
      character(len=127),    allocatable :: tmps_vec(:)
      type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
      integer                         :: msh_no_BCs
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      type(MSH_point_t), dimension(:), allocatable  :: msh_points
      type(MSH_curve_t), dimension(:), allocatable  :: msh_curves
      type(MSH_surf_t) , dimension(:), allocatable  :: msh_surfaces
      type(MSH_vol_t)  , dimension(:), allocatable  :: msh_volumes
      integer                         :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el

      type(MSH_node_block_t)  , dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t)  , dimension(:), allocatable  :: msh_element_blocks

      character(len=1024) :: msh_entity
      real(kind=RP), allocatable :: msh_entity_vec(:)
      integer, dimension(4)      :: check_eltype

      integer                         :: numberOfAllElements
      integer                         :: numberOfAllNodes
      integer, allocatable            :: globalToLocalNodeID(:)
      integer, allocatable            :: globalToLocalElementID(:)

      integer                         :: numberOfElements
      integer                         :: numberOfNodes
      integer                         :: numberOfBoundaryFaces
      integer                         :: numberOfFaces
       
      integer                         :: bFaceOrder, numBFacePoints
      integer                         :: i, j, k, l, pNode, pElement
      integer                         :: fUnit, fileStat
      integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
      real(kind=RP)                   :: x(NDIM)
      integer                         :: faceFlags(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
      real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
!-----Curved-patches------------------------------------------------------
      real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
      real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!-----Flat-patches--------------------------------------------------------
      real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
      real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
      real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
!  -----------------------------------------------------------------------

      success               = .TRUE.

      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         PRINT *, "Error opening file: ", fileName
         success = .FALSE.
         RETURN 
      END IF
!------------------------------------------------------------------------

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file." 
      read(fUnit,*) tmpd, tmpi, tmpi 
      read(fUnit,*) tmps
!------------------------------------------------------------------------

!-----Read-BC-info-------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined." 
      read(fUnit,*) tmpi

      allocate(tmpi_vec1(tmpi))
      allocate(tmpi_vec2(tmpi))
      allocate(tmps_vec(tmpi))

      msh_no_BCs = 0 ! only surfaces count!
      do i=1, tmpi
         read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
         if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
      end do ! tmpi
      if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

      allocate(msh_bcs(msh_no_BCs))

      j = 0
      do i=1, msh_no_BCs
         if(tmpi_vec1(i) .eq. 2) then
            j = j + 1
            msh_bcs(j)%dim  = tmpi_vec1(i)
            msh_bcs(j)%tag  = tmpi_vec2(i)
            msh_bcs(j)%name = tmps_vec(i)
         end if
      end do ! msh_no_BCs
 
      deallocate(tmpi_vec1)
      deallocate(tmpi_vec2)
      deallocate(tmps_vec)

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected." 
!------------------------------------------------------------------------

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found." 
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      allocate(msh_points(msh_no_points))
      allocate(msh_curves(msh_no_curves))
      allocate(msh_surfaces(msh_no_surfaces))
      allocate(msh_volumes(msh_no_volumes))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:6+msh_points(i)%no_ptags)) 
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_curves(i)%tag      = int(msh_entity_vec(1))
         msh_curves(i)%minX     = msh_entity_vec(2:4)
         msh_curves(i)%maxX     = msh_entity_vec(5:7)
         msh_curves(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_curves(i)%no_ptags .gt. 0) msh_curves(i)%ptags(1:msh_curves(i)%no_ptags) = int(msh_entity_vec(9:9+msh_curves(i)%no_ptags)) 
         msh_curves(i)%no_bps = msh_entity_vec(9+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = msh_entity_vec(10+msh_curves(i)%no_ptags:11+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = abs(msh_curves(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_curves
!------------------------------------------------------------------------

!-----Read-surfaces------------------------------------------------------
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_surfaces(i)%tag      = int(msh_entity_vec(1))
         msh_surfaces(i)%minX     = msh_entity_vec(2:4)
         msh_surfaces(i)%maxX     = msh_entity_vec(5:7)
         msh_surfaces(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_surfaces(i)%no_ptags .gt. 0) msh_surfaces(i)%ptags(1:msh_surfaces(i)%no_ptags) = int(msh_entity_vec(9:9+msh_surfaces(i)%no_ptags)) 
         msh_surfaces(i)%no_bps = msh_entity_vec(9+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = msh_entity_vec(10+msh_surfaces(i)%no_ptags:11+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = abs(msh_surfaces(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_surfaces
!------------------------------------------------------------------------

!-----Read-volumes-------------------------------------------------------
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_volumes(i)%tag      = int(msh_entity_vec(1))
         msh_volumes(i)%minX     = msh_entity_vec(2:4)
         msh_volumes(i)%maxX     = msh_entity_vec(5:7)
         msh_volumes(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_volumes(i)%no_ptags .gt. 0) msh_volumes(i)%ptags(1:msh_volumes(i)%no_ptags) = int(msh_entity_vec(9:9+msh_volumes(i)%no_ptags)) 
         msh_volumes(i)%no_bps = msh_entity_vec(9+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps = msh_entity_vec(10+msh_volumes(i)%no_ptags:11+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps = abs(msh_volumes(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_volumes
!------------------------------------------------------------------------

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected." 
!------------------------------------------------------------------------

!-----Read-nodes---------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found." 

      read(fUnit,*) msh_no_nodeblocks, numberOfNodes, tmpi1, tmpi2
      if (numberOfNodes .ne. tmpi2) error stop "READ_gmsh :: Incoherent node numbering."

      ! allocate nodes storage
      allocate(msh_node_blocks(msh_no_nodeblocks))

      msh_global_node = 0
      do msh_nodeblock=1, msh_no_nodeblocks
         read(fUnit,*) tmpi1, tmpi2, tmpi3, msh_nodes_per_block
         tmpb=tmpi3
         if (tmpb) error stop "READ_gmsh :: Parametric nodes not supported."

         call msh_node_blocks(msh_nodeblock) % Construct(tmpi1, tmpi2, tmpb, msh_nodes_per_block)

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % tags(msh_node)
         end do ! msh_nodes_per_block /for tags

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % cords(msh_node,:)
         end do ! msh_nodes_per_block /for coordinates

         msh_global_node = msh_global_node + msh_nodes_per_block

      end do ! msh_no_nodeblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected." 
      numberOfAllNodes = numberOfNodes
!------------------------------------------------------------------------

!-----Read-elements------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found." 

      read(fUnit,*) msh_no_elblocks, numberOfElements, tmpi1, tmpi2
      if (numberOfElements .ne. tmpi2) error stop "READ_gmsh :: Incoherent element numbering."

      allocate(msh_element_blocks(msh_no_elblocks))

      msh_global_el = 0
      do msh_elblock=1, msh_no_elblocks

         read(fUnit,*) tmpi1, tmpi2, element_type, msh_els_per_block

         call msh_element_blocks(msh_elblock) % Construct(element_type,msh_els_per_block)

         do msh_el=1, msh_els_per_block
            read(fUnit,*) msh_element_blocks(msh_elblock) % tags(msh_el), msh_element_blocks(msh_elblock) % nodes(msh_el,:)
         end do ! msh_els_per_block /for tags

         msh_global_el = msh_global_el + msh_els_per_block

      end do ! msh_no_elblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected." 
      close( fUnit )
!------------------------------------------------------------------------

!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype(1) = count(msh_element_blocks(:) % el_type .eq.  5)
      check_eltype(2) = count(msh_element_blocks(:) % el_type .eq. 12)
      check_eltype(3) = count(msh_element_blocks(:) % el_type .eq. 92)
      check_eltype(4) = count(msh_element_blocks(:) % el_type .eq. 93)
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      check_eltype = (/5,12,92,93/)
      org_element_type = check_eltype(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do
      numberOfAllElements = numberOfElements
!------------------------------------------------------------------------

!-----Reorder-nodes-in-elements------------------------------------------
      select case(org_element_type)
      case (5)
         allocate(tmpi_vec1(8))
      case (12)
         error stop "READ_GMSH :: TBC"
      case (92)
         error stop "READ_GMSH :: TBC"
      case (93)
         error stop "READ_GMSH :: TBC"
      case default 
         error stop "READ_GMSH :: No 3D elements detected in the mesh."
      end select
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els 
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               msh_element_blocks(msh_elblock) % nodes(j,1) = tmpi_vec1(4)
               msh_element_blocks(msh_elblock) % nodes(j,2:4) = tmpi_vec1(1:3)
               msh_element_blocks(msh_elblock) % nodes(j,5) = tmpi_vec1(8)
               msh_element_blocks(msh_elblock) % nodes(j,6:8) = tmpi_vec1(5:7)
            end do
         end if
      end do
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-BCs---------------------------------------------------------
      do i=1,msh_no_BCs
         ! find no surfaces
         msh_bcs(i) % no_of_surfaces = 0
         do j=1, msh_no_surfaces
            tmpi = findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
         end do ! msh_no_surfaces

         ! assign surfaces and curves
         allocate(tmpi_vec1(4 * msh_bcs(i) % no_of_surfaces))
         allocate(msh_bcs(i) % surface_tags(msh_bcs(i) % no_of_surfaces))
         k = 0
         do j=1, msh_no_surfaces
            tmpi = findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               msh_bcs(i) % surface_tags(k) = msh_surfaces(j) % tag
               tmpi_vec1(1 + 4*(k-1) : 4*k) = msh_surfaces(j) % bps 
            end if ! tmpi
         end do ! msh_no_surfaces
         call unique(tmpi_vec1,msh_bcs(i) % curve_tags)
         msh_bcs(i) % no_of_curves = size(msh_bcs(i) % curve_tags,1)
         deallocate(tmpi_vec1)

         ! assign points
         k = 0
         allocate(tmpi_vec1(2 * msh_bcs(i) % no_of_curves))
         do j=1, msh_no_curves
            tmpi = findloc(msh_bcs(i) % curve_tags, msh_curves(j)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               tmpi_vec1(1 + (k-1)*2: k*2) = msh_curves(j) % bps
            end if ! tmpi
         end do ! msh_no_points
         call unique(tmpi_vec1,msh_bcs(i) % point_tags)
         msh_bcs(i) % no_of_points = size(msh_bcs(i) % point_tags,1)
         deallocate(tmpi_vec1)

         ! assign nodes
         k = 0
         allocate(tmpi_vec1(numberOfNodes)) ! approximate upper bound
         tmpi_vec1 = 0
         do j=1, msh_no_nodeblocks

            select case (msh_node_blocks(j) % entity_dim)
            case (0) ! points 
               tmpi = findloc(msh_bcs(i) % point_tags, msh_node_blocks(j) % entity_tag, 1)
            case (1) ! curves 
               tmpi = findloc(msh_bcs(i) % curve_tags, msh_node_blocks(j) % entity_tag, 1)
            case (2) ! surfaces
               tmpi = findloc(msh_bcs(i) % surface_tags, msh_node_blocks(j) % entity_tag, 1)
            case (3) ! volume
               tmpi = 0
            end select ! msh_node_blocks(j) % entity_dim

            if (tmpi .gt. 0) then 
               tmpi_vec1(k+1 : k+msh_node_blocks(j) % no_nodes) = msh_node_blocks(j) % tags
               k = k + msh_node_blocks(j) % no_nodes
            end if ! tmpi

         end do ! msh_no_nodeblocks
         call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
         msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
         deallocate(tmpi_vec1)

      end do ! msh_no_BCs
!------------------------------------------------------------------------

!-----Assign-BC-to-nodes-------------------------------------------------
      allocate(tmpi_vec1(8)) ! help array for the walls
      do i=1, msh_no_BCs
         do msh_elblock=1, msh_no_elblocks
            if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
               do j=1, msh_element_blocks(msh_elblock) % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = findloc(msh_bcs(i) % node_tags,msh_element_blocks(msh_elblock) % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes 
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,4) = msh_bcs(i) % tag ! 4:front 
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_element_blocks(msh_elblock) % no_els 
            end if
         end do ! msh_no_elblocks
      end do ! msh_no_BCs
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-new-tags----------------------------------------------------
      tmp_eltag = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els 
               tmp_eltag = tmp_eltag + 1
               msh_element_blocks(msh_elblock) % tags(j) = tmp_eltag
            end do
         end if
      end do
      if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!------------------------------------------------------------------------

!-----Build-nodes--------------------------------------------------------
      self % nodeType = nodes
      self % no_of_elements = mpi_partition % no_of_elements
      self % no_of_allElements = numberOfAllElements
      ! print *, "MPI_Process % rank: ", MPI_Process % rank
      ! print *, "There are ", self % no_of_elements , " elements out of ", numberOfAllElements
      ! print *, "There are ", mpi_partition % no_of_nodes , " nodes out of ", numberOfAllNodes
!------------------------------------------------------------------------

!-----Set-up-face-patches------------------------------------------------
      numBFacePoints = bFaceOrder + 1
      allocate(uNodes(numBFacePoints))
      allocate(vNodes(numBFacePoints))
      allocate(values(3,numBFacePoints,numBFacePoints))
      do i = 1, numBFacePoints
         uNodes(i) = -cos((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         vNodes(i) = -cos((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
      end do
!------------------------------------------------------------------------

!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(mpi_partition % no_of_elements) )
      allocate( self % nodes(mpi_partition % no_of_nodes) )
      allocate ( self % Nx(self % no_of_elements) , self % Ny(self % no_of_elements) , self % Nz(self % no_of_elements) )
      allocate( globalToLocalNodeID(numberOfAllNodes) )
      allocate( globalToLocalElementID(numberOfAllElements) )
      
      globalToLocalNodeID = -1
      globalToLocalElementID = -1
!------------------------------------------------------------------------

!----Set-nodes-----------------------------------------------------------
      pNode = 1
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref

            if ( pNode .gt. mpi_partition % no_of_nodes ) cycle

            ! Construct only nodes that belong to the partition
            if ( msh_node_blocks(msh_nodeblock) % tags(msh_node) .eq. mpi_partition % nodeIDs(pNode) ) then
               call ConstructNode( self % nodes(pNode), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
               globalToLocalNodeID(msh_node_blocks(msh_nodeblock) % tags(msh_node)) = pNode
               pNode = pNode + 1 
            end if
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-elements-----------------------------------------------------------
      j = 0
      pElement = 1
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do msh_el = 1, msh_element_blocks(msh_elblock) % no_els 
               j = j + 1
               ! setting l'th element
               l =  msh_element_blocks(msh_elblock) % tags(msh_el)
               nodeIDs = msh_element_blocks(msh_elblock) % nodes(msh_el,1:8) ! only non-curved nodes
               nodeIDs = globalToLocalNodeID(nodeIDs)

               if ( pElement .gt. mpi_partition % no_of_elements ) then
                  cycle
               else if ( l .ne. mpi_partition % elementIDs(pElement) ) then
                  cycle
               end if

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(pElement) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(pElement) % SurfInfo % corners = corners
               else ! curved mesh

               end if

               call self % elements(pElement) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , pElement, l)

               self % Nx(pElement) = Nx(l)
               self % Ny(pElement) = Ny(l)
               self % Nz(pElement) = Nz(l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = findloc( msh_bcs % tag,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(pElement) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(pElement) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(pElement) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(pElement) % boundaryName(k) )
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(pElement) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(pElement) % boundaryName(k)), trim(self % elements(pElement) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k

               globalToLocalElementID(l) = pElement
               pElement = pElement + 1
            end do ! msh_element_blocks(msh_elblock) % no_els
         end if ! if el_type .eq. org_element_type
      end do ! msh_no_elblocks
!------------------------------------------------------------------------

!-----Deallocate-msh-vairables-------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         call msh_node_blocks(msh_nodeblock) % Destruct()
      end do
      deallocate(msh_node_blocks)

      do msh_elblock=1, msh_no_elblocks
         call msh_element_blocks(msh_elblock) % Destruct()
      end do
      deallocate(msh_element_blocks)

      do i=1,msh_no_BCs
         call msh_bcs(i) % Destruct()
      end do
      deallocate(msh_bcs)

      deallocate(msh_points)
      deallocate(msh_curves)
      deallocate(msh_surfaces)
      deallocate(msh_volumes)
      deallocate(msh_entity_vec)

      deallocate(uNodes) ! Check if we can do it! FIXME
      deallocate(vNodes) ! Check if we can do it! FIXME
      deallocate(values) ! Check if we can do it! FIXME
!------------------------------------------------------------------------
!
!     ---------------------------
!     Construct the element faces
!     ---------------------------
!
      numberOfFaces        = GetOriginalNumberOfFaces(self)
      self % numberOfFaces = numberOfFaces
      allocate( self % faces(self % numberOfFaces) )
      CALL ConstructFaces( self, success )
!        --------------------------------
!        Get actual mesh element face IDs
!        --------------------------------
!
      CALL getElementsFaceIDs(self)
!
!     --------------
!     Cast MPI faces
!     --------------
!
      call ConstructMPIFaces( self % MPIfaces )
      call self % UpdateFacesWithPartition(mpi_partition, &
                                           numberOfAllElements, &
                                           globalToLocalElementID)
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % ConstructZones()
!
!     ---------------------------
!     Construct periodic faces
!     ---------------------------
!
      CALL ConstructPeriodicFaces( self )
!
!     ---------------------------
!     Delete periodic- faces
!     ---------------------------
!
      CALL DeletePeriodicMinusFaces( self )
!
!     ---------------------------
!     Assign faces ID to elements
!     ---------------------------
!
      CALL getElementsFaceIDs(self)
!
!     ---------------------
!     Define boundary faces
!     ---------------------
!
      call self % DefineAsBoundaryFaces()
!
!     -----------------------------------
!     Check if this is a 2D extruded mesh
!     -----------------------------------
!
      call self % CheckIfMeshIs2D()
!
!     -------------------------------
!     Set the mesh as 2D if requested
!     -------------------------------
!
      if ( dir2D .ne. 0 ) then
         call SetMappingsToCrossProduct
         call self % CorrectOrderFor2DMesh(dir2D,1)
      end if
!
!     ------------------------------
!     Set the element connectivities
!     ------------------------------
!
      call self % SetConnectivitiesAndLinkFaces(nodes)
!
!     ---------------------------------------
!     Construct elements' and faces' geometry
!     ---------------------------------------
!
      call self % ConstructGeometry()

      if ( dir2D .ne. 0 ) then
         call self % CorrectOrderFor2DMesh(dir2D,0)
      end if
!
!     ---------
!     Finish up
!     ---------
!
      if (.not. self % child) then
         CALL self % Describe         ( trim(fileName), bFaceOrder )
         CALL self % DescribePartition( )
      end if
!
!     --------------------
!     Prepare mesh for I/O
!     --------------------
!
      call self % PrepareForIO

      deallocate(globalToLocalNodeID)
      deallocate(globalToLocalElementID)

   END SUBROUTINE ConstructMeshPartition_FromGMSHFile_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   subroutine ConstructSimplestMesh_FromGMSHFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success )
!  ---------------------------------------------------------
!  Build mesh from GMSH file. 
!  ---------------------------------------------------------
      USE Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                   :: self
      integer                         :: nodes
      character(len=*)                :: fileName
      integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      integer                         :: dir2D
      logical           , intent(out) :: success
!-----Local-Variables-----------------------------------------------------
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      logical                         :: tmpb
      integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2
      character(len=127),    allocatable :: tmps_vec(:)
      type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
      integer                         :: msh_no_BCs
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      type(MSH_point_t), dimension(:), allocatable  :: msh_points
      type(MSH_curve_t), dimension(:), allocatable  :: msh_curves
      type(MSH_surf_t) , dimension(:), allocatable  :: msh_surfaces
      type(MSH_vol_t)  , dimension(:), allocatable  :: msh_volumes
      integer                         :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el

      type(MSH_node_block_t)  , dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t)  , dimension(:), allocatable  :: msh_element_blocks

      character(len=1024) :: msh_entity
      real(kind=RP), allocatable :: msh_entity_vec(:)
      integer, dimension(4)      :: check_eltype

      integer                         :: numberOfElements
      integer                         :: numberOfNodes
      integer                         :: numberOfBoundaryFaces
      integer                         :: numberOfFaces
       
      integer                         :: bFaceOrder, numBFacePoints
      integer                         :: i, j, k, l
      integer                         :: fUnit, fileStat
      integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
      real(kind=RP)                   :: x(NDIM)
      integer                         :: faceFlags(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
      real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
!-----Curved-patches------------------------------------------------------
      real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
      real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!-----Flat-patches--------------------------------------------------------
      real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
      real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
      real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
!  -----------------------------------------------------------------------

      numberOfBoundaryFaces = 0
      success               = .TRUE.

      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         PRINT *, "Error opening file: ", fileName
         success = .FALSE.
         RETURN 
      END IF

!------------------------------------------------------------------------

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file." 
      read(fUnit,*) tmpd, tmpi, tmpi 
      read(fUnit,*) tmps
!------------------------------------------------------------------------

!-----Read-BC-info-------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined." 
      read(fUnit,*) tmpi

      allocate(tmpi_vec1(tmpi))
      allocate(tmpi_vec2(tmpi))
      allocate(tmps_vec(tmpi))

      msh_no_BCs = 0 ! only surfaces count!
      do i=1, tmpi
         read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
         if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
      end do ! tmpi
      if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

      allocate(msh_bcs(msh_no_BCs))

      j = 0
      do i=1, msh_no_BCs
         if(tmpi_vec1(i) .eq. 2) then
            j = j + 1
            msh_bcs(j)%dim  = tmpi_vec1(i)
            msh_bcs(j)%tag  = tmpi_vec2(i)
            msh_bcs(j)%name = tmps_vec(i)
         end if
      end do ! msh_no_BCs
 
      deallocate(tmpi_vec1)
      deallocate(tmpi_vec2)
      deallocate(tmps_vec)

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected." 
!------------------------------------------------------------------------

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found." 
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      allocate(msh_points(msh_no_points))
      allocate(msh_curves(msh_no_curves))
      allocate(msh_surfaces(msh_no_surfaces))
      allocate(msh_volumes(msh_no_volumes))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:6+msh_points(i)%no_ptags)) 
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_curves(i)%tag      = int(msh_entity_vec(1))
         msh_curves(i)%minX     = msh_entity_vec(2:4)
         msh_curves(i)%maxX     = msh_entity_vec(5:7)
         msh_curves(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_curves(i)%no_ptags .gt. 0) msh_curves(i)%ptags(1:msh_curves(i)%no_ptags) = int(msh_entity_vec(9:9+msh_curves(i)%no_ptags)) 
         msh_curves(i)%no_bps = msh_entity_vec(9+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = msh_entity_vec(10+msh_curves(i)%no_ptags:11+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = abs(msh_curves(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_curves
!------------------------------------------------------------------------

!-----Read-surfaces------------------------------------------------------
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_surfaces(i)%tag      = int(msh_entity_vec(1))
         msh_surfaces(i)%minX     = msh_entity_vec(2:4)
         msh_surfaces(i)%maxX     = msh_entity_vec(5:7)
         msh_surfaces(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_surfaces(i)%no_ptags .gt. 0) msh_surfaces(i)%ptags(1:msh_surfaces(i)%no_ptags) = int(msh_entity_vec(9:9+msh_surfaces(i)%no_ptags)) 
         msh_surfaces(i)%no_bps = msh_entity_vec(9+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = msh_entity_vec(10+msh_surfaces(i)%no_ptags:11+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = abs(msh_surfaces(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_surfaces
!------------------------------------------------------------------------

!-----Read-volumes-------------------------------------------------------
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(1024a)') msh_entity 
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_volumes(i)%tag      = int(msh_entity_vec(1))
         msh_volumes(i)%minX     = msh_entity_vec(2:4)
         msh_volumes(i)%maxX     = msh_entity_vec(5:7)
         msh_volumes(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_volumes(i)%no_ptags .gt. 0) msh_volumes(i)%ptags(1:msh_volumes(i)%no_ptags) = int(msh_entity_vec(9:9+msh_volumes(i)%no_ptags)) 
         msh_volumes(i)%no_bps = msh_entity_vec(9+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps = msh_entity_vec(10+msh_volumes(i)%no_ptags:11+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps = abs(msh_volumes(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_volumes
!------------------------------------------------------------------------

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected." 
!------------------------------------------------------------------------

!-----Read-nodes---------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found." 

      read(fUnit,*) msh_no_nodeblocks, numberOfNodes, tmpi1, tmpi2
      if (numberOfNodes .ne. tmpi2) error stop "READ_gmsh :: Incoherent node numbering."

      ! allocate nodes storage
      allocate(msh_node_blocks(msh_no_nodeblocks))

      msh_global_node = 0
      do msh_nodeblock=1, msh_no_nodeblocks
         read(fUnit,*) tmpi1, tmpi2, tmpi3, msh_nodes_per_block
         tmpb=tmpi3
         if (tmpb) error stop "READ_gmsh :: Parametric nodes not supported."

         call msh_node_blocks(msh_nodeblock) % Construct(tmpi1, tmpi2, tmpb, msh_nodes_per_block)

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % tags(msh_node)
         end do ! msh_nodes_per_block /for tags

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % cords(msh_node,:)
         end do ! msh_nodes_per_block /for coordinates

         msh_global_node = msh_global_node + msh_nodes_per_block

      end do ! msh_no_nodeblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected." 
!------------------------------------------------------------------------

!-----Read-elements------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found." 

      read(fUnit,*) msh_no_elblocks, numberOfElements, tmpi1, tmpi2
      if (numberOfElements .ne. tmpi2) error stop "READ_gmsh :: Incoherent element numbering."

      allocate(msh_element_blocks(msh_no_elblocks))

      msh_global_el = 0
      do msh_elblock=1, msh_no_elblocks

         read(fUnit,*) tmpi1, tmpi2, element_type, msh_els_per_block

         call msh_element_blocks(msh_elblock) % Construct(element_type,msh_els_per_block)

         do msh_el=1, msh_els_per_block
            read(fUnit,*) msh_element_blocks(msh_elblock) % tags(msh_el), msh_element_blocks(msh_elblock) % nodes(msh_el,:)
         end do ! msh_els_per_block /for tags

         msh_global_el = msh_global_el + msh_els_per_block

      end do ! msh_no_elblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected." 
      close( fUnit )
!------------------------------------------------------------------------

!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype(1) = count(msh_element_blocks(:) % el_type .eq.  5)
      check_eltype(2) = count(msh_element_blocks(:) % el_type .eq. 12)
      check_eltype(3) = count(msh_element_blocks(:) % el_type .eq. 92)
      check_eltype(4) = count(msh_element_blocks(:) % el_type .eq. 93)
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      check_eltype = (/5,12,92,93/)
      org_element_type = check_eltype(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do
!------------------------------------------------------------------------

!-----Reorder-nodes-in-elements------------------------------------------
      select case(org_element_type)
      case (5)
         allocate(tmpi_vec1(8))
      case (12)
         error stop "READ_GMSH :: TBC"
      case (92)
         error stop "READ_GMSH :: TBC"
      case (93)
         error stop "READ_GMSH :: TBC"
      case default 
         error stop "READ_GMSH :: No 3D elements detected in the mesh."
      end select
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els 
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               msh_element_blocks(msh_elblock) % nodes(j,1) = tmpi_vec1(4)
               msh_element_blocks(msh_elblock) % nodes(j,2:4) = tmpi_vec1(1:3)
               msh_element_blocks(msh_elblock) % nodes(j,5) = tmpi_vec1(8)
               msh_element_blocks(msh_elblock) % nodes(j,6:8) = tmpi_vec1(5:7)
            end do
         end if
      end do
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-BCs---------------------------------------------------------
      do i=1,msh_no_BCs
         ! find no surfaces
         msh_bcs(i) % no_of_surfaces = 0
         do j=1, msh_no_surfaces
            tmpi = findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
         end do ! msh_no_surfaces

         ! assign surfaces and curves
         allocate(tmpi_vec1(4 * msh_bcs(i) % no_of_surfaces))
         allocate(msh_bcs(i) % surface_tags(msh_bcs(i) % no_of_surfaces))
         k = 0
         do j=1, msh_no_surfaces
            tmpi = findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               msh_bcs(i) % surface_tags(k) = msh_surfaces(j) % tag
               tmpi_vec1(1 + 4*(k-1) : 4*k) = msh_surfaces(j) % bps 
            end if ! tmpi
         end do ! msh_no_surfaces
         call unique(tmpi_vec1,msh_bcs(i) % curve_tags)
         msh_bcs(i) % no_of_curves = size(msh_bcs(i) % curve_tags,1)
         deallocate(tmpi_vec1)

         ! assign points
         k = 0
         allocate(tmpi_vec1(2 * msh_bcs(i) % no_of_curves))
         do j=1, msh_no_curves
            tmpi = findloc(msh_bcs(i) % curve_tags, msh_curves(j)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               tmpi_vec1(1 + (k-1)*2: k*2) = msh_curves(j) % bps
            end if ! tmpi
         end do ! msh_no_points
         call unique(tmpi_vec1,msh_bcs(i) % point_tags)
         msh_bcs(i) % no_of_points = size(msh_bcs(i) % point_tags,1)
         deallocate(tmpi_vec1)

         ! assign nodes
         k = 0
         allocate(tmpi_vec1(numberOfNodes)) ! approximate upper bound
         tmpi_vec1 = 0
         do j=1, msh_no_nodeblocks

            select case (msh_node_blocks(j) % entity_dim)
            case (0) ! points 
               tmpi = findloc(msh_bcs(i) % point_tags, msh_node_blocks(j) % entity_tag, 1)
            case (1) ! curves 
               tmpi = findloc(msh_bcs(i) % curve_tags, msh_node_blocks(j) % entity_tag, 1)
            case (2) ! surfaces
               tmpi = findloc(msh_bcs(i) % surface_tags, msh_node_blocks(j) % entity_tag, 1)
            case (3) ! volume
               tmpi = 0
            end select ! msh_node_blocks(j) % entity_dim

            if (tmpi .gt. 0) then 
               tmpi_vec1(k+1 : k+msh_node_blocks(j) % no_nodes) = msh_node_blocks(j) % tags
               k = k + msh_node_blocks(j) % no_nodes
            end if ! tmpi

         end do ! msh_no_nodeblocks
         call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
         msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
         deallocate(tmpi_vec1)

      end do ! msh_no_BCs
!------------------------------------------------------------------------

!-----Assign-BC-to-nodes-------------------------------------------------
      allocate(tmpi_vec1(8)) ! help array for the walls
      do i=1, msh_no_BCs
         do msh_elblock=1, msh_no_elblocks
            if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
               do j=1, msh_element_blocks(msh_elblock) % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = findloc(msh_bcs(i) % node_tags,msh_element_blocks(msh_elblock) % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes 
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,4) = msh_bcs(i) % tag ! 4:front 
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_element_blocks(msh_elblock) % no_els 
            end if
         end do ! msh_no_elblocks
      end do ! msh_no_BCs
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-new-tags----------------------------------------------------
      tmp_eltag = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els 
               tmp_eltag = tmp_eltag + 1
               msh_element_blocks(msh_elblock) % tags(j) = tmp_eltag
            end do
         end if
      end do
      if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!------------------------------------------------------------------------

!-----Build-nodes--------------------------------------------------------
      self % nodeType = nodes
      self % no_of_elements = numberOfElements
!------------------------------------------------------------------------

!-----Set-up-face-patches------------------------------------------------
      numBFacePoints = bFaceOrder + 1
      allocate(uNodes(numBFacePoints))
      allocate(vNodes(numBFacePoints))
      allocate(values(3,numBFacePoints,numBFacePoints))
      do i = 1, numBFacePoints
         uNodes(i) = -cos((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         vNodes(i) = -cos((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
      end do
!------------------------------------------------------------------------

!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(numberOfelements) )
      allocate( self % nodes(numberOfNodes) )
      allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
!------------------------------------------------------------------------

!----Set-nodes-----------------------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref
            call ConstructNode( self % nodes(msh_node_blocks(msh_nodeblock) % tags(msh_node)), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-elements-----------------------------------------------------------
      j = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do msh_el = 1, msh_element_blocks(msh_elblock) % no_els 
               j = j + 1
               ! setting l'th element
               l =  msh_element_blocks(msh_elblock) % tags(msh_el)
               nodeIDs = msh_element_blocks(msh_elblock) % nodes(msh_el,1:8) ! only non-curved nodes

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(l) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(l) % SurfInfo % corners = corners
               else ! curved mesh

               end if

               call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = findloc( msh_bcs % tag,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(l) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(l) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(l) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(l) % boundaryName(k) )
                     numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(l) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(l) % boundaryName(k)), trim(self % elements(l) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k
            end do ! msh_element_blocks(msh_elblock) % no_els
         end if ! if el_type .eq. org_element_type
      end do ! msh_no_elblocks
      if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."
!------------------------------------------------------------------------

!-----Deallocate-msh-vairables-------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         call msh_node_blocks(msh_nodeblock) % Destruct()
      end do
      deallocate(msh_node_blocks)

      do msh_elblock=1, msh_no_elblocks
         call msh_element_blocks(msh_elblock) % Destruct()
      end do
      deallocate(msh_element_blocks)

      do i=1,msh_no_BCs
         call msh_bcs(i) % Destruct()
      end do
      deallocate(msh_bcs)

      deallocate(msh_points)
      deallocate(msh_curves)
      deallocate(msh_surfaces)
      deallocate(msh_volumes)
      deallocate(msh_entity_vec)

      deallocate(uNodes) ! Check if we can do it! FIXME
      deallocate(vNodes) ! Check if we can do it! FIXME
      deallocate(values) ! Check if we can do it! FIXME
!------------------------------------------------------------------------
!
!     ---------------------------
!     Construct the element faces
!     ---------------------------
!
      numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
      self % numberOfFaces = numberOfFaces
      allocate( self % faces(self % numberOfFaces) )
      CALL ConstructFaces( self, success )
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % ConstructZones()
!
!     ---------------------------
!     Construct periodic faces
!     ---------------------------
!
      CALL ConstructPeriodicFaces( self ) 
!
!     ---------------------------
!     Delete periodic- faces
!     ---------------------------
!
      CALL DeletePeriodicMinusFaces( self ) 
!
!     ---------------------------
!     Assign faces ID to elements
!     ---------------------------
!
      CALL getElementsFaceIDs(self)
!
!     ---------------------
!     Define boundary faces
!     ---------------------
!
      call self % DefineAsBoundaryFaces()
!
!     ------------------------------
!     Set the element connectivities
!     ------------------------------
!
      call self % SetConnectivitiesAndLinkFaces(nodes)
      
      call self % ExportBoundaryMesh (trim(fileName))
   end subroutine ConstructSimplestMesh_FromGMSHFile_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
    function NumOfElems_GMSH( fileName ) result(nelem)
        implicit none
        !----------------------------------
        CHARACTER(LEN=*), intent(in) :: fileName
        integer                      :: nelem
        !----------------------------------
        integer :: i, j, k, fUnit, nElBlocks, ElsInBlock, eltype
        integer :: ierr
        logical :: success=.false.
        character (len=100) :: cline
        character (len=10)  :: cword
        integer, dimension(4) :: check_eltype = (/5,12,92,93/)
        !----------------------------------
         
        open(newunit = fUnit, FILE = trim(fileName) )  

        do
            read(fUnit,'(A)',iostat=ierr) cline
            if (ierr /= 0) exit
            read (cline,*) cword ! read first word of line
            if (cword == "$Elements") then
               read(fUnit,*) nElBlocks, nelem, k, k           
               nelem = 0
               do i=1, nElBlocks
                  read(fUnit,*) k, k, eltype, ElsInBlock
                  if (findloc(check_eltype,eltype,1) ) nelem = nelem + ElsInBlock
                  do j=1, ElsInBlock
                     read(fUnit,*) k
                  end do
               end do
               read(fUnit,*) cline
               if (trim(cline) .ne. '$EndElements') error stop "NumOfElems_GMSH :: Wrong input file - not all elements read." 
               success=.true.
               exit
            end if
        end do

        if (success) then
        else
            print *, cline
            error stop "NumOfElems_GMSH :: Could not read number of elements."
        end if

        close(fUnit)

    end function NumOfElems_GMSH
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_DestructBCStorage(this)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_BCinfo_t)            , intent(inout)   :: this 
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      if(allocated(this % element_tags))  deallocate(this % element_tags)
      if(allocated(this % wall_side))     deallocate(this % wall_side)
      if(allocated(this % node_tags))     deallocate(this % node_tags)
      if(allocated(this % volume_tags))   deallocate(this % volume_tags)
      if(allocated(this % surface_tags))  deallocate(this % surface_tags)
      if(allocated(this % curve_tags))    deallocate(this % curve_tags)
      if(allocated(this % point_tags))    deallocate(this % point_tags)

   end subroutine MSH_DestructBCStorage
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_ConstructNodeBlock(this,edim,etag,par,no_nodes)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_node_block_t)            , intent(inout)   :: this ! element block to be constructed
      integer                , intent(in)      :: edim ! dimension of the entity
      integer                , intent(in)      :: etag ! entity tag
      logical                , intent(in)      :: par ! parametric or not (TODO: parametric not supported!)
      integer                , intent(in)      :: no_nodes ! number of nodes within the block
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      this % entity_dim = edim
      this % entity_tag = etag
      this % parametric = par
      this % no_nodes   = no_nodes
      allocate(this % tags(no_nodes))
      allocate(this % cords(no_nodes,3))

   end subroutine MSH_ConstructNodeBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_DestructNodeBlock(this)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_node_block_t)            , intent(inout)   :: this ! element block to be constructed
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      deallocate(this % tags)
      deallocate(this % cords)

   end subroutine MSH_DestructNodeBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_ConstructElementBlock(this,el_type,no_els)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_element_block_t)            , intent(inout)   :: this ! element block to be constructed
      integer                , intent(in)      :: el_type ! type of element
      integer                , intent(in)      :: no_els ! number of elements
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      this % el_type = el_type
      this % no_els = no_els
      allocate(this % tags(no_els))
      select case (el_type)               
      case (3)  ! 2D 1st order quadrangle
         allocate(this % nodes(no_els,4))
      case (10) ! 2D 2nd order quadrangle
         allocate(this % nodes(no_els,9))
      case (5)  ! 3D 1st order hexahedron
         allocate(this % nodes(no_els,8))
      case (12) ! 3D 2nd order hexahedron
         allocate(this % nodes(no_els,27))
      case (92) ! 3D 3rd order hexahedron
         allocate(this % nodes(no_els,64))
      case (93) ! 3D 4th order hexahedron
         allocate(this % nodes(no_els,125))
      case default
         print *, " READ_GMSH :: Wrong element type. Only proper tags will be stored."
         allocate(this % nodes(no_els,1))
      end select

      ! allocate BCs only if 3D element is detected
      if (findloc((/5,12,92,93/),el_type,1)) then
         allocate(this % BCs(no_els,6))
         this % BCs = 0
      end if
   end subroutine MSH_ConstructElementBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_DestructElementBlock(this)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_element_block_t)            , intent(inout)   :: this ! element block to be constructed
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      deallocate(this % tags)
      deallocate(this % nodes)
      if(allocated(this % BCs)) deallocate(this % BCs)

   end subroutine MSH_DestructElementBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine unique(vec,vec_unique)
      ! CREDIT: http://degenerateconic.com/unique/ email: jacob[at]degenerateconic[dot]com
      ! Return only the unique values from vec.
       
      implicit none
       
      integer,dimension(:),intent(in) :: vec
      integer,dimension(:),allocatable,intent(out) :: vec_unique
       
      integer :: i,num
      logical,dimension(size(vec)) :: mask
       
      mask = .false.
       
      do i=1,size(vec)
       
          !count the number of occurrences of this element:
          num = count( vec(i)==vec )
       
          if (num==1) then
              !there is only one, flag it:
              mask(i) = .true.
          else
              !flag this value only if it hasn't already been flagged:
              if (.not. any(vec(i)==vec .and. mask) ) mask(i) = .true.
          end if
       
      end do
       
      !return only flagged elements:
      allocate( vec_unique(count(mask)) )
      vec_unique = pack( vec, mask )
       
      !if you also need it sorted, then do so.
      ! For example, with slatec routine:
      !call ISORT (vec_unique, [0], size(vec_unique), 1)
       
   end subroutine unique
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module Read_GMSH
