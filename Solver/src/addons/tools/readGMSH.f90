#include "Includes.h"
module readGMSH
   use HexMeshClass
   use ElementClass
   use SMConstants
   use ElementConnectivityDefinitions
   use LocalRefinement
   use Utilities, only: UnusedUnit
   use NodeClass
   use Read_GMSH
   use FileReadingUtilities      , only: getFileName, getRealArrayFromStringNoCommas

   implicit none

   contains

   subroutine ConstructSimpleMesh_FromGMSHFile_v4_(self, fileName, locR, AllNx, AllNy, AllNz)
   
      !use PhysicsStorage
   
      implicit none
!        ---------------
!        Input variables
!        ---------------
!
      type(HexMesh)                    :: self
      character(len=*)                 :: fileName
      type(LocalRef_t), optional, intent(in)  :: locR
      integer, optional,intent(in)     :: AllNx(:), AllNy(:), AllNz(:)
!
!        ---------------
!        Local variables
!        ---------------
!
      integer                         :: Nx, Ny, Nz     !<  Polynomial orders for each element
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      integer                         :: i, j, k, l
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      integer                         :: funit, fileStat
      
      real(kind=RP), allocatable :: msh_entity_vec(:)
      
      logical                         :: tmpb
      
      integer, dimension(EL_MAX_ORDER)   :: check_eltype      
      integer, dimension(:), allocatable :: tmpi_vec1, el_types
      
      integer                    :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type, org_element_type_2D
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el
      integer                    :: numberOfElements, numberOfElements2D
      integer                    :: numberOfNodes
      integer                    :: bFaceOrder, no_nodes_i

      
      type(MSH_node_block_t),    dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t), dimension(:), allocatable  :: msh_element_blocks
      type(MSH_point_t),         dimension(:), allocatable  :: msh_points
                 
      real(kind=RP)          :: corners(NDIM,NODES_PER_ELEMENT)
      real(kind=RP)          :: x(NDIM), L_ref
      integer                :: nodeIDs(NODES_PER_ELEMENT)
      character(len=MSH_LEN) :: msh_entity

      L_ref = 1.0_RP
       
      funit = UnusedUnit()
      open( unit = fUnit, file = fileName, iostat = fileStat )
      if ( fileStat /= 0 ) then
         print *, "Error opening file: ", fileName
         return
      end if

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
      read(fUnit,*) tmpd, tmpi, tmpi
      read(fUnit,*) tmps

!-----Read-BC-info-------------------------------------------------------
      ! Bc not required, they are skipped
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
      read(fUnit,*) tmpi

      do i=1, tmpi
         read(fUnit,*) 
      end do

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found."
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      ! msh_no_curves, msh_no_surfaces, msh_no_volumes not required
      allocate(msh_points(msh_no_points))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:5+msh_points(i)%no_ptags))
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      ! Curves not needed
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') 
      end do ! msh_no_curves
!-----Read-surfaces------------------------------------------------------
      ! Surfaces not needed
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') 
      end do ! msh_no_surfaces
!-----Read-volumes-------------------------------------------------------
      ! Volumes not needed
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') 
      end do ! msh_no_volumes
      
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected."

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
      
!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype = 0
      do i=1, EL_MAX_ORDER
         do msh_elblock=1, msh_no_elblocks
            if( msh_element_blocks(msh_elblock)% el_type .eq. SUPPORTED_EL_TYPES(i) ) then !count(msh_element_blocks(:)% el_type .eq. SUPPORTED_EL_TYPES(i))
               check_eltype(i) =  check_eltype(i) + 1
            end if
         end do
      end do
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do         
         
!-----Reorder-nodes-in-elements------------------------------------------
      allocate(tmpi_vec1((bFaceOrder + 1)**3))
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_element_blocks(msh_elblock) % nodes(j,:) = tmpi_vec1
            end do
         end if
      end do
      deallocate(tmpi_vec1)
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

!-----Build-nodes--------------------------------------------------------
      self % no_of_elements = numberOfElements      
!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(numberOfelements) )
      allocate( self % nodes(numberOfNodes) )
      allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
      if (present(AllNx)) then
          self % Nx = AllNx
          self % Ny = AllNy
          self % Nz = AllNz
      end if 
      
!----Set-nodes-----------------------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)!/L_ref
            call ConstructNode( self % nodes(msh_node_blocks(msh_nodeblock) % tags(msh_node)), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
         end do 
      end do 

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
                   print*, "Curved mesh for IBM p-adaptation not supported"
                   error stop
                end if
                
            if( present(locR) ) then 
               call locR % getOrderOfPosition(corners, Nx, Ny, Nz)
               ! call self % elements(l) % Construct (Nx, Ny, Nz, falseNodeID , l, l) 
               call self % elements(l) % Construct (Nx, Ny, Nz, nodeIDs, l, l) 
            else
               call self % elements(l) % Construct (AllNx(l), AllNy(l), AllNz(l), nodeIDs , l, l)
            end if

             end do ! msh_element_blocks(msh_elblock) % no_els
          end if ! if el_type .eq. org_element_type
       end do ! msh_no_elblocks
       if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."      

   
!------Deallocate-msh-vairables-------------------------------------------
       do msh_nodeblock=1, msh_no_nodeblocks
          call msh_node_blocks(msh_nodeblock) % Destruct()
       end do
       deallocate(msh_node_blocks)

       do msh_elblock=1, msh_no_elblocks
          call msh_element_blocks(msh_elblock) % Destruct()
       end do
       deallocate(msh_element_blocks)

       deallocate(msh_points)
       deallocate(msh_entity_vec)
      
    end subroutine ConstructSimpleMesh_FromGMSHFile_v4_
    
    subroutine ConstructSimpleMesh_FromGMSHFile_v2_(self, filename, locR, AllNx, AllNy, AllNz)
    
                USE Physics
            use PartitionedMeshClass
            use MPI_Process_Info
    
       implicit none
!        ---------------
!        Input variables
!        ---------------
!
       type(HexMesh)                    :: self
       character(len=*)                 :: fileName
       type(LocalRef_t), optional, intent(in)  :: locR
       integer, optional,intent(in)     :: AllNx(:), AllNy(:), AllNz(:)
!
!        ---------------
!        Local variables
!        ---------------
!
      integer                         :: Nx, Ny, Nz     !<  Polynomial orders for each element
       character(len=1024)             :: tmps
       real(kind=RP)                   :: tmpd
       integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
       integer                         :: i, j, k, l
       integer                         :: msh_no_points
       integer                         :: msh_no_curves
       integer                         :: msh_no_surfaces
       integer                         :: msh_no_volumes
       integer                         :: funit, fileStat
      
       real(kind=RP),      allocatable :: msh_entity_vec(:)
      
       integer, dimension(EL_MAX_ORDER)   :: check_eltype      
       integer, dimension(:), allocatable :: tmpi_vec1, el_types
      
       integer                    :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_global_node
       integer                    :: element_type, org_element_type, org_element_type_2D
       integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el
       integer                    :: numberOfElements, numberOfElements2D
       integer                    :: numberOfNodes, numberOfElements_, numberOfElements2D_
       integer                    :: bFaceOrder, no_nodes_i, msh_node

       type(MSH_element_block_t)  :: msh_element_blocks
       type(MSH_point_t),         dimension(:), allocatable  :: msh_points
       
       type(MSH_node_block_t)     :: msh_nodes
       type(MSH_element_block_t)  :: msh_elements, msh_elements_3D, msh_elements_2D
                 
       real(kind=RP)          :: corners(NDIM,NODES_PER_ELEMENT)
       real(kind=RP)          :: x(NDIM)
       integer                :: nodeIDs(NODES_PER_ELEMENT)
       character(len=MSH_LEN) :: msh_entity

       funit = UnusedUnit()
       open( unit = fUnit, file = fileName, iostat = fileStat )
       if ( fileStat /= 0 ) then
          print *, "Error opening file: ", fileName
          return
       end if

!-----Read-header-info---------------------------------------------------
       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
       read(fUnit,*) tmpd, tmpi, tmpi
       read(fUnit,*) tmps
!-----Read-BC-info-------------------------------------------------------
     ! Bc not needed
       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
       read(fUnit,*) tmpi

       do i=1, tmpi
          read(fUnit,*) 
       end do ! tmpi

       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."

!-----Read-nodes---------------------------------------------------------
       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."
       read(fUnit,*) numberOfNodes

       call msh_nodes% Construct(0, 0, .false., numberOfNodes)

       do i=1, numberOfNodes
          read(fUnit,*) msh_nodes% tags(i), msh_nodes% cords(i,:)
       end do 

       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."

!-----Read-elements------------------------------------------------------
       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

       read(fUnit,*) numberOfElements
       call msh_elements% Construct(element_type,numberOfElements)
       allocate(msh_entity_vec(255)) ! arbitrary number
       allocate(el_types(numberOfElements)) ! arbitrary number
            
       do i=1, numberOfElements

          read(fUnit,'(4096a)') msh_entity ! read row

          msh_entity_vec=0.0_RP
               
          call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

          msh_elements% tags(i)     = int(msh_entity_vec(1))
          el_types(i)               = int(msh_entity_vec(2))
          msh_elements% no_ptags(i) = int(msh_entity_vec(3))
          if (msh_elements% no_ptags(i) .gt. 0) msh_elements% ptags(i,1:msh_elements% no_ptags(i)) = &
              int(msh_entity_vec(4:3+msh_elements % no_ptags(i)))

          select case (el_types(i))
             case (5) ! 3D - 1st order
                no_nodes_i = 8
             case (12) ! 3D - 2st order
                no_nodes_i = 27
             case (92) ! 3D - 3rd order
                no_nodes_i = 64
             case (93) ! 3D - 4th order
                no_nodes_i = 125
             case (94) ! 3D - 5th order
                no_nodes_i = 216
             case (3) ! 2D - 1st order
                no_nodes_i = 4
             case (10) ! 2D - 2st order
                no_nodes_i = 9
             case (36) ! 2D - 3rd order
                no_nodes_i = 16
             case (37) ! 2D - 4th order
                no_nodes_i = 25
             case (38) ! 2D - 5th order
                no_nodes_i = 36
             case default
                no_nodes_i = 0
          end select

          msh_elements% nodes(i,1:size(msh_entity_vec(4+msh_elements% no_ptags(i):3+msh_elements% no_ptags(i)+no_nodes_i) )) = &
          int(msh_entity_vec(4+msh_elements% no_ptags(i):3+msh_elements% no_ptags(i)+no_nodes_i))

       end do ! numberOfElements

       deallocate(msh_entity_vec)

       read(fUnit,*) tmps
       if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
       close( fUnit )

!-----Mesh-info----------------------------------------------------------
       ! find order of elements curvature
       check_eltype = 0
       do i=1, EL_MAX_ORDER
          do j = 1, numberOfElements
             if( el_types(j) .eq. SUPPORTED_EL_TYPES(i) )then
                check_eltype(i) = check_eltype(i) + 1
             end if
          end do
!~           check_eltype(i) = count(el_types .eq. SUPPORTED_EL_TYPES(i))
       end do
       if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
       if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
       bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
       org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
       msh_elements% el_type = org_element_type

       select case (bFaceOrder)
          case (1)
             org_element_type_2D = 3
          case (2)
             org_element_type_2D = 10
          case (3)
             org_element_type_2D = 36
          case (4)
             org_element_type_2D = 37
          case (5)
             org_element_type_2D = 38
          case default
       end select

!~        if (numberOfElements .ne. (count(el_types .eq. org_element_type) + count(el_types .eq. org_element_type_2D)) ) &
!~        error stop "READ_GMSH :: Too many different types of elements."
!~        numberOfElements = count(el_types .eq. org_element_type)
!~        numberOfElements2D = count(el_types .eq. org_element_type_2D)
       
       numberOfElements_ = 0
       numberOfElements2D_ = 0
       do i = 1, size(el_types)
          if( el_types(i) .eq. org_element_type ) then
             numberOfElements_ = numberOfElements_ + 1
          elseif( el_types(i) .eq. org_element_type_2D ) then
             numberOfElements2D_ = numberOfElements2D_ + 1
          end if
       end do
       
       if( numberOfElements .ne. (numberOfElements_ + numberOfElements2D_) ) &
       error stop "READ_GMSH :: Too many different types of elements."


       numberOfElements = numberOfElements_
       numberOfElements2D = numberOfElements2D_

       call msh_elements_3D % Construct(org_element_type,numberOfElements)
       call msh_elements_2D % Construct(org_element_type_2D,numberOfElements2D)

       j = 0
       k = 0
       do i=1, size(el_types)

          if ( el_types(i) .eq. org_element_type ) then
             j = j + 1
             msh_elements_3D % tags(j)     = msh_elements% tags(i)
             msh_elements_3D % no_ptags(j) = msh_elements% no_ptags(i)
             msh_elements_3D % ptags(j,:)  = msh_elements% ptags(i,1)
             msh_elements_3D % nodes(j,:)  = msh_elements% nodes(i,1:size(msh_elements_3D % nodes(j,:)))
          elseif ( el_types(i) .eq. org_element_type_2D ) then
             k = k + 1
             msh_elements_2D % tags(k)     = msh_elements% tags(i)
             msh_elements_2D % no_ptags(k) = msh_elements% no_ptags(i)
             msh_elements_2D % ptags(k,:)  = msh_elements% ptags(i,1)
             msh_elements_2D % nodes(k,:)  = msh_elements% nodes(i,1:size(msh_elements_2D % nodes(k,:)))
          else
             error stop "READ_GMSH :: Unknown element type in the mesh."
          end if

       end do

       deallocate(el_types)
       call msh_elements % Destruct()
!-----Reorder-nodes-in-elements------------------------------------------
       allocate(tmpi_vec1((bFaceOrder + 1)**3))
       do j=1, msh_elements_3D % no_els
          ! re-order nodes within element to match HORSES ordering
          tmpi_vec1 = msh_elements_3D % nodes(j,:)
          call ReorderElement(tmpi_vec1,bFaceOrder)
          msh_elements_3D % nodes(j,:) = tmpi_vec1
       end do
       deallocate(tmpi_vec1)
!-----Assign-new-tags----------------------------------------------------
       tmp_eltag = 0
       do j=1, msh_elements_3D % no_els
          tmp_eltag = tmp_eltag + 1
          msh_elements_3D % tags(j) = tmp_eltag
       end do
       if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!-----Build-nodes--------------------------------------------------------
       self % no_of_elements = numberOfElements
!-----Allocate-mem-for-elements-and-nodes--------------------------------
       allocate( self % elements(numberOfelements) )
       allocate( self % nodes(numberOfNodes) )
       allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
       self % Nx = AllNx
       self % Ny = AllNy
       self % Nz = AllNz
!~ !----Set-nodes-----------------------------------------------------------
       do msh_node=1, msh_nodes % no_nodes
          x = msh_nodes % cords(msh_node,1:NDIM)!/L_ref
          call ConstructNode( self % nodes(msh_nodes% tags(msh_node)), x, msh_nodes% tags(msh_node) )
       end do ! msh_nodes % no_nodes
!----Set-elements-----------------------------------------------------------
       j = 0
       do msh_el = 1, msh_elements_3D % no_els
          j = j + 1
          ! setting l'th element
          l =  msh_elements_3D % tags(msh_el)
          nodeIDs = msh_elements_3D % nodes(msh_el,1:8) ! only non-curved nodes

          if (bFaceOrder .eq. 1) then ! non-curved mesh
             do k = 1, NODES_PER_ELEMENT
                corners(:,k) = self % nodes(nodeIDs(k)) % x
             end do
             self % elements(l) % SurfInfo % IsHex8 = .TRUE.
             self % elements(l) % SurfInfo % corners = corners
          else ! curved mesh    
             print*, "Curved mesh for IBM p-adaptation not supported"
             error stop   
          end if

            if( present(locR) ) then 
               call locR % getOrderOfPosition(corners, Nx, Ny, Nz)
               ! call self % elements(l) % Construct (Nx, Ny, Nz, falseNodeID , l, l) 
               call self % elements(l) % Construct (Nx, Ny, Nz, nodeIDs, l, l) 
            else
               call self % elements(l) % Construct (AllNx(l), AllNy(l), AllNz(l), nodeIDs , l, l)
            end if

       end do ! msh_elements_3D % no_els
       if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."

!~ !-----Deallocate-msh-vairables-------------------------------------------
       call msh_nodes% Destruct()
       call msh_elements_3D % Destruct()
       call msh_elements_2D % Destruct()

   end subroutine ConstructSimpleMesh_FromGMSHFile_v2_

end module readGMSH
