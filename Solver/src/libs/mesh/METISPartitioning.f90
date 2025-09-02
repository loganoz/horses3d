subroutine GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain, useWeights, controlVariables, &
										nLevel, eID_Order, nElementLevel)
!
!     *********************************************************************
!        This subroutine performs the mesh partitioning using METIS
!        and then manual load balancing for mixed RK used in multiphase
!     *********************************************************************
!        
      use HexMeshClass
      use SMConstants
      use MPI_Process_Info
      use FTValueDictionaryClass
	  use Utilities                   , only : reindexIntegerList, combine_partitions
      implicit none
      type(HexMesh), intent(in)              :: mesh
      integer,       intent(in)              :: no_of_domains
      integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
      integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
      logical,       intent(in)              :: useWeights
      type(FTValueDictionary), intent(in)    :: controlVariables
	  integer,       intent(in)              :: nLevel 
	  integer,       intent(in)      	     :: eID_Order(mesh % no_of_elements)
	  integer,       intent(in)      	     :: nElementLevel(nLevel)
!
!     ---------------
!     Local Variables
!     ---------------
!
      integer                :: i, j, domain_id, level, starteID
      integer                :: ielem
      integer                :: ne                    ! # elements
      integer                :: nn                    ! # nodes
      integer                :: nvertex               ! # vertices per element
      integer, allocatable   :: eptr(:)               ! index in eind --> eptr(ne+1)
      integer, allocatable   :: eind(:)               ! vertices of each element   --> eind(nvertex*ne)
      integer, pointer       :: vwgt(:) => null()     ! vertices weights
      integer, allocatable   :: vsize(:) 
      integer, parameter     :: ncommon = 4           ! common faces for dual nesh (4 for hexahedrals)
      real(kind=RP), allocatable  :: tpwgt(:)         ! partitions' weights --> tpwgt(no_of_domains)
      integer, allocatable   :: opts(:)               ! metis options
      integer, parameter     :: metis_noptions = 39   ! number of metis options
      integer                :: objval                ! objective calculated value
      integer, parameter     :: MIN_EDGE_CUT = -1     ! option to minimize edge-cut
      integer, parameter     :: MIN_COMM_VOL = 1      ! option to minimize the communication volume
      integer, allocatable, target   :: weights(:)
      logical                :: is_water              ! Flag for water elements
      integer                :: water_count, air_count ! Counters for water and air elements
      integer, allocatable   :: domain_water_count(:), domain_air_count(:)
      integer, allocatable   :: water_elems(:), air_elems(:)
      integer                :: target_water, target_air, w_idx, a_idx
      logical, allocatable   :: elem_moved(:)
      logical                :: do_water_air_count    ! Flag to control water/air counting
	  integer, allocatable   :: elementsDomainLevel(:), nodesDomainLevel(:), mapToOld(:)
	  integer, allocatable   :: refEleDomain(:), inputMLRKDomain(:,:)
	  
#ifdef _HAS_METIS_
!
!     **************
!     Initialization    
!     **************
!
      ne = mesh % no_of_elements
      nn = size(mesh % nodes)
      nvertex = 8
	  
	  allocate(inputMLRKDomain(mesh % no_of_elements, nLevel))
	  inputMLRKDomain = 0
	
      ! Check if we should do water/air counting based on explicit method
      do_water_air_count = .false.
      if(trim(controlVariables % stringValueForKey('explicit method', requestedLength = LINE_LENGTH)) == 'mixed rk') then
          do_water_air_count = .true.
      end if

      ! Count water and air elements only if needed
      water_count = 0
      air_count = 0
      if (do_water_air_count) then
          do ielem=1,ne
              is_water = .false.
              if (associated(mesh % elements(ielem) % storage)) then
                    is_water = all(mesh % elements(ielem) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP)
              end if

              if (is_water) then
                  water_count = water_count + 1
              else
                  air_count = air_count + 1
              end if
          end do

          if (MPI_Process % isRoot) then
              write(*,*) "From inside METIS partitioning "
              write(*,'(A,I6,A,I6,A)') 'Identified ', water_count, ' water elements and ', air_count, ' air elements'
          end if
      end if
!
!     Construct the reference element partition from single level for MLRK
!     --------------------------------------------------------------------      
	  if (nLevel.gt.1) then
	    allocate(refEleDomain(ne))
		allocate(eptr(ne+1), eind(nvertex*ne))
		
		refEleDomain    =0
!
!     	Gather each element nodes: they are stored using C indexes (starting on 0)
!     	-------------------------
		i = 1
		do ielem=1,ne
!
!        	Save each element nodes
!        	-----------------------
			eind(i:i+nvertex-1) = mesh % elements(ielem) % nodeIDs - 1
!
!        	Save each element ID
!        	--------------------
			eptr(ielem) = i - 1
			i = i + nvertex
		end do
!
!     	Termination: set the last element position in the N+1 entry
!     	-----------
		eptr(ne+1) = i - 1
!
!     	*****************
!     	Set METIS options
!     	*****************
!
		allocate(opts(0:metis_noptions))
		call METIS_SetDefaultOptions(opts)
!
!     	First option chooses the method: edge-cut / communication volume      
!     	-------------------------------
		opts(1) = MIN_EDGE_CUT
!
!     	Disable verbosity
!     	-----------------
		opts(5) = 0
!
!     	*******************************
!     	Calculate weights based on NDOF
!     	*******************************
!
		if (useWeights) then
		  allocate(weights(ne))
		  do ielem=1,ne
			  ! Base weight from polynomial order
			  weights(ielem) = product(mesh % elements(ielem) % Nxyz + 1)
		  end do
		  if (maxval(weights) .ne. minval(weights)) then
			  vwgt => weights
		  endif
		end if 
	    
		call METIS_PartMeshDual(ne, nn,  eptr, eind,  vwgt, vsize, ncommon, no_of_domains, tpwgt,  opts,  objval, refEleDomain, nodesDomain)
!
!     	Recover FORTRAN displacements by adding 1
!     	-----------------------------------------
		refEleDomain = refEleDomain + 1
!
!       Free memory
!       -----------
!
		deallocate (eptr, eind, opts)
		if (associated(vwgt)) nullify(vwgt)       ! vwgt is a pointer to a target. nullify is enough
		if (allocated(weights)) deallocate(weights)
		if (allocated(tpwgt)) deallocate(tpwgt)
		if (allocated(vsize)) deallocate(vsize)
	  end if 
!
!     Perform METIS partitioning based on the element's level - if not MLRK then nLevel=1
!     -----------------------------------------------------------------------------------
	  starteID = 1								   ! counter of the elementID for multi level 
	  do level =1, nLevel
		ne = nElementLevel(level)
		allocate(eptr(ne+1))
		allocate(eind(nvertex*ne))
!
!     	Gather each element nodes: they are stored using C indexes (starting on 0)
!     	-------------------------------------------------------------------------
		i = 1
		j = 0
		do ielem=starteID, ne+starteID-1
		    j=j+1
!
!        	Save each element nodes
!        	-----------------------
			eind(i:i+nvertex-1) = mesh % elements(eID_Order(ielem)) % nodeIDs - 1
!
!        	Save each element ID
!        	--------------------
			eptr(j) = i - 1
			i = i + nvertex
		end do
!
!     	reindexIntegerList as such it is compactly renumbered from 0 to totalNodes-1
!     	----------------------------------------------------------------------------
        allocate(mapToOld(size(mesh % nodes)))
		mapToOld = [(i, i=1,size(mesh % nodes))]
		if (nLevel.gt.1) then
		    deallocate(mapToOld)
			call reindexIntegerList(eind,nn,mapToOld)
		end if
		allocate(elementsDomainLevel(1:ne), nodesDomainLevel(1:nn))
!
!     	Termination: set the last element position in the N+1 entry
!     	-----------------------------------------------------------
		eptr(ne+1) = i - 1
!
!     	*****************
!     	Set METIS options
!     	*****************
!
		allocate(opts(0:metis_noptions))
		call METIS_SetDefaultOptions(opts)
!
!     	First option chooses the method: edge-cut / communication volume      
!     	----------------------------------------------------------------
		opts(1) = MIN_EDGE_CUT
!
!     	Disable verbosity
!     	-----------------
		opts(5) = 0
!
!     	*******************************
!     	Calculate weights based on NDOF
!     	*******************************
!
		if (useWeights) then
		  allocate(weights(ne))
		  do ielem=1,ne
			  ! Base weight from polynomial order
			  weights(ielem) = product(mesh % elements(eID_Order(ielem+starteID-1)) % Nxyz + 1)
			  
			  ! Determine if this is a water element and apply weight factor only if needed
			  if (do_water_air_count) then
				  is_water = .false.
				  if (associated(mesh % elements(eID_Order(ielem+starteID-1)) % storage)) then
						  is_water = all(mesh % elements(eID_Order(ielem+starteID-1)) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP) 
				  end if
				  
				  ! Apply 4.666x weight factor for water elements
				  ! Not fully necessary because of the manual air-watter balancing done later in the function
				  ! but this is to protect against edge cases
				  if (is_water) then
					  weights(ielem) = weights(ielem) * 4.6666
				  end if
			  end if
		  end do

		  if (maxval(weights) .ne. minval(weights)) then
			  vwgt => weights
		  endif
		end if 
!     	**********************
!     	Perform the partitions
!     	**********************
!
		call METIS_PartMeshDual(ne, nn,  eptr, eind,  vwgt, vsize, ncommon, no_of_domains, tpwgt,  opts,  objval, elementsDomainLevel, nodesDomainLevel)
!
!     	Recover FORTRAN displacements by adding 1
!     	-----------------------------------------
		elementsDomainLevel = elementsDomainLevel + 1
		nodesDomainLevel    = nodesDomainLevel + 1
!
!     	Send information to inputMLRKDomain - relocate to true elementID
!     	----------------------------------------------------------------
        do i=1, ne
			inputMLRKDomain(eID_Order(starteID+i-1), level) = elementsDomainLevel(i)
		end do 
		starteID = starteID + ne
!
!       Free memory
!       -----------
!
		deallocate (nodesDomainLevel, elementsDomainLevel, mapToOld, eptr, eind, opts)
		if (associated(vwgt)) 	nullify(vwgt)       ! vwgt is a pointer to a target. nullify is enough
		if (allocated(weights)) deallocate(weights)
		if (allocated(tpwgt)) 	deallocate(tpwgt)
		if (allocated(vsize)) 	deallocate(vsize)
	 end do
	 
	 nodesDomain = 0 ! Not used 
	 do i=1, mesh % no_of_elements
		elementsDomain(i) = sum(inputMLRKDomain(i,:))
	 end do 
!
!    Combine_partitions between MLRK levels as such it has minimum MPI faces
!    -----------------------------------------------------------------------
	 if (nLevel.gt.1) then
		call combine_partitions(mesh % no_of_elements, nLevel, no_of_domains, refEleDomain, inputMLRKDomain, elementsDomain)
	 end if 
!    ************************************
!    Manual balancing of water and air. 
!    The aim is to have equal number of air and water elements in each core.
!    ************************************
     if (do_water_air_count .and. water_count > 0 .and. air_count > 0) then
          ! Allocate arrays for balancing
          allocate(domain_water_count(no_of_domains), domain_air_count(no_of_domains))
          allocate(water_elems(water_count), air_elems(air_count))
          allocate(elem_moved(ne))
          elem_moved = .false.
          
          ! Count water/air elements per domain and collect element lists
          domain_water_count = 0
          domain_air_count = 0
          w_idx = 0
          a_idx = 0
          
          do ielem = 1, ne
              is_water = .false.
              if (associated(mesh % elements(ielem) % storage)) then
                  is_water = all(mesh % elements(ielem) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP)
              end if
              
              if (is_water) then
                  w_idx = w_idx + 1
                  water_elems(w_idx) = ielem
                  domain_water_count(elementsDomain(ielem)) = domain_water_count(elementsDomain(ielem)) + 1
              else
                  a_idx = a_idx + 1
                  air_elems(a_idx) = ielem
                  domain_air_count(elementsDomain(ielem)) = domain_air_count(elementsDomain(ielem)) + 1
              end if
          end do
          
          ! Calculate target counts per domain
          target_water = water_count / no_of_domains
          target_air = air_count / no_of_domains
          
          if (MPI_Process % isRoot) then
              write(*,*) "Balancing water/air elements: target per domain =", target_water, target_air
          end if
          
          ! Redistribute elements to achieve balance while minimizing changes
          do domain_id = 1, no_of_domains
              ! Move water elements if needed
              do while (domain_water_count(domain_id) > target_water + 1)
                  ! Find a domain that needs more water elements
                  do j = 1, no_of_domains
                      if (domain_water_count(j) < target_water) then
                          ! Find an unmoved water element from this domain
                          do i = 1, water_count
                              if (elementsDomain(water_elems(i)) == domain_id .and. .not. elem_moved(water_elems(i))) then
                                  elementsDomain(water_elems(i)) = j
                                  domain_water_count(domain_id) = domain_water_count(domain_id) - 1
                                  domain_water_count(j) = domain_water_count(j) + 1
                                  elem_moved(water_elems(i)) = .true.
                                  exit
                              end if
                          end do
                          exit
                      end if
                  end do
                  
                  ! Break if we can't balance further
                  if (all(domain_water_count >= target_water)) exit
              end do
              
              ! Move air elements if needed
              do while (domain_air_count(domain_id) > target_air + 1)
                  ! Find a domain that needs more air elements
                  do j = 1, no_of_domains
                      if (domain_air_count(j) < target_air) then
                          ! Find an unmoved air element from this domain
                          do i = 1, air_count
                              if (elementsDomain(air_elems(i)) == domain_id .and. .not. elem_moved(air_elems(i))) then
                                  elementsDomain(air_elems(i)) = j
                                  domain_air_count(domain_id) = domain_air_count(domain_id) - 1
                                  domain_air_count(j) = domain_air_count(j) + 1
                                  elem_moved(air_elems(i)) = .true.
                                  exit
                              end if
                          end do
                          exit
                      end if
                  end do
                  
                  ! Break if we can't balance further
                  if (all(domain_air_count >= target_air)) exit
              end do
          end do
          
          ! Update node domains based on element assignments
          nodesDomain = 1
          do ielem = 1, ne
              do i = 1, nvertex
                  nodesDomain(mesh % elements(ielem) % nodeIDs(i)) = elementsDomain(ielem)
              end do
          end do
!
!         Free memory
!         -----------
!
          deallocate(domain_water_count, domain_air_count)
          deallocate(water_elems, air_elems)
          deallocate(elem_moved)
          
          if (MPI_Process % isRoot) then
              write(*,*) "Applied post-processing to balance water/air elements"
          end if
	  end if

!     **********************
!     Print domain statistics
!     **********************
	  if (do_water_air_count .and. MPI_Process % isRoot) then
          do domain_id = 1, no_of_domains
              water_count = 0
              air_count = 0
              do ielem = 1, ne
                  if (elementsDomain(ielem) == domain_id) then
                      ! Determine if this is a water element
                      is_water = .false.
                      if (associated(mesh % elements(ielem) % storage)) then
                            is_water = all(mesh % elements(ielem) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP)
                      end if
                      
                      if (is_water) then
                          water_count = water_count + 1
                      else
                          air_count = air_count + 1
                      end if
                  end if
              end do
              
              write(*,'(A,I3,A,I5,A,I5,A,I5)') 'Domain ', domain_id, &
                    ': Water elements = ', water_count, &
                    ', Air elements = ', air_count, &
                    ', Total = ', water_count + air_count
          end do
	  end if
!
!     ****
!     Free
!     ****
!
      if (allocated(weights)) deallocate(weights)
	  if (allocated(inputMLRKDomain)) deallocate(inputMLRKDomain)
	  if (allocated(refEleDomain)) deallocate(refEleDomain)
#endif
   end subroutine GetMETISElementsPartition
