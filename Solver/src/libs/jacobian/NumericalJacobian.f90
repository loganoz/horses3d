!
!//////////////////////////////////////////////////////
!
!      Routines for computing the Jacobian matrix numerically using the colorings technique
!  ! TODO1: Implement as a class with a destructor to prevent memory leaking
!////////////////////////////////////////////////////////////////////////
module NumericalJacobian
   use SMConstants
   use MatrixClass
   use ColorsClass            , only: Colors_t
   use HexMeshClass           , only: HexMesh, Neighbor_t, NUM_OF_NEIGHBORS
   use DGSEMClass             , only: DGSem, ComputeTimeDerivative_f
   use ElementClass
   use JacobianDefinitions    , only: JACEPS
   use JacobianComputerClass  , only: local2ijk, Look_for_neighbour, JacobianComputer_t
   use PhysicsStorage
   use Utilities              , only: Qsort, my_findloc
   use StorageClass           , only: SolutionStorage_t
   use IntegerDataLinkedList  , only: IntegerDataLinkedList_t
   use StopwatchClass         , only: StopWatch
   use BoundaryConditions     , only: NS_BC, C_BC, MU_BC
   use FTValueDictionaryClass
   use PartitionedMeshClass   , only: mpi_partition
   use ProgressBarsModule
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
   private
   public NumJacobian_t, GetRowsAndColsVector
   
!
!  *************************************************
!  Main type for the numerical Jacobian computations
!  *************************************************
   type, extends(JacobianComputer_t) :: NumJacobian_t
      
      contains
         procedure :: Construct           => NumJacobian_Construct
         procedure :: AssignColToJacobian => NumJacobian_AssignColToJacobian
         procedure :: Compute   => NumJacobian_Compute
   end type NumJacobian_t
   
!
!  Module variables
!  -> TODO: They will have to be moved to the class definition or to other types in the future
!  *******************
   type(Neighbor_t), allocatable :: nbr(:)  ! Neighbors information
   type(Neighbor_t), allocatable :: nbr_g(:)  ! Global neighbors array
   type(Colors_t)               :: ecolors
   type(Colors_t)               :: ecolors_g
   integer        , allocatable :: used(:)                  ! array containing index of elements whose contributions to Jacobian has already been considered (TODO: replace with integer linked list)
   integer                      :: usedctr                  ! counter to fill positions of used
   integer                      :: num_of_neighbor_levels   ! Number of neighboring levels that affect one element's column of the Jacobian
   integer                      :: max_num_of_neighbors     ! Maximum number of neighboring elements that affect one element's column of the Jacobian
   logical                      :: withMPI
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine NumJacobian_Construct(this, mesh, nEqn, controlVariables)
      implicit none
      !-arguments-----------------------------------------
      class(NumJacobian_t) , intent(inout) :: this
      type(HexMesh)        , intent(inout) :: mesh
      integer              , intent(in)    :: nEqn
      type(FTValueDictionary)  , intent(in)  :: controlVariables
      !---------------------------------------------------
      
!
!     Construct parent
!     ----------------
      call this % JacobianComputer_t % construct (mesh, nEqn, controlVariables)

      call SetNoNeighbours(this, controlVariables)
      
      call Stopwatch % CreateNewEvent("Numerical Jacobian construction")
      
      ! Big TODO: Move everything that is inside "if (isfirst)" to here!!
      
   end subroutine NumJacobian_Construct

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine NumJacobian_Compute(this, sem, nEqn, time, Matrix, TimeDerivative, TimeDerivativeIsolated, eps_in, BlockDiagonalized, mode )
      !-------------------------------------------------------------------
      class(NumJacobian_t)      , intent(inout)          :: this
      type(DGSem),                intent(inout)          :: sem
      integer,                    intent(in)             :: nEqn
      real(kind=RP),              intent(in)             :: time
      class(Matrix_t)          ,  intent(inout)          :: Matrix
      procedure(ComputeTimeDerivative_f), optional       :: TimeDerivative
      procedure(ComputeTimeDerivative_f), optional       :: TimeDerivativeIsolated
      real(kind=RP),   optional, intent(in)              :: eps_in
      logical, optional, intent(in)        :: BlockDiagonalized  !<? Construct only the block diagonal? (Only for AnJacobian_t)
      integer, optional, intent(in)                      :: mode
      !-------------------------------------------------------------------
      integer                                            :: nelm
      integer                                            :: thiscolor, thiselmidx, thiselm         ! specific counters
      integer                                            :: thisdof                           ! specific counters
      integer                                            :: ielm, felm                      
      integer, save                                      :: nnz, totalnnz
      integer                           , save           :: maxndofel
      integer, allocatable, dimension(:), save           :: Nx, Ny, Nz                             ! Polynomial orders
      integer, allocatable, dimension(:), save           :: ndofcol                                ! Maximum number of degrees of freedom in each color        
      integer, allocatable                               :: cols(:)
      integer, allocatable                               :: rows(:)
      integer, allocatable                               :: diag(:)
      real(kind=RP), allocatable, save                   :: Q0(:)
      real(kind=RP), allocatable, save                   :: QDot0(:)
      
      integer :: i, j, ii, jj, kk, eID ! General counters
      integer, dimension(4)                              :: ijkl                                   ! Indexes to locate certain degree of freedom i,j,k...l:equation number
      real(kind=RP), save                                :: eps                                    ! Perturbation magnitude
      
      logical, save                                      :: isfirst = .TRUE.
#if (!defined(NAVIERSTOKES))
      logical                                            :: computeGradients = .true.
#endif
      integer                         :: rank, ierror, mpisize
      integer, parameter :: faces_and_one = 7 ! hard coded number of faces for HEXA + 1
      integer :: Gloabl_nelm, thiselm_g
      real(kind=RP), pointer :: pbuffer(:)
      integer, allocatable, dimension(:) :: counts_recv, displacements
      type(TProgressBar) :: progress_bar
      integer, allocatable, dimension(:) :: el_reordering, el_reordering_idx
      !-------------------------------------------------------------------
      
      if(.not. present(TimeDerivative) ) error stop 'NumJacobian_Compute needs the time-derivative procedure'
!
!     --------------------------------------------------------------------
!     Initialize variables that will be used throughout all the simulation
!     --------------------------------------------------------------------
!
      
      call Stopwatch % Start("Numerical Jacobian construction")

#ifdef _HAS_MPI_
      withMPI = .true.
#else
      withMPI = .false.
#endif

      if (.NOT. isfirst) then 
         deallocate(nbr)
         deallocate(Nx)
         deallocate(Ny)
         deallocate(Nz)
         deallocate(used)
         deallocate(ndofcol)
         deallocate(QDot0)
         deallocate(Q0)
         deallocate(nbr_g)
      end if

      nelm = size(sem % mesh % elements)
      Gloabl_nelm = sem % mesh % no_of_allElements

!
!     Initialize the colorings structure
!     ----------------------------------
      allocate(nbr(nelm))
      CALL Look_for_neighbour(nbr, sem % mesh)
      allocate(nbr_g(Gloabl_nelm))

#ifdef _HAS_MPI_
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)
      call mpi_comm_size(MPI_COMM_WORLD, mpisize, ierror)
      allocate(counts_recv(mpisize))
      allocate(displacements(mpisize))

      call mpi_allgather(nelm, 1, MPI_INT, counts_recv, 1, MPI_INT, MPI_COMM_WORLD, ierror)

      counts_recv = counts_recv * faces_and_one
      displacements(1) = 0
      if (mpisize .ge. 2) then
         do i = 2, mpisize
            displacements(i) = displacements(i-1) + counts_recv(i-1)
         end do
      end if

      call mpi_allgatherv(nbr, nelm*faces_and_one, MPI_INT, nbr_g, counts_recv, displacements, MPI_INT, MPI_COMM_WORLD, ierror)

      deallocate(counts_recv)
      deallocate(displacements)
#else
      nbr_g = nbr
#endif

!
!     Re-order neighbours
!     -------------------
      allocate(el_reordering(Gloabl_nelm))
      allocate(el_reordering_idx(Gloabl_nelm))

      do i = 1, Gloabl_nelm
         el_reordering(i) = nbr_g(i) % elmnt(7)
      end do

      do i = 1, Gloabl_nelm
         el_reordering_idx(i) = my_findloc(el_reordering, i, 1)
      end do

      nbr_g = nbr_g(el_reordering_idx)

      deallocate(el_reordering)
      deallocate(el_reordering_idx)

!
!     Assemble colors
!     ---------------

      call ecolors % construct(nbr_g, num_of_neighbor_levels)       
!
!     Allocate storage
!     ----------------
      allocate(Nx(nelm), Ny(nelm), Nz(nelm))
      
      do i = 1, nelm
         Nx(i) = sem % mesh % elements(i) % Nxyz(1)
         Ny(i) = sem % mesh % elements(i) % Nxyz(2)
         Nz(i) = sem % mesh % elements(i) % Nxyz(3)
      end do         

      maxndofel = MAXVAL(this % ndofelm)                                             ! TODO: if there's p-adaptation, this value has to be recomputed         
!
!     ---------------------------------------------------------------------------------
!     Get the maximum number of neighbors ["of neighbors" * (num_of_neighbor_levels-1)] 
!     that are needed for the Jacobian computation (mesh dependent)
!     ---------------------------------------------------------------------------------
!
      max_num_of_neighbors = 0 ! Initialize to minimum possible value
      if (.not. withMPI) then
         do i=1, nelm
            max_num_of_neighbors = max (getNumOfNeighbors (i, num_of_neighbor_levels), max_num_of_neighbors)
         end do
      end if
      
!
!     ---------------------------------------------------------------
!     Allocate the used array that will contain the information about
!     which neighbor elements were already used in the numerical
!     computation of the Jacobian matrix entries
!     -> The neighbors (including itself) and a last entry that will be 0 always (boundary index)
!     ---------------------------------------------------------------
!
      allocate ( used(max_num_of_neighbors+1) )
!
!     -------------------------------------------------------------------------
!     Set max number of nonzero values expected in a row of the Jacobian matrix    TODO: if there's p-adaptation, this has to be recomputed
!           Assumes Legendre-Gauss quadrature and neglects zero values in each 
!              block (taken into account later when assembling)
!           For Legendre-Gauss-Lobatto much less entries are expected (a node on the
!              interface has more cols than an interior node)
!           IMPORTANT: These numbers assume conforming meshes!
!     -------------------------------------------------------------------------
!
      nnz = maxndofel * max_num_of_neighbors
!
!     --------------------------------------------------------------
!     Compute the maximum number of degrees of freedom in each color               TODO: if there's p-adaptation, this has to be recomputed
!     --------------------------------------------------------------
!    
      allocate(ndofcol(ecolors % num_of_colors))
      ndofcol = 0
      DO thiscolor = 1 , ecolors % num_of_colors
         ielm = ecolors%bounds(thiscolor)             
         felm = ecolors%bounds(thiscolor+1)
         DO thiselmidx = ielm, felm-1              !perturbs a dof in all elements within current color
            thiselm_g = ecolors%elmnts(thiselmidx) ! global eID
            thiselm = thiselm_g
            if (thiselm .gt. 0) ndofcol(thiscolor) = MAX(ndofcol(thiscolor),this % ndofelm(thiselm))
         END DO
      END DO
      
      allocate(QDot0(size(sem % mesh % storage % QDot)))
      allocate(Q0   (size(sem % mesh % storage % QDot)))
      
      ! All initializations done!
      isfirst = .FALSE.
!
!     ---------------------------------------------
!     Set value of eps (currently using Mettot et al. approach with L2 norm because it seems to work)
!        See:
!           > Mettot, Clément, Florent Renac, and Denis Sipp. "Computation of eigenvalue sensitivity to base flow modifications in a discrete framework: Application to open-loop control." Journal of Computational Physics 269 (2014): 234-258.
!           > Knoll, Dana A., and David E. Keyes. "Jacobian-free Newton–Krylov methods: a survey of approaches and applications." Journal of Computational Physics 193.2 (2004): 357-397.
!     --------------------------------------------
!
      if (present(eps_in)) then
         eps = eps_in
      else
         call sem % mesh % storage % local2GlobalQ (sem % NDOF)
         associate (Q => sem % mesh % storage % Q)
         eps = sqrt(EPSILON(eps))*(NORM2(Q)+1._RP) ! 1.e-8_RP: Sometimes gives better results
         end associate
      end if
!
!     *************************************************
!     If Jacobian matrix was not preallocated, allocate
!     *************************************************
!
      if (.not. this % preAllocate) then
         select type(Matrix_p => Matrix)
            type is(DenseBlockDiagMatrix_t)
               call Matrix_p % Preallocate(nnzs=this % ndofelm_l) ! Constructing with block size
            type is(SparseBlockDiagMatrix_t)
               call Matrix_p % Preallocate(nnzs=this % ndofelm_l) ! Constructing with block size
            type is(CSRMat_t)
!~                call GetRowsAndColsVector(sem, nbr, nEqn, Matrix_p % num_of_Rows, totalnnz, this % firstIdx, rows, cols, diag)
!~                call Matrix_p % PreAllocateWithStructure(totalnnz, rows, cols, diag) 
               call Matrix_p % Preallocate()
            class default ! Construct with nonzeros in each row
               call Matrix_p % Preallocate(nnz)
         end select
         
         call Matrix % SpecifyBlockInfo(this % firstIdx,this % ndofelm)
      end if
      
      call Matrix % Reset(ForceDiagonal = .TRUE.)
      
#if defined(CAHNHILLIARD)
      CALL TimeDerivative( sem % mesh, sem % particles, time, mode)
#else
      CALL TimeDerivative( sem % mesh, sem % particles, time, CTD_IGNORE_MODE )
#endif

#if defined(NAVIERSTOKES)
!$omp do schedule(runtime) private(ii,jj,kk)
      do eID = 1, sem % mesh % no_of_elements
         associate ( e => sem % mesh % elements(eID) )
         do kk = 0, e % Nxyz(3)   ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            e % storage % QDot(:,ii,jj,kk) = e % storage % QDot(:,ii,jj,kk) - e % storage % S_NS(:,ii,jj,kk)
         end do                  ; end do                ; end do
         end associate
      end do
!$omp end do
#elif defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
!$omp do schedule(runtime) private(ii,jj,kk)
      do eID = 1, sem % mesh % no_of_elements
         associate ( e => sem % mesh % elements(eID) )
         do kk = 0, e % Nxyz(3)   ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            e % storage % QDot(:,ii,jj,kk) = e % storage % QDot(:,ii,jj,kk) - e % storage % S_NS(:,ii,jj,kk)
         end do                  ; end do                ; end do
         end associate
      end do
!$omp end do
#endif

!
!     Save base state in Q0 and QDot0
!     -------------------------------
#if defined(CAHNHILLIARD)
      call sem % mesh % SetStorageToEqn(C_BC)
#endif
      
      call sem % mesh % storage % local2GlobalQdot (sem % NDOF)
      call sem % mesh % storage % local2GlobalQ    (sem % NDOF)
      QDot0 = sem % mesh % storage % QDot
      Q0    = sem % mesh % storage % Q
!
!     ------------------------------------------
!     Compute numerical Jacobian using colorings
!     ------------------------------------------
!
      call progress_bar % initialize(" Constructing numerical Jacobian...")
      if (withMPI) then

         select type(Matrix_p => Matrix)
            type is(DenseBlockDiagMatrix_t)
               ! all good
            type is(CSRMat_t)
               error stop "NumericalJacobian :: Full CSR Jacobian computation not compatible with MPI."
            class default
               error stop "NumericalJacobian :: Unknown matrix type."
         end select

!        Go through every color to obtain its elements' contribution to the Jacobian
!        ***************************************************************************
         do thiscolor = 1 , ecolors % num_of_colors
            if (this % verbose) write(STD_OUT,'(10X,A,1I6,A,1I6,A)') "Numerical Jacobian computing ", thiscolor , " out of ", ecolors % num_of_colors, " colors."

            ielm = ecolors%bounds(thiscolor)             ! Initial element of the color
            felm = ecolors%bounds(thiscolor+1)           ! Final element of the color! 

!           Iterate through the DOFs in thiscolor
!              ( Computes one column for each dof within an element )
!           ********************************************************
            do thisdof = 1, ndofcol(thiscolor)

!              Perturb the current degree of freedom in all elements within current color
!              --------------------------------------------------------------------------
               do thiselmidx = ielm, felm-1              
                  thiselm_g = ecolors%elmnts(thiselmidx) ! global eID
                  thiselm = mpi_partition % global2localeID(thiselm_g) ! local eID

                  if (thiselm .gt. 0) then
                  if (this % ndofelm(thiselm)<thisdof) cycle    ! Do nothing if the DOF exceeds the NDOF of thiselm

                     ijkl = local2ijk(thisdof,nEqn,Nx(thiselm),Ny(thiselm),Nz(thiselm))
                     sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                         sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 

                  end if
               end do ! thiselmidx = ielm, felm-1 
!  
!              Compute the time derivative
!              ---------------------------
#if defined(CAHNHILLIARD)
               CALL TimeDerivative( sem % mesh, sem % particles, time, mode)
#else
               CALL TimeDerivative( sem % mesh, sem % particles, time, CTD_IGNORE_MODE )
#endif

#if defined(NAVIERSTOKES)
!$omp do schedule(runtime) private(ii,jj,kk)
         do eID = 1, sem % mesh % no_of_elements
            associate ( e => sem % mesh % elements(eID) )
            do kk = 0, e % Nxyz(3)   ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
               e % storage % QDot(:,ii,jj,kk) = e % storage % QDot(:,ii,jj,kk) - e % storage % S_NS(:,ii,jj,kk)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
#elif defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
!$omp do schedule(runtime) private(ii,jj,kk)
         do eID = 1, sem % mesh % no_of_elements
            associate ( e => sem % mesh % elements(eID) )
            do kk = 0, e % Nxyz(3)   ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
               e % storage % QDot(:,ii,jj,kk) = e % storage % QDot(:,ii,jj,kk) - e % storage % S_NS(:,ii,jj,kk)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
#endif

               call sem % mesh % storage % local2GlobalQdot (sem %NDOF)
               sem % mesh % storage % QDot = (sem % mesh % storage % QDot - QDot0) / eps
               call sem % mesh % storage % global2LocalQdot
               
   !
   !           Add the contributions to the Jacobian
   !           -------------------------------------
               do thiselmidx = ielm, felm-1
                  thiselm_g = ecolors%elmnts(thiselmidx) ! global eID
                  thiselm = mpi_partition % global2localeID(thiselm_g) ! local eID

                  if (thiselm .gt. 0) then

                     IF (this % ndofelm(thiselm)<thisdof) CYCLE

                     pbuffer(1:this % ndofelm(thiselm)) => sem % mesh % storage % elements(thiselm) % QDot 

                     do j=1, this % ndofelm(thiselm)
                        call Matrix % AddToBlockEntry (sem % mesh % elements(thiselm) % GlobID, sem % mesh % elements(thiselm) % GlobID, &
                           j, thisdof, pbuffer(j) )
                     end do
                     
                  end if

               end do ! thiselmidx = ielm, felm-1      
   !
   !           Restore original values for Q (TODO: this can be improved)
   !           ----------------------------------------------------------
               sem % mesh % storage % Q = Q0
               call sem % mesh % storage % global2LocalQ

            ENDDO ! thisdof = 1, ndofcol(thiscolor)
         ENDDO ! thiscolor = 1 , ecolors % num_of_colors

      else ! NOMPI

!     Go through every color to obtain its elements' contribution to the Jacobian
!     ***************************************************************************
      do thiscolor = 1 , ecolors % num_of_colors
         ! if (this % verbose) call progress_bar % run(real(100 * thiscolor / ecolors % num_of_colors),5," Constructing numerical Jacobian...")
         if (this % verbose) write(STD_OUT,'(10X,A,1I6,A,1I6,A)') "Numerical Jacobian computing ", thiscolor , " out of ", ecolors % num_of_colors, " colors."
         ielm = ecolors%bounds(thiscolor)             ! Initial element of the color
         felm = ecolors%bounds(thiscolor+1)           ! Final element of the color
!         
!        Iterate through the DOFs in thiscolor
!           ( Computes one column for each dof within an element )
!        ********************************************************
         do thisdof = 1, ndofcol(thiscolor)
            
!           Perturb the current degree of freedom in all elements within current color
!           --------------------------------------------------------------------------
            DO thiselmidx = ielm, felm-1              
               thiselm = ecolors%elmnts(thiselmidx)
               IF (this % ndofelm(thiselm)<thisdof) CYCLE    ! Do nothing if the DOF exceeds the NDOF of thiselm
               
               ijkl = local2ijk(thisdof,nEqn,Nx(thiselm),Ny(thiselm),Nz(thiselm))
               
               sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                   sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 
            ENDDO
!
!           Compute the time derivative
!           ---------------------------
#if defined(CAHNHILLIARD)
            CALL TimeDerivative( sem % mesh, sem % particles, time, mode)
#else
            CALL TimeDerivative( sem % mesh, sem % particles, time, CTD_IGNORE_MODE )
#endif

#if defined(NAVIERSTOKES)
!$omp do schedule(runtime) private(ii,jj,kk)
      do eID = 1, sem % mesh % no_of_elements
         associate ( e => sem % mesh % elements(eID) )
         do kk = 0, e % Nxyz(3)   ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            e % storage % QDot(:,ii,jj,kk) = e % storage % QDot(:,ii,jj,kk) - e % storage % S_NS(:,ii,jj,kk)
         end do                  ; end do                ; end do
         end associate
      end do
!$omp end do
#elif defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
!$omp do schedule(runtime) private(ii,jj,kk)
      do eID = 1, sem % mesh % no_of_elements
         associate ( e => sem % mesh % elements(eID) )
         do kk = 0, e % Nxyz(3)   ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            e % storage % QDot(:,ii,jj,kk) = e % storage % QDot(:,ii,jj,kk) - e % storage % S_NS(:,ii,jj,kk)
         end do                  ; end do                ; end do
         end associate
      end do
!$omp end do
#endif

            call sem % mesh % storage % local2GlobalQdot (sem %NDOF)
            sem % mesh % storage % QDot = (sem % mesh % storage % QDot - QDot0) / eps
            call sem % mesh % storage % global2LocalQdot
            
!
!           Add the contributions to the Jacobian
!           -------------------------------------
            do thiselmidx = ielm, felm-1
               thiselm = ecolors%elmnts(thiselmidx)
               IF (this % ndofelm(thiselm)<thisdof) CYCLE
               ! Redefine used array and counter
               used    = 0
               usedctr = 1
               
               call this % AssignColToJacobian(Matrix, sem % mesh, thiselm, thiselm, thisdof, num_of_neighbor_levels)
               
            END DO           
!
!           Restore original values for Q (TODO: this can be improved)
!           ----------------------------------------------------------
            sem % mesh % storage % Q = Q0
            call sem % mesh % storage % global2LocalQ
         ENDDO
      ENDDO

      end if ! MPI

      CALL Matrix % Assembly()                             ! Matrix A needs to be assembled before being used
      
      call Stopwatch % Pause("Numerical Jacobian construction")
      
      IF (this % verbose) PRINT*, "Numerical Jacobian construction: ", Stopwatch % lastElapsedTime("Numerical Jacobian construction"), "seconds"

!
!     --------------------
!     Return storage to NS
!     --------------------
!
#if defined(FLOW) && defined(CAHNHILLIARD)
      call sem % mesh % SetStorageToEqn(NS_BC)
#endif
      
      ! call Matrix % Visualize('Jacobian.txt')
      ! error stop "TBC"
   end subroutine NumJacobian_Compute
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------------------------
!  Assign to Jacobian the column that corresponds to the element "eID" and the degree of freedom "thisdof" 
!  taking into account the contribution of the neighbors [(depth-1) * "of neighbors"]
!  -------------------------------------------------------------------------------------------------------
   recursive subroutine NumJacobian_AssignColToJacobian(this, Matrix, mesh, eID, eIDn, thisdof, depth)
      implicit none
      !-arguments---------------------------------------
      class(NumJacobian_t)           , intent(inout) :: this
      class(Matrix_t)                , intent(inout) :: Matrix  !<> Jacobian Matrix
      type(HexMesh)                  , intent(in)    :: mesh
      integer                        , intent(in)    :: eID     !<  Element ID 
      integer                        , intent(in)    :: eIDn    !<  ID of the element, whose neighbors' contributions are added
      integer                        , intent(in)    :: thisdof !<  Current degree of freedom
      integer                        , intent(in)    :: depth   !<  Amount of neighbors to visit
      !-local-variables---------------------------------
      integer :: elmnbr                    ! Neighbor element index
      integer :: i,j                       ! Counter
      integer :: ndof                      ! Number of degrees of freedom of element
      real(kind=RP), pointer :: pbuffer(:) ! Buffer to point to an element's Qdot
      !-------------------------------------------------
      
      if ( (eID  == 0) .or. (eIDn == 0) ) return

!
!     Go through all the neighbors
!     ----------------------------
      do i = 1, size(nbr(eIDn) % elmnt)
         elmnbr = nbr(eIDn) % elmnt(i) 
      
         if (.NOT. any(used == elmnbr)) THEN
            ndof   = this % ndofelm(elmnbr)
            pbuffer(1:ndof) => mesh % storage % elements(elmnbr) % QDot  !maps Qdot array into a 1D pointer
            
            do j=1, this % ndofelm(elmnbr)
               call Matrix % AddToBlockEntry (mesh % elements(elmnbr) % GlobID, mesh % elements(eID) % GlobID, j, thisdof, pbuffer(j) )
            end do
            
            used(usedctr) = elmnbr
            usedctr = usedctr + 1
         end if
         
         if (depth > 1) call this % AssignColToJacobian(Matrix, mesh, eID, elmnbr, thisdof, depth-1)
         
      end do
      
   end subroutine NumJacobian_AssignColToJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------------
!  Returns the number of neighbors [(depth-1) * "of neighbors"] for a specific element (counting itself)
!  -----------------------------------------------------------------------------------------------------
   function getNumOfNeighbors (eID, depth) result(num)
      implicit none
      !-arguments---------------------------------------
      integer                      , intent(in)              :: eID     !<  Element ID 
      integer                      , intent(in)              :: depth   !<  Amount of neighbors to visit
      integer           :: num
      !-local-variables---------------------------------
      type(IntegerDataLinkedList_t) :: neighborsList
      !-------------------------------------------------
      
      num = 0
      if (eID == 0) return 
      
!
!     Create list of already counted elements
!     ---------------------------------------
      neighborsList = IntegerDataLinkedList_t(.FALSE.)
      
!
!     Add neighbors to list
!     ---------------------
      call addNeighborsToList(eID,depth,neighborsList)
      num = neighborsList % no_of_entries
      
!
!     destruct list of already counted elements
!     -----------------------------------------
      call neighborsList % destruct
      
   end function getNumOfNeighbors
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   recursive subroutine addNeighborsToList(eID, depth, neighborsList)
      implicit none
      !-arguments---------------------------------------
      integer                      , intent(in)    :: eID     !<  Element ID 
      integer                      , intent(in)    :: depth   !<  Amount of neighbors to visit
      type(IntegerDataLinkedList_t), intent(inout) :: neighborsList
      !-local-variables---------------------------------
      integer :: elmnbr                    ! Neighbor element index
      integer :: i                         ! Counter
      !-------------------------------------------------
      
      do i = 1, NUM_OF_NEIGHBORS + 1
         elmnbr = nbr(eID) % elmnt(i)
         
         if (elmnbr == 0) cycle
            
         call neighborsList % add(elmnbr)
         if (depth > 1) call addNeighborsToList (elmnbr, depth - 1, neighborsList)
         
      end do
   end subroutine addNeighborsToList
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
!  ! To be deprecated(?)
   subroutine GetRowsAndColsVector(sem, nbr, nEqn, nRows, nnz, firstIdx, rows, cols, diag)
      implicit none
      class(DGSEM)         :: sem
      type(Neighbor_t)     :: nbr(:)
      integer, intent(in)  :: nEqn, nRows
      integer, intent(in)  :: firstIdx(1:sem % mesh % no_of_elements)
      integer, intent(out)  :: nnz
      integer, allocatable, intent(out) :: rows(:)
      integer, allocatable, intent(out) :: cols(:)
      integer, allocatable, intent(out) :: diag(:)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                 :: eID, i, j, k, counter, csr_pos, ieq, nID
      integer                 :: ii, jj, kk, iieq
      integer                 :: pos_i, pos_j
      integer                 :: lb, ub
      integer, dimension(26)  :: neighbours
  
!
!     *********************
!     First loop to get nnz: TODO could I already set rows here?
!     *********************
!
!!$omp parallel private(counter, neighbours, i, j)
!!$omp do reduction(+:nnz)
      nnz = 0
      do eID = 1, sem % mesh % no_of_elements
         associate(e => sem % mesh % elements(eID))
!
!     1/ For each element, get its neighbours
!        ------------------------------------
         counter = 0
         neighbours = -1

         do i = 1, size(nbr(eID) % elmnt)
            if ( nbr(eID) % elmnt(i) .le. 0 ) cycle
            do j = 1, size(nbr(nbr(eID) % elmnt(i)) % elmnt)
               if (nbr(nbr(eID) % elmnt(i)) % elmnt(j) .le. 0) cycle

               if ( .not. any(neighbours .eq. nbr(nbr(eID) % elmnt(i)) % elmnt(j)) ) then
                  counter = counter + 1 
                  neighbours(counter) = nbr(nbr(eID) % elmnt(i)) % elmnt(j)
               end if
            end do
         end do

         do i = 1, counter
            associate(eL => sem % mesh % elements(neighbours(i)))
            nnz = nnz + nEqn*(eL % Nxyz(1)+1)*(eL % Nxyz(2)+1)*(eL % Nxyz(3)+1)*nEqn*(e % Nxyz(1)+1)*(e % Nxyz(2)+1)*(e % Nxyz(3)+1)
            end associate
         end do
         end associate
      end do
!!$omp end do
!!$omp end parallel
      allocate(rows(1:nRows+1))
      allocate(cols(1:nnz))
      allocate(diag(1:nRows))
!
!     ****************************
!     We need to set rows and cols
!     ****************************
!
      csr_pos = 1
      do eID = 1, sem % mesh % no_of_elements
         associate(e => sem % mesh % elements(eID))
!
!     1/ For each element, get its neighbours
!        ------------------------------------
         counter = 0
         neighbours = -1

         do i = 1, size(nbr(eID) % elmnt)
            if ( nbr(eID) % elmnt(i) .le. 0 ) cycle
            do j = 1, size(nbr(nbr(eID) % elmnt(i)) % elmnt)          ! For only neighbors comment this
               if (nbr(nbr(eID) % elmnt(i)) % elmnt(j) .le. 0) cycle  ! For only neighbors comment this

               if ( .not. any(neighbours .eq. nbr(nbr(eID) % elmnt(i)) % elmnt(j)) ) then ! For neighbors change by if ( .not. any(neighbours .eq. nbr(eID) % elmnt(i)) ) then
                  counter = counter + 1 
                  neighbours(counter) = nbr(eID) % elmnt(i)                               ! For neighbors change by neighbours(counter) = nbr(eID) % elmnt(i)
               end if
            end do                                                    ! For only neighbors comment this
         end do
         call Qsort(neighbours(1:counter))

         pos_i = firstIdx(eID)

         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) ; do ieq = 1, nEqn
            rows(pos_i) = csr_pos
            pos_i = pos_i + 1 
            do nID = 1, counter
               associate(neigh => sem % mesh % elements(neighbours(nID)))

               pos_j = firstIdx(neighbours(nID))

               do kk = 0, neigh % Nxyz(3) ; do jj = 0, neigh % Nxyz(2) ; do ii = 0, neigh % Nxyz(1) ; do iieq = 1, nEqn
                    
                  cols(csr_pos) = pos_j

                  pos_j = pos_j + 1
                  csr_pos = csr_pos + 1 
               end do                  ; end do                ; end do                 ; end do
               end associate
            end do   
         end do                ; end do                ; end do                ; end do

         end associate
      end do

      rows(nRows+1) = nnz+1
!
!     **************************
!     Get the diagonal positions
!     **************************
!
      do i = 1, nRows
         lb = rows(i)
         ub = rows(i+1)-1

         pos_i = minloc(abs(cols(lb:ub)-i),dim=1)
         diag(i) = lb + pos_i - 1
      end do

   end subroutine GetRowsAndColsVector
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------
!     Destruct the JacobianInfo variables
!     ------------------------------------
   subroutine SetNoNeighbours(this, controlVariables)
      implicit none
      !-arguments-----------------------------------------
      class(NumJacobian_t), intent(inout)   :: this
      type(FTValueDictionary)  , intent(in) :: controlVariables
      !---------------------------------------------------
      character(len=LINE_LENGTH)                     :: tmpc
      !---------------------------------------------------

#if defined(CAHNHILLIARD)
      num_of_neighbor_levels = 4
#elif defined(NAVIERSTOKES)
      if (flowIsNavierStokes) then
         if (controlVariables % containsKey("viscous discretization")) then

            tmpc = controlVariables % StringValueForKey("viscous discretization",LINE_LENGTH)
            select case (tmpc)
            case('BR1')
               num_of_neighbor_levels = 2
            case('BR2')
               num_of_neighbor_levels = 1
            case('IP')
               num_of_neighbor_levels = 1
            case default 
               if (MPI_Process % isRoot) error stop 'JacobianComputerClass :: Viscous discretization not recognized.'
            end select
         else
            if (MPI_Process % isRoot) write(STD_OUT,*) 'JacobianComputerClass :: Viscous discretization not defined. Jacobian assumes BR1.'
            num_of_neighbor_levels = 2
         end if 
      else
         num_of_neighbor_levels = 1
      end if
#else
      num_of_neighbor_levels = 2
#endif

   end subroutine SetNoNeighbours
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module NumericalJacobian