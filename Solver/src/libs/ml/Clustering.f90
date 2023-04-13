!**************************************************************************************************!
!
! Clustering algorithms implemented:
!
!   * k-means (kMeans_t)
!   * Gaussian Mixture Model (GMM_t)
!
! The algorithms are abstracted as classes with similar APIs. An instance of any of these types
! always has the attributes:
!
!   * ndims: dimension of the feature space
!   * nclusters: number of clusters
!
! It also has the following methods:
!
!   * init: initialize algorithm
!   * fit: fit clusters to data
!   * predict: cluster nodes
!
!**************************************************************************************************!
!   k-means
!**************************************************************************************************!
!
! Attributes holding the results of fit and predict methods:
!
!   * centroids(ndims, nclusters)
!   * clusters(npts)
!
! Methods:
!
!   * init(ndims, nclusters, maxiters)
!
!      [IN]
!      + ndims
!         dimension of the feature space
!      + nclusters
!         number of clusters
!      + maxiters [optional] = 100
!         maximum number of iterations
!
!   * fit(x, info, centroids, reset)
!
!      [IN]
!      + x(ndims, npts)
!         points to cluster
!      + centroids(ndims, nclusters) [optional]
!         starting centroids of the clusters
!      + reset [optional] = .false.
!         only used when no centroids are given. Restart centroids with random values when .true.
!         and use centroids from the previous step otherwise
!
!      [OUT]
!      + info [optional]
!         status of the algorithm: number of iterations or -1 if no convergence
!
!   * predict(x, info)
!
!      [IN]
!      + x(ndims, npts)
!         points to cluster
!      + info [optional]
!         status of the algorithm: <0 if the computation failed
!
! Note that calling fit also updates the clusters, so predict must be called only when the points x
! are different from the ones used for fitting.
!
!**************************************************************************************************!
!   GMM
!**************************************************************************************************!
!
! This adapted version of the GMM implements an optional adaptation step to remove overlapping
! clusters. When using this, the number of clusters will be <= nclusters, and the allowed maximum
! number of clusters is stored in max_nclusters.
!
! Attributes holding the results of fit and predict methods:
!
!   * centroids(ndims, max_nclusters)
!   * prob(npts, max_nclusters)
!
! Methods:
!
!   * init(ndims, nclusters, maxiters, collapse_tol, logl_tol, zero_tol)
!
!      [IN]
!      + ndims
!         dimension of the feature space
!      + nclusters
!         number of clusters
!      + maxiters [optional] = 100
!         maximum number of iterations
!      + collapse_tol [optional] = 2e-5
!         maximum distance (non-dimensional) between the centroids of overlapping clusters
!      + logl_tol [optional] = 1e-3
!         relative tolerance between succesive computations of log L to assume convergence
!      + zero_tol [optional] = 1e-6
!         regularization value added to the diagonal of the covariance matrices
!
!   * fit(x, info, centroids, reset, adapt, from_kmeans)
!
!      [IN]
!      + x(ndims, npts)
!         points to cluster
!      + centroids(ndims, nclusters) [optional]
!         starting centroids of deleted clusters (all if restart = .true.)
!      + reset [optional] = .false.
!         restart all centroids if .true.
!      + adapt [optional] = .false.
!         remove overlapping clusters when .true. Deleted clusters are filled with centroids or
!         random values
!      + from_kmeans [optional] = .false.
!         if .true. update the centroids with 10 iterations of k-means after setting their values
!         according to the other options
!
!      [OUT]
!      + info [optional]
!         status of the algorithm: number of iterations, or -1 if dimensions do not match,
!         or -2 if no convergence after maxiters of GMM and maxiters of k-means
!
!   * predict(x, info)
!
!      [IN]
!      + x(ndims, npts)
!         points to cluster
!      + info [optional]
!         status of the algorithm: <0 if the computation failed
!
! Note that fit does not compute the probabilities, so predict must always be called before using
! them. Also note that the GMM only provide probabilities, not clusters. These can be computed by,
! for example, taking the cluster with the maximum probability for each point.
!
!**************************************************************************************************!
module Clustering

   use SMConstants,       only: RP, PI
   use Utilities,         only: AlmostEqual, BubblesortWithFriend
   use MPI_Process_Info,  only: MPI_Process
   use MPI_Utilities,     only: MPI_SumAll, MPI_MinMax
#if defined(_HAS_MPI_)
   use mpi
#endif

   implicit none

   private

   public :: kMeans_t
   public :: GMM_t
   public :: standardize
   public :: rescale

   type :: kMeans_t
      logical               :: initialized = .false.
      integer               :: ndims
      integer               :: nclusters
      integer               :: maxiters
      logical               :: centroids_set
      integer,  allocatable :: prevClusters(:)  ! [p]: clusters at the previous iteration
      integer,  allocatable :: clusters(:)      ! [p]: clusters
      real(RP), allocatable :: centroids(:,:)   ! [d,n]: centroids of the clusters
   contains
      procedure :: init    => kMeans_init
      procedure :: fit     => kMeans_fit
      procedure :: predict => kMeans_predict
      final     :: kMeans_final
   end type kMeans_t

   type :: GaussianList_t
      real(RP), pointer, contiguous :: storage(:)    => null()  ! contiguous storage for MPI communication
      real(RP), pointer, contiguous :: logtau(:)     => null()  ! [n]: log of weights
      real(RP), pointer, contiguous :: mu(:,:)       => null()  ! [d,n]: centroids
      real(RP), pointer, contiguous :: cov(:,:,:)    => null()  ! [d,d,n]: covariance matrices
      real(RP), pointer, contiguous :: covinv(:,:,:) => null()  ! [d,d,n]: inverse covariance matrices
      real(RP), pointer, contiguous :: logdet(:)     => null()  ! [n]: log-determinant of covariance matrices
      real(RP), allocatable         :: tcov(:,:,:)              ! [d,d,t]: thread-local covariance matrix
   end type GaussianList_t

   type :: GMM_t
      logical                       :: initialized = .false.
      logical                       :: with_kmeans
      integer                       :: ndims
      integer                       :: nclusters
      integer                       :: max_nclusters
      integer                       :: maxiters
      real(RP)                      :: mutol           ! Relative tolerance to collapse clusters
      real(RP)                      :: lltol           ! Relative change in log(L) to stop iterations
      real(RP)                      :: zerotol         ! Almost-zero tolerance
      real(RP)                      :: logL
      real(RP), allocatable         :: prob(:,:)       ! [p,n]: probability matrix
      real(RP), pointer, contiguous :: centroids(:,:)  ! [d,n]: centroids of the clusters
      type(GaussianList_t)          :: g
      type(kMeans_t)                :: kmeans
   contains
      procedure :: init    => GMM_init
      procedure :: fit     => GMM_fit
      procedure :: predict => GMM_predict
      final     :: GMM_final
   end type GMM_t

   real(RP), parameter :: LOG2PI = log(2.0_RP*PI)
!
!  ========
   contains
!  ========
!
   subroutine kMeans_init(self, ndims, nclusters, maxiters)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t),   intent(inout) :: self
      integer,           intent(in)    :: ndims
      integer,           intent(in)    :: nclusters
      integer, optional, intent(in)    :: maxiters
!
!     ---------------
!     Local variables
!     ---------------
      logical :: needs_alloc


      if (present(maxiters)) then
         self % maxiters = maxiters
      else
         self % maxiters = 100
      end if

      self % ndims     = ndims
      self % nclusters = nclusters

      needs_alloc = .true.
      if (allocated(self % centroids)) then
         if (size(self % centroids, 1) == ndims .and. size(self % centroids, 2) == nclusters) then
            needs_alloc = .false.
         else
            deallocate(self % centroids)
         end if
      end if

      if (needs_alloc) then
         allocate(self % centroids(ndims, nclusters))
      end if

      self % initialized   = .true.
      self % centroids_set = .false.

   end subroutine kMeans_init
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_reset(self)
!
!     ---------
!     Interface
!     ---------
      type(kMeans_t), intent(inout) :: self
!
!     ---------------
!     Local variables
!     ---------------
      integer :: ierr


      if (self % initialized) then

         if (MPI_Process % isRoot) then
            call random_number(self % centroids)
         end if

#if defined(_HAS_MPI_)
         call MPI_Bcast(self % centroids, size(self % centroids), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
#endif
         self % centroids_set = .true.

      end if

   end subroutine kMeans_reset
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_final(self)
!
!     ---------
!     Interface
!     ---------
      type(kMeans_t), intent(inout) :: self


      ! Required for gfortran?
      if (self % initialized) then
         if (allocated(self % prevClusters)) deallocate(self % prevClusters)
         if (allocated(self % clusters))     deallocate(self % clusters)
         if (allocated(self % centroids))    deallocate(self % centroids)
      end if

      self % initialized   = .false.
      self % centroids_set = .false.

   end subroutine kMeans_final
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_fit(self, x, info, centroids, reset)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t),    intent(inout) :: self
      real(RP),           intent(in)    :: x(:,:)
      integer,  optional, intent(out)   :: info
      real(RP), optional, intent(in)    :: centroids(self % ndims, self % nclusters)
      logical,  optional, intent(in)    :: reset
!
!     ---------------
!     Local variables
!     ---------------
      integer :: ndims, npts, nclusters
      integer :: i, ierr
      logical :: need_alloc
      logical :: reset_centroids
      logical :: breakFlag

!
!     Checks
!     ------
      if (.not. self % initialized .or. size(x, dim=1) /= self % ndims) then
         if (present(info)) info = -1
         return
      end if

      if (present(centroids)) then
         if (size(centroids, dim=1) /= self % ndims .or. &
               size(centroids, dim=2) /= self % nclusters) then
            if (present(info)) info = -1
            return
         end if
      end if
!
!     Initial clusters
!     ----------------
      if (present(centroids)) then
         self % centroids     = centroids
         self % centroids_set = .true.
         reset_centroids      = .false.
      elseif (present(reset)) then
         reset_centroids = reset
      else
         reset_centroids = .not. self % centroids_set
      end if

      if (reset_centroids) then
         call kMeans_reset(self)
      end if
!
!     Initialization
!     --------------
      npts = size(x, dim=2)
      ndims = self % ndims
      nclusters = self % nclusters

      ! This is a bit unsafe since we only check `clusters`, but should be fine as long as we
      ! always use `clusters` and `prevClusters` together
      need_alloc = .false.
      if (allocated(self % clusters)) then
         if (size(self % clusters, dim=1) /= npts) then
            deallocate(self % clusters)
            deallocate(self % prevClusters)
            need_alloc = .true.
         end if
      else
         need_alloc = .true.
      end if

      if (need_alloc) then
         allocate(self % clusters(npts))
         allocate(self % prevClusters(npts))
      end if

      call kMeans_compute_clusters(nclusters, ndims, npts, x, self % centroids, self % clusters)
!
!     Loop until convergence
!     ----------------------
      do i = 1, self % maxiters

         self % prevClusters = self % clusters
         call kMeans_compute_centroids(nclusters, ndims, npts, x, self % centroids, self % clusters)
         call kMeans_compute_clusters(nclusters, ndims, npts, x, self % centroids, self % clusters)

         breakFlag = all(self % prevClusters == self % clusters)
#if defined(_HAS_MPI_)
         call MPI_AllReduce(MPI_IN_PLACE, breakFlag, 1, MPI_LOGICAL, MPI_LAND, &
                            MPI_COMM_WORLD, ierr)
#endif
         if (breakFlag) exit

      end do
!
!     Check convergence
!     -----------------
      if (present(info)) info = merge(-1, i, i > self % maxiters)
!
!     Sort clusters
!     -------------
      call kMeans_sort_clusters(nclusters, ndims, npts, self % centroids, self % clusters)

   end subroutine kMeans_fit
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_predict(self, x, info)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t),   intent(inout) :: self
      real(RP),          intent(in)    :: x(:,:)  ! [ndims, npts]
      integer, optional, intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      integer :: ndims
      integer :: npts
      integer :: nclusters

!
!     Checks
!     ------
      ndims = size(x, dim=1)
      if (self % ndims /= ndims) then
         if (present(info)) info = -1
         return
      end if
!
!     Compute clusters and output status
!     ----------------------------------
      npts = size(x, dim=2)
      nclusters = self % nclusters
      call kMeans_compute_clusters(nclusters, ndims, npts, x, self % centroids, self % clusters)

      if (present(info)) info = 1

   end subroutine kMeans_predict
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_compute_clusters(nclusters, ndims, npts, x, xavg, clusters)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)  :: nclusters
      integer,  intent(in)  :: ndims
      integer,  intent(in)  :: npts
      real(RP), intent(in)  :: x(ndims, npts)
      real(RP), intent(in)  :: xavg(ndims, nclusters)
      integer,  intent(out) :: clusters(npts)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j
      real(RP) :: xi(ndims)
      real(RP) :: dist, minDist


!$omp parallel do default(private) firstprivate(npts, nclusters) shared(x, xavg, clusters)
      do i = 1, npts
         xi = x(:,i)
         minDist = huge(1.0_RP)
         do j = 1, nclusters
            dist = norm2(xi - xavg(:,j))
            if (dist < minDist) then
               minDist = dist
               clusters(i) = j
            end if
         end do
      end do
!$omp end parallel do

   end subroutine kMeans_compute_clusters
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_compute_centroids(nclusters, ndims, npts, x, xavg, clusters)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)  :: nclusters
      integer,  intent(in)  :: ndims
      integer,  intent(in)  :: npts
      real(RP), intent(in)  :: x(ndims, npts)
      real(RP), intent(out) :: xavg(ndims, nclusters)
      integer,  intent(in)  :: clusters(npts)
!
!     ---------------
!     Local variables
!     ---------------
      integer :: i
      integer :: ptsInCluster(nclusters)

!
!     Summation over the points
!     -------------------------
      xavg = 0.0_RP
      ptsInCluster = 0
      do i = 1, npts
         xavg(:,clusters(i)) = xavg(:,clusters(i)) + x(:,i)
         ptsInCluster(clusters(i)) = ptsInCluster(clusters(i)) + 1
      end do
!
!     Syncronization between MPI processes
!     ------------------------------------
      call MPI_SumAll(xavg)
      call MPI_SumAll(ptsInCluster)
!
!     Final computation of the average
!     --------------------------------
      do i = 1, nclusters
         if (ptsInCluster(i) == 0) cycle
         xavg(:,i) = xavg(:,i) / real(ptsInCluster(i), kind=RP)
      end do

   end subroutine kMeans_compute_centroids
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_sort_clusters(nclusters, ndims, npts, xavg, clusters)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)    :: nclusters
      integer,  intent(in)    :: ndims
      integer,  intent(in)    :: npts
      real(RP), intent(inout) :: xavg(ndims, nclusters)
      integer,  intent(inout) :: clusters(npts)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, ierr
      real(RP) :: normvec(nclusters)
      integer  :: indices(nclusters)
      integer  :: clustermap(nclusters)

!
!     Sort the clusters
!     -----------------
      call sort_clusters(ndims, nclusters, xavg, clustermap=clustermap)
!
!     Reorder clusters
!     ----------------
      do i = 1, npts
         clusters(i) = clustermap(clusters(i))
      end do

   end subroutine kMeans_sort_clusters
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_init(self, ndims, nclusters, maxiters, collapse_tol, logl_tol, zero_tol)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t),       intent(inout) :: self
      integer,            intent(in)    :: ndims
      integer,            intent(in)    :: nclusters
      integer,  optional, intent(in)    :: maxiters
      real(RP), optional, intent(in)    :: collapse_tol
      real(RP), optional, intent(in)    :: logl_tol
      real(RP), optional, intent(in)    :: zero_tol
!
!     ---------------
!     Local variables
!     ---------------
      integer :: tsize, msize, csize, dsize, ssize
      integer :: ind1, ind2
      integer :: nthreads
      logical :: needs_allocation

!
!     Tolerances
!     ----------
      if (present(maxiters)) then
         self % maxiters = maxiters
      else
         self % maxiters = 100
      end if

      if (present(collapse_tol)) then
         self % mutol = collapse_tol
      else
         self % mutol = 2e-5_RP
      end if

      if (present(logl_tol)) then
         self % lltol = logl_tol
      else
         self % lltol = 1e-3_RP
      end if

      if (present(zero_tol)) then
         self % zerotol = zero_tol
      else
         self % zerotol = 1e-6_RP
      end if
!
!     Reset all pointers
!     ------------------
      nullify(self % centroids)
      nullify(self % g % logtau)
      nullify(self % g % mu)
      nullify(self % g % cov)
      nullify(self % g % covinv)
      nullify(self % g % logdet)
!
!     Allocate contiguous storage and set pointers
!     --------------------------------------------
      tsize = nclusters
      msize = ndims * nclusters
      csize = ndims * ndims * nclusters
      dsize = nclusters
      ssize = tsize + msize + 2 * csize + dsize

      needs_allocation = .true.
      if (self % initialized) then  ! Required for gfortran?
         if (associated(self % g % storage)) then
            if (size(self % g % storage) == ssize) then
               needs_allocation = .false.
            else
               deallocate(self % g % storage)
            end if
         end if
      end if

      if (needs_allocation) then
         allocate(self % g % storage(ssize))
      end if

      nthreads = get_num_threads()
      needs_allocation = .true.
      if (self % initialized) then  ! Required for gfortran?
         if (allocated(self % g % tcov)) then
            if (all(size(self % g % tcov) == [ndims, ndims, nthreads])) then
               needs_allocation = .false.
            else
               deallocate(self % g % tcov)
            end if
         end if
      end if

      if (needs_allocation) then
         allocate(self % g % tcov(ndims, ndims, nthreads))
      end if

      ind1 = 1
      ind2 = tsize
      self % g % logtau(1:nclusters) => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + msize
      self % g % mu(1:ndims, 1:nclusters)    => self % g % storage(ind1:ind2)
      self % centroids(1:ndims, 1:nclusters) => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + csize
      self % g % cov(1:ndims, 1:ndims, 1:nclusters) => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + csize
      self % g % covinv(1:ndims, 1:ndims, 1:nclusters) => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + dsize
      self % g % logdet(1:nclusters) => self % g % storage(ind1:ind2)
!
!     Set initial values
!     ------------------
      self % ndims         = ndims
      self % nclusters     = 0
      self % max_nclusters = nclusters
      self % logL          = huge(1.0_RP)
      self % initialized   = .true.

   end subroutine GMM_init
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_reset(self, centroids, use_kmeans, x)
!
!     ---------
!     Interface
!     ---------
      type(GMM_t),        intent(inout) :: self
      real(RP), optional, intent(in)    :: centroids(self % ndims, self % max_nclusters)
      logical,  optional, intent(in)    :: use_kmeans ! reset from kmeans
      real(RP), optional, intent(in)    :: x(:,:)
!
!     ---------------
!     Local variables
!     ---------------
      integer :: i1, i2
      integer :: i, j
      integer :: ierr

!
!     Keep the last state if all the clusters are used
!     ------------------------------------------------
      if (self % nclusters /= self % max_nclusters) then

         i1 = self % nclusters + 1
         i2 = self % max_nclusters
!
!        Set default values for the new clusters
!        ---------------------------------------
         self % logL       = huge(1.0_RP)
         self % g % logtau = log(1.0_RP / self % max_nclusters)
         do i = i1, i2
            do j = 1, self % ndims
               self % g % cov(:,j,i)    = 0.0_RP
               self % g % covinv(:,j,i) = 0.0_RP
               self % g % cov(j,j,i)    = 1.0_RP
               self % g % covinv(j,j,i) = 1.0_RP
            end do
         end do
         self % g % logdet(i1:i2) = 0.0_RP
!
!        Set the centroids
!        -----------------
         if (present(centroids)) then
            self % centroids(:,i1:i2) = centroids(:,i1:i2)

         else
            if (MPI_Process % isRoot) then
               call random_number(self % centroids(:,i1:i2))
            end if

#if defined(_HAS_MPI_)
            call MPI_Bcast(self % centroids(:,i1:i2), size(self % centroids(:,i1:i2)), MPI_DOUBLE, &
                           0, MPI_COMM_WORLD, ierr)
#endif
         end if

      end if

      if (present(use_kmeans)) then
         if (use_kmeans) then
            call self % kmeans % init(self % ndims, self % max_nclusters, 10)
            call self % kmeans % fit(x, centroids=self % centroids)
            self % centroids = self % kmeans % centroids
         end if
      end if

      self % nclusters = self % max_nclusters

   end subroutine GMM_reset
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_final(self)
!
!     ---------
!     Interface
!     ---------
      type(GMM_t), intent(inout) :: self


      nullify(self % centroids)
      nullify(self % g % logtau)
      nullify(self % g % mu)
      nullify(self % g % cov)
      nullify(self % g % covinv)
      nullify(self % g % logdet)

      ! Required for gfortran?
      if (self % initialized) then
         if (associated(self % g % storage)) deallocate(self % g % storage)
         if (allocated(self % g % tcov))     deallocate(self % g % tcov)
         if (allocated(self % prob))         deallocate(self % prob)
      end if

      call kMeans_final(self % kmeans)

      self % initialized = .false.

   end subroutine GMM_final
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_fit(self, x, info, centroids, reset, adapt, from_kmeans)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t),       intent(inout) :: self
      real(RP),           intent(in)    :: x(:,:)
      integer,  optional, intent(out)   :: info
      real(RP), optional, intent(in)    :: centroids(self % ndims, self % max_nclusters)
      logical,  optional, intent(in)    :: reset
      logical,  optional, intent(in)    :: adapt
      logical,  optional, intent(in)    :: from_kmeans
!
!     ---------------
!     Local variables
!     ---------------
      logical               :: kmeans_
      integer               :: info_
      integer               :: npts
      integer               :: iter
      logical               :: breakFlag
      real(RP)              :: llprev
      real(RP)              :: mudiff
      real(RP), allocatable :: minimum(:), maximum(:)
      integer               :: ierr

!
!     Checks
!     ------
      if (.not. self % initialized .or. size(x, dim=1) /= self % ndims) then
         if (present(info)) info = -1
         return
      end if

      if (present(centroids)) then
         if (size(centroids, dim=1) /= self % ndims .or. &
               size(centroids, dim=2) /= self % nclusters) then
            if (present(info)) info = -1
            return
         end if
      end if
!
!     Reset "empty" clusters
!     ----------------------
      if (present(reset)) then
         if (reset) self % nclusters = 0
      end if

      if (present(from_kmeans)) then
         kmeans_ = from_kmeans
      else
         kmeans_ = .false.
      end if

      if (present(centroids)) then
         call GMM_reset(self, centroids=centroids, use_kmeans=kmeans_, x=x)
      else
         call GMM_reset(self, use_kmeans=kmeans_, x=x)
      end if
!
!     Minimum distance between cluster centroids
!     ------------------------------------------
      minimum = minval(x, dim=2)
      maximum = maxval(x, dim=2)

      call MPI_MinMax(minimum, maximum)

      mudiff  = self % mutol * maxval(maximum - minimum)
!
!     Initialization
!     --------------
      npts = size(x, dim=2)
      self % with_kmeans = .false.

      if (allocated(self % prob)) then
         if (size(self % prob, dim=1) /= npts) then
            deallocate(self % prob)
            allocate(self % prob(npts, self % nclusters))
         end if
      else
         allocate(self % prob(npts, self % nclusters))
      end if
!
!     EM algorithm
!     ------------
      do iter = 1, self % maxiters
      associate(ndims     => self % ndims,     &
                nclusters => self % nclusters, &
                ll        => self % logL,      &
                prob      => self % prob       )

         llprev = ll
         call GMM_Estep(ndims, npts, nclusters, self % g, x, prob(:,1:nclusters), ll)
         call GMM_Mstep(ndims, npts, nclusters, self % g, x, prob(:,1:nclusters), self % zerotol)

         if (present(adapt)) then
            if (adapt) call GMM_adapt_clusters(ndims, nclusters, self % g, mudiff, info_)
         end if

         if (MPI_Process % isRoot) then
            breakFlag = abs((ll - llprev) / ll) <= self % lltol .or. nclusters < 1
         end if
#if defined(_HAS_MPI_)
         call MPI_Bcast(breakFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
         if (breakFlag) exit

      end associate
      end do
!
!     Check convergence
!     -----------------
      if (iter > self % maxiters) then
         associate(ndims => self % ndims, nclusters => self % nclusters)
            call self % kmeans % init(ndims, nclusters, self % maxiters)
            call self % kmeans % fit(x, info=info_, centroids=self % centroids(:,1:nclusters))
            self % centroids(:,1:nclusters) = self % kmeans % centroids
            self % with_kmeans = .true.
         end associate

         if (info_ == -1) then
            info_ = -2
         else
            info_ = info_ + self % maxiters
         end if

      else
!
!        Sort clusters
!        -------------
         call sort_clusters(self % ndims, self % nclusters, self % centroids(:,1:self % nclusters))
         info_ = iter

      end if

      if (present(info)) info = info_

   end subroutine GMM_fit
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_predict(self, x, info)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t),      intent(inout) :: self
      real(RP),          intent(in)    :: x(:,:)  ! [ndims, npts]
      integer, optional, intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j
      integer  :: ndims
      integer  :: npts
      integer  :: nclusters
      real(RP) :: logtau
      real(RP) :: mu(size(x, dim=1))
      real(RP) :: covinv(size(x, dim=1), size(x, dim=1))
      real(RP) :: logdet
      real(RP) :: lognorm

!
!     Checks
!     ------
      ndims = size(x, dim=1)
      if (self % ndims /= ndims) then
         if (present(info)) info = -1
         return
      end if
!
!     Compute probabilities and output status
!     ---------------------------------------
      npts = size(x, dim=2)
      nclusters = self % nclusters

      ! The probabilities collapse when using k-means
      if (self % with_kmeans) then
!$omp parallel do private(i, j) firstprivate(npts) shared(self)
         do i = 1, npts
            j = self % kmeans % clusters(i)
            self % prob(i,:) = 0.0_RP
            self % prob(i,j) = 1.0_RP
         end do
!$omp end parallel do

      else

!$omp parallel default(private) firstprivate(ndims, npts, nclusters) shared(self, x)
!$omp do collapse(2)
         do j = 1, nclusters
            do i = 1, npts
               logtau = self % g % logtau(j)
               mu = self % g % mu(:,j)
               covinv = self % g % covinv(:,:,j)
               logdet = self % g % logdet(j)
               self % prob(i,j) = GMM_logpdf(ndims, x(:,i), logtau, mu, covinv, logdet)
            end do
         end do
!$omp end do

!$omp do
         do i = 1, npts
            lognorm = logsumexp(self % prob(i,:))
            self % prob(i,:) = exp(self % prob(i,:) - lognorm)
         end do
!$omp end do
!$omp end parallel

      end if

      if (present(info)) info = 1

   end subroutine GMM_predict
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function GMM_logpdf(ndims, x, logtau, mu, Sinv, logdet)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in) :: ndims
      real(RP), intent(in) :: x(ndims)
      real(RP), intent(in) :: logtau
      real(RP), intent(in) :: mu(ndims)
      real(RP), intent(in) :: Sinv(ndims, ndims)
      real(RP), intent(in) :: logdet
      real(RP)             :: GMM_logpdf
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: xmu(ndims)


      xmu = x - mu
      GMM_logpdf = logtau - 0.5_RP * (ndims * LOG2PI + logdet + &
                   dot_product(xmu, matmul(Sinv, xmu)))

   end function GMM_logpdf
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_Estep(ndims, npts, nclusters, g, x, R, logL)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(in)    :: npts
      integer,              intent(in)    :: nclusters
      type(GaussianList_t), intent(in)    :: g
      real(RP),             intent(in)    :: x(ndims, npts)
      real(RP),             intent(inout) :: R(npts, nclusters)
      real(RP),             intent(out)   :: logL
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j
      integer  :: ierr
      real(RP) :: logtau
      real(RP) :: mu(ndims)
      real(RP) :: covinv(ndims, ndims)
      real(RP) :: logdet
      real(RP) :: lognorm

!
!     E-step loop
!     -----------
      logL = 0.0_RP
!$omp parallel default(private) firstprivate(ndims, npts, nclusters) shared(R, g, x) reduction(+:logL)
!$omp do collapse(2)
      do j = 1, nclusters
         do i = 1, npts
            logtau = g % logtau(j)
            mu = g % mu(:,j)
            covinv = g % covinv(:,:,j)
            logdet = g % logdet(j)
            R(i,j) = GMM_logpdf(ndims, x(:,i), logtau, mu, covinv, logdet)
         end do
      end do
!$omp end do

!$omp do
      do i = 1, npts
         lognorm = logsumexp(R(i,:))
         R(i,:) = exp(R(i,:) - lognorm)
         logL = logL + lognorm
      end do
!$omp end do
!$omp end parallel
!
!     Global log-likelihood (only in root)
!     ------------------------------------
#if defined(_HAS_MPI_)
      if (MPI_Process % isRoot) then
         call MPI_Reduce(MPI_IN_PLACE, logL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call MPI_Reduce(logL, logL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      end if
#endif

   end subroutine GMM_Estep
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_Mstep(ndims, npts, nclusters, g, x, R, covtol)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(in)    :: npts
      integer,              intent(in)    :: nclusters
      type(GaussianList_t), intent(inout) :: g
      real(RP),             intent(in)    :: x(ndims, npts)
      real(RP),             intent(in)    :: R(:,:)  ! [npts, nclusters]
      real(RP),             intent(in)    :: covtol
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, l, m
      integer  :: tid
      integer  :: info
      real(RP) :: ns(nclusters)
      real(RP) :: inv_ns(nclusters)
      real(RP) :: xmu(ndims)
      real(RP) :: Rij
      real(RP) :: inv_sumtau
      real(RP) :: det

!
!     "Number of points" in each cluster
!     ----------------------------------
      ns = sum(R, dim=1)
      call MPI_SumAll(ns)
      ns = ns + 10.0_RP * epsilon(1.0_RP)  ! Avoid div by 0. Factor taken from scikit
      inv_ns = 1.0_RP / ns
!
!     M-step loop
!     -----------
      do j = 1, nclusters
         g % logtau(j) = ns(j)   ! No scaling yet
         g % mu(:,j) = matmul(x, R(:,j)) * inv_ns(j)
      end do

      ! Compute the global centroids
      call MPI_SumAll(g % mu(:,1:nclusters))

      ! Compute the local covariance matrices
!$omp parallel default(private) firstprivate(ndims, npts, nclusters) shared(R, g, x, inv_ns)
      tid = get_thread_id()
      do j = 1, nclusters

         g % tcov(:,:,tid) = 0.0_RP
!$omp do
         do i = 1, npts
            xmu = x(:,i) - g % mu(:,j)
            Rij = R(i,j)
            do m = 1, ndims
               do l = 1, ndims
                  g % tcov(l,m,tid) = g % tcov(l,m,tid) + Rij * xmu(l) * xmu(m)
               end do
            end do
         end do
!$omp end do

!$omp single
         ! Reduction over threads
         g % cov(:,:,j) = g % tcov(:,:,1)
         do i = 2, size(g % tcov, dim=3)
            g % cov(:,:,j) = g % cov(:,:,j) + g % tcov(:,:,i)
         end do

         ! Scaling
         g % cov(:,:,j) = g % cov(:,:,j) * inv_ns(j)
!$omp end single

      end do
!$omp end parallel

      ! Compute the global covariance matrices
      call MPI_SumAll(g % cov(:,:,1:nclusters))

      inv_sumtau = 1.0_RP / sum(g % logtau(1:nclusters))
      do j = 1, nclusters
         ! Regularize the covariance matrices
         do l = 1, ndims
            g % cov(l,l,j) = g % cov(l,l,j) + covtol
         end do

         ! Compute their inverse
         call matinv(ndims, g % cov(:,:,j), g % covinv(:,:,j), det)

         ! Update the "log" values
         g % logdet(j) = log(det)
         g % logtau(j) = log(g % logtau(j) * inv_sumtau)
      end do

   end subroutine GMM_Mstep
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_adapt_clusters(ndims, nclusters, g, mudiff, info)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(inout) :: nclusters
      type(GaussianList_t), intent(inout) :: g
      real(RP),             intent(in)    :: mudiff
      integer, optional,    intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: j, k, l
      integer  :: init_nclusters
      integer  :: bcast_req(2)
      logical  :: deleteit
      real(RP) :: tau_scale

!
!     Make sure the clusters do not overlap
!     -------------------------------------
      if (MPI_Process % isRoot) then
         k = 1
         init_nclusters = nclusters
         do j = 1, init_nclusters

            deleteit = .false.
            do l = 1, k-1
               if (all(g % mu(:,k) - g % mu(:,l) < mudiff)) then
                  deleteit = .true.
                  exit
               end if
            end do

            ! Delete the marked clusters
            if (deleteit) then
               do l = k+1, nclusters
                  g % logtau(l-1) = g % logtau(l)
                  g % mu(:,l-1) = g % mu(:,l)
                  g % cov(:,:,l-1) = g % cov(:,:,l)
                  g % covinv(:,:,l-1) = g % covinv(:,:,l)
                  g % logdet(l-1) = g % logdet(l)
               end do
               nclusters = nclusters - 1
            else
               k = k + 1
            end if

         end do

         ! Rescale tau to make sure it adds to 1
         tau_scale = logsumexp(g % logtau(1:nclusters))
         g % logtau(1:nclusters) = g % logtau(1:nclusters) - tau_scale
      end if
!
!     Syncronize
!     ----------
#if defined(_HAS_MPI_)
      call MPI_IBcast(nclusters, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_req(1), info)
      call MPI_IBcast(g % storage, size(g % storage), MPI_DOUBLE, 0, MPI_COMM_WORLD, bcast_req(2), info)
      call MPI_Waitall(2, bcast_req, MPI_STATUSES_IGNORE, info)
#endif

   end subroutine GMM_adapt_clusters
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine sort_clusters(ndims, nclusters, xavg, clustermap)
!
!     ---------
!     Interface
!     ---------
      integer,           intent(in)  :: ndims
      integer,           intent(in)  :: nclusters
      real(RP),          intent(out) :: xavg(ndims, nclusters)
      integer, optional, intent(out) :: clustermap(nclusters)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, ierr
      integer  :: indices(nclusters)
      integer  :: clustermap_(nclusters)
      real(RP) :: normvec(nclusters)
      real(RP) :: xavg_tmp(ndims, nclusters)

!
!     Sort the clusters by their distance to the origin
!     -------------------------------------------------
      if (MPI_Process % isRoot) then
         normvec = norm2(xavg, dim=1)
         indices = [(i, i = 1, nclusters)]
         call BubblesortWithFriend(normvec, indices)
         do i = 1, nclusters
            clustermap_(indices(i)) = i
         end do
      end if
#if defined(_HAS_MPI_)
      call MPI_Bcast(clustermap_, size(clustermap_), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

      xavg_tmp = xavg
      do i = 1, nclusters
         j = clustermap_(i)
         xavg(:,j) = xavg_tmp(:,i)
      end do

      if (present(clustermap)) then
         clustermap = clustermap_
      end if

   end subroutine sort_clusters
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function logsumexp(x) result(lse)
!
!     ---------
!     Interface
!     ---------
      real(RP), intent(in) :: x(:)
      real(RP)             :: lse
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      real(RP) :: m


      m = maxval(x)
      lse = 0.0_RP
      do i = 1, size(x, dim=1)
         lse = lse + exp(x(i) - m)
      end do
      lse = m + log(lse)

   end function logsumexp
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine matinv(n, A, Ainv, det)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)  :: n
      real(RP), intent(in)  :: A(n,n)
      real(RP), intent(out) :: Ainv(n,n)
      real(RP), intent(out) :: det
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: invdet


      if (n == 1) then
         ! Determinant
         det = A(1,1)

         ! Inverse
         Ainv(1,1) = 1.0_RP / A(1,1)

      elseif (n == 2) then
         ! Determinant
         det = A(1,1) * A(2,2) - A(1,2) * A(2,1)
         invdet = 1.0_RP / det

         ! Inverse
         Ainv(1,1) =  A(2,2) * invdet
         Ainv(2,1) = -A(2,1) * invdet
         Ainv(1,2) = -A(1,2) * invdet
         Ainv(2,2) =  A(1,1) * invdet

      elseif (n == 3) then
         ! Determinant
         det = A(1,1) * A(2,2) * A(3,3) + &
               A(1,2) * A(2,3) * A(3,1) + &
               A(1,3) * A(2,1) * A(3,2) - &
               A(1,3) * A(2,2) * A(3,1) - &
               A(1,2) * A(2,1) * A(3,3) - &
               A(1,1) * A(2,3) * A(3,2)
         invdet = 1.0_RP / det

         ! Inverse
         Ainv(1,1) = (A(2,2) * A(3,3) - A(2,3) * A(3,2)) * invdet
         Ainv(2,1) = (A(2,3) * A(3,1) - A(2,1) * A(3,3)) * invdet
         Ainv(3,1) = (A(2,1) * A(3,2) - A(2,2) * A(3,1)) * invdet
         Ainv(1,2) = (A(1,3) * A(3,2) - A(1,2) * A(3,3)) * invdet
         Ainv(2,2) = (A(1,1) * A(3,3) - A(1,3) * A(3,1)) * invdet
         Ainv(3,2) = (A(1,2) * A(3,1) - A(1,1) * A(3,2)) * invdet
         Ainv(1,3) = (A(1,2) * A(2,3) - A(1,3) * A(2,2)) * invdet
         Ainv(2,3) = (A(1,3) * A(2,1) - A(1,1) * A(2,3)) * invdet
         Ainv(3,3) = (A(1,1) * A(2,2) - A(1,2) * A(2,1)) * invdet

      else
         call LDL_factor(n, A, Ainv, det)
         call LDL_inverse(n, Ainv)

      end if

   end subroutine matinv
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine LDL_factor(n, A, LD, det)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)  :: n
      real(RP), intent(in)  :: A(n,n)
      real(RP), intent(out) :: LD(n,n)
      real(RP), intent(out) :: det
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, k
      real(RP) :: s


      do j = 1, n
!
!        Diagonal at `j`
!        ---------------
         LD(j,j) = A(j,j)
         do k = 1, j - 1
            LD(j,j) = LD(j,j) - LD(j,k)**2 * LD(k,k)
         end do
!
!        Column `j` with `i > j`
!        -----------------------
         do i = j + 1, n
            LD(i,j) = A(i,j)
            do k = 1, j - 1
               LD(i,j) = LD(i,j) - LD(i,k) * LD(j,k) * LD(k,k)
            end do
            LD(i,j) = LD(i,j) / LD(j,j)
         end do

      end do
!
!     Determinant
!     -----------
      det = LD(1,1)
      do i = 2, n
         det = det * LD(i,i)
      end do

   end subroutine LDL_factor

   pure subroutine LDL_inverse(n, A)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)    :: n
      real(RP), intent(inout) :: A(n,n)
!
!     ---------------
!     Local variables
!     ---------------
      integer :: i
      integer :: j
      integer :: k

!
!     Inverse following: https://arxiv.org/pdf/1111.4144.pdf
!     ------------------------------------------------------
      do j = n, 1, -1

         ! Diagonal term
         A(j,j) = 1.0_RP / A(j,j)
         do k = j + 1, n
            A(j,j) = A(j,j) - A(k,j) * A(j,k)
         end do

         ! Upper triangular terms
         do i = j - 1, 1, -1
            A(i,j) = 0.0_RP
            do k = i + 1, j
               A(i,j) = A(i,j) - A(k,i) * A(k,j)
            end do
            do k = j + 1, n
               A(i,j) = A(i,j) - A(k,i) * A(j,k)
            end do
         end do

      end do
!
!     Fill the lower triangle
!     -----------------------
      do i = 1, n - 1
         A(i + 1:,i) = A(i,i + 1:)
      end do

   end subroutine LDL_inverse
!
!///////////////////////////////////////////////////////////////////////////////
!
   ! TODO: move to OpenMP module
   function get_num_threads() result(n)
!
!     -------
!     Modules
!     -------
      use omp_lib
!
!     ---------
!     Interface
!     ---------
      integer :: n


#if defined(_OPENMP)
!$omp parallel
!$omp single
      n = omp_get_num_threads()
!$omp end single
!$omp end parallel
#else
      n = 1
#endif

   end function get_num_threads
!
!///////////////////////////////////////////////////////////////////////////////
!
   ! TODO: move to OpenMP module
   function get_thread_id() result(id)
!
!     -------
!     Modules
!     -------
      use omp_lib
!
!     ---------
!     Interface
!     ---------
      integer :: id


#if defined(_OPENMP)
      id = omp_get_thread_num() + 1
#else
      id = 1
#endif

   end function get_thread_id
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine standardize(x)
!
!     ---------
!     Interface
!     ---------
      real(RP), intent(inout) :: x(:,:)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: ndims
      integer  :: npts


      ndims = size(x, dim=1)
      npts  = size(x, dim=2)

      call standardize_(ndims, npts, x)

   end subroutine standardize

   subroutine standardize_(ndims, npts, x)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)    :: ndims
      integer,  intent(in)    :: npts
      real(RP), intent(inout) :: x(:,:)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: npts_total
      real(RP) :: inv_npts
      real(RP) :: mean(ndims)
      real(RP) :: var(ndims)

#if defined(_HAS_MPI_)
      integer  :: ierr
      integer  :: req(2)
#endif

!
!     Total number of points
!     ----------------------
#if defined(_HAS_MPI_)
      call MPI_IAllReduce(npts, npts_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, req(1), ierr)
#else
      npts_total = npts
#endif
!
!     Compute the average
!     -------------------
      mean = sum(x, dim=2)

#if defined(_HAS_MPI_)
      call MPI_IAllReduce(MPI_IN_PLACE, mean, ndims, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(size(req), req, MPI_STATUSES_IGNORE, ierr)
#endif

      inv_npts = 1.0_RP / npts_total
      mean = mean * inv_npts
!
!     Compute the variance
!     --------------------
      do i = 1, ndims
         var(i) = sum((x(i, :) - mean(i))**2) * inv_npts
      end do

#if defined(_HAS_MPI_)
      call MPI_AllReduce(MPI_IN_PLACE, var, ndims, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
!
!     Standardize data
!     ----------------
      do i = 1, ndims
         x(i,:) = (x(i,:) - mean(i)) / sqrt(var(i))
      end do

   end subroutine standardize_
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine rescale(x, s, o)
!
!     ---------
!     Interface
!     ---------
      real(RP),           intent(inout) :: x(:,:)
      real(RP), optional, intent(in)    :: s
      real(RP), optional, intent(in)    :: o
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: ndims
      real(RP) :: factor
      real(RP) :: offset


      ndims = size(x, dim=1)

      if (present(s)) then
         factor = s
      else
         factor = 1.0_RP
      end if

      if (present(o)) then
         offset = o
      else
         offset = 0.0_RP
      end if

      call rescale_(ndims, x, factor, offset)

   end subroutine rescale

   subroutine rescale_(ndims, x, scale, offset)
!
!     -------
!     Modules
!     -------
      use Utilities, only: AlmostEqual
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)    :: ndims
      real(RP), intent(inout) :: x(:,:)
      real(RP), intent(in)    :: scale
      real(RP), intent(in)    :: offset
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      real(RP) :: minimum(ndims)
      real(RP) :: maximum(ndims)


      minimum = minval(x, dim=2)
      maximum = maxval(x, dim=2)

      call MPI_MinMax(minimum, maximum)

      do i = 1, ndims
         if (AlmostEqual(maximum(i), minimum(i))) then
            if (maximum(i) > 0.0_RP) then
               x(i,:) = scale + offset
            else
               x(i,:) = offset
            end if
         else
            x(i,:) = (x(i,:) - minimum(i)) / (maximum(i) - minimum(i)) * scale + offset
         end if
      end do

   end subroutine rescale_

end module Clustering
