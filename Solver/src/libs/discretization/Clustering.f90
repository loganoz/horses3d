#include "Includes.h"
module Clustering

   use SMConstants,       only: RP, PI
   use Utilities,         only: AlmostEqual, BubblesortWithFriend
   use MPI_Process_Info,  only: MPI_Process
   use MPI_Utilities,     only: MPI_SumAll, MPI_MinMax
#ifdef _HAS_MPI_
   use mpi
#endif

   implicit none

   private

   public :: kMeans_t
   public :: GMM_t

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


      if (self % initialized) then
         call random_number(self % centroids)
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
         if (MPI_Process % doMPIAction) then
            call MPI_AllReduce(MPI_IN_PLACE, breakFlag, 1, MPI_LOGICAL, MPI_LAND, &
                               MPI_COMM_WORLD, ierr)
         end if
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
      integer  :: i1, i2
      integer  :: i, j

!
!     Keep the last state if all the clusters are used
!     ------------------------------------------------
      if (self % nclusters == self % max_nclusters) then
         return
      end if

      i1 = self % nclusters + 1
      i2 = self % max_nclusters
!
!     Set default values for the new clusters
!     ---------------------------------------
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
!     Set the centroids
!     -----------------
      if (present(centroids)) then
         self % centroids(:,i1:i2) = centroids(:,i1:i2)
      else
         call random_number(self % centroids(:,i1:i2))
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
      mudiff  = self % mutol * maxval(maximum - minimum)
!
!     Trivial case
!     ------------
      if (self % nclusters == 1) then
         call MPI_MinMax(minimum, maximum)
         self % centroids(:,1) = (minimum + maximum) * 0.5_RP
         if (present(info)) info = 0
         return
      end if
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
         if (MPI_Process % doMPIAction) then
            call MPI_Bcast(breakFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
         end if
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
      associate(prob => self % prob)
      if (self % with_kmeans) then
         prob = 0.0_RP
         do i = 1, npts
            j = self % kmeans % clusters(i)
            prob(i,j) = 1.0_RP
         end do

      else
         do j = 1, nclusters
            logtau = self % g % logtau(j)
            mu = self% g % mu(:,j)
            covinv = self % g % covinv(:,:,j)
            logdet = self % g % logdet(j)
            do i = 1, npts
               prob(i,j) = GMM_logpdf(ndims, x(:,i), logtau, mu, covinv, logdet)
            end do
         end do
         do i = 1, npts
            lognorm = logsumexp(prob(i,:))
            prob(i,:) = exp(prob(i,:) - lognorm)
         end do

      end if
      end associate

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
      do j = 1, nclusters
         logtau = g % logtau(j)
         mu = g % mu(:,j)
         covinv = g % covinv(:,:,j)
         logdet = g % logdet(j)
         do i = 1, npts
            R(i,j) = GMM_logpdf(ndims, x(:,i), logtau, mu, covinv, logdet)
         end do
      end do

      logL = 0.0_RP
      do i = 1, npts
         lognorm = logsumexp(R(i,:))
         R(i,:) = exp(R(i,:) - lognorm)
         logL = logL + lognorm
      end do
!
!     Global log-likelihood (only in root)
!     ------------------------------------
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIRootAction) then
         call MPI_Reduce(MPI_IN_PLACE, logL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      elseif (MPI_Process % doMPIAction) then
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
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIRootAction) then
         call MPI_Reduce(MPI_IN_PLACE, g % mu, size(g % mu), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      elseif (MPI_Process % doMPIAction) then
         call MPI_Reduce(g % mu, g % mu, size(g % mu), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      end if
#endif

      ! Compute the local covariance matrices
      do j = 1, nclusters
         g % cov(:,:,j) = 0.0_RP
         do i = 1, npts
            xmu = x(:,i) - g % mu(:,j)
            Rij = R(i,j)
            do m = 1, ndims
               do l = 1, ndims
                  g % cov(l,m,j) = g % cov(l,m,j) + Rij * xmu(l) * xmu(m)
               end do
            end do
         end do
         g % cov(:,:,j) = g % cov(:,:,j) * inv_ns(j)
      end do

      ! Compute the global covariance matrices
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIRootAction) then
         call MPI_Reduce(MPI_IN_PLACE, g % cov, size(g % cov), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      elseif (MPI_Process % doMPIAction) then
         call MPI_Reduce(g % cov, g % cov, size(g % cov), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      end if
#endif

      inv_sumtau = 1.0_RP / sum(g % logtau)
      do j = 1, nclusters
         ! Regularize the covariance matrices
         do l = 1, ndims
            g % cov(l,l,j) = g % cov(l,l,j) + covtol
         end do

         ! Compute their inverse
         det = matdet(ndims, g % cov(:,:,j))
         call matinv(ndims, det, g % cov(:,:,j), g % covinv(:,:,j))

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
         tau_scale = 1.0_RP / sum(exp(g % logtau(1:nclusters)))
         g % logtau(1:nclusters) = g % logtau(1:nclusters) + log(tau_scale)
      end if
!
!     Syncronize
!     ----------
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         call MPI_IBcast(nclusters, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_req(1), info)
         call MPI_IBcast(g % storage, size(g % storage), MPI_DOUBLE, 0, MPI_COMM_WORLD, bcast_req(2), info)
         call MPI_Waitall(2, bcast_req, MPI_STATUSES_IGNORE, info)
      end if
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
      if (MPI_Process % doMPIAction) then
         call MPI_Bcast(clustermap_, size(clustermap_), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
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
   pure function matdet(n, A)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in) :: n
      real(RP), intent(in) :: A(:,:)
      real(RP)             :: matdet

      if (n == 1) then
          matdet = A(1,1)

      elseif (n == 2) then
          matdet = A(1,1) * A(2,2) - A(1,2) * A(2,1)

      elseif (n == 3) then
          matdet = A(1,1) * A(2,2) * A(3,3) + &
                   A(1,2) * A(2,3) * A(3,1) + &
                   A(1,3) * A(2,1) * A(3,2) - &
                   A(1,3) * A(2,2) * A(3,1) - &
                   A(1,2) * A(2,1) * A(3,3) - &
                   A(1,1) * A(2,3) * A(3,2)

      else
          matdet = 0.0_RP

      end if

   end function matdet
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine matinv(n, det, A, Ainv)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)  :: n
      real(RP), intent(in)  :: det
      real(RP), intent(in)  :: A(n,n)
      real(RP), intent(out) :: Ainv(n,n)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: invd


      if (n == 1) then
         Ainv(1,1) = 1.0_RP / A(1,1)

      elseif (n == 2) then
         invd = 1.0_RP / det
         Ainv(1,1) =  A(2,2) * invd
         Ainv(2,1) = -A(2,1) * invd
         Ainv(1,2) = -A(1,2) * invd
         Ainv(2,2) =  A(1,1) * invd

      elseif (n == 3) then
         invd = 1.0_RP / det
         Ainv(1,1) = (A(2,2) * A(3,3) - A(2,3) * A(3,2)) * invd
         Ainv(2,1) = (A(2,3) * A(3,1) - A(2,1) * A(3,3)) * invd
         Ainv(3,1) = (A(2,1) * A(3,2) - A(2,2) * A(3,1)) * invd
         Ainv(1,2) = (A(1,3) * A(3,2) - A(1,2) * A(3,3)) * invd
         Ainv(2,2) = (A(1,1) * A(3,3) - A(1,3) * A(3,1)) * invd
         Ainv(3,2) = (A(1,2) * A(3,1) - A(1,1) * A(3,2)) * invd
         Ainv(1,3) = (A(1,2) * A(2,3) - A(1,3) * A(2,2)) * invd
         Ainv(2,3) = (A(1,3) * A(2,1) - A(1,1) * A(2,3)) * invd
         Ainv(3,3) = (A(1,1) * A(2,2) - A(1,2) * A(2,1)) * invd

      else
         ! TODO: Add generic method for symmetric matrices

      end if

   end subroutine matinv

end module Clustering
