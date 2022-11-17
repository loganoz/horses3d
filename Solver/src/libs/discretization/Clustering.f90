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
      integer,  allocatable :: prevClusters(:)
      real(RP), allocatable :: centroids(:,:)
   contains
      procedure :: init  => kMeans_init
      procedure :: reset => kMeans_reset
      procedure :: fit   => kMeans_fit
      final     :: kMeans_final
   end type kMeans_t

   type :: GaussianList_t
      real(RP), pointer, contiguous :: storage(:)    => null()  ! contiguous storage for MPI communication
      real(RP), pointer, contiguous :: taumu_st(:)   => null()  ! contiguous storage for MPI communication
      real(RP), pointer, contiguous :: logtau(:)     => null()  ! [n]: log of weights
      real(RP), pointer, contiguous :: mu(:,:)       => null()  ! [d,n]: centroids
      real(RP), pointer, contiguous :: cov(:,:,:)    => null()  ! [d,d,n]: covariance matrices
      real(RP), pointer, contiguous :: covinv(:,:,:) => null()  ! [d,d,n]: inverse covariance matrices
      real(RP), pointer, contiguous :: logdet(:)     => null()  ! [n]: log-determinant of covariance matrices
   end type GaussianList_t

   type :: GMM_t
      logical                       :: initialized = .false.
      integer                       :: ndims
      integer                       :: nclusters
      integer                       :: max_nclusters
      real(RP), allocatable         :: R(:,:)
      real(RP), pointer, contiguous :: centroids(:,:)
      type(GaussianList_t)          :: g
      type(kMeans_t)                :: kmeans
   contains
      procedure :: init  => GMM_init
      procedure :: reset => GMM_reset
      procedure :: fit   => GMM_fit
      final     :: GMM_final
   end type GMM_t

   real(RP), parameter :: LOG2PI = log(2.0_RP*PI)
!
!  ========
   contains
!  ========
!
   subroutine kMeans_init(self, ndims, nclusters)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t), intent(inout) :: self
      integer,         intent(in)    :: ndims
      integer,         intent(in)    :: nclusters


      self % ndims     = ndims
      self % nclusters = nclusters

      if (allocated(self % centroids)) deallocate(self % centroids)
      allocate(self % centroids(ndims, nclusters))

      call self % reset()
      self % initialized = .true.

   end subroutine kMeans_init
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_reset(self)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t), intent(inout) :: self


      call random_number(self % centroids)

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


      if (allocated(self % prevClusters)) deallocate(self % prevClusters)
      if (allocated(self % centroids))    deallocate(self % centroids)

      self % initialized = .false.

   end subroutine kMeans_final
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_fit(self, x, clusters, info, centroids)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t),    intent(inout) :: self
      real(RP),           intent(in)    :: x(:,:)
      integer,            intent(out)   :: clusters(:)
      integer,  optional, intent(out)   :: info
      real(RP), optional, intent(in)    :: centroids(self % ndims, self % nclusters)
!
!     ---------------
!     Local variables
!     ---------------
      integer, parameter :: maxIters = 50
      integer            :: npts
      integer            :: i, ierr
      logical            :: breakFlag


      if (.not. self % initialized) then
         if (present(info)) info = -1
         return
      end if
!
!     Initial clusters
!     ----------------
      npts = size(x, dim=2)
      if (present(centroids)) then
         self % centroids = centroids
      end if
      call kMeans_compute_clusters(self % nclusters, self % ndims, npts, x, &
                                   self % centroids, clusters)
!
!     Loop until convergence
!     ----------------------
      do i = 1, maxIters
      associate(ndims => self % ndims, nclusters => self % nclusters)

         self % prevClusters = clusters
         call kMeans_compute_centroids(nclusters, ndims, npts, x, self % centroids, clusters)
         call kMeans_compute_clusters(nclusters, ndims, npts, x, self % centroids, clusters)

         breakFlag = all(self % prevClusters == clusters)
#if defined(_HAS_MPI_)
         if (MPI_Process % doMPIAction) then
            call MPI_AllReduce(MPI_IN_PLACE, breakFlag, 1, MPI_LOGICAL, MPI_LAND, &
                               MPI_COMM_WORLD, ierr)
         end if
#endif
         if (breakFlag) exit

      end associate
      end do
!
!     Check convergence
!     -----------------
      if (present(info)) info = merge(-1, i, i > maxIters)
!
!     Sort clusters
!     -------------
      call kMeans_sort_clusters(self % nclusters, self % ndims, npts, self % centroids, clusters)

   end subroutine kMeans_fit
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
      real(RP) :: xavg_tmp(ndims,nclusters)
      integer  :: indices(nclusters)
      integer  :: clustermap(nclusters)

!
!     Sort the clusters by their distance to the origin
!     -------------------------------------------------
      if (MPI_Process % isRoot) then
         normvec = norm2(xavg, dim=1)
         indices = [(i, i = 1, nclusters)]
         call BubblesortWithFriend(normvec, indices)
         do i = 1, nclusters
            clustermap(indices(i)) = i
         end do
      end if
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         call MPI_Bcast(clustermap, size(clustermap), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
!
!     Assign nodes to clusters
!     ------------------------
      do i = 1, npts
         clusters(i) = clustermap(clusters(i))
      end do

      xavg_tmp = xavg
      do i = 1, nclusters
         j = clustermap(i)
         xavg(:,j) = xavg_tmp(:,i)
      end do

   end subroutine kMeans_sort_clusters
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_init(self, ndims, nclusters)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t), intent(inout) :: self
      integer,      intent(in)    :: ndims
      integer,      intent(in)    :: nclusters
!
!     ---------------
!     Local variables
!     ---------------
      integer :: tsize, msize, csize, dsize
      integer :: ind1, ind2
      integer :: i, j

!
!     Reset all pointers
!     ------------------
      nullify(self % centroids)
      nullify(self % g % taumu_st)
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

      ! Required for gfortran?
      if (self % initialized) then
         if (associated(self % g % storage)) deallocate(self % g % storage)
      end if
      allocate(self % g % storage(tsize + msize + 2 * csize + dsize))

      self % g % taumu_st(1:tsize + msize) => self % g % storage(1:tsize + msize)

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
      self % initialized   = .true.
      self % ndims         = ndims
      self % max_nclusters = nclusters
      call self % reset()

   end subroutine GMM_init
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_reset(self)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t), intent(inout) :: self
!
!     ---------------
!     Local variables
!     ---------------
      integer :: i, j


      if (self % initialized) then
         self % nclusters = 0
         call GMM_reset_empty(self)
      end if

   end subroutine GMM_reset
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_reset_empty(self)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t), intent(inout) :: self
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i1, i2
      integer  :: i, j


      if (self % nclusters == self % max_nclusters) then
         return
      end if

      i1 = self % nclusters + 1
      i2 = self % max_nclusters
!
!     Set default values for the new clusters
!     ---------------------------------------
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

      call random_number(self % g % mu(:, i1:i2))

      self % nclusters = self % max_nclusters

   end subroutine GMM_reset_empty
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
      nullify(self % g % taumu_st)
      nullify(self % g % logtau)
      nullify(self % g % mu)
      nullify(self % g % cov)
      nullify(self % g % covinv)
      nullify(self % g % logdet)

      ! Required for gfortran?
      if (self % initialized) then
         if (associated(self % g % storage)) deallocate(self % g % storage)
      end if

      if (allocated(self % R)) deallocate(self % R)

      self % initialized = .false.

   end subroutine GMM_final
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_fit(self, x, clusters, info, centroids)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t),       intent(inout) :: self
      real(RP),           intent(in)    :: x(:,:)
      integer,            intent(out)   :: clusters(:)
      integer,  optional, intent(out)   :: info
      real(RP), optional, intent(in)    :: centroids(self % ndims, self % max_nclusters)
!
!     ---------------
!     Local variables
!     ---------------
      integer,  parameter :: maxIters = 100                    ! Magic (scikit)!!
      real(RP), parameter :: mutol    = 2e-5_RP                ! Magic
      real(RP), parameter :: lltol    = 1e-3_RP                ! Magic (scikit)!!
      real(RP), parameter :: zerotol  = 1e-10_RP               ! Magic (scikit)!!

      integer               :: info_
      integer               :: npts
      integer               :: iter
      logical               :: breakFlag
      real(RP), save        :: ll = huge(1.0_RP)
      real(RP)              :: llprev
      real(RP)              :: mudiff
      real(RP), allocatable :: minimum(:), maximum(:)
      integer               :: i
      integer               :: ierr

!
!     Reset "empty" clusters
!     ----------------------
      if (present(centroids)) then
         call self % reset()
         self % centroids = centroids
      else
         if (self % nclusters < self % max_nclusters) ll = huge(1.0_RP)
         call GMM_reset_empty(self)
      end if
!
!     Minimum distance between cluster centroids
!     ------------------------------------------
      minimum = minval(x, dim=2)
      maximum = maxval(x, dim=2)
      mudiff  = mutol * maxval(maximum - minimum)
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

      if (allocated(self % R)) then
         if (size(self % R, dim=1) /= npts) then
            deallocate(self % R)
            allocate(self % R(npts, self % nclusters))
         end if
      else
         allocate(self % R(npts, self % nclusters))
      end if
!
!     EM algorithm
!     ------------
      do iter = 1, maxIters
      associate(ndims => self % ndims, nclusters => self % nclusters)

         llprev = ll
         call GMM_Estep(ndims, npts, nclusters, self % g, x, self % R(:,1:nclusters), zerotol, ll)
         call GMM_Mstep(ndims, npts, nclusters, self % g, x, self % R(:,1:nclusters), zerotol, mudiff)

         if (MPI_Process % isRoot) then
            breakFlag = abs((ll - llprev) / ll) <= lltol .or. self % nclusters < 1
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
      if (iter > maxIters) then

         call self % kmeans % init(self % ndims, self % nclusters)
         call self % kmeans % fit(x, clusters, info_, self % centroids(:,1:self % nclusters))

         if (info_ == -1) then
            info_ = -2
         end if

      else
!
!        Sort clusters
!        -------------
         call GMM_sort_clusters(self % ndims, npts, self % nclusters, self % g, x, &
                                self % centroids(:,1:self % nclusters), clusters)
         info_ = iter

      end if

      if (present(info)) info = info_

   end subroutine GMM_fit
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
   subroutine GMM_Estep(ndims, npts, nclusters, g, x, R, small, logL)
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
      real(RP),             intent(in)    :: small
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
      real(RP) :: den

!
!     E-step loop
!     -----------
      logL = 0.0_RP
      do j = 1, nclusters
         logtau = g % logtau(j)
         mu = g % mu(:,j)
         covinv = g % covinv(:,:,j)
         logdet = g % logdet(j)
         do i = 1, npts
            R(i,j) = exp(GMM_logpdf(ndims, x(:,i), logtau, mu, covinv, logdet))
         end do
      end do

      do i = 1, npts
         den = sum(R(i,:))
         if (AlmostEqual(den, 0.0_RP)) then  ! The point is far from all clusters
             R(i,:) = 1.0_RP / nclusters
             den = small
         else
             R(i,:) = R(i,:) / den
         end if
         logL = logL + log(den)
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
   subroutine GMM_Mstep(ndims, npts, nclusters, g, x, R, covtol, mudiff)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(in)    :: npts
      integer,              intent(inout) :: nclusters
      type(GaussianList_t), intent(inout) :: g
      real(RP),             intent(in)    :: x(ndims, npts)
      real(RP),             intent(in)    :: R(:,:)
      real(RP),             intent(in)    :: covtol
      real(RP),             intent(in)    :: mudiff
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, k, l, m
      integer  :: init_nclusters
      integer  :: info
      integer  :: bcast_req(2)
      logical  :: deleteit
      real(RP) :: ns(nclusters)
      real(RP) :: inv_nsAll(nclusters)
      real(RP) :: inv_nptsAll
      real(RP) :: xmu(ndims)
      real(RP) :: Rij
      real(RP) :: tau_scale

!
!     Total number of points
!     ----------------------
      inv_nptsAll = real(npts, kind=RP)
      call MPI_SumAll(inv_nptsAll)
      inv_nptsAll = 1.0_RP / inv_nptsAll
!
!     "Number of points" in the cluster
!     ---------------------------------
      ns = sum(R, dim=1)
      inv_nsAll = ns
      call MPI_SumAll(inv_nsAll)
      inv_nsAll = 1.0_RP / inv_nsAll
!
!     M-step loop
!     -----------
      do j = 1, nclusters
         g % logtau(j) = ns(j) * inv_nptsAll
         g % mu(:,j) = matmul(x, R(:,j)) * inv_nsAll(j)
      end do

      ! The covariance matrix needs the new centroids
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIRootAction) then
         call MPI_Reduce(MPI_IN_PLACE, g % taumu_st, size(g % taumu_st), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      elseif (MPI_Process % doMPIAction) then
         call MPI_Reduce(g % taumu_st, g % taumu_st, size(g % taumu_st), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      end if
#endif

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
         g % cov(:,:,j) = g % cov(:,:,j) * inv_nsAll(j)
      end do

#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIRootAction) then
         call MPI_Reduce(MPI_IN_PLACE, g % cov, size(g % cov), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      elseif (MPI_Process % doMPIAction) then
         call MPI_Reduce(g % cov, g % cov, size(g % cov), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      end if
#endif
!
!     Make sure the clusters are OK
!     -----------------------------
      if (MPI_Process % isRoot) then
         k = 1
         init_nclusters = nclusters
         do j = 1, init_nclusters

            ! Compute the inverse of the covariance matrix
            call matinv(ndims, g % cov(:,:,k), g % covinv(:,:,k), tol=covtol, info=info)

            ! If a cluster collapses, limit its size
            if (info < 0) then
               do l = 1, ndims
                  g % cov(:,l,k) = 0.0_RP
                  g % cov(l,l,k) = covtol**(1.0_RP/ndims)
                  g % covinv(:,l,k) = 0.0_RP
                  g % covinv(l,l,k) = 1.0_RP / g % cov(l,l,k)
               end do
            end if

            ! Update the "log" values
            g % logtau(k) = log(g % logtau(k))
            g % logdet(k) = log(matdet(ndims, g % cov(:,:,k)))

            ! Also look for overlapping clusters
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

   end subroutine GMM_Mstep
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_sort_clusters(ndims, npts, nclusters, g, x, xavg, clusters)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(in)    :: npts
      integer,              intent(in)    :: nclusters
      type(GaussianList_t), intent(in)    :: g
      real(RP),             intent(in)    :: x(ndims, npts)
      real(RP),             intent(out)   :: xavg(ndims, nclusters)
      integer,              intent(inout) :: clusters(npts)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, ierr
      real(RP) :: pdf, maxpdf
      real(RP) :: xi(ndims)
      real(RP) :: normvec(nclusters)
      integer  :: indices(nclusters)
      integer  :: clustermap(nclusters)

!
!     Sort the clusters by their distance to the origin
!     -------------------------------------------------
      if (MPI_Process % isRoot) then
         normvec = norm2(g % mu(:,1:nclusters), dim=1)
         indices = [(i, i = 1, nclusters)]
         call BubblesortWithFriend(normvec, indices)
         do i = 1, nclusters
            clustermap(indices(i)) = i
         end do
      end if
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         call MPI_Bcast(clustermap, size(clustermap), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
!
!     Assign nodes to clusters
!     ------------------------
      do i = 1, npts
         clusters(i) = 0
         maxpdf = 0.0_RP
         pdf = 0.0_RP
         xi = x(:,i)
         do j = 1, nclusters
            pdf = exp(GMM_logpdf(ndims, xi, g % logtau(j), g % mu(:,j), &
               g % covinv(:,:,j), g % logdet(j)))
            if (pdf > maxpdf) then
               maxpdf = pdf
               clusters(i) = clustermap(j)
            end if
         end do
      end do

      do i = 1, nclusters
         j = clustermap(i)
         xavg(:,j) = g % mu(:,i)
      end do

   end subroutine GMM_sort_clusters
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
   pure subroutine matinv(n, A, Ainv, tol, info)
!
!     ---------
!     Interface
!     ---------
      integer,            intent(in)  :: n
      real(RP),           intent(in)  :: A(n,n)
      real(RP),           intent(out) :: Ainv(n,n)
      real(RP), optional, intent(in)  :: tol
      integer,  optional, intent(out) :: info
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: d, invd
      logical  :: iszero
      integer  :: info_


      d = matdet(n, A)
      if (present(tol)) then
         iszero = AlmostEqual(d, 0.0_RP, tol)
      else
         iszero = AlmostEqual(d, 0.0_RP)
      end if

      if (.not. iszero) then
         if (n == 1) then
            info_ = 1
            Ainv(1,1) = 1.0_RP / A(1,1)

         elseif (n == 2) then
            info_ = 1
            invd = 1.0_RP / d
            Ainv(1,1) =  A(2,2) * invd
            Ainv(2,1) = -A(2,1) * invd
            Ainv(1,2) = -A(1,2) * invd
            Ainv(2,2) =  A(1,1) * invd

         elseif (n == 3) then
            info_ = 1
            invd = 1.0_RP / d
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
            info_ = -2

         end if

      else
          info_ = -1

      end if

      if (present(info)) info = info_

   end subroutine matinv

end module Clustering
