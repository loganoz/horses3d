module Clustering

   use SMConstants,       only: RP, PI
   use Utilities,         only: AlmostEqual, BubblesortWithFriend
   use MPI_Process_Info,  only: MPI_Process
   use MPI_Utilities,     only: MPI_OpAll, MPI_SumAll, MPI_MinMax
#ifdef _HAS_MPI_
   use mpi
#endif

   implicit none

   private

   public :: KMeans_t
   public :: GMM_t

   type :: KMeans_t
      integer               :: ndims
      integer               :: nclusters
      integer,  allocatable :: prevClusters(:)
   contains
      procedure :: init => KMeans_init
      procedure :: fit  => kMeans_fit
   end type KMeans_t

   type :: GaussianList_t
      real(RP), pointer, contiguous :: storage(:)      ! contiguous storage for MPI communication
      real(RP), pointer, contiguous :: tau(:)          ! [n]: weights
      real(RP), pointer, contiguous :: mu(:,:)         ! [d,n]: centroids
      real(RP), pointer, contiguous :: cov(:,:,:)      ! [d,d,n]: covariance matrices
      real(RP), pointer, contiguous :: covinv(:,:,:)   ! [d,d,n]: inverse covariance matrices
   end type GaussianList_t

   type :: GMM_t
      integer                       :: ndims
      integer                       :: nclusters
      type(GaussianList_t)          :: g
      real(RP),         allocatable :: R(:,:)
      type(KMeans_t)                :: kmeans
   contains
      procedure :: init => GMM_init
      procedure :: fit  => GMM_fit
      final     :: GMM_final
   end type GMM_t
!
!  ========
   contains
!  ========
!
   subroutine kMeans_init(self)
!
!     ---------
!     Interface
!     ---------
      class(kMeans_t), intent(inout) :: self


      ! Allocate in fit
      if (allocated(self % prevClusters)) deallocate(self % prevClusters)

   end subroutine kMeans_init
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine kMeans_fit(self, x, clusters, centroids, info)
!
!     ---------
!     Interface
!     ---------
      class(KMeans_t),   intent(inout) :: self
      real(RP),          intent(in)    :: x(:,:)
      integer,           intent(out)   :: clusters(:)
      real(RP),          intent(out)   :: centroids(:,:)
      integer, optional, intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      integer, parameter   :: maxIters = 50
      integer              :: npts
      integer              :: i, ierr
      logical              :: breakFlag

!
!     Initial clusters
!     ----------------
      npts = size(x, dim=2)
      call kMeans_compute_clusters(self % nclusters, self % ndims, npts, x, centroids, clusters)
!
!     Loop until convergence
!     ----------------------
      do i = 1, maxIters

         self % prevClusters = clusters

         call kMeans_compute_centroids(self % nclusters, self % ndims, npts, x, centroids, clusters)
         call kMeans_compute_clusters(self % nclusters, self % ndims, npts, x, centroids, clusters)

         breakFlag = all(self % prevClusters == clusters)
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
      if (present(info)) info = merge(-1, i, i > maxIters)
!
!     Sort clusters
!     -------------
      call kMeans_sort_clusters(self % nclusters, self % ndims, npts, centroids, clusters)

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
      real(RP), intent(in)  :: x(ndims,npts)
      real(RP), intent(out) :: xavg(ndims,nclusters)
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
   subroutine GMM_init(self)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t), intent(inout) :: self
!
!     ---------------
!     Local variables
!     ---------------
      integer :: ndims, nclusters
      integer :: tsize, msize, csize
      integer :: ind1, ind2
      integer :: i, j


      ndims = self % ndims
      nclusters = self % nclusters

      nullify(self % g % tau)
      nullify(self % g % mu)
      nullify(self % g % cov)
      nullify(self % g % covinv)

      if (associated(self % g % storage)) deallocate(self % g % storage)
      if (allocated(self % R))            deallocate(self % R)    ! Allocate in fit

      tsize = nclusters
      msize = ndims * nclusters
      csize = ndims * ndims * nclusters

      allocate(self % g % storage(tsize + msize + 2 * csize))

      ind1 = 1
      ind2 = tsize
      self % g % tau => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + msize
      self % g % mu(1:ndims, 1:nclusters) => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + csize
      self % g % cov(1:ndims, 1:ndims, 1:nclusters) => self % g % storage(ind1:ind2)

      ind1 = ind2 + 1
      ind2 = ind2 + csize
      self % g % covinv(1:ndims, 1:ndims, 1:nclusters) => self % g % storage(ind1:ind2)

      self % g % tau = 1.0_RP / nclusters
      self % g % cov = 0.0_RP
      self % g % covinv = 0.0_RP
      do i = 1, nclusters
         do j = 1, ndims
            self % g % cov(j,j,i) = 1.0_RP
            self % g % covinv(j,j,i) = 1.0_RP
         end do
      end do

      self % kmeans % ndims = ndims
      self % kmeans % nclusters = nclusters
      call self % kmeans % init()

   end subroutine GMM_init
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_final(self)
!
!     ---------
!     Interface
!     ---------
      type(GMM_t), intent(inout) :: self
!
!     ---------------
!     Local variables
!     ---------------
      integer :: ierr


      nullify(self % g % tau)
      nullify(self % g % mu)
      nullify(self % g % cov)
      nullify(self % g % covinv)

      if (associated(self % g % storage)) deallocate(self % g % storage)
      if (allocated(self % R))            deallocate(self % R)

   end subroutine GMM_final
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_fit(self, x, clusters, centroids, info)
!
!     ---------
!     Interface
!     ---------
      class(GMM_t),      intent(inout) :: self
      real(RP),          intent(in)    :: x(:,:)
      integer,           intent(out)   :: clusters(:)
      real(RP),          intent(out)   :: centroids(:,:)
      integer, optional, intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      real(RP), parameter :: tol = 1e-2_RP                     ! Magic!!
      real(RP), parameter :: small = 1e4_RP * epsilon(1.0_RP)  ! Magic!!
      integer,  parameter :: maxIters = 50

      integer               :: info_
      integer               :: npts
      integer               :: iter
      logical               :: breakFlag
      real(RP)              :: ll
      real(RP)              :: llprev
      real(RP), allocatable :: minimum(:), maximum(:)
      integer               :: i
      integer               :: ierr

!
!     Trivial case
!     ------------
      if (self % nclusters == 1) then
         minimum = minval(x, dim=2)
         maximum = maxval(x, dim=2)
         call MPI_MinMax(minimum, maximum)
         centroids(:,1) = (minimum + maximum) * 0.5_RP
         if (present(info)) info = 0
         return
      end if
!
!     Initialization
!     --------------
      npts = size(x, dim=2)
      self % g % mu = centroids
      if (size(self % R, dim=1) /= npts) then
         if (allocated(self % R)) deallocate(self % R)
         allocate(self % R(npts, self % nclusters))
      end if
!
!     EM algorithm
!     ------------
      ll = huge(1.0_RP)
      do iter = 1, maxIters

         llprev = ll
         call GMM_Estep(self % ndims, npts, self % nclusters, self % g, x, self % R, small, ll)
         call GMM_Mstep(self % ndims, npts, self % nclusters, self % g, x, self % R, small)

         if (MPI_Process % isRoot) then
            breakFlag = abs((ll - llprev) / ll) <= tol .or. self % nclusters < 1
         end if
#if defined(_HAS_MPI_)
         if (MPI_Process % doMPIAction) then
            call MPI_Bcast(breakFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
         end if
#endif
         if (breakFlag) exit

      end do
!
!     Check convergence
!     -----------------
      if (iter > maxIters) then
         call self % kmeans % fit(x, clusters, centroids, info_)
         if (info_ == -1) then
            info_ = -2
         end if
      else
         info_ = iter
      end if
      if (present(info)) info = info_
!
!     Sort clusters
!     -------------
      call GMM_sort_clusters(self % ndims, npts, self % nclusters, self % g, x, centroids, clusters)

   end subroutine GMM_fit
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function GMM_pdf(ndims, x, mu, S, Sinv)
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in) :: ndims
      real(RP), intent(in) :: x(ndims)
      real(RP), intent(in) :: mu(ndims)
      real(RP), intent(in) :: S(ndims, ndims)
      real(RP), intent(in) :: Sinv(ndims, ndims)
      real(RP)             :: GMM_pdf
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: pi2
      real(RP) :: xmu(ndims)
      real(RP) :: det


      pi2 = 2.0_RP * PI
      xmu = x - mu
      det = matdet(ndims, S)
      GMM_pdf = exp(-0.5_RP * dot_product(xmu, matmul(Sinv, xmu))) / &
                    sqrt(pi2 ** ndims * det)

   end function GMM_pdf
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
      real(RP),             intent(in)    :: x(:,:)
      real(RP),             intent(inout) :: R(:,:)
      real(RP),             intent(in)    :: small
      real(RP),             intent(out)   :: logL
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j
      integer  :: ierr
      real(RP) :: tau
      real(RP) :: mu(ndims)
      real(RP) :: cov(ndims,ndims)
      real(RP) :: covinv(ndims,ndims)
      real(RP) :: den

!
!     E-step loop
!     -----------
      logL = 0.0_RP
      do j = 1, nclusters
         tau = g % tau(j)
         mu = g % mu(:,j)
         cov = g % cov(:,:,j)
         covinv = g % covinv(:,:,j)
         do i = 1, npts
            R(i,j) = tau * GMM_pdf(ndims, x(:,i), mu, cov, covinv)
         end do
      end do

      do i = 1, npts
         den = sum(R(i,1:nclusters))
         if (AlmostEqual(den, 0.0_RP)) then  ! The point is far from all clusters
             R(i,1:nclusters) = 1.0_RP / nclusters
             den = small
         else
             R(i,1:nclusters) = R(i,1:nclusters) / den
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
   subroutine GMM_Mstep(ndims, npts, nclusters, g, x, R, tol)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(in)    :: npts
      integer,              intent(inout) :: nclusters
      type(GaussianList_t), intent(inout) :: g
      real(RP),             intent(in)    :: x(:,:)
      real(RP),             intent(in)    :: R(:,:)
      real(RP),             intent(in)    :: tol
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, k, l, m
      integer  :: nptsAll
      integer  :: init_nclusters
      integer  :: info
      logical  :: deleteit
      real(RP) :: ns(nclusters)
      real(RP) :: xmu(ndims)

!
!     Total number of points
!     ----------------------
      nptsAll = npts
      call MPI_SumAll(nptsAll)
!
!     "Number of points" in the cluster
!     ---------------------------------
      ns = sum(R, dim=1)
      call MPI_SumAll(ns)
!
!     M-step loop
!     -----------
      do j = 1, nclusters
         g % tau(j) = ns(j) / nptsAll
         g % mu(:,j) = matmul(x, R(:,j)) / ns(j)
         g % cov(:,:,j) = 0.0_RP
         do i = 1, npts
            xmu = x(:,i) - g % mu(:,j)
            do m = 1, ndims
               do l = 1, ndims
                  g % cov(l,m,j) = g % cov(l,m,j) + R(i,j) * xmu(l) * xmu(m)
               end do
            end do
         end do
         g % cov(:,:,j) = g % cov(:,:,j) / ns(j)
      end do
!
!     Make sure the clusters are OK
!     -----------------------------
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIRootAction) then
         call MPI_Reduce(MPI_IN_PLACE, g % storage, size(g % storage), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      elseif (MPI_Process % doMPIAction) then
         call MPI_Reduce(g % storage, g % storage, size(g % storage), MPI_DOUBLE, &
                         MPI_SUM, 0, MPI_COMM_WORLD, info)
      end if
#endif

      if (MPI_Process % isRoot) then
         k = 1
         init_nclusters = nclusters
         do j = 1, init_nclusters

            ! Compute the inverse of the covariance matrix
            call matinv(ndims, g % cov(:,:,k), g % covinv(:,:,k), tol=tol, info=info)

            ! If a cluster collapses, limit its size
            if (info < 0) then
               do l = 1, ndims
                  g % cov(:,l,k) = 0.0_RP
                  g % cov(l,l,k) = tol**(1.0_RP/ndims) * 1e1_RP
               end do
               call matinv(ndims, g % cov(:,:,k), g % covinv(:,:,k), tol=tol, info=info)
            end if

            ! Also look for overlapping clusters
            deleteit = .false.
            do l = 1, k-1
               if (all(AlmostEqual(g % mu(:,k), g % mu(:,l), tol))) then
                  deleteit = .true.
                  exit
               end if
            end do

            ! Delete the marked clusters
            if (deleteit) then
               do l = k+1, nclusters
                  g % tau(l-1) = g % tau(l)
                  g % mu(:,l-1) = g % mu(:,l)
                  g % cov(:,:,l-1) = g % cov(:,:,l)
                  g % covinv(:,:,l-1) = g % covinv(:,:,l)
               end do
               nclusters = nclusters - 1
            else
               k = k + 1
            end if

         end do
      end if
!
!     Syncronize
!     ----------
#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         call MPI_Bcast(g % storage, size(g % storage), MPI_DOUBLE, 0, MPI_COMM_WORLD, info)
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
      real(RP),             intent(in)    :: x(:,:)
      real(RP),             intent(out)   :: xavg(:,:)
      integer,              intent(inout) :: clusters(:)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, ierr
      real(RP) :: val, tmp
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
         val = 0.0_RP
         tmp = 0.0_RP
         do j = 1, nclusters
            tmp = GMM_pdf(ndims, x(:,i), g % mu(:,j), g % cov(:,:,j), g % covinv(:,:,j))
            if (tmp > val) then
               val = tmp
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
      real(RP) :: d
      logical  :: iszero
      integer  :: info_


      if (n == 1) then
          info_ = 1
          Ainv = 1.0_RP / A

      elseif (n == 2) then
          d = matdet(2, A)
          if (present(tol)) then
              iszero = AlmostEqual(d, 0.0_RP, tol)
          else
              iszero = AlmostEqual(d, 0.0_RP)
          end if
          if (iszero) then
              info_ = -1
          else
              info_ = 1
              Ainv(1,1) =  A(2,2) / d
              Ainv(2,1) = -A(2,1) / d
              Ainv(1,2) = -A(1,2) / d
              Ainv(2,2) =  A(1,1) / d
          end if

      else
          info_ = -2

      end if

      if (present(info)) info = info_

   end subroutine matinv

end module Clustering
