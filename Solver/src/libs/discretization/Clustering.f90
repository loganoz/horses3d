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

   type GaussianList_t
      integer               :: d              ! Space dimension
      integer               :: n              ! # of components
      real(RP), allocatable :: tau(:)         ! [n]: weights
      real(RP), allocatable :: mu(:,:)        ! [d,n]: centroids
      real(RP), allocatable :: cov(:,:,:)     ! [d,d,n]: covariance matrices
      real(RP), allocatable :: covinv(:,:,:)  ! [d,d,n]: inverse covariance matrices
   end type GaussianList_t

   public :: kMeans
   public :: GMM
!
!  ========
   contains
!  ========
!
   subroutine kMeans(nclusters, x, xavg, clusters, info)
!
!     ---------
!     Interface
!     ---------
      integer,           intent(in)    :: nclusters
      real(RP),          intent(in)    :: x(:,:)
      real(RP),          intent(inout) :: xavg(:,:)
      integer,           intent(out)   :: clusters(:)
      integer, optional, intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      integer, parameter   :: maxIters = 50
      integer, allocatable :: prevClusters(:)
      integer              :: ndims
      integer              :: npts
      integer              :: i


      ndims = size(x, dim=1)
      npts  = size(x, dim=2)
!
!     Initial clusters
!     ----------------
      call kMeans_compute_clusters(nclusters, ndims, npts, x, xavg, clusters)
!
!     Loop until convergence
!     ----------------------
      do i = 1, maxIters
         prevClusters = clusters
         call kMeans_compute_centroids(nclusters, ndims, npts, x, xavg, clusters)
         call kMeans_compute_clusters(nclusters, ndims, npts, x, xavg, clusters)
         if (all(prevClusters == clusters)) exit
      end do
!
!     Check convergence
!     -----------------
      if (present(info)) info = merge(-1, i, i > maxIters)
!
!     Sort clusters
!     -------------
      call kMeans_sort_clusters(nclusters, ndims, npts, xavg, clusters)

   end subroutine kMeans
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine kMeans_compute_clusters(nclusters, ndims, npts, x, xavg, clusters)
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
      real(RP) :: dist, minDist


      do i = 1, npts
         minDist = huge(1.0_RP)
         do j = 1, nclusters
            dist = norm2(x(:,i) - xavg(:,j))
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
      integer  :: i, j
      real(RP) :: normvec(nclusters)
      real(RP) :: xavg_tmp(ndims,nclusters)
      integer  :: indices(nclusters)
      integer  :: clustermap(nclusters)

!
!     Sort the clusters by their distance to the origin
!     -------------------------------------------------
      normvec = norm2(xavg, dim=1)
      indices = [(i, i = 1, nclusters)]
      call BubblesortWithFriend(normvec, indices)
      do i = 1, nclusters
         clustermap(indices(i)) = i
      end do
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
   subroutine GMM(init_nclusters, nclusters, x, xavg, clusters, info)
!
!     ---------
!     Interface
!     ---------
      integer,           intent(in)    :: init_nclusters
      integer,           intent(out)   :: nclusters
      real(RP),          intent(in)    :: x(:,:)
      real(RP),          intent(inout) :: xavg(:,:)
      integer,           intent(out)   :: clusters(:)
      integer, optional, intent(out)   :: info
!
!     ---------------
!     Local variables
!     ---------------
      real(RP), parameter :: tol = 1e-2_RP                     ! Magic!!
      real(RP), parameter :: small = 1e4_RP * epsilon(1.0_RP)  ! Magic!!
      integer,  parameter :: maxIters = 50

      type(GaussianList_t)  :: g, local_g
      integer               :: info_
      integer               :: ndims
      integer               :: npts
      integer               :: iter
      real(RP)              :: ll
      real(RP)              :: llprev
      real(RP), allocatable :: R(:,:)
      real(RP), allocatable :: minimum(:), maximum(:)


      nclusters = init_nclusters
!
!     Trivial case
!     ------------
      if (nclusters == 1) then
         clusters = 1
         minimum = minval(x, dim=2)
         maximum = maxval(x, dim=2)
         call MPI_MinMax(minimum, maximum)
         xavg(:,1) = (minimum + maximum) * 0.5_RP
         return
      end if

      ndims = size(x, dim=1)
      npts  = size(x, dim=2)

      call GMM_init(g, ndims, nclusters, xavg)
!
!     EM algorithm
!     ------------
      allocate(R(npts, nclusters))
      ll = huge(1.0_RP)
      do iter = 1, maxIters
         llprev = ll
         call GMM_Estep(ndims, npts, nclusters, g, x, R, small, ll)
         call GMM_Mstep(ndims, npts, nclusters, g, local_g, x, R, small)
         nclusters = g % n
         if (abs((ll - llprev) / ll) <= tol .or. nclusters < 1) exit
      end do
!
!     Check convergence
!     -----------------
      if (iter > maxIters) then
         call kMeans(nclusters, x, xavg, clusters, info_)
         if (info_ == -1) then
            info_ = 2
         end if
      else
         info_ = 0
      end if
      if (present(info)) info = info_
!
!     Sort clusters
!     -------------
      call GMM_sort_clusters(ndims, npts, nclusters, g, x, xavg, clusters)

   end subroutine GMM
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine GMM_init(g, ndims, nclusters, xavg)
!
!     ---------
!     Interface
!     ---------
      type(GaussianList_t), intent(out) :: g
      integer,              intent(in)  :: ndims
      integer,              intent(in)  :: nclusters
      real(RP),             intent(in)  :: xavg(ndims, nclusters)
!
!     ---------------
!     Local variables
!     ---------------
      integer :: i, j


      g % d = ndims
      g % n = nclusters
      allocate(g % tau(nclusters))
      allocate(g % mu(ndims, nclusters))
      allocate(g % cov(ndims, ndims, nclusters)); g % cov = 0.0_RP
      allocate(g % covinv, source=g % cov)

      g % tau = 1.0_RP / nclusters
      g % mu = xavg
      do i = 1, nclusters
         do j = 1, ndims
            g % cov(j,j,i) = 1.0_RP
            g % covinv(j,j,i) = 1.0_RP
         end do
      end do

   end subroutine GMM_init
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
      real(RP) :: den

!
!     E-step loop
!     -----------
      logL = 0.0_RP
      do i = 1, npts
         den = 0.0_RP
         do j = 1, nclusters
            R(i,j) = g % tau(j) * GMM_pdf(                                   &
               ndims, x(:,i), g % mu(:,j), g % cov(:,:,j), g % covinv(:,:,j) &
            )
            den = den + R(i,j)
         end do

         ! If the point is far from all clusters
         if (AlmostEqual(den, 0.0_RP)) then
             R(i,1:nclusters) = 1.0_RP / nclusters
             den = small
         else
             R(i,1:nclusters) = R(i,1:nclusters) / den
         end if
         logL = logL + log(den)
      end do
!
!     Global log-likelihood
!     ---------------------
      call MPI_SumAll(logL)

   end subroutine GMM_Estep
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine GMM_Mstep(ndims, npts, nclusters, g, local_g, x, R, tol)
!
!     ---------
!     Interface
!     ---------
      integer,              intent(in)    :: ndims
      integer,              intent(in)    :: npts
      integer,              intent(inout) :: nclusters
      type(GaussianList_t), intent(inout) :: g
      type(GaussianList_t), intent(inout) :: local_g
      real(RP),             intent(in)    :: x(:,:)
      real(RP),             intent(in)    :: R(:,:)
      real(RP),             intent(in)    :: tol
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, k, l, m
      integer  :: nptsAll
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
      if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
         local_g = g
         call MPI_reduce(local_g % tau, g % tau, size(g % tau), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, info)
         call MPI_reduce(local_g % mu, g % mu, size(g % mu), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, info)
         call MPI_reduce(local_g % cov, g % cov, size(g % cov), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, info)
#endif
      end if

      if (MPI_Process % isRoot) then
         k = 1
         do j = 1, nclusters

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
               do l = k+1, g % n
                  g % tau(l-1) = g % tau(l)
                  g % mu(:,l-1) = g % mu(:,l)
                  g % cov(:,:,l-1) = g % cov(:,:,l)
                  g % covinv(:,:,l-1) = g % covinv(:,:,l)
               end do
               g % n = g % n - 1
            else
               k = k + 1
            end if

         end do
      end if
!
!     Syncronize
!     ----------
      if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
         call MPI_Bcast(g % tau, size(g % tau), MPI_DOUBLE, 0, MPI_COMM_WORLD, info)
         call MPI_Bcast(g % mu, size(g % mu), MPI_DOUBLE, 0, MPI_COMM_WORLD, info)
         call MPI_Bcast(g % cov, size(g % cov), MPI_DOUBLE, 0, MPI_COMM_WORLD, info)
         call MPI_Bcast(g % covinv, size(g % covinv), MPI_DOUBLE, 0, MPI_COMM_WORLD, info)
#endif
      end if

   end subroutine GMM_Mstep
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine GMM_sort_clusters(ndims, npts, nclusters, g, x, xavg, clusters)
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
      integer  :: i, j
      real(RP) :: val, tmp
      real(RP) :: normvec(nclusters)
      integer  :: indices(nclusters)
      integer  :: clustermap(nclusters)

!
!     Sort the clusters by their distance to the origin
!     -------------------------------------------------
      normvec = norm2(g % mu(:,1:nclusters), dim=1)
      indices = [(i, i = 1, nclusters)]
      call BubblesortWithFriend(normvec, indices)
      do i = 1, nclusters
         clustermap(indices(i)) = i
      end do
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
          Ainv = A

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
