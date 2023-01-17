
!//////////////////////////////////////////////////////
!
! Small program to test that the NS Jacobians are working
program TestEulerJacobians
   use SMConstants           , only: RP, NDIM
   use PhysicsStorage_NS     , only: ConstructPhysicsStorage_NS, DestructPhysicsStorage_NS, NCONS, STANDARD_SPLIT, RIEMANN_LXF
   use Physics_NS            , only: ViscousJacobian, ViscousFlux0D
   use FluidData_NS          , only: dimensionless, thermodynamics
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities             , only: AlmostEqual
   implicit none
   !
   real(kind=RP) :: Q(NCONS)
   real(kind=RP) :: Q_x(NCONS)
   real(kind=RP) :: Q_y(NCONS)
   real(kind=RP) :: Q_z(NCONS)
   real(kind=RP) :: Lref, timeref
   type(FTValueDictionary) :: controlVariables
   logical :: success
   
   integer :: i,j
   real(kind=RP) :: dFdQ_an    (NCONS,NCONS,NDIM)
   real(kind=RP) :: dFdgradQ_an(NCONS,NCONS,NDIM,NDIM)
   real(kind=RP) :: dFdQ_nu    (NCONS,NCONS,NDIM)
   real(kind=RP) :: dFdgradQ_nu(NCONS,NCONS,NDIM,NDIM)
   
   !
   
!
!  Default values
!  --------------
   Lref = 1.0_RP
   timeref = 1.0_RP
   
   call controlVariables % initWithSize(16)
   call controlVariables % addValueForKey(6._RP,'mach number')
   call controlVariables % addValueForKey(1._RP,'reynolds number')
   call controlVariables % addValueForKey('NS','flow equations')
   
!
!  Construct physics
!  -----------------
   call ConstructPhysicsStorage_NS( controlVariables, Lref, timeref, success )
!
!  
!  -----------------
   Q   = [0.95_RP, 1._RP, 1.2_RP, -0.5_RP , 44._RP]
   Q_x = [0.1_RP, 0.5_RP, 0.3_RP,  0.7_RP , 1._RP]
   Q_y = [0.8_RP,-0.2_RP, 0.5_RP,  0.4_RP , 3.2_RP]
   Q_z = [0.2_RP, 0.4_RP, 0.6_RP,  0.6_RP , -2.05_RP]

!
!  Test the flux Jacobians
!  -----------------------
   
   call ViscousJacobian (Q, Q_x, Q_y, Q_z, dFdgradQ_an, dFdQ_an)
   call FluxNumJacobians(Q, Q_x, Q_y, Q_z, dFdgradQ_nu, dFdQ_nu)
   
   print*, '* dF1/dQ an'
   call PrintMat(dFdQ_an(:,:,1))
   print*, '* dF1/dQ nu'
   call PrintMat(dFdQ_nu(:,:,1))
   
   print*, '* dF2/dQ an'
   call PrintMat(dFdQ_an(:,:,2))
   print*, '* dF2/dQ nu'
   call PrintMat(dFdQ_nu(:,:,2))
   
   print*, '* dF3/dQ an'
   call PrintMat(dFdQ_an(:,:,3))
   print*, '* dF3/dQ nu'
   call PrintMat(dFdQ_nu(:,:,3))
   
   
   print*, '-----------'
   do i=1, NDIM
      print*, '* dF',i,'/dQ'
      print*,  maxval( abs( dFdQ_an(:,:,i) - dFdQ_nu(:,:,i) ) ), maxloc( abs( dFdQ_an(:,:,i) - dFdQ_nu(:,:,i) ) )
   end do
   print*, '-----------'
   do j=1,NDIM ; do i=1, NDIM
      print*, '* dF',j,'/dQ_',i
      print*,  maxval( abs( dFdgradQ_an(:,:,i,j) - dFdgradQ_nu(:,:,i,j) ) ), &
               maxloc( abs( dFdgradQ_an(:,:,i,j) - dFdgradQ_nu(:,:,i,j) ) )
   end do ; end do
   
!
!  Test the numerical flux Jacobians
!  ---------------------------------
   
   
   
!
!  Finish up
!  ---------
   
   call DestructPhysicsStorage_NS
   
contains

!~   subroutine TestFluxJacobians(Q)
!~      implicit none
!~      real(kind=RP), intent(in) :: Q(NCONS)
!~      !
!~      real(kind=RP) :: dFdQnum(NCONS,NCONS), dFdQan(NCONS,NCONS)
!~      real(kind=RP) :: dGdQnum(NCONS,NCONS), dGdQan(NCONS,NCONS)
!~      real(kind=RP) :: dHdQnum(NCONS,NCONS), dHdQan(NCONS,NCONS)
!~      !
      
!~      call FluxNumJacobians   (Q,dFdQnum,dGdQnum,dHdQnum)
!~      call InviscidJacobian   (Q,dFdQan ,dGdQan ,dHdQan )
      
!~      print*, maxval(abs(dFdQnum-dFdQan)), maxloc(abs(dFdQnum-dFdQan))
!~      print*, maxval(abs(dGdQnum-dGdQan)), maxloc(abs(dGdQnum-dGdQan))
!~      print*, maxval(abs(dHdQnum-dHdQan)), maxloc(abs(dHdQnum-dHdQan))
      
!~   end subroutine TestFluxJacobians
   
   
   subroutine FluxNumJacobians(Q_in,Q_x_in,Q_y_in,Q_z_in, dFdgradQ, dFdQ)
      implicit none
      !
      real(kind=RP), intent(in) :: Q_in(NCONS)
      real(kind=RP), intent(in) :: Q_x_in(NCONS)
      real(kind=RP), intent(in) :: Q_y_in(NCONS)
      real(kind=RP), intent(in) :: Q_z_in(NCONS)
      real(kind=RP) :: dFdQ    (NCONS,NCONS,NDIM)
      real(kind=RP) :: dFdgradQ(NCONS,NCONS,NDIM,NDIM)
      !
      integer       :: i
      real(kind=RP) :: buffer
      real(kind=RP) :: q (NCONS)
      real(kind=RP) :: F    (1:NCONS , 1:NDIM)
      real(kind=RP) :: Fbase(1:NCONS , 1:NDIM)
      real(kind=RP), parameter :: eps = 1.e-8_RP
      real(kind=RP) :: mu, kappa, beta
      !
      
      mu    = dimensionless % mu
      kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * dimensionless % Mach**2 * dimensionless % Pr ) * dimensionless % mu
      beta  = 0._RP
      
      call ViscousFlux0D(5, 5, Q_in, Q_x_in, Q_y_in, Q_z_in, mu, beta, kappa, Fbase)
      
!     dF/dQ
!     -----
      q = Q_in
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call ViscousFlux0D(5, 5, Q, Q_x_in, Q_y_in, Q_z_in, mu, beta, kappa, F)
         
         dFdQ(:,i,1) = (F(:,1)-Fbase(:,1))/eps
         dFdQ(:,i,2) = (F(:,2)-Fbase(:,2))/eps
         dFdQ(:,i,3) = (F(:,3)-Fbase(:,3))/eps
         
         q(i) = buffer
      end do
      
!     dF/d Qx
!     -------
      q = Q_x_in
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call ViscousFlux0D(5, 5, Q_in, q, Q_y_in, Q_z_in, mu, beta, kappa, F)
         
         dFdgradQ(:,i,1,1) = (F(:,1)-Fbase(:,1))/eps
         dFdgradQ(:,i,1,2) = (F(:,2)-Fbase(:,2))/eps
         dFdgradQ(:,i,1,3) = (F(:,3)-Fbase(:,3))/eps
         
         q(i) = buffer
      end do
!     dF/d Qy
!     -------
      q = Q_y_in
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call ViscousFlux0D(5, 5, Q_in, Q_x_in, q, Q_z_in, mu, beta, kappa, F)
         
         dFdgradQ(:,i,2,1) = (F(:,1)-Fbase(:,1))/eps
         dFdgradQ(:,i,2,2) = (F(:,2)-Fbase(:,2))/eps
         dFdgradQ(:,i,2,3) = (F(:,3)-Fbase(:,3))/eps
         
         q(i) = buffer
      end do
      
!     dF/d Qz
!     -------
      q = Q_z_in
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call ViscousFlux0D(5, 5, Q_in, Q_x_in, Q_y_in, q, mu, beta, kappa, F)
         
         dFdgradQ(:,i,3,1) = (F(:,1)-Fbase(:,1))/eps
         dFdgradQ(:,i,3,2) = (F(:,2)-Fbase(:,2))/eps
         dFdgradQ(:,i,3,3) = (F(:,3)-Fbase(:,3))/eps
         
         q(i) = buffer
      end do
   end subroutine FluxNumJacobians
   
   
   subroutine PrintMat(A)
      implicit none
      real(kind=RP) :: A(:,:)
      
      integer :: i
      
      do i=1, size(A,1)
         write(*,'(5ES16.8)') A(i,:)
      end do
   end subroutine PrintMat
   
   ! Get a set of two valid tangent vectors
   subroutine GetTangentVectors(nHat,t1,t2)
      implicit none
      real(kind=RP), intent(in)  :: nHat(NDIM)
      real(kind=RP), intent(out) :: t1(NDIM)
      real(kind=RP), intent(out) :: t2(NDIM)
      !
      real(kind=RP) :: t1_ini(NDIM)
      real(kind=RP) :: t2_ini(NDIM)
      !
      
      ! Create set of vectors
      t1_ini = [1._RP,1._RP,-(nHat(1)+nHat(2))/nHat(3)]
      t2_ini = [1._RP,0._RP,-(nHat(1))/nHat(3)]
      
      !Normalize
      t1 = t1_ini / norm2(t1_ini)
      
      t2_ini = t2_ini - t1 * dot_product(t2_ini, t1)
      t2 = t2_ini / norm2(t2_ini)
      
   end subroutine GetTangentVectors
end program TestEulerJacobians
