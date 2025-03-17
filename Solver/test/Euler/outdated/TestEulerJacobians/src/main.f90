!
!//////////////////////////////////////////////////////
!
! Small program to test that the Euler Jacobians are working
program TestEulerJacobians
   use SMConstants           , only: RP, NDIM
   use PhysicsStorage_NS     , only: ConstructPhysicsStorage_NS, DestructPhysicsStorage_NS, NCONS, STANDARD_SPLIT, RIEMANN_LXF
   use Physics_NS            , only: EulerFlux0D, InviscidJacobian
   use RiemannSolvers_NS     , only: RiemannSolver, RiemannSolver_dFdQ, SetRiemannSolver
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities             , only: AlmostEqual
   implicit none
   !
   real(kind=RP) :: QL(NCONS)
   real(kind=RP) :: QR(NCONS)
   real(kind=RP) :: nHat(NDIM), t1(NDIM), t2(NDIM)
   real(kind=RP) :: Lref, timeref
   type(FTValueDictionary) :: controlVariables
   logical :: success
   
   
   real(kind=RP) :: dFdQL(NCONS,NCONS)
   real(kind=RP) :: dFdQR(NCONS,NCONS)
   
   !
   
!
!  Default values
!  --------------
   Lref = 1.0_RP
   timeref = 1.0_RP
   
   call controlVariables % initWithSize(16)
   call controlVariables % addValueForKey(0.7_RP,'mach number')
   call controlVariables % addValueForKey('euler','flow equations')
!~   call controlVariables % addValueForKey('lax-friedrichs','riemann solver')
   
!
!  Construct physics
!  -----------------
   call ConstructPhysicsStorage_NS( controlVariables, Lref, timeref, success )
   call SetRiemannSolver(RIEMANN_LXF, STANDARD_SPLIT)
!
!  
!  -----------------
   QL = [1.5_RP, 1._RP, -0.3_RP, -0.5_RP , 44._RP]
   QR = [5._RP , 0.8_RP,-0.2_RP, -0.55_RP, 33._RP]
   
   nHat = [1._RP,1._RP,1._RP]
   nHat = nHat / norm2(nHat)
   
   call GetTangentVectors(nHat,t1,t2)

!
!  Test the flux Jacobians
!  -----------------------
   
   call TestFluxJacobians(QL)
!
!  Test the numerical flux Jacobians
!  ---------------------------------
   
   call TestNumFluxJacobians(QL,QR,nHat,t1,t2)
   
   print*, 'getting numerical'
   call NumFluxNumJacobians(QL,QR,nHat,t1,t2,dFdQL,dFdQR)
   
   print*, 'Numerical dFdQL'
   call PrintMat(dFdQL)
   print*, 'Numerical dFdQR'
   call PrintMat(dFdQR)
   
   call RiemannSolver_dFdQ(QL,QR,nHat,dFdQL,1)
   call RiemannSolver_dFdQ(QL,QR,nHat,dFdQR,2)
   print*, 'Analytical dFdQL'
   call PrintMat(dFdQL)
   print*, 'Analytical dFdQR'
   call PrintMat(dFdQR)
   
   
!
!  Finish up
!  ---------
   
   call DestructPhysicsStorage_NS
   
contains

   subroutine TestFluxJacobians(Q)
      implicit none
      real(kind=RP), intent(in) :: Q(NCONS)
      !
      real(kind=RP) :: dFdQnum(NCONS,NCONS), dFdQan(NCONS,NCONS)
      real(kind=RP) :: dGdQnum(NCONS,NCONS), dGdQan(NCONS,NCONS)
      real(kind=RP) :: dHdQnum(NCONS,NCONS), dHdQan(NCONS,NCONS)
      !
      
      call FluxNumJacobians   (Q,dFdQnum,dGdQnum,dHdQnum)
      call InviscidJacobian   (Q,dFdQan ,dGdQan ,dHdQan )
      
      print*, maxval(abs(dFdQnum-dFdQan)), maxloc(abs(dFdQnum-dFdQan))
      print*, maxval(abs(dGdQnum-dGdQan)), maxloc(abs(dGdQnum-dGdQan))
      print*, maxval(abs(dHdQnum-dHdQan)), maxloc(abs(dHdQnum-dHdQan))
      
   end subroutine TestFluxJacobians
   
   subroutine TestNumFluxJacobians(QL,QR,nHat,t1,t2)
      implicit none
      !
      real(kind=RP), intent(in) :: QL(NCONS)
      real(kind=RP), intent(in) :: QR(NCONS)
      real(kind=RP), intent(in) :: nHat(NDIM),t1(NDIM),t2(NDIM)
      !
      real(kind=RP) :: dFdQLnum(NCONS,NCONS), dFdQLan(NCONS,NCONS)
      real(kind=RP) :: dFdQRnum(NCONS,NCONS), dFdQRan(NCONS,NCONS)
      !
      
      call NumFluxNumJacobians(QL,QR,nHat,t1,t2,dFdQLnum,dFdQRnum)
      call RiemannSolver_dFdQ(QL,QR,nHat,dFdQLan,1)
      call RiemannSolver_dFdQ(QL,QR,nHat,dFdQRan,2)
      
      print*, maxval(abs(dFdQLnum-dFdQLan)), maxloc(abs(dFdQLnum-dFdQLan))
      print*, maxval(abs(dFdQRnum-dFdQRan)), maxloc(abs(dFdQRnum-dFdQRan))
      
   end subroutine TestNumFluxJacobians
   
   subroutine FluxNumJacobians(Q_in,dFdQ,dGdQ,dHdQ)
      implicit none
      !
      real(kind=RP), intent(in) :: Q_in(NCONS)
      real(kind=RP) :: dFdQ(NCONS,NCONS)
      real(kind=RP) :: dGdQ(NCONS,NCONS)
      real(kind=RP) :: dHdQ(NCONS,NCONS)
      !
      integer       :: i
      real(kind=RP) :: buffer
      real(kind=RP) :: q (NCONS)
      real(kind=RP) :: F    (1:NCONS , 1:NDIM)
      real(kind=RP) :: Fbase(1:NCONS , 1:NDIM)
      real(kind=RP), parameter :: eps = 1.e-8_RP
      !
      
      q = Q_in
      call EulerFlux0D(Q, Fbase)
      
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call EulerFlux0D(Q, F)
         
         dFdQ(:,i) = (F(:,1)-Fbase(:,1))/eps
         dGdQ(:,i) = (F(:,2)-Fbase(:,2))/eps
         dHdQ(:,i) = (F(:,3)-Fbase(:,3))/eps
         
         q(i) = buffer
      end do
   end subroutine FluxNumJacobians
   
   subroutine NumFluxNumJacobians(QL_in,QR_in,nHat,t1,t2,dFdQL,dFdQR)
      implicit none
      !
      real(kind=RP), intent(in) :: QL_in(NCONS)
      real(kind=RP), intent(in) :: QR_in(NCONS)
      real(kind=RP), intent(in) :: nHat(NDIM),t1(NDIM),t2(NDIM)
      real(kind=RP) :: dFdQL(NCONS,NCONS)
      real(kind=RP) :: dFdQR(NCONS,NCONS)
      !
      integer       :: i
      real(kind=RP) :: buffer
      real(kind=RP) :: q (NCONS)
      real(kind=RP) :: F    (1:NCONS)
      real(kind=RP) :: Fbase(1:NCONS)
      real(kind=RP), parameter :: eps = 1.e-8_RP
      !
      
      call RiemannSolver( QL_in, QR_in, nHat, t1, t2, Fbase )
      
      ! Left
      q = QL_in
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call RiemannSolver( q, QR_in, nHat, t1, t2, F ) 
         
         dFdQL(:,i) = (F-Fbase)/eps
         
         q(i) = buffer
      end do
      
      !Right
      q = QR_in
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call RiemannSolver( QL_in, q, nHat, t1, t2, F ) 
         
         dFdQR(:,i) = (F-Fbase)/eps
         
         q(i) = buffer
      end do
   end subroutine NumFluxNumJacobians
   
   subroutine PrintMat(A)
      implicit none
      real(kind=RP) :: A(:,:)
      
      integer :: i
      
      do i=1, size(A,1)
         print*, A(i,:)
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