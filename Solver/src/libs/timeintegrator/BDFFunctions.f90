module BDFFunctions
   use SMConstants
   implicit none
   private
   public BDF_SetOrder, BDF_SetPreviousSolution, BDF_MatrixShift, BDF_GetRHS, bdf_order
   
!
!  BDF coefficients for constant time-step
!  ---------------------------------------
   real(kind=RP), parameter :: BDFCoeff(6,5) = &
!                    a_1             a_2     a_3           a_4             a_5          a_6
         reshape( (/ 1.0_RP        , -1._RP, 0._RP       , 0._RP         , 0._RP      , 0._RP        ,  &   ! BDF1
                     1.5_RP        , -2._RP, 0.5_RP      , 0._RP         , 0._RP      , 0._RP        ,  &   ! BDF2
                     11._RP/6_RP   , -3._RP, 3._RP/2._RP , -1._RP/3._RP  , 0._RP      , 0._RP        ,  &   ! BDF3
                     25._RP/12_RP  , -4._RP, 3._RP       , -4._RP/3._RP  , 1._RP/4._RP, 0._RP        ,  &   ! BDF4
                     137._RP/60_RP , -5._RP, 5._RP       , -10._RP/3._RP , 5._RP/4._RP, -1._RP/5._RP /) &   ! BDF5
                                                                                                      , (/6,5/) )
   integer, parameter :: MAX_ORDER_CONS_DT = 5
!
!  Module variables
!  ----------------

   integer :: bdf_order       ! BDF order specified by user
   integer :: order           ! BDF order to be used
   integer :: StepsTaken = 0
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine BDF_SetOrder(order)
      implicit none
      !------------------------------------------------------
      integer, intent(in) :: order
      !------------------------------------------------------
      
      if (order > MAX_ORDER_CONS_DT) then
         write(STD_OUT,*) 'WARNING :: Maximum BDF order for constant time-step is 5. Using 1 by default.'
         bdf_order = 1
      else
         bdf_order = order
      end if
      
   end subroutine BDF_SetOrder
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine BDF_SetPreviousSolution(PrevQ,Q)
      implicit none
      !------------------------------------------------------
      real(kind=RP) :: PrevQ(:,:)
      real(kind=RP) :: Q(:)
      !------------------------------------------------------
      integer :: i      ! Counter
      !------------------------------------------------------
      
      StepsTaken = StepsTaken + 1
    
      order = min(StepsTaken, bdf_order)
      
      do i=Order, 2, -1
         PrevQ(:,i) = PrevQ(:,i-1)
      end do
      
      PrevQ(:,1) = Q
      
   end subroutine BDF_SetPreviousSolution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function BDF_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = -BDFCoeff(1,order)/dt
      
   end function BDF_MatrixShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function BDF_GetRHS(Q,PrevQ,Qdot,dt) result(RHS)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in)  :: Q(:)
      real(kind=RP), intent(in)  :: PrevQ(:,:)
      real(kind=RP), intent(in)  :: Qdot(:)
      real(kind=RP), intent(in)  :: dt
      real(kind=RP)              :: RHS(size(Q))
      !------------------------------------------------------
      integer :: k
      real(kind=RP) :: invdt
      !------------------------------------------------------
      
      invdt = 1._RP/dt
      
      RHS = Q * BDFCoeff(1,order)*invdt - Qdot
      
      do k=1, order
         RHS = RHS + BDFCoeff(k+1,order) * PrevQ(:,k) * invdt
      end do
      
   end function BDF_GetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!


!~   subroutine BDFCoefficientsVariable_dt(order,dt)
!~      implicit none
!~      integer      , intent(in) :: order
!~      real(kind=RP), intent(in) :: dt(order)
      
!~      select case(order)
!~         case (1)
!~            a(1) = 1.0_RP / dt(1)
!~            a(2) = -1.0_RP / dt(1)
!~         case (2)
!~            a(1) = a(1) + 1.0_RP / (dt(1)+dt(2)) 
!~            a(2) = a(2) - (1.0_RP + dt(1)/dt(2)) / (dt(1)+dt(2)) 
!~            a(3) = (dt(1)/dt(2)) / (dt(1)+dt(2)) 
!~         case (3)
!~            a(1) = a(1) + 1.0_RP / (dt(1)+dt(2)+dt(3)) 
!~            a(2) = a(2) - (1.0_RP + dt(1)/dt(2)*(1.0+(dt(1)+dt(2))/(dt(2)+dt(3)))) / (dt(1)+dt(2)+dt(3)) 
!~            a(3) = a(3) + (dt(1)/dt(2)*(1.0+(dt(1)+dt(2))/(dt(2)+dt(3))) + &
!~                     dt(1)/dt(3)*(dt(1)+dt(2))/(dt(2)+dt(3)) ) / (dt(1)+dt(2)+dt(3)) 
!~            a(4) = -(dt(1)/dt(3))*(dt(1)+dt(2))/(dt(2)+dt(3)) / (dt(1)+dt(2)+dt(3))
!~         case default
!~            ERROR stop ':: variable time-step BDF only up to order 3'
!~      end select
!~   end subroutine
end module BDFFunctions
