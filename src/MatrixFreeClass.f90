!
!//////////////////////////////////////////////////////
!
!
!   Matrix-free operations.
!
!//////////////////////////////////////////////////////
!
module MatrixFreeClass
   USE SMConstants             
   use GenericMatrixClass   , only: Matrix_t, DenseBlock_t
   use LinkedListMatrixClass
   use JacobianDefinitions  , only: JACEPS
   use DGSEMClass
   use PhysicsStorage
#include "Includes.h"
   implicit none
   !-----------------------------------------------------------------------------   
   private
   public  :: MF_JacVecMul, MF_p_F
   
!
!========
 contains
!========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !  ----------------------------------------------------
   !  MF_JacVecMul:
   !  Matrix vector product (v = Au) 
   !  ----------------------------------------------------
   subroutine MF_JacVecMul(p_sem, DimPrb, Ur, F_Ur, x, Ax, dt, timesolve, shift, ComputeTimeDerivative)
      implicit none
      !-arguments-----------------------------------------------
      type(DGSem),               intent(inout) :: p_sem
      integer,                   intent(in)    :: DimPrb
      real(kind=RP),             intent(in)    :: x    (DimPrb)
      real(kind=RP),             intent(in)    :: Ur   (DimPrb)
      real(kind=RP),             intent(in)    :: F_Ur (DimPrb)
      real(kind=RP),             intent(in)    :: timesolve
      real(kind=RP),             intent(in)    :: dt
      real(kind=RP),             intent(in)    :: shift
      real(kind=RP),             intent(out)   :: Ax   (DimPrb)
      procedure(ComputeTimeDerivative_f)   :: ComputeTimeDerivative
      !-local-variables-----------------------------------------
      ! real(kind=RP) :: eps, shift
      real(kind=RP) :: eps
      real(kind=RP) :: L2x
      !---------------------------------------------------------
      
      L2x = norm2(x)

      if (L2x < 1e-12) then 
         Ax = 0.0_RP
      else 
         !eps = 1e-11_RP * (1._RP + norm2(p_sem % mesh % storage % Q) ) / norm2(x)
         eps = 1e-8_RP * (1._RP + norm2(x) )
         ! eps = 1e-2_RP * (1._RP + norm2(x) )

         Ax = ( MF_p_F(p_sem, DimPrb, Ur + x * eps, dt + timesolve, ComputeTimeDerivative) - F_Ur)/ eps + shift * x  ! First Order 

         !Ax = ( MF_p_F(p_sem, DimPrb, Ur + x * eps, dt + timesolve, ComputeTimeDerivative) - MF_p_F(p_sem, DimPrb, Ur, dt + timesolve, ComputeTimeDerivative))/ eps + shift * x  ! First Order 

         !Ax = ( MF_p_F(p_sem,DimPrb,Ur + x * eps, dt+timesolve,ComputeTimeDerivative) - MF_p_F(p_sem,DimPrb,Ur - x * eps,dt+timesolve, ComputeTimeDerivative))  /(2._RP * eps)  + shift*x   !Second order

      end if
   end subroutine MF_JacVecMul
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Function to return the time derivative using a specific solution vector
!     ---------------------------------------------------
      function MF_p_F(p_sem, DimPrb, u, CTD_time, computeTimeDerivative) result(F)
         implicit none
         !-arguments-----------------------------------------------
         type(DGSem),           intent(inout) :: p_sem
         integer,               intent(in)    :: DimPrb
         real(kind=rp),         intent(in)    :: u(DimPrb)
         real(kind=RP),         intent(in)    :: CTD_time
         procedure(ComputeTimeDerivative_f)   :: ComputeTimeDerivative
         real(kind=rp)                        :: F(DimPrb)
         !-local-variables-----------------------------------------
         real(kind=rp)                        :: u_p(DimPrb)
         !---------------------------------------------------------
         
         ! Save original Q
         u_p = p_sem % mesh % storage % Q
         
         ! Obtain derivative with new Q
         p_sem % mesh % storage % Q = u
         call p_sem % mesh % storage % global2LocalQ
         call ComputeTimeDerivative(p_sem % mesh, p_sem % particles, CTD_time, CTD_IGNORE_MODE)
         call p_sem % mesh % storage % local2GlobalQdot(p_sem % NDOF)
         
         F = p_sem % mesh % storage % Qdot

         ! Restore original Q
         p_sem % mesh % storage % Q = u_p   ! TODO: this step can be avoided if Ur is not read in the "child" GMRES (the preconditioner)
         call p_sem % mesh % storage % global2LocalQ
      end function MF_p_F
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MatrixFreeClass