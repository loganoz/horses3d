module variableConversion
   use SMConstants
   use FluidData
   implicit none   

   private
   public GradientValuesForQ, pressure, temperature

   interface GradientValuesForQ
       module procedure GradientValuesForQ_0D , GradientValuesForQ_3D
   end interface GradientValuesForQ

   contains
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      SUBROUTINE GradientValuesForQ_0D( Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(5)     , INTENT(IN)  :: Q
      REAL(KIND=RP), DIMENSION(4), INTENT(OUT) :: U
!
!     ---------------
!     Local Variables
!     ---------------
!     
      U(1) = Q(1)
      U(2) = Q(2)
      U(3) = Q(3)
      U(4) = Q(4)

      END SUBROUTINE GradientValuesForQ_0D

      SUBROUTINE GradientValuesForQ_3D( Nx, Ny, Nz, Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,       intent(in)  :: Nx, Ny, Nz
      REAL(KIND=RP), INTENT(IN)  :: Q(1:5,0:Nx,0:Ny,0:Nz)
      REAL(KIND=RP), INTENT(OUT) :: U(1:4,0:Nx,0:Ny,0:Nz)
!
!     ---------------
!     Local Variables
!     ---------------
!     
      U(1,:,:,:) = Q(1,:,:,:)
      U(2,:,:,:) = Q(2,:,:,:)
      U(3,:,:,:) = Q(3,:,:,:)
      U(4,:,:,:) = Q(4,:,:,:)

      END SUBROUTINE GradientValuesForQ_3D
!
!
! /////////////////////////////////////////////////////////////////////
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
      PURE FUNCTION Pressure(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(5), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = 1.0_RP

      END FUNCTION Pressure

      PURE FUNCTION Temperature(Q) RESULT(T)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(5), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: T
!
      T = 1.0_RP

      END FUNCTION Temperature

end module variableConversion
