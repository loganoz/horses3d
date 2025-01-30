#include "Includes.h"
module FluidData_MU
   use SMConstants
   implicit none

   private
   public Thermodynamics_t , thermodynamics , SetThermodynamics
   public Dimensionless_t  , dimensionless  , SetDimensionless
   public RefValues_t      , refValues      , SetRefValues

   integer, parameter   :: STR_LEN_FLUIDDATA = 128
!
!  ----------------
!  Type definitions
!  ----------------
!
   type Thermodynamics_t
      real(kind=RP) :: rho(2)
      real(kind=RP) :: mu(2)
      real(kind=RP) :: c02
   end type Thermodynamics_t

   type RefValues_t
      real(kind=RP) :: rho
      real(kind=RP) :: V
      real(kind=RP) :: p
      real(kind=RP) :: mu
      real(kind=RP) :: g0
      real(kind=RP) :: AoAPhi 
      real(kind=RP) :: AoATheta
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP) :: rho(2)
      real(kind=RP) :: mu(2)
      real(kind=RP) :: Re(2)
      real(kind=RP) :: Fr
      real(kind=RP) :: gravity_dir(NDIM)
      real(kind=RP) :: invFr2
      real(kind=RP) :: vel_dir(NDIM)
      real(kind=RP) :: rho_max
      real(kind=RP) :: rho_min
      real(kind=RP) :: Ma2
      real(kind=RP) :: invMa2
   end type Dimensionless_t
!
!  ---------
!  Instances: they are public and can be accessed from the outside.
!  ---------
!
   type(Thermodynamics_t), protected   :: thermodynamics
   !$acc declare create(thermodynamics)
   type(RefValues_t),      protected   :: refValues
   !$acc declare create(refValues)
   type(Dimensionless_t),  protected   :: dimensionless
   !$acc declare create(dimensionless)

   contains
!
!///////////////////////////////////////////////////////////////////////
!
!        This Set* subroutines are required since the instances defined
!     here are protected. The protected argument implies that they
!     can not be modified outside this module (these quantities are
!     really important to risk accidental modification).
!
!///////////////////////////////////////////////////////////////////////
!
      subroutine SetThermodynamics( thermodynamics_ )
         implicit none
         type(Thermodynamics_t), intent(in)  :: thermodynamics_

         thermodynamics % rho = thermodynamics_ % rho
         thermodynamics % mu  = thermodynamics_ % mu
         thermodynamics % c02 = thermodynamics_ % c02

         !$acc update device(thermodynamics, thermodynamics % rho, thermodynamics % mu, thermodynamics % c02)

      end subroutine SetThermodynamics

      subroutine SetRefValues( refValues_ )
         implicit none
         type(RefValues_t),   intent(in)     :: refValues_

         refValues % p   = refValues_ % p
         refValues % rho = refValues_ % rho
         refValues % V   = refValues_ % V
         refValues % mu  = refValues_ % mu
         refValues % g0  = refValues_ % g0
         refValues % AoAPhi  = refValues_ % AoAPhi
         refValues % AoATheta  = refValues_ % AoATheta

         !$acc update device(refValues, refValues % p, refValues % rho, refValues % V, refValues % mu, refValues % g0)
         !$acc update device(refValues % AoAPhi, refValues % AoATheta)

      end subroutine SetRefValues

      subroutine SetDimensionless( dimensionless_ )
         implicit none
         type(Dimensionless_t),  intent(in)     :: dimensionless_

         dimensionless % rho         = dimensionless_ % rho
         dimensionless % mu          = dimensionless_ % mu
         dimensionless % Re          = dimensionless_ % Re
         dimensionless % Fr          = dimensionless_ % Fr
         dimensionless % invFr2      = dimensionless_ % invFr2
         dimensionless % gravity_dir = dimensionless_ % gravity_dir
         dimensionless % vel_dir     = dimensionless_ % vel_dir
         dimensionless % rho_max     = dimensionless_ % rho_max
         dimensionless % rho_min     = dimensionless_ % rho_min
         dimensionless % Ma2         = dimensionless_ % Ma2
         dimensionless % invMa2      = dimensionless_ % invMa2

         !$acc update device(dimensionless, dimensionless % rho, dimensionless % mu, dimensionless % Re)
         !$acc update device(dimensionless % Fr, dimensionless % invFr2, dimensionless % gravity_dir)
         !$acc update device(dimensionless % vel_dir, dimensionless % rho_max, dimensionless % rho_min)
         !$acc update device(dimensionless % Ma2, dimensionless % invMa2)
         
      end subroutine SetDimensionless
end module FluidData_MU