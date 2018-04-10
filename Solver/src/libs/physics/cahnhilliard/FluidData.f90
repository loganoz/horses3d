#include "Includes.h"
module FluidData
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
      real(kind=RP)  :: M        ! Mobility
      real(kind=RP)  :: rhoS     ! Double-well function height
      real(kind=RP)  :: kappa    ! Gradient energy coefficient
      real(kind=RP)  :: c_alpha  ! Alpha equilibrium concentration
      real(kind=RP)  :: c_beta   ! Beta equilibrium concentration
      real(kind=RP)  :: thetaw   ! Wall angle
   end type Thermodynamics_t

   type RefValues_t
      real(kind=RP)        :: L     ! Reference length
      real(kind=RP)        :: time  ! Reference time
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP)  :: eps         ! Coefficient in the dimensionless CH equation
      real(kind=RP)  :: w           ! Dimensionless interface width: 7.071 sqrt(kappa/rhoS)/Lref
      real(kind=RP)  :: sigma       ! Interface energy: 0.01508 sqrt(kappa*rhoS)/Lref
   end type Dimensionless_t
!
!  ---------
!  Instances: they are public and can be accessed from the outside.
!  ---------
!
   type(Thermodynamics_t), protected   :: thermodynamics
   type(RefValues_t),      protected   :: refValues
   type(Dimensionless_t),  protected   :: dimensionless

   contains
!
!///////////////////////////////////////////////////////////////////////
!
!        This Set* subroutines are required since the intances defined
!     here are protected. The protected argument implies that they
!     can not be modified outside this module (these quantities are
!     really important to risk accidental modification).
!
!///////////////////////////////////////////////////////////////////////
!
      subroutine SetThermodynamics( thermodynamics_ )
         implicit none
         type(Thermodynamics_t), intent(in)  :: thermodynamics_

         thermodynamics % M       = thermodynamics_ % M
         thermodynamics % rhoS    = thermodynamics_ % rhoS
         thermodynamics % kappa   = thermodynamics_ % kappa
         thermodynamics % c_alpha = thermodynamics_ % c_alpha
         thermodynamics % c_beta  = thermodynamics_ % c_beta
         thermodynamics % thetaw  = thermodynamics_ % thetaw

      end subroutine SetThermodynamics

      subroutine SetRefValues( refValues_ )
         implicit none
         type(RefValues_t),   intent(in)     :: refValues_

         refValues % L        = refValues_ % L
         refValues % time     = refValues_ % time

      end subroutine SetRefValues

      subroutine SetDimensionless( dimensionless_ )
         implicit none
         type(Dimensionless_t),  intent(in)     :: dimensionless_

         dimensionless % w     = dimensionless_ % w
         dimensionless % eps   = dimensionless_ % eps
         dimensionless % sigma = dimensionless_ % sigma

      end subroutine SetDimensionless

end module FluidData
