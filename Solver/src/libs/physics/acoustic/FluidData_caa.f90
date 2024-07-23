#include "Includes.h"
module FluidData_NS
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
      character(len=STR_LEN_FLUIDDATA) :: fluidName   ! Fluid name
      real(kind=RP)                    :: R           ! Gas constant
      real(kind=RP)                    :: gamma       ! Heat ratio
      real(kind=RP)                    :: sqrtGamma
      real(kind=RP)                    :: gammaMinus1
      real(kind=RP)                    :: gammaMinus1Div2
      real(kind=RP)                    :: gammaPlus1Div2
      real(kind=RP)                    :: gammaMinus1Div2sg
      real(kind=RP)                    :: gammaMinus1Div2g
      real(kind=RP)                    :: InvGammaPlus1Div2
      real(kind=RP)                    :: InvGammaMinus1
      real(kind=RP)                    :: InvGamma
      real(kind=RP)                    :: gammaDivGammaMinus1
      real(kind=RP)                    :: cp          ! R * gogm1
      real(kind=RP)                    :: cv          ! R * invgm1
   end type Thermodynamics_t


   type RefValues_t
      real(kind=RP)        :: V
      real(kind=RP)        :: T
      real(kind=RP)        :: p
      real(kind=RP)        :: rho
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP)        :: cp
      real(kind=RP)        :: cv
      real(kind=RP)        :: Mach
      real(kind=RP)        :: gammaM2
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

         thermodynamics % fluidName           = trim(thermodynamics_ % fluidName)
         thermodynamics % R                   = thermodynamics_ % R
         thermodynamics % gamma               = thermodynamics_ % gamma
         thermodynamics % sqrtGamma           = thermodynamics_ % sqrtGamma
         thermodynamics % gammaMinus1         = thermodynamics_ % gammaMinus1
         thermodynamics % gammaMinus1Div2     = thermodynamics_ % gammaMinus1Div2
         thermodynamics % gammaPlus1Div2      = thermodynamics_ % gammaPlus1Div2
         thermodynamics % gammaMinus1Div2sg   = thermodynamics_ % gammaMinus1Div2sg
         thermodynamics % gammaMinus1Div2g    = thermodynamics_ % gammaMinus1Div2g
         thermodynamics % InvGammaPlus1Div2   = thermodynamics_ % InvGammaPlus1Div2
         thermodynamics % InvGammaMinus1      = thermodynamics_ % InvGammaMinus1
         thermodynamics % InvGamma            = thermodynamics_ % InvGamma
         thermodynamics % gammaDivGammaMinus1 = thermodynamics_ % gammaDivGammaMinus1
         thermodynamics % cp                  = thermodynamics_ % cp
         thermodynamics % cv                  = thermodynamics_ % cv

      end subroutine SetThermodynamics

      subroutine SetRefValues( refValues_ )
         implicit none
         type(RefValues_t),   intent(in)     :: refValues_

         refValues % T        = refValues_ % T
         refValues % p        = refValues_ % p
         refValues % rho      = refValues_ % rho
         refValues % V        = refValues_ % V

      end subroutine SetRefValues

      subroutine SetDimensionless( dimensionless_ )
         implicit none
         type(Dimensionless_t),  intent(in)     :: dimensionless_

         dimensionless % cp              = dimensionless_ % cp
         dimensionless % cv              = dimensionless_ % cv
         dimensionless % Mach            = dimensionless_ % Mach
         dimensionless % gammaM2         = dimensionless_ % gammaM2

      end subroutine SetDimensionless
!
end module FluidData_NS
