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
   end type Thermodynamics_t


   type RefValues_t
      real(kind=RP)        :: L
      real(kind=RP)        :: time
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP)  :: eps = 0.01_RP
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

      end subroutine SetDimensionless

end module FluidData
