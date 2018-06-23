!
!//////////////////////////////////////////////////////
!
!   @File:    FluidData_iNS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:25 2018
!   @Last revision date: Sat Jun 23 10:20:34 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module FluidData_iNS
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
      real(kind=RP)                    :: rho0c02     ! Artificial compressibility const
   end type Thermodynamics_t


   type RefValues_t
      real(kind=RP)        :: V
      real(kind=RP)        :: p
      real(kind=RP)        :: rho
      real(kind=RP)        :: mu
      real(kind=RP)        :: AoATheta
      real(kind=RP)        :: AoAPhi
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP)        :: Re
      real(kind=RP)        :: Fr
      real(kind=RP)        :: mu
      real(kind=RP)        :: gammaM2
      real(kind=RP)        :: invFroudeSquare
      real(kind=RP)        :: gravity_dir(NDIM)
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

         thermodynamics % fluidName           = trim(thermodynamics_ % fluidName)
         thermodynamics % rho0c02             = thermodynamics_ % rho0c02

      end subroutine SetThermodynamics

      subroutine SetRefValues( refValues_ )
         implicit none
         type(RefValues_t),   intent(in)     :: refValues_

         refValues % p        = refValues_ % p
         refValues % rho      = refValues_ % rho
         refValues % V        = refValues_ % V
         refValues % mu       = refValues_ % mu
         refValues % AoATheta = refValues_ % AoATheta
         refValues % AoAPhi   = refValues_ % AoAPhi

      end subroutine SetRefValues

      subroutine SetDimensionless( dimensionless_ )
         implicit none
         type(Dimensionless_t),  intent(in)     :: dimensionless_

         dimensionless % Re              = dimensionless_ % Re
         dimensionless % Fr              = dimensionless_ % Fr
         dimensionless % mu              = dimensionless_ % mu
         dimensionless % invFroudeSquare = dimensionless_ % invFroudeSquare
         dimensionless % gravity_dir     = dimensionless_ % gravity_dir

      end subroutine SetDimensionless
end module FluidData_iNS
