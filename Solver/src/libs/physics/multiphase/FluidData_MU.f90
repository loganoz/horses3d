!
!//////////////////////////////////////////////////////
!
!   @File:    FluidData_MU.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:25 2018
!   @Last revision date: Mon Jul  2 14:17:28 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7af1f42fb2bc9ea3a0103412145f2a925b4fac5e
!
!//////////////////////////////////////////////////////
!
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
      integer                    :: number_of_fluids
      real(kind=RP)              :: rho0c02   ! Artificial compressibility constant
      real(kind=RP), allocatable :: rho(:)
      real(kind=RP), allocatable :: mu(:)
   end type Thermodynamics_t

   type RefValues_t
      real(kind=RP)        :: V
      real(kind=RP)        :: p
      real(kind=RP)        :: rho
      real(kind=RP)        :: mu
      real(kind=RP)        :: g0
      real(kind=RP)        :: AoATheta
      real(kind=RP)        :: AoAPhi
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP), allocatable :: rho(:)
      real(kind=RP), allocatable :: mu(:)
      real(kind=RP)              :: Re
      real(kind=RP)              :: Fr
      real(kind=RP)              :: invFr2
      real(kind=RP)              :: gravity_dir(NDIM)
      real(kind=RP)              :: rho_max
      real(kind=RP)              :: rho_min
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

         thermodynamics % number_of_fluids = thermodynamics_ % number_of_fluids  
         thermodynamics % rho0c02          = thermodynamics_ % rho0c02              

         if(allocated(thermodynamics % rho)) deallocate(thermodynamics % rho)
         allocate(thermodynamics % rho(thermodynamics % number_of_fluids))
         thermodynamics % rho              = thermodynamics_ % rho              

         if(allocated(thermodynamics % mu )) deallocate(thermodynamics % mu )
         allocate(thermodynamics % mu (thermodynamics % number_of_fluids))
         thermodynamics % mu               = thermodynamics_ % mu               


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

         if(allocated(dimensionless % rho)) deallocate(dimensionless % rho)
         allocate(dimensionless % rho(size(dimensionless_ % rho)))
         dimensionless % rho              = dimensionless_ % rho              

         if(allocated(dimensionless % mu )) deallocate(dimensionless % mu )
         allocate(dimensionless % mu (size(dimensionless_ % mu)))
         dimensionless % mu               = dimensionless_ % mu               

         dimensionless % Re              = dimensionless_ % Re
         dimensionless % Fr              = dimensionless_ % Fr
         dimensionless % invFr2 = dimensionless_ % invFr2
         dimensionless % gravity_dir     = dimensionless_ % gravity_dir

         dimensionless % rho_max          = dimensionless_ % rho_max          
         dimensionless % rho_min          = dimensionless_ % rho_min          

      end subroutine SetDimensionless
end module FluidData_MU
