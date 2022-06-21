!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Mon May 14 19:03:28 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
module PhysicsStorage
   use SMConstants, only: RP
#if defined(NAVIERSTOKES)
   use PhysicsStorage_NS
#endif
#if defined(CAHNHILLIARD)
   use PhysicsStorage_CH
#endif
   implicit none

#if (defined(NAVIERSTOKES) && !defined(CAHNHILLIARD))
   integer, parameter   :: NCONS = NCONS
   integer, parameter   :: NTOTALGRADS = NGRAD
#elif (!defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
   integer, parameter   :: NCONS = NCOMP
   integer, parameter   :: NTOTALGRADS = NCOMP
#elif (defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
   integer, parameter   :: NCONS = NCONS + NCOMP
   integer, parameter   :: NTOTALGRADS = NGRAD + NCOMP
#endif

   real(kind=RP)     :: Lref
   real(kind=RP)     :: timeref
   
   contains
      SUBROUTINE ConstructPhysicsStorage( machArg, REArg, PRArg, flowIsNavierStokesArg )
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: machArg, REArg, PRArg
      LOGICAL       :: flowIsNavierStokesArg
!
!        Default values
!        --------------
         Lref = 1.0_RP
         timeref = 1.0_RP
!
!        Construct NSE physics
!        ---------------------
         call ConstructPhysicsStorage_NS( machArg, REArg, PRArg, flowIsNavierStokesArg )
   
      end subroutine ConstructPhysicsStorage
end module PhysicsStorage
