!
!////////////////////////////////////////////////////////////////////////
!
!      ExplicitMethods.f90
!      Created: 2007-10-23 09:25:32 -0400 
!      By: David Kopriva  
!
!      RK integrators for DG approximation to conservation
!      laws in 3D
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE ExplicitMethods
   USE SMConstants
   use HexMeshClass
   use TimeIntegratorDefinitions
   use DGSEMClass, only: ComputeTimeDerivative_f
   use ParticlesClass
   use PhysicsStorage, only: CTD_IGNORE_MODE
   IMPLICIT NONE

   private
   public   TakeRK3Step, TakeRK5Step, TakeExplicitEulerStep, Enable_CTD_AFTER_STEPS

   integer, protected :: eBDF_order = 3
   logical, protected :: CTD_AFTER_STEPS = .false.
!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a RK3 step.
!  ------------------------------
   SUBROUTINE TakeRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative )
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      type(HexMesh)      :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)
      
      INTEGER :: k, id
      
      DO k = 1,3
         
         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do
         
      END DO
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

!$omp parallel do schedule(runtime)
      do k=1, mesh % no_of_elements
         if ( any(isnan(mesh % elements(k) % storage % Q))) then
            print*, "Numerical divergence obtained in solver."
            call exit(99)
         endif
      end do
!$omp end parallel do
      
   END SUBROUTINE TakeRK3Step

   SUBROUTINE TakeRK5Step( mesh, particles, t, deltaT, ComputeTimeDerivative )
!  
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, k
      integer, parameter         :: N_STAGES = 5
      real(kind=RP), parameter  :: a(N_STAGES) = [0.0_RP , -0.4178904745_RP, -1.192151694643_RP ,     -1.697784692471_RP , -1.514183444257_RP ]
      real(kind=RP), parameter  :: b(N_STAGES) = [0.0_RP , 0.1496590219993_RP , 0.3704009573644_RP , 0.6222557631345_RP , 0.9582821306748_RP ]
      real(kind=RP), parameter  :: c(N_STAGES) = [0.1496590219993_RP , 0.3792103129999_RP , 0.8229550293869_RP , 0.6994504559488_RP , 0.1530572479681_RP]

      DO k = 1, N_STAGES
         
         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do
         
      END DO
      
!$omp parallel do schedule(runtime)
      do k=1, mesh % no_of_elements
         if ( any(isnan(mesh % elements(k) % storage % Q))) then
            print*, "Numerical divergence obtained in solver."
            call exit(99)
         endif
      end do
!$omp end parallel do

!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
      
   end subroutine TakeRK5Step

   SUBROUTINE TakeExplicitEulerStep( mesh, particles, t, deltaT, ComputeTimeDerivative )
!  
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      !
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, k

      CALL ComputeTimeDerivative( mesh, particles, t, CTD_IGNORE_MODE)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
         END DO
!$omp end parallel do

!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         
   end subroutine TakeExplicitEulerStep

   subroutine TakeExplicitBDFStep(mesh, particles, t, deltaT, ComputeTimeDerivative)
      implicit none
      type(HexMesh)                      :: mesh
#ifdef FLOW
      type(Particles_t)                  :: particles
#else
      logical                            :: particles
#endif
      REAL(KIND=RP)                      :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                  :: id
      real(kind=RP), parameter :: invGamma1 = 1.0_RP, invGamma2 = 2.0_RP/3.0_RP, invGamma3 = 6.0_RP / 11.0_RP
      logical, save            :: isFirst = .true., isSecond = .false., isThird = .false.

      if (isThird) then      
!
!        Perform the third order stages
!        ------------------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
!           Set y^{*,n+1} in Q and downgrade y^n and y^{n-1}
!           ------------------------------------------------
            mesh % elements(id) % storage % QDot = mesh % elements(id) % storage % prevQ(2) % Q
            mesh % elements(id) % storage % prevQ(2) % Q = mesh % elements(id) % storage % prevQ(1) % Q
            mesh % elements(id) % storage % prevQ(1) % Q = mesh % elements(id) % storage % Q
            mesh % elements(id) % storage % Q =   3.0_RP * mesh % elements(id) % storage % prevQ(1) % Q &
                                                  - 3.0_RP * mesh % elements(id) % storage % prevQ(2) % Q & 
                                                  + mesh % elements(id) % storage % QDot
         end do
!$omp end parallel do
!
!        Compute QDot
!        ------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
!
!        Perform the time-step
!        ---------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q =   2.0_RP * mesh % elements(id) % storage % prevQ(1) % Q &
                                                - 0.5_RP * mesh % elements(id) % storage % prevQ(2) % Q &
                                                + (1.0_RP/3.0_RP) * mesh % elements(id) % storage % Q & 
                                                + deltaT * mesh % elements(id) % storage % QDot  
            mesh % elements(id) % storage % Q = invGamma3 * mesh % elements(id) % storage % Q
         end do
!$omp end parallel do


      elseif (isSecond) then
!
!        Perform the second order stages
!        -------------------------------
         if (eBDF_ORDER > 2) then
!$omp parallel do schedule(runtime)
            do id = 1, size(mesh % elements)          
!
!              Set for the previous solution
!              -----------------------------
               mesh % elements(id) % storage % prevQ(2) % Q = mesh % elements(id) % storage % prevQ(1) % Q
            end do
!$omp end parallel do
         end if

!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)          
!
!           Set y^{*,n+1} in Q
!           ------------------
            mesh % elements(id) % storage % Q = 2.0_RP*mesh % elements(id) % storage % Q -mesh % elements(id) % storage % prevQ(1) % Q
!
!           Set y^{n} in prevQ
!           --------------------
            mesh % elements(id) % storage % prevQ(1) % Q = 0.5_RP*(mesh % elements(id) % storage % Q + mesh % elements(id) % storage % prevQ(1) % Q)

         end do
!$omp end parallel do
!
!        Compute QDot
!        ------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
!
!        Perform the time-step
!        ---------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)          
            mesh % elements(id) % storage % Q =   mesh % elements(id) % storage % prevQ(1) % Q + 0.5_RP*mesh % elements(id) % storage % Q  &
                                                + deltaT * mesh % elements(id) % storage % QDot
            mesh % elements(id) % storage % Q = invGamma2 * mesh % elements(id) % storage % Q

         end do
!$omp end parallel do

         if (eBDF_ORDER > 1 ) then
!
!           Move to third order
!           -------------------
            isFirst = .false.
            isSecond = .false.
            isThird = .true.
         end if

      elseif (isFirst) then
!
!        Perform the first order stages
!        ------------------------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            if (eBDF_ORDER > 1 ) then
!
!              Set for the previous solution
!              -----------------------------
               mesh % elements(id) % storage % prevQ(1) % Q = mesh % elements(id) % storage % Q

            end if

            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot

         END DO
!$omp end parallel do

            if (eBDF_ORDER > 1 ) then
!
!              Move to second order
!              --------------------
               isFirst = .false.
               isSecond = .true.
            end if

      end if




   end subroutine TakeExplicitBDFStep

   subroutine Enable_CTD_AFTER_STEPS()
      implicit none
      CTD_AFTER_STEPS = .true.
   end subroutine Enable_CTD_AFTER_STEPS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE ExplicitMethods
