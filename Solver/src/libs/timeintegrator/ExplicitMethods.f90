!
!////////////////////////////////////////////////////////////////////////
!
!   @File:    ExplicitMethods.f90
!   @Author:  2007-10-23 09:25:32 -0400
!   @Created: David Kopriva
!   @Last revision date: Wed May 5 16:30:01 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: a699bf7e073bc5d10666b5a6a373dc4e8a629897
!
!////////////////////////////////////////////////////////////////////////
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
   use IBMClass
   IMPLICIT NONE

   private
   public   TakeRK3Step, TakeRK5Step, TakeExplicitEulerStep, Enable_CTD_AFTER_STEPS
   public   TakeRKOptStep

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
   SUBROUTINE TakeRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt )
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
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)

      
      INTEGER :: i, j, k, id

      if (present(dt_vec)) then   
         if( mesh% IBM% active ) then
            if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  if( mesh% IBM% TimePenal ) mesh% IBM% penalization(id) = 0.5_RP*dt_vec(id)
                  do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                     if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                        associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                        call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                        end associate
                     end if
                  end do; end do; end do
               end do
!$omp end parallel do
            end if
         end if 
         
         do k = 1,3
            tk = t + b(k)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do
   
         end do ! k

         if( mesh% IBM% active ) then
            if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                     if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                        associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                        call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q )
                        end associate
                     end if
                  end do; end do; end do
               end do
!$omp end parallel do
            end if
         end if 

      else

         if( mesh% IBM% active ) then
            if( mesh% IBM% semiImplicit ) then
               if( mesh% IBM% TimePenal ) mesh% IBM% penalization = 0.5_RP*deltaT
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                     if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                        associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                        call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q )
                        end associate
                     end if
                  end do; end do; end do
               end do
!$omp end parallel do
            end if
         end if 

         do k = 1,3
            tk = t + b(k)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do
       
         end do ! k

         if( mesh% IBM% active ) then
            if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                     if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                        associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                        call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                        end associate
                     end if
                  end do; end do; end do
               end do
!$omp end parallel do
            end if
         end if 
      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeRK3Step

   SUBROUTINE TakeRK5Step( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt )
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
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, i, j, k
      integer, parameter         :: N_STAGES = 5
      real(kind=RP), parameter  :: a(N_STAGES) = [0.0_RP , -0.4178904745_RP, -1.192151694643_RP ,     -1.697784692471_RP , -1.514183444257_RP ]
      real(kind=RP), parameter  :: b(N_STAGES) = [0.0_RP , 0.1496590219993_RP , 0.3704009573644_RP , 0.6222557631345_RP , 0.9582821306748_RP ]
      real(kind=RP), parameter  :: c(N_STAGES) = [0.1496590219993_RP , 0.3792103129999_RP , 0.8229550293869_RP , 0.6994504559488_RP , 0.1530572479681_RP]


      if (present(dt_vec)) then 
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               if( mesh% IBM% TimePenal ) mesh% IBM% penalization(id) = 0.5_RP*dt_vec(id)
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q )
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if 

      DO k = 1, N_STAGES

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO
      
      if( mesh% IBM% active ) then 
         if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if

      else
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
            if( mesh% IBM% TimePenal ) mesh% IBM% penalization = 0.5_RP*deltaT
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if

      DO k = 1, N_STAGES

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

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
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                    end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeRK5Step

   SUBROUTINE TakeExplicitEulerStep( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt )
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
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      !
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, k

      CALL ComputeTimeDerivative( mesh, particles, t, CTD_IGNORE_MODE)
      if ( present(dts) ) then
         if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
      end if

      if (present(dt_vec)) then
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + dt_vec(id)*mesh % elements(id) % storage % QDot
         END DO
!$omp end parallel do
      else
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
         END DO
!$omp end parallel do
      end if

      call checkForNan(mesh, t)
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
   SUBROUTINE TakeRKOptStep( mesh, particles, t, deltaT, ComputeTimeDerivative , N_STAGES, dt_vec, dts, global_dt )
!
!        *****************************************************************************************
!       Optimal RK coefficients from Bassi2009
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
      integer, intent(in)         :: N_STAGES
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, i, j, k
      real(kind=RP), dimension(6,7) :: Am, Bm
      real(kind=RP) :: a(N_STAGES), b(N_STAGES)

      Bm(1,:) = (/ 0.9998596842476678, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Bm(2,:) = (/ 0.528003175664866, 0.5193233361621609, 0.3209132144853066, 0.0, 0.0, 0.0, 0.0 /)
      Bm(3,:) = (/ 0.4062766523561424, 0.3590274668186006, 0.2782786562184366, 0.3031000737788218, 0.0, 0.0, 0.0 /)
      Bm(4,:) = (/ 0.32451232607547, 0.2850381110111294, 0.2299189950459751, 0.3245118697485674, 0.1925289659886354, 0.0, 0.0 /)
      Bm(5,:) = (/ 0.2712469842000987, 0.2506886071464794, 0.1571621623659122, 0.2281484031198761, 0.2523511867383585, 0.1918765315676819, 0.0 /)
      Bm(6,:) = (/ 0.2328811281838825, 0.2007437175473038, 0.1577268986576998, 0.2052094755795549, 0.2138853585222901, 0.192146045128098, 0.1428487872093191 /)

      Am(1,:) = (/ 0.0 / 1.0, -0.4998596842476678 / 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Am(2,:) = (/ 0.0 / 1.0, -0.1645541151603977 / 0.315637725010225, -0.20368561115193584 / 0.3209132144853066, 0.0, 0.0, 0.0, 0.0 /)
      Am(3,:) = (/ 0.0 / 1.0, -0.11076800268308828 / 0.1937832814991212, -0.1652441853194794 / 0.207607995049003, &
            -0.0706706611694336 / 0.3031000737788218, 0.0, 0.0, 0.0 /)
      Am(4,:) = (/ 0.0 / 1.0, -0.11517541944711737 / 0.1896693024156952, -0.0953688085954342 / 0.2159361297115306, &
            -0.01398286533444451 / 0.1925286952557861, -0.1319831744927813 / 0.1925289659886354, 0.0, 0.0 /)
      Am(5,:) = (/ 0.0 / 1.0, -0.056511008554826186 / 0.1229547684000589, -0.1277338387464205 / 0.1094339401033201, &
            -0.047728222262592115 / 0.1824888900957738, -0.045659513024102316 / 0.1785098941878928,-0.0738412925504657 / 0.1918765315676819, 0.0 /)
      Am(6,:) = (/ 0.0 / 1.0, -0.044038551801784204 / 0.1267393084745009, -0.0740044090728029 / 0.1061061571001407, &
            -0.051620741557559094 / 0.1525740388611983, -0.052635436718356604 / 0.1650271578073516,-0.04885820071493849 / 0.1178619741653911, &
             -0.0742840709627069 / 0.1428487872093191 /)

      a = Am(N_STAGES-1,1:N_STAGES)
      b = Bm(N_STAGES-1,1:N_STAGES)

      tk = t + deltaT

      if (present(dt_vec)) then 
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               if( mesh% IBM% TimePenal ) mesh% IBM% penalization(id) = 0.5_RP*dt_vec(id)
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if 

      DO k = 1, N_STAGES

         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +  dt_vec(id) * mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q   +  b(k) * mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH  +  dt_vec(id) * mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c          +  b(k) * mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if 

      else
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
            if( mesh% IBM% TimePenal ) mesh% IBM% penalization = 0.5_RP*deltaT
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q ) 
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if

      DO k = 1, N_STAGES

         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  + deltaT * mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + b(k) * mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + deltaT * mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + b(k) * mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO
      
      if( mesh% IBM% active ) then
         if( mesh% IBM% semiImplicit ) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
               do i = 0, mesh% elements(id)% Nxyz(1); do j = 0, mesh% elements(id)% Nxyz(2); do k = 0, mesh% elements(id)% Nxyz(3)
                  if( mesh% elements(id)% isInsideBody(i,j,k) ) then
                     associate( Q => mesh% elements(id)% storage% Q(:,i,j,k) )
                     call mesh% IBM% GetSemiImplicitStep( id, 0.5_RP*deltaT, Q )
                     end associate
                  end if
               end do; end do; end do
            end do
!$omp end parallel do
         end if
      end if

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeRKOptStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   Subroutine checkForNan(mesh, t)
!
!        **************************************************************************************************************
!        Look if there is a nan in the solution, if at least one if found, stops and creates a hsol for user inspection
!        **************************************************************************************************************
!
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none

      type(HexMesh), intent(in)     :: mesh
      real(kind=RP), intent(in)     :: t

      !local variables
      integer                       :: eID, i, j, k
      CHARACTER(len=LINE_LENGTH)    :: FinalName      !  Final name for particular restart file
      logical                       :: NanNotFound, allNan
      integer                       :: ierr

      ! use not found instead of found as OMP reduction initialized the private value as true
      ! this is redundant for the OMP but left for non parallel compilations
      NanNotFound = .TRUE.

!$omp parallel do reduction(.AND.:NanNotFound) schedule(runtime)
      do eID=1, mesh % no_of_elements
         if ( any(isnan(mesh % elements(eID) % storage % Q))) then
            NanNotFound = .FALSE.
         endif
      end do
!$omp end parallel do

#ifdef _HAS_MPI_
      call mpi_allreduce(NanNotFound, allNan, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      NanNotFound = allNan
#endif

      if (.not. NanNotFound) then
          if ( MPI_Process % isRoot ) print*, "Numerical divergence obtained in solver."
          WRITE(FinalName,'(A,ES11.5,A)')  'RESULTS/horses_divergence_',t,'.hsol'
          if ( MPI_Process % isRoot ) print *, "Writing failed solution: ", FinalName
          call mesh % SaveSolution(0, t, FinalName, .FALSE.)
      end if

#ifdef _HAS_MPI_
     call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

      if (.not. NanNotFound) call exit(99)


   End Subroutine checkForNan
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE ExplicitMethods
