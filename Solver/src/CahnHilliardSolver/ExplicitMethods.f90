!
!//////////////////////////////////////////////////////
!
!   @File:    ExplicitMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:37 2018
!   @Last revision date: Wed Jan 31 18:27:04 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 1181c365aba00e78739d327d06901d6d8ca99e02
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      EcplicitMethods.f90
!      Created: 2007-10-23 09:25:32 -0400 
!      By: David Kopriva  
!
!      RK integrators for DG approximation to conservation
!      laws in 3D
!
!////////////////////////////////////////////////////////////////////////
!
MODULE ExplicitMethods
   USE SMConstants
   USE DGSEMClass

   IMPLICIT NONE

   private
   public TakeRK3Step, TakeRK5Step, TakeExplicitEulerStep
!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a RK3 step.
!  ------------------------------
   SUBROUTINE TakeRK3Step( sem, t, deltaT )
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
      TYPE(DGSem)     :: sem
      REAL(KIND=RP)   :: t, deltaT, tk
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
         CALL ComputeTimeDerivative( sem, tk )
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % storage % G = a(k)*sem % mesh % elements(id) % storage % G  +             sem % mesh % elements(id) % storage % QDot
            sem % mesh % elements(id) % storage % Q =      sem % mesh % elements(id) % storage % Q  + c(k)*deltaT*sem % mesh % elements(id) % storage % G
         END DO
!$omp end parallel do
         
      END DO
      
   END SUBROUTINE TakeRK3Step

   subroutine TakeRK5Step(sem, t, deltaT)
!  
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
      implicit none
      TYPE(DGSem)     :: sem
      REAL(KIND=RP)   :: t, deltaT, tk
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
         CALL ComputeTimeDerivative( sem, tk )
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % storage % G = a(k)*sem % mesh % elements(id) % storage % G  +             sem % mesh % elements(id) % storage % QDot
            sem % mesh % elements(id) % storage % Q =      sem % mesh % elements(id) % storage % Q  + c(k)*deltaT*sem % mesh % elements(id) % storage % G
         END DO
!$omp end parallel do
         
      END DO

   end subroutine TakeRK5Step

   subroutine TakeExplicitEulerStep(sem, t, deltaT)
!  
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
      implicit none
      TYPE(DGSem)     :: sem
      REAL(KIND=RP)   :: t, deltaT, tk
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, k

      CALL ComputeTimeDerivative( sem, tk )
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % storage % Q = sem % mesh % elements(id) % storage % Q  + deltaT*sem % mesh % elements(id) % storage % QDot
         END DO
!$omp end parallel do
         
   end subroutine TakeExplicitEulerStep


!
!///////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE ExplicitMethods
