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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE ExplicitMethods
