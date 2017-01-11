!
!////////////////////////////////////////////////////////////////////////
!
!      InitialFlowState.f90
!      Created: June 11, 2015 at 5:43 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initialFlowState(x,t,state)
         USE PhysicsStorage
         IMPLICIT NONE  
         REAL(KIND=RP) :: x(3)
         REAL(KIND=RP) :: t
         REAL(KIND=RP) :: state(N_EQN)
         
         state(1) = x(1)
         state(2) = x(2)
         state(3) = x(3)
         state(4) = 1.0_RP
         state(5) = SUM(x)
         
!         state(1) = x(1)
!         state(2) = x(2)
!         state(3) = x(3)
!         state(4) = SQRT(x(1)**2 + x(2)**2 + x(3)**2)
!         state(5) = SUM(x)

      END SUBROUTINE initialFlowState
!
!////////////////////////////////////////////////////////////////////////
!
         SUBROUTINE externalBoundaryState(x,t,nHat,Q,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: Q(:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalBoundaryState
!
!////////////////////////////////////////////////////////////////////////
!
         SUBROUTINE externalGradientState(x,t,nHat,gradU,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: gradU(:,:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalGradientState
