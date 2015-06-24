!
!////////////////////////////////////////////////////////////////////////
!
!      BoundaryConditions.f90
!      Created: June 19, 2015 at 8:43 AM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeBoundaryFlux( e, faceID, t, externalState)  
         USE ElementClass
         USE Physics
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(element) :: e
         INTEGER       :: faceID
         REAL(KIND=RP) :: t
         
         INTEGER       :: i, j
         INTEGER       :: N
         REAL(KIND=RP) :: nHat(3)
         EXTERNAL      :: externalState
         
         N = e % N
         
         DO j = 0, N
            DO i = 0, N
               nHat = e % geom % normal(:,i,j,faceID)
               e % FStarb(:,i,j,faceID) = e % Qb(:,i,j,faceID) * e % geom % scal(i,j,faceID) * (Nhat(1) + nHat(2) + nHat(3))
            END DO   
         END DO  
         
      END SUBROUTINE computeBoundaryFlux
