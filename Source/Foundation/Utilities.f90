!
! /////////////////////////////////////////////////////////////////////
!
!
!     Utilities.F
!
!!
!!     Modification History:
!!        version 0.0 June 1, 2005 David A. Kopriva
!!
!!     The main entry here is AlmostEqual, which returns true if the arguments are 
!!     within rounding error of each other. No adjustments are made for scaling;
!!     We assume numbers are in [-1,1] since this routine is meant to be used
!!     by the Gauss point routines.
!
!      PUBLIC DATA:
!          None
!      PUBLIC METHODS:
!          ALGORITHM 139: AlmostEqual( a, b )
!
!!     @author David A. Kopriva
!
! /////////////////////////////////////////////////////////////////////
!
!-----------------------------------------------------------------------
!! Returns .TRUE. if two numbers are within rounding error of each other
!-----------------------------------------------------------------------
!
      LOGICAL FUNCTION AlmostEqual( a, b ) 
      USE SMConstants
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: a, b
!
      IF ( a == 0.0_RP .OR. b == 0.0_RP )     THEN
         IF ( ABS(a-b) <= 2*EPSILON(b) )     THEN
            AlmostEqual = .TRUE.
         ELSE
            AlmostEqual = .FALSE.
         END IF
      ELSE
         IF( ABS( b - a ) <= 2*EPSILON(b)*MAX(ABS(a), ABS(b)) )     THEN
            AlmostEqual = .TRUE.
         ELSE
            AlmostEqual = .FALSE.
         END IF
      END IF

      END FUNCTION AlmostEqual
