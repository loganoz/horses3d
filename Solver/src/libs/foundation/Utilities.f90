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
!
! /////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!>    This subroutine returns an unused
!!    unit number. It is assumed that the unit number will be used
!!    immediately after it is assigned. It is also assumed that the
!!    compiler has access to units 1-99. The function returns "0"
!!    if it cannot find an unused unit number.
!     ----------------------------------------------------------------
!
      INTEGER FUNCTION UnusedUnit()
!
         IMPLICIT NONE
         INTEGER :: j
         LOGICAL :: unitIsOpened, unitDoesExist
         
         UnusedUnit = 0
         DO j = 1, 99
           INQUIRE( UNIT = j, OPENED = unitIsOpened, EXIST = unitDoesExist )
           IF ( .NOT.unitIsOpened ) EXIT
         END DO
         
         IF (  j <= 99 .AND. unitDoesExist )     THEN
            UnusedUnit = j
         ELSE
            UnusedUnit = 0
         END IF
!
      END FUNCTION UnusedUnit
!
! /////////////////////////////////////////////////////////////////////
!
      function SolveThreeEquationLinearSystem(A,b)
         use SMConstants
         implicit none
         real(kind=RP), intent(in)  :: A(3,3)
         real(kind=RP), intent(in)  :: b(3)
         real(kind=RP)     :: SolveThreeEquationLinearSystem(3)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: detA
         real(kind=RP)  :: r1, r2, r3

         detA =   A(1,1) * A(2,2) * A(3,3) &
                + A(2,1) * A(3,2) * A(1,3) &
                + A(1,2) * A(2,3) * A(3,1) &
                - A(2,2) * A(1,3) * A(3,1) &
                - A(1,1) * A(2,3) * A(3,2) &
                - A(3,3) * A(1,2) * A(2,1)

         r1   =   b(1) * A(2,2) * A(3,3) &
                + b(2) * A(3,2) * A(1,3) &
                + A(1,2) * A(2,3) * b(3) &
                - A(2,2) * A(1,3) * b(3) &
                - b(1) * A(2,3) * A(3,2) &
                - A(3,3) * A(1,2) * b(2)

         r2   =   A(1,1) * b(2) * A(3,3) &
                + A(2,1) * b(3) * A(1,3) &
                + b(1) * A(2,3) * A(3,1) &
                - b(2) * A(1,3) * A(3,1) &
                - A(1,1) * A(2,3) * b(3) &
                - A(3,3) * b(1) * A(2,1)

         r3   =   A(1,1) * A(2,2) * b(3) &
                + A(2,1) * A(3,2) * b(1) &
                + A(1,2) * b(2) * A(3,1) &
                - A(2,2) * b(1) * A(3,1) &
                - A(1,1) * b(2) * A(3,2) &
                - b(3) * A(1,2) * A(2,1)

         detA = 1.0_RP / detA

         SolveThreeEquationLinearSystem(1) = r1 * detA
         SolveThreeEquationLinearSystem(2) = r2 * detA
         SolveThreeEquationLinearSystem(3) = r3 * detA


      end function SolveThreeEquationLinearSystem
!
! /////////////////////////////////////////////////////////////////////
!
   subroutine toLower(str)
!
!  ----------------
!  From ResettaCode
!  ----------------
!
     character(*), intent(in out) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do  
   end subroutine toLower
