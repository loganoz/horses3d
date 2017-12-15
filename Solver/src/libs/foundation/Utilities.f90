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
#include "Includes.h"
module Utilities
   use SMConstants
   implicit none

   private
   public   AlmostEqual, UnusedUnit, SolveThreeEquationLinearSystem
   public   toLower, Qsort
   public   logarithmicMean


   contains
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

   pure subroutine logarithmicMean(aL, aR, lnMean)
!
!     ************************************************
!        Algorithm by Ismail and Roe
!
!     ************************************************
!
      implicit none
      real(kind=RP), intent(in)  :: aL, aR
      real(kind=RP), intent(out) :: lnMean
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP), parameter   :: eps = 0.01_RP
      real(kind=RP)              :: xi, f, u, FF

      xi = aL / aR
      f = (xi - 1.0_RP) / (xi + 1.0_RP)
      u = f * f

      if ( u .lt. eps ) then
         FF = 1.0_RP + (1.0_RP / 3.0_RP) * u + (1.0_RP / 5.0_RP) * POW2(u) + (1.0_RP / 7.0_RP) * POW3(u)
   
      else
         FF = log(xi) / ( 2.0_RP * f )

      end if

      lnMean = 0.5_RP * (aL + aR) / FF
         
   end subroutine logarithmicMean
!
!////////////////////////////////////////////////////////////////////////
!
!      Qsort.f90
!      Created: 2016-08-25 14:30 (GMT+0)
!      By: Juli Rew (juliana@ucar.edu)
! 		Recursive Fortran 95 quicksort rOUTINe
! 			sorts INTEGER numbers INto ascENDing numerical order
!  		Based on algorithm from Cormen et al., Introduction to Algorithms,
! 			1997 prINting
!
!////////////////////////////////////////////////////////////////////////////////////////
!
RECURSIVE SUBROUTINE Qsort(A)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER                               :: iq

  IF(size(A) > 1) THEN
     call Partition(A, iq)
     call Qsort(A(:iq-1))
     call Qsort(A(iq:))
  ENDIF
END SUBROUTINE Qsort

SUBROUTINE Partition(A, marker)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(OUT)                  :: marker
  INTEGER                               :: i, j
  INTEGER                               :: temp
  INTEGER                               :: x      ! pivot poINt
  x = A(1)
  i = 0
  j = size(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseIF (i == j) THEN
        marker = i+1
        return
     else
        marker = i
        return
     ENDIF
  END DO

END SUBROUTINE Partition

end module Utilities

