!
! /////////////////////////////////////////////////////////////////////
!
!     The main entry here is AlmostEqual, which returns true if the arguments are 
!     within rounding error of each other. No adjustments are made for scaling;
!     We assume numbers are in [-1,1] since this routine is meant to be used
!     by the Gauss point routines.
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
   public   AlmostEqual, UnusedUnit, SolveThreeEquationLinearSystem, GreatestCommonDivisor, outer_product, AlmostEqualRelax
   public   toLower, Qsort, QsortWithFriend, BubblesortWithFriend, my_findloc
   public   logarithmicMean, dot_product
   public   LeastSquaresLinRegression
   
   interface dot_product
      module procedure dot_product_3Tensor_Vec
   end interface
   contains
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------
!     Function to compute the GreatestCommonDivisor. Modified from:
!     -> https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap04/gcd.html
!     -------------------------------------------------------------------
      pure function GreatestCommonDivisor(a,b) result(c)
         implicit none
         !-arguments-------------------------------------
         integer, intent(in) :: a, b
         integer             :: c
         !-local-variables-------------------------------
         integer          :: a1, b1, c1
         !-----------------------------------------------
         
         if (a < b) then
            c1 = a
            a1 = b
            b1 = a
         else
            a1 = a
            b1 = b
         end if

         do                    ! now we have a <= b
            c1 = MOD(a1, b1)   !    compute c, the reminder
            IF (c1 == 0) exit  !    if c is zero, we are done.  GCD = b
            a1 = b1            !    otherwise, b becomes a
            b1 = c1            !    and c becomes b
         END DO                !    go back
         
         c = b1
      end function GreatestCommonDivisor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure function outer_product(A,B) result(C)
         implicit none
         !-arguments----------------------------------------------
         real(kind=RP), intent(in) :: A(:)
         real(kind=RP), intent(in) :: B(:)
         real(kind=RP)             :: C(size(A),size(B))
         !-local-variables----------------------------------------
         real(kind=RP) :: Amat(size(A),size(B))
         real(kind=RP) :: Bmat(size(A),size(B))
         !--------------------------------------------------------
         
         Amat = spread( A, dim = 2, ncopies = size(B) )
         Bmat = spread( B, dim = 1, ncopies = size(A) )
         
         C = Amat * Bmat
         
      end function outer_product
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!     dot_product_3Tensor_Vec:
!        Inner product of a third-order tensor with a vector
!     ------------------------------------------------------
      pure function dot_product_3Tensor_Vec(X,Y) result(Z)
         implicit none
         !-arguments----------------------------------------------
         real(kind=RP), intent(in) :: X(:,:,:)
         real(kind=RP), intent(in) :: Y(size(X,3))
         real(kind=RP)             :: Z(size(X,1),size(X,2))
         !-local-variables----------------------------------------
         integer :: i
         !--------------------------------------------------------
         
         Z = 0._RP
         do i=1, size(Y)
            Z = Z + X(:,:,i) * Y(i)
         end do
      end function dot_product_3Tensor_Vec
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      LOGICAL elemental FUNCTION AlmostEqual( a, b, tol )
      USE SMConstants
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), intent(in) :: a, b
      REAL(KIND=RP), intent(in), optional :: tol
      REAL(KIND=RP) :: mytol
!
      mytol = 2*epsilon(b)
      if (present(tol)) mytol = tol

      IF ( a == 0.0_RP .OR. b == 0.0_RP )     THEN
         IF ( ABS(a-b) <= mytol )     THEN
            AlmostEqual = .TRUE.
         ELSE
            AlmostEqual = .FALSE.
         END IF
      ELSE
         IF( ABS( b - a ) <= mytol*MAX(ABS(a), ABS(b)) )     THEN
            AlmostEqual = .TRUE.
         ELSE
            AlmostEqual = .FALSE.
         END IF
      END IF

      END FUNCTION AlmostEqual
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      LOGICAL pure FUNCTION AlmostEqualRelax( a, b, reference_length, tol) 
      USE SMConstants
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), intent(in) :: a, b, reference_length
      REAL(KIND=RP), intent(in), optional :: tol
!
      REAL(KIND=RP)             :: mytol
!
      mytol = 1.0E-4_RP
      if (present(tol)) mytol = tol

      IF( ABS( b - a ) <= mytol * reference_length) THEN
          AlmostEqualRelax = .TRUE.
      ELSE
          AlmostEqualRelax = .FALSE.
      END IF

      END FUNCTION AlmostEqualRelax
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
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for performing least square regression and giving up the coefficients
!  -----------------------------------------------------------------------
   pure subroutine LeastSquaresLinRegression(N,x,y,C,eta,r)
      implicit none
      !--------------------------------------
      integer      , intent(in)  :: N        !< Number of points
      real(kind=RP), intent(in)  :: x(N) 
      real(kind=RP), intent(in)  :: y(N)
      real(kind=RP), intent(out) :: C
      real(kind=RP), intent(out) :: eta
      real(kind=RP), intent(out) :: r
      !--------------------------------------
      real(kind=RP)              :: sumx,sumy,sumsqx,sumxy,deno,sumsqy
      !--------------------------------------
      
      sumx = sum(x)
      sumy = sum(y)
      sumsqx = sum(x*x)
      sumsqy = sum(y*y)
      sumxy  = sum(x*y)
      
      deno=n*sumsqx-sumx*sumx
      
      eta = (n*sumxy-sumx*sumy)/deno
      C   = (sumsqx*sumy-sumx*sumxy)/deno
      r   = (  sumxy-sumx*sumy/n) / sqrt((sumsqx - sumx*sumx/n) * (sumsqy - sumy*sumy/n))
      
   end subroutine LeastSquaresLinRegression
!
!////////////////////////////////////////////////////////////////////////
!
!      Qsort.f90
!      Created: 2016-08-25 14:30 (GMT+0)
!      By: Juli Rew (juliana@ucar.edu)
!     Recursive Fortran 95 quicksort rOUTINe
!        sorts INTEGER numbers INto ascENDing numerical order
!        Based on algorithm from Cormen et al., Introduction to Algorithms,
!        1997 prINting
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
!
!////////////////////////////////////////////////////////////////////////
!
!     A routine to perform quicksort for array A and reorder array B with the same operations 
!
!////////////////////////////////////////////////////////////////////////////////////////
!
RECURSIVE SUBROUTINE QsortWithFriend(A,B)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(IN OUT), DIMENSION(size(A)) :: B
  INTEGER                               :: iq
  
  IF(size(A) > 1) THEN
     call PartitionWithFriend(A, B, iq)
     call QsortWithFriend(A(:iq-1),B(:iq-1))
     call QsortWithFriend(A(iq:),B(iq:))
  ENDIF
END SUBROUTINE QsortWithFriend

SUBROUTINE PartitionWithFriend(A, B, marker)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: B
  INTEGER, INTENT(OUT)                  :: marker
  INTEGER                               :: i, j
  INTEGER                               :: temp
  INTEGER                               :: tempB
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
        
        tempB = B(i)
        B(i) = B(j)
        B(j) = tempB
     elseIF (i == j) THEN
        marker = i+1
        return
     else
        marker = i
        return
     ENDIF
  END DO

END SUBROUTINE PartitionWithFriend
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
pure subroutine BubblesortWithFriend(A, B)
    implicit none
    !--------------------------------------
    real(RP), intent(inout) :: A(:)
    integer,  intent(inout) :: B(size(A))
    !--------------------------------------
    integer :: n, newn
    integer :: i, temp
    !--------------------------------------

    n = size(A)
    do
        newn = 0
        do i = 1, n-1
            if (A(i) > A(i+1)) then
                ! Swap A
                temp = A(i)
                A(i) = A(i+1)
                A(i+1) = temp
                ! Swap B
                temp = B(i)
                B(i) = B(i+1)
                B(i+1) = temp
                ! Update number of unsorted elements
                newn = i
            end if
        end do
        n = newn
        if (n <= 1) exit
    end do

end subroutine BubblesortWithFriend
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
integer function my_findloc(arr, val, dim)
!  ---------------------------------------------------------
!  Vanilla routine to find an index of matching value in the array.
!  For INTEL or GNU v>9.0 just switch to 'findloc'.
!  ---------------------------------------------------------
      implicit none
      !-----Arguments---------------------------------------------------
      integer,dimension(:),intent(in) :: arr
      integer,             intent(in) :: val, dim
      !-----Local-Variables---------------------------------------------
      integer :: i
      !  -----------------------------------------------------------------------

      ! my_findloc = findloc(arr,val,dim) ! intrinsic function
      my_findloc = 0
      do i = 1, size(arr,1)
         if (arr(i) .eq. val) then
         my_findloc = i
         exit
         end if
      end do
       
   end function my_findloc
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module Utilities
