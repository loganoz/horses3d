!
!////////////////////////////////////////////////////////////////////////
!
!      Assert.f90
!      Created: February 21, 2013 10:02 AM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
!> Defines procedures that test equality of different kinds of arguments.
!> Procedures defined here are USEd by the FTAssertions Module.
!
      Module ComparisonsModule
      USE ISO_FORTRAN_ENV
      USE FTOLConstants
      IMPLICIT NONE
      PRIVATE
      
      INTEGER, PARAMETER, PUBLIC :: ASSERT_SUCCESS       = 0, ASSERT_SIZE_DIFFERS = 1
      INTEGER, PARAMETER, PUBLIC :: ASSERT_VALUES_DIFFER = 2
      
      CHARACTER(LEN=21), PARAMETER :: ASSERT_SIZE_DIFFERS_NAME   = "Array sizes differ"
      CHARACTER(LEN=21), PARAMETER :: ASSERT_VALUES_DIFFERS_NAME = "Array elements differ"
      CHARACTER(LEN=21), PARAMETER :: ASSERT_VALUES_OK_NAME      = "Arrays match"
      
      CHARACTER(LEN=21), PARAMETER, PUBLIC :: compareCodeStrings(0:2) = [ASSERT_VALUES_OK_NAME,   &
                                                                         ASSERT_SIZE_DIFFERS_NAME,&
                                                                         ASSERT_VALUES_DIFFERS_NAME]
      
      INTERFACE isEqual
         MODULE PROCEDURE isEqualTwoIntegers
         MODULE PROCEDURE isEqualTwoIntegerArrays1D
         MODULE PROCEDURE isEqualTwoIntegerArrays2D
         MODULE PROCEDURE isWithinToleranceTwoReal
         MODULE PROCEDURE isWithinToleranceTwoRealArrays1D
         MODULE PROCEDURE isWithinToleranceTwoRealArrays2D
         MODULE PROCEDURE isWithinToleranceTwoDouble
         MODULE PROCEDURE isWithinToleranceTwoDoubleArrays1D
         MODULE PROCEDURE isWithinToleranceTwoDoubleArrays2D
         MODULE PROCEDURE isEqualString
#ifdef _has_Quad
         MODULE PROCEDURE isWithinToleranceTwoQuad
#endif
      END INTERFACE isEqual
      
      TYPE assertInfoArray1D
         CHARACTER(LEN=128)                 :: failureName
         INTEGER                            :: failureType
         LOGICAL, DIMENSION(:), ALLOCATABLE :: locations
      END TYPE assertInfoArray1D
      
      TYPE assertInfoArray2D
         CHARACTER(LEN=128)                   :: failureName
         INTEGER                              :: failureType
         LOGICAL, DIMENSION(:,:), ALLOCATABLE :: locations
      END TYPE assertInfoArray2D
      
      PUBLIC :: isEqual,assertInfoArray1D,assertInfoArray2D
!
!     ========
      CONTAINS
!     ========
!
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isTrue(condition)  
         IMPLICIT NONE  
         LOGICAL :: condition
         IF ( condition )     THEN
            isTrue = .true. 
         ELSE 
            isTrue = .false. 
         END IF 
      END FUNCTION isTrue
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isFalse(condition)  
         IMPLICIT NONE  
         LOGICAL :: condition
         IF ( .NOT.condition )     THEN
            isFalse = .true. 
         ELSE 
            isFalse = .false. 
         END IF 
      END FUNCTION isFalse
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isEqualTwoIntegers(i,j)  
         IMPLICIT NONE  
         INTEGER, INTENT(in)  :: i, j
         
         IF ( i == j )     THEN
            isEqualTwoIntegers = .true.
         ELSE
            isEqualTwoIntegers = .false.
         END IF 
         
      END FUNCTION isEqualTwoIntegers    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isEqualTwoIntegerArrays1D(a,b,info)  
         IMPLICIT NONE  
         INTEGER, INTENT(in)    , DIMENSION(:)            :: a, b
         TYPE(assertInfoArray1D), INTENT(INOUT), OPTIONAL :: info
         
         isEqualTwoIntegerArrays1D = .true.
         
         IF(PRESENT(info)) THEN
            info % failureType = ASSERT_SUCCESS
            info % failureName = ASSERT_VALUES_OK_NAME
         END IF
        
         IF ( SIZE(a) /= SIZE(b) )     THEN
            isEqualTwoIntegerArrays1D = .false.
            IF(PRESENT(info))     THEN 
               info % failureType = ASSERT_SIZE_DIFFERS
               info % failureName = ASSERT_SIZE_DIFFERS_NAME
            END IF
         ELSE IF(ANY(a /= b))     THEN
            isEqualTwoIntegerArrays1D = .false.
             IF(PRESENT(info))     THEN 
               info % failureType = ASSERT_VALUES_DIFFER
               info % failureName = ASSERT_VALUES_DIFFERS_NAME
               ALLOCATE(info % locations(SIZE(a)))
               info % locations = .true.
               WHERE(a /= b) info % locations = .false.
            END IF
         END IF 
         
      END FUNCTION isEqualTwoIntegerArrays1D
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isEqualTwoIntegerArrays2D(a,b,info)  
         IMPLICIT NONE  
         INTEGER, INTENT(in)    , DIMENSION(:,:)          :: a, b
         TYPE(assertInfoArray2D), INTENT(INOUT), OPTIONAL :: info
         
         isEqualTwoIntegerArrays2D = .true.
         
         IF(PRESENT(info)) THEN
            info % failureType = ASSERT_SUCCESS
            info % failureName = ASSERT_VALUES_OK_NAME
         END IF
         
         IF ( SIZE(a) /= SIZE(b) )     THEN
            isEqualTwoIntegerArrays2D = .false.
            IF(PRESENT(info))     THEN 
               info % failureType = ASSERT_SIZE_DIFFERS
               info % failureName = ASSERT_SIZE_DIFFERS_NAME
            END IF
         ELSE IF(ANY(a /= b))     THEN
            isEqualTwoIntegerArrays2D = .false.
             IF(PRESENT(info))     THEN 
               info % failureType = ASSERT_VALUES_DIFFER
               info % failureName = ASSERT_VALUES_DIFFERS_NAME
               ALLOCATE(info % locations(SIZE(a,1),SIZE(a,2)))
               info % locations = .true.
               WHERE(a /= b) info % locations = .false.
            END IF
         END IF 
         
      END FUNCTION isEqualTwoIntegerArrays2D
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoReal(x,y,tol)  
         IMPLICIT NONE  
         REAL, INTENT(in)  :: x,y,tol
         LOGICAL           :: test
         
         IF ( x <= tol )     THEN
            test = ABS(x-y) <= tol
         ELSE
            test = ABS(x-y) <= tol*MAX(ABS(x),ABS(y))
         END IF 
         
         IF ( test )     THEN
            isWithinToleranceTwoReal = .true.
         ELSE
            isWithinToleranceTwoReal = .false.
         END IF 
         
      END FUNCTION isWithinToleranceTwoReal    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoRealArrays1D(a,b,tol,code)  
         IMPLICIT NONE  
         REAL, INTENT(IN), DIMENSION(:) :: a, b
         REAL, INTENT(IN)               :: tol
         INTEGER, INTENT(OUT), OPTIONAL :: code
         
         isWithinToleranceTwoRealArrays1D = .true.
         IF(PRESENT(code)) code = ASSERT_SUCCESS
         
         IF ( SIZE(a) /= SIZE(b) )     THEN
            isWithinToleranceTwoRealArrays1D = .false.
            IF(PRESENT(code)) code = ASSERT_SIZE_DIFFERS
         ELSE IF(ANY(ABS(a-b) > tol*MAX(ABS(a),ABS(b))))     THEN
            isWithinToleranceTwoRealArrays1D = .false.
            IF(PRESENT(code)) code = ASSERT_VALUES_DIFFER
         END IF 
         
      END FUNCTION isWithinToleranceTwoRealArrays1D
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoRealArrays2D(a,b,tol,code)  
         IMPLICIT NONE  
         REAL, INTENT(IN), DIMENSION(:,:) :: a, b
         REAL, INTENT(IN)                 :: tol
         INTEGER, INTENT(OUT), OPTIONAL   :: code
         
         isWithinToleranceTwoRealArrays2D = .true.
         IF(PRESENT(code)) code = ASSERT_SUCCESS
         
         IF ( SIZE(a) /= SIZE(b) )     THEN
            isWithinToleranceTwoRealArrays2D = .false.
            IF(PRESENT(code)) code = ASSERT_SIZE_DIFFERS
         ELSE IF(ANY(ABS(a-b)> tol*MAX(ABS(a),ABS(b))))     THEN
            isWithinToleranceTwoRealArrays2D = .false.
            IF(PRESENT(code)) code = ASSERT_VALUES_DIFFER
         END IF 
         
      END FUNCTION isWithinToleranceTwoRealArrays2D
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoDouble(x,y,tol)  
         IMPLICIT NONE  
         DOUBLE PRECISION, INTENT(in) :: x,y,tol
         LOGICAL                      :: test
         
         IF ( x <= tol )     THEN
            test = ABS(x-y) <= tol
         ELSE
            test = ABS(x-y) <= tol*MAX(ABS(x),ABS(y))
         END IF 
         
         IF ( test )     THEN
            isWithinToleranceTwoDouble = .true.
         ELSE
            isWithinToleranceTwoDouble = .false.
         END IF 
         
      END FUNCTION isWithinToleranceTwoDouble    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoDoubleArrays1D(a,b,tol,code)
         IMPLICIT NONE  
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: a, b
         DOUBLE PRECISION, INTENT(IN)               :: tol
         INTEGER, INTENT(OUT), OPTIONAL             :: code
         
         isWithinToleranceTwoDoubleArrays1D = .true.
         IF(PRESENT(code)) code = ASSERT_SUCCESS
         
         IF ( SIZE(a) /= SIZE(b) )     THEN
            isWithinToleranceTwoDoubleArrays1D = .false.
            IF(PRESENT(code)) code = ASSERT_SIZE_DIFFERS
         ELSE IF(ANY(ABS(a-b) > tol*MAX(ABS(a),ABS(b))))     THEN
            isWithinToleranceTwoDoubleArrays1D = .false.
            IF(PRESENT(code)) code = ASSERT_VALUES_DIFFER
         END IF 
         
      END FUNCTION isWithinToleranceTwoDoubleArrays1D
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoDoubleArrays2D(a,b,tol,code)  
         IMPLICIT NONE  
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: a, b
         DOUBLE PRECISION, INTENT(IN)                 :: tol
         INTEGER, INTENT(OUT), OPTIONAL               :: code
         
         isWithinToleranceTwoDoubleArrays2D = .true.
         code = ASSERT_SUCCESS
         
         IF ( SIZE(a) /= SIZE(b) )     THEN
            isWithinToleranceTwoDoubleArrays2D = .false.
            IF(PRESENT(code)) code = ASSERT_SIZE_DIFFERS
         ELSE IF(ANY(ABS(a-b) > tol*MAX(ABS(a),ABS(b))))     THEN
            isWithinToleranceTwoDoubleArrays2D = .false.
            IF(PRESENT(code)) code = ASSERT_VALUES_DIFFER
         END IF 
         
      END FUNCTION isWithinToleranceTwoDoubleArrays2D
!@mark -
#ifdef _has_Quad
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isWithinToleranceTwoQuad(x,y,tol)  
         IMPLICIT NONE  
         REAL(KIND=SELECTED_REAL_KIND(QUAD_DIGITS)), INTENT(in)  :: x,y,tol
         LOGICAL                         :: test
         
         IF ( x == 0.0d0 )     THEN
            test = ABS(x-y) <= tol
         ELSE
            test = ABS(x-y) <= tol*MAX(ABS(x),ABS(y))
         END IF 
         
         IF ( test )     THEN
            isWithinToleranceTwoQuad = .true.
         ELSE
            isWithinToleranceTwoQuad = .false.
         END IF 
         
      END FUNCTION isWithinToleranceTwoQuad    
#endif
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isEqualString(s1,s2)
         IMPLICIT NONE
         CHARACTER(LEN=*) :: s1,s2
         
         isEqualString = .true.
         IF ( TRIM(s1) /= TRIM(s2) )     THEN
            isEqualString = .false. 
         END IF 
         
      END FUNCTION isEqualString
      
      END Module ComparisonsModule    