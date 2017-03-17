!
! ////////////////////////////////////////////////////////////////////
!
!     SMConstants.F
!
!!
!!     Modification History:
!!       version 0.0 August 10, 2005 David A. Kopriva
!
!     MODULE SMConstants
!
!!        Defines constants for use by the spectral demonstaration
!!        routines, including precision definitions. 
!
!!    @author David A. Kopriva
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      MODULE SMConstants
         INTEGER      , PARAMETER, PRIVATE:: DIGITS        = 14                       ! # of desired digits
         INTEGER      , PARAMETER, PRIVATE:: SINGLE_DIGITS = 6                        ! # of desired digits
         INTEGER      , PARAMETER         :: RP = SELECTED_REAL_KIND( DIGITS )        ! Real Kind
         INTEGER      , PARAMETER         :: SP = SELECTED_REAL_KIND( SINGLE_DIGITS ) ! Single Real Kind
         INTEGER      , PARAMETER         :: CP = SELECTED_REAL_KIND( DIGITS )        ! Complex Kind
         REAL(KIND=RP), PARAMETER         :: PI = 3.141592653589793238462643_RP
         
         INTEGER, PARAMETER               :: FORWARD  = +1
         INTEGER, PARAMETER               :: BACKWARD = -1
         
         INTEGER, PARAMETER               :: STD_OUT = 6
         INTEGER, PARAMETER               :: STD_IN  = 5
         INTEGER, PARAMETER               :: LINE_LENGTH = 132
         INTEGER, PARAMETER               :: STRING_CONSTANT_LENGTH = 64
         
         COMPLEX(KIND=CP)                 :: ImgI = ( 0.0_RP, 1.0_RP) ! = SQRT(-1.0_RP)

#define errorMessage(UNIT) write(UNIT,'(A,A,A,I0,A)')   "Error in file ", __FILE__ , ", in line " , __LINE__ ,"."
#define stopMessage(UNIT)  write(UNIT,'(A,A,A,I0,A)') "Stopped in file ", __FILE__ , ", in line " , __LINE__ ,"."
         
      END MODULE SMConstants
