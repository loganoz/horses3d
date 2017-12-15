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

         character(len=1), parameter, private   :: SINGLE_CHARACTER = "a"
         integer      , parameter         :: SIZEOF_INT  = sizeof(integer)
         integer      , parameter         :: SIZEOF_RP   = RP
         integer      , parameter         :: SIZEOF_CHAR = sizeof(SINGLE_CHARACTER)

#if defined(ARCH_32_BITS)
         INTEGER      , PARAMETER         :: AddrInt = SELECTED_INT_KIND(9)
#else
         INTEGER      , PARAMETER         :: AddrInt = SELECTED_INT_KIND(18)
#endif
         
         INTEGER, PARAMETER               :: FORWARD  = +1
         INTEGER, PARAMETER               :: BACKWARD = -1
         
         INTEGER, PARAMETER               :: STD_OUT = 6
         INTEGER, PARAMETER               :: STD_IN  = 5
         INTEGER, PARAMETER               :: LINE_LENGTH = 132
         INTEGER, PARAMETER               :: STRING_CONSTANT_LENGTH = 64
         
         COMPLEX(KIND=CP), parameter      :: ImgI = ( 0.0_RP, 1.0_RP) ! = SQRT(-1.0_RP)

         INTEGER, PARAMETER :: LEFT   = 1, RIGHT  = 2, TOP  = 2, BOTTOM  = 1
         INTEGER, PARAMETER :: FRONT  = 1, BACK   = 2

         INTEGER, PARAMETER :: BC_STRING_LENGTH = 32

         CHARACTER(len=*), parameter   :: VERSION = "Development- v0.6.1"
         
      END MODULE SMConstants
