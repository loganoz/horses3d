!
! ////////////////////////////////////////////////////////////////////
!
! Defines constants for use by the spectral solver, including precision definitions. 
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
         REAL(KIND=RP), PARAMETER         :: DEG2RAD = PI / 180.0_RP

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

         integer, parameter               :: NDIM = 3, IX = 1, IY = 2, IZ = 3, IXY = 4, IXZ = 5, IYZ = 6 , IXYZ = 7
         !$acc declare copyin(NDIM, IX, IY, IZ, IXY, IXZ, IYZ, IXYZ)
         
         integer, parameter               :: DT_FIXED = 0
         integer, parameter               :: DT_DIFF  = 1
         integer, parameter               :: DT_CONV  = 2

         INTEGER, PARAMETER :: LEFT   = 1, RIGHT  = 2, TOP  = 2, BOTTOM  = 1
         INTEGER, PARAMETER :: FRONT  = 1, BACK   = 2
         !$acc declare copyin(LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)
         
         INTEGER, PARAMETER :: BC_STRING_LENGTH = 64

         CHARACTER(len=*), parameter   :: VERSION = "v0.8.9: Compressible Navier-Stokes physics now using conservative gradients."
         integer, protected            :: solver
   

         enum, bind(C)
            enumerator :: NAVIERSTOKES_SOLVER, INCNS_SOLVER, CAHNHILLIARD_SOLVER
            enumerator :: MULTIPHASE_SOLVER, NAVIERSTOKESSA_SOLVER, ACOUSTIC_SOLVER
            enumerator :: NO_OF_SOLVERS
            enumerator :: UNKNOWN_SOLVER = -1
         end enum

         enum, bind(C)
            enumerator :: INFLOW_BC = 1 , OUTFLOW_BC
            enumerator :: NOSLIPWALL_BC , FREESLIPWALL_BC
            enumerator :: PERIODIC_BC   , USERDEFINED_BC
         end enum
   
         character(len=BC_STRING_LENGTH), dimension(8)  :: implementedBCNames = [&
                   "inflow              ",  &
                   "outflow             ",  &
                   "noslipwall          ",  &
                   "freeslipwall        ",  &
                   "periodic            ",  &
                   "user-defined        ",  &
                   "manufacturedsol     ",  &
                   "msoutflowspecifyp   "]

         contains
            subroutine SetSolver(which)
               implicit none
               integer, intent(in)  :: which

               if ((solver < 0) .or. (solver >= NO_OF_SOLVERS)) then
                  print*, "Solver not recognized"
                  solver = UNKNOWN_SOLVER
               else
                  solver = which
               end if
   
            

            end subroutine SetSolver
         
      END MODULE SMConstants
