#include "Includes.h"
module UserDefinedBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities, only: toLower, almostEqual
   use GenericBoundaryConditionClass
   use FluidData
   use VariableConversion, only: GetGradientValues_f
#if defined(NAVIERSTOKES)
   use VariableConversion, only: pressure
#endif
   use PolynomialInterpAndDerivsModule 
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
!
!  ****************
!  Public variables
!  ****************
!
!
!  ******************
!  Public definitions
!  ******************
!
   public UserDefinedBC_t
!
!  ****************************
!  Static variables definitions
!  ****************************
!
!
!  ****************
!  Class definition
!  ****************
!
   type, extends(GenericBC_t) ::  UserDefinedBC_t
      integer     :: udf_no
      contains
         procedure         :: Destruct          => UserDefinedBC_Destruct
#ifdef FLOW
         procedure         :: FlowState           => UserDefinedBC_FlowState
         procedure         :: FlowGradVars        => UserDefinedBC_FlowGradVars
         procedure         :: FlowNeumann         => UserDefinedBC_FlowNeumann
         procedure         :: FlowState_HOIBM     => UserDefinedBC_FlowState_HOIBM
         procedure         :: FlowStateWeak_HOIBM => UserDefinedBC_FlowStateWeak_HOIBM
         procedure         :: FlowGradVars_HOIBM  => UserDefinedBC_FlowGradVars_HOIBM
         procedure         :: FlowNeumann_HOIBM   => UserDefinedBC_FlowNeumann_HOIBM
         procedure         :: FlowStateMoving_IBM => UserDefinedBC_FlowStateMoving_IBM
         procedure         :: PositionMoving_IBM  => UserDefinedBC_PositionMoving_IBM
#endif
#if defined(CAHNHILLIARD)
         procedure         :: PhaseFieldState   => UserDefinedBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => UserDefinedBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => UserDefinedBC_ChemPotState
         procedure         :: ChemPotNeumann    => UserDefinedBC_ChemPotNeumann
#endif
   end type UserDefinedBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface UserDefinedBC_t
      module procedure ConstructUserDefinedBC
   end interface UserDefinedBC_t
!
!  *******************
!  Function prototypes
!  *******************
!
   abstract interface
   end interface
!
!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      function ConstructUserDefinedBC(bname,isIBM)
!
!        ********************************************************************
!        · Definition of the user-defined boundary condition in the control file:
!              #define boundary bname
!                 type             = user-defined
!                 velocity         = #value        (only in incompressible NS)
!                 Mach number      = #value        (only in compressible NS)
!                 AoAPhi           = #value
!                 AoATheta         = #value
!                 density          = #value        (only in monophase)
!                 pressure         = #value        (only in compressible NS)
!                 multiphase type  = mixed/layered
!                 phase 1 layer x  > #value
!                 phase 1 layer y  > #value
!                 phase 1 layer z  > #value
!                 phase 1 velocity = #value
!                 phase 2 velocity = #value
!              #end
!        ********************************************************************
!
         implicit none
         type(UserDefinedBC_t)             :: ConstructUserDefinedBC
         character(len=*), intent(in) :: bname
         logical, optional, intent(in) :: isIBM
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: fid, io
         character(len=LINE_LENGTH) :: boundaryHeader
         character(len=LINE_LENGTH) :: currentLine
         character(len=LINE_LENGTH) :: keyword, keyval
         logical                    :: inside
         type(FTValueDIctionary)    :: bcdict

         if( present(isIBM) ) return 

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructUserDefinedBC % BCType = "user-defined"
         ConstructUserDefinedBC % bname  = bname
         call toLower(ConstructUserDefinedBC % bname)

         write(boundaryHeader,'(A,A)') "#define boundary ",trim(bname)
         call toLower(boundaryHeader)
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do 
            read(fid, '(A)', iostat=io) currentLine

            IF(io .ne. 0 ) EXIT

            call PreprocessInputLine(currentLine)
            call toLower(currentLine)

            if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructUserDefinedBC % bname)) 
            end if
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit

               call bcdict % InitWithSize(16)

               keyword  = ADJUSTL(GetKeyword(currentLine))
               keyval   = ADJUSTL(GetValueAsString(currentLine))
               call ToLower(keyword)
      
               call bcdict % AddValueForKey(keyval, trim(keyword))

            end if

         end do
!
!        Analyze the gathered data
!        -------------------------
         if ( bcdict % ContainsKey("udf number") ) then
            ConstructUserDefinedBC % udf_no = bcdict % IntegerValueForKey("udf number")
         else
            ConstructUserDefinedBC % udf_no = 1
         end if

         close(fid)
         call bcdict % Destruct
   
      end function ConstructUserDefinedBC

      subroutine UserDefinedBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(UserDefinedBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "UserDefined"
         write(STD_OUT,'(30X,A,A28,I0)') "->", " UD Function number: ", self % udf_no
         
      end subroutine UserDefinedBC_Describe



!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine UserDefinedBC_Destruct(self)
         implicit none
         class(UserDefinedBC_t)    :: self

      end subroutine UserDefinedBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Navier--Stokes equations
!        ----------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef FLOW
      subroutine UserDefinedBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(UserDefinedBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         interface
            subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage
               use FluidData
               implicit none
               real(kind=RP)  :: x(NDIM)
               real(kind=RP)  :: t
               real(kind=RP)  :: nHat(NDIM)
               real(kind=RP)  :: Q(NCONS)
               type(Thermodynamics_t), intent(in)  :: thermodynamics_
               type(Dimensionless_t),  intent(in)  :: dimensionless_
               type(RefValues_t),      intent(in)  :: refValues_
            end subroutine UserDefinedState1
         end interface

         select case(self % udf_no)
         case(1)
            call UserDefinedState1(x, t, nHat, Q, thermodynamics, dimensionless, refValues)
         case default
            print*, "Unrecognized UDF number for boundary", self % bname
         end select
   
      end subroutine UserDefinedBC_FlowState

      subroutine UserDefinedBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
         implicit none
         class(UserDefinedBC_t),   intent(in) :: self
         real(kind=RP),       intent(in)      :: x(NDIM)
         real(kind=RP),       intent(in)      :: t
         real(kind=RP),       intent(in)      :: nHat(NDIM)
         real(kind=RP),       intent(in)      :: Q(NCONS)
         real(kind=RP),       intent(inout)   :: U(NGRAD)
         procedure(GetGradientValues_f)       :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         interface
            subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage
               use FluidData
               use VariableConversion, only: GetGradientValues_f
               implicit none
               real(kind=RP), intent(in)          :: x(NDIM)
               real(kind=RP), intent(in)          :: t
               real(kind=RP), intent(in)          :: nHat(NDIM)
               real(kind=RP), intent(in)          :: Q(NCONS)
               real(kind=RP), intent(inout)       :: U(NGRAD)
               procedure(GetGradientValues_f)     :: GetGradients
               type(Thermodynamics_t), intent(in) :: thermodynamics_
               type(Dimensionless_t),  intent(in) :: dimensionless_
               type(RefValues_t),      intent(in) :: refValues_
            end subroutine UserDefinedGradVars1
         end interface

         select case(self % udf_no)
         case(1)
            call UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics, dimensionless, refValues)
         case default
            print*, "Unrecognized UDF number for boundary", self % bname
         end select
   
      end subroutine UserDefinedBC_FlowGradVars

      subroutine UserDefinedBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(UserDefinedBC_t),   intent(in) :: self
         real(kind=RP),       intent(in)      :: x(NDIM)
         real(kind=RP),       intent(in)      :: t
         real(kind=RP),       intent(in)      :: nHat(NDIM)
         real(kind=RP),       intent(in)      :: Q(NCONS)
         real(kind=RP),       intent(in)      :: U_x(NGRAD)
         real(kind=RP),       intent(in)      :: U_y(NGRAD)
         real(kind=RP),       intent(in)      :: U_z(NGRAD)
         real(kind=RP),       intent(inout)   :: flux(NCONS)
         interface
            subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
            end subroutine UserDefinedNeumann1
         end interface

         select case(self % udf_no)
         case(1)
            call UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics, dimensionless, refValues)
         case default
            print*, "Unrecognized UDF number for boundary", self % bname
         end select

      end subroutine UserDefinedBC_FlowNeumann

      subroutine UserDefinedBC_FlowStateMoving_IBM( self, Q, x, dt, cL, cD, Qsb )  
         implicit none 
         class(UserDefinedBC_t),   intent(in)    :: self
         real(kind=RP),            intent(in)    :: Q(NCONS)
         real(kind=RP),            intent(inout) :: x(NDIM)
         real(kind=RP),            intent(in)    :: dt
         real(kind=RP),            intent(in)    :: cL, cD 
         real(kind=RP),            intent(inout) :: Qsb(NCONS)

         logical       :: updatePosition, GetVelocity
         real(kind=RP) :: V(NDIM), P 
         interface
            subroutine UserDefinedIBMKinematicsNS( x, V, cL, cD, dt, refValues_, updatePosition, GetVelocity )
               use SMConstants
               use FluidData
               use PhysicsStorage
               IMPLICIT NONE
               real(kind=RP),           intent(inout) :: x(NDIM), V(NDIM)
               real(kind=RP),           intent(in)    :: dt
               real(kind=RP),           intent(in)    :: cL, cD 
               type(RefValues_t),       intent(in)    :: refValues_
               logical,                 intent(in)    :: GetVelocity, UpdatePosition
            end subroutine UserDefinedIBMKinematicsNS
         end interface
#if defined(NAVIERSTOKES)
         associate(gammaMinus1 => thermodynamics% gammaMinus1 )
         
         call UserDefinedIBMKinematicsNS( x=x, V=V, cL=cL, cD=cD, dt=dt, refValues_=refValues, updatePosition=.false., GetVelocity=.true. )

         P                = pressure(Q)
         Qsb(IRHO)        = Q(IRHO) 
         Qsb(IRHOU:IRHOW) = Q(IRHO)*V
         Qsb(IRHOE)       = P/gammaMinus1 + 0.5_RP * Q(IRHO) * sum(V*V)

         end associate 
#endif
      end subroutine UserDefinedBC_FlowStateMoving_IBM

      subroutine UserDefinedBC_PositionMoving_IBM( self, x, dt, cL, cD )  
         implicit none 
         class(UserDefinedBC_t),  intent(in)    :: self
         real(kind=RP),           intent(inout) :: x(NDIM)
         real(kind=RP),           intent(in)    :: dt
         real(kind=RP),           intent(in)    :: cL, cD 

         logical       :: updatePosition, GetVelocity
         real(kind=RP) :: V(NDIM)
         interface
            subroutine UserDefinedIBMKinematicsNS( x, V, cL, cD, dt, refValues_, updatePosition, GetVelocity )
               use SMConstants
               use FluidData
               use PhysicsStorage
               IMPLICIT NONE
               real(kind=RP),           intent(inout) :: x(NDIM), V(NDIM)
               real(kind=RP),           intent(in)    :: dt
               real(kind=RP), optional, intent(in)    :: cL, cD 
               type(RefValues_t),       intent(in)    :: refValues_
               logical,                 intent(in)    :: GetVelocity, UpdatePosition
            end subroutine UserDefinedIBMKinematicsNS
         end interface
         
         call UserDefinedIBMKinematicsNS( x=x, V=V, cL=cL, cD=cD, dt=dt, refValues_=refValues, updatePosition=.true., GetVelocity=.false. )

      end subroutine UserDefinedBC_PositionMoving_IBM

      subroutine UserDefinedBC_FlowState_HOIBM( self, Q, xb, xsb, nodes, N, x, time, nHat, Qsb )   
         implicit none
         class(UserDefinedBC_t), intent(in)    :: self
         real(kind=RP),          intent(inout) :: Q(NCONS,0:N)
         real(kind=RP),          intent(in)    :: xb, xsb, nodes(0:N), x(NDIM), time, nHat(NDIM)
         integer,                intent(in)    :: N
         real(kind=RP),          intent(inout) :: Qsb(NCONS) 

         real(kind=RP) :: Q_(NCONS), lj(0:N), w(0:N), den, dQ(NCONS)
         integer       :: i     
         interface
            subroutine UserDefinedState1(x, time, nHat, Q_, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage
               use FluidData
               implicit none
               real(kind=RP)  :: x(NDIM)
               real(kind=RP)  :: time
               real(kind=RP)  :: nHat(NDIM)
               real(kind=RP)  :: Q_(NCONS)
               type(Thermodynamics_t), intent(in)  :: thermodynamics_
               type(Dimensionless_t),  intent(in)  :: dimensionless_
               type(RefValues_t),      intent(in)  :: refValues_
            end subroutine UserDefinedState1
         end interface

         Qsb = 0.0_RP 

         do i = 0, N 
            lj(i) = LagrangeInterpolatingPolynomial( i, xsb, N, nodes )
         end do 

         call UserDefinedState1(x, time, nHat, Q(:,0), thermodynamics, dimensionless, refValues)

         do i = 0, N 
            Qsb = Qsb + lj(i) * Q(:,i)
         end do 

      end subroutine UserDefinedBC_FlowState_HOIBM

      subroutine UserDefinedBC_FlowStateWeak_HOIBM( self, QIn, Qsb, Qsb_weak )    
         implicit none
         class(UserDefinedBC_t), intent(in)    :: self
         real(kind=RP),          intent(in)    :: QIn(NCONS), Qsb(NCONS)  
         real(kind=RP),          intent(inout) :: Qsb_weak(NCONS)  

         real(kind=RP) :: rho, invRho, Vsb(NDIM), V(NDIM)

         Qsb_weak = Qsb

      end subroutine UserDefinedBC_FlowStateWeak_HOIBM

      subroutine UserDefinedBC_FlowGradVars_HOIBM( self, QIn, Qsb, x, time, nHat, Usb, GetGradients )    
         implicit none
         class(UserDefinedBC_t),  intent(in)    :: self
         real(kind=RP),           intent(in)    :: QIn(NCONS), Qsb(NCONS), x(NDIM), time, nHat(NDIM)
         real(kind=RP),           intent(inout) :: Usb(NGRAD)
         procedure(GetGradientValues_f)         :: GetGradients 

         interface
            subroutine UserDefinedGradVars1(x, time, nHat, Q, Usb, GetGradients, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage
               use FluidData
               use VariableConversion, only: GetGradientValues_f
               implicit none
               real(kind=RP), intent(in)          :: x(NDIM)
               real(kind=RP), intent(in)          :: time
               real(kind=RP), intent(in)          :: nHat(NDIM)
               real(kind=RP), intent(in)          :: Q(NCONS)
               real(kind=RP), intent(inout)       :: Usb(NGRAD)
               procedure(GetGradientValues_f)     :: GetGradients
               type(Thermodynamics_t), intent(in) :: thermodynamics_
               type(Dimensionless_t),  intent(in) :: dimensionless_
               type(RefValues_t),      intent(in) :: refValues_
            end subroutine UserDefinedGradVars1
         end interface
      end subroutine UserDefinedBC_FlowGradVars_HOIBM

      subroutine UserDefinedBC_FlowNeumann_HOIBM( self, QIn, UIn_x, UIn_y, UIn_z, Qsb, x, time, nHat, flux )    
         implicit none
         class(UserDefinedBC_t), intent(in)    :: self
         real(kind=RP),         intent(in)    :: QIn(NCONS), Qsb(NCONS), x(NDIM), time, nHat(NDIM)
         real(kind=RP),         intent(in)    :: UIn_x(NGRAD), UIn_y(NGRAD), UIn_z(NGRAD)
         real(kind=RP),         intent(inout) :: flux(NCONS)

      end subroutine UserDefinedBC_FlowNeumann_HOIBM
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard: all do--nothing in UserDefineds
!        ----------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)
      subroutine UserDefinedBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(UserDefinedBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine UserDefinedBC_PhaseFieldState

      subroutine UserDefinedBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(UserDefinedBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP

      end subroutine UserDefinedBC_PhaseFieldNeumann

      subroutine UserDefinedBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(UserDefinedBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine UserDefinedBC_ChemPotState

      subroutine UserDefinedBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(UserDefinedBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP

      end subroutine UserDefinedBC_ChemPotNeumann
#endif
end module UserDefinedBCClass