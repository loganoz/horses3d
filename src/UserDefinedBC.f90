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
         procedure         :: FlowState         => UserDefinedBC_FlowState
         procedure         :: FlowGradVars      => UserDefinedBC_FlowGradVars
         procedure         :: FlowNeumann       => UserDefinedBC_FlowNeumann
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
      function ConstructUserDefinedBC(bname)
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