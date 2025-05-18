#include "Includes.h"
module GenericBoundaryConditionClass
   use SMConstants
   use PhysicsStorage
   use FTValueDictionaryClass, only: FTValueDictionary
   use VariableConversion, only: GetGradientValues_f
   use FluidData
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
   public NS_BC, C_BC, MU_BC
!
!  ******************
!  Public definitions
!  ******************
!
   public GenericBC_t, GetValueWithDefault, CheckIfBoundaryNameIsContained, CheckIfStlNameIsContained
!
!  ****************************
!  Static variables definitions
!  ****************************
!
   enum, bind(C)
      enumerator :: NONE_BC
      enumerator :: NS_BC
      enumerator :: C_BC, MU_BC
   end enum
!
!  ****************
!  Class definition
!  ****************
!
   type GenericBC_t
      logical                    :: constructed = .false.
      character(len=LINE_LENGTH) :: bname
      character(len=LINE_LENGTH) :: BCType
      integer                    :: currentEqn = 1
      contains
         procedure         :: Destruct          => GenericBC_Destruct
         procedure         :: Describe          => GenericBC_Describe
         procedure         :: GetPeriodicPair   => GenericBC_GetPeriodicPair
#ifdef FLOW
         procedure         :: FlowState           => GenericBC_FlowState
         procedure         :: FlowGradVars        => GenericBC_FlowGradVars
         procedure         :: FlowNeumann         => GenericBC_FlowNeumann
         procedure         :: FlowState_VPIBM     => GenericBC_FlowState_VPIBM
         procedure         :: FlowState_HOIBM     => GenericBC_FlowState_HOIBM
         procedure         :: FlowStateWeak_HOIBM => GenericBC_FlowStateWeak_HOIBM
         procedure         :: FlowGradVars_HOIBM  => GenericBC_FlowGradVars_HOIBM
         procedure         :: FlowNeumann_HOIBM   => GenericBC_FlowNeumann_HOIBM
         procedure         :: FlowStateMoving_IBM => GenericBC_FlowStateMoving_IBM
         procedure         :: PositionMoving_IBM  => GenericBC_PositionMoving_IBM
#endif
#ifdef CAHNHILLIARD
         procedure         :: PhaseFieldState   => GenericBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => GenericBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => GenericBC_ChemPotState
         procedure         :: ChemPotNeumann    => GenericBC_ChemPotNeumann
#endif
         procedure         :: StateForEqn
         procedure         :: GradVarsForEqn
         procedure         :: NeumannForEqn
   end type GenericBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface GenericBC_t
      module procedure ConstructGenericBC
   end interface GenericBC_t
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
      function ConstructGenericBC()
!
!        ************************************************************
!        This is a default constructor (i.e. without input arguments)
!        ************************************************************
!
         implicit none
         type(GenericBC_t)  :: ConstructGenericBC
   
         print*, "Default Boundary condition class can not be constructed."
         errorMessage(STD_OUT)
         error stop

      end function ConstructGenericBC
!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine GenericBC_Destruct(self)
         implicit none
         class(GenericBC_t)    :: self

         if (.not. self % constructed) return
      
      end subroutine GenericBC_Destruct
!
!//////////////////////////////////////////////////////////////////////////////////
!
!        This subroutine allows to call the State/Neumann for generic equations
!        ----------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine StateForEqn(self, nEqn, x, t, nHat, Q)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         integer,             intent(in)    :: nEqn
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(nEqn)

#ifndef CAHNHILLIARD
         call self % FlowState(x, t, nHat, Q)

#else
         select case(self % currentEqn)
#ifdef FLOW
         case(NS_BC)
            call self % FlowState(x, t, nHat, Q)
#endif
         case(C_BC)
            call self % PhaseFieldState(x, t, nHat, Q)
   
         case(MU_BC)
            call self % ChemPotState(x, t, nHat, Q)

         end select
#endif

      end subroutine StateForEqn

      subroutine GradVarsForEqn(self, nEqn, nGradEqn, x, t, nHat, Q, U, GetGradients)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         integer,             intent(in)    :: nEqn
         integer,             intent(in)    :: nGradEqn
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(nEqn)
         real(kind=RP),       intent(inout) :: U(nGradEqn)
         procedure(GetGradientValues_f)     :: GetGradients

#ifndef CAHNHILLIARD
         call self % FlowGradVars(x, t, nHat, Q, U, GetGradients)

#else
         select case(self % currentEqn)
#ifdef FLOW
         case(NS_BC)
            call self % FlowGradVars(x, t, nHat, Q, U, GetGradients)
#endif
         case(C_BC)
            call self % PhaseFieldState(x, t, nHat, U)
   
         case(MU_BC)
            call self % ChemPotState(x, t, nHat, U)

         end select
#endif

      end subroutine GradVarsForEqn

      subroutine NeumannForEqn(self, nEqn, nGradEqn, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         integer,             intent(in)    :: nEqn
         integer,             intent(in)    :: nGradEqn
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(nEqn)
         real(kind=RP),       intent(in)    :: U_x(nGradEqn)
         real(kind=RP),       intent(in)    :: U_y(nGradEqn)
         real(kind=RP),       intent(in)    :: U_z(nGradEqn)
         real(kind=RP),       intent(inout) :: flux(nEqn)

#ifdef CAHNHILLIARD
         select case(self % currentEqn)
         case(C_BC)
            call self % PhaseFieldNeumann(x, t, nHat, Q, U_x, U_y, U_z, flux)
   
         case(MU_BC)
            call self % ChemPotNeumann(x, t, nHat, Q, U_x, U_y, U_z, flux)

         case default
            print*, "Unexpected equation choice"
            errorMessage(STD_OUT)
            error stop

         end select
#else
         print*, "This function is only supported for Cahn-Hilliard"
         errorMessage(STD_OUT)
         error stop
#endif
      end subroutine NeumannForEqn
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Navier--Stokes equations
!        ----------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef FLOW
      subroutine GenericBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
      end subroutine GenericBC_FlowState

      subroutine GenericBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)     :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: Q_aux(NCONS), U_aux(NGRAD), rho

         Q_aux = Q
         U_aux = U

         call self % FlowState(x,t,nHat,Q_aux)

#ifdef MULTIPHASE
!
!        Set the chemical potential to the interior
!        ------------------------------------------
         rho = dimensionless % rho(1) * Q_aux(IMC) + dimensionless % rho(2) * (1.0_RP-Q_aux(IMC))
         rho = min(max(rho, dimensionless % rho_min),dimensionless % rho_max)
         call GetGradients(NCONS,NGRAD,Q_aux, U_aux, rho)
         U_aux(IGMU) = U(IGMU)
#else
         call GetGradients(NCONS,NGRAD,Q_aux, U_aux)
#endif

         U = 0.5_RP * (U_aux + U)

      end subroutine GenericBC_FlowGradVars

      subroutine GenericBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NGRAD)
         real(kind=RP),       intent(in)    :: U_y(NGRAD)
         real(kind=RP),       intent(in)    :: U_z(NGRAD)
         real(kind=RP),       intent(inout) :: flux(NCONS)
      end subroutine GenericBC_FlowNeumann

      subroutine GenericBC_FlowState_VPIBM( self, Q, x, Qsb )   
         implicit none
         class(GenericBC_t),    intent(in)    :: self
         real(kind=RP),         intent(in)    :: Q(NCONS)
         real(kind=RP),         intent(in)    :: x(NDIM)
         real(kind=RP),         intent(inout) :: Qsb(NCONS)  
      end subroutine GenericBC_FlowState_VPIBM

      subroutine GenericBC_FlowStateMoving_IBM( self, Q, x, t, dt, cL, cD, Qsb )  
         implicit none 
         class(GenericBC_t),       intent(in)    :: self
         real(kind=RP),            intent(in)    :: Q(NCONS)
         real(kind=RP),            intent(inout) :: x(NDIM)
         real(kind=RP),            intent(in)    :: t, dt
         real(kind=RP),            intent(in)    :: cL, cD 
         real(kind=RP),            intent(inout) :: Qsb(NCONS)
      end subroutine GenericBC_FlowStateMoving_IBM

      subroutine GenericBC_PositionMoving_IBM( self, x, t, dt, cL, cD )  
         implicit none 
         class(GenericBC_t),      intent(in)    :: self
         real(kind=RP),           intent(inout) :: x(NDIM)
         real(kind=RP),           intent(in)    :: t, dt
         real(kind=RP),           intent(in)    :: cL, cD 
      end subroutine GenericBC_PositionMoving_IBM

      subroutine GenericBC_FlowState_HOIBM( self, Q, xb, xsb, nodes, N, x, time, nHat, Qsb )  
         implicit none
         class(GenericBC_t),    intent(in)    :: self
         real(kind=RP),         intent(inout) :: Q(NCONS,0:N)
         real(kind=RP),         intent(in)    :: xb, xsb, nodes(0:N), x(NDIM), time, nHat(NDIM)
         integer,               intent(in)    :: N
         real(kind=RP),         intent(inout) :: Qsb(NCONS)
      end subroutine GenericBC_FlowState_HOIBM

      subroutine GenericBC_FlowStateWeak_HOIBM( self, QIn, Qsb, nHat, Qsb_weak )    
         implicit none
         class(GenericBC_t),    intent(in)    :: self
         real(kind=RP),         intent(in)    :: QIn(NCONS), Qsb(NCONS)  
         real(kind=RP),         intent(in)    :: nHat(NDIM)
         real(kind=RP),         intent(inout) :: Qsb_weak(NCONS)
      end subroutine GenericBC_FlowStateWeak_HOIBM

      subroutine GenericBC_FlowGradVars_HOIBM( self, QIn, Qsb, x, time, nHat, Usb, GetGradients )    
         implicit none
         class(GenericBC_t),     intent(in)    :: self
         real(kind=RP),          intent(in)    :: QIn(NCONS), Qsb(NCONS), x(NDIM), time, nHat(NDIM)
         real(kind=RP),          intent(inout) :: Usb(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
      end subroutine GenericBC_FlowGradVars_HOIBM

      subroutine GenericBC_FlowNeumann_HOIBM( self, QIn, UIn_x, UIn_y, UIn_z, Qsb, x, time, nHat, flux )    
         implicit none
         class(GenericBC_t),    intent(in)    :: self
         real(kind=RP),         intent(in)    :: QIn(NCONS), Qsb(NCONS), x(NDIM), time, nHat(NDIM)
         real(kind=RP),         intent(in)    :: UIn_x(NGRAD), UIn_y(NGRAD), UIn_z(NGRAD)
         real(kind=RP),         intent(inout) :: flux(NCONS)
      end subroutine GenericBC_FlowNeumann_HOIBM
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard
!        ------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef CAHNHILLIARD
      subroutine GenericBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine GenericBC_PhaseFieldState

      subroutine GenericBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)
      end subroutine GenericBC_PhaseFieldNeumann

      subroutine GenericBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine GenericBC_ChemPotState

      subroutine GenericBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)
      end subroutine GenericBC_ChemPotNeumann
#endif

      subroutine GenericBC_Describe(self)
         implicit none
         class(GenericBC_t),  intent(in)  :: self
      end subroutine GenericBC_Describe

      subroutine GenericBC_GetPeriodicPair(self, bname)
!
!        *****************************************
!        Only for periodic BCs, empty for the rest
!        *****************************************
!
         implicit none
         class(GenericBC_t),  intent(in)  :: self
         character(len=*), intent(out) :: bname
      end subroutine GenericBC_GetPeriodicPair
         
!
!////////////////////////////////////////////////////////////////////////////
!
!        Some utilities
!        --------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine GetValueWithDefault(controlVariables, keyword, default, val)
         implicit none
         type(FTValueDictionary),    intent(in)  :: controlVariables
         character(len=*),           intent(in)  :: keyword
         real(kind=RP),              intent(in)  :: default
         real(kind=RP),              intent(out) :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         if ( controlVariables % ContainsKey(trim(keyword)) ) then
            val = controlVariables % DoublePrecisionValueForKey(trim(keyword))
         else
            val = default
         end if
   
      end subroutine GetValueWithDefault

      logical function CheckIfBoundaryNameIsContained(line, bname)
!
!        **********************************************************************
!        We will allow several boundary names defined within a single region.
!           They should be separated by two underscores (e.g. bname1__bname2)
!        **********************************************************************
!
         implicit none
         character(len=*), intent(in)  :: line
         character(len=*), intent(in)  :: bname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: pos1, pos2
         character(len=LINE_LENGTH) :: str, curName
         logical     :: found
         

         str = "#define boundary"
!
!        Exit if not a #define boundary sentinel
!        ---------------------------------------
         if ( line(1:len_trim(str)) .ne. trim(str) ) then
            CheckIfBoundaryNameIsContained = .false.
            return
         end if
!
!        Get only the boundary names
!        ---------------------------
         str = adjustl(line(len_trim(str)+1:len_trim(line)))
!
!        Exit if empty
!        -------------
         if (len_trim(str) .eq. 0) then
            CheckIfBoundaryNameIsContained = .false.
            return
         end if
!
!        Check the entries
!        -----------------
         pos1 = 1
         pos2 = len_trim(str)
         found = .false.
         do while(pos2 .ne. 0)
            pos2 = index(str(pos1:len_trim(str)),"__") + pos1 - 1
            if ( pos2 .eq. pos1-1 ) then
!
!              Last iteration
!              --------------
               curName = str(pos1:len_trim(str))
               pos2 = 0
            else
               curName = str(pos1:pos2-1)
            end if

            if ( trim(curName) .eq. trim(bname) ) then
               found = .true.
               CheckIfBoundaryNameIsContained = .true.
               return      
            else
               pos1 = pos2+2
            end if
         end do

         CheckIfBoundarynameIsContained = .false.

      end function CheckIfBoundaryNameIsContained

      logical function CheckIfStlNameIsContained(line, stlname)

         implicit none
         character(len=*), intent(in)  :: line
         character(len=*), intent(in)  :: stlname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: pos1, pos2
         character(len=LINE_LENGTH) :: str, curName
         logical     :: found
               
      
         str = "#define stl"
!
!        Exit if not a #define boundary sentinel
!        ---------------------------------------
         if ( line(1:len_trim(str)) .ne. trim(str) ) then
            CheckIfStlNameIsContained = .false.
            return
         end if
!
!        Get only the boundary names
!        ---------------------------
         str = adjustl(line(len_trim(str)+1:len_trim(line)))
!
!        Exit if empty
!        -------------
         if (len_trim(str) .eq. 0) then
            CheckIfStlNameIsContained = .false.
            return
         end if
!
!        Check the entries
!        -----------------
         pos1 = 1
         pos2 = len_trim(str)
         found = .false.
         do while(pos2 .ne. 0)
            pos2 = index(str(pos1:len_trim(str)),"__") + pos1 - 1
            if ( pos2 .eq. pos1-1 ) then
!
!              Last iteration
!              --------------
               curName = str(pos1:len_trim(str))
               pos2 = 0
            else
               curName = str(pos1:pos2-1)
            end if

            if ( trim(curName) .eq. trim(stlname) ) then
               found = .true.
               CheckIfStlNameIsContained = .true.
               return      
            else
               pos1 = pos2+2
            end if
         end do

         CheckIfStlNameIsContained = .false.
      
      end function CheckIfStlNameIsContained

end module GenericBoundaryConditionClass