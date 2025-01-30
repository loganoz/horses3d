#include "Includes.h"
module GenericBoundaryConditionClass
   use SMConstants
   use PhysicsStorage
   use FTValueDictionaryClass, only: FTValueDictionary
   use VariableConversion
   use FluidData
   use HexMeshClass
   use ZoneClass
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
   public GenericBC_t, GetValueWithDefault
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
         procedure         :: FlowState         => GenericBC_FlowState
         procedure         :: CreateDeviceData  => GenericBC_CreateDeviceData
         procedure         :: FlowGradVars      => GenericBC_FlowGradVars
         procedure         :: FlowNeumann       => GenericBC_FlowNeumann
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
      subroutine StateForEqn(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t), intent(in)            :: zone

#ifndef CAHNHILLIARD
        call self % FlowState(mesh, zone)

#else
         select case(self % currentEqn)
#ifdef FLOW
         case(NS_BC)
            call self % FlowState(mesh, zone)
#endif
         case(C_BC)
            call self % PhaseFieldState(mesh, zone)
   
         case(MU_BC)
            call self % ChemPotState(mesh, zone)

         end select
#endif

      end subroutine StateForEqn

      subroutine GradVarsForEqn(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t), intent(in)            :: zone

#ifndef CAHNHILLIARD
         call self % FlowGradVars(mesh, zone)

#else
         select case(self % currentEqn)
#ifdef FLOW
         case(NS_BC)
            call self % FlowGradVars(mesh, zone)
#endif
         case(C_BC)
            call self % PhaseFieldState(mesh, zone)
   
         case(MU_BC)
            call self % ChemPotState(mesh, zone)

         end select
#endif

      end subroutine GradVarsForEqn

      subroutine NeumannForEqn(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t), intent(in)            :: zone

#ifdef CAHNHILLIARD
         select case(self % currentEqn)
         case(C_BC)
            call self % PhaseFieldNeumann(mesh, zone)
   
         case(MU_BC)
            call self % ChemPotNeumann(mesh, zone)

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

      subroutine GenericBC_CreateDeviceData(self)
         implicit none 
         class(GenericBC_t), intent(in)    :: self
         !$acc enter data copyin(self)
      end subroutine GenericBC_CreateDeviceData

      subroutine GenericBC_FlowState(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t), intent(in)            :: zone
!         integer,             intent(in)     :: zoneID 
         
      end subroutine GenericBC_FlowState

      subroutine GenericBC_FlowGradVars(self, mesh, zone)
!
!        **************************************************************
!              Computes the set of gradient variables U* at the wall
!        **************************************************************
!
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         type(HexMesh), intent(inout)       :: mesh
         type(Zone_t), intent(in)           :: zone

!         integer,                 intent(in)    :: zoneID 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: Q_aux(NCONS),Q(NCONS)
         real(kind=RP)  :: u_int(NGRAD), u_star(NGRAD)
         real(kind=RP)  :: e_int, U1
         real(kind=RP)  :: invRho, rho
         integer        :: i,j,zonefID,fID
   
   
         !$acc parallel loop gang present(mesh, zone) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, Q_aux, u_star, u_int)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q     = mesh % faces(fID) % storage(1) % Q(:,i,j)
               Q_aux = mesh % faces(fID) % storage(2) % Q(:,i,j)

#ifdef MULTIPHASE

               ! Set the chemical potential to the interior
               ! ------------------------------------------

               call mGradientVariables(NCONS, NGRAD, Q, u_int, mesh % faces(fID) % storage(1) % rho(i,j))
               u_int(IGMU) = Q(IGMU) 

               rho = dimensionless % rho(1) * Q_aux(IMC) + dimensionless % rho(2) * (1.0_RP-Q_aux(IMC))
               rho = min(max(rho, dimensionless % rho_min),dimensionless % rho_max)
               call mGradientVariables(NCONS, NGRAD, Q_aux, u_star, rho)
               u_star(IGMU) = u_int(IGMU)           
#else
               call NSGradientVariables_STATE(NCONS, NGRAD, Q    , u_int)
               call NSGradientVariables_STATE(NCONS, NGRAD, Q_aux, u_star)
#endif

               u_star = 0.5_RP* (u_star + u_int)
               
               mesh % faces(fID) % storage(1) % unStar(:,1,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(1,i,j) * mesh % faces(fID) % geom % jacobian(i,j)
               mesh % faces(fID) % storage(1) % unStar(:,2,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(2,i,j) * mesh % faces(fID) % geom % jacobian(i,j)    
               mesh % faces(fID) % storage(1) % unStar(:,3,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(3,i,j) * mesh % faces(fID) % geom % jacobian(i,j)

            enddo ; enddo

         enddo
         !$acc end parallel loop
      end subroutine GenericBC_FlowGradVars

!      subroutine GenericBC_FlowGradVars(self, mesh, zoneID)
!         implicit none
!         class(GenericBC_t),  intent(in)    :: self
!         type(HexMesh),       intent(in)    :: mesh
!         integer,             intent(in)    :: zoneID 
!
!        ---------------
!        Local variables
!        ---------------
!
!         real(kind=RP)  :: Q_aux(NCONS), U_aux(NGRAD), rho
!         Q_aux = Q
!         U_aux = U
!         call self % FlowState(x,t,nHat,Q_aux)
!#ifdef MULTIPHASE
!
!        Set the chemical potential to the interior
!        ------------------------------------------
!         rho = dimensionless % rho(1) * Q_aux(IMC) + dimensionless % rho(2) * (1.0_RP-Q_aux(IMC))
!         rho = min(max(rho, dimensionless % rho_min),dimensionless % rho_max)
!         call GetGradients(NCONS,NGRAD,Q_aux, U_aux, rho)
!         U_aux(IGMU) = U(IGMU)
!#else
!         call GetGradients(NCONS,NGRAD,Q_aux, U_aux)
!#endif
!         U = 0.5_RP * (U_aux + U)
!      end subroutine GenericBC_FlowGradVars

      subroutine GenericBC_FlowNeumann(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         type(HexMesh),       intent(inout)    :: mesh
         type(Zone_t), intent(in)               :: zone
      end subroutine GenericBC_FlowNeumann
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
      subroutine GenericBC_PhaseFieldState(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)    :: self
         type(HexMesh),       intent(inout)    :: mesh
         type(Zone_t), intent(in)               :: zone
      end subroutine GenericBC_PhaseFieldState

      subroutine GenericBC_PhaseFieldNeumann(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t),        intent(in)     :: zone
      end subroutine GenericBC_PhaseFieldNeumann

      subroutine GenericBC_ChemPotState(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t),        intent(in)     :: zone
      end subroutine GenericBC_ChemPotState

      subroutine GenericBC_ChemPotNeumann(self, mesh, zone)
         implicit none
         class(GenericBC_t),  intent(in)     :: self
         type(HexMesh),       intent(inout)  :: mesh
         type(Zone_t),        intent(in)     :: zone
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



end module GenericBoundaryConditionClass