!
!//////////////////////////////////////////////////////
!
!   @File:    NoSlipWallBC.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jul 25 15:26:42 2018
!   @Last revision date: Thu Oct 18 16:09:47 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: f0ca5b23053e717fbb5fcc06b6de56d366b37b53
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module NoSlipWallBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString ,PreprocessInputLine
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   use Utilities, only: toLower, almostEqual
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
   public NoSlipWallBC_t
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
   type, extends(GenericBC_t) ::  NoSlipWallBC_t
#if defined(NAVIERSTOKES)
      logical           :: isAdiabatic
      real(kind=RP)     :: Twall
      real(kind=RP)     :: ewall       ! Wall internal energy
      real(kind=RP)     :: kWallType
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
      real(kind=RP)     :: vWall(NDIM)
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP)     :: thetaw
#endif
      contains
         procedure         :: Destruct          => NoSlipWallBC_Destruct
         procedure         :: Describe          => NoSlipWallBC_Describe
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure         :: FlowState         => NoSlipWallBC_FlowState
         procedure         :: FlowNeumann       => NoSlipWallBC_FlowNeumann
#endif
#if defined(CAHNHILLIARD)
         procedure         :: PhaseFieldState   => NoSlipWallBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => NoSlipWallBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => NoSlipWallBC_ChemPotState
         procedure         :: ChemPotNeumann    => NoSlipWallBC_ChemPotNeumann
#endif
   end type NoSlipWallBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface NoSlipWallBC_t
      module procedure ConstructNoSlipWallBC
   end interface NoSlipWallBC_t
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
      function ConstructNoSlipWallBC(bname)
!
!        ********************************************************************
!        · Definition of the noSlipWall boundary condition in the control file:
!              #define boundary bname
!                 type             = noSlipWall
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
         type(NoSlipWallBC_t)             :: ConstructNoSlipWallBC
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

         ConstructNoSlipWallBC % BCType = "noslipwall"
         ConstructNoSlipWallBC % bname  = bname
         call toLower(ConstructNoSlipWallBC % bname)

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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructNoSlipWallBC % bname)) 
            end if
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit

               keyword  = ADJUSTL(GetKeyword(currentLine))
               keyval   = ADJUSTL(GetValueAsString(currentLine))
               call ToLower(keyword)
      
               call bcdict % AddValueForKey(keyval, trim(keyword))

            end if

         end do
!
!        Analyze the gathered data
!        -------------------------
#if defined(NAVIERSTOKES)
         if ( bcdict % ContainsKey("wall type (adiabatic/isothermal)") ) then
            keyval = bcdict % StringValueForKey("wall type (adiabatic/isothermal)", LINE_LENGTH)
            call tolower(keyval)
         else
            keyval = "adiabatic"
         end if

         if ( trim(keyval) .eq. "adiabatic" ) then
            ConstructNoSlipWallBC % isAdiabatic = .true.
            ConstructNoSlipWallBC % kWallType      = 0.0_RP   ! This is to avoid an "if" when setting the wall temp 
            ConstructNoSlipWallBC % ewall = 0.0_RP
            ConstructNoSlipWallBC % Twall = 0.0_RP
         else
            ConstructNoSlipWallBC % isAdiabatic = .false.
            call GetValueWithDefault(bcdict, "wall temperature" , refValues % T, ConstructNoSlipWallBC % Twall     )
            ConstructNoSlipWallBC % ewall = ConstructNoSlipWallBC % Twall / (refValues % T*thermodynamics % gammaMinus1*dimensionless % gammaM2)
            ConstructNoSlipWallBC % kWallType = 1.0_RP
         end if
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
         if ( bcdict % ContainsKey("wall velocity") ) then
            ConstructNoSlipWallBC % vWall = getRealArrayFromString( bcdict % StringValueForKey("wall velocity",&
                                                                                           LINE_LENGTH))    
            
         else
            ConstructNoSlipWallBC % vWall = 0.0_RP
         end if
#endif
#if defined(CAHNHILLIARD)
         call GetValueWithDefault(bcdict, "contact angle", 0.0_RP, ConstructNoSlipWallBC % thetaw)
#endif

         close(fid)
         call bcdict % Destruct
   
      end function ConstructNoSlipWallBC

      subroutine NoSlipWallBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(NoSlipWallBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "NoSlipWall"
#if defined(NAVIERSTOKES) || defined(INCNS)
         write(STD_OUT,'(30X,A,A28,A,F10.2,A,F10.2,A,F10.2,A)') "->", ' Wall velocity: ',"[",self % vWall(1),",",self % vWall(2),",",self % vWall(3),"]"
#endif
#if defined(NAVIERSTOKES) 
         if ( self % isAdiabatic ) then
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Adiabatic"


         else
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Isothermal"
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall temperature: ', self % Twall * refValues % T
         end if
#endif
#if defined(CAHNHILLIARD)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall contact angle coef: ', self % thetaw
#endif
         
      end subroutine NoSlipWallBC_Describe

!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine NoSlipWallBC_Destruct(self)
         implicit none
         class(NoSlipWallBC_t)    :: self

      end subroutine NoSlipWallBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine NoSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Internal energy is either the interior state for 
!              adiabatic walls, or the imposed for isothermal.
!              eBC = eInt + kWallType (eIso - eInt)
!              where kWallType = 0 for adiabatic and 1 for isothermal.        
!        *************************************************************
!
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(inout) :: Q(NCONS)

         Q(IRHOU:IRHOW) = 2.0_RP * self % vWall - Q(IRHOU:IRHOW)
         Q(IRHOE) = Q(IRHOE) + self % kWallType*(Q(IRHO)*self % ewall-Q(IRHOE))

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
!
!        ***********************************************************
!           Change the sign of the temperature gradient only
!           if the adiabatic wall BC is specified.
!        ***********************************************************
!
         implicit none
         class(NoSlipWallBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
         real(kind=RP),       intent(inout) :: U_x(NGRAD)
         real(kind=RP),       intent(inout) :: U_y(NGRAD)
         real(kind=RP),       intent(inout) :: U_z(NGRAD)
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP) :: dTdn

         dTdn = (U_x(IGT)*nHat(IX) + U_y(IGT)*nHat(IY) + U_z(IGT)*nHat(IZ))
         U_x(IGT) = U_x(IGT) - (1.0_RP - self % kWallType) * 2.0_RP * dTdn * nHat(IX)
         U_y(IGT) = U_y(IGT) - (1.0_RP - self % kWallType) * 2.0_RP * dTdn * nHat(IY) 
         U_z(IGT) = U_z(IGT) - (1.0_RP - self % kWallType) * 2.0_RP * dTdn * nHat(IZ) 

      end subroutine NoSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for incompressible Navier--Stokes equations
!        -------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(INCNS)
      subroutine NoSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!

         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NINC)

         Q(INSRHOU:INSRHOW) = 2.0_RP * self % vWall - Q(INSRHOU:INSRHOW)

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NINC)
         real(kind=RP),       intent(inout) :: U_x(NINC)
         real(kind=RP),       intent(inout) :: U_y(NINC)
         real(kind=RP),       intent(inout) :: U_z(NINC)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
!
!
!        Remove the normal component of the density gradient
!        ---------------------------------------------------
         gradUNorm =  nHat(1)*U_x(INSRHO) + nHat(2)*U_y(INSRHO)+ nHat(3)*U_z(INSRHO)
         UTanx = U_x(INSRHO) - gradUNorm*nHat(1)
         UTany = U_y(INSRHO) - gradUNorm*nHat(2)
         UTanz = U_z(INSRHO) - gradUNorm*nHat(3)
   
         U_x(INSRHO) = UTanx - gradUNorm*nHat(1)
         U_y(INSRHO) = UTany - gradUNorm*nHat(2)
         U_z(INSRHO) = UTanz - gradUNorm*nHat(3)
!
!        Remove the normal component of the pressure gradient
!        ----------------------------------------------------
         gradUNorm =  nHat(1)*U_x(INSP) + nHat(2)*U_y(INSP)+ nHat(3)*U_z(INSP)
         UTanx = U_x(INSP) - gradUNorm*nHat(1)
         UTany = U_y(INSP) - gradUNorm*nHat(2)
         UTanz = U_z(INSP) - gradUNorm*nHat(3)
   
         U_x(INSP) = UTanx - gradUNorm*nHat(1)
         U_y(INSP) = UTany - gradUNorm*nHat(2)
         U_z(INSP) = UTanz - gradUNorm*nHat(3)


      end subroutine NoSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard
!        ------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)
      subroutine NoSlipWallBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine NoSlipWallBC_PhaseFieldState

      subroutine NoSlipWallBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
         real(kind=RP),       intent(inout) :: U_x(NCOMP)
         real(kind=RP),       intent(inout) :: U_y(NCOMP)
         real(kind=RP),       intent(inout) :: U_z(NCOMP)

         U_x = self % thetaw * nHat(1)
         U_y = self % thetaw * nHat(2)
         U_z = self % thetaw * nHat(3)

      end subroutine NoSlipWallBC_PhaseFieldNeumann

      subroutine NoSlipWallBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine NoSlipWallBC_ChemPotState

      subroutine NoSlipWallBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
         real(kind=RP),       intent(inout) :: U_x(NCOMP)
         real(kind=RP),       intent(inout) :: U_y(NCOMP)
         real(kind=RP),       intent(inout) :: U_z(NCOMP)

         U_x = 0.0_RP
         U_y = 0.0_RP
         U_z = 0.0_RP

      end subroutine NoSlipWallBC_ChemPotNeumann
#endif
end module NoSlipWallBCClass
