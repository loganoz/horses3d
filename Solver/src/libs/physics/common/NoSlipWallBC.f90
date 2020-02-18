!
!//////////////////////////////////////////////////////
!
!   @File:    NoSlipWallBC.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jul 25 15:26:42 2018
!   @Last revision date: Mon Feb 25 16:07:51 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 17d60e4e57235a57aa406023ebe4c26157bc211a
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
#ifdef FLOW
      real(kind=RP)     :: vWall(NDIM)
#endif
#ifdef CAHNHILLIARD
      real(kind=RP)     :: thetaw
#endif
      contains
         procedure         :: Destruct          => NoSlipWallBC_Destruct
         procedure         :: Describe          => NoSlipWallBC_Describe
#ifdef FLOW
         procedure         :: FlowState         => NoSlipWallBC_FlowState
         procedure         :: FlowNeumann       => NoSlipWallBC_FlowNeumann
#endif
#ifdef CAHNHILLIARD
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
#ifdef FLOW
         if ( bcdict % ContainsKey("wall velocity") ) then
            ConstructNoSlipWallBC % vWall = getRealArrayFromString( bcdict % StringValueForKey("wall velocity",&
                                                                                           LINE_LENGTH))    
            
         else
            ConstructNoSlipWallBC % vWall = 0.0_RP
         end if
#endif
#ifdef CAHNHILLIARD
         call GetValueWithDefault(bcdict, "contact angle", 90.0_RP, ConstructNoSlipWallBC % thetaw)

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
#ifdef FLOW
         write(STD_OUT,'(30X,A,A28,A,F10.2,A,F10.2,A,F10.2,A)') "->", ' Wall velocity: ',"[",self % vWall(1),",",self % vWall(2),",",self % vWall(3),"]"
#endif
#ifdef NAVIERSTOKES
         if ( self % isAdiabatic ) then
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Adiabatic"


         else
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Isothermal"
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall temperature: ', self % Twall
         end if
#endif
#ifdef CAHNHILLIARD
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall contact angle: ', self % thetaw
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

         Q(IRHOU:IRHOW) = 2.0_RP * Q(IRHO) * self % vWall - Q(IRHOU:IRHOW)
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
         real(kind=RP) :: invRho, invRho2, nablaT(NDIM), uDivRho(NDIM), vel_x(NDIM), vel_y(NDIM), vel_z(NDIM), u, v, w
         real(kind=RP) :: constA, constB, constC, constD
         REAL(KIND=RP) :: dTdn
         
         if (self % isAdiabatic) then
            
            invRho = 1._RP / Q(IRHO)
            invRho2 = invRho * invRho
            u = Q(IRHOU) * invRho
            v = Q(IRHOV) * invRho
            w = Q(IRHOW) * invRho
            
            uDivRho = [Q(IRHOU) , Q(IRHOV) , Q(IRHOW) ] * invRho2
            
            vel_x = invRho * U_x(IRHOU:IRHOW) - uDivRho * U_x(IRHO)
            vel_y = invRho * U_y(IRHOU:IRHOW) - uDivRho * U_y(IRHO)
            vel_z = invRho * U_z(IRHOU:IRHOW) - uDivRho * U_z(IRHO)
            
            constA = thermodynamics % gammaMinus1 * dimensionless % gammaM2
            constB = Q(IRHOE)*invRho2*U_x(IRHO) + u*vel_x(IX) + v*vel_x(IY) + w*vel_x(IZ)
            constC = Q(IRHOE)*invRho2*U_y(IRHO) + u*vel_y(IX) + v*vel_y(IY) + w*vel_y(IZ)
            constD = Q(IRHOE)*invRho2*U_z(IRHO) + u*vel_z(IX) + v*vel_z(IY) + w*vel_z(IZ)
            
            nablaT(IX) = constA * (invRho*U_x(IRHOE) - constB ) ! Inner dT/dx
            nablaT(IY) = constA * (invRho*U_y(IRHOE) - constC ) ! Inner dT/dy
            nablaT(IZ) = constA * (invRho*U_z(IRHOE) - constD ) ! Inner dT/dz
            
            dTdn = ( nablaT(IX)*nHat(IX) + nablaT(IY)*nHat(IY) + nablaT(IZ)*nHat(IZ) )
            
            nablaT(IX) = nablaT(IX) - 2.0_RP * dTdn * nHat(IX) ! Adiabatic reflection
            nablaT(IY) = nablaT(IY) - 2.0_RP * dTdn * nHat(IY) ! Adiabatic reflection
            nablaT(IZ) = nablaT(IZ) - 2.0_RP * dTdn * nHat(IZ) ! Adiabatic reflection
            
            U_x(IRHOE) = Q(IRHO) * ( nablaT(IX)/constA + constB )
            U_y(IRHOE) = Q(IRHO) * ( nablaT(IY)/constA + constC )
            U_z(IRHOE) = Q(IRHO) * ( nablaT(IZ)/constA + constD )
         end if

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
         real(kind=RP),       intent(inout) :: Q(NCONS)

         Q(INSRHOU:INSRHOW) = 2.0_RP * self % vWall - Q(INSRHOU:INSRHOW)

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
         real(kind=RP),       intent(inout) :: U_x(NCONS)
         real(kind=RP),       intent(inout) :: U_y(NCONS)
         real(kind=RP),       intent(inout) :: U_z(NCONS)
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
!        Subroutines for the multiphase solver
!        -------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef MULTIPHASE
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
         real(kind=RP),       intent(inout) :: Q(NCONS)

         Q(IMSQRHOU:IMSQRHOW) = 2.0_RP * self % vWall - Q(IMSQRHOU:IMSQRHOW)

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
         real(kind=RP),       intent(inout) :: U_x(NCONS)
         real(kind=RP),       intent(inout) :: U_y(NCONS)
         real(kind=RP),       intent(inout) :: U_z(NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
!
!
!        Remove the normal component of the chemical potential gradient
!        --------------------------------------------------------------
         gradUNorm =  nHat(1)*U_x(IGMU) + nHat(2)*U_y(IGMU)+ nHat(3)*U_z(IGMU)
         UTanx = U_x(IGMU) - gradUNorm*nHat(1)
         UTany = U_y(IGMU) - gradUNorm*nHat(2)
         UTanz = U_z(IGMU) - gradUNorm*nHat(3)
   
         U_x(IGMU) = UTanx - gradUNorm*nHat(1)
         U_y(IGMU) = UTany - gradUNorm*nHat(2)
         U_z(IGMU) = UTanz - gradUNorm*nHat(3)
!
!        Remove the normal component of the pressure gradient
!        ----------------------------------------------------
         gradUNorm =  nHat(1)*U_x(IMP) + nHat(2)*U_y(IMP)+ nHat(3)*U_z(IMP)
         UTanx = U_x(IMP) - gradUNorm*nHat(1)
         UTany = U_y(IMP) - gradUNorm*nHat(2)
         UTanz = U_z(IMP) - gradUNorm*nHat(3)
   
         U_x(IMP) = UTanx - gradUNorm*nHat(1)
         U_y(IMP) = UTany - gradUNorm*nHat(2)
         U_z(IMP) = UTanz - gradUNorm*nHat(3)

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
!
!        ---------------
!        Local variables
!        ---------------
!  
         real(kind=RP), parameter   :: MIN_ = 1.0e-1_RP
         real(kind=RP)              :: prod
         real(kind=RP)              :: prod12, prod13, prod23, c3

         prod = Q(1) * (1.0_RP - Q(1))

         if ( prod .le. MIN_ ) then
            prod = 0.0_RP
         end if

         U_x = -4.0_RP * multiphase % invEps * cos(DEG2RAD*self % thetaw) * nHat(1) * prod 
         U_y = -4.0_RP * multiphase % invEps * cos(DEG2RAD*self % thetaw) * nHat(2) * prod
         U_z = -4.0_RP * multiphase % invEps * cos(DEG2RAD*self % thetaw) * nHat(3) * prod

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
