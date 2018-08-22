!
!//////////////////////////////////////////////////////
!
!   @File:    InflowBC.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jul 25 15:26:42 2018
!   @Last revision date: Tue Aug 21 18:19:11 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: f09a57fa71226958a52c9f81e1bbaeb68e07bd37
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module InflowBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities, only: toLower, almostEqual
   use GenericBoundaryConditionClass
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
!
!  ******************
!  Public definitions
!  ******************
!
   public InflowBC_t
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
   type, extends(GenericBC_t) ::  InflowBC_t
#if defined(NAVIERSTOKES)
      real(kind=RP)              :: AoAPhi
      real(kind=RP)              :: AoATheta
      real(kind=RP)              :: v
      real(kind=RP)              :: rho
      real(kind=RP)              :: p
#endif
#if defined(INCNS)
      real(kind=RP)              :: AoAPhi
      real(kind=RP)              :: AoATheta
      real(kind=RP)              :: v
#if defined(CAHNHILLIARD)
      logical                    :: isLayered  = .false.
      logical                    :: isXLimited = .false.
      logical                    :: isYLimited = .false.
      logical                    :: isZLimited = .false.
      real(kind=RP)              :: xLim, yLim, zLim
      real(kind=RP)              :: phase1Vel
      real(kind=RP)              :: phase2Vel
#else
      real(kind=RP)              :: rho
#endif
#endif
      contains
         procedure         :: Destruct          => InflowBC_Destruct
         procedure         :: Describe          => InflowBC_Describe
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure         :: FlowState         => InflowBC_FlowState
         procedure         :: FlowNeumann       => InflowBC_FlowNeumann
#endif
#if defined(CAHNHILLIARD)
         procedure         :: PhaseFieldState   => InflowBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => InflowBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => InflowBC_ChemPotState
         procedure         :: ChemPotNeumann    => InflowBC_ChemPotNeumann
#endif
   end type InflowBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface InflowBC_t
      module procedure ConstructInflowBC
   end interface InflowBC_t
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
      function ConstructInflowBC(bname)
!
!        ********************************************************************
!        · Definition of the inflow boundary condition in the control file:
!              #define boundary bname
!                 type             = inflow
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
         type(InflowBC_t)             :: ConstructInflowBC
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
         interface
            subroutine PreprocessInputLine(line)
               implicit none
               character(len=*), intent(inout) :: line
            end subroutine PreprocessInputLine
         end interface

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructInflowBC % bname  = bname
         ConstructInflowBC % BCType = "inflow"

         call toLower(ConstructInflowBC % bname)
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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructInflowBC % bname)) 
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
         call GetValueWithDefault(bcdict, "pressure", refValues % p / dimensionless % gammaM2, ConstructInflowBC % p       )
         call GetValueWithDefault(bcdict, "density" , refValues % rho                        , ConstructInflowBC % rho     )
         call GetValueWithDefault(bcdict, "mach"    , dimensionless % Mach                   , ConstructInflowBC % v       )
         call GetValueWithDefault(bcdict, "aoaphi"  , refValues % AoAPhi                     , ConstructInflowBC % AoAPhi  )
         call GetValueWithDefault(bcdict, "aoatheta", refValues % AoATheta                   , ConstructInflowBC % AoATheta)

         ConstructInflowBC % p        = ConstructInflowBC % p / refValues % p
         ConstructInflowBC % rho      = ConstructInflowBC % rho / refValues % rho
         ConstructInflowBC % v        = ConstructInflowBC % v * sqrt(thermodynamics % gamma * ConstructInflowBC % p / ConstructInflowBC % rho)
         ConstructInflowBC % AoATheta = ConstructInflowBC % AoATheta * PI / 180.0_RP
         ConstructInflowBC % AoAPhi   = ConstructInflowBC % AoAPhi * PI / 180.0_RP

#elif defined(INCNS)
         call GetValueWithDefault(bcdict, "velocity", refValues % v, ConstructInflowBC % v)
         call GetValueWithDefault(bcdict, "aoaphi"  , refValues % AoAPhi                     , ConstructInflowBC % AoAPhi  )
         call GetValueWithDefault(bcdict, "aoatheta", refValues % AoATheta                   , ConstructInflowBC % AoATheta)

         ConstructInflowBC % v = ConstructInflowBC % v / refValues % v
         ConstructInflowBC % AoATheta = ConstructInflowBC % AoATheta * PI / 180.0_RP
         ConstructInflowBC % AoAPhi   = ConstructInflowBC % AoAPhi * PI / 180.0_RP

#if (!defined(CAHNHILLIARD))
         call GetValueWithDefault ( bcdict , "density" , refValues % rho , ConstructInflowBC % rho ) 
         ConstructInflowBC % rho = ConstructInflowBC % rho / refValues % rho
#else
!
!        *********************
!        Multiphase input data
!        *********************
!
         if (bcdict % ContainsKey("multiphase type (mixed/layered)") ) then   
            keyval = bcdict % StringValueForKey("multiphase type (mixed/layered)", LINE_LENGTH) 
            call tolower(keyval)
            if ( trim(keyval) .eq. "layered" ) then
               ConstructInflowBC % isLayered = .true.
            else
               ConstructInflowBC % isLayered = .false.
            end if

         else
            ConstructInflowBC % isLayered = .false.

         end if

         if (ConstructInflowBC % isLayered) then
!
!           Require for x, y, and/or z limits
!           ---------------------------------
            if ( bcdict % ContainsKey("interface x > x0") ) then
               ConstructInflowBC % isXLimited = .true.
               ConstructInflowBC % xLim = bcdict % DoublePrecisionValueForKey("interface x > x0")
            else
               ConstructInflowBC % isXLimited = .false.
            end if

            if ( bcdict % ContainsKey("interface y > y0") ) then
               ConstructInflowBC % isYLimited = .true.
               ConstructInflowBC % yLim = bcdict % DoublePrecisionValueForKey("interface y > y0")
            else
               ConstructInflowBC % isYLimited = .false.
            end if

            if ( bcdict % ContainsKey("interface z > z0") ) then
               ConstructInflowBC % isZLimited = .true.
               ConstructInflowBC % zLim = bcdict % DoublePrecisionValueForKey("interface z > z0")
            else
               ConstructInflowBC % isZLimited = .false.
            end if

            call GetValueWithDefault(bcdict, "phase 1 velocity", refValues % v, ConstructInflowBC % phase1Vel) 
            call GetValueWithDefault(bcdict, "phase 2 velocity", refValues % v, ConstructInflowBC % phase2Vel) 
   
            ConstructInflowBC % phase1Vel = ConstructInflowBC % phase1Vel / refValues % v
            ConstructInflowBC % phase2Vel = ConstructInflowBC % phase2Vel / refValues % v

         end if
#endif
#endif
         call bcdict % Destruct
         close(fid)
   
      end function ConstructInflowBC

      subroutine InflowBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(InflowBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "Inflow"
#if defined(NAVIERSTOKES)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Velocity: ', self % v * refValues % v
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Mach number: ', self % v / sqrt(thermodynamics % gamma * self % p / self % rho)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Pressure: ', self % p * refValues % p
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Density: ', self % rho * refValues % rho
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' AoaPhi: ', self % AoAPhi * 180.0_RP / PI
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' AoaTheta: ', self % AoATheta * 180.0_RP / PI
#elif defined(INCNS)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Velocity: ', self % v * refValues % v
#if (!defined(CAHNHILLIARD))
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Density: ', self % rho * refValues % rho
#endif
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' AoaPhi: ', self % AoAPhi * 180.0_RP / PI
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' AoaTheta: ', self % AoATheta * 180.0_RP / PI
#if defined(CAHNHILLIARD)
         if ( self % isLayered ) then
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Multiphase type: '," Layered"
         else
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Multiphase type: '," Mixed"
         end if
         
         if ( self % isLayered ) then
            if ( self % isXLimited ) write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Interface limits x>x0: ", self % xLim
            if ( self % isYLimited ) write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Interface limits y>y0: ", self % yLim
            if ( self % isZLimited ) write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Interface limits z>z0: ", self % zLim
            
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Phase 1 velocity: ", self % phase1Vel * refValues % v
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Phase 2 velocity: ", self % phase2Vel * refValues % v
         end if
#endif
#endif
         
      end subroutine InflowBC_Describe
!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine InflowBC_Destruct(self)
         implicit none
         class(InflowBC_t)    :: self

      end subroutine InflowBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine InflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(InflowBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: qq, u, v, w
         associate ( gammaM2 => dimensionless % gammaM2, &
                     gamma => thermodynamics % gamma )
         
         qq = self % v
         u  = qq*cos(self % AoAtheta)*COS(self % AoAphi)
         v  = qq*sin(self % AoAtheta)*COS(self % AoAphi)
         w  = qq*SIN(self % AoAphi)
         
         Q(1) = self % rho
         Q(2) = Q(1)*u
         Q(3) = Q(1)*v
         Q(4) = Q(1)*w
         Q(5) = self % p/(gamma - 1._RP) + 0.5_RP*Q(1)*(u**2 + v**2 + w**2)

         end associate

      end subroutine InflowBC_FlowState

      subroutine InflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(InflowBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
         real(kind=RP),       intent(inout) :: U_x(NGRAD)
         real(kind=RP),       intent(inout) :: U_y(NGRAD)
         real(kind=RP),       intent(inout) :: U_z(NGRAD)

         INTEGER :: k
         REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz

         DO k = 1, NGRAD
            gradUNorm =  nHat(1)*U_x(k) + nHat(2)*U_y(k) + nHat(3)*U_z(k)
            UTanx = U_x(k) - gradUNorm*nHat(1)
            UTany = U_y(k) - gradUNorm*nHat(2)
            UTanz = U_z(k) - gradUNorm*nHat(3)
      
            U_x(k) = UTanx - gradUNorm*nHat(1)
            U_y(k) = UTany - gradUNorm*nHat(2)
            U_z(k) = UTanz - gradUNorm*nHat(3)
         END DO

      end subroutine InflowBC_FlowNeumann
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
      subroutine InflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rho, u, v, w, vel

#if (!defined(CAHNHILLIARD))
         u = self % v * cos(self % AoAtheta) * cos(self % AoAphi)
         v = self % v * sin(self % AoAtheta) * cos(self % AoAphi)
         w = self % v * sin(self % AoAphi)

         Q(INSRHO)  = self % rho
         Q(INSRHOU) = Q(INSRHO)*u
         Q(INSRHOV) = Q(INSRHO)*v
         Q(INSRHOW) = Q(INSRHO)*w
#else
         if ( self % isLayered ) then
            if (self % isXLimited) then
               if ( x(IX) > self % xLim ) then
                  rho = dimensionless % rho(1)
                  vel = self % phase1Vel 
               else
                  rho = dimensionless % rho(2)
                  vel = self % phase2Vel 
               end if
            elseif ( self % isYLimited) then
               if ( x(IY) > self % yLim ) then
                  rho = dimensionless % rho(1)
                  vel = self % phase1Vel 
               else
                  rho = dimensionless % rho(2)
                  vel = self % phase2Vel 
               end if
            elseif ( self % isZLimited) then
               if ( x(IZ) > self % zLim ) then
                  rho = dimensionless % rho(1)
                  vel = self % phase1Vel 
               else
                  rho = dimensionless % rho(2)
                  vel = self % phase2Vel 
               end if
            end if
         else
            rho = 0.5_RP * (dimensionless % rho(1) + dimensionless % rho(2))
            vel = self % v
         end if

         u = vel * cos(self % AoAtheta) * cos(self % AoAphi)
         v = vel * sin(self % AoAtheta) * cos(self % AoAphi)
         w = vel * sin(self % AoAphi)

         Q(INSRHO)  = rho
         Q(INSRHOU) = Q(INSRHO)*u
         Q(INSRHOV) = Q(INSRHO)*v
         Q(INSRHOW) = Q(INSRHO)*w
#endif

      end subroutine InflowBC_FlowState

      subroutine InflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
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

      end subroutine InflowBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard: all do--nothing in Inflows
!        ----------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)
      subroutine InflowBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine InflowBC_PhaseFieldState

      subroutine InflowBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
         real(kind=RP),       intent(inout) :: U_x(NCOMP)
         real(kind=RP),       intent(inout) :: U_y(NCOMP)
         real(kind=RP),       intent(inout) :: U_z(NCOMP)
      end subroutine InflowBC_PhaseFieldNeumann

      subroutine InflowBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine InflowBC_ChemPotState

      subroutine InflowBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
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

      end subroutine InflowBC_ChemPotNeumann
#endif
end module InflowBCClass
