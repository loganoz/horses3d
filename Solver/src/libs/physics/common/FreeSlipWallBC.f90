!
!//////////////////////////////////////////////////////
!
!   @File:    FreeSlipWallBC.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jul 25 15:26:41 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module FreeSlipWallBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use FluidData
   use FileReadingUtilities, only: getArrayFromString
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
   public FreeSlipWallBC_t
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
   type, extends(GenericBC_t) ::  FreeSlipWallBC_t
#if defined(NAVIERSTOKES)
      logical           :: isAdiabatic
      real(kind=RP)     :: Twall
      real(kind=RP)     :: ewall       ! Wall internal energy
      real(kind=RP)     :: kWallType
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP)     :: thetaw
#endif
      contains
         procedure         :: Destruct          => FreeSlipWallBC_Destruct
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure         :: FlowState         => FreeSlipWallBC_FlowState
         procedure         :: FlowNeumann       => FreeSlipWallBC_FlowNeumann
#endif
#if defined(CAHNHILLIARD)
         procedure         :: PhaseFieldState   => FreeSlipWallBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => FreeSlipWallBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => FreeSlipWallBC_ChemPotState
         procedure         :: ChemPotNeumann    => FreeSlipWallBC_ChemPotNeumann
#endif
   end type FreeSlipWallBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface FreeSlipWallBC_t
      module procedure ConstructFreeSlipWallBC
   end interface FreeSlipWallBC_t
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
      function ConstructFreeSlipWallBC(bname)
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
         type(FreeSlipWallBC_t)             :: ConstructFreeSlipWallBC
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

         ConstructFreeSlipWallBC % BCType = "freeslipwall"
         ConstructFreeSlipWallBC % bname  = bname

         write(boundaryHeader,'(A,A)') "#define boundary ",trim(bname)
         call toLower(boundaryHeader)
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do 
            write(fid, '(A)', iostat=io) currentLine

            IF(io .ne. 0 ) EXIT

            call PreprocessInputLine(currentLine)

            if ( trim(currentLine) .eq. trim(boundaryHeader) ) then
               inside = .true.
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
#if defined(NAVIERSTOKES)
         if ( bcdict % ContainsKey("wall type (adiabatic/isothermal)") ) then
            keyval = bcdict % StringValueForKey("wall type (adiabatic/isothermal)", LINE_LENGTH)
            call tolower(keyval)
         else
            keyval = "adiabatic"
         end if

         if ( trim(keyval) .eq. "adiabatic" ) then
            ConstructFreeSlipWallBC % isAdiabatic = .true.
            ConstructFreeSlipWallBC % kWallType      = 0.0_RP   ! This is to avoid an "if" when setting the wall temp 
         else
            ConstructFreeSlipWallBC % isAdiabatic = .false.
            call GetValueWithDefault(bcdict, "wall temperature" , refValues % T, ConstructFreeSlipWallBC % Twall     )
            ConstructFreeSlipWallBC % ewall = ConstructFreeSlipWallBC % Twall / (refValues % T*thermodynamics % gammaMinus1*dimensionless % gammaM2)
            ConstructFreeSlipWallBC % kWallType = 1.0_RP
         end if
#endif
#if defined(CAHNHILLIARD)
         call GetValueWithDefault(bcdict, "contact angle", 0.0_RP, ConstructFreeSlipWallBC % thetaw)
#endif

         close(fid)
         call bcdict % Destruct
   
      end function ConstructFreeSlipWallBC
!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine FreeSlipWallBC_Destruct(self)
         implicit none
         class(FreeSlipWallBC_t)    :: self

      end subroutine FreeSlipWallBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to v_interior - 2v_normal:
!                 -> Maintains tangential speed.
!                 -> Removes normal speed (weakly).
!           · Internal energy is either the interior state for 
!              adiabatic walls, or the imposed for isothermal.
!              eBC = eInt + kWallType (eIso - eInt)
!              where kWallType = 0 for adiabatic and 1 for isothermal.        
!        *************************************************************
!
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         real(kind=RP),           intent(in)    :: x(NDIM)
         real(kind=RP),           intent(in)    :: t
         real(kind=RP),           intent(in)    :: nHat(NDIM)
         real(kind=RP),           intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: qNorm

         qNorm = nHat(IX) * Q(IRHOU) + nHat(IY) * Q(IRHOV) + nHat(IZ) * Q(IRHOW)

         Q(IRHOU:IRHOW) = Q(IRHOU:IRHOW) - 2.0_RP * qNorm * nHat

         Q(IRHOE) = Q(IRHOE) + self % kWallType*(Q(IRHO)*self % ewall-Q(IRHOE))

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
!
!        ***********************************************************
!           Remove all normal gradients
!        ***********************************************************
!
         implicit none
         class(FreeSlipWallBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
         real(kind=RP),       intent(inout) :: U_x(NGRAD)
         real(kind=RP),       intent(inout) :: U_y(NGRAD)
         real(kind=RP),       intent(inout) :: U_z(NGRAD)
!
!        ---------------
!        Local Variables
!        ---------------
!   
         REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
         INTEGER       :: k
!   
         DO k = 1, NGRAD
            gradUNorm =  nHat(1)*U_x(k) + nHat(2)*U_y(k) + nHat(3)*U_z(k)
            UTanx = U_x(k) - gradUNorm*nHat(1)
            UTany = U_y(k) - gradUNorm*nHat(2)
            UTanz = U_z(k) - gradUNorm*nHat(3)
   
            U_x(k) = UTanx - gradUNorm*nHat(1)
            U_y(k) = UTany - gradUNorm*nHat(2)
            U_z(k) = UTanz - gradUNorm*nHat(3)
         END DO

      end subroutine FreeSlipWallBC_FlowNeumann
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
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
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
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: vn
!
!        -----------------------------------------------
!        Generate the external flow along the face, that
!        represents a solid wall.
!        -----------------------------------------------
!
         vn = sum(Q(INSRHOU:INSRHOW)*nHat)

         Q(INSRHO)  = Q(INSRHO)
         Q(INSRHOU) = Q(INSRHOU) - 2.0_RP * vn * nHat(IX)
         Q(INSRHOV) = Q(INSRHOV) - 2.0_RP * vn * nHat(IY)
         Q(INSRHOW) = Q(INSRHOW) - 2.0_RP * vn * nHat(IZ)
         Q(INSP)    = Q(INSP)

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
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


      end subroutine FreeSlipWallBC_FlowNeumann
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
      subroutine FreeSlipWallBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine FreeSlipWallBC_PhaseFieldState

      subroutine FreeSlipWallBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
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

      end subroutine FreeSlipWallBC_PhaseFieldNeumann

      subroutine FreeSlipWallBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine FreeSlipWallBC_ChemPotState

      subroutine FreeSlipWallBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
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

      end subroutine FreeSlipWallBC_ChemPotNeumann
#endif
end module FreeSlipWallBCClass
