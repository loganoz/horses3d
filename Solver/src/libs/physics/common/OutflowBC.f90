!
!//////////////////////////////////////////////////////
!
!   @File:    OutflowBC.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jul 25 15:26:43 2018
!   @Last revision date: Fri Jul 27 18:59:33 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 52dc2d6e64f5c93de205ce8bb5d283128f3c5e11
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module OutflowBCClass
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
   public OutflowBC_t
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
   type, extends(GenericBC_t) ::  OutflowBC_t
#if defined(NAVIERSTOKES) || defined(INCNS)
      real(kind=RP)     :: pExt
#endif
      contains
         procedure         :: Destruct          => OutflowBC_Destruct
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure         :: FlowState         => OutflowBC_FlowState
         procedure         :: FlowNeumann       => OutflowBC_FlowNeumann
#endif
#if defined(CAHNHILLIARD)
         procedure         :: PhaseFieldState   => OutflowBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => OutflowBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => OutflowBC_ChemPotState
         procedure         :: ChemPotNeumann    => OutflowBC_ChemPotNeumann
#endif
   end type OutflowBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface OutflowBC_t
      module procedure ConstructOutflowBC
   end interface OutflowBC_t
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
      function ConstructOutflowBC(bname)
!
!        ********************************************************************
!        · Definition of the outflow boundary condition in the control file:
!              #define boundary bname
!                 type             = outflow
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
         type(OutflowBC_t)             :: ConstructOutflowBC
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

         ConstructOutflowBC % BCType = "outflow"
         ConstructOutflowBC % bname  = bname
         call toLower(ConstructOutflowBC % bname)

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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructOutflowBC % bname)) 
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
#if defined(NAVIERSTOKES) || defined(INCNS)
         if ( bcdict % ContainsKey("pressure") ) then
            ConstructOutflowBC % pExt = bcdict % DoublePrecisionValueForKey("pressure")
         else if ( bcdict % ContainsKey("dp") ) then
#if defined(NAVIERSTOKES)
            ConstructOutflowBC % pExt = refValues % p / dimensionless % gammaM2 - bcdict % DoublePrecisionValueForKey("dp")
#elif defined(INCNS)
            ConstructOutflowBC % pExt = 0.0_RP - bcdict % DoublePrecisionValueForKey("dp")
#endif
         else
#if defined(NAVIERSTOKES)
            ConstructOutflowBC % pExt = refValues % p / dimensionless % gammaM2 
#elif defined(INCNS)
            ConstructOutflowBC % pExt = 0.0_RP
#endif
         end if

         ConstructOutflowBC % pExt = ConstructOutflowBC % pExt / refValues % p
#endif

         close(fid)
         call bcdict % Destruct
   
      end function ConstructOutflowBC
!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine OutflowBC_Destruct(self)
         implicit none
         class(OutflowBC_t)    :: self

      end subroutine OutflowBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine OutflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!   
!        ---------------
!        Local Variables
!        ---------------
!   
         REAL(KIND=RP) :: qDotN, qTanx, qTany, qTanz, p, a, a2
         REAL(KIND=RP) :: rPlus, entropyConstant, u, v, w, rho, normalMachNo
!         
         associate ( gammaMinus1 => thermodynamics % gammaMinus1, &
                     gamma => thermodynamics % gamma )
         
         qDotN = (nHat(1)*Q(2) + nHat(2)*Q(3) + nHat(3)*Q(4))/Q(1)
         qTanx = Q(2)/Q(1) - qDotN*nHat(1)
         qTany = Q(3)/Q(1) - qDotN*nHat(2)
         qTanz = Q(4)/Q(1) - qDotN*nHat(3)
         
         p            = gammaMinus1*( Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1) )
         a2           = gamma*p/Q(1)
         a            = SQRT(a2)
         normalMachNo = ABS(qDotN/a)
         
         IF ( normalMachNo <= 1.0_RP )     THEN
!   
!           -------------------------------
!           Quantities coming from upstream
!           -------------------------------
!   
            rPlus           = qDotN + 2.0_RP*a/gammaMinus1
            entropyConstant = p - a2*Q(1)
!   
!           ----------------
!           Resolve solution
!           ----------------
!   
            rho   = -(entropyConstant - self % pExt)/a2
            a     = SQRT(gamma*self % pExt/rho)
            qDotN = rPlus - 2.0_RP*a/gammaMinus1
            u     = qTanx + qDotN*nHat(1)
            v     = qTany + qDotN*nHat(2)
            w     = qTanz + qDotN*nHat(3)
            
            Q(1) = rho
            Q(2) = rho*u
            Q(3) = rho*v
            Q(4) = rho*w
            Q(5) = self % pExt/gammaMinus1 + 0.5_RP*rho*(u*u + v*v + w*w)
           
         END IF
   
      end associate

      end subroutine OutflowBC_FlowState

      subroutine OutflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(OutflowBC_t),   intent(in)    :: self
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

      end subroutine OutflowBC_FlowNeumann
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
      subroutine OutflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: u, v, w, theta, phi, un

         un = Q(INSRHOU)*nHat(IX) + Q(INSRHOV)*nHat(IY) + Q(INSRHOW)*nHat(IZ)
         
         if ( un .ge. -1.0e-4_RP ) then
         
            Q(INSP) = self % pExt

         else

            theta = refValues % AoATheta * PI / 180.0_RP
            phi   = refValues % AoAPhi   * PI / 180.0_RP
   
            u = cos(theta) * cos(phi)
            v = sin(theta) * cos(phi)
            w = sin(phi)

            Q(INSRHOU) = Q(INSRHO)*u
            Q(INSRHOV) = Q(INSRHO)*v
            Q(INSRHOW) = Q(INSRHO)*w

         end if

      end subroutine OutflowBC_FlowState

      subroutine OutflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NINC)
         real(kind=RP),       intent(inout) :: U_x(NINC)
         real(kind=RP),       intent(inout) :: U_y(NINC)
         real(kind=RP),       intent(inout) :: U_z(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: drhodn, dudn, dvdn, dwdn

         drhodn = (U_x(INSRHO)*nHat(IX) + U_y(INSRHO)*nHat(IY) + U_z(INSRHO)*nHat(IZ))
         U_x(INSRHO) = U_x(INSRHO) - 2.0_RP * drhodn * nHat(IX)
         U_y(INSRHO) = U_y(INSRHO) - 2.0_RP * drhodn * nHat(IY) 
         U_z(INSRHO) = U_z(INSRHO) - 2.0_RP * drhodn * nHat(IZ) 

         dudn = (U_x(INSRHOU)*nHat(IX) + U_y(INSRHOU)*nHat(IY) + U_z(INSRHOU)*nHat(IZ))
         U_x(INSRHOU) = U_x(INSRHOU) - 2.0_RP * dudn * nHat(IX)
         U_y(INSRHOU) = U_y(INSRHOU) - 2.0_RP * dudn * nHat(IY) 
         U_z(INSRHOU) = U_z(INSRHOU) - 2.0_RP * dudn * nHat(IZ) 

         dvdn = (U_x(INSRHOV)*nHat(IX) + U_y(INSRHOV)*nHat(IY) + U_z(INSRHOV)*nHat(IZ))
         U_x(INSRHOV) = U_x(INSRHOV) - 2.0_RP * dvdn * nHat(IX)
         U_y(INSRHOV) = U_y(INSRHOV) - 2.0_RP * dvdn * nHat(IY) 
         U_z(INSRHOV) = U_z(INSRHOV) - 2.0_RP * dvdn * nHat(IZ) 

         dwdn = (U_x(INSRHOW)*nHat(IX) + U_y(INSRHOW)*nHat(IY) + U_z(INSRHOW)*nHat(IZ))
         U_x(INSRHOW) = U_x(INSRHOW) - 2.0_RP * dwdn * nHat(IX)
         U_y(INSRHOW) = U_y(INSRHOW) - 2.0_RP * dwdn * nHat(IY) 
         U_z(INSRHOW) = U_z(INSRHOW) - 2.0_RP * dwdn * nHat(IZ)

      end subroutine OutflowBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard: all do--nothing in Outflows
!        ----------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)
      subroutine OutflowBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine OutflowBC_PhaseFieldState

      subroutine OutflowBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
         real(kind=RP),       intent(inout) :: U_x(NCOMP)
         real(kind=RP),       intent(inout) :: U_y(NCOMP)
         real(kind=RP),       intent(inout) :: U_z(NCOMP)
      end subroutine OutflowBC_PhaseFieldNeumann

      subroutine OutflowBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine OutflowBC_ChemPotState

      subroutine OutflowBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
         real(kind=RP),       intent(inout) :: U_x(NCOMP)
         real(kind=RP),       intent(inout) :: U_y(NCOMP)
         real(kind=RP),       intent(inout) :: U_z(NCOMP)
      end subroutine OutflowBC_ChemPotNeumann
#endif
end module OutflowBCClass
