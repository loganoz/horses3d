#include "Includes.h"
module InflowBCClass
   use SMConstants
   use PhysicsStorage
   use VariableConversion, only: GetGradientValues_f
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine, CheckIfBoundaryNameIsContained
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities, only: toLower, almostEqual
   use GenericBoundaryConditionClass
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
      real(kind=RP)              :: TurbIntensity
      real(kind=RP)              :: eddy_theta
#endif
#if defined(INCNS)
      real(kind=RP)              :: AoAPhi
      real(kind=RP)              :: AoATheta
      real(kind=RP)              :: v
      real(kind=RP)              :: rho
#endif
#if defined(MULTIPHASE)
      real(kind=RP)              :: AoAPhi
      real(kind=RP)              :: AoATheta
      logical                    :: isLayered  = .false.
      logical                    :: isXLimited = .false.
      logical                    :: isYLimited = .false.
      logical                    :: isZLimited = .false.
      real(kind=RP)              :: xLim, yLim, zLim
      real(kind=RP)              :: phase1Vel
      real(kind=RP)              :: phase2Vel
      real(kind=RP)              :: c
      real(kind=RP)              :: v
#endif
#if defined(ACOUSTIC)
      real(kind=RP)              :: v
      real(kind=RP)              :: rho
      real(kind=RP)              :: p
#endif
      contains
         procedure         :: Destruct          => InflowBC_Destruct
         procedure         :: Describe          => InflowBC_Describe
#if defined(FLOW)
         procedure         :: FlowState         => InflowBC_FlowState
         procedure         :: FlowNeumann       => InflowBC_FlowNeumann
         procedure         :: CreateDeviceData  => InflowBC_CreateDeviceData
         procedure         :: ExitDeviceData    => InflowBC_ExitDeviceData
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
!                 TurbIntensity    = #value        (only in compressible NS)
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
         call GetValueWithDefault(bcdict, "TurbIntensity", 0.0_RP                            , ConstructInflowBC % TurbIntensity)
#if defined(SPALARTALMARAS)
         call GetValueWithDefault(bcdict, "Turbulence parameter theta", refValues % mu    , ConstructInflowBC % eddy_theta)
#endif
         ConstructInflowBC % p        = ConstructInflowBC % p / refValues % p
         ConstructInflowBC % rho      = ConstructInflowBC % rho / refValues % rho
         ConstructInflowBC % v        = ConstructInflowBC % v * sqrt(thermodynamics % gamma * ConstructInflowBC % p / ConstructInflowBC % rho)
         ConstructInflowBC % AoATheta = ConstructInflowBC % AoATheta * PI / 180.0_RP
         ConstructInflowBC % AoAPhi   = ConstructInflowBC % AoAPhi * PI / 180.0_RP
         ConstructInflowBC % TurbIntensity   = ConstructInflowBC % TurbIntensity / 100.0_RP
#if defined(SPALARTALMARAS)
         ConstructInflowBC % eddy_theta = ConstructInflowBC % eddy_theta / refValues % mu
#endif

#elif defined(INCNS)
         call GetValueWithDefault(bcdict, "velocity", refValues % v, ConstructInflowBC % v)
         call GetValueWithDefault(bcdict, "aoaphi"  , refValues % AoAPhi                     , ConstructInflowBC % AoAPhi  )
         call GetValueWithDefault(bcdict, "aoatheta", refValues % AoATheta                   , ConstructInflowBC % AoATheta)
         call GetValueWithDefault ( bcdict , "density" , refValues % rho , ConstructInflowBC % rho ) 

         ConstructInflowBC % v = ConstructInflowBC % v / refValues % v
         ConstructInflowBC % AoATheta = ConstructInflowBC % AoATheta * PI / 180.0_RP
         ConstructInflowBC % AoAPhi   = ConstructInflowBC % AoAPhi * PI / 180.0_RP
         ConstructInflowBC % rho = ConstructInflowBC % rho / refValues % rho
#endif
#if defined(MULTIPHASE)
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

         else
             ! standard one phase bc
             call GetValueWithDefault(bcdict, "velocity", refValues % v, ConstructInflowBC % v)
             call GetValueWithDefault(bcdict, "aoaphi"  , refValues % AoAPhi                     , ConstructInflowBC % AoAPhi  )
             call GetValueWithDefault(bcdict, "aoatheta", refValues % AoATheta                   , ConstructInflowBC % AoATheta)
             call GetValueWithDefault ( bcdict , "concentration" , 0.0_RP , ConstructInflowBC % c ) 

             ConstructInflowBC % v = ConstructInflowBC % v / refValues % v
             ConstructInflowBC % AoATheta = ConstructInflowBC % AoATheta * PI / 180.0_RP
             ConstructInflowBC % AoAPhi   = ConstructInflowBC % AoAPhi * PI / 180.0_RP
         end if
#endif

#if defined(ACOUSTIC)
         ! by default Acoustic perturbations are assumed to be 0, is a far field BC rather than a real inflow
         call GetValueWithDefault(bcdict, "acoustic pressure", 0.0_RP, ConstructInflowBC % p   )
         call GetValueWithDefault(bcdict, "acoustic density" , 0.0_RP, ConstructInflowBC % rho )
         call GetValueWithDefault(bcdict, "acoustic velocity", 0.0_RP, ConstructInflowBC % v   )
         ConstructInflowBC % p     = ConstructInflowBC % p / refValues % p
         ConstructInflowBC % rho   = ConstructInflowBC % rho / refValues % rho
         ConstructInflowBC % v     = ConstructInflowBC % v / refValues % v
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

         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Max. Vel. Fluct. in % (from TurbIntensity): ', (self % TurbIntensity)
#if defined(SPALARTALMARAS)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Initial Value of Turbulence variable: ', (3.0_RP * self % eddy_theta)
#endif
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
         else
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Concentration: ", self % c
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", "All Phases velocity: ", self % V * refValues % v
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' AoaPhi: ', self % AoAPhi * 180.0_RP / PI
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' AoaTheta: ', self % AoATheta * 180.0_RP / PI
         end if
#endif
#elif defined(ACOUSTIC)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Velocity: ', self % v * refValues % v
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Pressure: ', self % p * refValues % p
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Density: ', self % rho * refValues % rho
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
!        Subroutines for all flow equations
!        -----------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(FLOW)

      subroutine InflowBC_CreateDeviceData(self)
         implicit none 
         class(InflowBC_t), intent(in)    :: self

         !$acc enter data copyin(self)

      end subroutine InflowBC_CreateDeviceData

      subroutine InflowBC_ExitDeviceData(self)
         implicit none 
         class(InflowBC_t), intent(in)    :: self

         !$acc exit data delete(self)

      end subroutine InflowBC_ExitDeviceData
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)

      subroutine InflowBC_FlowState(self, mesh, zone)
         implicit none
         class(InflowBC_t),   intent(in)    :: self
         type(HexMesh),       intent(inout)    :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: qq, u, v, w
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j
         integer       :: fID
         integer       :: zonefID
      
         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) independent private(Q)  
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)

               qq = self % v
               u  = qq*cos(self % AoAtheta)*COS(self % AoAphi)
               v  = qq*sin(self % AoAtheta)*COS(self % AoAphi)
               w  = qq*SIN(self % AoAphi)

               Q(1) = self % rho
               Q(2) = Q(1)*u
               Q(3) = Q(1)*v
               Q(4) = Q(1)*w
               Q(5) = self % p/(thermodynamics % gamma - 1._RP) + 0.5_RP*Q(1)*(u**2 + v**2 + w**2)
#if defined(SPALARTALMARAS)
               Q(6) = Q(1) * 3.0_RP*self % eddy_theta
#endif
               mesh % faces(fID) % storage(2) % Q(:,i,j) = Q 
            enddo 
          enddo
         enddo
         !$acc end parallel loop

      end subroutine InflowBC_FlowState

      subroutine InflowBC_FlowNeumann(self, mesh, zone)
         implicit none
         class(InflowBC_t),   intent(in)    :: self
         type(HexMesh),       intent(inout) :: mesh
         type(Zone_t),        intent(in)    :: zone

!         integer,             intent(in)    :: zoneID 
!
!        *******************************************
!        Cancel out the viscous flux at the inlet
!        *******************************************
!
         real(kind=RP)  :: flux(NCONS)
         integer       :: i,j
         integer       :: fID
         integer       :: zonefID
         
         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) independent private(flux)  
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = 0.0_RP
            enddo 
          enddo
         enddo
         !$acc end parallel loop

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
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rho, u, v, w, vel

         u = self % v * cos(self % AoAtheta) * cos(self % AoAphi)
         v = self % v * sin(self % AoAtheta) * cos(self % AoAphi)
         w = self % v * sin(self % AoAphi)

         Q(INSRHO)  = self % rho
         Q(INSRHOU) = Q(INSRHO)*u
         Q(INSRHOV) = Q(INSRHO)*v
         Q(INSRHOW) = Q(INSRHO)*w

      end subroutine InflowBC_FlowState

      subroutine InflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(InflowBC_t),  intent(in)     :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)

         flux = 0.0_RP

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

      subroutine InflowBC_PhaseFieldState(self, mesh, zone)
         use HexMeshClass
         implicit none
         class(InflowBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j,zonefID,fID
         
         !$acc parallel loop gang present(mesh, self, zone) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)

               mesh % faces(fID) % storage(2) % Q(:,i,j) = Q(:)
            
            enddo ; enddo
         enddo
         !$acc end parallel loop
      end subroutine InflowBC_PhaseFieldState

      subroutine InflowBC_PhaseFieldGradVars(self, mesh, zone)
         implicit none
         class(InflowBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!
!        ---------------
!        Local variables
!        ---------------
!        
         integer        :: i,j,zonefID,fID

         !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1) 
         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               mesh % faces(fID) % storage(1) % unStar(:,1,i,j) = 0.0
               mesh % faces(fID) % storage(1) % unStar(:,2,i,j) = 0.0    
               mesh % faces(fID) % storage(1) % unStar(:,3,i,j) = 0.0
               
            enddo ; enddo
         enddo
         !$acc end parallel loop
         
      end subroutine InflowBC_PhaseFieldGradVars
      
      subroutine InflowBC_PhaseFieldNeumann(self, mesh, zone)
         use HexMeshClass
         implicit none
         class(InflowBC_t),      intent(in)    :: self
         type(HexMesh),           intent(inout) :: mesh
         type(Zone_t), intent(in)               :: zone

         !local variables
         integer        :: i,j,zonefID,fID

         !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) independent 
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = 0.0_RP
            enddo ; enddo
         enddo
         !$acc end parallel loop

      end subroutine InflowBC_PhaseFieldNeumann

      subroutine InflowBC_ChemPotState(self, mesh, zone)
         use HexMeshClass
         implicit none
         class(InflowBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j,zonefID,fID
         
         !$acc parallel loop gang present(mesh, self, zone)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)

               mesh % faces(fID) % storage(2) % Q(:,i,j) = Q(:)
            
            enddo ; enddo
         enddo
         !$acc end parallel loop
      end subroutine InflowBC_ChemPotState

      subroutine InflowBC_ChemPotGradVars(self, mesh, zone)
         implicit none
         class(InflowBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!
!        ---------------
!        Local variables
!        ---------------
!        
         integer        :: i,j,zonefID,fID

         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
   
               mesh % faces(fID) % storage(1) % unStar(:,1,i,j) = 0.0_RP
               mesh % faces(fID) % storage(1) % unStar(:,2,i,j) = 0.0_RP
               mesh % faces(fID) % storage(1) % unStar(:,3,i,j) = 0.0_RP
               
            enddo ; enddo
         enddo
         !$acc end parallel loop
         
      end subroutine InflowBC_ChemPotGradVars

      subroutine InflowBC_ChemPotNeumann(self, mesh, zone)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
         type(HexMesh),           intent(inout) :: mesh
         type(Zone_t), intent(in)               :: zone

         !local variables
         integer        :: i,j,zonefID,fID

         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) independent  
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = 0.0_RP
            enddo ; enddo
         enddo
         !$acc end parallel loop

      end subroutine InflowBC_ChemPotNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Multiphase Solver only
!        --------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(MULTIPHASE)
      subroutine InflowBC_FlowState(self, mesh, zone)
         implicit none
         class(InflowBC_t),  intent(in)        :: self
         type(HexMesh),       intent(inout)    :: mesh
         type(Zone_t), intent(in)              :: zone
!
!        ---------------
!        Local variables
!        ---------------

         real(kind=RP) :: qq, u, v, w, c, rho, sqrtRho
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j
         integer       :: fID
         integer       :: zonefID
!
         if (self % isLayered) then

             ! flow always in the x direction unless is the interphase is normal to x, in that case is in z direction
            direction_cond:if (self % isXLimited) then
                !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
                !$acc parallel loop gang present(mesh, self, zone) private(fID)
                do zonefID = 1, zone % no_of_faces
                   fID = zone % faces(zonefID)
                   !$acc loop vector collapse(2) independent private(Q)  
                   do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                      c  = 1.0 - 0.5*(1.0+tanh(2.0*(( mesh % faces(fID) % geom % x(IX,i,j) - self % xLim))/multiphase % eps))
                      rho = dimensionless % rho(1)*c + dimensionless % rho(2)*(1.0_RP - c)
                      sqrtRho = sqrt(rho)
                      u  = 0.0_RP
                      v  = 0.0_RP
                      w = self % phase1Vel * c + self % phase2Vel * (1.0_RP - c)

                      Q(IMC) = c
                      Q(IMSQRHOU) = sqrtRho*u
                      Q(IMSQRHOV) = sqrtRho*v
                      Q(IMSQRHOW) = sqrtRho*w
                      Q(IMP) = mesh % faces(fID) % storage(1) % Q(IMP,i,j)

                      mesh % faces(fID) % storage(2) % Q(:,i,j) = Q 
                    enddo ; enddo
                 enddo
                !$acc end parallel loop
            else if(self % isYLimited) then
                !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
                !$acc parallel loop gang present(mesh, self, zone) private(fID)
                do zonefID = 1, zone % no_of_faces
                   fID = zone % faces(zonefID)
                   !$acc loop vector collapse(2) independent private(Q)  
                   do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                      c  = 1.0 - 0.5*(1.0+tanh(2.0*(( mesh % faces(fID) % geom % x(IY,i,j) - self % yLim))/multiphase % eps))
                      rho = dimensionless % rho(1)*c + dimensionless % rho(2)*(1.0_RP - c)
                      sqrtRho = sqrt(rho)
                      u = self % phase1Vel * c + self % phase2Vel * (1.0_RP - c)
                      v  = 0.0_RP
                      w  = 0.0_RP

                      Q(IMC) = c
                      Q(IMSQRHOU) = sqrtRho*u
                      Q(IMSQRHOV) = sqrtRho*v
                      Q(IMSQRHOW) = sqrtRho*w
                      Q(IMP) = mesh % faces(fID) % storage(1) % Q(IMP,i,j)

                      mesh % faces(fID) % storage(2) % Q(:,i,j) = Q 
                    enddo ; enddo
                enddo
                !$acc end parallel loop
            else if(self % isZLimited) then
                !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
                !$acc parallel loop gang present(mesh, self, zone) private(fID)
                do zonefID = 1, zone % no_of_faces
                   fID = zone % faces(zonefID)
                   !$acc loop vector collapse(2) independent private(Q)  
                   do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                      c  = 1.0 - 0.5*(1.0+tanh(2.0*(( mesh % faces(fID) % geom % x(IZ,i,j) - self % zLim))/multiphase % eps))
                      rho = dimensionless % rho(1)*c + dimensionless % rho(2)*(1.0_RP - c)
                      sqrtRho = sqrt(rho)
                      u = self % phase1Vel * c + self % phase2Vel * (1.0_RP - c)
                      v  = 0.0_RP
                      w  = 0.0_RP

                      Q(IMC) = c
                      Q(IMSQRHOU) = sqrtRho*u
                      Q(IMSQRHOV) = sqrtRho*v
                      Q(IMSQRHOW) = sqrtRho*w
                      Q(IMP) = mesh % faces(fID) % storage(1) % Q(IMP,i,j)

                      mesh % faces(fID) % storage(2) % Q(:,i,j) = Q 
                   enddo ; enddo
                enddo
                !$acc end parallel loop
            end if direction_cond

         else
            !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
            !$acc parallel loop gang present(mesh, self, zone) private(fID)
            do zonefID = 1, zone % no_of_faces
               fID = zone % faces(zonefID)
               !$acc loop vector collapse(2) independent private(Q)  
               do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)

                  qq = self % v
                  u  = qq*cos(self % AoAtheta)*COS(self % AoAphi)
                  v  = qq*sin(self % AoAtheta)*COS(self % AoAphi)
                  w  = qq*SIN(self % AoAphi)
                  rho = dimensionless % rho(1)*self % c + dimensionless % rho(2)*(1.0_RP - self % c)
                  sqrtRho = sqrt(rho)

                  Q(IMC) = self % c
                  Q(IMSQRHOU) = sqrtRho*u
                  Q(IMSQRHOV) = sqrtRho*v
                  Q(IMSQRHOW) = sqrtRho*w
                  Q(IMP) = mesh % faces(fID) % storage(1) % Q(IMP,i,j)

                  mesh % faces(fID) % storage(2) % Q(:,i,j) = Q 
               enddo ; enddo
            enddo
            !$acc end parallel loop
         end if 

      end subroutine InflowBC_FlowState

      subroutine InflowBC_FlowNeumann(self, mesh, zone)
         implicit none
         class(InflowBC_t),   intent(in)    :: self
         type(HexMesh),       intent(inout) :: mesh
         type(Zone_t),        intent(in)    :: zone

         integer       :: i,j
         integer       :: fID
         integer       :: zonefID
         
         ! flux = 0.0_RP directly stored in Fstar

         !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = 0.0_RP
            enddo 
          enddo
         enddo
         !$acc end parallel loop

      end subroutine InflowBC_FlowNeumann

#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Acoustic APE equations
!        ---------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(ACOUSTIC)
      subroutine InflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(InflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: u, v, w

         ! not use aoa here
         u = self % v
         v = self % v
         w = self % v

         Q(ICAARHO) = self % rho
         Q(ICAAU) = u
         Q(ICAAV) = v
         Q(ICAAW) = w
         Q(ICAAP) = self % p

      end subroutine InflowBC_FlowState

      subroutine InflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(InflowBC_t),  intent(in)     :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)

         flux = 0.0_RP

      end subroutine InflowBC_FlowNeumann
#endif
!
end module InflowBCClass
