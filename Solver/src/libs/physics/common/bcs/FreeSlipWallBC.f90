#include "Includes.h"
module FreeSlipWallBCClass
   use SMConstants
   use PhysicsStorage
   use VariableConversion
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine, CheckIfBoundaryNameIsContained
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   use Utilities, only: toLower, almostEqual
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
#ifdef NAVIERSTOKES
      logical           :: isAdiabatic
      real(kind=RP)     :: Twall
      real(kind=RP)     :: ewall       ! Wall internal energy
      real(kind=RP)     :: invTwall
      real(kind=RP)     :: wallType    ! 0/Adia 1/Isothermal
#endif
#ifdef CAHNHILLIARD
      real(kind=RP)     :: thetaw
#endif
      contains
         procedure         :: Destruct          => FreeSlipWallBC_Destruct
         procedure         :: Describe          => FreeSlipWallBC_Describe
#ifdef FLOW
         procedure         :: FlowState         => FreeSlipWallBC_FlowState
         procedure         :: CreateDeviceData  => FreeSlipWallBC_CreateDeviceData
         procedure         :: ExitDeviceData    => FreeSlipWallBC_ExitDeviceData
         procedure         :: FlowGradVars      => FreeSlipWallBC_FlowGradVars
         procedure         :: FlowNeumann       => FreeSlipWallBC_FlowNeumann
#endif
#ifdef CAHNHILLIARD
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

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructFreeSlipWallBC % BCType = "freeslipwall"
         ConstructFreeSlipWallBC % bname  = bname
         call toLower(ConstructFreeSlipWallBC % bname)

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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructFreeSlipWallBC % bname)) 
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
#ifdef NAVIERSTOKES
         if ( bcdict % ContainsKey("wall type (adiabatic/isothermal)") ) then
            keyval = bcdict % StringValueForKey("wall type (adiabatic/isothermal)", LINE_LENGTH)
            call tolower(keyval)
         else
            keyval = "adiabatic"
         end if

         if ( trim(keyval) .eq. "adiabatic" ) then
            ConstructFreeSlipWallBC % isAdiabatic = .true.
            ConstructFreeSlipWallBC % ewall = 0.0_RP
            ConstructFreeSlipWallBC % Twall = 0.0_RP
            ConstructFreeSlipWallBC % invTwall = 0.0_RP
            ConstructFreeSlipWallBC % wallType = 0.0_RP
         else
            ConstructFreeSlipWallBC % isAdiabatic = .false.
            call GetValueWithDefault(bcdict, "wall temperature" , refValues % T, ConstructFreeSlipWallBC % Twall     )
            ConstructFreeSlipWallBC % ewall = ConstructFreeSlipWallBC % Twall / (refValues % T*thermodynamics % gammaMinus1*dimensionless % gammaM2)
            ConstructFreeSlipWallBC % invTwall = dimensionless % gammaM2*refValues % T / ConstructFreeSlipWallBC % Twall
            ConstructFreeSlipWallBC % wallType = 1.0_RP
         end if
#endif
#ifdef CAHNHILLIARD
         call GetValueWithDefault(bcdict, "contact angle", 90.0_RP, ConstructFreeSlipWallBC % thetaw)
#endif

         close(fid)
         call bcdict % Destruct

      end function ConstructFreeSlipWallBC

      subroutine FreeSlipWallBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "FreeSlipWall"
#ifdef NAVIERSTOKES
         if ( self % isAdiabatic ) then
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Adiabatic"


         else
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Isothermal"
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall temperature: ', self % Twall * refValues % T
         end if
#endif
#ifdef CAHNHILLIARD
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall contact angle: ', self % thetaw
#endif
         
      end subroutine FreeSlipWallBC_Describe

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
#ifdef NAVIERSTOKES

      subroutine FreeSlipWallBC_CreateDeviceData(self)
         implicit none 
         class(FreeSlipWallBC_t), intent(in)    :: self

         !$acc enter data copyin(self)
      end subroutine FreeSlipWallBC_CreateDeviceData

      subroutine FreeSlipWallBC_ExitDeviceData(self)
         implicit none 
         class(FreeSlipWallBC_t), intent(in)    :: self

         !$acc exit data delete(self)

      end subroutine FreeSlipWallBC_ExitDeviceData

      subroutine FreeSlipWallBC_FlowState(self, mesh, zone)
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
         use HexMeshClass
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!         integer,                 intent(in)    :: zoneID                              
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: qNorm, pressure_aux
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j,zonefID,fID
         
         !$acc parallel loop gang present(mesh, self, zone) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)

               qNorm = mesh % faces(fID) % geom % normal(IX,i,j) * Q(IRHOU) + &
                       mesh % faces(fID) % geom % normal(IY,i,j) * Q(IRHOV) + &
                       mesh % faces(fID) % geom % normal(IZ,i,j) * Q(IRHOW) 
         
               Q(IRHOU:IRHOW) = Q(IRHOU:IRHOW) - 2.0_RP * qNorm * mesh % faces(fID) % geom % normal(:,i,j)
               
               mesh % faces(fID) % storage(2) % Q(IRHO:IRHOW,i,j) = Q(IRHO:IRHOW)
         
               !Isothermal BC
               pressure_aux = Q(IRHO) * self % Twall / (refValues % T * dimensionless % gammaM2)
               mesh % faces(fID) % storage(2) % Q(IRHOE,i,j) = Q(IRHOE) + self % wallType*(pressure_aux/thermodynamics % gammaMinus1 + & 
                                                           0.5_RP*(POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW)))/Q(IRHO) - Q(IRHOE))
               
            enddo ; enddo
         enddo
         !$acc end parallel loop

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowGradVars(self, mesh, zone)
!
!        *****************************************************************
!           Only set the temperature, velocity is Neumann, use interior!
!        *****************************************************************
!
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone

!         integer,                 intent(in)    :: zoneID 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhou_n
         real(kind=RP)  :: Q_aux(NCONS),Q(NCONS)
         real(kind=RP)  :: u_int(NGRAD), u_star(NGRAD)
         integer        :: i,j,zonefID,fID
         
         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, Q_aux, u_star, u_int)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)

               call NSGradientVariables_STATE(NCONS, NGRAD, Q, u_int)

               Q_aux(IRHO) = Q(IRHO)
               Q_aux(IRHOU:IRHOW) = Q(IRHOU:IRHOW)
               Q_aux(IRHOE) = Q(IRHOE) + self % wallType*(Q(IRHO)*self % eWall+0.5_RP*(POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW)))/Q(IRHO)-Q(IRHOE))
#if defined(SPALARTALMARAS)
               Q_aux(IRHOTHETA)= Q(IRHOTHETA)
#endif
               call NSGradientVariables_STATE(NCONS, NGRAD, Q_aux, u_star)

               mesh % faces(fID) % storage(1) % unStar(:,1,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(1,i,j) * mesh % faces(fID) % geom % jacobian(i,j)
               mesh % faces(fID) % storage(1) % unStar(:,2,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(2,i,j) * mesh % faces(fID) % geom % jacobian(i,j)    
               mesh % faces(fID) % storage(1) % unStar(:,3,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(3,i,j) * mesh % faces(fID) % geom % jacobian(i,j)

            enddo ; enddo
         enddo
         !$acc end parallel loop  
      end subroutine FreeSlipWallBC_FlowGradVars

      subroutine FreeSlipWallBC_FlowNeumann(self, mesh, zone)
!
!        ***********************************************************
!           In momentum, free slip is Neumann. In temperature, 
!           depends on the adiabatic/isothermal choice
!        ***********************************************************
!
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         type(Zone_t), intent(in)               :: zone
!
!        ---------------
!        Local Variables
!        ---------------
!   
         integer        :: i,j,zonefID,fID
         real(kind=RP)  :: viscWork, heatFlux
         real(kind=RP)  :: flux(NCONS),Q(NCONS)

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, flux, viscWork, heatFlux)     
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)

               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)
               flux = mesh % faces(fID) % storage(2) % FStar(:,i,j)
               
               viscWork = (flux(IRHOU)*Q(IRHOU)+flux(IRHOV)*Q(IRHOV)+flux(IRHOW)*Q(IRHOW))/Q(IRHO)
               heatFlux = flux(IRHOE) - viscWork
               flux(IRHO:IRHOW) = 0.0_RP
               flux(IRHOE) = self % wallType * heatFlux  ! 0 (Adiabatic)/ heatFlux (Isothermal)
#if defined(SPALARTALMARAS)
               flux(IRHOTHETA) = 0.0_RP
#endif
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = flux(:)
            enddo ; enddo
         enddo
         !$acc end parallel loop  

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
#ifdef INCNS
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
         real(kind=RP),       intent(inout) :: Q(NCONS)
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

         Q(INSRHO)          = Q(INSRHO)
         Q(INSRHOU:INSRHOW) = Q(INSRHOU:INSRHOW) - 2.0_RP * vn * nHat
         Q(INSP)            = Q(INSP)

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
!        **************************************************************
!           Use the interior velocity: Neumann BC!
!        **************************************************************
!
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients

      end subroutine FreeSlipWallBC_FlowGradVars

      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        ***************************************************************
!           Set homogeneous Neumann BCs everywhere
!        ***************************************************************
!        
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)      :: x(NDIM)
         real(kind=RP),       intent(in)      :: t
         real(kind=RP),       intent(in)      :: nHat(NDIM)
         real(kind=RP),       intent(in)      :: Q(NCONS)
         real(kind=RP),       intent(in)      :: U_x(NCONS)
         real(kind=RP),       intent(in)      :: U_y(NCONS)
         real(kind=RP),       intent(in)      :: U_z(NCONS)
         real(kind=RP),       intent(inout)   :: flux(NCONS)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for multiphase solver
!        ---------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef MULTIPHASE

      subroutine FreeSlipWallBC_CreateDeviceData(self)
         implicit none 
         class(FreeSlipWallBC_t), intent(in)    :: self

         !$acc enter data copyin(self)
      end subroutine FreeSlipWallBC_CreateDeviceData

      subroutine FreeSlipWallBC_ExitDeviceData(self)
         implicit none 
         class(FreeSlipWallBC_t), intent(in)    :: self

         !$acc exit data delete(self)

      end subroutine FreeSlipWallBC_ExitDeviceData

      subroutine FreeSlipWallBC_FlowState(self, mesh, zoneID)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!

         use HexMeshClass
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)              :: mesh
         integer,                 intent(in)    :: zoneID 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: vn
         real(kind=RP) :: Q(NCONS)
         integer       :: i,j,zonefID,fID
         !!$acc parallel loop gang present(mesh, self, zone) async(1)
         !$acc parallel loop gang present(mesh, self, zone)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID = mesh % zones(zoneID) % faces(zonefID)
            !$acc loop vector collapse(2) private(Q)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)
               
!
!              -----------------------------------------------
!              Generate the external flow along the face, that
!              represents a solid wall.
!              -----------------------------------------------
!
               vn =  mesh % faces(fID) % geom % normal(IX,i,j) * Q(IMSQRHOU) + &
                     mesh % faces(fID) % geom % normal(IY,i,j) * Q(IMSQRHOV) + &
                     mesh % faces(fID) % geom % normal(IZ,i,j) * Q(IMSQRHOW) 

               Q(IMC)          = Q(IMC)
               Q(IMSQRHOU:IMSQRHOW) = Q(IMSQRHOU:IMSQRHOW) - 2.0_RP * vn * mesh % faces(fID) % geom % normal(:,i,j)
               Q(IMP)            = Q(IMP)

               mesh % faces(fID) % storage(2) % Q(:,i,j) = Q
            enddo ; enddo
         enddo
         !$acc end parallel loop

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowGradVars(self, mesh, zoneID)
         !TODO slightly unsure about this one 
!
!        **************************************************************
!           Use the interior velocity: Neumann BC!
!        **************************************************************
!
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)           :: mesh
         integer,                 intent(in)    :: zoneID 

!
!        ---------------
!        Local variables
!        ---------------
!        
         real(kind=RP)  :: Q(NCONS)
         integer        :: i,j,zonefID,fID
         real(kind=RP)  :: u_int(NGRAD), u_star(NGRAD)

         !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1) 
         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID = mesh % zones(zoneID) % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, u_int)            
            do j = 0, mesh % faces(fID) % Nf(2)  ; do i = 0, mesh % faces(fID) % Nf(1)
               
               Q = mesh % faces(fID) % storage(1) % Q(:,i,j)

               rho = dimensionless % rho(1)*Q(IMC) + dimensionless % rho(2)*(1.0_RP - Q(IMC))
               call mGradientVariables(NCONS, NGRAD, Q, u_int,rho )
               u_star = u_int

               mesh % faces(fID) % storage(1) % unStar(:,1,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(1,i,j) * mesh % faces(fID) % geom % jacobian(i,j)
               mesh % faces(fID) % storage(1) % unStar(:,2,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(2,i,j) * mesh % faces(fID) % geom % jacobian(i,j)    
               mesh % faces(fID) % storage(1) % unStar(:,3,i,j) = (u_star-u_int) * mesh % faces(fID) % geom % normal(3,i,j) * mesh % faces(fID) % geom % jacobian(i,j)
               
            enddo ; enddo
         enddo
         !$acc end parallel loop
         
      end subroutine FreeSlipWallBC_FlowGradVars

      subroutine FreeSlipWallBC_FlowNeumann(self, mesh, zoneID)
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         type(HexMesh), intent(inout)              :: mesh
         integer,                 intent(in)    :: zoneID 

         integer        :: i,j,zonefID,fID
         real(kind=RP)  :: flux(NCONS),Q(NCONS)

         !!$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         !$acc parallel loop gang present(mesh, self, zone) private(fID)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID =  mesh % zones(zoneID) % faces(zonefID)
            !$acc loop vector collapse(2) independent private(Q,flux)  
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = 0.0_RP
            enddo 
          enddo
         enddo
         !$acc end parallel loop

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

      subroutine FreeSlipWallBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)
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

         flux = -4.0_RP * multiphase % invEps * cos(DEG2RAD*self % thetaw) * prod 

      end subroutine FreeSlipWallBC_PhaseFieldNeumann

      subroutine FreeSlipWallBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine FreeSlipWallBC_ChemPotState

      subroutine FreeSlipWallBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_ChemPotNeumann

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
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),       intent(in)       :: x(NDIM)
         real(kind=RP),       intent(in)       :: t
         real(kind=RP),       intent(in)       :: nHat(NDIM)
         real(kind=RP),       intent(inout)    :: Q(NCONS)
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
         vn = sum(Q(ICAAU:ICAAW)*nHat)

         Q(ICAARHO) = Q(ICAARHO)
         Q(ICAAU:ICAAW) = Q(ICAAU:ICAAW) - 2.0_RP * vn * nHat
         Q(ICAAP) = Q(ICAAP)

      end subroutine FreeSlipWallBC_FlowState
!
!     not use for acoustic
      subroutine FreeSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
      end subroutine FreeSlipWallBC_FlowGradVars
!
      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),            intent(in)    :: x(NDIM)
         real(kind=RP),            intent(in)    :: t
         real(kind=RP),            intent(in)    :: nHat(NDIM)
         real(kind=RP),            intent(in)    :: Q(NCONS)
         real(kind=RP),            intent(in)    :: U_x(NCONS)
         real(kind=RP),            intent(in)    :: U_y(NCONS)
         real(kind=RP),            intent(in)    :: U_z(NCONS)
         real(kind=RP),            intent(inout) :: flux(NCONS)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_FlowNeumann
#endif
!
end module FreeSlipWallBCClass
