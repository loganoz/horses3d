module SCsensorClass

   use SMConstants,  only: RP, PI
   use ElementClass, only: Element
   use HexMeshClass, only: HexMesh

   use ShockCapturingKeywords

   implicit none

   private

   public :: SCsensor
   public :: Set_SCsensor

   type :: SCsensor_t

         real(RP) :: s0    !< Centre of the sensor scaling
         real(RP) :: ds    !< Bandwith of the sensor scaling
         real(RP) :: ds2   !< Half-bandwidth
         real(RP) :: low   !< Lower threshold
         real(RP) :: high  !< Upper threshold
         real(RP) :: s1    !< Threshold one (from 0 to 1)
         real(RP) :: s2    !< Threshold two (from 0 to 1), higher that s1
         integer  :: sVar  !< Variable used as input for the sensor

         procedure(Compute_Int), pointer :: Compute => null()
         procedure(Rescale_Int), pointer :: Rescale => SinRamp

   end type SCsensor_t
!
!  Interfaces
!  ----------
   abstract interface
      pure function Compute_Int(this, mesh, elem) result(val)
         import RP, SCsensor_t, HexMesh, Element
         class(SCsensor_t), intent(in) :: this
         type(HexMesh),     intent(in) :: mesh
         type(Element),     intent(in) :: elem
         real(RP)                      :: val
      end function Compute_Int

      pure function Rescale_Int(this, sensVal) result(scaled)
         import RP, SCsensor_t
         class(SCsensor_t), intent(in) :: this
         real(RP),          intent(in) :: sensVal
         real(RP)                      :: scaled
      end function Rescale_Int
   end interface

   type(SCsensor_t), protected :: SCsensor
!
!  ========
   contains
!  ========
!
   subroutine Set_SCsensor(sensor, sensorType, variable,        &
                           sensorLow, sensorHigh, thres1, thres2)
      implicit none
      type(SCsensor_t), intent(inout) :: sensor
      integer,          intent(in)    :: sensorType
      integer,          intent(in)    :: variable
      real(RP),         intent(in)    :: sensorLow
      real(RP),         intent(in)    :: sensorHigh
      real(RP),         intent(in)    :: thres1
      real(RP),         intent(in)    :: thres2

!
!     Sensor type
!     -----------
      select case (sensorType)
      case (SC_MODAL)
         sensor % Compute => Sensor_modal

      end select
!
!     Sensor parameters
!     -----------------
      sensor % sVar = variable
      sensor % low  = sensorLow
      sensor % high = sensorHigh
      sensor % s1   = thres1
      sensor % s2   = thres2
      sensor % s0   = (sensorHigh + sensorLow) / 2.0_RP
      sensor % ds   = (sensorHigh - sensorLow)
      sensor % ds2  = sensor % ds / 2.0_RP

   end subroutine Set_SCsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Sensors
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function Sensor_modal(sensor, mesh, elem) result(val)
      use NodalStorageClass, only: NodalStorage
      implicit none
      class(SCsensor_t), intent(in) :: sensor
      type(HexMesh),     intent(in) :: mesh
      type(Element),     intent(in) :: elem
      real(RP)                      :: val
!
!     Local variables
!     ---------------
      integer               :: i, j, k, r
      real(RP), allocatable :: sVar(:,:,:)
      real(RP), allocatable :: sVarMod(:,:,:)


      associate(Nx => elem % Nxyz(1), Ny => elem % Nxyz(2), Nz => elem % Nxyz(3))
!
!     Compute the sensed variable
!     ---------------------------
      allocate(sVar(0:Nx, 0:Ny, 0:Nz))
      do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
         sVar(i,j,k) = GetSensedVariable(sensor % sVar, elem % storage % Q(:,i,j,k))
      end do       ; end do       ; end do
!
!     Switch to modal space
!     ---------------------
      allocate(sVarMod(0:Nx, 0:Ny, 0:Nz), source=0.0_RP)

      ! Xi direction
      do k = 0, Nz ; do j = 0, Ny ; do r = 0, Nx ; do i = 0, Nx
         sVarMod(i,j,k) = sVarMod(i,j,k) + NodalStorage(Nx) % Fwd(i,r) * sVar(r,j,k)
      end do       ; end do       ; end do       ; end do

      ! Eta direction
      do k = 0, Nz ; do r = 0, Ny ; do j = 0, Ny ; do i = 0, Nx
         sVarMod(i,j,k) = sVarMod(i,j,k) + NodalStorage(Ny) % Fwd(j,r) * sVar(i,r,k)
      end do       ; end do       ; end do       ; end do

      ! Zeta direction
      do r = 0, Nz ; do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
         sVarMod(i,j,k) = sVarMod(i,j,k) + NodalStorage(Nz) % Fwd(k,r) * sVar(i,j,r)
      end do       ; end do       ; end do       ; end do
!
!     Ratio of higher modes vs all the modes
!     --------------------------------------
      val = log10( abs(sVarMod(Nx,Ny,Nz)) / norm2(sVarMod) )

      end associate

   end function Sensor_modal
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Scaling functions
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function SinRamp(sensor, sensedValue) result(scaled)
      implicit none
      class(SCsensor_t), intent(in) :: sensor
      real(RP),          intent(in) :: sensedValue
      real(RP)                      :: scaled
!
!     Local variables
!     ---------------
      real(RP) :: dev

      ! Deviation of the input value wrt the centre of the sensor
      dev = sensedValue - sensor % s0

      ! Sin-ramp activation function
      if (dev <= -sensor % ds2) then
         scaled = 0.0_RP

      else if (dev >= sensor % ds2) then
         scaled = 1.0_RP

      else
         scaled = ( 1.0_RP + sin(PI * dev / sensor % ds) ) / 2.0_RP

      end if

   end function SinRamp
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Utilities
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function GetSensedVariable(varType, Q) result(s)
!
!     Modules
!     -------
      use PhysicsStorage_NS,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
      use VariableConversion_NS, only: Pressure
      implicit none
!
!     Interface
!     ---------
      integer,  intent(in) :: varType
      real(RP), intent(in) :: Q(:)
      real(RP)             :: s


      select case (varType)
      case (SC_RHO)
         s = Q(IRHO)

      case (SC_RHOU)
         s = Q(IRHOU)

      case (SC_RHOV)
         s = Q(IRHOV)

      case (SC_RHOW)
         s = Q(IRHOW)

      case (SC_RHOE)
         s = Q(IRHOE)

      case (SC_U)
         s = Q(IRHOU) / Q(IRHO)

      case (SC_V)
         s = Q(IRHOV) / Q(IRHO)

      case (SC_W)
         s = Q(IRHOW) / Q(IRHO)

      case (SC_P)
         s = Pressure(Q)

      case (SC_RHOP)
         s = Pressure(Q) * Q(IRHO)

      end select

   end function GetSensedVariable

end module SCsensorClass
