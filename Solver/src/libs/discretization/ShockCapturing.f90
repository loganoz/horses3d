!
!//////////////////////////////////////////////////////
!
!   @File:    ShockCapturing.f90
!   @Author:  Andrés Mateo (andres.mgabin@upm.es)
!   @Created: Thu Jun 17 2021
!   @Last revision date: Thu Jun 17 2021
!   @Last revision author: Andrés Mateo (andres.mgabin@upm.es)
!
!//////////////////////////////////////////////////////
!
module ShockCapturing

   use SMConstants,   only: RP, STD_OUT, LINE_LENGTH
   use Utilities,     only: toLower
   use HexMeshClass,  only: HexMesh
   use ElementClass,  only: Element

   use SCsensorClass, only: SCsensor_t, SetSCsensor, DestructSCsensor
   use ShockCapturingKeywords

   implicit none

   public :: Initialize_ShockCapturing
   public :: Destruct_ShockCapturing
   public :: ShockCapturingDriver

   type SCdriver_t

         logical                       :: isActive          = .false.
         logical                       :: hasHyperbolicTerm = .false.
         character(len=:), allocatable :: method

         type(SCsensor_t),     private :: sensor

         procedure(Apply_Int), pointer :: Apply => null()

      contains

         procedure :: Detect => SC_detect

         final :: Destruct_ShockCapturing

   end type SCdriver_t

   type(SCdriver_t) :: ShockCapturingDriver
!
!  Interfaces
!  ----------
   abstract interface
      pure subroutine Apply_Int(self, mesh, elem, SVVflux)
         import RP, SCdriver_t, HexMesh, Element
         class(SCdriver_t), intent(in)    :: self
         type(HexMesh),     intent(inout) :: mesh
         type(Element),     intent(inout) :: elem
         real(RP),          intent(out)   :: SVVflux(:,:,:,:,:)
      end subroutine Apply_Int
   end interface
!
!  ========
   contains
!  ========
!
   subroutine Initialize_ShockCapturing(self, controlVariables)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t),        intent(inout) :: self
      class(FTValueDictionary), intent(in)    :: controlVariables
!
!     ---------------
!     Local variables
!     ---------------
      real(RP)                      :: lowThr
      real(RP)                      :: highThr
      real(RP)                      :: thr1
      real(RP)                      :: thr2
      character(len=:), allocatable :: sType
      character(len=:), allocatable :: sVar

!
!     Check if shock-capturing is requested
!     -------------------------------------
      if (controlVariables % containsKey(SC_KEY)) then
         self % isActive = controlVariables%logicalValueForKey(SC_KEY)
      else
         self % isActive = .false.
      end if

      if (.not. self % isActive) return
!
!     Sensor thresholds
!     --------------
      if (controlVariables % containsKey(SC_LOW_THRES_KEY)) then
         lowThr = controlVariables % doublePrecisionValueForKey(SC_LOW_THRES_KEY)
      else
         write(STD_OUT,*) 'ERROR. Lower threshold of the sensor must be specified.'
         stop
      end if

      if (controlVariables % containsKey(SC_HIGH_THRES_KEY)) then
         highThr = controlVariables % doublePrecisionValueForKey(SC_HIGH_THRES_KEY)
      else
         write(STD_OUT,*) 'ERROR. Higher threshold of the sensor must be specified.'
         stop
      end if

      if (controlVariables % containsKey(SC_THRES_1_KEY)) then
         thr1 = controlVariables % doublePrecisionValueForKey(SC_THRES_1_KEY)
      else
         thr1 = (highThr-lowThr) / 3.0_RP
      end if

      if (controlVariables % containsKey(SC_THRES_2_KEY)) then
         thr2 = controlVariables % doublePrecisionValueForKey(SC_THRES_2_KEY)
      else
         thr2 = (highThr-lowThr) * 2.0_RP / 3.0_RP
      end if
!
!     Sensor type
!     -----------
      if (controlVariables % containsKey(SC_SENSOR_KEY)) then
         sType = controlVariables % stringValueForKey(SC_SENSOR_KEY, LINE_LENGTH)
      else
         sType = SC_MODAL_KEY
      end if
      call toLower(sType)
!
!     Sensed variable
!     ---------------
      if (controlVariables % containsKey(SC_VARIABLE_KEY)) then
         sVar = controlVariables % stringValueForKey(SC_VARIABLE_KEY, LINE_LENGTH)
      else
         sVar = SC_RHOP_KEY
      end if
      call toLower(sVar)
!
!     Shock-capturing method
!     ----------------------
      if (controlVariables % containsKey(SC_METHOD_KEY)) then
         self % method = controlVariables % stringValueForKey(SC_METHOD_KEY, LINE_LENGTH)
      else
         self % method = SC_SSFV_KEY
      end if
      call toLower(self % method)

      select case (trim(self % method))
      case (SC_NOSVV_KEY)
         !self % Apply => SC_non_filtered
         self % hasHyperbolicTerm = .false.

      case (SC_SVV_KEY)
         !self % Apply => SC_SVV_filtered
         self % hasHyperbolicTerm = .false.

      case (SC_SSFV_KEY)
         !self % Apply => SC_entropy_stable_FV
         self % hasHyperbolicTerm = .true.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable Shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_SSFV_KEY
         stop

      end select
!
!     Construct sensor
!     ----------------
      call SetSCsensor(self % sensor, sType, sVar, lowThr, highThr, thr1, thr2)

   end subroutine Initialize_ShockCapturing
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine Destruct_ShockCapturing(self)
!
!     ---------
!     Interface
!     ---------
      type(SCdriver_t), intent(inout) :: self


      self%isActive = .false.
      self%hasHyperbolicTerm = .false.

      if (allocated(self % method)) deallocate(self % method)
      if (associated(self % Apply)) nullify(self % Apply)

      call DestructSCsensor(self % sensor)

   end subroutine Destruct_ShockCapturing
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine SC_detect(self, mesh, elem)
!
!     ---------
!     Interface
!     ---------
      class(SCdriver_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: elem

   end subroutine SC_detect

end module ShockCapturing
