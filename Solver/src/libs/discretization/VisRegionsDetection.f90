module VisRegionsDetection

    use SMConstants,                only: RP, STD_OUT, LINE_LENGTH, NDIM, IX, IY, IZ
    use PhysicsStorage,             only: NCONS, NGRAD
    use FluidData,                  only: dimensionless
    use Utilities,                  only: toLower
    use DGSEMClass,                 only: ComputeTimeDerivative_f, DGSem
    use HexMeshClass,               only: HexMesh
    use ElementClass,               only: Element
    use RegionsDetectionKeywords
    use VisSensorClass,             only: Set_VISsensor, VISsensor_t,Destruct_VISsensor
    implicit none

    public :: Initialize_ViscousRegionDetection
    public :: ViscousRegionDetectionDriver
    type VISdriver_t

      logical  :: isActive = .false.  !< On/Off flag
      integer  :: IterJump            !< Interval of iterations to activate the sensor each IterJump iteration
      integer  :: IterMin             !< Minimum number of iterations to start activating the sensor
      logical  :: toAdapt  = .false.  !< If true perform P-adaptation with the viscous sensor
      logical  :: toHybrid = .false.  !< if true set off viscous fluxes in the inviscid region
      type(VISsensor_t),private :: sensor   !< Sensor
      contains

      procedure :: Detect           => VIS_detect
      procedure :: Describe         => VIS_describe
      final :: VIS_destruct
 
    end type 
    type(VISdriver_t), allocatable :: ViscousRegionDetectionDriver 
    ! ========
      contains
   !  ========
   

subroutine Initialize_ViscousRegionDetection(self, controlVariables, sem, &
        TimeDerivative, TimeDerivativeIsolated)
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
type(VISdriver_t), allocatable, intent(inout) :: self
class(FTValueDictionary),      intent(in)    :: controlVariables
class(DGSem),                  intent(inout) :: sem
procedure(ComputeTimeDerivative_f)           :: TimeDerivative
procedure(ComputeTimeDerivative_f)           :: TimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
character(len=:), allocatable :: method
integer                       :: minSteps

!
!     Check if viscous regions detection is requested
!     -------------------------------------
allocate(self)
if (controlVariables % containsKey(VIS_KEY)) then
self % isActive = controlVariables % logicalValueForKey(VIS_KEY)
else
self % isActive = .false.
end if

if (.not. self % isActive) then
return
end if
!     Sensor 'inertia'
!     ----------------
if (controlVariables % containsKey(VIS_SENSOR_INERTIA_KEY)) then
    minSteps = controlVariables % realValueForKey(VIS_SENSOR_INERTIA_KEY)
    if (minSteps < 1) then
       write(STD_OUT,*) 'ERROR. Sensor inertia must be at least 1.'
       stop
    end if
 else
    minSteps = 1
 end if
!--------- Minimum numbero of iterations to activate the sensor -----------------------------------
if (controlVariables % containsKey( VIS_ITER_MIN)) then
   self % IterMin = controlVariables % integerValueForKey(VIS_ITER_MIN)
else
   self % IterMin = 1
end if 
!--------------------------------------------------------------------------------------------------
!--------- Interval of number of iteration to activate the sensor --------------------------------- 
if (controlVariables % containsKey(VIS_ITER_JUMP)) then 
   self % IterJump = controlVariables % integerValueForKey(VIS_ITER_JUMP)
else
   self % IterJump = 1
end if
!---------- Set off the viscous fluxes computations in the inviscid region if requested -----------
if (controlVariables % containsKey(VIS_TO_HYBRID)) then 
   self % toHybrid = controlVariables % logicalValueForKey(VIS_TO_HYBRID)
else
   self % toHybrid = .false.   
end if
!--------------------------- Perform p-adaptation with the viscous sensor if requested -------------
if (controlVariables % containsKey(VIS_TO_ADAPT)) then 
   self % toAdapt = controlVariables % logicalValueForKey(VIS_TO_ADAPT)
else
   self % toAdapt = .false.
end if
!
!     Construct sensor
!     ----------------
 call Set_VISsensor(self % sensor, controlVariables, sem, minSteps, &
                   TimeDerivative, TimeDerivativeIsolated)

end subroutine Initialize_ViscousRegionDetection
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Base class
!
!///////////////////////////////////////////////////////////////////////////////
!
subroutine VIS_detect(self, sem, t)
    !
    !     ---------
    !     Interface
    !     ---------
          implicit none
          class(VISdriver_t), intent(inout) :: self
          type(DGSem),       intent(inout) :: sem
          real(RP),          intent(in)    :: t
    
          
          call self % sensor % Compute(sem, t)
          
end subroutine VIS_detect
subroutine VIS_describe(self)
        !
        !     -------
        !     Modules
        !     -------
!use MPI_Process_Info, only: MPI_Process
        use Headers,          only: Subsection_Header
        !
        !     ---------
        !     Interface
        !     ---------
        implicit none
        class(VISdriver_t), intent(in) :: self
if (.not. self % isActive) return 
write(STD_OUT, "(/)")
      call Subsection_Header("Viscous regions detection")

      call self % sensor % Describe()
      write(STD_OUT,"(30X,A,A52,L)") "->", " Set off the viscous fluxes in the inviscid region: ", self % toHybrid
      write(STD_OUT,"(30X,A,A52,L)") "->", " Perform p-adaptation with the viscous sensor: ", self % toAdapt
      write(STD_OUT,"(30X,A,A52,I0,A)") "->", " minimum number of iterations to enable the sensor: ", self % IterMin, " iteration(s)"
      write(STD_OUT,"(30X,A,A52,I0,A)") "->", " interval of iterations to activate the sensor: ", self % IterJump, " iteration(s)"
end subroutine VIS_describe
subroutine VIS_destruct(self)
    !
    !     ---------
    !     Interface
    !     ---------
          implicit none
          type(VISdriver_t), intent(inout) :: self
    
    
          self%isActive = .false.
          
    
          call Destruct_VISsensor(self % sensor)
end subroutine VIS_destruct
end module VisRegionsDetection      


!      
                 
