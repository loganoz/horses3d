module VisSensorClass
   use SMConstants,       only: RP, PI, STD_OUT, NDIM, IX, IY, IZ, LINE_LENGTH
   use DGSEMClass,        only: ComputeTimeDerivative_f, DGSem
   use PhysicsStorage,    only: grad_vars, GRADVARS_STATE, GRADVARS_ENERGY, GRADVARS_ENTROPY, NCONS
   use NodalStorageClass, only: NodalStorage, NodalStorage_t
   use ElementClass,      only: Element
   use Utilities,         only: toLower
   use Clustering,        only: GMM_t,standardize
   use StopwatchClass                  , only: Stopwatch
   use RegionsDetectionKeywords

   implicit none

   private
   public :: VISsensor_t!, VIS_GMM
   public :: Set_VISsensor
   public :: Destruct_VISsensor

   type :: VISsensor_t
   integer  :: sens_type  ! ID of the sensor type
   integer  :: sVar       !< Variable used as input for the sensor
   integer  :: min_steps  !< Minimum number of steps before the sensor is deactivated
   integer  :: e_clust    !< method used to cluster elements to different regions after fitting GMM (either 0 for max method (default method) or 1 for mean method)

         type(GMM_t)           :: gmm          !< Gaussian mixture model
         real(RP), allocatable :: x(:,:)       !< Feature space for each point
         integer,  allocatable :: clusters(:)  !< Cluster ID for each point
    
    procedure(Compute_Int), pointer :: Compute_Raw => null() 

         contains
         procedure :: Describe => Describe_VISsensor
         procedure :: Compute  => Compute_VISsensor
         final     :: Destruct_VISsensor
   end type VISsensor_t
   !
!  Interfaces
!  ----------
   abstract interface
      subroutine Compute_Int(this, sem, t)
         import RP, VISsensor_t, DGsem
         class(VISsensor_t), target, intent(inout) :: this
         type(DGSem),       target, intent(inout) :: sem
         real(RP),                  intent(in)    :: t
      end subroutine Compute_Int
   end interface
!  ========
   contains
!  ========
subroutine Set_VISsensor(sensor, controlVariables, sem, minSteps, &
    ComputeTimeDerivative, ComputeTimeDerivativeIsolated)   
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
type(VISsensor_t),                  intent(inout) :: sensor
type(FTValueDictionary),            intent(in)    :: controlVariables
type(DGSem),                        intent(in)    :: sem
integer,                            intent(in)    :: minSteps
procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivative
procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivativeIsolated
integer                                           :: nclusters
character(len=:), allocatable                     :: sensorType

sensor % min_steps = minSteps

if (controlVariables % containsKey(VIS_SENSOR_KEY)) then
    sensorType = controlVariables % stringValueForKey(VIS_SENSOR_KEY, LINE_LENGTH)
end if 
select case (sensorType)
case (VIS_GMM_VAL)
    sensor % sens_type = VIS_GMM_ID
    sensor % Compute_Raw => SensorViscous_GMM
    allocate(sensor % x(3, sem % NDOF))
    allocate(sensor % clusters(sem % NDOF))
    if (controlVariables % containsKey(VIS_NUM_CLUSTERS_KEY)) then
       nclusters = controlVariables % doublePrecisionValueForKey(VIS_NUM_CLUSTERS_KEY)
    else
       nclusters = 2
    end if
    if (controlVariables % containsKey(VIS_SENSOR_ELEM_CLUST)) then
      sensor % e_clust = controlVariables % integerValueForKey(VIS_SENSOR_ELEM_CLUST)
    else
      sensor % e_clust = 0
    end if  
     call Stopwatch % CreateNewEvent("Clustering")
    call sensor % gmm % init(3, nclusters, logl_tol=1e-4_RP,zero_tol= 1e-11_RP)
end select
if (sensor % sens_type == VIS_GMM_ID) then
      sensor % sVar = VIS_INV_ID
end if
end subroutine Set_VISsensor
subroutine Describe_VISsensor(sensor)
    !
    !     ---------
    !     Interface
    !     ---------
          implicit none
          class(VISsensor_t), intent(in) :: sensor
    
    
          write(STD_OUT,"(30X,A,A52)", advance="no") "->", "Sensor type: "
          select case (sensor % sens_type)
          case (VIS_GMM_ID);          write(STD_OUT,"(A)") VIS_GMM_VAL
          end select
          write(STD_OUT,"(30X,A,A52)", advance="no") "->", "Sensed variable: "
          select case (sensor % sVar)
          case (VIS_INV_ID ); write(STD_OUT,"(A)") VIS_INV_VAL   
          end select
          write(STD_OUT,"(30X,A,A52,I0,A)") "->", "Sensor inertia: ", sensor % min_steps, " timesteps"
end subroutine Describe_VISsensor
subroutine Destruct_VISsensor(sensor)
    !
    !     ---------
    !     Interface
    !     ---------
          type(VISsensor_t), intent(inout) :: sensor
    

          if (allocated(sensor % x))            deallocate(sensor % x)
          if (allocated(sensor % clusters))     deallocate(sensor % clusters)
          if (associated(sensor % Compute_Raw)) nullify(sensor % Compute_Raw)
    
end subroutine Destruct_VISsensor
subroutine Compute_VISsensor(sensor, sem, t)
        !
        !     ---------
        !     Interface
        !     ---------
              implicit none
              class(VISsensor_t), target, intent(inout) :: sensor
              type(DGsem),       target, intent(inout) :: sem
              real(RP),                  intent(in)    :: t
        !
        !     ---------------
        !     Local variables
        !     ---------------
        !
              type(Element), pointer :: e
              integer                :: eID
              real(RP)               :: s
        
        !
        !     Compute the sensor
        !     ------------------
              call sensor % Compute_Raw(sem, t)
        !
        !     Add 'inertia' to the scaled value
        !     ---------------------------------
              if (sensor % min_steps > 1) then   ! Enter the loop only if necessary
        !$omp parallel do default(private) shared(sem)
                 do eID = 1, sem % mesh % no_of_elements
                    e => sem % mesh % elements(eID)
                    s = e % storage % sensor
                    if (s > 0.0_RP) then
                       if (e % storage % prev_sensor <= 0.0_RP) then
                          e % storage % first_sensed = 0
                          e % storage % prev_sensor = s
                       else
                          e % storage % first_sensed = e % storage % first_sensed + 1
                          e % storage % prev_sensor = s
                       end if
                    elseif (e % storage % first_sensed < sensor % min_steps) then
                       e % storage % first_sensed = e % storage % first_sensed + 1
                       e % storage % sensor = e % storage % prev_sensor
                    end if
                 end do
        !$omp end parallel do
              end if
        
              nullify(e)
        
end subroutine Compute_VISsensor
!
!///////////////////////////////////////////////////////////////////////////////
!      Subroutine for GMM based detection of viscous dominated rotational
!      regions
!///////////////////////////////////////////////////////////////////////////////
!
subroutine SensorViscous_GMM(sensor,sem,t)
    use VariableConversion, only: getVelocityGradients
 
    use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
    use FluidData,          only: thermodynamics, dimensionless
 
    !
 !     ---------
 !     Interface
 !     ---------
    implicit none
    class(VISsensor_t), target, intent(inout) :: sensor
    type(DGSem),       target, intent(inout) :: sem
    real(RP),                  intent(in)    :: t
 !
 !     ---------------
 !     Local variables
 !     ---------------
    type(Element), pointer :: e
    integer                :: eID
    integer                :: i, j, k, l
    integer                :: cnt
    integer                :: n
    integer                :: cluster
    real(RP)               :: ux(3), uy(3), uz(3)
    real(RP)               :: P1,P2
 !
 !     Compute the clustering variables and store them in a global array
 !     -----------------------------------------------------------------

cnt = 0
! 
    do eID = 1, sem % mesh % no_of_elements
 
       e => sem % mesh % elements(eID)
 
       do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
 
          cnt = cnt + 1
 
          
          call getVelocityGradients(     &
             e % storage % Q(:,i,j,k),   &
             e % storage % U_x(:,i,j,k), &
             e % storage % U_y(:,i,j,k), &
             e % storage % U_z(:,i,j,k), &
             ux, uy, uz                  &
          )
     sensor % x(1,cnt)= ux(1)*uy(2)+uy(2)*uz(3)+ux(1)*uz(3)-(0.5*(ux(2)+uy(1)))**2-(0.5*(uz(3)+uz(2)))**2-(0.5*(uz(1)+ux(3)))**2  ! 2nd invariant of strain rate tensor
     sensor % x(2,cnt)= -0.25*uy(2)*(ux(3)+uz(1))**2+0.25*(ux(2)+uy(1))*(uy(3)+uz(2))*(uz(1)+ux(1))-0.25*ux(1)*(uz(2)+uy(3)) ! 3rd invariant of strain rate tensor
     sensor % x(3,cnt)= -0.25*(uy(1)-ux(2))**2-0.25*(uz(1)-ux(3))**2-0.25*(uz(2)-uy(3))**2 ! 2nd invariant of rotational rate tensor

           
       end do   ; end do   ; end do 
    end do     
    !
 !     Rescale the values
 !     ------------------
    
    
    !     Compute the GMM clusters
    call Stopwatch % Start("Clustering")
    call sensor % gmm % fit(sensor % x, reset= .false.,adapt= .false. , from_kmeans= .true.)
    call sensor % gmm % predict(sensor % x)
    call Stopwatch % Pause("Clustering")
     
    

if (sensor % e_clust == 0) then    ! the element is assigned to the cluster that has maximum numbers in that element  
   cnt = 0
   do eID = 1, sem % mesh % no_of_elements
       e => sem % mesh % elements(eID)
       if (sensor % gmm % nclusters <= 1) then
          e % storage % sensor = 0.0_RP
       else
         n = product(e % Nxyz + 1)
         cluster = maxval(maxloc(sensor % gmm % prob(cnt+1:cnt+n,:), dim=2))
         e % storage % sensor = real(cluster - 1, RP) / (sensor % gmm % nclusters - 1)
       end if
   cnt = cnt + n    
   end do       
else if (sensor % e_clust == 1) then ! the mean method is used to decide whether the element is viscous or inviscid
   cnt = 0
   do eID = 1, sem % mesh % no_of_elements
      e => sem % mesh % elements(eID)
      if (sensor % gmm % nclusters <= 1) then
         e % storage % sensor = 0.0_RP
      else
         n = product(e % Nxyz + 1)
         P1=sum(sensor % gmm % prob(cnt+1:cnt+n,1))/n ! mean of probability memberships of being in Cluster 1
         P2=sum(sensor % gmm % prob(cnt+1:cnt+n,2))/n ! mean of probability memberships of being in Cluster 2   
            if (P2 < P1) then   ! the element is assigned to the cluster with maximum of mean probability membership       
            e % storage % sensor = 0.0_RP
            else if (P2 >= P1)  then 
            e % storage % sensor = 1.0_RP
            end if
      end if       
   cnt = cnt + n
   end do
end if   

    nullify(e)    
        
 
  end subroutine SensorViscous_GMM
  !///////////////////////////////////////////////////////////////////////////////
end module VisSensorClass
    

!

    
