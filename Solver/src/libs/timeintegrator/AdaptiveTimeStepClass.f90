!
!//////////////////////////////////////////////////////
!
!      Class cotaining routines for adapting time step based on Reinforcement Learning.
!        -> The adaptation procedure is performed with a RL agent trined with a Value Iteration algorithm.
!        -> The current implementation is compatible with OpenMP and MPI.
!
!////////////////////////////////////////////////////////////////////////
!
module AdaptiveTimeStepClass
   use SMConstants
#ifdef NAVIERSTOKES
   use PhysicsStorage                  , only: CTD_IGNORE_MODE, flowIsNavierStokes
   use FluidData                       , only: thermodynamics
#elif defined(MULTIPHASE)
   use PhysicsStorage                  , only: CTD_IGNORE_MODE, IMP
#else
   use PhysicsStorage                  , only: CTD_IGNORE_MODE
#endif
   use ElementClass
   use HexMeshClass
   USE DGSEMClass
   use FTValueDictionaryClass          , only: FTValueDictionary
   use StorageClass
   use MPI_Process_Info
   
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
#include "Includes.h"
   private
   public adaptiveTimeStep_t

   type :: adaptiveTimeStep_t
      real(kind=RP) :: min_dt = 1e-8_RP ! Minimum time step
      real(kind=RP) :: max_dt = 1e-2_RP ! Maximum time step
      real(kind=RP) :: b1     = 0.0 !PID 1st parameter
      real(kind=RP) :: b2     = -0.2 !PID 2nd parameter
      real(kind=RP) :: b3     = 0.2 !PID 3rd parameter
      real(kind=RP) :: k      = 3.0 !Order of the RK method (only RK3 is implemented)
      real(kind=RP) :: dt_eps(3)    !Epsilon for the PID controller
      logical       :: constructed = .FALSE. !If the adaptive time step is constructed
      real(kind=RP) :: dtAdaptationStep = 1e-3_RP !Adaptation step
      real(kind=RP) :: lastAdaptationTime = 0.0_RP !Last adaptation time
      integer       :: error_variable !1:u, 2:v, 3:w, 4:rho*u, 5:rho*v, 6:rho*w, 7:p (only for Navier-Stokes)
      logical       :: pAdapted = .FALSE.

      contains 
         procedure :: construct      => adaptiveTimeStep_Construct
         procedure :: destruct       => adaptiveTimeStep_Destruct
         procedure :: hasToAdapt     => adaptiveTimeStep_HasToAdapt
         procedure :: update         => adaptiveTimeStep_Update
   end type adaptiveTimeStep_t

contains 

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------
!  Routine for constructing the adaptive time step
!  -----------------------------------------------
   subroutine adaptiveTimeStep_Construct(this, controlVariables, t0)
      implicit none
      class(adaptiveTimeStep_t), intent(inout) :: this
      type(FTValueDictionary)  , intent(in)    :: controlVariables
      real(kind=RP)            , intent(in)    :: t0
!!     ---------------
      this % dt_eps(1) = 1.0_RP
      this % dt_eps(2) = 1.0_RP
      this % dt_eps(3) = 1.0_RP
      this % constructed = .TRUE.
      this % lastAdaptationTime = t0

      if (controlVariables % containsKey("dt adaptation step")) then
         this % dtAdaptationStep = controlVariables % doublePrecisionValueForKey("dt adaptation step")
      end if

      if (controlVariables % containsKey("minimum dt")) then
         this % min_dt = controlVariables % doublePrecisionValueForKey("minimum dt")
      end if

      if (controlVariables % containsKey("maximum dt")) then
         this % max_dt = controlVariables % doublePrecisionValueForKey("maximum dt")
      end if
   end subroutine adaptiveTimeStep_Construct
!
!
!  -----------------------------------------------
!  Subroutine for destructing the adaptive time step
!  -----------------------------------------------
   subroutine adaptiveTimeStep_Destruct(this)
      implicit none
      class(adaptiveTimeStep_t), intent(inout) :: this
!!     ---------------
   end subroutine adaptiveTimeStep_Destruct
!
!  -----------------------------------------------
!  Subroutine for updating the adaptive time step
!  -----------------------------------------------
   subroutine adaptiveTimeStep_Update(this, mesh, t, dt)
      implicit none
      class(adaptiveTimeStep_t) , intent(inout)    :: this
      class(HexMesh)            , intent(inout)    :: mesh
      real(kind=RP)             , intent(in)       :: t
      real(kind=RP)             , intent(inout)    :: dt
!!     ---------------
!     Local variables
!!     ---------------
      integer                                :: eID, ierr, i
      real(kind=RP)                          :: dt_weight, sum_dt_weight, avg_sum_dt_weight, dt_lim, min_dt_lim, max_dt_lim

      min_dt_lim = 0.5_RP ! Minimum limit for the time step
      max_dt_lim = 1.3_RP ! Maximum limit for the time step

      this % lastAdaptationTime = t
      dt_weight = 0.0_RP
      sum_dt_weight = 0.0_RP
      avg_sum_dt_weight = 0.0_RP
!$omp parallel shared(dt_weight, mesh, this)
!$omp do reduction(+:dt_weight) schedule(runtime)
      do eID = 1, mesh % no_of_elements
         dt_weight = dt_weight + adaptiveTimeStep_ComputeWeights(this, mesh % elements(eID))
      end do
!$omp end do
!$omp end parallel

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         call mpi_allreduce ( dt_weight, sum_dt_weight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
      else
         sum_dt_weight = dt_weight
      end if  

      avg_sum_dt_weight = sum_dt_weight / mesh % no_of_allElements

!$omp parallel do schedule(runtime)
      do eID = 1, mesh % no_of_elements
         mesh % elements(eID) % ML_error_ratio = mesh % elements(eID) % ML_error_ratio / avg_sum_dt_weight
      end do
!$omp end parallel do

      if (isnan(sum_dt_weight)) then
         this % dt_eps(3) = 1e-10_RP
      else
         this % dt_eps(3) = 1.0_RP / max(sqrt(avg_sum_dt_weight), 1e-10_RP)
      end if

      if ((this % dt_eps(1) == this % dt_eps(2)) .and. (this % dt_eps(2) == 1.0_RP)) then
         this % dt_eps(1) = this % dt_eps(3) ! Initialize the first epsilon value
         this % dt_eps(2) = this % dt_eps(3) ! Initialize the second epsilon value
      else
         this % dt_eps(3) = min(this % dt_eps(3), this % dt_eps(2) * 10.0_RP) ! Ensure that the third epsilon is not too large
      end if

      if (this % dt_eps(2) > this % dt_eps(3)) then
         if (this % pAdapted) then !Right after p-Adaptation
            this % b3 = 0.15_RP 
            this % b2 = -0.15_RP 
            this % pAdapted = .FALSE.
         else
            this % b3 = 0.2_RP 
            this % b2 = -0.2_RP 
         end if
      else
         this % b3 = 0.2_RP 
         if (this % dt_eps(2) < 1.0_RP) then
            this % b2 = -0.1_RP
         else if (this % dt_eps(2) / this % dt_eps(3) < 0.5_RP) then
            if (this % pAdapted) then
               this % b3 = 0.15_RP 
               this % b2 = max(-0.15_RP * (1.0_RP + 0.01_RP * this % dt_eps(3) / this % dt_eps(2)), -0.2_RP)
               this % pAdapted = .FALSE.
            else
               this % b2 = max(-0.2_RP * (1.0_RP + 0.01_RP * this % dt_eps(3) / this % dt_eps(2)), -0.25_RP) 
            end if
         else
            this % b2 = -0.2_RP 
         end if
      end if

      dt_lim = dt_limiter(this % dt_eps(3) ** (this % b3/this % k) * &
                        this % dt_eps(2) ** (this % b2/this % k) * &
                        this % dt_eps(1) ** (this % b1/this % k))

      if (dt_lim > 1.05_RP .or. dt_lim < 0.95_RP) then

         this % dt_eps(1) = this % dt_eps(2)
         this % dt_eps(2) = this % dt_eps(3)

         if (dt_lim < 0.85_RP) then
            dt_lim = dt_lim * 0.9_RP
         else if (dt_lim < 0.95_RP) then
            dt_lim = dt_lim ** 1.2_RP
         else if (dt_lim < 1.0_RP) then
            dt_lim = dt_lim
         else
            if (this % dt_eps(3) > 100.0_RP) then 
               dt_lim = dt_lim ** (1.1_RP + this % dt_eps(3) / 10000.0_RP)
            else
               dt_lim = dt_lim * (1.0_RP - 0.2_RP / this % dt_eps(3))
            end if
         end if

         dt_lim = min(max(dt_lim, min_dt_lim), max_dt_lim) ! Limit the time step
         dt = dt * dt_lim

         dt = min(max(dt, this % min_dt), this % max_dt) ! Ensure dt is not below the minimum value or above the maximum value

         this % dtAdaptationStep = dt

      end if !end if dt_lim > 1.05_RP .or. dt_lim < 0.95_RP

   end subroutine AdaptiveTimeStep_Update

   function adaptiveTimeStep_ComputeWeights(this, e) result(dt_weight)
      implicit none 
      class(adaptiveTimeStep_t), intent(in)    :: this
      type(Element)         , intent(inout)    :: e
      real(kind=RP)                            :: dt_weight
!!     ---------------
      integer               :: i, j, k, dir, Ndir
      integer               :: Pxyz(3)
      real(kind=RP)         :: Q_error(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)), QLowRK_error(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))

#ifdef FLOW
!     Initialization of P
      Pxyz = e % Nxyz

      Ndir = 3

      dt_weight = 0.0_RP

      !Select the error variable
      do i = 0, Pxyz(3)
         do j = 0, Pxyz(2)
            do k = 0, Pxyz(1)
               do dir = 1, Ndir
                  if(mod(this % error_variable, 3) == mod(dir, 3) .and. dir <= 3) then
                     if (this % error_variable <= 6) then
                        Q_error(k, j, i) = e % storage % Q(dir+1, k, j, i)
                        QLowRK_error(k, j, i) = e % storage % QLowRK(dir+1, k, j, i)
                        if (this % error_variable <= 3) then
                           Q_error(k, j, i) = Q_error(k, j, i) / e % storage % Q(1, k, j, i)
                           QLowRK_error(k, j, i) = QLowRK_error(k, j, i) / e % storage % QLowRK(1, k, j, i)
                        end if
                     else if (this % error_variable == 7) then
#ifdef NAVIERSTOKES
                        !Pressure
                        Q_error(k, j, i) = thermodynamics % gammaMinus1*(e % storage % Q(5,k,j,i) - 0.5_RP*(e % storage % Q(2,k,j,i)**2 + e % storage % Q(3,k,j,i)**2 + e % storage % Q(4,k,j,i)**2)/e % storage % Q(1,k,j,i))
                        QLowRK_error(k, j, i) = thermodynamics % gammaMinus1*(e % storage % QLowRK(5,k,j,i) - 0.5_RP*(e % storage % QLowRK(2,k,j,i)**2 + e % storage % QLowRK(3,k,j,i)**2 + e % storage % QLowRK(4,k,j,i)**2)/e % storage % QLowRK(1,k,j,i))
#elif defined(MULTIPHASE)
                        !Pressure
                        Q_error(k, j, i) = e % storage % QNS(IMP,k,j,i)
                        QLowRK_error(k, j, i) = e % storage % QLowRK(IMP,k,j,i)
#endif
                     else if (this % error_variable == 8) then
                        !Density
                        Q_error(k, j, i) = e % storage % Q(1,k,j,i)
                        QLowRK_error(k, j, i) = e % storage % QLowRK(1,k,j,i)
                     end if
                  end if
               end do
               dt_weight = dt_weight + (Q_error(k, j, i) - QLowRK_error(k, j, i))**2.0_RP
            end do
         end do
      end do

      dt_weight = dt_weight / (e % storage % sensor)**2.0_RP
      dt_weight = dt_weight / ((Pxyz(1)+1) * (Pxyz(2)+1) * (Pxyz(3)+1)) !Average over all Gauss points
      ! MLRK correction
      dt_weight = dt_weight / (3.0_RP ** (e % MLevel - 1))**2.0_RP
      e % ML_error_ratio = dt_weight
#endif
   end function adaptiveTimeStep_ComputeWeights

   subroutine adaptiveTimeStep_HasToAdapt(this, t, dt, dtHasToAdapt)
      implicit none
      class(adaptiveTimeStep_t), intent(in) :: this
      real(kind=RP), intent(in) :: t, dt
      logical, intent(out) :: dtHasToAdapt
!!       ---------------
      dtHasToAdapt = (t - this % lastAdaptationTime + dt) >= this % dtAdaptationStep

   end subroutine adaptiveTimeStep_HasToAdapt

   function dt_limiter(a) result(dt)
      implicit none
      real(kind=RP), intent(in) :: a
      real(kind=RP) :: dt
!!       ---------------
      dt = 1.0_RP + atan(a-1.0_RP)
   end function dt_limiter
end module AdaptiveTimeStepClass