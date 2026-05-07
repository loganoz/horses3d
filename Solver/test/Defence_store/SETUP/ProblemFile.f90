#include "Includes.h"

module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f

      subroutine UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use HexMeshClass
         use FluidData
         implicit none
         class(HexMesh)                     :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in) :: multiphase_
#endif
      end subroutine UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use SMConstants
         use PhysicsStorage
         use HexMeshClass
         use FluidData
         implicit none
         class(HexMesh)                     :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in) :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f

#ifdef FLOW
      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(inout)       :: Q(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         use VariableConversion, only: GetGradientValues_f
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(in)          :: Q(NCONS)
         real(kind=RP), intent(inout)       :: U(NGRAD)
         procedure(GetGradientValues_f)     :: GetGradients
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedGradVars_f

      subroutine UserDefinedNeumann_f(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(in)          :: Q(NCONS)
         real(kind=RP), intent(in)          :: U_x(NGRAD)
         real(kind=RP), intent(in)          :: U_y(NGRAD)
         real(kind=RP), intent(in)          :: U_z(NGRAD)
         real(kind=RP), intent(inout)       :: flux(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedNeumann_f
#endif

      subroutine UserDefinedPeriodicOperation_f(mesh, time, dt, monitors)
         use SMConstants
         use HexMeshClass
         use MonitorsClass
         implicit none
         class(HexMesh)               :: mesh
         real(kind=RP)                :: time
         real(kind=RP)                :: dt
         type(Monitor_t), intent(in)  :: monitors
      end subroutine UserDefinedPeriodicOperation_f

#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
                                           , multiphase_ &
#endif
                                           )
         use SMConstants
         use HexMeshClass
         use FluidData
         use PhysicsStorage
         implicit none
         real(kind=RP),             intent(in)    :: x(NDIM)
         real(kind=RP),             intent(in)    :: Q(NCONS)
         real(kind=RP),             intent(in)    :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t),    intent(in)    :: thermodynamics_
         type(Dimensionless_t),     intent(in)    :: dimensionless_
         type(RefValues_t),         intent(in)    :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),        intent(in)    :: multiphase_
#endif
      end subroutine UserDefinedSourceTermNS_f

      subroutine UserDefinedIBMKinematicsNS_f(x, V, cL, cD, t, dt, refValues_, UpdatePosition, GetVelocity)
         use SMConstants
         use FluidData
         use PhysicsStorage
         implicit none
         real(kind=RP),             intent(inout) :: x(NDIM), V(NDIM)
         real(kind=RP),             intent(in)    :: t, dt
         real(kind=RP),             intent(in)    :: cL, cD
         type(RefValues_t),         intent(in)    :: refValues_
         logical,                   intent(in)    :: UpdatePosition, GetVelocity
      end subroutine UserDefinedIBMKinematicsNS_f
#endif

      subroutine UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#ifdef FLOW
                                      , thermodynamics_ &
                                      , dimensionless_  &
                                      , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                      , multiphase_ &
#endif
                                      , monitors, elapsedTime, CPUTime)
         use SMConstants
         use HexMeshClass
         use FluidData
         use MonitorsClass
         implicit none
         class(HexMesh)                     :: mesh
         real(kind=RP)                      :: time
         integer                            :: iter
         real(kind=RP)                      :: maxResidual
#ifdef FLOW
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in) :: multiphase_
#endif
         type(Monitor_t),        intent(in) :: monitors
         real(kind=RP),          intent(in) :: elapsedTime
         real(kind=RP),          intent(in) :: CPUTime
      end subroutine UserDefinedFinalize_f

      subroutine UserDefinedTermination_f
         implicit none
      end subroutine UserDefinedTermination_f
   end interface
end module ProblemFileFunctions

subroutine UserDefinedStartup
   implicit none
end subroutine UserDefinedStartup

subroutine UserDefinedFinalSetup(mesh &
#ifdef FLOW
                                , thermodynamics_ &
                                , dimensionless_  &
                                , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                , multiphase_ &
#endif
                                )
   use HexMeshClass
   use PhysicsStorage
   use FluidData
   implicit none
   class(HexMesh)                     :: mesh
#ifdef FLOW
   type(Thermodynamics_t), intent(in) :: thermodynamics_
   type(Dimensionless_t),  intent(in) :: dimensionless_
   type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
   type(Multiphase_t),     intent(in) :: multiphase_
#endif
end subroutine UserDefinedFinalSetup

subroutine UserDefinedInitialCondition(mesh &
#ifdef FLOW
                                      , thermodynamics_ &
                                      , dimensionless_  &
                                      , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                      , multiphase_ &
#endif
                                      )
   use SMConstants
   use PhysicsStorage
   use HexMeshClass
   use FluidData
   implicit none
   class(HexMesh)                     :: mesh
#ifdef FLOW
   type(Thermodynamics_t), intent(in) :: thermodynamics_
   type(Dimensionless_t),  intent(in) :: dimensionless_
   type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
   type(Multiphase_t),     intent(in) :: multiphase_
#endif

!! integer       :: eID, i, j, k
!! real(kind=RP) :: qq, u, v, w, p
!! #if defined(NAVIERSTOKES)
!!    real(kind=RP) :: Q(NCONS), phi, theta
!! #endif
!! 
!! #if defined(NAVIERSTOKES)
!!    associate(gammaM2 => dimensionless_%gammaM2, gamma => thermodynamics_%gamma)
!!       theta = refValues_%AOAtheta * (PI / 180.0_RP)
!!       phi   = refValues_%AOAphi   * (PI / 180.0_RP)
!! 
!!       do eID = 1, mesh%no_of_elements
!!          associate(Nx => mesh%elements(eID)%Nxyz(1), Ny => mesh%elements(eID)%Nxyz(2), Nz => mesh%elements(eID)%Nxyz(3))
!!             do k = 0, Nz
!!                do j = 0, Ny
!!                   do i = 0, Nx
!!                      qq = 1.0_RP
!!                      u  = qq * cos(theta) * cos(phi)
!!                      v  = qq * sin(theta) * cos(phi)
!!                      w  = qq * sin(phi)
!! 
!!                      Q(1) = 1.0_RP
!!                      p    = 1.0_RP / gammaM2
!!                      Q(2) = Q(1) * u
!!                      Q(3) = Q(1) * v
!!                      Q(4) = Q(1) * w
!!                      Q(5) = p / (gamma - 1.0_RP) + 0.5_RP * Q(1) * (u**2 + v**2 + w**2)
!! 
!!                      mesh%elements(eID)%storage%Q(:, i, j, k) = Q
!!                   end do
!!                end do
!!             end do
!!          end associate
!!       end do
!!    end associate
!! #endif
!! ##################### NEW ############################## !!
   integer       :: eID, i, j, k
   integer       :: loaded_iter
   logical       :: with_gradients
   real(kind=RP) :: loaded_time
   real(kind=RP) :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
   logical, parameter :: use_solution_ic = .true.
   character(len=LINE_LENGTH), parameter :: solution_ic_file = "RESULTS/WPFS_store_IBM_p1_0000430000.hsol"
   real(kind=RP) :: Q(NCONS), phi, theta

   if (use_solution_ic) then
      call mesh % LoadSolution(trim(solution_ic_file), loaded_iter, loaded_time, with_gradients)
      return
   end if

   associate(gammaM2 => dimensionless_%gammaM2, gamma => thermodynamics_%gamma)
      theta = refValues_%AOAtheta * (PI / 180.0_RP)
      phi   = refValues_%AOAphi   * (PI / 180.0_RP)

      do eID = 1, mesh%no_of_elements
         associate(Nx => mesh%elements(eID)%Nxyz(1), Ny => mesh%elements(eID)%Nxyz(2), Nz => mesh%elements(eID)%Nxyz(3))
            do k = 0, Nz
               do j = 0, Ny
                  do i = 0, Nx
                     qq = 1.0_RP
                     u  = qq * cos(theta) * cos(phi)
                     v  = qq * sin(theta) * cos(phi)
                     w  = qq * sin(phi)

                     Q(1) = 1.0_RP
                     p    = 1.0_RP / gammaM2
                     Q(2) = Q(1) * u
                     Q(3) = Q(1) * v
                     Q(4) = Q(1) * w
                     Q(5) = p / (gamma - 1.0_RP) + 0.5_RP * Q(1) * (u**2 + v**2 + w**2)

                     mesh%elements(eID)%storage%Q(:, i, j, k) = Q
                  end do
               end do
            end do
         end associate
      end do
   end associate
#endif

#if defined(INCNS)
   do eID = 1, mesh%no_of_elements
      associate(Nx => mesh%elements(eID)%Nxyz(1), Ny => mesh%elements(eID)%Nxyz(2), Nz => mesh%elements(eID)%Nxyz(3))
         do k = 0, Nz
            do j = 0, Ny
               do i = 0, Nx
                  mesh%elements(eID)%storage%Q(:, i, j, k) = [1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP]
               end do
            end do
         end do
      end associate
   end do
#endif

#ifdef CAHNHILLIARD
   call random_seed()

   do eID = 1, mesh%no_of_elements
      associate(e => mesh%elements(eID)%storage)
         call random_number(e%c)
         e%c = 2.0_RP * (e%c - 0.5_RP)
      end associate
   end do
#endif
end subroutine UserDefinedInitialCondition

#ifdef FLOW
subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
   use SMConstants
   use PhysicsStorage
   use FluidData
   implicit none
   real(kind=RP), intent(in)          :: x(NDIM)
   real(kind=RP), intent(in)          :: t
   real(kind=RP), intent(in)          :: nHat(NDIM)
   real(kind=RP), intent(inout)       :: Q(NCONS)
   type(Thermodynamics_t), intent(in) :: thermodynamics_
   type(Dimensionless_t),  intent(in) :: dimensionless_
   type(RefValues_t),      intent(in) :: refValues_
end subroutine UserDefinedState1

subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
   use SMConstants
   use PhysicsStorage
   use FluidData
   use VariableConversion, only: GetGradientValues_f
   implicit none
   real(kind=RP), intent(in)          :: x(NDIM)
   real(kind=RP), intent(in)          :: t
   real(kind=RP), intent(in)          :: nHat(NDIM)
   real(kind=RP), intent(in)          :: Q(NCONS)
   real(kind=RP), intent(inout)       :: U(NGRAD)
   procedure(GetGradientValues_f)     :: GetGradients
   type(Thermodynamics_t), intent(in) :: thermodynamics_
   type(Dimensionless_t),  intent(in) :: dimensionless_
   type(RefValues_t),      intent(in) :: refValues_
end subroutine UserDefinedGradVars1

subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
   use SMConstants
   use PhysicsStorage
   use FluidData
   implicit none
   real(kind=RP), intent(in)          :: x(NDIM)
   real(kind=RP), intent(in)          :: t
   real(kind=RP), intent(in)          :: nHat(NDIM)
   real(kind=RP), intent(in)          :: Q(NCONS)
   real(kind=RP), intent(in)          :: U_x(NGRAD)
   real(kind=RP), intent(in)          :: U_y(NGRAD)
   real(kind=RP), intent(in)          :: U_z(NGRAD)
   real(kind=RP), intent(inout)       :: flux(NCONS)
   type(Thermodynamics_t), intent(in) :: thermodynamics_
   type(Dimensionless_t),  intent(in) :: dimensionless_
   type(RefValues_t),      intent(in) :: refValues_
end subroutine UserDefinedNeumann1
#endif

subroutine UserDefinedPeriodicOperation(mesh, time, dt, monitors)
   use SMConstants
   use HexMeshClass
   use MonitorsClass
   implicit none
   class(HexMesh)              :: mesh
   real(kind=RP)               :: time
   real(kind=RP)               :: dt
   type(Monitor_t), intent(in) :: monitors
end subroutine UserDefinedPeriodicOperation

#ifdef FLOW
subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
                                  , multiphase_ &
#endif
                                  )
   use SMConstants
   use HexMeshClass
   use PhysicsStorage
   use FluidData
   implicit none
   real(kind=RP),             intent(in)    :: x(NDIM)
   real(kind=RP),             intent(in)    :: Q(NCONS)
   real(kind=RP),             intent(in)    :: time
   real(kind=RP),             intent(inout) :: S(NCONS)
   type(Thermodynamics_t),    intent(in)    :: thermodynamics_
   type(Dimensionless_t),     intent(in)    :: dimensionless_
   type(RefValues_t),         intent(in)    :: refValues_
#ifdef CAHNHILLIARD
   type(Multiphase_t),        intent(in)    :: multiphase_
#endif
end subroutine UserDefinedSourceTermNS
#endif

#ifdef FLOW
subroutine UserDefinedIBMKinematicsNS(x, V, cL, cD, t, dt, refValues_, UpdatePosition, GetVelocity)
   use SMConstants
   use FluidData
   use PhysicsStorage
   implicit none
   real(kind=RP),             intent(inout) :: x(NDIM), V(NDIM)
   real(kind=RP),             intent(in)    :: t, dt
   real(kind=RP),             intent(in)    :: cL, cD
   type(RefValues_t),         intent(in)    :: refValues_
   logical,                   intent(in)    :: UpdatePosition, GetVelocity

   ! Motion parameters in solver units. Edit these three values for the store motion.
   real(kind=RP), parameter :: store_translation_velocity(NDIM) = [0.0_RP, -1.0_RP, 0.0_RP]
   real(kind=RP), parameter :: store_rotation_center(NDIM)      = [4.2618_RP, -0.894_RP, 3.3018_RP]
   real(kind=RP), parameter :: store_omega_z                    = 0.0_RP

   real(kind=RP) :: center_old(NDIM), center_new(NDIM), relative_pos(NDIM)
   real(kind=RP) :: x_new(NDIM), cos_dtheta, sin_dtheta, dtheta

#if defined(NAVIERSTOKES)
   center_old = store_rotation_center + store_translation_velocity * (t - dt)
   center_new = store_rotation_center + store_translation_velocity * t
   dtheta     = store_omega_z * dt

   if (UpdatePosition) then
      relative_pos = x - center_old
      cos_dtheta   = cos(dtheta)
      sin_dtheta   = sin(dtheta)

      x_new(IX) =  center_new(IX) + cos_dtheta * relative_pos(IX) - sin_dtheta * relative_pos(IY)
      x_new(IY) =  center_new(IY) + sin_dtheta * relative_pos(IX) + cos_dtheta * relative_pos(IY)
      x_new(IZ) =  center_new(IZ) + relative_pos(IZ)

      x = x_new
   else if (GetVelocity) then
      relative_pos = x - center_new
      V(IX) = store_translation_velocity(IX) - store_omega_z * relative_pos(IY)
      V(IY) = store_translation_velocity(IY) + store_omega_z * relative_pos(IX)
      V(IZ) = store_translation_velocity(IZ)
   end if
#endif
end subroutine UserDefinedIBMKinematicsNS
#endif

subroutine UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef FLOW
                              , thermodynamics_ &
                              , dimensionless_  &
                              , refValues_ &
#endif
#ifdef CAHNHILLIARD
                              , multiphase_ &
#endif
                              , monitors, elapsedTime, CPUTime)
   use SMConstants
   use HexMeshClass
   use PhysicsStorage
   use FluidData
   use MonitorsClass
   implicit none
   class(HexMesh)                     :: mesh
   real(kind=RP)                      :: time
   integer                            :: iter
   real(kind=RP)                      :: maxResidual
#ifdef FLOW
   type(Thermodynamics_t), intent(in) :: thermodynamics_
   type(Dimensionless_t),  intent(in) :: dimensionless_
   type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
   type(Multiphase_t),     intent(in) :: multiphase_
#endif
   type(Monitor_t),        intent(in) :: monitors
   real(kind=RP),          intent(in) :: elapsedTime
   real(kind=RP),          intent(in) :: CPUTime
end subroutine UserDefinedFinalize

subroutine UserDefinedTermination
   implicit none
end subroutine UserDefinedTermination
