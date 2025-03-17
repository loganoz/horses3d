!
!//////////////////////////////////////////////////////
!
!  Variables and utilities for mulgridrid solvers.
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "Includes.h"
module MultigridTypes
   use SMConstants
   use HexMeshClass              , only: HexMesh
   use InterpolationMatrices     , only: Tset
   use TimeIntegratorDefinitions
   use DGSEMClass                , only: ComputeMaxResiduals
   use MPI_Process_Info          , only: MPI_Process
   implicit none

   public
!
!  Multigrid solution storage
!  --------------------------
   type :: MGSolStorage_t
      real(kind=RP), dimension(:,:,:,:), allocatable :: Q     ! Solution of the conservative level (before the smoothing)
      real(kind=RP), dimension(:,:,:,:), allocatable :: E     ! Error (for correction)
      real(kind=RP), dimension(:,:,:,:), allocatable :: S     ! Source term interpolated from next finer grid
      real(kind=RP), dimension(:,:,:,:), allocatable :: Scase ! Source term from the specific case that is running (this is actually not necessary for the MG scheme, but it's needed to estimate the truncation error) .. Currently, it only considers the source term from manufactured solutions (state of the code when this module was written)
      real(kind=RP), dimension(:,:,:,:), allocatable :: R     !
      real(kind=RP), dimension(:,:,:,:), allocatable :: Q0    !
   end type MGSolStorage_t

!
!  Interface for the smoother
!  --------------------------
   abstract interface
      subroutine SmoothIt_t( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt)
         use SMConstants, only: RP
         use HexMeshClass, only: HexMesh
         use DGSEMClass, only: ComputeTimeDerivative_f
#ifdef FLOW
         use ParticlesClass, only: Particles_t
#endif
         IMPLICIT NONE
         type(HexMesh)      :: mesh
#ifdef FLOW
         type(Particles_t)  :: particles
#else
         logical            :: particles
#endif
         REAL(KIND=RP)              :: t, deltaT
         procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
         real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
         ! Dual (pseudo) time stepping arguments (also optional):
         logical, intent(in), optional :: dts
         real(kind=RP), intent(in), optional :: global_dt
      end subroutine SmoothIt_t
   end interface

!  Parameters
!  ----------

   ! Node type
   integer, parameter :: NMIN_GAUSS        = 1 ! Minimum polynomial order when using Gauss nodes .... The threshold should actually be zero... Using 1 because code doesn't support 0
   integer, parameter :: NMIN_GAUSSLOBATTO = 1 ! Minimum polynomial order when using Gauss-Lobatto nodes

   ! smoothers
   integer, parameter :: Euler_SMOOTHER   = 0 !
   integer, parameter :: RK3_SMOOTHER     = 1 ! Williamson's 3rd order low-storage Runge-Kutta (only for steady state cases)
   integer, parameter :: RK5_SMOOTHER     = 2 !
   integer, parameter :: RKOPT_SMOOTHER   = 3 !
   integer, parameter :: IMPLICIT_SMOOTHER_IDX = 4 ! All smoothers with index >= IMPLICIT_SMOOTHER_IDX are implicit
   integer, parameter :: IRK_SMOOTHER = 4     ! Implicit Euler smoother (full matrix assembly)
   integer, parameter :: BIRK5_SMOOTHER = 5   ! Semi-implicit RK smoother (5-stage Block-Jacobi)
   integer, parameter :: SGS_SMOOTHER = 6     ! symmetric Gauss-Seidel smoother
   integer, parameter :: ILU_SMOOTHER = 7     ! symmetric Gauss-Seidel smoother

   ! Variables for IO
   integer        :: ThisTimeStep   ! Current time step
   integer        :: plotInterval   ! Read to display output

   ! preconditioners
   integer, parameter :: PRECONDIIONER_NONE = 0
   integer, parameter :: PRECONDIIONER_LTS = 1 ! Local time stepping

!========
 contains
!========

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Creates the restriction and prolongation operators for a certain element
!  for multigrid. Takes into account order anisotropy, but the coarse grid
!  is constructed by reducing the polynomial order uniformly.
!  ------------------------------------------------------------------------
   subroutine CreateInterpolationOperators(N1,N2,NMAX,MGlevels,lvl,DeltaN,nodes)
      use NodalStorageClass
      implicit none
      !-----------------------------------------------------
      integer, intent(in)     :: N1           !<  Fine grid order(anisotropic) of the element
      integer, intent(out)    :: N2           !>  Coarse grid order(anisotropic) of the element
      integer, intent(in)     :: NMAX         !<  Maximum polynomial order of multigrid representation
      integer, intent(in)     :: MGlevels     !<  Number of multigrid levels of representation
      integer, intent(in)     :: lvl          !<  Current multigrid level
      integer, intent(in)     :: DeltaN       !<  Interval of reduction of polynomial order for coarser level
      integer, intent(in)     :: nodes        !<  Is the quadrature a Legendre-Gauss-Lobatto representation?
      !-----------------------------------------------------
!
!     --------------------------------------
!     Compute order of coarse representation
!     --------------------------------------
!
      ! Uniform coarsening (not used currently)
!~      N2 = N1 - DeltaN

      ! High order coarsening (see paper p-adaptation + Multigrid)
      N2 = N1
      if ( N2 > (NMAX - deltaN*(MGlevels - lvl)) ) N2 = NMAX - deltaN*(MGlevels - lvl)

      ! Check the minimum polynomial order
      if (nodes == GAUSSLOBATTO) then
         if (N2 < NMIN_GAUSSLOBATTO) N2 = NMIN_GAUSSLOBATTO
      else
         if (N2 < NMIN_GAUSS)        N2 = NMIN_GAUSS
      end if

      call NodalStorage(N2) % construct( nodes, N2)

      call Tset(N1, N2) % construct(N1, N2)
      call Tset(N2, N1) % construct(N2, N1)
   end subroutine CreateInterpolationOperators

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!       Multigrid IO
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Internal subroutine to print the residuals
!  ------------------------------------------
   subroutine PlotResiduals( lvl, sweeps , mesh, white )
      use PhysicsStorage
      implicit none
      !--------------------------------------------------------
      integer    , intent(in)           :: lvl
      integer    , intent(in)           :: sweeps
      type(HexMesh), intent(in)           :: mesh
      logical    , intent(in), optional :: white
      !--------------------------------------------------------
      real(kind=RP)             :: maxResiduals(NCONS)
      character(len=5)          :: color1
      character(len=5)          :: color2
      integer                   :: eqn
      !--------------------------------------------------------

      if (present(white)) then
         if (white) then
            color1 = achar(27)//'[00m'
         else
            color1 = achar(27)//'[34m'
         end if
      else
         color1 = achar(27)//'[34m'
      end if
      color2 = achar(27)//'[00m'

      if( (MOD( ThisTimeStep+1, plotInterval) == 0) .or. (ThisTimeStep .eq. 0) ) then
         maxResiduals = ComputeMaxResiduals(mesh)

         if (MPI_Process % isRoot) then
            write(STD_OUT,'(A,A,I3,X,A,X,A,I8)',advance="no") color1,'FAS lvl', lvl ,"|","it",sweeps

            do eqn = 1, NCONS
               write(STD_OUT ,'(X,A,X,ES10.3)',advance="no") "|", maxResiduals(eqn)
            end do

            write(STD_OUT,'(A)') color2
         end if
      end if

   end subroutine PlotResiduals
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CFLRamp(cfl_max,cfl,cflboost_rate,CFLboost)
      IMPLICIT NONE
!
!     ------------------------------------
!      CLF ramping strategies.
!     ------------------------------------
!
!     ----------------------------------------------
      real(kind=rp), intent(in)              :: cfl_max
      real(kind=rp), intent(inout)           :: cfl
      real(kind=rp), intent(in)              :: cflboost_rate
      character(len=LINE_LENGTH), intent(in) :: CFLboost
!     ----------------------------------------------
      real(kind=rp)                :: conv
!     ----------------------------------------------


      if ( trim(CFLboost) .eq. "linear") then
         cfl = cfl + cfl * cflboost_rate
         if (cfl .ge. cfl_max) then
            cfl = cfl_max
        end if
      elseif ( trim(CFLboost) .eq. "exponential") then
         cfl = cfl * cflboost_rate
         if (cfl .ge. cfl_max) then
            cfl = cfl_max
        end if
      end if ! CLFBoost

         ! real(kind=rp), intent(in)       :: res0     ! RES before
         ! real(kind=rp), intent(in)       :: res1     ! RES after
         !  conv = log10(res0/res1)
         !  if (conv .le. 0.0_RP) then
         !  elseif ( (conv .gt. 0.0_RP) .and. (conv .le. 1.0_RP )  ) then
         !      cfl = cfl + cfl * conv * cflboost_rate
         !  elseif (conv .gt. 1.0_RP ) then
         !      cfl = cfl + cfl * cflboost_rate
         !  end if

   end subroutine CFLRamp
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MultigridTypes