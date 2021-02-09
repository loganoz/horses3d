!
!//////////////////////////////////////////////////////
!
!   @File:    MultigridTypes.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Sun Apr 27 12:57:00 2017
!   @Last revision date: Wed Jan 27 16:23:12 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: e199f09aa7589b8bf0cca843e5f1caf3e59586af
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
   end type MGSolStorage_t
   
!
!  Interface for the smoother
!  --------------------------
   abstract interface
      subroutine SmoothIt_t( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec)
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
   integer, parameter :: RK5OPT_SMOOTHER  = 3 ! 
   integer, parameter :: IMPLICIT_SMOOTHER_IDX = 4 ! All smoothers with index >= IMPLICIT_SMOOTHER_IDX are implicit
   integer, parameter :: BJ_SMOOTHER      = 4 ! Block Jacobi smoother
   integer, parameter :: JFGMRES_SMOOTHER = 5 ! Jacobian-Free GMRES
   
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
      
      if (present(white) .AND. white) then
         color1 = achar(27)//'[00m'
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
   subroutine CFLRamp(cfl_ini,cfl,n,res0,res1,CFLboost)
      IMPLICIT NONE
!
!     ------------------------------------
!     Variation of CFL ramping according to Jiang et al. 2015 (version for FAS). 
!     ------------------------------------
!
!     ----------------------------------------------
      real(kind=rp), intent(in)       :: cfl_ini ! initial, user-defined CFL
      real(kind=rp), intent(inout)    :: cfl     ! previous iteration CFL
      integer, intent(in)             :: n ! outer iteration
      real(kind=rp), intent(in)       :: res0     ! RES before
      real(kind=rp), intent(in)       :: res1     ! RES after
      logical, intent(in)             :: CFLboost
!     ----------------------------------------------
      real(kind=rp)                :: eta = 1.01d0   ! Ideally 1.0 < eta < 1.05
      real(kind=rp)                :: eps = 1e-6     !
!     ----------------------------------------------
      if (CFLboost) then
         if ( (res1 / res0 .gt. 1.0d0) .and. (res0 .gt. eps) ) then
            cfl = cfl / 2.0d0
         else
            ! cfl = cfl_ini * eta**n ! original work
            cfl = cfl_ini  * (1.d0 + (eta-1.d0)*n) ! linear variation 
         end if
      end if ! CLFBoost 
   
   end subroutine CFLRamp
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MultigridTypes