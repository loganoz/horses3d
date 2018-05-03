module MultigridTypes
   use SMConstants
   use HexMeshClass
   use InterpolationMatrices
   use TimeIntegratorDefinitions
   use DGSEMClass
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
      subroutine SmoothIt_t( mesh, particles, t, BCFunctions, deltaT, ComputeTimeDerivative )
         use SMConstants, only: RP
         use HexMeshClass, only: HexMesh
         use DGSEMClass, only: ComputeQDot_FCN, BCFunctions_t, no_of_BCsets
#if defined(NAVIERSTOKES)
         use ParticlesClass, only: Particles_t
#endif
         IMPLICIT NONE
         type(HexMesh)              :: mesh
#if defined(NAVIERSTOKES)
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
         REAL(KIND=RP)              :: t, deltaT
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
         procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      end subroutine SmoothIt_t
   end interface
   
!  Parameters
!  ----------
   
   ! Node type
   integer, parameter :: NMIN_GAUSS        = 1 ! Minimum polynomial order when using Gauss nodes .... The threshold should actually be zero... Using 1 because code doesn't support 0
   integer, parameter :: NMIN_GAUSSLOBATTO = 1 ! Minimum polynomial order when using Gauss-Lobatto nodes
   
   ! smoothers
   integer, parameter :: RK3_SMOOTHER     = 0 ! Williamson's 3rd order low-storage Runge-Kutta (only for steady state cases)
   integer, parameter :: BJ_SMOOTHER      = 1 ! Block Jacobi smoother
   integer, parameter :: JFGMRES_SMOOTHER = 2 ! Jacobian-Free GMRES
   integer, parameter :: IMPLICIT_SMOOTHER_IDX = 1 ! All smoothers with index >= IMPLICIT_SMOOTHER_IDX are implicit
   
   ! Variables for IO
   integer        :: ThisTimeStep   ! Current time step
   integer        :: plotInterval   ! Read to display output
   
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
      
      call ConstructInterpolationMatrices(N1, N2)
      
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
      real(kind=RP)             :: maxResiduals(N_EQN)
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

         write(STD_OUT,'(A,A,I3,X,A,X,A,I8)',advance="no") color1,'FAS lvl', lvl ,"|","it",sweeps

         do eqn = 1, N_EQN
            write(STD_OUT ,'(X,A,X,ES10.3)',advance="no") "|", maxResiduals(eqn)
         end do
   
         write(STD_OUT,'(A)') color2
      end if
      
   end subroutine PlotResiduals
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MultigridTypes
