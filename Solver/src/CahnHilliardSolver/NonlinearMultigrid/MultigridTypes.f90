!
!//////////////////////////////////////////////////////
!
!   @File:    MultigridTypes.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:43 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module MultigridTypes
   use SMConstants
   use DGSEMClass
   use InterpolationMatrices
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
      subroutine SmoothIt_t(sem, t, dt)
         use DGSEMClass
         type(DGSem)   :: sem
         real(kind=RP) :: t, dt
      end subroutine SmoothIt_t
   end interface
   
   ! Parameters
   integer, parameter :: NMIN_GAUSS        = 1 ! Minimum polynomial order when using Gauss nodes .... The threshold should actually be zero... Using 1 because code doesn't support 0
   integer, parameter :: NMIN_GAUSSLOBATTO = 1 ! Minimum polynomial order when using Gauss-Lobatto nodes
   
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
      if (N2 > (NMAX - deltaN*(MGlevels - lvl))) N2 = N2 - deltaN
      
      ! The order must be greater or equal to 0 (Legendre-Gauss quadrature) or 1 (Legendre-Gauss-Lobatto)
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
   subroutine PlotResiduals( lvl, sweeps , sem, white )
      use PhysicsStorage
      implicit none
      !--------------------------------------------------------
      integer    , intent(in)           :: lvl
      integer    , intent(in)           :: sweeps
      type(DGSem), intent(in)           :: sem
      logical    , intent(in), optional :: white
      !--------------------------------------------------------
      real(kind=RP)             :: maxResiduals(N_EQN)
      character(len=5)          :: color1
      character(len=5)          :: color2
      !--------------------------------------------------------
      
      if (present(white) .AND. white) then
         color1 = achar(27)//'[00m'
      else
         color1 = achar(27)//'[34m'
      end if
      color2 = achar(27)//'[00m'
      
      if( (MOD( ThisTimeStep+1, plotInterval) == 0) .or. (ThisTimeStep .eq. 0) ) then
         maxResiduals = ComputeMaxResidual(sem)
         write(STD_OUT , 110) color1,'FAS lvl', lvl ,"|","it",sweeps,"|", maxResiduals(IRHO) , "|" , maxResiduals(IRHOU) , &
                                 "|", maxResiduals(IRHOV) , "|" , maxResiduals(IRHOW) , "|" , maxResiduals(IRHOE),color2
      end if
      
      110 format (A,A,I3,X,A,X,A,I8,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,A)
      
   end subroutine PlotResiduals
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MultigridTypes
