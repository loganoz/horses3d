#include "Includes.h"

#if defined(NAVIERSTOKES) || defined(INCNS)
module HyperbolicSubcell
   use SMConstants
   use HyperbolicDiscretizationClass
   use FluidData
#if defined(SPALARTALMARAS)
   use RiemannSolvers_NSSA
#elif defined(NAVIERSTOKES)
   use RiemannSolvers_NS
#elif defined(INCNS)
   use RiemannSolvers_iNS
#endif

   implicit none

   private
   public SubcellDG_t

   type, extends(HyperbolicDiscretization_t) :: SubcellDG_t
         real(RP) :: max_blend
      contains
         procedure :: Initialize               => SubcellDG_Initialize
         procedure :: ComputeSplitFormFluxes   => SubcellDG_ComputeSplitFormFluxes
         procedure :: ComputeSubcellDivergence => SubcellDG_ComputeSubcellDivergence
   end type SubcellDG_t
!
!  ========
   contains
!  ========
!
      subroutine SubcellDG_Initialize(self, controlVariables)
!
!        -------
!        Modules
!        -------
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
!
!        ---------
!        Interface
!        ---------
         class(SubcellDG_t)                   :: self
         class(FTValueDictionary), intent(in) :: controlVariables

!
!        Set the maximum blending coefficient (0.5 by default)
!        -----------------------------------------------------
         if (controlVariables % containsKey(maxBlendKey)) then
            self % max_blend = controlVariables % realValueForKey(maxBlendKey)
         else
            self % max_blend = 0.5_RP
         end if
!
!        Setup the Riemann solver and two-point flux
!        -------------------------------------------
         call SetRiemannSolver(controlVariables)
!
!        Describe
!        --------
         if (.not. MPI_Process % isRoot ) return

         call Subsection_Header("Hyperbolic discretization")

         write(STD_OUT,'(30X,A,A30,A)') "->", "Numerical scheme: ", "Split-Form (subcell)"
         write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Maximum DGFV blending:", self % max_blend

         call DescribeRiemannSolver

         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->", "Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->", "Gradients computation: ", "Disabled."
         end if

      end subroutine SubcellDG_Initialize

      subroutine SubcellDG_ComputeSplitFormFluxes(self, e, contravariantFlux, fSharp, gSharp, hSharp)
!
!        -------
!        Modules
!        -------
         use ElementClass
         use PhysicsStorage
!
!        ---------
!        Interface
!        ---------
         class(SubcellDG_t), intent(in)  :: self
         type(Element),      intent(in)  :: e
         real(kind=RP),      intent(in)  :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP),      intent(out) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3))
         real(kind=RP),      intent(out) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3))
         real(kind=RP),      intent(out) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
         integer     :: i, j, k, l

         associate ( Q => e % storage % Q )
!
!        First, diagonal results are introduced directly (consistency property)
!        ----------------------------------------------------------------------
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            fSharp(:,i,i,j,k) = contravariantFlux(:,i,j,k,IX)
            gSharp(:,j,i,j,k) = contravariantFlux(:,i,j,k,IY)
            hSharp(:,k,i,j,k) = contravariantFlux(:,i,j,k,IZ)
         end do               ; end do             ; end do
!
!        Then, terms out of the diagonal are computed
!        --------------------------------------------
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = i+1, e%Nxyz(1)
               call TwoPointFlux(Q(:,i,j,k), Q(:,l,j,k), e % geom % jGradXi(:,i,j,k), e % geom % jGradXi(:,l,j,k), fSharp(:,l,i,j,k))
               fSharp(:,i,l,j,k) = fSharp(:,l,i,j,k)
            end do
         end do               ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = j+1, e%Nxyz(2)
               call TwoPointFlux(Q(:,i,j,k), Q(:,i,l,k), e % geom % jGradEta(:,i,j,k), e % geom % jGradEta(:,i,l,k), gSharp(:,l,i,j,k))
               gSharp(:,j,i,l,k) = gSharp(:,l,i,j,k)
            end do
         end do               ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = k+1, e%Nxyz(3)
               call TwoPointFlux(Q(:,i,j,k), Q(:,i,j,l), e % geom % jGradZeta(:,i,j,k), e % geom % jGradZeta(:,i,j,l), hSharp(:,l,i,j,k))
               hSharp(:,k,i,j,l) = hSharp(:,l,i,j,k)
            end do
         end do               ; end do             ; end do

         end associate

      end subroutine SubcellDG_ComputeSplitFormFluxes

      function SubcellDG_ComputeSubcellDivergence(self, e, fSharp, gSharp, hSharp, Fv) result (div)
!
!        -------
!        Modules
!        -------
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
!
!        ---------
!        Interface
!        ---------
         class(SubcellDG_t), intent(in) :: self
         type(Element),      intent(in) :: e
         real(kind=RP),      intent(in) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),      intent(in) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),      intent(in) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),      intent(in) :: Fv(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP)                  :: div(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
         integer  :: i, j, k, l
         real(RP) :: blending
         real(RP) :: Fcx(NCONS, 0:e % Nxyz(1)+1)
         real(RP) :: Fcy(NCONS, 0:e % Nxyz(2)+1)
         real(RP) :: Fcz(NCONS, 0:e % Nxyz(3)+1)

!
!        Blending between DG and FV fluxes
!        ---------------------------------
         blending = min(e % storage % sensor, self % max_blend)

         associate(Nx => e % Nxyz(1), &
                   Ny => e % Nxyz(2), &
                   Nz => e % Nxyz(3)  )
         associate(spAxi   => NodalStorage(Nx), &
                   spAeta  => NodalStorage(Ny), &
                   spAzeta => NodalStorage(Nz)  )
!
!        ------------
!        Direction Xi
!        ------------
         do k = 0, Nz ; do j = 0, Ny
!
!           Hyperbolic term
!           ---------------
            call ComputeSubcellFluxes_1D(e, blending, e % storage % Q(:,:,j,k),                &
                                         spAxi % w, spAxi % D, e % geom % ncXi(:,:,j,k),       &
                                         e % geom % t1cXi(:,:,j,k), e % geom % t2cXi(:,:,j,k), &
                                         e % geom % JfcXi(:,j,k), fSharp(:,:,:,j,k), Fcx)
            do i = 0, Nx
               ! Note that we initialize div here
               div(:,i,j,k) = (Fcx(:,i) - Fcx(:,i+1)) / spAxi % w(i)
            end do
!
!           Viscous term
!           ------------
            do l = 0, Nx ; do i = 0, Nx
               div(:,i,j,k) = div(:,i,j,k) + spAxi % hatD(i,l) * Fv(:,l,j,k,IX)
            end do       ; end do
         end do       ; end do
!
!        -------------
!        Direction Eta
!        -------------
         do k = 0, Nz ; do i = 0, Nx
!
!           Hyperbolic term
!           ---------------
            call ComputeSubcellFluxes_1D(e, blending, e % storage % Q(:,i,:,k),                  &
                                         spAeta % w, spAeta % D, e % geom % ncEta(:,i,:,k),      &
                                         e % geom % t1cEta(:,i,:,k), e % geom % t2cEta(:,i,:,k), &
                                         e % geom % JfcEta(i,:,k), gSharp(:,:,i,:,k), Fcy)
            do j = 0, Ny
               div(:,i,j,k) = div(:,i,j,k) + (Fcy(:,j) - Fcy(:,j+1)) / spAeta % w(j)
            end do
!
!           Viscous term
!           ------------
            do l = 0, Ny ; do j = 0, Ny
               div(:,i,j,k) = div(:,i,j,k) + spAeta % hatD(j,l) * Fv(:,i,l,k,IY)
            end do       ; end do
         end do       ; end do
!
!        --------------
!        Direction Zeta
!        --------------
         do j = 0, Ny ; do i = 0, Nx
!
!           Hyperbolic term
!           ---------------
            call ComputeSubcellFluxes_1D(e, blending, e % storage % Q(:,i,j,:),                    &
                                         spAzeta % w, spAzeta % D, e % geom % ncZeta(:,i,j,:),     &
                                         e % geom % t1cZeta(:,i,j,:), e % geom % t2cZeta(:,i,j,:), &
                                         e % geom % JfcZeta(i,j,:), hSharp(:,:,i,j,:), Fcz)
            do k = 0, Nz
               div(:,i,j,k) = div(:,i,j,k) + (Fcz(:,k) - Fcz(:,k+1)) / spAzeta % w(k)
            end do
!
!           Viscous term
!           ------------
            do l = 0, Nz ; do k = 0, Nz
               div(:,i,j,k) = div(:,i,j,k) + spAzeta % hatD(k,l) * Fv(:,i,j,l,IZ)
            end do       ; end do
         end do       ; end do

         end associate
         end associate

      end function SubcellDG_ComputeSubcellDivergence

      subroutine ComputeSubcellFluxes_1D(e, blending, Q, w, D, n, t1, t2, Jf, Fs, Fc)
!
!        -------
!        Modules
!        -------
         use ElementClass
!
!        ---------
!        Interface
!        ---------
         type(Element), intent(in)  :: e
         real(RP),      intent(in)  :: blending
         real(RP),      intent(in)  :: Q(1:, 0:)       ! (NCONS, 0:N)
         real(RP),      intent(in)  :: w(0:)           ! (0:N)
         real(RP),      intent(in)  :: D(0:, 0:)       ! (0:N, 0:N)
         real(RP),      intent(in)  :: n(1:, 0:)       ! (NDIMS, 0:N+1)
         real(RP),      intent(in)  :: t1(1:, 0:)      ! (NDIMS, 0:N+1)
         real(RP),      intent(in)  :: t2(1:, 0:)      ! (NDIMS, 0:N+1)
         real(RP),      intent(in)  :: Jf(0:)          ! (0:N+1)
         real(RP),      intent(in)  :: Fs(1:, 0:, 0:)  ! (NCONS, 0:N, 0:N)
         real(RP),      intent(out) :: Fc(1:, 0:)      ! (NCONS, 0:N+1)
!
!        ---------------
!        Local variables
!        ---------------
         integer  :: nx, nc
         integer  :: i
         integer  :: r, s
         real(RP) :: ws, wv
         real(RP) :: Fv(1:size(Fc, 1))


         nx = size(Fs, 2) - 1
         nc = size(Fc, 2) - 1  ! n + 1
!
!        Boundary fluxes are zero here, its value is set in the surface subroutine
!        -------------------------------------------------------------------------
         Fc(:,0)  = 0.0_RP
         Fc(:,nc) = 0.0_RP
!
!        Compute subcell, DG flux
!        ------------------------
         do i = 1, nc-1
            Fc(:,i) = 0.0_RP
            do s = i, nx ; do r = 0, i-1
               Fc(:,i) = Fc(:,i) + w(r) * D(r,s) * Fs(:,r,s)
            end do                  ; end do
            Fc(:,i) = 2.0_RP * Fc(:,i)
         end do
!
!        Compute subcell, FV flux
!        ------------------------
         if (blending > 0.0_RP) then
            ws = 1.0_RP - blending
            wv = blending
            do i = 1, nc-1
               call RiemannSolver(Q(:,i-1), Q(:,i), n(:,i), t1(:,i), t2(:,i), Fv)
               Fv = Fv * Jf(i)
               Fc(:,i) = ws * Fc(:,i) + wv * Fv
            end do
         end if

      end subroutine ComputeSubcellFluxes_1D

end module HyperbolicSubcell
#endif
