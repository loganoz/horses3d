!
!//////////////////////////////////////////////////////
!
!   @File:    SVV.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sat Jan  6 11:47:48 2018
!   @Last revision date: Tue Feb 26 18:17:37 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 2bf19fe8b635bd13e65827bd1023dd656a765b42
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
module SpectralVanishingViscosity
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use HexMeshClass
   use NodalStorageClass
   use GaussQuadrature
   use FluidData
   use MPI_Process_Info             , only: MPI_Process
   use Headers                      , only: Subsection_Header
   use LESModels                    , only: Smagorinsky_t
   use Utilities                    , only: toLower
   implicit none

   private
   public   SVV, InitializeSVV

   integer,          parameter  :: Nmax = 20
   character(len=*), parameter  :: SVV_KEY        = "enable svv"
   character(len=*), parameter  :: SVV_MU_KEY     = "svv viscosity"
   character(len=*), parameter  :: SVV_CUTOFF_KEY = "svv filter cutoff"

   type FilterMatrices_t
      logical                    :: constructed = .false.
      integer                    :: N
      real(kind=RP), allocatable :: Q(:,:)
   end type FilterMatrices_t

   type  SVV_t
      logical                :: enabled
      logical                :: muIsSmagorinsky = .FALSE.
      real(kind=RP)          :: muSVV
      real(kind=RP)          :: Psvv
      type(FilterMatrices_t) :: filters(0:Nmax)
      contains
         procedure      :: ConstructFilter    => SVV_ConstructFilter
         procedure      :: ComputeInnerFluxes => SVV_ComputeInnerFluxes
         procedure      :: RiemannSolver      => SVV_RiemannSolver
         procedure      :: describe           => SVV_Describe
         procedure      :: destruct           => SVV_destruct
   end type SVV_t

   type(SVV_t), protected    :: SVV
   type(Smagorinsky_t)       :: Smagorinsky
!
!  ========
   contains
!  ========
!
      subroutine InitializeSVV(self, controlVariables, mesh)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use PhysicsStorage
         implicit none
         class(SVV_t)                          :: self
         class(FTValueDictionary),  intent(in) :: controlVariables
         class(HexMesh),            intent(in) :: mesh
!
!        ---------------
!        Local variables         
!        ---------------
!
         integer     :: eID
         character(len=LINE_LENGTH) :: mu
!
!        -------------------------
!        Check if SVV is requested
!        -------------------------
!
         if ( controlVariables % containsKey(SVV_KEY) ) then
            if ( controlVariables % logicalValueForKey(SVV_KEY) ) then
               self % enabled = .true.

            else
               self % enabled = .false.

            end if

         else
            self % enabled = .false.
   
         end if

         if ( .not. self % enabled ) return
!
!        ---------------------
!        Get the SVV viscosity: the viscosity is later multiplied by 1/N
!        ---------------------
!
         if ( controlVariables % containsKey(SVV_MU_KEY) ) then
            
            mu = trim(controlVariables % StringValueForKey(SVV_MU_KEY,LINE_LENGTH) )
            call ToLower(mu)
            
            select case ( mu )
               case ('smagorinsky')
                  self % muIsSmagorinsky = .TRUE.
                  call Smagorinsky % Initialize (controlVariables)
               case default
                  self % muSVV = controlVariables % doublePrecisionValueForKey(SVV_MU_KEY)
            end select

         else
            self % muSVV = 0.1_RP

         end if 
!
!        -------------------------
!        Get the SVV kernel cutoff: the cutoff is later multiplied by sqrt(N)
!        -------------------------
!
         if ( controlVariables % containsKey(SVV_CUTOFF_KEY) ) then
            
            
            self % Psvv = controlVariables % doublePrecisionValueForKey(SVV_CUTOFF_KEY)

         else
            self % Psvv = 1.0_RP

         end if
!
!        Display the configuration
!        -------------------------
         call self % describe
!
!        Construct the filters
!        ---------------------
         do eID = 1, mesh % no_of_elements
            associate(Nxyz => mesh % elements(eID) % Nxyz)
            call self % ConstructFilter(Nxyz(1))
            call self % ConstructFilter(Nxyz(2))
            call self % ConstructFilter(Nxyz(3))
            end associate
         end do

      end subroutine InitializeSVV
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_Describe(this)
         implicit none
         !-arguments----------------------------------------
         class(SVV_t), intent (in) :: this
         !--------------------------------------------------
         
         if (.not. MPI_Process % isRoot) return
         
         write(STD_OUT,'(/)')
         call Subsection_Header("Spectral vanishing viscosity (SVV)")
         
         if (this % muIsSmagorinsky) then
            write(STD_OUT,'(30X,A,A30,A,F4.2,A)') "->","viscosity: ", "Smagorinsky (Cs = ", Smagorinsky % CS,  ")"
         else
            write(STD_OUT,'(30X,A,A30,F10.3)') "->","viscosity: ", this % muSVV
         end if
         write(STD_OUT,'(30X,A,A30,F10.3)') "->","filter cutoff: ", this % Psvv
         
      end subroutine SVV_Describe
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_ComputeInnerFluxes( self , e , EllipticFlux, contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         use LESModels
         implicit none
         class(SVV_t) ,     intent (in)         :: self
         type(Element)                          :: e
         procedure(EllipticFlux3D_f)            :: EllipticFlux
         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: Uxf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Uyf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Uzf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: cartesianFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k, ii, jj, kk
         real(kind=RP)       :: Q3D, delta
!
!        -------------------------
!        Compute the SVV viscosity
!        -------------------------
!
         if (self % muIsSmagorinsky) then
!
!           (1+Psvv) * muSmag
!           -----------------
            delta = (e % geom % Volume / product(e % Nxyz + 1) ) ** (1._RP / 3._RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call Smagorinsky % ComputeViscosity ( delta, e % geom % dWall(i,j,k), e % storage % Q  (:,i,j,k) &
                                                                                   , e % storage % U_x(:,i,j,k) &
                                                                                   , e % storage % U_y(:,i,j,k) &
                                                                                   , e % storage % U_z(:,i,j,k) &
                                                                                   , mu(i,j,k) )
               mu(i,j,k) = mu(i,j,k) * (1._RP + self % Psvv)
            end do                ; end do                ; end do
         else
!
!           Fixed value
!           -----------
            mu    = self % muSVV !/ maxval(e % Nxyz+1)
         end if

         beta  = 0.0_RP
         kappa = mu / ( thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Pr)
!
!        --------------------
!        Filter the gradients
!        --------------------
!
         associate(Qx => self % filters(e % Nxyz(1)) % Q, & 
                   Qy => self % filters(e % Nxyz(2)) % Q, &
                   Qz => self % filters(e % Nxyz(3)) % Q    )

         Uxf = 0.0_RP   ; Uyf = 0.0_RP    ; Uzf = 0.0_RP
         do kk = 0, e % Nxyz(3)  ; do jj = 0, e % Nxyz(2)   ; do ii = 0, e % Nxyz(1)
            do k = 0, e % Nxyz(3)  ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
               Q3D = Qx(ii,i) * Qy(jj,j) * Qz(kk,k)
               Uxf(:,ii,jj,kk) = Uxf(:,ii,jj,kk) + Q3D * e % storage % U_x(:,i,j,k)
               Uyf(:,ii,jj,kk) = Uyf(:,ii,jj,kk) + Q3D * e % storage % U_y(:,i,j,k)
               Uzf(:,ii,jj,kk) = Uzf(:,ii,jj,kk) + Q3D * e % storage % U_z(:,i,j,k)
            end do                 ; end do                  ; end do
         end do                  ; end do                   ; end do

         end associate

         call EllipticFlux( NCONS, NGRAD, e%Nxyz, e % storage % Q , Uxf, Uyf, Uzf, mu, beta, kappa, cartesianFlux )

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,i,j,k,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,i,j,k,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,i,j,k,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do            ; end do

      end subroutine SVV_ComputeInnerFluxes

      subroutine SVV_RiemannSolver ( self, f, EllipticFlux, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, flux)
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         use LESModels
         implicit none
         class(SVV_t)                :: self
         class(Face),   intent(in)   :: f
         procedure(EllipticFlux2D_f) :: EllipticFlux
         real(kind=RP), intent(in)   :: QLeft(NCONS, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: QRight (NCONS, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: U_xLeft(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: U_yLeft(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: U_zLeft(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: U_xRight(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: U_yRight(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(in)   :: U_zRight(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP), intent(out)  :: flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer           :: i, j, ii, jj
         real(kind=RP)     :: Q(NCONS, 0:f % Nf(1), 0:f % Nf(2)) 
         real(kind=RP)     :: U_x(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: U_y(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: U_z(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: Uxf(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: Uyf(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: Uzf(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: flux_vec(NCONS,NDIM, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: mu(0:f % Nf(1), 0:f % Nf(2)), kappa(0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: beta(0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)     :: delta, Q2D
!
!        Interface averages
!        ------------------
         Q   = 0.5_RP * ( QLeft + QRight)
         U_x = 0.5_RP * ( U_xLeft + U_xRight)
         U_y = 0.5_RP * ( U_yLeft + U_yRight)
         U_z = 0.5_RP * ( U_zLeft + U_zRight)
!
!        -------------------------
!        Compute the SVV viscosity
!        -------------------------
!
         if (self % muIsSmagorinsky) then
!
!           (1+Psvv) * muSmag
!           -----------------
            delta = sqrt(f % geom % surface / product(f % Nf + 1) )
            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
               call Smagorinsky % ComputeViscosity ( delta, f % geom % dWall(i,j), Q  (:,i,j) &
                                                                                 , U_x(:,i,j) &
                                                                                 , U_y(:,i,j) &
                                                                                 , U_z(:,i,j) &
                                                                                 , mu(i,j) )
               mu(i,j) = mu(i,j) * (1._RP + self % Psvv)
            end do                ; end do
         else
!
!           Fixed value
!           -----------
            mu    = self % muSVV !/ maxval(f % Nf+1)
         end if
         
         beta  = 0.0_RP
         kappa = mu / ( thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Pr)
!
!        --------------------
!        Filter the gradients
!        --------------------
!
         associate(Qx => self % filters(f % Nf(1)) % Q, & 
                   Qy => self % filters(f % Nf(2)) % Q   )

         Uxf = 0.0_RP   ; Uyf = 0.0_RP    ; Uzf = 0.0_RP
         do jj = 0, f % Nf(2)   ; do ii = 0, f % Nf(1)
            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
               Q2D = Qx(ii,i) * Qy(jj,j) 
               Uxf(:,ii,jj) = Uxf(:,ii,jj) + Q2D * U_x(:,i,j)
               Uyf(:,ii,jj) = Uyf(:,ii,jj) + Q2D * U_y(:,i,j)
               Uzf(:,ii,jj) = Uzf(:,ii,jj) + Q2D * U_z(:,i,j)
            end do                  ; end do
         end do                   ; end do

         end associate

         call EllipticFlux(NCONS, NGRAD, f % Nf, Q,U_x,U_y,U_z, mu, beta, kappa, flux_vec)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            flux(:,i,j) =   flux_vec(:,IX,i,j) * f % geom % normal(IX,i,j) &
                          + flux_vec(:,IY,i,j) * f % geom % normal(IY,i,j) &
                          + flux_vec(:,IZ,i,j) * f % geom % normal(IZ,i,j) 
         end do               ; end do

      end subroutine SVV_RiemannSolver
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           Auxiliar subroutines
!           --------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_constructFilter(self, N)
         implicit none
         class(SVV_t)         :: self
         integer, intent(in) :: N
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: Nodal2Modal(0:N,0:N)
         real(kind=RP)  :: Modal2Nodal(0:N,0:N)
         real(kind=RP)  :: filterCoefficients(0:N)
         real(kind=RP)  :: Lkj(0:N,0:N), dLk_dummy
         real(kind=RP)  :: normLk(0:N)
         real(kind=RP)  :: filterExp

         if ( self % filters(N) % Constructed ) return
!
!        Get the evaluation of Legendre polynomials at the interpolation nodes
!        ---------------------------------------------------------------------
         do j = 0 , N ;    do k = 0 , N
            call LegendrePolyAndDerivative(k, NodalStorage(N) % x(j), Lkj(k,j), dLk_dummy)
         end do       ;    end do
!
!        Get the norm of Legendre polynomials
!        ------------------------------------
         normLk = 0.0_RP
         do k = 0 , N   ; do j = 0 , N
            normLk(k) = normLk(k) + NodalStorage(N) % w(j) * Lkj(k,j) * Lkj(k,j)
         end do         ; end do
!
!        Get the transformation from Nodal to Modal and viceversa matrices
!        -----------------------------------------------------------------
         do k = 0 , N   ; do i = 0 , N
            Nodal2Modal(k,i) = NodalStorage(N) % w(i) * Lkj(k,i) / normLk(k)
            Modal2Nodal(i,k) = Lkj(k,i)
         end do         ; end do
!
!        Get the filter coefficients
!        ---------------------------
         filterExp = self % Psvv !* sqrt( real(N, kind=RP) )

         do k = 0, N
            filterCoefficients(k) = (real(k, kind=RP) / N + 1.0e-12_RP) ** filterExp
         end do
!
!        Compute the filtering matrix
!        ----------------------------
         self % filters(N) % N = N
         allocate(self % filters(N) % Q(0:N,0:N))

         self % filters(N) % Q = 0.0_RP
         do k = 0, N ; do j = 0, N  ; do i = 0, N
            self % filters(N) % Q(i,j) = self % filters(N) % Q(i,j) + Modal2Nodal(i,k) * filterCoefficients(k) * Nodal2Modal(k,j)
         end do      ; end do       ; end do

         self % filters(N) % constructed = .true.

      end subroutine SVV_constructFilter
      
      subroutine SVV_destruct(this)
         implicit none
         class(SVV_t) :: this
         integer :: i
         
         do i = 0, Nmax
            if ( this % filters(i) % constructed ) deallocate(this % filters(i) % Q)
         end do
         
         !if (this % muIsSmagorinsky) call Smagorinsky % destruct
         
      end subroutine SVV_destruct
      
end module SpectralVanishingViscosity
#endif
