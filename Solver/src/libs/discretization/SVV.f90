!
!//////////////////////////////////////////////////////
!
!   @File:    SVV.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sat Jan  6 11:47:48 2018
!   @Last revision date: Tue Mar 12 15:50:35 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: e91212aadaa211fa6b91bd7bc18009c1e482533d
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
!
!  Keywords
!  --------
   character(len=*), parameter  :: SVV_KEY          = "enable svv"
   character(len=*), parameter  :: SVV_MU_KEY       = "svv viscosity"
   character(len=*), parameter  :: SVV_CUTOFF_KEY   = "svv filter cutoff"
   character(len=*), parameter  :: FILTER_SHAPE_KEY = "svv filter shape"
   character(len=*), parameter  :: FILTER_TYPE_KEY  = "svv filter type"
   character(len=*), parameter  :: POST_FILTER_KEY  = "svv post filter"
!
!  Filter types
!  ------------
   enum, bind(C)
      enumerator :: HPASS_FILTER, LPASS_FILTER
   end enum
!
!  Filter shapes
!  -------------
   enum, bind(C)
      enumerator :: POW_FILTER, SHARP_FILTER, EXP_FILTER
   end enum
   
   type FilterMatrices_t
      logical                    :: constructed = .false.
      integer                    :: N
      real(kind=RP), allocatable :: Q(:,:)
   end type FilterMatrices_t

   type  SVV_t
      logical                :: enabled
      logical                :: muIsSmagorinsky = .FALSE.
      logical                :: postFiltering   = .FALSE.
      integer                :: filterType
      integer                :: filterShape
      real(kind=RP)          :: muSVV
      real(kind=RP)          :: Psvv
      type(FilterMatrices_t) :: filters(0:Nmax)
      contains
         procedure      :: ConstructFilter    => SVV_ConstructFilter
         procedure      :: ComputeInnerFluxes => SVV_ComputeInnerFluxes
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
!        ------------
!        Filter order
!        ------------
!
         if ( controlVariables % containsKey(POST_FILTER_KEY) ) then
            self % postFiltering = controlVariables % logicalValueForKey(POST_FILTER_KEY)
         end if
!
!        --------------
!        Type of filter
!        --------------
!
         if ( controlVariables % containsKey(FILTER_TYPE_KEY) ) then
            select case ( trim(controlVariables % stringValueForKey(FILTER_TYPE_KEY,LINE_LENGTH)) )
               case ("low-pass" ) ; self % filterType = LPASS_FILTER
               case ("high-pass") ; self % filterType = HPASS_FILTER
               case default
                  write(STD_OUT,*) 'ERROR. SVV filter type not recognized. Options are:'
                  write(STD_OUT,*) '   * low-pass'
                  write(STD_OUT,*) '   * high-pass'
                  stop
            end select
         else
            self % filterType = HPASS_FILTER
         end if
!
!        ---------------
!        Shape of filter
!        ---------------
!
         if ( controlVariables % containsKey(FILTER_SHAPE_KEY) ) then
            select case ( trim(controlVariables % stringValueForKey(FILTER_SHAPE_KEY,LINE_LENGTH)) )
               case ("power")       ; self % filterShape = POW_FILTER
               case ("sharp")       ; self % filterShape = SHARP_FILTER
               case ("exponential") ; self % filterShape = EXP_FILTER
               case default
                  write(STD_OUT,*) 'ERROR. SVV filter shape not recognized. Options are:'
                  write(STD_OUT,*) '   * power'
                  write(STD_OUT,*) '   * sharp'
                  write(STD_OUT,*) '   * exponential'
                  stop
            end select
         else
            self % filterShape = POW_FILTER
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
         
         write(STD_OUT,'(30X,A,A30,L)') "->","post-filtering: ", this % postFiltering
         
         select case (this % filterType)
            case (LPASS_FILTER)
               write(STD_OUT,'(30X,A,A30,A)') "->","filter type: ", 'low-pass'
            case (HPASS_FILTER)
               write(STD_OUT,'(30X,A,A30,A)') "->","filter type: ", 'high-pass'
         end select
         
         select case (this % filterShape)
            case (POW_FILTER)
               write(STD_OUT,'(30X,A,A30,A)') "->","filter shape: ", 'power kernel'
            case (EXP_FILTER)
               write(STD_OUT,'(30X,A,A30,A)') "->","filter shape: ", 'exponential kernel'
            case (SHARP_FILTER)
               write(STD_OUT,'(30X,A,A30,A)') "->","filter shape: ", 'sharp kernel'
         end select
         
      end subroutine SVV_Describe
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_ComputeInnerFluxes(self, mesh, e, contravariantFlux)
         use ElementClass
         use PhysicsStorage
         use Physics
         use LESModels
         implicit none
         class(SVV_t)                           :: self
         type(HexMesh)                          :: mesh
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer             :: i, j, k, l, fIDs(6)
         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hx(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hy(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hz(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hxf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hyf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hzf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hxf_aux(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hyf_aux(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hzf_aux(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: cartesianFlux(1:NCONS, 1:NDIM)
         real(kind=RP)       :: SVV_diss

         mu    = self % muSVV

         associate(Qx => self % filters(e % Nxyz(1)) % Q, & 
                   Qy => self % filters(e % Nxyz(2)) % Q, &
                   Qz => self % filters(e % Nxyz(3)) % Q    )

         associate(spA_xi   => NodalStorage(e % Nxyz(1)), &
                   spA_eta  => NodalStorage(e % Nxyz(2)), &
                   spA_zeta => NodalStorage(e % Nxyz(3))) 

!
!        -----------------
!        Compute the Hflux
!        -----------------
!
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
           Hx(:,i,j,k) = sqrt(e % geom % jacobian(i,j,k)*mu(i,j,k)) * e % storage % U_x(:,i,j,k)
           Hy(:,i,j,k) = sqrt(e % geom % jacobian(i,j,k)*mu(i,j,k)) * e % storage % U_y(:,i,j,k)
           Hz(:,i,j,k) = sqrt(e % geom % jacobian(i,j,k)*mu(i,j,k)) * e % storage % U_z(:,i,j,k)
         end do                ; end do                ; end do
!
!        ----------------
!        Filter the Hflux
!        ----------------
!
!        Perform filtering in xi Hf_aux -> Hf
!        -----------------------
         Hxf_aux = Hx     ; Hyf_aux = Hy     ; Hzf_aux = Hz
         Hxf     = 0.0_RP ; Hyf     = 0.0_RP ; Hzf     = 0.0_RP
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            do l = 0, e % Nxyz(1)
               Hxf(:,i,j,k) = Hxf(:,i,j,k) + Qx(i,l) * Hxf_aux(:,l,j,k)
               Hyf(:,i,j,k) = Hyf(:,i,j,k) + Qx(i,l) * Hyf_aux(:,l,j,k)
               Hzf(:,i,j,k) = Hzf(:,i,j,k) + Qx(i,l) * Hzf_aux(:,l,j,k)
            end do
         end do                ; end do                ; end do
!
!        Perform filtering in eta Hf -> Hf_aux
!        ------------------------
         Hxf_aux = 0.0_RP  ; Hyf_aux = 0.0_RP  ; Hzf_aux = 0.0_RP
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do l = 0, e % Nxyz(1) ; do i = 0, e % Nxyz(1)
            Hxf_aux(:,i,j,k) = Hxf_aux(:,i,j,k) + Qy(j,l) * Hxf(:,i,l,k)
            Hyf_aux(:,i,j,k) = Hyf_aux(:,i,j,k) + Qy(j,l) * Hyf(:,i,l,k)
            Hzf_aux(:,i,j,k) = Hzf_aux(:,i,j,k) + Qy(j,l) * Hzf(:,i,l,k)
         end do                ; end do                ; end do                ; end do
!
!        Perform filtering in zeta Hf_aux -> Hf
!        -------------------------
         Hxf = 0.0_RP  ; Hyf = 0.0_RP  ; Hzf = 0.0_RP
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do l = 0, e % Nxyz(1) ; do i = 0, e % Nxyz(1)
            Hxf(:,i,j,k) = Hxf(:,i,j,k) + Qz(k,l) * Hxf_aux(:,i,j,l)
            Hyf(:,i,j,k) = Hyf(:,i,j,k) + Qz(k,l) * Hyf_aux(:,i,j,l)
            Hzf(:,i,j,k) = Hzf(:,i,j,k) + Qz(k,l) * Hzf_aux(:,i,j,l)
         end do                ; end do                ; end do                ; end do
         
         if (self % filterType == LPASS_FILTER) then
            Hxf = Hx - Hxf
            Hyf = Hy - Hyf
            Hzf = Hz - Hzf
         end if
!
!        ----------------
!        Compute the flux
!        ----------------
!
         SVV_diss = 0.0_RP
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call Hflux_physical_dissipation_ENERGY( NCONS, NGRAD, e % storage % Q(:,i,j,k), Hxf(:,i,j,k), Hyf(:,i,j,k), &
                               Hzf(:,i,j,k), 1.0_RP, 0.0_RP, cartesianFlux )
   
            SVV_diss = SVV_diss + spA_xi % w(i) * spA_eta % w(j) * spA_zeta % w(k) * (sum(cartesianFlux(2:4,IX)*Hx(2:4,i,j,k) + &
                                                                                          cartesianFlux(2:4,IY)*Hy(2:4,i,j,k) + & 
                                                                                          cartesianFlux(2:4,IZ)*Hz(2:4,i,j,k)))

            cartesianFlux = sqrt(mu(i,j,k) * e % geom % invJacobian(i,j,k)) * cartesianFlux
            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do            ; end do

         e % storage % SVV_diss = SVV_diss
!
!        --------------------
!        Prolong to the faces
!        --------------------
!
         fIDs = e % faceIDs
         call e % ProlongHfluxToFaces(NCONS, contravariantFlux, &
                                   mesh % faces(fIDs(1)),&
                                   mesh % faces(fIDs(2)),&
                                   mesh % faces(fIDs(3)),&
                                   mesh % faces(fIDs(4)),&
                                   mesh % faces(fIDs(5)),&
                                   mesh % faces(fIDs(6))   )


         end associate
         end associate
         
      end subroutine SVV_ComputeInnerFluxes
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        Library with H fluxes
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Hflux_diagonal_dissipation(NCONS, NGRAD, Q, Hx, Hy, Hz, mu, alpha, F)
         implicit none
         integer, intent(in)        :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Hx(NGRAD)
         real(kind=RP), intent(in)  :: Hy(NGRAD)
         real(kind=RP), intent(in)  :: Hz(NGRAD)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: alpha
         real(kind=RP), intent(out) :: F(NCONS, NDIM)

         F(:,IX) = Hx
         F(:,IY) = Hy
         F(:,IZ) = Hz

      end subroutine Hflux_diagonal_dissipation

      subroutine Hflux_physical_dissipation_ENERGY(NCONS, NGRAD, Q, Hx, Hy, Hz, mu, alpha, F)
         implicit none
         integer, intent(in)        :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Hx(NGRAD)
         real(kind=RP), intent(in)  :: Hy(NGRAD)
         real(kind=RP), intent(in)  :: Hz(NGRAD)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: alpha
         real(kind=RP), intent(out) :: F(NCONS, NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invRho, divV, u(NDIM)
         real(kind=RP)  :: kappa

         kappa = mu * dimensionless % kappa * dimensionless % Re 
  
         invRho  = 1.0_RP / Q(IRHO)
         u = Q(IRHOU:IRHOW)*invRho

         divV = Hx(IX) + Hy(IY) + Hz(IZ)

         F(IRHO,IX)  = 0.0_RP
         F(IRHOU,IX) = mu * (2.0_RP * Hx(IRHOU) - 2.0_RP/3.0_RP * divV ) 
         F(IRHOV,IX) = mu * ( Hx(IRHOV) + Hy(IRHOU) ) 
         F(IRHOW,IX) = mu * ( Hx(IRHOW) + Hz(IRHOU) ) 
         F(IRHOE,IX) = F(IRHOU,IX) * u(IRHOU) + F(IRHOV,IX) * u(IRHOV) + F(IRHOW,IX) * u(IRHOW) + kappa * Hx(IRHOE) 

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX) 
         F(IRHOV,IY) = mu * (2.0_RP * Hy(IRHOV) - 2.0_RP / 3.0_RP * divV )
         F(IRHOW,IY) = mu * ( Hy(IRHOW) + Hz(IRHOV) ) 
         F(IRHOE,IY) = F(IRHOU,IY) * u(IRHOU) + F(IRHOV,IY) * u(IRHOV) + F(IRHOW,IY) * u(IRHOW) + kappa * Hy(IRHOE)

         F(IRHO,IZ) = 0.0_RP
         F(IRHOU,IZ) = F(IRHOW,IX) 
         F(IRHOV,IZ) = F(IRHOW,IY) 
         F(IRHOW,IZ) = mu * ( 2.0_RP * Hz(IRHOW) - 2.0_RP / 3.0_RP * divV ) 
         F(IRHOE,IZ) = F(IRHOU,IZ) * u(IRHOU) + F(IRHOV,IZ) * u(IRHOV) + F(IRHOW,IZ) * u(IRHOW) + kappa * Hz(IRHOE)

      end subroutine Hflux_physical_dissipation_ENERGY

!      subroutine SVV_ComputeInnerFluxes( self , e , EllipticFlux, contravariantFlux )
!         use ElementClass
!         use PhysicsStorage
!         use Physics
!         use LESModels
!         implicit none
!         class(SVV_t) ,     intent (in)         :: self
!         type(Element)                          :: e
!         procedure(EllipticFlux_f)              :: EllipticFlux
!         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!!
!!        ---------------
!!        Local variables
!!        ---------------
!!
!         real(kind=RP)       :: Uxf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!         real(kind=RP)       :: Uyf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!         real(kind=RP)       :: Uzf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!         real(kind=RP)       :: cartesianFlux      (1:NCONS, 1:NDIM)
!         real(kind=RP)       :: contravariantFluxF (1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!         real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!         integer             :: i, j, k, ii, jj, kk
!         real(kind=RP)       :: Q3D, delta
!!
!!        -------------------------
!!        Compute the SVV viscosity
!!        -------------------------
!!
!         if (self % muIsSmagorinsky) then
!!
!!           (1+Psvv) * muSmag
!!           -----------------
!            delta = (e % geom % Volume / product(e % Nxyz + 1) ) ** (1._RP / 3._RP)
!            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!               call Smagorinsky % ComputeViscosity ( delta, e % geom % dWall(i,j,k), e % storage % Q  (:,i,j,k) &
!                                                                                   , e % storage % U_x(:,i,j,k) &
!                                                                                   , e % storage % U_y(:,i,j,k) &
!                                                                                   , e % storage % U_z(:,i,j,k) &
!                                                                                   , mu(i,j,k) )
!               mu(i,j,k) = mu(i,j,k)! * (1._RP + self % Psvv)
!            end do                ; end do                ; end do
!         else
!!
!!           Fixed value
!!           -----------
!            mu    = self % muSVV !/ maxval(e % Nxyz+1)
!         end if
!
!         beta  = 0.0_RP
!         kappa = mu / ( thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Prt)
!!
!!        --------------------
!!        Filter the gradients
!!        --------------------
!!
!         if (.not. self % postFiltering) then
!            associate(Qx => self % filters(e % Nxyz(1)) % Q, & 
!                      Qy => self % filters(e % Nxyz(2)) % Q, &
!                      Qz => self % filters(e % Nxyz(3)) % Q    )
!
!            Uxf = 0.0_RP   ; Uyf = 0.0_RP    ; Uzf = 0.0_RP
!            do k = 0, e % Nxyz(3)  ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
!               do kk = 0, e % Nxyz(3)  ; do jj = 0, e % Nxyz(2)   ; do ii = 0, e % Nxyz(1)
!                  Q3D = Qx(ii,i) * Qy(jj,j) * Qz(kk,k)
!                  Uxf(:,ii,jj,kk) = Uxf(:,ii,jj,kk) + Q3D * e % storage % U_x(:,i,j,k)
!                  Uyf(:,ii,jj,kk) = Uyf(:,ii,jj,kk) + Q3D * e % storage % U_y(:,i,j,k)
!                  Uzf(:,ii,jj,kk) = Uzf(:,ii,jj,kk) + Q3D * e % storage % U_z(:,i,j,k)
!               end do                 ; end do                  ; end do
!            end do                  ; end do                   ; end do
!
!            end associate
!            
!            if (self % filterType == LPASS_FILTER) then
!               Uxf = e % storage % U_x - Uxf
!               Uyf = e % storage % U_y - Uyf
!               Uzf = e % storage % U_z - Uzf
!            end if
!         end if
!!
!!        ----------------
!!        Compute the flux
!!        ----------------
!!
!!
!!        ----------------------
!!        Get contravariant flux
!!        ----------------------
!!
!         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
!            call EllipticFlux( NCONS, NGRAD, e % storage % Q(:,i,j,k), Uxf(:,i,j,k), Uyf(:,i,j,k), &
!                               Uzf(:,i,j,k), mu(i,j,k), beta(i,j,k), kappa(i,j,k), cartesianFlux )
!            contravariantFluxF(:,i,j,k,IX) =    cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
!                                             +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
!                                             +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)
!
!
!            contravariantFluxF(:,i,j,k,IY) =    cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
!                                             +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
!                                             +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)
!
!
!            contravariantFluxF(:,i,j,k,IZ) =    cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
!                                             +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
!                                             +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)
!
!         end do               ; end do            ; end do
!!
!!        ----------------------
!!        Post-filtering ?
!!        ----------------------
!!
!         if (self % postFiltering) then
!            associate(Qx => self % filters(e % Nxyz(1)) % Q, & 
!                      Qy => self % filters(e % Nxyz(2)) % Q, &
!                      Qz => self % filters(e % Nxyz(3)) % Q    )
!
!            contravariantFlux = 0._RP
!            do k = 0, e % Nxyz(3)  ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
!               do kk = 0, e % Nxyz(3)  ; do jj = 0, e % Nxyz(2)   ; do ii = 0, e % Nxyz(1)
!                  Q3D = Qx(ii,i) * Qy(jj,j) * Qz(kk,k)
!                  contravariantFlux(:,ii,jj,kk,IX) = contravariantFlux(:,ii,jj,kk,IX) + Q3D * contravariantFluxF(:,i,j,k,IX)
!                  contravariantFlux(:,ii,jj,kk,IY) = contravariantFlux(:,ii,jj,kk,IY) + Q3D * contravariantFluxF(:,i,j,k,IY)
!                  contravariantFlux(:,ii,jj,kk,IZ) = contravariantFlux(:,ii,jj,kk,IZ) + Q3D * contravariantFluxF(:,i,j,k,IZ)
!               end do                 ; end do                  ; end do
!            end do                  ; end do                   ; end do
!
!            end associate
!            
!            if (self % filterType == LPASS_FILTER) then
!               contravariantFlux = contravariantFluxF - contravariantFlux
!            end if
!         else
!            contravariantFlux = contravariantFluxF
!         end if
!         
!      end subroutine SVV_ComputeInnerFluxes

!      subroutine SVV_RiemannSolver ( self, f, EllipticFlux, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, flux)
!         use SMConstants
!         use PhysicsStorage
!         use Physics
!         use FaceClass
!         use LESModels
!         implicit none
!         class(SVV_t)                :: self
!         class(Face),   intent(in)   :: f
!         procedure(EllipticFlux_f)   :: EllipticFlux
!         real(kind=RP), intent(in)   :: QLeft(NCONS, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: QRight (NCONS, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: U_xLeft(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: U_yLeft(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: U_zLeft(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: U_xRight(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: U_yRight(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(in)   :: U_zRight(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP), intent(out)  :: flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
!!
!!        ---------------
!!        Local variables
!!        ---------------
!!
!         integer           :: i, j, ii, jj
!         real(kind=RP)     :: Q(NCONS, 0:f % Nf(1), 0:f % Nf(2)) 
!         real(kind=RP)     :: U_x(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: U_y(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: U_z(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: Uxf(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: Uyf(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: Uzf(NGRAD, 0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: flux_vec(NCONS,NDIM)
!         real(kind=RP)     :: mu(0:f % Nf(1), 0:f % Nf(2)), kappa(0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: beta(0:f % Nf(1), 0:f % Nf(2))
!         real(kind=RP)     :: delta, Q2D
!         real(kind=RP)     :: fluxF(NCONS, 0:f % Nf(1), 0:f % Nf(2))
!!
!!        Interface averages
!!        ------------------
!         Q   = 0.5_RP * ( QLeft + QRight)
!         U_x = 0.5_RP * ( U_xLeft + U_xRight)
!         U_y = 0.5_RP * ( U_yLeft + U_yRight)
!         U_z = 0.5_RP * ( U_zLeft + U_zRight)
!!
!!        -------------------------
!!        Compute the SVV viscosity
!!        -------------------------
!!
!         if (self % muIsSmagorinsky) then
!!
!!           (1+Psvv) * muSmag
!!           -----------------
!            delta = sqrt(f % geom % surface / product(f % Nf + 1) )
!            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!               call Smagorinsky % ComputeViscosity ( delta, f % geom % dWall(i,j), Q  (:,i,j) &
!                                                                                 , U_x(:,i,j) &
!                                                                                 , U_y(:,i,j) &
!                                                                                 , U_z(:,i,j) &
!                                                                                 , mu(i,j) )
!               mu(i,j) = mu(i,j) !* (1._RP + self % Psvv)
!            end do                ; end do
!         else
!!
!!           Fixed value
!!           -----------
!            mu    = self % muSVV !/ maxval(f % Nf+1)
!         end if
!         
!         beta  = 0.0_RP
!         kappa = mu / ( thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Prt)
!!
!!        --------------------
!!        Filter the gradients
!!        --------------------
!!
!         if (.not. self % postFiltering) then
!            associate(Qx => self % filters(f % Nf(1)) % Q, & 
!                      Qy => self % filters(f % Nf(2)) % Q   )
!
!            Uxf = 0.0_RP   ; Uyf = 0.0_RP    ; Uzf = 0.0_RP
!            
!            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
!               do jj = 0, f % Nf(2)   ; do ii = 0, f % Nf(1)
!                  Q2D = Qx(ii,i) * Qy(jj,j) 
!                  Uxf(:,ii,jj) = Uxf(:,ii,jj) + Q2D * U_x(:,i,j)
!                  Uyf(:,ii,jj) = Uyf(:,ii,jj) + Q2D * U_y(:,i,j)
!                  Uzf(:,ii,jj) = Uzf(:,ii,jj) + Q2D * U_z(:,i,j)
!               end do                  ; end do
!            end do                   ; end do
!            end associate
!            
!            if (self % filterType == LPASS_FILTER) then
!               Uxf = U_x - Uxf
!               Uyf = U_y - Uyf
!               Uzf = U_z - Uzf
!            end if
!         end if
!!
!!        ----------------------------------
!!        Compute flux and project to normal
!!        ----------------------------------
!!
!         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
!            call EllipticFlux(NCONS, NGRAD, Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), mu(i,j), beta(i,j), kappa(i,j), flux_vec)
!            fluxF(:,i,j) =  flux_vec(:,IX) * f % geom % normal(IX,i,j) &
!                          + flux_vec(:,IY) * f % geom % normal(IY,i,j) &
!                          + flux_vec(:,IZ) * f % geom % normal(IZ,i,j) 
!         end do               ; end do
!!
!!        ---------------
!!        Filter the flux
!!        ---------------
!!
!         if (self % postFiltering) then
!            associate(Qx => self % filters(f % Nf(1)) % Q, & 
!                      Qy => self % filters(f % Nf(2)) % Q   )
!
!            flux = 0.0_RP
!            
!            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
!               do jj = 0, f % Nf(2)   ; do ii = 0, f % Nf(1)
!                  Q2D = Qx(ii,i) * Qy(jj,j) 
!                  flux(:,ii,jj) = flux(:,ii,jj) + Q2D * fluxF(:,i,j)
!               end do                  ; end do
!            end do                   ; end do
!            end associate
!            
!            if (self % filterType == LPASS_FILTER) then
!               flux = fluxF - flux
!            end if
!         else
!            flux = fluxF
!         end if
!
!      end subroutine SVV_RiemannSolver
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
         integer        :: sharpCutOff
         real(kind=RP)  :: Nodal2Modal(0:N,0:N)
         real(kind=RP)  :: Modal2Nodal(0:N,0:N)
         real(kind=RP)  :: filterCoefficients(0:N)
         real(kind=RP)  :: Lkj(0:N,0:N), dLk_dummy
         real(kind=RP)  :: normLk(0:N)

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
         select case (self % filterShape)
            case (POW_FILTER)
               do k = 0, N
                  filterCoefficients(k) = (real(k, kind=RP) / N + 1.0e-12_RP) ** self % Psvv
               end do
            
            case (SHARP_FILTER)
               sharpCutOff = nint(self % Psvv)
               if (sharpCutOff >= N) then
                  write(STD_OUT) 'ERROR :: sharp cut-off must be lower than N'
                  stop
               end if
               filterCoefficients = 0._RP
               filterCoefficients(sharpCutOff+1:N) = 1._RP
               
            case (EXP_FILTER)
               filterCoefficients = 0._RP
               do k = 0, N
                  if (k > self % Psvv) filterCoefficients(k) = exp( -real( (k-N)**2 , kind=RP) / (k - self % Psvv) ** 2 )
               end do
               
         end select
         
         if (self % filterType == LPASS_FILTER) then
            filterCoefficients = 1._RP - filterCoefficients
         end if
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
