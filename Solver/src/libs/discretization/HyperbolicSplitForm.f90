#include "Includes.h"
#if defined(NAVIERSTOKES)
#define NNS NCONS
#define NGRADNS NGRAD
#elif defined(INCNS)
#define NNS NINC
#define NGRADNS NINC
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
module HyperbolicSplitForm
   use SMConstants
#if defined(NAVIERSTOKES)
   use RiemannSolvers_NS
#elif defined(INCNS)
   use RiemannSolvers_iNS
#endif
   use HyperbolicDiscretizationClass
   use FluidData
   implicit none

   private
   public   SplitDG_t

   type, extends(HyperbolicDiscretization_t)  :: SplitDG_t
      contains 
         procedure   :: Initialize             => SplitDG_Initialize
         procedure   :: ComputeSplitFormFluxes => SplitDG_ComputeSplitFormFluxes
   end type SplitDG_t
!
!  ========
   contains
!  ========
!
      subroutine SplitDG_Initialize(self, controlVariables)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(SplitDG_t) :: self
         class(FTValueDictionary),  intent(in)   :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)    :: splitForm
         integer                       :: splitType
!
!        Get the split form from control variables
!        -----------------------------------------
         splitForm = controlVariables % stringValueForKey(splitFormKey, requestedLength = LINE_LENGTH)
         call toLower(splitForm)

         select case ( trim(splitForm) )
#if defined(NAVIERSTOKES)
         case ( "standard" )
!
!           Skew-symmetric version of the standard DG. Useful for testing
!           -------------------------------------------------------------
            self % ComputeVolumetricSharpFlux => StandardDG_VolumetricSharpFlux
            splitType = STANDARD_SPLIT

         case ( "morinishi" )
            self % ComputeVolumetricSharpFlux => Morinishi_VolumetricSharpFlux
            splitType = MORINISHI_SPLIT

         case ( "ducros" )
            self % ComputeVolumetricSharpFlux => Ducros_VolumetricSharpFlux
            splitType = DUCROS_SPLIT

         case ( "kennedy-gruber" )
            self % ComputeVolumetricSharpFlux => KennedyGruber_VolumetricSharpFlux
            splitType = KENNEDYGRUBER_SPLIT

         case ( "pirozzoli" )
            self % ComputeVolumetricSharpFlux => Pirozzoli_VolumetricSharpFlux
            splitType = PIROZZOLI_SPLIT

         case ( "entropy conserving")
            self % ComputeVolumetricSharpFlux => EntropyConserving_VolumetricSharpFlux
            splitType = ENTROPYCONS_SPLIT

         case ( "entropy and energy conserving")
            self % ComputeVolumetricSharpFlux => EntropyAndEnergyConserving_VolumetricSharpFlux
            splitType = ENTROPYANDENERGYCONS_SPLIT

         case default
            if ( MPI_Process % isRoot ) then   
            write(STD_OUT,'(A,A,A)') 'Requested split form "',trim(splitForm),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented split forms are:"
            write(STD_OUT,'(A)') "  * Standard"
            write(STD_OUT,'(A)') "  * Morinishi"
            write(STD_OUT,'(A)') "  * Ducros"
            write(STD_OUT,'(A)') "  * Kennedy-Gruber"
            write(STD_OUT,'(A)') "  * Pirozzoli"
            write(STD_OUT,'(A)') "  * Entropy conserving"
            write(STD_OUT,'(A)') "  * Entropy and energy conserving"
            errorMessage(STD_OUT)
            stop 
            end if
#endif
         end select

         call SetRiemannSolver( whichRiemannSolver, splitType )
!
!        ********
!        Describe
!        ********
!
         call Subsection_Header("Hyperbolic discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Split-Form"
#if defined(NAVIERSTOKES)
         select case ( splitType )
         case (STANDARD_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Standard"

         case (MORINISHI_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Morinishi"

         case (DUCROS_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Ducros"

         case (KENNEDYGRUBER_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Kennedy-Gruber"

         case (PIROZZOLI_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Pirozzoli"

         case (ENTROPYCONS_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Entropy conserving"

         case (ENTROPYANDENERGYCONS_SPLIT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Split form scheme: ","Entropy and energy conserving"

         end select

         select case (whichRiemannSolver)
         case (RIEMANN_ROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Roe"

         case (RIEMANN_LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"

         case (RIEMANN_RUSANOV)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Rusanov"

         case (RIEMANN_STDROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Standard Roe"

         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

         case (RIEMANN_ROEPIKE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Roe-Pike"
         
         case (RIEMANN_LOWDISSROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Low dissipation Roe"
         
         case (RIEMANN_MATRIXDISS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Matrix dissipation"
         
         case (RIEMANN_VISCOUSNS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Viscous NS"
         
         end select
#endif

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab
         
         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Disabled."
         end if

      end subroutine SplitDG_Initialize

      subroutine SplitDG_ComputeSplitFormFluxes(self, e, contravariantFlux, fSharp, gSharp, hSharp)
         use ElementClass
         use PhysicsStorage
         implicit none
         class(SplitDG_t), intent(in)  :: self
         type(Element),    intent(in)  :: e
         real(kind=RP),    intent(in)  :: contravariantFlux(1:NNS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP),    intent(out) :: fSharp(1:NNS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),    intent(out) :: gSharp(1:NNS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),    intent(out) :: hSharp(1:NNS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
!
!        ---------------
!        Local variables
!        ---------------
!
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
               call self % ComputeVolumetricSharpFlux(Q(:,i,j,k), Q(:,l,j,k), e % geom % jGradXi(:,i,j,k), e % geom % jGradXi(:,l,j,k), fSharp(:,l,i,j,k)) 
               fSharp(:,i,l,j,k) = fSharp(:,l,i,j,k)
            end do            
         end do               ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = j+1, e%Nxyz(2)
               call self % ComputeVolumetricSharpFlux(Q(:,i,j,k), Q(:,i,l,k), e % geom % jGradEta(:,i,j,k), e % geom % jGradEta(:,i,l,k), gSharp(:,l,i,j,k)) 
               gSharp(:,j,i,l,k) = gSharp(:,l,i,j,k)
            end do            
         end do               ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = k+1, e%Nxyz(3)
               call self % ComputeVolumetricSharpFlux(Q(:,i,j,k), Q(:,i,j,l), e % geom % jGradZeta(:,i,j,k), e % geom % jGradZeta(:,i,j,l), hSharp(:,l,i,j,k))
               hSharp(:,k,i,j,l) = hSharp(:,l,i,j,k)
            end do            
         end do               ; end do             ; end do

         end associate

      end subroutine SplitDG_ComputeSplitFormFluxes
!
!///////////////////////////////////////////////////////////////////////
!
!        Volumetric sharp flux library
!        -----------------------------
!///////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine StandardDG_VolumetricSharpFlux(QL,QR,JaL,JaR, fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))
!
!        Average metrics: (Note: Here all average (1/2)s are accounted later)
!        ---------------
         Ja = (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = ( QL(IRHOU) + QR(IRHOU) )
         f(IRHOU) = ( QL(IRHOU) * uL + QR(IRHOU) * uR + pL + pR )
         f(IRHOV) = ( QL(IRHOU) * vL + QR(IRHOU) * vR )
         f(IRHOW) = ( QL(IRHOU) * wL + QR(IRHOU) * wR )
         f(IRHOE) = ( uL*(QL(IRHOE) + pL) + uR*(QR(IRHOE) + pR) )

         g(IRHO)  = (QL(IRHOV) + QR(IRHOV))
         g(IRHOU) = ( QL(IRHOV) * uL + QR(IRHOV) * uR )
         g(IRHOV) = ( QL(IRHOV) * vL + QR(IRHOV) * vR + pL + pR )
         g(IRHOW) = ( QL(IRHOV) * wL + QR(IRHOV) * wR )
         g(IRHOE) = ( vL*(QL(IRHOE) + pL) + vR*(QR(IRHOE) + pR) )

         h(IRHO)  = (QL(IRHOW) + QR(IRHOW))
         h(IRHOU) = (QL(IRHOW)*uL + QR(IRHOW)*uR)
         h(IRHOV) = (QL(IRHOW)*vL + QR(IRHOW)*vR)
         h(IRHOW) = (QL(IRHOW)*wL + QR(IRHOW)*wR + pL + pR)
         h(IRHOE) = ( wL*(QL(IRHOE) + pL) + wR*(QR(IRHOE) + pR) )
!
!        Compute the sharp flux: (And account for the (1/2)^2)
!        ----------------------         
         fSharp = 0.25_RP * ( f*Ja(IX) + g*Ja(IY) + h*Ja(IZ) )

      end subroutine StandardDG_VolumetricSharpFlux

      subroutine Morinishi_VolumetricSharpFlux(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL, hL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR, hR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)
   
         associate(gm1 => thermodynamics % GammaMinus1)
   
         pL = gm1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                             + QL(IRHOV) * vL &
                                             + QL(IRHOW) * wL ))

         pR = gm1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                             + QR(IRHOV) * vR &
                                             + QR(IRHOW) * wR ))
         end associate
!
!        Here the enthalpy does not contain the kinetic energy
!        -----------------------------------------------------
         hL = dimensionless % cp * pL  ; hR = dimensionless % cp * pR 
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = 0.5_RP * ( QL(IRHOU) + QR(IRHOU) )
         f(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * ( uL + uR ) + 0.5_RP * ( pL + pR )
         f(IRHOV) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * ( vL + vR )
         f(IRHOW) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * ( wL + wR )
         f(IRHOE) = 0.5_RP * ( uL*hL + uR*hR) + 0.25_RP * ( QL(IRHOU)*uL + QR(IRHOU)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QL(IRHOU)*vL + QR(IRHOU)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QL(IRHOU)*wL + QR(IRHOU)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QL(IRHOU)*POW2(uL) + QR(IRHOU)*POW2(uR)   ) &
                                              - 0.25_RP * ( QL(IRHOU)*POW2(vL) + QR(IRHOU)*POW2(vR)   ) &
                                              - 0.25_RP * ( QL(IRHOU)*POW2(wL) + QR(IRHOU)*POW2(wR)   ) 

         g(IRHO)  = 0.5_RP * ( QL(IRHOV) + QR(IRHOV) )
         g(IRHOU) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * ( uL + uR )
         g(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * ( vL + vR ) + 0.5_RP * ( pL + pR )
         g(IRHOW) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * ( wL + wR )
         g(IRHOE) = 0.5_RP * ( vL*hL + vR*hR) + 0.25_RP * ( QL(IRHOV)*uL + QR(IRHOV)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QL(IRHOV)*vL + QR(IRHOV)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QL(IRHOV)*wL + QR(IRHOV)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QL(IRHOV)*POW2(uL) + QR(IRHOV)*POW2(uR)   ) &
                                              - 0.25_RP * ( QL(IRHOV)*POW2(vL) + QR(IRHOV)*POW2(vR)   ) &
                                              - 0.25_RP * ( QL(IRHOV)*POW2(wL) + QR(IRHOV)*POW2(wR)   ) 

         h(IRHO)  = 0.5_RP * ( QL(IRHOW) + QR(IRHOW) )
         h(IRHOU) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * ( uL + uR )
         h(IRHOV) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * ( vL + vR )
         h(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * ( wL + wR ) + 0.5_RP * ( pL + pR )
         h(IRHOE) = 0.5_RP * ( wL*hL + wR*hR) + 0.25_RP * ( QL(IRHOW)*uL + QR(IRHOW)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QL(IRHOW)*vL + QR(IRHOW)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QL(IRHOW)*wL + QR(IRHOW)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QL(IRHOW)*POW2(uL) + QR(IRHOW)*POW2(uR)   ) &
                                              - 0.25_RP * ( QL(IRHOW)*POW2(vL) + QR(IRHOW)*POW2(vR)   ) &
                                              - 0.25_RP * ( QL(IRHOW)*POW2(wL) + QR(IRHOW)*POW2(wR)   ) 
!
!        Compute the sharp flux
!        ----------------------         
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ) 

      end subroutine Morinishi_VolumetricSharpFlux

      subroutine Ducros_VolumetricSharpFlux(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)
         
         associate(gm1 => thermodynamics % GammaMinus1)
   
         pL = gm1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                             + QL(IRHOV) * vL &
                                             + QL(IRHOW) * wL ))

         pR = gm1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                             + QR(IRHOV) * vR &
                                             + QR(IRHOW) * wR ))
         end associate
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = 0.25_RP * ( QL(IRHO) + QR(IRHO) ) * (uL + uR)
         f(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * (uL + uR) + 0.5_RP * (pL + pR)
         f(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * (uL + uR)
         f(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * (uL + uR)
         f(IRHOE) = 0.25_RP * ( QL(IRHOE) + pL + QR(IRHOE) + pR ) * (uL + uR)

         g(IRHO)  = 0.25_RP * ( QL(IRHO) + QR(IRHO) ) * (vL + vR)
         g(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * (vL + vR)
         g(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * (vL + vR) + 0.5_RP * (pL + pR)
         g(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * (vL + vR)
         g(IRHOE) = 0.25_RP * ( QL(IRHOE) + pL + QR(IRHOE) + pR ) * (vL + vR)

         h(IRHO)  = 0.25_RP * ( QL(IRHO) + QR(IRHO) ) * (wL + wR)
         h(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * (wL + wR)
         h(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * (wL + wR)
         h(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * (wL + wR) + 0.5_RP * (pL + pR)
         h(IRHOE) = 0.25_RP * ( QL(IRHOE) + pL + QR(IRHOE) + pR ) * (wL + wR)
!
!        Compute the sharp flux
!        ----------------------         
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine Ducros_VolumetricSharpFlux
 
      subroutine KennedyGruber_VolumetricSharpFlux(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: rho, u, v, w, e, p
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         associate(gm1 => thermodynamics % GammaMinus1)
   
         pL = gm1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                             + QL(IRHOV) * vL &
                                             + QL(IRHOW) * wL ))

         pR = gm1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                             + QR(IRHOV) * vR &
                                             + QR(IRHOW) * wR ))
         end associate

         rho = 0.5_RP * (QL(IRHO) + QR(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         e   = 0.5_RP * (QL(IRHOE)*invRhoL + QR(IRHOE)*invRhoR)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = rho * u
         f(IRHOU) = rho * u * u + p
         f(IRHOV) = rho * u * v
         f(IRHOW) = rho * u * w
         f(IRHOE) = rho * u * e + p * u
         
         g(IRHO)  = rho * v
         g(IRHOU) = rho * v * u
         g(IRHOV) = rho * v * v + p
         g(IRHOW) = rho * v * w
         g(IRHOE) = rho * v * e + p * v

         h(IRHO)  = rho * w
         h(IRHOU) = rho * w * u
         h(IRHOV) = rho * w * v
         h(IRHOW) = rho * w * w + p
         h(IRHOE) = rho * w * e + p * w
!
!        Compute the sharp flux
!        ----------------------         
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine KennedyGruber_VolumetricSharpFlux

      subroutine Pirozzoli_VolumetricSharpFlux(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: rho, u, v, w, h, p
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))

         rho = 0.5_RP * (QL(IRHO) + QR(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         h   = 0.5_RP * ((QL(IRHOE)+pL)*invRhoL + (QR(IRHOE)+pR)*invRhoR)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h
         
         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h 
!
!        Compute the sharp flux
!        ----------------------         
         fSharp = ff*Ja(IX) + gg*Ja(IY) + hh*Ja(IZ)

      end subroutine Pirozzoli_VolumetricSharpFlux

      subroutine EntropyConserving_VolumetricSharpFlux(QL,QR,JaL,JaR,fSharp) 
!
!        *******************************************************************
!           Entropy conserving split form by Ismail and Roe. 
!           Parameter vector definition:
!
!           z1 = sqrt(rho / p)
!           z2 = z1 * u
!           z3 = z1 * v
!           z4 = z1 * w
!           z5 = sqrt(p rho)
!
!           Averaged states:
!           
!           rho = {{z1}} {{z5:ln}}
!           u   = {{z2}} / {{z1}}
!           v   = {{z3}} / {{z1}}
!           w   = {{z4}} / {{z1}}
!           p   = {{z5}} / {{z1}}
!           h   = g*p2/(rho(g-1)) + 0.5*(u^2 + v^2 + w^2)
!
!           where
!           p2   = (g+1)/(2g) {{z5:ln}}/{{z1:ln}} + (g-1)/(2g){{z5}}/{{z1}}
!           
!        *******************************************************************
!               
         use SMConstants
         use PhysicsStorage
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, rhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, rhoR, uR, vR, wR, pR
         real(kind=RP)     :: rho, u, v, w, h, p, p2
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: zL(NCONS), zR(NCONS), zAv(NCONS), invZ1Sum
         real(kind=RP)     :: z5Log, z1Log, invZ1Av
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         associate ( gammaPlus1Div2      => thermodynamics % gammaPlus1Div2, &
                     gammaMinus1Div2     => thermodynamics % gammaMinus1Div2, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1, &
                     invGamma            => thermodynamics % invGamma ) 

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         rhoL = QL(IRHO)               ; rhoR = QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))

!
!        Compute Ismail and Roe parameter vector
!        ---------------------------------------
         zL(5) = sqrt(rhoL*pL)      ; zR(5) = sqrt(rhoR*pR)
         zL(1) = rhoL / zL(5)       ; zR(1) = rhoR / zR(5)
         zL(2) = zL(1) * uL         ; zR(2) = zR(1) * uR
         zL(3) = zL(1) * vL         ; zR(3) = zR(1) * vR
         zL(4) = zL(1) * wL         ; zR(4) = zR(1) * wR

         zAv = 0.5_RP * (zL + zR)
         invZ1Av = 1.0_RP / zAv(1)

         call logarithmicMean(zL(1),zR(1), z1Log)
         call logarithmicMean(zL(5),zR(5), z5Log)

         rho = zAv(1) * z5Log
         u   = zAv(2) * invZ1Av
         v   = zAv(3) * invZ1Av
         w   = zAv(4) * invZ1Av
         p   = zAv(5) * invZ1Av
         p2  = (gammaPlus1Div2 * z5Log / z1Log + gammaMinus1Div2 * p) * invGamma
         h   = gammaDivGammaMinus1 * p2 / rho + 0.5_RP*(POW2(u) + POW2(v) + POW2(w))
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h
         
         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h 
!
!        Compute the sharp flux
!        ----------------------         
         fSharp = ff*Ja(IX) + gg*Ja(IY) + hh*Ja(IZ)

         end associate

      end subroutine EntropyConserving_VolumetricSharpFlux

      subroutine EntropyAndEnergyConserving_VolumetricSharpFlux(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, rhoL, uL, vL, wL, pL, betaL
         real(kind=RP)     :: invRhoR, rhoR, uR, vR, wR, pR, betaR
         real(kind=RP)     :: rho, u, v, w, h, p, betaLog
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         rhoL = QL(IRHO)               ; rhoR = QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))

!
!        Compute Chandrasekar's variables
!        --------------------------------
         betaL = 0.5_RP * rhoL / pL    ; betaR = 0.5_RP * rhoR / pR
         call logarithmicMean(betaL, betaR, betaLog)

         call logarithmicMean(rhoL,rhoR,rho)
         u   = AVERAGE(uL, uR)
         v   = AVERAGE(vL, vR)
         w   = AVERAGE(wL, wR)
         p   = 0.5_RP * (rhoL + rhoR) / (betaL + betaR)
         h   =   0.5_RP/(betaLog*(gammaMinus1)) &
               - 0.5_RP*AVERAGE(POW2(uL)+POW2(vL)+POW2(wL), POW2(uR)+POW2(vR)+POW2(wR)) &
               + p/rho + POW2(u) + POW2(v) + POW2(w)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h
         
         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h 
!
!        Compute the sharp flux
!        ----------------------         
         fSharp = ff*Ja(IX) + gg*Ja(IY) + hh*Ja(IZ)

         end associate

      end subroutine EntropyAndEnergyConserving_VolumetricSharpFlux
#endif
end module HyperbolicSplitForm
#endif
