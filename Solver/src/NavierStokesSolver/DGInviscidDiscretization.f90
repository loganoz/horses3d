#include "Includes.h"
module DGInviscidDiscretization
   use SMConstants
   implicit none

   private
   public   InviscidMethod_t , StandardDG_t , SplitDG_t, InviscidMethod

   integer,    parameter   :: STANDARD_DG = 1
   integer,    parameter   :: SPLIT_DG = 2

   integer,    parameter   :: STANDARD_SPLIT = 1
   integer,    parameter   :: MORINISHI_SPLIT = 2
   integer,    parameter   :: DUCROS_SPLIT = 3
   integer,    parameter   :: KENNEDYGRUBER_SPLIT = 4
   integer,    parameter   :: PIROZZOLI_SPLIT = 5

   type InviscidMethod_t
      procedure(VolumetricSharpFlux_FCN), nopass, pointer  :: ComputeVolumetricSharpFlux => NULL()
      contains
         procedure   :: Initialize           => InitializeInviscidMethod
         procedure   :: ComputeInnerFluxes   => BaseClass_ComputeInnerFluxes
         procedure   :: ComputeSplitFormFluxes => BaseClass_ComputeSplitFormFluxes
   end type InviscidMethod_t

   type, extends(InviscidMethod_t)  :: StandardDG_t
   end type StandardDG_t

   type, extends(InviscidMethod_t)  :: SplitDG_t
      contains 
         procedure   :: ComputeSplitFormFluxes => SplitDG_ComputeSplitFormFluxes
   end type SplitDG_t

   class(InviscidMethod_t), allocatable         :: InviscidMethod

   abstract interface
      function VolumetricSharpFlux_FCN(QL,QR,JaL,JaR) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), dimension(NCONS) :: VolumetricSharpFlux_FCN
      end function VolumetricSharpFlux_FCN
   end interface
!
!  ========
   contains
!  ========
!
      subroutine InitializeInviscidMethod(self, controlVariables)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(InviscidMethod_t) :: self
         class(FTValueDictionary),  intent(in)   :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)    :: splitForm
         integer                       :: splitType
         interface
            subroutine toLower(str)
               character(*), intent(in out) :: str
            end subroutine toLower
         end interface


         select type ( self ) 
         type is (SplitDG_t)
            splitForm = controlVariables % stringValueForKey(splitFormKey, requestedLength = LINE_LENGTH)
   
            call toLower(splitForm)
   
            select case ( trim(splitForm) )
            case ( "standard" )
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
   
            case default
if ( MPI_Process % isRoot ) then   
               write(STD_OUT,'(A,A,A)') 'Requested split form "',trim(splitForm),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented split forms are:"
               write(STD_OUT,'(A)') "  * Standard"
               write(STD_OUT,'(A)') "  * Morinishi"
               write(STD_OUT,'(A)') "  * Ducros"
               write(STD_OUT,'(A)') "  * Kennedy-Gruber"
               write(STD_OUT,'(A)') "  * Pirozzoli"
               errorMessage(STD_OUT)
               stop 
end if
   
            end select
         end select
!
!        Describe
!        --------
         call Subsection_Header("Inviscid discretization")

if (MPI_Process % isRoot ) then
         select type ( self ) 
   
         type is (StandardDG_t)
            write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Standard"

         type is (SplitDG_t)
            write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Split-Form"

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

            end select
         end select

         select case (riemannSolverChoice)
         case (ROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Roe"

         case (LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"

         case (RUSANOV)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Rusanov"

         case (DUCROS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Ducros"

         case (MORINISHI)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Morinishi"

         case (KENNEDYGRUBER)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Kennedy-Gruber"

         case (PIROZZOLI)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Pirozzoli"

         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab
         
         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Disabled."
         end if
end if

      end subroutine InitializeInviscidMethod
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           BaseClass Procedures
!           --------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeInnerFluxes( self , e , contravariantFlux )
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         class(InviscidMethod_t), intent(in)  :: self
         type(Element),           intent(in)  :: e
         real(kind=RP),           intent(out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: i, j, k
         real(kind=RP)      :: cartesianFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)

         cartesianFlux = InviscidFlux( e%Nxyz(1) , e%Nxyz(2) , e%Nxyz(3) , e % storage % Q )

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
         
            contravariantFlux(:,i,j,k,IX) =    cartesianFlux(:,i,j,k,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             + cartesianFlux(:,i,j,k,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             + cartesianFlux(:,i,j,k,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =   cartesianFlux(:,i,j,k,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             + cartesianFlux(:,i,j,k,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             + cartesianFlux(:,i,j,k,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =   cartesianFlux(:,i,j,k,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             + cartesianFlux(:,i,j,k,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             + cartesianFlux(:,i,j,k,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do                ; end do

      end subroutine BaseClass_ComputeInnerFluxes

      subroutine BaseClass_ComputeSplitFormFluxes(self, e, contravariantFlux, fSharp, gSharp, hSharp)
         use ElementClass
         use PhysicsStorage
         implicit none
         class(InviscidMethod_t), intent(in)  :: self
         type(Element),           intent(in)  :: e
         real(kind=RP),           intent(in)  :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP),           intent(out) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),           intent(out) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),           intent(out) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )

      end subroutine BaseClass_ComputeSplitFormFluxes
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           StandardDG Procedures
!           ---------------------
!///////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           SplitDG Procedures
!           ------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine SplitDG_ComputeSplitFormFluxes(self, e, contravariantFlux, fSharp, gSharp, hSharp)
         use ElementClass
         use PhysicsStorage
         implicit none
         class(SplitDG_t), intent(in)  :: self
         type(Element),           intent(in)  :: e
         real(kind=RP),           intent(in)  :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP),           intent(out) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),           intent(out) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),           intent(out) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
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
               fSharp(:,l,i,j,k) = self % ComputeVolumetricSharpFlux(Q(:,i,j,k), Q(:,l,j,k), e % geom % jGradXi(:,i,j,k), e % geom % jGradXi(:,l,j,k))
               fSharp(:,i,l,j,k) = fSharp(:,l,i,j,k)
            end do            
         end do               ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = j+1, e%Nxyz(2)
               gSharp(:,l,i,j,k) = self % ComputeVolumetricSharpFlux(Q(:,i,j,k), Q(:,i,l,k), e % geom % jGradEta(:,i,j,k), e % geom % jGradEta(:,i,l,k))
               gSharp(:,j,i,l,k) = gSharp(:,l,i,j,k)
            end do            
         end do               ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            do l = k+1, e%Nxyz(3)
               hSharp(:,l,i,j,k) = self % ComputeVolumetricSharpFlux(Q(:,i,j,k), Q(:,i,j,l), e % geom % jGradZeta(:,i,j,k), e % geom % jGradZeta(:,i,j,l))
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
      function StandardDG_VolumetricSharpFlux(QL,QR,JaL,JaR) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), dimension(NCONS) :: StandardDG_VolumetricSharpFlux
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
         associate ( fSharp => StandardDG_VolumetricSharpFlux )

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
         end associate

      end function StandardDG_VolumetricSharpFlux

      function Morinishi_VolumetricSharpFlux(QL,QR,JaL,JaR) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), dimension(NCONS) :: Morinishi_VolumetricSharpFlux
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
   
         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))
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
         associate ( fSharp => Morinishi_VolumetricSharpFlux )

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
         end associate

      end function Morinishi_VolumetricSharpFlux

      function Ducros_VolumetricSharpFlux(QL,QR,JaL,JaR) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), dimension(NCONS) :: Ducros_VolumetricSharpFlux
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
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         associate ( fSharp => Ducros_VolumetricSharpFlux )

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

         end associate

      end function Ducros_VolumetricSharpFlux
 
      function KennedyGruber_VolumetricSharpFlux(QL,QR,JaL,JaR) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), dimension(NCONS) :: KennedyGruber_VolumetricSharpFlux
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
         e   = 0.5_RP * (QL(IRHOE)*invRhoL + QR(IRHOE)*invRhoR)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         associate ( fSharp => KennedyGruber_VolumetricSharpFlux )

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

         end associate

      end function KennedyGruber_VolumetricSharpFlux

      function Pirozzoli_VolumetricSharpFlux(QL,QR,JaL,JaR) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), dimension(NCONS) :: Pirozzoli_VolumetricSharpFlux
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
         associate ( fSharp => Pirozzoli_VolumetricSharpFlux )

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

      end function Pirozzoli_VolumetricSharpFlux

end module DGInviscidDiscretization
