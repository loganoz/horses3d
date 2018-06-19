!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicDiscretizationClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:30 2017
!   @Last revision date: Wed Jun 20 18:14:32 2018
!   @Last revision author: Juan Manzanero (j.manzanero1992@gmail.com)
!   @Last revision commit: 9c8ed8b6306ad0912cb55b510aa73d1610bb1cb5
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
#define NNS NCONS
#define NGRADNS NGRAD
#elif defined(INCNS)
#define NNS NINC
#define NGRADNS NINC
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
module HyperbolicDiscretizationClass
   use SMConstants
#if defined(NAVIERSTOKES)
   use RiemannSolvers_NS
#elif defined(INCNS)
   use RiemannSolvers_iNS
#endif
   implicit none

   private
   public   HyperbolicDiscretization_t, volumetricSharpFlux_FCN

   integer,    parameter   :: STANDARD_DG = 1
   integer,    parameter   :: SPLIT_DG = 2

   type HyperbolicDiscretization_t
      procedure(VolumetricSharpFlux_FCN), nopass, pointer  :: ComputeVolumetricSharpFlux => NULL()
      contains
         procedure   :: Initialize               => BaseClass_Initialize
         procedure   :: ComputeInnerFluxes       => BaseClass_ComputeInnerFluxes
         procedure   :: ComputeSplitFormFluxes   => BaseClass_ComputeSplitFormFluxes
#if defined(NAVIERSTOKES)
         procedure   :: ComputeInnerFluxJacobian => BaseClass_ComputeInnerFluxJacobian
#endif
   end type HyperbolicDiscretization_t

   abstract interface
      pure subroutine HyperbolicFlux0D_f(Q, F)
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)  :: Q   (1:NNS     )
         real(kind=RP), intent(out) :: F(1:NNS, 1:NDIM)
      end subroutine HyperbolicFlux0D_f

      pure subroutine HyperbolicFlux2D_f( N, Q, F)
         use SMConstants
         use PhysicsStorage
         implicit none
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:NNS, 0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:NNS, 1:NDIM, 0:N(1), 0:N(2))
      end subroutine HyperbolicFlux2D_f

      pure subroutine HyperbolicFlux3D_f( N, Q, F)
         use SMConstants
         use PhysicsStorage
         implicit none
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:NNS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:NNS, 0:N(1), 0:N(2), 0:N(3), 1:NDIM )
      end subroutine HyperbolicFlux3D_f

      subroutine VolumetricSharpFlux_FCN(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NNS)
         real(kind=RP), intent(in)       :: QR(1:NNS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(1:NNS)
      end subroutine VolumetricSharpFlux_FCN
   end interface
!
!  ========
   contains
!  ========
!
      subroutine BaseClass_Initialize(self, controlVariables)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(HyperbolicDiscretization_t) :: self
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
!
!        Set up the Hyperbolic Discretization
!        ----------------------------------
         splitType = STANDARD_SPLIT

         call SetRiemannSolver( whichRiemannSolver, splitType )
!
!        ********
!        Describe
!        ********
!
         call Subsection_Header("Hyperbolic discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Standard"

         select case (whichRiemannSolver)
#if defined(NAVIERSTOKES)
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
         
         case (RIEMANN_MATRIXDISS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Matrix dissipation"
         
         case (RIEMANN_LOWDISSROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Low dissipation Roe"

#elif defined(INCNS)
         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

         case (RIEMANN_LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"

         case (RIEMANN_EXACT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Exact"

#endif         
         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab
         
         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Disabled."
         end if

      end subroutine BaseClass_Initialize

      subroutine BaseClass_ComputeInnerFluxes( self , e , HyperbolicFlux, contravariantFlux )
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         class(HyperbolicDiscretization_t), intent(in)  :: self
         type(Element),           intent(in)  :: e
         procedure(HyperbolicFlux3D_f)        :: HyperbolicFlux
         real(kind=RP),           intent(out) :: contravariantFlux(1:NNS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: i, j, k
         real(kind=RP)      :: cartesianFlux(1:NNS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)

         call HyperbolicFlux( e%Nxyz, e % storage % Q, cartesianFlux)

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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!     Subroutine to get the Jacobian of the contravariant fluxes
!     -> dFdQ (:,:,i,j,k,dim)
!                 |     |
!              jac|coord|flux in cartesian direction dim 
!     ----------------------------------------------------------
#if defined(NAVIERSTOKES)
      subroutine BaseClass_ComputeInnerFluxJacobian( self, e, dFdQ) 
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(HyperbolicDiscretization_t), intent(in)  :: self
         type(Element)          , intent(in)  :: e
         real(kind=RP)          , intent(out) :: dFdQ( NNS, NNS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3), NDIM )
         !--------------------------------------------
         real(kind=RP), DIMENSION(NNS,NNS)  :: dfdq_,dgdq_,dhdq_
         integer                                :: i,j,k
         !--------------------------------------------
         
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            
            call InviscidJacobian(e % storage % Q(:,i,j,k),dfdq_,dgdq_,dhdq_)
            
            
            dFdQ(:,:,i,j,k,IX) = e % geom % jGradXi  (1,i,j,k) * dfdq_ + &
                                 e % geom % jGradXi  (2,i,j,k) * dgdq_ + &
                                 e % geom % jGradXi  (3,i,j,k) * dhdq_ 

            dFdQ(:,:,i,j,k,IY) = e % geom % jGradEta (1,i,j,k) * dfdq_ + &
                                 e % geom % jGradEta (2,i,j,k) * dgdq_ + &
                                 e % geom % jGradEta (3,i,j,k) * dhdq_ 

            dFdQ(:,:,i,j,k,IZ) = e % geom % jGradZeta(1,i,j,k) * dfdq_ + &
                                 e % geom % jGradZeta(2,i,j,k) * dgdq_ + &
                                 e % geom % jGradZeta(3,i,j,k) * dhdq_ 
         end do                ; end do                ; end do
         
      end subroutine BaseClass_ComputeInnerFluxJacobian
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeSplitFormFluxes(self, e, contravariantFlux, fSharp, gSharp, hSharp)
         use ElementClass
         use PhysicsStorage
         implicit none
         class(HyperbolicDiscretization_t), intent(in)  :: self
         type(Element),           intent(in)  :: e
         real(kind=RP),           intent(in)  :: contravariantFlux(1:NNS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP),           intent(out) :: fSharp(1:NNS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),           intent(out) :: gSharp(1:NNS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),           intent(out) :: hSharp(1:NNS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )

      end subroutine BaseClass_ComputeSplitFormFluxes
end module HyperbolicDiscretizationClass
#endif
