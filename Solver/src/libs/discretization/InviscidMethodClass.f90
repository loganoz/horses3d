!
!//////////////////////////////////////////////////////
!
!   @File:    InviscidMethodClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:30 2017
!   @Last revision date: Sat Dec 16 11:09:33 2017
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: ddec1961effb920e19c9b28e302f01273406797c
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module InviscidMethodClass
   use SMConstants
   use RiemannSolvers
   implicit none

   private
   public   InviscidMethod_t, volumetricSharpFlux_FCN

   integer,    parameter   :: STANDARD_DG = 1
   integer,    parameter   :: SPLIT_DG = 2

   type InviscidMethod_t
      procedure(VolumetricSharpFlux_FCN), nopass, pointer  :: ComputeVolumetricSharpFlux => NULL()
      contains
         procedure   :: Initialize             => BaseClass_Initialize
         procedure   :: ComputeInnerFluxes     => BaseClass_ComputeInnerFluxes
         procedure   :: ComputeSplitFormFluxes => BaseClass_ComputeSplitFormFluxes
   end type InviscidMethod_t

   abstract interface
      subroutine VolumetricSharpFlux_FCN(QL,QR,JaL,JaR,fSharp) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(1:NCONS)
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
!
!        Set up the Inviscid Discretization
!        ----------------------------------
         splitType = STANDARD_SPLIT

         call SetRiemannSolver( whichRiemannSolver, splitType )
!
!        ********
!        Describe
!        ********
!
         call Subsection_Header("Inviscid discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Standard"

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
         
         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab
         
         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Disabled."
         end if

      end subroutine BaseClass_Initialize

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

         call InviscidFlux( e%Nxyz, e % storage % Q, cartesianFlux)

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
end module InviscidMethodClass
