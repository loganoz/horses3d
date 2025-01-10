#include "Includes.h"

#if defined(NAVIERSTOKES) || defined(INCNS)
module HyperbolicSplitForm
   use SMConstants
#if defined(SPALARTALMARAS)
   use RiemannSolvers_NSSA
#elif defined(NAVIERSTOKES)
   use RiemannSolvers_NS
#elif defined(INCNS)
   use RiemannSolvers_iNS
#endif
   use HyperbolicDiscretizationClass
   use FluidData
   implicit none

   private
   public SplitDG_t, SplitDG_ComputeSplitFormFluxes

   type, extends(HyperbolicDiscretization_t) :: SplitDG_t
      contains
         procedure   :: Initialize             => SplitDG_Initialize
         !procedure   :: ComputeSplitFormFluxes => SplitDG_ComputeSplitFormFluxes
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
!        Setup the Riemann solver and two-point flux
!        -------------------------------------------
         call SetRiemannSolver(controlVariables)
!
!        Describe
!        --------
         if (.not. MPI_Process % isRoot ) return

         call Subsection_Header("Hyperbolic discretization")

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Split-Form"

         call DescribeRiemannSolver

         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Disabled."
         end if

      end subroutine SplitDG_Initialize

      subroutine SplitDG_ComputeSplitFormFluxes(e, contravariantFlux, fSharp, gSharp, hSharp)
         !$acc routine vector
         use ElementClass
         use PhysicsStorage
         implicit none
         type(Element),    intent(in)  :: e
         real(kind=RP),    intent(in)  :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP),    intent(out) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),    intent(out) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
         real(kind=RP),    intent(out) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0: e%Nxyz(3) )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l
!
!        First, diagonal results are introduced directly (consistency property)
!        ----------------------------------------------------------------------
         !$acc loop vector collapse(3)
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            fSharp(:,i,i,j,k) = contravariantFlux(:,i,j,k,IX)
            gSharp(:,j,i,j,k) = contravariantFlux(:,i,j,k,IY)
            hSharp(:,k,i,j,k) = contravariantFlux(:,i,j,k,IZ)
         end do               ; end do             ; end do
!
!        Then, terms out of the diagonal are computed
!        --------------------------------------------
         !$acc loop vector collapse(3)
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            !$acc loop seq
            do l = i+1, e%Nxyz(1)
               call TwoPointFlux_Selector(e % storage % Q(:,i,j,k), e % storage % Q(:,l,j,k), e % geom % jGradXi(:,i,j,k), &
                                          e % geom % jGradXi(:,l,j,k), fSharp(:,l,i,j,k))
               fSharp(:,i,l,j,k) = fSharp(:,l,i,j,k)
            end do
         end do               ; end do             ; end do
         
         !$acc loop vector collapse(3)
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            !$acc loop seq
            do l = j+1, e%Nxyz(2)
               call TwoPointFlux_Selector(e % storage % Q(:,i,j,k), e % storage % Q(:,i,l,k), e % geom % jGradEta(:,i,j,k), &
                                          e % geom % jGradEta(:,i,l,k), gSharp(:,l,i,j,k))
               gSharp(:,j,i,l,k) = gSharp(:,l,i,j,k)
            end do
         end do               ; end do             ; end do
         
         !$acc loop vector collapse(3)
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            !$acc loop seq
            do l = k+1, e%Nxyz(3)
               call TwoPointFlux_Selector(e % storage % Q(:,i,j,k), e % storage % Q(:,i,j,l), e % geom % jGradZeta(:,i,j,k), &
                                          e % geom % jGradZeta(:,i,j,l), hSharp(:,l,i,j,k))
               hSharp(:,k,i,j,l) = hSharp(:,l,i,j,k)
            end do
         end do               ; end do             ; end do

      end subroutine SplitDG_ComputeSplitFormFluxes
end module HyperbolicSplitForm
#endif