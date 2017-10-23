module DGInviscidDiscretization
   use SMConstants
   implicit none

   private
   public   InviscidMethod_t , StandardDG_t , InviscidMethod


   type InviscidMethod_t
      contains
         procedure   :: ComputeInnerFluxes   => BaseClass_ComputeInnerFluxes
   end type InviscidMethod_t

   type, extends(InviscidMethod_t)  :: StandardDG_t
      contains
         procedure   :: ComputeInnerFluxes   => StandardDG_ComputeInnerFluxes
   end type StandardDG_t

   class(InviscidMethod_t), allocatable         :: InviscidMethod
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           BaseClass Procedures
!           --------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeInnerFluxes( self , e  , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         class(InviscidMethod_t) :: self
         type(Element)           :: e
         type(NodalStorage)      :: spA
         real(kind=RP)           :: contravariantFlux(1:NCONS, 0:spA % Nx , 0:spA % Ny , 0:spA % Nz, 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end subroutine BaseClass_ComputeInnerFluxes
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           StandardDG Procedures
!           ---------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine StandardDG_ComputeInnerFluxes( self , e , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(StandardDG_t) :: self
         type(Element)       :: e
         type(NodalStorage) :: spA
         real(kind=RP)      :: contravariantFlux(1:NCONS, 0:spA % Nx , 0:spA % Ny , 0:spA % Nz, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: i, j, k
         real(kind=RP)      :: cartesianFlux(1:NCONS, 0:spA % Nx , 0:spA % Ny , 0:spA % Nz, 1:NDIM)

         cartesianFlux = InviscidFlux( spA % Nx , spA % Ny , spA % Nz , e % storage % Q )

         do k = 0, spA % Nz   ; do j = 0, spA % Ny    ; do i = 0, spA % Nx
         
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

      end subroutine StandardDG_ComputeInnerFluxes

end module DGInviscidDiscretization
