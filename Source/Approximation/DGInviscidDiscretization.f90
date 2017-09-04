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
         real(kind=RP)           :: contravariantFlux(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1:N_EQN , 1:NDIM)
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
         real(kind=RP)      :: contravariantFlux(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1:N_EQN , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)      :: cartesianFlux(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1:N_EQN , 1:NDIM)
         integer             :: nv

         cartesianFlux = InviscidFlux( spA % Nx , spA % Ny , spA % Nz , e % Q )

         do nv = 1 , N_EQN
         
            contravariantFlux(:,:,:,nv,IX) =    cartesianFlux(:,:,:,nv,IX) * e % geom % jGradXi(IX,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IY) * e % geom % jGradXi(IY,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IZ) * e % geom % jGradXi(IZ,:,:,:)


            contravariantFlux(:,:,:,nv,IY) =    cartesianFlux(:,:,:,nv,IX) * e % geom % jGradEta(IX,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IY) * e % geom % jGradEta(IY,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IZ) * e % geom % jGradEta(IZ,:,:,:)


            contravariantFlux(:,:,:,nv,IZ) =    cartesianFlux(:,:,:,nv,IX) * e % geom % jGradZeta(IX,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IY) * e % geom % jGradZeta(IY,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IZ) * e % geom % jGradZeta(IZ,:,:,:)

         end do

      end subroutine StandardDG_ComputeInnerFluxes

end module DGInviscidDiscretization
