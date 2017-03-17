module DGViscousDiscretization
!
!  *******
   private
!  *******
!
!
!  *****************************
!  Generic viscous methods class
!  *****************************
!
   type ViscousMethod_t
      contains
!         procedure      :: ComputeGradients      => BaseClass_ComputeGradients
         procedure      :: ComputeInnerFluxes    => BaseClass_ComputeInnerFluxes
   end type ViscousMethod_t

   type, extends(ViscousMethod_t)   :: BassiRebay1_t
      contains
!         procedure      :: ComputeGradients      => BassiRebay1_ComputeGradients
         procedure      :: ComputeInnerFluxes    => BassiRebay1_ComputeInnerFluxes

   end type BassiRebay1_t

   type, extends(ViscousMethod_t)   :: InteriorPenalty_t
      contains
!         procedure      :: ComputeGradients     => IP_ComputeGradients
!         procedure      :: ComputeInnerFluxes   => IP_ComputeInnerFluxes
   end type InteriorPenalty_t
!
!  ****************************************
   public   ViscousMethod_t , BassiRebay1_t
!  ****************************************
!

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
      subroutine BaseClass_ComputeInnerFluxes( self , e , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t) ,  intent (in)   :: self
         type(Element)                           :: e
         type(NodalStorage)      ,  intent (in)  :: spA
         real(kind=RP)           ,  intent (out) :: contravariantFlux(0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN, 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end subroutine BaseClass_ComputeInnerFluxes

!///////////////////////////////////////////////////////////////////////////////////
!
!           Bassi-Rebay 1 Procedures
!           ------------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine BassiRebay1_ComputeInnerFluxes( self , e , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         type(Element)                          :: e
         type(NodalStorage)      , intent (in)  :: spA
         real(kind=RP)           , intent (out) :: contravariantFlux(0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: cartesianFlux(0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN , 1:NDIM)
         integer             :: nv

         cartesianFlux = ViscousFlux( spA % N , e % Q , e % U_x , e % U_y , e % U_z )

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

      end subroutine BassiRebay1_ComputeInnerFluxes

end module DGViscousDiscretization
