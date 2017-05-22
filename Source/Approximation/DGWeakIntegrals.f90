module DGWeakIntegrals
   use SMConstants
   use ElementClass
   use NodalStorageClass
   use PhysicsStorage
   use Physics
   implicit none


   private

   
   type  ScalarWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => ScalarWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => ScalarWeakIntegrals_StdFace
   end type ScalarWeakIntegrals_t

   type  VectorWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => VectorWeakIntegrals_StdVolumeGreen
!         procedure, nopass    :: StdFace => VectorWeakIntegrals_StdFace
   end type VectorWeakIntegrals_t

   public ScalarWeakIntegrals_t , VectorWeakIntegrals_t
   public ScalarWeakIntegrals   , VectorWeakIntegrals
   
   type(ScalarWeakIntegrals_t)      :: ScalarWeakIntegrals
   type(VectorWeakIntegrals_t)      :: VectorWeakIntegrals
!
!  ========
   contains
!  ========
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with scalar test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      function ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , F ) result ( volInt )
         use MatrixOperations
         implicit none
         integer, intent(in)     :: NEQ
         class(Element)          :: e
         class(NodalStorage)     :: spA
         real(kind=RP)           :: F      (0:spA % N , 0:spA % N , 0:spA % N , 1:NEQ , 1:NDIM )
         real(kind=RP)           :: volInt (0:spA % N , 0:spA % N , 0:spA % N , 1:NEQ          )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: N 

         N = spA % N

         volInt =  MatrixMultiplyInIndex_F(F(:,:,:,:,IX) , spA % hatD , N+1 , N+1 , N+1 , NEQ , IX) &
                 + MatrixMultiplyInIndex_F(F(:,:,:,:,IY) , spA % hatD , N+1 , N+1 , N+1 , NEQ , IY) &
                 + MatrixMultiplyInIndex_F(F(:,:,:,:,IZ) , spA % hatD , N+1 , N+1 , N+1 , NEQ , IZ)

      end function ScalarWeakIntegrals_StdVolumeGreen

      function ScalarWeakIntegrals_StdFace( e , spA , loc , F ) result ( faceInt )
         use MatrixOperations
         implicit none
         class(Element)          :: e
         class(NodalStorage)     :: spA 
         integer                 :: loc
         real(kind=RP)           :: F( 1:N_EQN , 0:spA % N , 0:spA % N )
         real(kind=RP)           :: faceInt(0:spA % N,0:spA % N,0:spA % N,1:N_EQN)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: iXi , iEta , iZeta , iVar
         INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2

         select case (loc)
!
!     ----------------
! >   Xi-contributions
!     ----------------
!
            case (ELEFT)
         
               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % N ; do iEta = 0 , spA % N ; do iXi = 0 , spa % N
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iEta , iZeta) * spA % b(iXi , LEFT)
               end do              ; end do                 ; end do                ; end do
            
            case (ERIGHT)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % N ; do iEta = 0 , spA % N ; do iXi = 0 , spa % N
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iEta , iZeta) * spA % b(iXi , RIGHT)
               end do              ; end do                 ; end do                ; end do
!
!     -----------------
! >   Eta-contributions
!     -----------------
!
            case (EFRONT)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % N ; do iEta = 0 , spA % N ; do iXi = 0 , spa % N
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iZeta) * spA % b(iEta , LEFT)
               end do              ; end do                 ; end do                ; end do

            case (EBACK)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % N ; do iEta = 0 , spA % N ; do iXi = 0 , spa % N
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iZeta) * spA % b(iEta , RIGHT)
               end do              ; end do                 ; end do                ; end do
!
!     ------------------
!  >  Zeta-contributions
!     ------------------
!
            case (EBOTTOM)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % N ; do iEta = 0 , spA % N ; do iXi = 0 , spa % N
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iEta) * spA % b(iZeta , LEFT)
               end do              ; end do                 ; end do                ; end do

            case (ETOP)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % N ; do iEta = 0 , spA % N ; do iXi = 0 , spa % N
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iEta) * spA % b(iZeta , RIGHT)
               end do              ; end do                 ; end do                ; end do

         end select

      end function ScalarWeakIntegrals_StdFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with vector test function
!        ----------------------------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdVolumeGreen( NEQ , e , spA , U , volInt_x , volInt_y , volInt_z )
         use ElementClass
         use NodalStorageClass
         use Physics
         use PhysicsStorage
         implicit none
         integer, intent(in)        :: NEQ
         class(Element)             :: e
         class(NodalStorage)        :: spA
         real(kind=RP), intent(in)  :: U        ( 0 : spA % N , 0 : spA % N , 0 : spA % N , NEQ ) 
         real(kind=RP), intent(out) :: volInt_x ( 0 : spA % N , 0 : spA % N , 0 : spA % N , NEQ ) 
         real(kind=RP), intent(out) :: volInt_y ( 0 : spA % N , 0 : spA % N , 0 : spA % N , NEQ ) 
         real(kind=RP), intent(out) :: volInt_z ( 0 : spA % N , 0 : spA % N , 0 : spA % N , NEQ ) 
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: eq , dimID
         real(kind=RP)        :: contravariantU(0:spA % N,0:spA % N,0:spA % N , NEQ , NDIM )

         do eq = 1 , NEQ
            contravariantU ( :,:,:,eq,IX ) = U ( :,:,:,eq ) * e % geom % jGradXi   ( IX,:,:,: ) 
            contravariantU ( :,:,:,eq,IY ) = U ( :,:,:,eq ) * e % geom % jGradEta  ( IX,:,:,: ) 
            contravariantU ( :,:,:,eq,IZ ) = U ( :,:,:,eq ) * e % geom % jGradZeta ( IX,:,:,: ) 
         end do

         volInt_x = ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , contravariantU )

         do eq = 1 , NEQ
            contravariantU ( :,:,:,eq,IX ) = U ( :,:,:,eq ) * e % geom % jGradXi   ( IY,:,:,: ) 
            contravariantU ( :,:,:,eq,IY ) = U ( :,:,:,eq ) * e % geom % jGradEta  ( IY,:,:,: ) 
            contravariantU ( :,:,:,eq,IZ ) = U ( :,:,:,eq ) * e % geom % jGradZeta ( IY,:,:,: ) 
         end do

         volInt_y = ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , contravariantU )

         do eq = 1 , NEQ
            contravariantU ( :,:,:,eq,IX ) = U ( :,:,:,eq ) * e % geom % jGradXi   ( IZ,:,:,: ) 
            contravariantU ( :,:,:,eq,IY ) = U ( :,:,:,eq ) * e % geom % jGradEta  ( IZ,:,:,: ) 
            contravariantU ( :,:,:,eq,IZ ) = U ( :,:,:,eq ) * e % geom % jGradZeta ( IZ,:,:,: ) 
         end do

         volInt_z = ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , contravariantU )
   
      end subroutine VectorWeakIntegrals_StdVolumeGreen

end module DGWeakIntegrals
