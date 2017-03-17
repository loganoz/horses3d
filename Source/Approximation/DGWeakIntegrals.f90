module DGWeakIntegrals
   use SMConstants
   use ElementClass
   use NodalStorageClass
   use PhysicsStorage
   use Physics
   implicit none


   private

   
   type  WeakIntegrals_t
      contains
         procedure, nopass    :: VolumeGreen  => DGWeakIntegrals_VolumeGreen
         procedure, nopass    :: FaceIntegral => DGWeakIntegrals_FaceIntegral
   end type WeakIntegrals_t

   public WeakIntegrals_t
!
!  ========
   contains
!  ========
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      function DGWeakIntegrals_VolumeGreen( e , spA , F ) result ( volInt )
         use MatrixOperations
         implicit none
         class(Element)          :: e
         class(NodalStorage)     :: spA
         real(kind=RP)           :: F      (0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN , 1:NDIM )
         real(kind=RP)           :: volInt (0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN          )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: N 

         N = spA % N

         volInt =  MatrixMultiplyInIndex_F(F(:,:,:,:,IX) , spA % hatD , N+1 , N+1 , N+1 , N_EQN , IX) &
                 + MatrixMultiplyInIndex_F(F(:,:,:,:,IY) , spA % hatD , N+1 , N+1 , N+1 , N_EQN , IY) &
                 + MatrixMultiplyInIndex_F(F(:,:,:,:,IZ) , spA % hatD , N+1 , N+1 , N+1 , N_EQN , IZ)

      end function DGWeakIntegrals_VolumeGreen

      function DGWeakIntegrals_FaceIntegral( e , spA , loc , F ) result ( faceInt )
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

      end function DGWeakIntegrals_FaceIntegral

end module DGWeakIntegrals
