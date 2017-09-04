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
         procedure, nopass    :: StdFace => VectorWeakIntegrals_StdFace
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
         integer, intent(in)              :: NEQ
         class(Element),      intent(in)  :: e
         class(NodalStorage), intent(in)  :: spA
         real(kind=RP),       intent(in)  :: F      ( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , 1:NEQ , 1:NDIM )
         real(kind=RP)                    :: volInt ( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , 1:NEQ          )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: Nx, Ny, Nz

         Nx = spA % Nx
         Ny = spA % Ny
         Nz = spA % Nz

         volInt =  MatrixMultiplyInIndex_F(F(:,:,:,:,IX) , spA % hatDx , Nx+1 , Ny+1 , Nz+1 , NEQ , IX) &
                 + MatrixMultiplyInIndex_F(F(:,:,:,:,IY) , spA % hatDy , Nx+1 , Ny+1 , Nz+1 , NEQ , IY) &
                 + MatrixMultiplyInIndex_F(F(:,:,:,:,IZ) , spA % hatDz , Nx+1 , Ny+1 , Nz+1 , NEQ , IZ)

      end function ScalarWeakIntegrals_StdVolumeGreen
!
!///////////////////////////////////////////////////////////////
!
      pure function ScalarWeakIntegrals_StdFace( e , spA , loc , F ) result ( faceInt )
         use MatrixOperations
         implicit none
         class(Element),      intent(in)     :: e
         class(NodalStorage), intent(in)     :: spA 
         integer,             intent(in)     :: loc
         real(kind=RP),       intent(in)     :: F(1:,0:,0:)
         real(kind=RP)                       :: faceInt(0:spA % Nx,0:spA % Ny,0:spA % Nz,1:N_EQN)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: iXi , iEta , iZeta , iVar

         select case (loc)
!
!     ----------------
! >   Xi-contributions
!     ----------------
!
            case (ELEFT)
         
               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iEta , iZeta) * spA % bx(iXi , LEFT)
               end do              ; end do                 ; end do                ; end do
            
            case (ERIGHT)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iEta , iZeta) * spA % bx(iXi , RIGHT)
               end do              ; end do                 ; end do                ; end do
!
!     -----------------
! >   Eta-contributions
!     -----------------
!
            case (EFRONT)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iZeta) * spA % by(iEta , LEFT)
               end do              ; end do                 ; end do                ; end do

            case (EBACK)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iZeta) * spA % by(iEta , RIGHT)
               end do              ; end do                 ; end do                ; end do
!
!     ------------------
!  >  Zeta-contributions
!     ------------------
!
            case (EBOTTOM)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % Nx ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iEta) * spA % bz(iZeta , LEFT)
               end do              ; end do                 ; end do                ; end do

            case (ETOP)

               do iVar = 1 , N_EQN ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt(iXi,iEta,iZeta,iVar) = F(iVar , iXi , iEta) * spA % bz(iZeta , RIGHT)
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
         integer,             intent(in)  :: NEQ
         class(Element),      intent(in)  :: e
         class(NodalStorage), intent(in)  :: spA
         real(kind=RP),       intent(in)  :: U        ( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , NEQ )
         real(kind=RP),       intent(out) :: volInt_x ( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , NEQ )
         real(kind=RP),       intent(out) :: volInt_y ( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , NEQ )
         real(kind=RP),       intent(out) :: volInt_z ( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , NEQ )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: eq , dimID
         real(kind=RP)               :: contravariantU( 0 : spA % Nx , 0 : spA % Ny , 0 : spA % Nz , NEQ , NDIM )
!!
!!      ╔═════════╗
!!      ║ Get U_x ║
!!      ╚═════════╝
!!
         do eq = 1 , NEQ
            contravariantU ( :,:,:,eq,IX ) = U ( :,:,:,eq ) * e % geom % jGradXi   ( IX,:,:,: ) 
            contravariantU ( :,:,:,eq,IY ) = U ( :,:,:,eq ) * e % geom % jGradEta  ( IX,:,:,: ) 
            contravariantU ( :,:,:,eq,IZ ) = U ( :,:,:,eq ) * e % geom % jGradZeta ( IX,:,:,: ) 
         end do

         volInt_x = ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , contravariantU )
!!
!!      ╔═════════╗
!!      ║ Get U_y ║
!!      ╚═════════╝
!!
         do eq = 1 , NEQ
            contravariantU ( :,:,:,eq,IX ) = U ( :,:,:,eq ) * e % geom % jGradXi   ( IY,:,:,: ) 
            contravariantU ( :,:,:,eq,IY ) = U ( :,:,:,eq ) * e % geom % jGradEta  ( IY,:,:,: ) 
            contravariantU ( :,:,:,eq,IZ ) = U ( :,:,:,eq ) * e % geom % jGradZeta ( IY,:,:,: ) 
         end do
!!
!!      ╔═════════╗
!!      ║ Get U_z ║
!!      ╚═════════╝
!!
         volInt_y = ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , contravariantU )

         do eq = 1 , NEQ
            contravariantU ( :,:,:,eq,IX ) = U ( :,:,:,eq ) * e % geom % jGradXi   ( IZ,:,:,: ) 
            contravariantU ( :,:,:,eq,IY ) = U ( :,:,:,eq ) * e % geom % jGradEta  ( IZ,:,:,: ) 
            contravariantU ( :,:,:,eq,IZ ) = U ( :,:,:,eq ) * e % geom % jGradZeta ( IZ,:,:,: ) 
         end do

         volInt_z = ScalarWeakIntegrals_StdVolumeGreen( NEQ , e , spA , contravariantU )
   
      end subroutine VectorWeakIntegrals_StdVolumeGreen
!
!/////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdFace( NEQ , e , spA , loc , U , faceInt_x , faceInt_y , faceInt_z )
         use ElementClass
         use NodalStorageClass
         use Physics
         use PhysicsStorage
         implicit none
         integer,             intent(in)  :: NEQ
         class(Element),      intent(in)  :: e
         class(NodalStorage), intent(in)  :: spA
         integer,             intent(in)  :: loc
         real(kind=RP),       intent(in)  :: U( 1: , 0: , 0: , 1: )                                   !<  
         real(kind=RP),       intent(out) :: faceInt_x(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , NEQ )   !>
         real(kind=RP),       intent(out) :: faceInt_y(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , NEQ )   !>
         real(kind=RP),       intent(out) :: faceInt_z(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , NEQ )   !>
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iVar , iXi , iEta , iZeta

         select case ( loc ) 
!
!     ----------------
! >   Xi-contributions
!     ----------------
!
            case (ELEFT)
         
               do iVar = 1 , NEQ ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt_x(iXi,iEta,iZeta,iVar) =   U(iVar , iEta , iZeta , ELEFT) * spA % bx(iXi , LEFT) &
                                                   * e % geom % normal(IX,iEta,iZeta,ELEFT) * e % geom % scal(iEta,iZeta,ELEFT)

                  faceInt_y(iXi,iEta,iZeta,iVar) =   U(iVar , iEta , iZeta , ELEFT) * spA % bx(iXi , LEFT) &
                                                   * e % geom % normal(IY,iEta,iZeta,ELEFT) * e % geom % scal(iEta,iZeta,ELEFT)

                  faceInt_z(iXi,iEta,iZeta,iVar) =   U(iVar , iEta , iZeta , ELEFT) * spA % bx(iXi , LEFT) &
                                                   * e % geom % normal(IZ,iEta,iZeta,ELEFT) * e % geom % scal(iEta,iZeta,ELEFT)
               end do              ; end do                 ; end do                ; end do
            
            case (ERIGHT)

               do iVar = 1 , NEQ ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt_x(iXi,iEta,iZeta,iVar) =   U(iVar , iEta , iZeta , ERIGHT) * spA % bx(iXi , RIGHT) &
                                                   * e % geom % normal(IX,iEta,iZeta,ERIGHT) * e % geom % scal(iEta,iZeta,ERIGHT)

                  faceInt_y(iXi,iEta,iZeta,iVar) =   U(iVar , iEta , iZeta , ERIGHT) * spA % bx(iXi , RIGHT) &
                                                   * e % geom % normal(IY,iEta,iZeta,ERIGHT) * e % geom % scal(iEta,iZeta,ERIGHT)

                  faceInt_z(iXi,iEta,iZeta,iVar) =   U(iVar , iEta , iZeta , ERIGHT) * spA % bx(iXi , RIGHT) &
                                                   * e % geom % normal(IZ,iEta,iZeta,ERIGHT) * e % geom % scal(iEta,iZeta,ERIGHT)
               end do              ; end do                 ; end do                ; end do
!
!     -----------------
! >   Eta-contributions
!     -----------------
!
            case (EFRONT)

               do iVar = 1 , NEQ ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt_x(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iZeta , EFRONT) * spA % by(iEta , LEFT) &
                                                   * e % geom % normal(IX,iXi,iZeta,EFRONT) * e % geom % scal(iXi,iZeta,EFRONT)

                  faceInt_y(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iZeta , EFRONT) * spA % by(iEta , LEFT) &
                                                   * e % geom % normal(IY,iXi,iZeta,EFRONT) * e % geom % scal(iXi,iZeta,EFRONT)

                  faceInt_z(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iZeta , EFRONT) * spA % by(iEta , LEFT) &
                                                   * e % geom % normal(IZ,iXi,iZeta,EFRONT) * e % geom % scal(iXi,iZeta,EFRONT)
               end do              ; end do                 ; end do                ; end do

            case (EBACK)

               do iVar = 1 , NEQ ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt_x(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iZeta , EBACK) * spA % by(iEta , RIGHT) &
                                                   * e % geom % normal(IX,iXi,iZeta,EBACK) * e % geom % scal(iXi,iZeta,EBACK)

                  faceInt_y(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iZeta , EBACK) * spA % by(iEta , RIGHT) &
                                                   * e % geom % normal(IY,iXi,iZeta,EBACK) * e % geom % scal(iXi,iZeta,EBACK)

                  faceInt_z(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iZeta , EBACK) * spA % by(iEta , RIGHT) &
                                                   * e % geom % normal(IZ,iXi,iZeta,EBACK) * e % geom % scal(iXi,iZeta,EBACK)
               end do              ; end do                 ; end do                ; end do
!
!     ------------------
!  >  Zeta-contributions
!     ------------------
!
            case (EBOTTOM)

               do iVar = 1 , NEQ ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt_x(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iEta , EBOTTOM) * spA % bz(iZeta , LEFT) &
                                                   * e % geom % normal(IX,iXi,iEta,EBOTTOM) * e % geom % scal(iXi,iEta,EBOTTOM)

                  faceInt_y(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iEta , EBOTTOM) * spA % bz(iZeta , LEFT) &
                                                   * e % geom % normal(IY,iXi,iEta,EBOTTOM) * e % geom % scal(iXi,iEta,EBOTTOM)

                  faceInt_z(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iEta , EBOTTOM) * spA % bz(iZeta , LEFT) &
                                                   * e % geom % normal(IZ,iXi,iEta,EBOTTOM) * e % geom % scal(iXi,iEta,EBOTTOM)
               end do              ; end do                 ; end do                ; end do

            case (ETOP)

               do iVar = 1 , NEQ ; do iZeta = 0 , spA % Nz ; do iEta = 0 , spA % Ny ; do iXi = 0 , spa % Nx
                  faceInt_x(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iEta , ETOP) * spA % bz(iZeta , RIGHT) &
                                                   * e % geom % normal(IX,iXi,iEta,ETOP) * e % geom % scal(iXi,iEta,ETOP)

                  faceInt_y(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iEta , ETOP) * spA % bz(iZeta , RIGHT) &
                                                   * e % geom % normal(IY,iXi,iEta,ETOP) * e % geom % scal(iXi,iEta,ETOP)

                  faceInt_z(iXi,iEta,iZeta,iVar) =   U(iVar , iXi , iEta , ETOP) * spA % bz(iZeta , RIGHT) &
                                                   * e % geom % normal(IZ,iXi,iEta,ETOP) * e % geom % scal(iXi,iEta,ETOP)
               end do              ; end do                 ; end do                ; end do

         end select

      end subroutine VectorWeakIntegrals_StdFace

end module DGWeakIntegrals
