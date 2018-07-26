!
!//////////////////////////////////////////////////////
!
!   @File:
!   @Author:  David Kopriva
!   @Created: Tue Jun 04 15:34:44 2008
!   @Last revision date: Wed Jul 25 17:15:33 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: d886ff7a7d37081df645692157131f3ecc98f761
!
!//////////////////////////////////////////////////////
!
!      Implements Algorithms:
!         Algorithm: 124: ElementClass (QuadElementClass)
!
!       The Quad Element class, Alg. 124. See Sec. 8.2.1.2. self has
!       been modified to add the association of a boundary name to an element
!       edge so that different boundary conditions can be applied to different
!       elements. The names of the boundaries (not necessarily the names of the
!       *boundary conditions* to be applied) are of length BC_STRING_LENGTH.
!       One will associate boundary conditions to boundaries in the routine
!       "ExternalState".
!
!       Modified 2D Code to move solution into element class. 5/14/15, 5:36 PM
!
!////////////////////////////////////////////////////////////////////////
!
      Module ElementClass
      USE SMConstants
      USE PolynomialInterpAndDerivsModule
      USE GaussQuadrature
      USE TransfiniteMapClass
      USE MappedGeometryClass
      USE MeshTypes
      USE ElementConnectivityDefinitions
      USE ConnectivityClass
      use StorageClass
      USE NodalStorageClass
      use PhysicsStorage
      use Physics
      use VariableConversion
      IMPLICIT NONE

      private
      public   Element, axisMap 
      public   DestructElement, PrintElement, SetElementBoundaryNames
      
      TYPE Element
         logical                         :: hasSharedFaces
         integer                         :: dir2D
         integer                         :: globDir(3)        ! If the global coordinate (GLOBAL) is aligned with the local coordinate (REFERENCE): globDir(GLOBAL) = REFERENCE 
         integer                         :: eID               ! ID of this element
         integer                         :: globID            ! globalID of the element
         integer                         :: offsetIO          ! Offset from the first element for IO
         INTEGER                         :: nodeIDs(8)
         integer                         :: faceIDs(6)
         integer                         :: faceSide(6)
         INTEGER, DIMENSION(3)           :: Nxyz              ! Polynomial orders in every direction (Nx,Ny,Nz)
         TYPE(MappedGeometry)            :: geom
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName(6)
         INTEGER                         :: NumberOfConnections(6)
         TYPE(Connectivity)              :: Connection(6)
         type(ElementStorage_t)          :: storage
         type(NodalStorage_t), pointer   :: spAxi
         type(NodalStorage_t), pointer   :: spAeta
         type(NodalStorage_t), pointer   :: spAzeta
         type(TransfiniteHexMap)         :: hexMap            ! High-order mapper
         contains
            procedure   :: Construct               => HexElement_Construct
            procedure   :: ConstructGeometry       => HexElement_ConstructGeometry
            procedure   :: FindPointWithCoords     => HexElement_FindPointWithCoords
            procedure   :: EvaluateSolutionAtPoint => HexElement_EvaluateSolutionAtPoint
            procedure   :: ProlongSolutionToFaces  => HexElement_ProlongSolutionToFaces
            procedure   :: ProlongGradientsToFaces => HexElement_ProlongGradientsToFaces
            procedure   :: ComputeLocalGradient    => HexElement_ComputeLocalGradient
#if defined(NAVIERSTOKES)
            procedure   :: ComputeRhoGradient      => HexElement_ComputeRhoGradient
#endif
            procedure   :: InterpolateSolution     => HexElement_InterpolateSolution
      END TYPE Element 
      
!
!     -------------------------------------------------------------------------
!!    axisMap gives the element local coordinate number for the two directions
!!    on each face. The coordinate numbers are given by (xi,eta,zeta) = (1,2,3).
!!    For instance, the two coordinate directions on Face 1 are (xi,zeta).
!     -------------------------------------------------------------------------
!
      INTEGER, DIMENSION(2,6) :: axisMap =                        &
                                 RESHAPE( (/1, 3,                 & ! Face 1 (x,z)
                                            1, 3,                 & ! Face 2 (x,z)
                                            1, 2,                 & ! Face 3 (x,y)
                                            2, 3,                 & ! Face 4 (y,z)
                                            1, 2,                 & ! Face 5 (x,y)
                                            2, 3/)                & ! Face 6 (y,z)
                                 ,(/2,6/))
            
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE HexElement_Construct( self, Nx, Ny, Nz, nodeIDs, eID, globID)
         IMPLICIT NONE
         class(Element)      :: self
         integer, intent(in) :: Nx, Ny, Nz  !<  Polynomial orders         
         integer, intent(in) :: nodeIDs(8)
         integer, intent(in) :: eID, globID
         
         self % eID                 = eID
         self % dir2D               = 0
         self % globDir             = 0
         self % globID              = globID
         self % nodeIDs             = nodeIDs
         self % Nxyz(1)             = Nx
         self % Nxyz(2)             = Ny
         self % Nxyz(3)             = Nz
         self % boundaryName        = emptyBCName
         self % hasSharedFaces      = .false.
         self % NumberOfConnections = 0
         self % spAxi   => NodalStorage(Nx)
         self % spAeta  => NodalStorage(Ny)
         self % spAzeta => NodalStorage(Nz)
!
!        ----------------------------------------
!        Solution Storage is allocated separately
!        ----------------------------------------
!
      END SUBROUTINE HexElement_Construct
!
!////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------
!     Constructs the mapped geometry of the element (metric terms)
!     ------------------------------------------------------------
      subroutine HexElement_ConstructGeometry( self, hexMap)
         implicit none
         !--------------------------------------
         class(Element)         , intent(inout) :: self
         TYPE(TransfiniteHexMap), intent(in)    :: hexMap
         !--------------------------------------
         
         self % hexMap = hexMap
         CALL self % geom % Construct( self % spAxi, self % spAeta, self % spAzeta, hexMap )
         
      end subroutine HexElement_ConstructGeometry
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetElementBoundaryNames( self, names ) 
         use Utilities, only: toLower
         IMPLICIT NONE
         TYPE(Element)                   :: self
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(6)
         INTEGER                         :: j
         
         DO j = 1, 6
            CALL toLower(names(j)) 
            self % boundaryName(j) = names(j)
         END DO  
      END SUBROUTINE SetElementBoundaryNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructElement( self )
         IMPLICIT NONE
         TYPE(Element) :: self
         
         CALL self % geom % Destruct
         call self % Storage % Destruct   
         
         nullify( self % spAxi   )
         nullify( self % spAeta  )
         nullify( self % spAzeta )     

      END SUBROUTINE DestructElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintElement( self, id )
         IMPLICIT NONE 
         TYPE(Element) :: self
         INTEGER      :: id
         PRINT *, id, self % nodeIDs
         PRINT *, "   ",self % boundaryName
      END SUBROUTINE PrintElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionStorageToUnit( self, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: self
         INTEGER       :: fUnit
         
         WRITE(funit) self % storage % Q
      
      END SUBROUTINE SaveSolutionStorageToUnit
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionFromUnit( self, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: self
         INTEGER       :: fUnit
         
         READ(funit) self % storage % Q
      
      END SUBROUTINE LoadSolutionFromUnit
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexElement_ProlongSolutionToFaces(self, nEqn, fFR, fBK, fBOT, fR, fT, fL)
         use FaceClass
         implicit none
         class(Element),   intent(in)  :: self
         integer,          intent(in)  :: nEqn
         class(Face),      intent(inout) :: fFR, fBK, fBOT, fR, fT, fL
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l, N(3)
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: QFR, QBK
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: QBOT, QT
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: QL, QR

         N = self % Nxyz
!
!        *************************
!        Prolong solution to faces
!        *************************
!
         QL   = 0.0_RP     ; QR   = 0.0_RP
         QFR  = 0.0_RP     ; QBK  = 0.0_RP
         QBOT = 0.0_RP     ; QT   = 0.0_RP
         
         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            QL  (:,j,k)= QL  (:,j,k)+ self % storage % Q(:,i,j,k)* self % spAxi % v  (i,LEFT  )
            QR  (:,j,k)= QR  (:,j,k)+ self % storage % Q(:,i,j,k)* self % spAxi % v  (i,RIGHT )
            QFR (:,i,k)= QFR (:,i,k)+ self % storage % Q(:,i,j,k)* self % spAeta % v (j,FRONT )
            QBK (:,i,k)= QBK (:,i,k)+ self % storage % Q(:,i,j,k)* self % spAeta % v (j,BACK  )
            QBOT(:,i,j)= QBOT(:,i,j)+ self % storage % Q(:,i,j,k)* self % spAzeta % v(k,BOTTOM)
            QT  (:,i,j)= QT  (:,i,j)+ self % storage % Q(:,i,j,k)* self % spAzeta % v(k,TOP   )
         end do                   ; end do                   ; end do

         
         call fL   % AdaptSolutionToFace(nEqn, N(2), N(3), QL   , self % faceSide(ELEFT  ))
         call fR   % AdaptSolutionToFace(nEqn, N(2), N(3), QR   , self % faceSide(ERIGHT ))
         call fFR  % AdaptSolutionToFace(nEqn, N(1), N(3), QFR  , self % faceSide(EFRONT ))
         call fBK  % AdaptSolutionToFace(nEqn, N(1), N(3), QBK  , self % faceSide(EBACK  ))
         call fBOT % AdaptSolutionToFace(nEqn, N(1), N(2), QBOT , self % faceSide(EBOTTOM))
         call fT   % AdaptSolutionToFace(nEqn, N(1), N(2), QT   , self % faceSide(ETOP   ))

      end subroutine HexElement_ProlongSolutionToFaces

      subroutine HexElement_ProlongGradientsToFaces(self, nGradEqn, fFR, fBK, fBOT, fR, fT, fL)
         use FaceClass
         implicit none
         class(Element),   intent(in)  :: self
         integer,          intent(in)  :: nGradEqn
         class(Face),      intent(inout) :: fFR, fBK, fBOT, fR, fT, fL
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l, N(3)
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: UxFR, UyFR, UzFR
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: UxBK, UyBK, UzBK
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: UxBT, UyBT, UzBT
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: UxT, UyT, UzT
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: UxL, UyL, UzL
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: UxR, UyR, UzR

         N = self % Nxyz
!
!        *************************
!        Prolong solution to faces
!        *************************
!
         UxL  = 0.0_RP ; UyL  = 0.0_RP ; UzL  = 0.0_RP
         UxR  = 0.0_RP ; UyR  = 0.0_RP ; UzR  = 0.0_RP
         UxFR = 0.0_RP ; UyFR = 0.0_RP ; UzFR = 0.0_RP
         UxBK = 0.0_RP ; UyBK = 0.0_RP ; UzBK = 0.0_RP
         UxBT = 0.0_RP ; UyBT = 0.0_RP ; UzBT = 0.0_RP
         UxT  = 0.0_RP ; UyT  = 0.0_RP ; UzT  = 0.0_RP
         
         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            UxL (:,j,k) = UxL (:,j,k) + self % storage % U_x(:,i,j,k)* self % spAxi   % v (i,LEFT  )
            UxR (:,j,k) = UxR (:,j,k) + self % storage % U_x(:,i,j,k)* self % spAxi   % v (i,RIGHT )
            UxFR(:,i,k) = UxFR(:,i,k) + self % storage % U_x(:,i,j,k)* self % spAeta  % v (j,FRONT )
            UxBK(:,i,k) = UxBK(:,i,k) + self % storage % U_x(:,i,j,k)* self % spAeta  % v (j,BACK  )
            UxBT(:,i,j) = UxBT(:,i,j) + self % storage % U_x(:,i,j,k)* self % spAzeta % v (k,BOTTOM)
            UxT (:,i,j) = UxT (:,i,j) + self % storage % U_x(:,i,j,k)* self % spAzeta % v (k,TOP   )

            UyL (:,j,k) = UyL (:,j,k) + self % storage % U_y(:,i,j,k)* self % spAxi   % v (i,LEFT  )
            UyR (:,j,k) = UyR (:,j,k) + self % storage % U_y(:,i,j,k)* self % spAxi   % v (i,RIGHT )
            UyFR(:,i,k) = UyFR(:,i,k) + self % storage % U_y(:,i,j,k)* self % spAeta  % v (j,FRONT )
            UyBK(:,i,k) = UyBK(:,i,k) + self % storage % U_y(:,i,j,k)* self % spAeta  % v (j,BACK  )
            UyBT(:,i,j) = UyBT(:,i,j) + self % storage % U_y(:,i,j,k)* self % spAzeta % v (k,BOTTOM)
            UyT (:,i,j) = UyT (:,i,j) + self % storage % U_y(:,i,j,k)* self % spAzeta % v (k,TOP   )

            UzL (:,j,k) = UzL (:,j,k) + self % storage % U_z(:,i,j,k)* self % spAxi   % v (i,LEFT  )
            UzR (:,j,k) = UzR (:,j,k) + self % storage % U_z(:,i,j,k)* self % spAxi   % v (i,RIGHT )
            UzFR(:,i,k) = UzFR(:,i,k) + self % storage % U_z(:,i,j,k)* self % spAeta  % v (j,FRONT )
            UzBK(:,i,k) = UzBK(:,i,k) + self % storage % U_z(:,i,j,k)* self % spAeta  % v (j,BACK  )
            UzBT(:,i,j) = UzBT(:,i,j) + self % storage % U_z(:,i,j,k)* self % spAzeta % v (k,BOTTOM)
            UzT (:,i,j) = UzT (:,i,j) + self % storage % U_z(:,i,j,k)* self % spAzeta % v (k,TOP   )

         end do                   ; end do                   ; end do
         
         call fL   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxL , UyL , UzL , self % faceSide(ELEFT  ))
         call fR   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxR , UyR , UzR , self % faceSide(ERIGHT ))
         call fFR  % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxFR, UyFR, UzFR, self % faceSide(EFRONT ))
         call fBK  % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxBK, UyBK, UzBK, self % faceSide(EBACK  ))
         call fBOT % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxBT, UyBT, UzBT, self % faceSide(EBOTTOM))
         call fT   % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxT , UyT , UzT , self % faceSide(ETOP   ))

      end subroutine HexElement_ProlongGradientsToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexElement_ComputeLocalGradient(self, nEqn, nGradEqn, GetGradientValues3D)
!
!        ****************************************************************
!           This subroutine computes local gradients as:
!  
!              nabla U = (1/J) \sum_i Ja^i \cdot \partial(u)/ \partial \xi^i
!
!        ****************************************************************
!  
         implicit none
         class(Element),   intent(inout)  :: self
         integer,          intent(in)     :: nEqn
         integer,          intent(in)     :: nGradEqn
         procedure(GetGradientValues3D_f) :: GetGradientValues3D
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l
         real(kind=RP)  :: U(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: U_xi(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: U_eta(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: U_zeta(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))

         associate( N => self % Nxyz )
!
!        **********************
!        Get gradient variables
!        **********************
!
         call GetGradientValues3D(nEqn, nGradEqn, N(1), N(2), N(3), self % storage % Q, U )
!
!        ************
!        Compute U_xi
!        ************
!
         U_xi = 0.0_RP
         do k = 0, N(3)   ; do j = 0, N(2) ; do l = 0, N(1) ; do i = 0, N(1)
            U_xi(:,i,j,k) = U_xi(:,i,j,k) + U(:,l,j,k) * self % spAxi % D(i,l)
         end do           ; end do         ; end do         ; end do
!
!        *************
!        Compute U_eta
!        *************
!
         U_eta = 0.0_RP
         do k = 0, N(3)   ; do l = 0, N(2) ; do j = 0, N(2) ; do i = 0, N(1)
            U_eta(:,i,j,k) = U_eta(:,i,j,k) + U(:,i,l,k) * self % spAeta % D(j,l)
         end do           ; end do         ; end do         ; end do
!
!        **************
!        Compute U_zeta
!        **************
!
         U_zeta = 0.0_RP
         do l = 0, N(3)   ; do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            U_zeta(:,i,j,k) = U_zeta(:,i,j,k) + U(:,i,j,l) * self % spAzeta % D(k,l)
         end do           ; end do         ; end do         ; end do
!
!        ******************************
!        Project to the cartesian basis
!        ******************************
!
         do k = 0, N(3) ; do j = 0, N(2)  ; do i = 0, N(1)
            self % storage % U_x(:,i,j,k) = (   U_xi(:,i,j,k) * self % geom % jGradXi(1,i,j,k) &
                                              + U_eta(:,i,j,k) * self % geom % jGradEta(1,i,j,k) & 
                                              + U_zeta(:,i,j,k) * self % geom % jGradZeta(1,i,j,k))&
                                             * self % geom % InvJacobian(i,j,k)

            self % storage % U_y(:,i,j,k) = (    U_xi(:,i,j,k) * self % geom % jGradXi(2,i,j,k) &
                                              + U_eta(:,i,j,k) * self % geom % jGradEta(2,i,j,k) & 
                                              + U_zeta(:,i,j,k) * self % geom % jGradZeta(2,i,j,k))& 
                                            * self % geom % InvJacobian(i,j,k)

            self % storage % U_z(:,i,j,k) = (    U_xi(:,i,j,k) * self % geom % jGradXi(3,i,j,k) &
                                              + U_eta(:,i,j,k) * self % geom % jGradEta(3,i,j,k) & 
                                              + U_zeta(:,i,j,k) * self % geom % jGradZeta(3,i,j,k))& 
                                            * self % geom % InvJacobian(i,j,k)
         end do         ; end do          ; end do
         end associate

      end subroutine HexElement_ComputeLocalGradient
!
!////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)  
      subroutine HexElement_ComputeRhoGradient(self)
!
!        ****************************************************************
!           This subroutine computes local gradients as:
!  
!              nabla U = (1/J) \sum_i Ja^i \cdot \partial(u)/ \partial \xi^i
!
!        ****************************************************************
!  
         implicit none
         class(Element),   intent(inout)  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l
         real(kind=RP)  :: rho      (0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: rho_xi   (0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: rho_eta  (0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: rho_zeta (0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))

         associate( N => self % Nxyz )
!
!        **********************
!        Get density
!        **********************
!
         rho = self % storage % Q (1,:,:,:)
!
!        **************
!        Compute rho_xi
!        **************
!
         rho_xi = 0.0_RP
         do k = 0, N(3)   ; do j = 0, N(2) ; do l = 0, N(1) ; do i = 0, N(1)
            rho_xi(i,j,k) = rho_xi(i,j,k) + rho(l,j,k) * self % spAxi % D(i,l)
         end do           ; end do         ; end do         ; end do
!
!        ***************
!        Compute rho_eta
!        ***************
!
         rho_eta = 0.0_RP
         do k = 0, N(3)   ; do l = 0, N(2) ; do j = 0, N(2) ; do i = 0, N(1)
            rho_eta(i,j,k) = rho_eta(i,j,k) + rho(i,l,k) * self % spAeta % D(j,l)
         end do           ; end do         ; end do         ; end do

!        ****************
!        Compute rho_zeta
!        ****************
!
         rho_zeta = 0.0_RP
         do l = 0, N(3)   ; do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            rho_zeta(i,j,k) = rho_zeta(i,j,k) + rho(i,j,l) * self % spAzeta % D(k,l)
         end do           ; end do         ; end do         ; end do
!
!        ******************************
!        Project to the cartesian basis
!        ******************************
!
         do k = 0, N(3) ; do j = 0, N(2)  ; do i = 0, N(1)
            self % storage % gradRho(1,i,j,k) = (  rho_xi(i,j,k) * self % geom % jGradXi(1,i,j,k) &
                                                 + rho_eta(i,j,k) * self % geom % jGradEta(1,i,j,k) & 
                                                 + rho_zeta(i,j,k) * self % geom % jGradZeta(1,i,j,k))&
                                                * self % geom % InvJacobian(i,j,k)

            self % storage % gradRho(2,i,j,k) = (  rho_xi(i,j,k) * self % geom % jGradXi(2,i,j,k) &
                                                 + rho_eta(i,j,k) * self % geom % jGradEta(2,i,j,k) & 
                                                 + rho_zeta(i,j,k) * self % geom % jGradZeta(2,i,j,k))& 
                                                * self % geom % InvJacobian(i,j,k)
   
            self % storage % gradRho(3,i,j,k) = (  rho_xi(i,j,k) * self % geom % jGradXi(3,i,j,k) &
                                                 + rho_eta(i,j,k) * self % geom % jGradEta(3,i,j,k) & 
                                                 + rho_zeta(i,j,k) * self % geom % jGradZeta(3,i,j,k))& 
                                               * self % geom % InvJacobian(i,j,k)
         end do         ; end do          ; end do
         end associate

      end subroutine HexElement_ComputeRhoGradient
#endif
!
!////////////////////////////////////////////////////////////////////////
!
      logical function HexElement_FindPointWithCoords(self, x, xi)
!
!        **********************************************************
!          
!           This function finds whether a point is inside or not 
!           of the element. This is done solving
!           the mapping non-linear system
!
!        **********************************************************
!          
!
         use Utilities, only: SolveThreeEquationLinearSystem
         implicit none
         class(Element),      intent(in)  :: self
         real(kind=RP),       intent(in)  :: x(NDIM)
         real(kind=RP),       intent(out) :: xi(NDIM)
!
!        ----------------------------------
!        Newton iterative solver parameters
!        ----------------------------------
!
         integer,       parameter   :: N_MAX_ITER = 50
         real(kind=RP), parameter   :: TOL = 1.0e-12_RP
         integer,       parameter   :: STEP = 1.0_RP
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j, k, iter
         real(kind=RP), parameter      :: INSIDE_TOL = 1.0e-08_RP
         real(kind=RP)                 :: lxi   (0:self % Nxyz(1)) 
         real(kind=RP)                 :: leta  (0:self % Nxyz(2)) 
         real(kind=RP)                 :: lzeta (0:self % Nxyz(3)) 
         real(kind=RP)                 :: dlxi   (0:self % Nxyz(1)) 
         real(kind=RP)                 :: dleta  (0:self % Nxyz(2)) 
         real(kind=RP)                 :: dlzeta (0:self % Nxyz(3)) 
         real(kind=RP)                 :: F(NDIM)
         real(kind=RP)                 :: Jac(NDIM,NDIM)
         real(kind=RP)                 :: dx(NDIM)
!
!        Initial seed
!        ------------      
         xi = 0.0_RP    

         do iter = 1 , N_MAX_ITER
!
!           Get Lagrange polynomials and derivatives
!           ----------------------------------------
            lxi     = self %spAxi % lj   (xi(1))
            leta    = self %spAeta % lj  (xi(2))
            lzeta   = self %spAzeta % lj (xi(3))
  
            F = 0.0_RP
            do k = 0, self %spAzeta % N   ; do j = 0, self %spAeta % N ; do i = 0, self %spAxi % N
               F = F + self % geom % x(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do
   
            F = F - x
!
!           Stopping criteria: there are several
!           ------------------------------------
            if ( maxval(abs(F)) .lt. TOL ) exit
            if ( abs(xi(1)) .ge. 2.5_RP ) exit
            if ( abs(xi(2)) .ge. 2.5_RP ) exit
            if ( abs(xi(3)) .ge. 2.5_RP ) exit
!
!           Perform a step
!           --------------
            dlxi    = self %spAxi % dlj  (xi(1))
            dleta   = self %spAeta % dlj (xi(2))
            dlzeta  = self %spAzeta % dlj(xi(3))

            Jac = 0.0_RP
            do k = 0, self %spAzeta % N   ; do j = 0, self %spAeta % N ; do i = 0, self %spAxi % N
               Jac(:,1) = Jac(:,1) + self % geom % x(:,i,j,k) * dlxi(i) * leta(j) * lzeta(k) 
               Jac(:,2) = Jac(:,2) + self % geom % x(:,i,j,k) * lxi(i) * dleta(j) * lzeta(k) 
               Jac(:,3) = Jac(:,3) + self % geom % x(:,i,j,k) * lxi(i) * leta(j) * dlzeta(k) 
            end do               ; end do             ; end do

            dx = solveThreeEquationLinearSystem( Jac , -F )
            xi = xi + STEP * dx
   
         end do

         if ( (abs(xi(1)) .lt. 1.0_RP + INSIDE_TOL) .and. &
              (abs(xi(2)) .lt. 1.0_RP + INSIDE_TOL) .and. &
              (abs(xi(3)) .lt. 1.0_RP + INSIDE_TOL)          ) then
!
!           Solution is valid
!           -----------------
            HexElement_FindPointWithCoords = .true.
   
         else
!
!           Solution is not valid
!           ---------------------
            HexElement_FindPointWithCoords = .false.
         
         end if

      end function HexElement_FindPointWithCoords

      function HexElement_EvaluateSolutionAtPoint(self, nEqn, xi)
         implicit none
         class(Element),   intent(in)    :: self
         integer,          intent(in)    :: nEqn
         real(kind=RP),    intent(in)    :: xi(NDIM)
         real(kind=RP)                   :: HexElement_EvaluateSolutionAtPoint(nEqn)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: lxi(0:self % Nxyz(1))
         real(kind=RP)  :: leta(0:self % Nxyz(2))
         real(kind=RP)  :: lzeta(0:self % Nxyz(3))
         real(kind=RP)  :: Q(nEqn)
!
!        Compute Lagrange basis
!        ----------------------
         lxi   = self %spAxi % lj(xi(1))
         leta  = self %spAeta % lj(xi(2))
         lzeta = self %spAzeta % lj(xi(3))
!
!        Compute the tensor product
!        --------------------------
         Q = 0.0_RP
      
         do k = 0, self % spAzeta % N   ; do j = 0, self % spAeta % N ; do i = 0, self % spAxi % N
            Q = Q + self % storage % Q(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
         end do               ; end do             ; end do   

         HexElement_EvaluateSolutionAtPoint = Q

      end function HexElement_EvaluateSolutionAtPoint
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------
!     Interpolate solution to another element
!            this % storage % Q  ->  e % storage % Q
!     ----------------------------------------------
      impure elemental subroutine HexElement_InterpolateSolution(this,e,nodes,with_gradients)
         use InterpolationMatrices, only: ConstructInterpolationMatrices, Interp3DArrays
         implicit none
         !-arguments----------------------------------------------
         class(Element), intent(in)    :: this
         type(Element) , intent(inout) :: e
         integer       , intent(in)    :: nodes
         logical, optional, intent(in) :: with_gradients
         !-local-variables----------------------------------------
         logical                       :: gradients
         !--------------------------------------------------------
         
         if ( present(with_gradients) ) then
            gradients = with_gradients
         else
            gradients = .FALSE.
         end if
         
         ! Copy the solution if the polynomial orders are the same, if not, interpolate
         
         
         if (all(this % Nxyz == e % Nxyz)) then
            e % storage % Q = this % storage % Q
         else
            
            call NodalStorage(this % Nxyz(1)) % construct(nodes,this % Nxyz(1))
            call NodalStorage(this % Nxyz(2)) % construct(nodes,this % Nxyz(2))
            call NodalStorage(this % Nxyz(3)) % construct(nodes,this % Nxyz(3))
            call NodalStorage(e % Nxyz(1)) % construct(nodes,e % Nxyz(1))
            call NodalStorage(e % Nxyz(2)) % construct(nodes,e % Nxyz(2))
            call NodalStorage(e % Nxyz(3)) % construct(nodes,e % Nxyz(3))
            
            !------------------------------------------------------------------
            ! Construct the interpolation matrices in every direction if needed
            !------------------------------------------------------------------
            call ConstructInterpolationMatrices( this % Nxyz(1), e % Nxyz(1) )  ! Xi
            call ConstructInterpolationMatrices( this % Nxyz(2), e % Nxyz(2) )  ! Eta
            call ConstructInterpolationMatrices( this % Nxyz(3), e % Nxyz(3) )  ! Zeta
            
            !---------------------------------------------
            ! Interpolate solution to new solution storage
            !---------------------------------------------
            
            call Interp3DArrays  (Nvars      = NTOTALVARS         , &
                                  Nin        = this % Nxyz        , &
                                  inArray    = this % storage % Q , &
                                  Nout       = e % Nxyz           , &
                                  outArray   = e % storage % Q          )
            
            if (gradients) then
               call Interp3DArrays  (Nvars      = NTOTALGRADS     , &
                                  Nin        = this % Nxyz        , &
                                  inArray    = this % storage % U_x , &
                                  Nout       = e % Nxyz           , &
                                  outArray   = e % storage % U_x          )
               
               call Interp3DArrays  (Nvars      = NTOTALGRADS              , &
                                  Nin        = this % Nxyz        , &
                                  inArray    = this % storage % U_y , &
                                  Nout       = e % Nxyz           , &
                                  outArray   = e % storage % U_y        )
                                  
               call Interp3DArrays  (Nvars      = NTOTALGRADS              , &
                                  Nin        = this % Nxyz        , &
                                  inArray    = this % storage % U_z , &
                                  Nout       = e % Nxyz           , &
                                  outArray   = e % storage % U_z        )
            end if
            
         end if
      end subroutine HexElement_InterpolateSolution
      
      END Module ElementClass
