!!
!////////////////////////////////////////////////////////////////////////
!
!      Modification history:
!           Modified to 3D             5/27/15, 11:13 AM: David A. Kopriva
!           Added isotropic mortars    4/26/17, 11:12 AM: Andrés Rueda
!           Added anisotropic mortars  5/16/17, 11:11 AM: Andrés Rueda
!           Embedded mortars          11/13/17, 05:37 PM: Juan Manzanero
!
!////////////////////////////////////////////////////////////////////////
!
      Module FaceClass
      use SMConstants
      use MeshTypes
      use PolynomialInterpAndDerivsModule
      use GaussQuadrature
      use MappedGeometryClass
      use StorageClass                    , only: FaceStorage_t
      use PhysicsStorage
      use NodalStorageClass
      use InterpolationMatrices           , only: Tset
      IMPLICIT NONE 

      private
      public   Face, stencil_t

      type stencil_t

         real(kind=RP)              :: x(NDIM), normal(NDIM), dist, xiI, xiB, u, v, w, rho, p
         integer                    :: N, partition
         logical                    :: state = .false.
         real(kind=RP), allocatable :: x_s(:,:), xi_s(:,:), Q(:,:)
         integer,       allocatable :: eIDs(:)

         contains
            procedure :: build        => stencil_build
            procedure :: ComputeState => stencil_ComputeState 
            procedure :: destroy      => stencil_destroy

      end type stencil_t

      type IBM_HO_faces_t

         type(face), allocatable :: faces(:)
         integer                 :: NumOfFaces, NumOfObjs, NumOfDoFs

      end type IBM_HO_faces_t

      type(IBM_HO_faces_t), allocatable, public :: IBM_HO_faces(:), IBM_HO_mpifaces(:)
!
!     ************************************************************************************
!
!           Face derived type. The face connects two elements,
!           =================  specified in "elementIDs". 
!
!        The face handles two problems:
!
!           1) Rotation: faces are not oriented in the same direction.
!              --------  The direction adopted is always that of the left
!                        element.
!           2) Different polynomial order: elements at each sidecan have different
!              --------------------------  polynomial orders.
!
!        We have defined the following quantities:
!
!           -> NelLeft/Right: The face original order in the element, without
!              ~~~~~~~~~~~~~  rotation or projections.
!
!           -> NfLeft/Right: The face original order rotated to match each other,
!              ~~~~~~~~~~~~  without projection.
!                       ** Thus, always NfLeft = NelLeft, and
!                          NfRight = rotation(NelRight)                           
!
!           -> Nf: The face polynomial order. This is the maximum along 
!              ~~  each side of both directions.
!
!        The process from the element storage to face storage is:
!
!          (0:Nxyz)           (0:NelLeft)                (0:Nf)
!        Left element -- Interpolation to Boundary -- Projection to face
!
!          (0:Nxyz)           (0:NelRight)             (0:NfRight)       (0:Nf)
!        Right element -- Interpolation to Boundary --  Rotation -- Projection to face
!
!     ************************************************************************************
!
      type Face
         logical                         :: flat
         integer                         :: ID                       ! face ID
         integer                         :: FaceType                 ! Type of face: 1 = HMESH_INTERIOR, 2 = HMESH_BOUNDARY, 3 = 2 = HMESH_MPI
         integer                         :: zone                     ! In the case of HMESH_BOUNDARY, which zone it belongs to
         integer                         :: rotation                 ! Relative orientation between faces
         integer                         :: NelLeft(2)               ! Left element face polynomial order
         integer                         :: NelRight(2)              ! Right element face polynomial order
         integer                         :: NfLeft(2)                ! Left face polynomial order
         integer                         :: NfRight(2)               ! Right face polynomial order
         integer                         :: Nf(2)                    ! Face polynomial order
         integer                         :: nodeIDs(4)               ! Face nodes ID
         integer                         :: elementIDs(2)            ! Convention is: 1 = Left, 2 = Right
         integer                         :: elementSide(2)
         integer                         :: projectionType(2)
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
         type(MappedGeometryFace)        :: geom
         type(FaceStorage_t)             :: storage(2)
         logical                         :: HO_IBM = .false., mpiHO_IBM = .false., corrGrad = .false., fmpi =.false.
         integer                         :: HOSIDE, HO_ID, domain
         type(stencil_t), allocatable    :: stencil(:,:)
         contains
            procedure   :: Construct                     => ConstructFace
            procedure   :: Destruct                      => DestructFace
            procedure   :: Print                         => PrintFace
            procedure   :: LinkWithElements              => Face_LinkWithElements
            procedure   :: AdaptSolutionToFace           => Face_AdaptSolutionToFace
            procedure   :: AdaptGradientsToFace          => Face_AdaptGradientsToFace
            procedure   :: AdaptAviscFluxToFace          => Face_AdaptAviscFluxToFace
            procedure   :: ProjectFluxToElements         => Face_ProjectFluxToElements
            procedure   :: ProjectGradientFluxToElements => Face_ProjectGradientFluxToElements
#if defined(NAVIERSTOKES)
            procedure   :: ProjectFluxJacobianToElements => Face_ProjectFluxJacobianToElements
            procedure   :: ProjectGradJacobianToElements => Face_ProjectGradJacobianToElements
            procedure   :: ProjectBCJacobianToElements   => Face_ProjectBCJacobianToElements
#endif
            procedure   :: copy           => Face_Assign
            generic     :: assignment(=)  => copy
            procedure   :: HO_IBM_correction             => Face_HO_IBM_correction
            procedure   :: HO_IBM_Gradcorrection         => Face_HO_IBM_Gradcorrection
            procedure   :: StencilConstruct              => face_StencilConstruct
      end type Face
!
!     ========
      CONTAINS
!     ========
!
      SUBROUTINE ConstructFace( self, ID, nodeIDs, elementID, side )
!
!        *******************************************************
!
!           This is not the full face construction, but a 
!         variable initialization. Just nodes are introduced and
!         the left element.
!
!        *******************************************************
!
         IMPLICIT NONE 
         class(Face) :: self
         integer     :: ID, nodeIDs(4), elementID, side
!
!        --------------------------
!        Set nodes and left element 
!        --------------------------
!
         self % ID             = ID
         self % nodeIDS        = nodeIDs
         self % elementIDs(1)  = elementID
         self % elementIDs(2)  = -1
         self % elementSide(1) = side
!
!        ------------
!        Set defaults
!        ------------
!
         self % FaceType       = HMESH_UNDEFINED
         self % elementIDs(2)  = HMESH_NONE
         self % elementSide(2) = HMESH_NONE
         self % boundaryName   = emptyBCName
         self % rotation       = 0
         self % zone           = 0

      end SUBROUTINE ConstructFace
!
!////////////////////////////////////////////////////////////////////////
!
      elemental SUBROUTINE DestructFace( self )
         IMPLICIT NONE 
         class(Face), intent(inout) :: self
         
         self % ID = -1
         self % FaceType = HMESH_NONE
         self % rotation = 0
         self % NelLeft = -1
         self % NelRight = -1       
         self % NfLeft = -1       
         self % NfRight = -1       
         self % Nf = -1         
         self % nodeIDs = -1             
         self % elementIDs = -1
         self % elementSide = -1
         self % projectionType = -1
         self % boundaryName = ""
         
         call self % geom % Destruct
         call self % storage % Destruct

      end SUBROUTINE DestructFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintFace( self ) 
      IMPLICIT NONE
      class(Face) :: self
      PRINT *, "Face type = "   , self % FaceType
      PRINT *, "Element IDs: "  , self % elementIDs
      PRINT *, "Element Sides: ", self % elementSide
      IF ( self % FaceType == HMESH_INTERIOR )     THEN
         PRINT *, "Neighbor rotation: ", self  %  rotation
      ELSE
         PRINT *, "Boundary name = ", self % boundaryName
      end IF
      PRINT *, "-----------------------------------"
      end SUBROUTINE PrintFace
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Face_LinkWithElements( self, NelLeft, NelRight, nodeType)
      IMPLICIT NONE
      class(Face)        ,     intent(INOUT) :: self        ! Current face
      integer,                 intent(in)    :: NelLeft(2)  ! Left element face polynomial order
      integer,                 intent(in)    :: NelRight(2) ! Right element face polynomial order
      integer,                 intent(in)    :: nodeType    ! Either Gauss or Gauss-Lobatto
#if (!defined(NAVIERSTOKES)) && (!defined(INCNS))
      logical  :: computeGradients = .true.
#endif
!
!     -------------------------------------------------------------
!     First, get face elements polynomial orders (without rotation)
!     -------------------------------------------------------------
!
      self % NelLeft = NelLeft
      self % NelRight = NelRight
!
!     ---------------------------------------------------------
!     Second, get face polynomial orders (considering rotation)
!     ---------------------------------------------------------
!
      self % NfLeft = self % NelLeft     ! Left elements are always oriented.
      
      SELECT CASE ( self % rotation )
      CASE ( 0, 2, 5, 7 ) ! Local x and y axis are parallel or antiparallel
         self % NfRight = self % NelRight
      CASE ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
         self % NfRight(2)  = self % NelRight(1)
         self % NfRight(1)  = self % NelRight(2)
      end SELECT
!
!     ------------------------------------------------------------------
!     Third, the face polynomial order (the maximum for both directions)
!     ------------------------------------------------------------------
!
      self % Nf(1) = max(self % NfLeft(1), self % NfRight(1))
      self % Nf(2) = max(self % NfLeft(2), self % NfRight(2))
!
!     ------------------------------
!     Construct needed nodal storage
!     ------------------------------
!
      if ( .not. NodalStorage(self % Nf(1)) % Constructed ) then
         call NodalStorage(self % Nf(1)) % Construct(nodeType, self % Nf(1))
      end if

      if ( .not. NodalStorage(self % Nf(2)) % Constructed ) then
         call NodalStorage(self % Nf(2)) % Construct(nodeType, self % Nf(2))
      end if
!
!     -----------------------------------------------------------------------
!     Construction of the projection matrices (simple Lagrange interpolation)
!     -----------------------------------------------------------------------
!
      call Tset(self % NfLeft(1), self % Nf(1)) % construct(self % NfLeft(1), self % Nf(1))
      call Tset(self % Nf(1), self % NfLeft(1)) % construct(self % Nf(1), self % NfLeft(1))
      
      call Tset(self % NfLeft(2), self % Nf(2)) % construct(self % NfLeft(2), self % Nf(2))
      call Tset(self % Nf(2), self % NfLeft(2)) % construct(self % Nf(2), self % NfLeft(2))
      
      call Tset(self % NfRight(1), self % Nf(1)) % construct(self % NfRight(1), self % Nf(1))
      call Tset(self % Nf(1), self % NfRight(1)) % construct(self % Nf(1), self % NfRight(1))
      
      call Tset(self % NfRight(2), self % Nf(2)) % construct(self % NfRight(2), self % Nf(2))
      call Tset(self % Nf(2), self % NfRight(2)) % construct(self % Nf(2), self % NfRight(2))
!
!     -----------------------  0- no projection
!     Set the projection type: 1- x needs projection
!     -----------------------  2- y needs projection
!                              3- both x and y need projection
      if ( (self % NfLeft(1) .eq. self % Nf(1)) .and. (self % NfLeft(2) .eq. self % Nf(2)) ) then
         self % projectionType(1) = 0
      elseif ( (self % NfLeft(1) .ne. self % Nf(1)) .and. (self % NfLeft(2) .eq. self % Nf(2)) ) then
         self % projectionType(1) = 1
      elseif ( (self % NfLeft(1) .eq. self % Nf(1)) .and. (self % NfLeft(2) .ne. self % Nf(2)) ) then
         self % projectionType(1) = 2
      elseif ( (self % NfLeft(1) .ne. self % Nf(1)) .and. (self % NfLeft(2) .ne. self % Nf(2)) ) then
         self % projectionType(1) = 3
      end if

      if ( (self % NfRight(1) .eq. self % Nf(1)) .and. (self % NfRight(2) .eq. self % Nf(2)) ) then
         self % projectionType(2) = 0
      elseif ( (self % NfRight(1) .ne. self % Nf(1)) .and. (self % NfRight(2) .eq. self % Nf(2)) ) then
         self % projectionType(2) = 1
      elseif ( (self % NfRight(1) .eq. self % Nf(1)) .and. (self % NfRight(2) .ne. self % Nf(2)) ) then
         self % projectionType(2) = 2
      elseif ( (self % NfRight(1) .ne. self % Nf(1)) .and. (self % NfRight(2) .ne. self % Nf(2)) ) then
         self % projectionType(2) = 3
      end if

   end SUBROUTINE Face_LinkWithElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Face_AdaptSolutionToFace(self, nEqn, Nelx, Nely, Qe, side, QdotE, computeQdot)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      integer,       intent(in)                 :: Nelx, Nely
      real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)                 :: side
      real(kind=RP), intent(in), optional       :: QdotE(1:nEqn, 0:Nelx, 0:Nely)
      logical,       intent(in), optional       :: computeQdot
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: Qe_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      ! real(kind=RP) :: QdotE_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      logical :: prolongQdot

      ! prolongQdot = present(QdotE)
      if (present(computeQdot)) then
          prolongQdot = computeQdot
      else
          prolongQdot = .FALSE.
      end if

      ! if (prolongQdot) then
      !     print *, "side: ", side
      !     ! print *, "projectionType 1: ",  self % projectionType(1)
      !     ! print *, "projectionType 2: ",  self % projectionType(2)
      !     print *, "projectionType side: ",  self % projectionType(side)
      ! end if

      select case (side)
      case(1)
         associate(Qf => self % storage(1) % Q)
         select case ( self % projectionType(1) )
         case (0)
            Qf = Qe
            if (prolongQdot) self % storage(1) % Qdot = QdotE

         case (1)
            Qf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Qe(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Qf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Qe(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Qf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * Qe(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2)
         associate( Qf => self % storage(2) % Q )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            Qe_rot(:,i,j) = Qe(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            Qf = Qe_rot
         case (1)
            Qf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Qe_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Qf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Qe_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Qf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * Qe_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      end select

   end subroutine Face_AdaptSolutionToFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Face_AdaptGradientsToFace(self, nEqn, Nelx, Nely, Uxe, Uye, Uze, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn 
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: Uxe(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP), intent(in)     :: Uye(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP), intent(in)     :: Uze(nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: Uxe_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      real(kind=RP) :: Uye_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      real(kind=RP) :: Uze_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      select case (side)
      case(1)
         associate(Uxf => self % storage(1) % U_x, &
                   Uyf => self % storage(1) % U_y, &
                   Uzf => self % storage(1) % U_z   )
         select case ( self % projectionType(1) )
         case (0)
            Uxf = Uxe
            Uyf = Uye
            Uzf = Uze
         case (1)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uxe(:,l,j)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uye(:,l,j)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uze(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uxe(:,i,l)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uye(:,i,l)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uze(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * Uxe(:,m,l)
                  Uyf(:,i,j) = Uyf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * Uye(:,m,l)
                  Uzf(:,i,j) = Uzf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * Uze(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2)
         associate(Uxf => self % storage(2) % U_x, &
                   Uyf => self % storage(2) % U_y, &
                   Uzf => self % storage(2) % U_z   )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            Uxe_rot(:,i,j) = Uxe(:,ii,jj) 
            Uye_rot(:,i,j) = Uye(:,ii,jj) 
            Uze_rot(:,i,j) = Uze(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            Uxf = Uxe_rot
            Uyf = Uye_rot
            Uzf = Uze_rot
         case (1)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uxe_rot(:,l,j)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uye_rot(:,l,j)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uze_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uxe_rot(:,i,l)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uye_rot(:,i,l)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uze_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * Uxe_rot(:,m,l)
                  Uyf(:,i,j) = Uyf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * Uye_rot(:,m,l)
                  Uzf(:,i,j) = Uzf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * Uze_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      end select

   end subroutine Face_AdaptGradientsToFace

   subroutine Face_AdaptAviscFluxToFace(self, nEqn, Nelx, Nely, AVn_e, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: AVn_e(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: AVn_e_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      select case (side)
      case(1)
         associate(AVf => self % storage(1) % AviscFlux)
         select case ( self % projectionType(1) )
         case (0)
            AVf = AVn_e
         case (1)
            AVf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * AVn_e(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            AVf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * AVn_e(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            AVf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * AVn_e(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2)
         associate( AVf => self % storage(2) % AviscFlux )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            AVn_e_rot(:,i,j) = AVn_e(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            AVf = AVn_e_rot
         case (1)
            AVf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * AVn_e_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            AVf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * AVn_e_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            AVf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * AVn_e_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select

         AVf = -AVf
         end associate
      end select

   end subroutine Face_AdaptAviscFluxToFace

   subroutine Face_ProjectFluxToElements(self, nEqn, flux, whichElements)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(in)  :: flux(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
      integer,       intent(in)  :: whichElements(2)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer           :: i, j, ii, jj, l, m, side
      real(kind=RP)     :: fStarAux(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      do side = 1, 2
         select case ( whichElements(side) )
         case (1)    ! Prolong to left element
            associate(fStar => self % storage(1) % Fstar)
            select case ( self % projectionType(1) )
            case (0)
               fStar(1:nEqn,:,:) = flux
            case (1)
               fStar(1:nEqn,:,:) = 0.0
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * flux(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStar(1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * flux(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStar(1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) +   Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                             * flux(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate

         case (2)    ! Prolong to right element
!      
!           *********
!           1st stage: Projection
!           *********
!      
            select case ( self % projectionType(2) )
            case (0)
               fStarAux(1:nEqn,:,:) = flux
            case (1)
               fStarAux(1:nEqn,:,:) = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * flux(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStarAux(1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * flux(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStarAux(1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                             * flux(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
            associate(fStar => self % storage(2) % Fstar)
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               fStar(1:nEqn,ii,jj) = fStarAux(1:nEqn,i,j) 
            end do                        ; end do
!
!           *********
!           3rd stage: Inversion
!           *********
!
            fStar = -fStar

            end associate
         end select
      end do

   end subroutine Face_ProjectFluxToElements
#if defined(NAVIERSTOKES) 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  Face_ProjectBCJacobianToElements:
!  Routine to project the Jacobian of the boundary condition to the element
!  --
   subroutine Face_ProjectBCJacobianToElements(self,nEqn,BCjacF)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face), target        :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(in)  :: BCjacF(nEqn,nEqn,0:self % Nf(1), 0:self % Nf(2))
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                :: i, j, l, m
      
!     Always project to left element (this is only for boundaries)
!     ------------------------------------------------------------
      associate(BCJac => self % storage(1) % BCJac)
      
      select case ( self % projectionType(1) )
         case (0)
            BCJac(1:nEqn,1:nEqn,:,:) = BCJacF
         case (1)
            BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
            do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
               BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                            + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * BCJacF(:,:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
            do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
               BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                            + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * BCJacF(:,:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
            do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
               do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                  BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                           +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                            * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                            * BCJacF(:,:,m,l)
               end do                 ; end do
            end do                  ; end do
      end select
      
      end associate
   end subroutine Face_ProjectBCJacobianToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!            
   subroutine Face_ProjectFluxJacobianToElements(self, nEqn, whichElement,whichderiv)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face), target        :: self
      integer,       intent(in)  :: nEqn
      integer,       intent(in)  :: whichElement
      integer,       intent(in)  :: whichderiv           !<  One can either transfer the derivative with respect to qL (LEFT) or to qR (RIGHT)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                :: i, j, ii, jj, l, m, side
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:)
      real(kind=RP)          :: fStarAux(nEqn,nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      
      fluxDeriv(1:,1:,0:,0:) => self % storage(whichderiv) % dFStar_dqF
      
      select case ( whichElement )
         case (LEFT)    ! Prolong to left element
            associate(dFStar_dq => self % storage(1) % dFStar_dqEl)
            select case ( self % projectionType(1) )
            case (0)
               dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = fluxDeriv
            case (1)
               dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                               + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                               + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                              +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                               * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                               * fluxDeriv(:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate

         case (RIGHT)    ! Prolong to right element
!      
!           *********
!           1st stage: Projection
!           *********
!      
            select case ( self % projectionType(2) )
            case (0)
               fStarAux(1:nEqn,1:nEqn,:,:) = fluxDeriv
            case (1)
               fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                             * fluxDeriv(:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
            associate(dFStar_dq => self % storage(2) % dFStar_dqEl)
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               dFStar_dq(1:nEqn,1:nEqn,ii,jj,whichderiv) = fStarAux(1:nEqn,1:nEqn,i,j) 
            end do                        ; end do
!
!           *********
!           3rd stage: Inversion
!           *********
!
            dFStar_dq = -dFStar_dq
            end associate
      end select
      

   end subroutine Face_ProjectFluxJacobianToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#endif
   subroutine Face_ProjectGradientFluxToElements(self, nEqn, Hflux, whichElements, factor)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(in)  :: Hflux(nEqn, NDIM, 0:self % Nf(1), 0:self % Nf(2))
      integer,       intent(in)  :: whichElements(2)
      integer,       intent(in)  :: factor               ! A factor that relates LEFT and RIGHT fluxes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer           :: i, j, ii, jj, l, m, side
      real(kind=RP)     :: hStarAux(nEqn, NDIM, 0:self % NfRight(1), 0:self % NfRight(2))

      do side = 1, 2
         select case ( whichElements(side) )
         case (1)    ! Prolong from left element
            associate(unStar => self % storage(1) % unStar)
            select case ( self % projectionType(1) )
            case (0)
               unStar(:,:,:,:) = Hflux
            case (1)
               unStar(:,:,:,:) = 0.0
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  unStar(:,:,i,j) = unStar(:,:,i,j) + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * Hflux(:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               unStar(:,:,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  unStar(:,:,i,j) = unStar(:,:,i,j) + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * Hflux(:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               unStar(:,:,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     unStar(:,:,i,j) = unStar(:,:,i,j) +   Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                             * Hflux(:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate

         case (2)    ! Prolong from right element
!      
!           *********
!           1st stage: Projection
!           *********
!      
            select case ( self % projectionType(2) )
            case (0)
               HstarAux = Hflux
            case (1)
               HstarAux = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  HstarAux(:,:,i,j) = HstarAux(:,:,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * Hflux(:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               HstarAux = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  HstarAux(:,:,i,j) = HstarAux(:,:,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * Hflux(:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               HstarAux = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     HstarAux(:,:,i,j) = HstarAux(:,:,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                             * Hflux(:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
            associate(unStar => self % storage(2) % unStar)
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               unStar(:,:,ii,jj) = HstarAux(:,:,i,j) 
            end do                        ; end do
!
!           *********
!           3rd stage: Multiplication by a factor (inversion usually)
!           *********
!
            unStar = factor * unStar

            end associate
         end select
      end do

   end subroutine Face_ProjectGradientFluxToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
   subroutine Face_ProjectGradJacobianToElements(self, whichElement, whichderiv)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      !---------------------------------------------------------
      class(Face), target        :: self
      integer,       intent(in)  :: whichElement
      integer,       intent(in)  :: whichderiv
      !---------------------------------------------------------
      integer                :: i, j, ii, jj, l, m, side
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:,:)
      real(kind=RP)          :: fStarAux(NCONS,NCONS,1:NDIM, 0:self % NfRight(1), 0:self % NfRight(2))
      
      !---------------------------------------------------------
      
      fluxDeriv(1:,1:,1:,0:,0:) => self % storage(whichderiv) % dFv_dGradQF
      
      select case ( whichElement )
         case (LEFT)    ! Prolong to left element
            associate(dFv_dGradQEl => self % storage(1) % dFv_dGradQEl)
            
            select case ( self % projectionType(1) )
            case (0)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = fluxDeriv
            case (1)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                          + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                          + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                                   +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                                    * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                                    * fluxDeriv(:,:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            
            end associate

         case (RIGHT)    ! Prolong to right element
!      
!           *********
!           1st stage: Projection
!           *********
!      
            select case ( self % projectionType(2) )
            case (0)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = fluxDeriv
            case (1)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                             * fluxDeriv(:,:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
            associate(dFv_dGradQEl => self % storage(2) % dFv_dGradQEl)
            
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,ii,jj,whichderiv) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) 
            end do                        ; end do
!
!           *********
!           3rd stage: Inversion
!           *********
!
            dFv_dGradQEl = -dFv_dGradQEl
            
            end associate
      end select
      
      
   end subroutine Face_ProjectGradJacobianToElements
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE InterpolateToBoundary( u, v, Nx, Ny, Nz, which_dim , bValue , NEQ)
      use SMConstants
      use PhysicsStorage
      IMPLICIT NONE
      integer                      , INTENT(IN)    :: Nx, Ny, Nz
      real(kind=RP)                , intent(in)    :: u(1:NEQ,0:Nx,0:Ny,0:Nz) , v(0:)
      integer                      , intent(in)    :: which_dim
      REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
      integer                      , intent(in)    :: NEQ
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                                    :: i , j , k , eq

      select case (which_dim)
      case (IX)

         do j = 0 , Nz ; do i = 0 , Ny ; do k = 0 , Nx
            bValue(:,i,j) = bValue(:,i,j) + u(:,k,i,j) * v(k)
         end do        ; end do        ; end do

      case (IY)

         do j = 0 , Nz ; do k = 0 , Ny ; do i = 0 , Nx
            bValue(:,i,j) = bValue(:,i,j) + u(:,i,k,j) * v(k)
         end do        ; end do        ; end do

      case (IZ)

         do k = 0 , Nz ; do j = 0 , Ny ; do i = 0 , Nx
            bValue(:,i,j) = bValue(:,i,j) + u(:,i,j,k) * v(k)
         end do        ; end do        ; end do

      end select
      
      end SUBROUTINE InterpolateToBoundary
!
!////////////////////////////////////////////////////////////////////////
!
      elemental subroutine Face_Assign(to,from)
         implicit none
         class(Face), intent(inout) :: to
         class(Face), intent(in)    :: from
         
         
         to % flat = from % flat
         to % ID = from % ID
         to % FaceType = from % FaceType
         to % zone = from % zone
         to % rotation = from % rotation
         to % NelLeft = from % NelLeft
         to % NelRight = from % NelRight
         to % NfLeft = from % NfLeft
         to % NfRight = from % NfRight
         to % Nf = from % Nf
         to % nodeIDs = from % nodeIDs
         to % elementIDs = from % elementIDs
         to % elementSide = from % elementSide
         to % projectionType = from % projectionType
         to % boundaryName = from % boundaryName
         to % geom = from % geom
         to % storage = from % storage
      end subroutine Face_Assign

      subroutine Face_HO_IBM_correction( self, nEqn, Nelx, Nely, sideOut, Q, eID )
         use MPI_Process_Info
         use FluidData
         implicit none 

         class(face),   intent(inout) :: self 
         real(kind=RP), intent(inout) :: Q(nEqn,0:Nelx,0:Nely)
         integer,       intent(in)    :: nEqn, eID, Nelx, Nely, sideOut

         real(kind=RP) :: rho, p, u2, v2, w2 
         integer       :: i, j, side, domain 

         side = self% HOSIDE 

         if( side .eq. sideOut ) then 
            domain = self% domain 
            associate( f => IBM_HO_faces(domain)% faces(self% HO_ID) )
            do i = 0, Nelx; do j = 0, Nely
#if defined(NAVIERSTOKES)
               Q(IRHO,i,j)  = f% stencil(i,j)% rho 
               Q(IRHOU,i,j) = f% stencil(i,j)% u 
               Q(IRHOV,i,j) = f% stencil(i,j)% v 
               Q(IRHOW,i,j) = f% stencil(i,j)% w 
               u2 = f% stencil(i,j)% u**2
               v2 = f% stencil(i,j)% v**2
               w2 = f% stencil(i,j)% w**2
               Q(IRHOE,i,j) = f% stencil(i,j)% p/(thermodynamics% gamma - 1._RP) + 0.5_RP * 1._RP/f% stencil(i,j)% rho * &
                              ( u2 + v2 + w2 ) 
#endif
            end do; end do    
            end associate      
         end if 
         
         end subroutine Face_HO_IBM_correction

         subroutine Face_HO_IBM_Gradcorrection( self, sideIn, sideOut )

            implicit none 

            class(face), intent(inout) :: self 
            integer,     intent(in)    :: sideIn, sideOut

            integer :: i, j 

            do i = 0, self% Nf(1); do j = 0, self% Nf(2) 
               self% storage(sideOut)% U_x(:,i,j) = self% storage(sideIn)% U_x(:,i,j) 
               self% storage(sideOut)% U_y(:,i,j) = self% storage(sideIn)% U_y(:,i,j) 
               self% storage(sideOut)% U_z(:,i,j) = self% storage(sideIn)% U_z(:,i,j) 
            end do; end do    

         end subroutine Face_HO_IBM_Gradcorrection
! STENCIL PROCEDURES 

      subroutine face_StencilConstruct( this )

         implicit none 

         class(face), intent(inout) :: this 

         integer :: i, j, N 

         allocate( this% stencil(0:this% Nf(1),0:this% Nf(2)) )

         N = max(this% Nf(1),this% Nf(2))

         do i = 0, this% Nf(1); do j = 0, this% Nf(2) 
            this% stencil(i,j)% N = N 
            allocate( this% stencil(i,j)% Q(NCONS,0:N) )
            if( .not. allocated(this% stencil(i,j)% x_s)  ) allocate(this% stencil(i,j)% x_s(NDIM,0:N))
            if( .not. allocated(this% stencil(i,j)% xi_s) ) allocate(this% stencil(i,j)% xi_s(NDIM,0:N))
            if( .not. allocated(this% stencil(i,j)% eIDs) ) allocate(this% stencil(i,j)% eIDs(0:N))
         end do; end do 
         
      end subroutine face_StencilConstruct

      subroutine stencil_build( this, x0, normal, L, N )

         implicit none 

         class(stencil_t), intent(inout) :: this 
         real(kind=RP),    intent(in)    :: x0(NDIM), normal(NDIM), L 
         integer,          intent(in)    :: N 

         real(kind=RP) :: xi 
         integer       :: i 

         allocate( this% x_s(NDIM,0:N),  &
                   this% xi_s(NDIM,0:N), &
                   this% eIDs(0:N)       )

         do i = 0, N 
            xi = NodalStorage(N)% x(i)
            this% x_s(:,i) = x0 + 0.5_RP * (1.0_RP + xi) * L * normal
         end do 

      end subroutine stencil_build

      subroutine stencil_ComputeState( this )
         use VariableConversion, only:Pressure
         implicit none 

         class(stencil_t), intent(inout) :: this 

         real(kind=RP)                 :: u_s, v_s, w_s, rho_s, p_s
         real(kind=RP)                 :: uB, vB, wB, drhoB, dpB, p(0:this% N)
         real(kind=RP)                 :: lj(0:this% N), dlj(0:this% N)
         integer                       :: i, N 
         type(NodalStorage_t), pointer :: spA
#if defined(NAVIERSTOKES)     
         N = this% N 

         uB    = 0.0_RP
         vB    = 0.0_RP
         wB    = 0.0_RP
         drhoB = 0.0_RP
         dpB   = 0.0_RP

         spA => NodalStorage(N)
      
         lj  = spA% lj(this% xiB)
         dlj = spA% dlj(this% xiB)

         this% Q(IRHO ,0) = drhoB
         this% Q(IRHOU,0) = uB 
         this% Q(IRHOV,0) = vB 
         this% Q(IRHOW,0) = wB 
         p(0)             = dpB 
         
         do i = 1, N 
            p(i) = Pressure(this% Q(:,i))
         end do

         do i = 1, N
            this% Q(IRHO ,0) = this% Q(IRHO ,0) - this% Q(IRHO ,i) * dlj(i) 
            this% Q(IRHOU,0) = this% Q(IRHOU,0) - this% Q(IRHOU,i) * lj(i)
            this% Q(IRHOV,0) = this% Q(IRHOV,0) - this% Q(IRHOV,i) * lj(i)
            this% Q(IRHOW,0) = this% Q(IRHOW,0) - this% Q(IRHOW,i) * lj(i)
            p(0)             = p(0) - p(i) * dlj(i)
         end do
         this% Q(IRHO ,0) = this% Q(IRHO, 0)/dlj(0)
         this% Q(IRHOU,0) = this% Q(IRHOU,0)/lj(0)
         this% Q(IRHOV,0) = this% Q(IRHOV,0)/lj(0)
         this% Q(IRHOW,0) = this% Q(IRHOW,0)/lj(0)
         p(0)             = p(0)/dlj(0)

         lj = spA% lj(this% xiI)

         u_s   = 0.0_RP 
         v_s   = 0.0_RP 
         w_s   = 0.0_RP 
         rho_s = 0.0_RP 
         p_s   = 0.0_RP

         do i = 0, N 
            u_s   = u_s   + this% Q(IRHOU,i) * lj(i)
            v_s   = v_s   + this% Q(IRHOV,i) * lj(i)
            w_s   = w_s   + this% Q(IRHOW,i) * lj(i)
            rho_s = rho_s + this% Q(IRHO ,i) * lj(i) 
            p_s   = p_s   + p(i)             * lj(i) 
         end do   

         this% u   = u_s
         this% v   = v_s 
         this% w   = w_s 
         this% rho = rho_s 
         this% p   = p_s 
#endif   
      end subroutine stencil_ComputeState

      subroutine stencil_destroy( this )

         implicit none 

         class(stencil_t), intent(inout) :: this

         deallocate( this% x_s,  &
                     this% xi_s, &
                     this% eIDs  )

      end subroutine stencil_destroy



end Module FaceClass