!!
!////////////////////////////////////////////////////////////////////////
!
!      FaceClass.f
!      Created: 2008-06-05 14:12:52 -0400 
!      By: David Kopriva  
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
      use StorageClass
      use PhysicsStorage
      use NodalStorageClass
      use InterpolationMatrices
      IMPLICIT NONE 

      private
      public   Face
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
         class(NodalStorage_t), pointer    :: spAxi
         class(NodalStorage_t), pointer    :: spAeta
         contains
            procedure   :: Construct => ConstructFace
            procedure   :: Destruct  => DestructFace
            procedure   :: Print     => PrintFace
            procedure   :: LinkWithElements      => Face_LinkWithElements
            procedure   :: AdaptSolutionToFace   => Face_AdaptSolutionToFace
            procedure   :: AdaptGradientsToFace   => Face_AdaptGradientsToFace
            procedure   :: ProjectFluxToElements => Face_ProjectFluxToElements
            procedure   :: ProjectGradientFluxToElements => Face_ProjectGradientFluxToElements
#if defined(NAVIERSTOKES)
            procedure   :: ProjectFluxJacobianToElements => Face_ProjectFluxJacobianToElements
            procedure   :: ProjectGradJacobianToElements => Face_ProjectGradJacobianToElements
            procedure   :: AdaptDensityGradientToFace => Face_AdaptDensityGradientToFace
#endif
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
               
         self % spAxi => NULL()
         self % spAeta => NULL()

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
      self % NfLeft = self % NelLeft     ! Left elements is always oriented.
      
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
!     ----------------------
!     Construct face storage
!     ----------------------
!
      if ( .not. NodalStorage(self % Nf(1)) % Constructed ) then
         call NodalStorage(self % Nf(1)) % Construct(nodeType, self % Nf(1))
      end if

      if ( .not. NodalStorage(self % Nf(2)) % Constructed ) then
         call NodalStorage(self % Nf(2)) % Construct(nodeType, self % Nf(2))
      end if

      self % spAxi => NodalStorage(self % Nf(1))
      self % spAeta => NodalStorage(self % Nf(2))

      call self % storage(1) % Construct(NDIM, self % Nf, self % NelLeft, computeGradients)
      call self % storage(2) % Construct(NDIM, self % Nf, self % NelRight, computeGradients)
!
!     -----------------------------------------------------------------------
!     Construction of the projection matrices (simple Lagrange interpolation)
!     -----------------------------------------------------------------------
!
      call ConstructInterpolationMatrices(self % NfLeft(1), self % Nf(1))
      call ConstructInterpolationMatrices(self % NfLeft(2), self % Nf(2))

      call ConstructInterpolationMatrices(self % NfRight(1), self % Nf(1))
      call ConstructInterpolationMatrices(self % NfRight(2), self % Nf(2))
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
   subroutine Face_AdaptSolutionToFace(self, nEqn, Nelx, Nely, Qe, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: Qe(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: Qe_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      select case (side)
      case(1)
         associate(Qf => self % storage(1) % Q)
         select case ( self % projectionType(1) )
         case (0)
            Qf = Qe
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
            call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
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
!  -------------------------------
!  Only needed by implicit methods
!  -------------------------------
#if defined(NAVIERSTOKES)
   subroutine Face_AdaptDensityGradientToFace(self, Nelx, Nely, Qe, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: Qe(1:NDIM, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: Qe_rot(1:NDIM, 0:self % NfRight(1), 0:self % NfRight(2))

      select case (side)
      case(1)
         associate(Qf => self % storage(1) % gradRho)
         select case ( self % projectionType(1) )
         case (0)
            Qf = Qe
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
         associate( Qf => self % storage(2) % gradRho )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
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

   end subroutine Face_AdaptDensityGradientToFace
#endif
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
            call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
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
         case (1)    ! Prolong from left element
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

         case (2)    ! Prolong from right element
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
               call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES) 
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
      ! real(kind=RP), intent(in)  :: flux(1:nEqn,1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
      integer                :: i, j, ii, jj, l, m, side
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:)
      real(kind=RP)     :: fStarAux(nEqn,nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      
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
               call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
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
               call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
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
   subroutine Face_ProjectGradJacobianToElements(self, whichElement)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      !---------------------------------------------------------
      class(Face), target        :: self
      integer,       intent(in)  :: whichElement
      !---------------------------------------------------------
      integer                :: i, j, ii, jj, l, m, side
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:,:,:)
      real(kind=RP)          :: fStarAux(NCONS,NCONS,1:NDIM,2, 0:self % NfRight(1), 0:self % NfRight(2))
      integer                :: whichderiv
      !---------------------------------------------------------
      
      whichderiv = whichElement  ! Hardcoded: df/d∇q⁺ goes to left element and df/d∇q⁻ goes to right element
      
      fluxDeriv(1:,1:,1:,1:,0:,0:) => self % storage(whichderiv) % dFv_dGradQF
      
      select case ( whichElement )
         case (LEFT)    ! Prolong to left element
            associate(dFv_dGradQEl => self % storage(1) % dFv_dGradQEl)
            
            select case ( self % projectionType(1) )
            case (0)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = fluxDeriv
            case (1)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = 0._RP
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) &
                                                               + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) &
                                                               + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) &
                                                                                   +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                                    * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                                    * fluxDeriv(:,:,:,:,m,l)
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
               fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = fluxDeriv
            case (1)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                             * fluxDeriv(:,:,:,:,m,l)
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
               call iijjIndexes(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,1:2,ii,jj) = fStarAux(1:NCONS,1:NCONS,1:NDIM,1:2,i,j) 
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
end Module FaceClass
