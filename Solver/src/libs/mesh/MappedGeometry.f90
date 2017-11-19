!
!////////////////////////////////////////////////////////////////////////
!
!      MappedGeometry.f95
!      Created: 2008-06-19 15:58:02 -0400 
!      By: David Kopriva  
!
!      Modification history:
!        2008-06-19: Created by David Kopriva
!        XXXX-XX-XX: Gonzalo Rubio implemented cross-product metrics
!        2017-05-05: Andrés Rueda implemented polynomial anisotropy
!      Contains:
!         ALGORITHM 101: MappedGeometryClass
!         ALGORITHM 102: ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
Module MappedGeometryClass 
   USE SMConstants
   USE TransfiniteMapClass
   USE NodalStorageClass
   IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      integer, parameter :: EFRONT = 1, EBACK = 2, EBOTTOM = 3
      integer, parameter :: ERIGHT = 4, ETOP = 5, ELEFT = 6
!
!     -----
!     Class
!     -----
!
      TYPE MappedGeometry
            INTEGER                                         :: Nx, Ny, Nz                    ! Polynomial order
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: jGradXi, jGradEta, jGradZeta  ! 
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: x                             ! Position of points in absolute coordinates
            REAL(KIND=RP), DIMENSION(:,:,:)   , ALLOCATABLE :: jacobian 
            
            CONTAINS
            
            PROCEDURE :: construct => ConstructMappedGeometry
            PROCEDURE :: destruct  => DestructMappedGeometry
      END TYPE MappedGeometry
      
      type MappedGeometryFace
         real(kind=RP), dimension(:,:,:),   allocatable :: x
         real(kind=RP), dimension(:,:)  , allocatable :: scal   ! |ja^i|: Normalization term of the normal vectors on a face
         real(kind=RP), dimension(:,:,:), allocatable :: normal ! normal vector on a face
         contains
            procedure :: construct => ConstructMappedGeometryFace
            procedure :: destruct  => DestructMappedGeometryFace
      end type MappedGeometryFace
      
      LOGICAL       :: useCrossProductMetrics = .false. ! A switch for debugging purposes. Cross product metrics are fine (and more precise) for 2D geometries... But not for 3D.
                                                        ! Before changing, read: Kopriva, David A. "Metric identities and the discontinuous spectral element method on curvilinear meshes." Journal of Scientific Computing 26.3 (2006): 301-327.
!
!  ========
   CONTAINS 
!  ========
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ConstructMappedGeometry( self, spAxi, spAeta, spAzeta, mapper )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
      CLASS(MappedGeometry)  , intent(inout) :: self
      TYPE(TransfiniteHexMap), intent(in)    :: mapper
      TYPE(NodalStorage)     , intent(in)    :: spAxi
      TYPE(NodalStorage)     , intent(in)    :: spAeta
      TYPE(NodalStorage)     , intent(in)    :: spAzeta
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: Nx, Ny, Nz, Nmax
      INTEGER       :: i, j, k
      REAL(KIND=RP) :: nrm
      REAL(KIND=RP) :: grad_x(3,3), jGrad(3)
!
!     -----------
!     Allocations
!     -----------
!
      Nx        = spAxi   % N
      Ny        = spAeta  % N
      Nz        = spAzeta % N
      Nmax      = MAX(Nx,Ny,Nz)
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
      
      ALLOCATE( self % JGradXi  (3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % JGradEta (3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % JGradZeta(3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % jacobian   (0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % x        (3,0:Nx,0:Ny,0:Nz)    )
!
!     --------------------------
!     Compute interior locations
!     --------------------------
!
!     TODO: this could be inconsistent if curveOrder > N. Otherwise I think is ok.
      DO k = 0, Nz
         DO j= 0, Ny       
            DO i = 0,Nx 
               self % x(:,i,j,k) = mapper %  transfiniteMapAt([spAxi % x(i), spAeta % x(j), spAzeta % x(k)])
            END DO
         END DO
      END DO
!
!     ------------
!     Metric terms
!     ------------
!
!
!     ------------------------------------------------------------
!     If the faces are straight, the CrossProductForm is OK
!     If there are curved faces, the Conservative form is required
!     see Kopriva 2006
!     ------------------------------------------------------------
!
      IF ( useCrossProductMetrics) THEN 
      
         CALL computeMetricTermsCrossProductForm(self, spAxi, spAeta, spAzeta, mapper)
         
      ELSE
         
         CALL computeMetricTermsConservativeForm(self, spAxi, spAeta, spAzeta, mapper)
      
      ENDIF
      
   END SUBROUTINE ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructMappedGeometry(self)
         IMPLICIT NONE 
         CLASS(MappedGeometry) :: self
         DEALLOCATE( self % jGradXi, self % jGradEta, self % jGradZeta, self % jacobian )
         DEALLOCATE( self % x)
      END SUBROUTINE DestructMappedGeometry
!
!//////////////////////////////////////////////////////////////////////// 
!
!  -----------------------------------------------------------------------------------
!  Computation of the metric terms on a face: TODO only the Left element (rotation 0)
!  -----------------------------------------------------------------------------------
   subroutine ConstructMappedGeometryFace(self, Nf, spA, geom, hexMap, side)
      use PhysicsStorage
      implicit none
      class(MappedGeometryFace), intent(inout)  :: self
      integer,                   intent(in)     :: Nf(2)
      type(NodalStorage),        intent(in)     :: spA(2)
      type(MappedGeometry),      intent(in)     :: geom
      type(TransfiniteHexMap),   intent(in)     :: hexMap
      integer,                   intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i, j, l
      real(kind=RP)  :: x_xi(NDIM,0:Nf(1),0:Nf(2))
      real(kind=RP)  :: x_eta(NDIM,0:Nf(1),0:Nf(2))
      real(kind=RP)  :: crossTerm1(0:Nf(1),0:Nf(2))
      real(kind=RP)  :: crossTerm2(0:Nf(1),0:Nf(2))
      real(kind=RP)  :: crossTerm1_eta(0:Nf(1),0:Nf(2))
      real(kind=RP)  :: crossTerm2_xi(0:Nf(1),0:Nf(2))
      
      
      allocate( self % x(NDIM, 0:Nf(1), 0:Nf(2)))
      allocate( self % scal(0:Nf(1), 0:Nf(2)))
      allocate( self % normal(NDIM, 0:Nf(1), 0:Nf(2)))
!
!     First step: get the surface coordinates directly from the mapping
!     -----------------------------------------------------------------
      select case(side)
         case(ELEFT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([-1.0_RP      , spA(1) % x(i), spA(2) % x(j) ])
            end do ; end do
         
         case(ERIGHT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([ 1.0_RP      , spA(1) % x(i), spA(2) % x(j) ])
            end do ; end do
         
         case(EBOTTOM)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), spA(2) % x(j),    -1.0_RP    ])
            end do ; end do
            
         case(ETOP)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), spA(2) % x(j),     1.0_RP    ])
            end do ; end do
            
         case(EFRONT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), -1.0_RP      , spA(2) % x(j) ])
            end do ; end do
            
         case(EBACK)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i),  1.0_RP      , spA(2) % x(j) ])
            end do ; end do
            
      end select
!
!     Get the mappings interpolant derivatives
!     ----------------------------------------
      x_xi = 0.0_RP
      do j = 0, Nf(2) ; do i = 0, Nf(1) ; do l = 0, Nf(1)  
         x_xi(:,i,j) = x_xi(:,i,j) + spA(1) % D(i,l) * self % x(:,l,j)
      end do              ; end do              ; end do

      x_eta = 0.0_RP
      do j = 0, Nf(2) ; do i = 0, Nf(1) ; do l = 0, Nf(2)  
         x_eta(:,i,j) = x_eta(:,i,j) + spA(2) % D(j,l) * self % x(:,i,l)
      end do              ; end do              ; end do
!
!     Compute the metric terms
!     ------------------------
      if ( useCrossProductMetrics ) then

         do j = 0, Nf(2) ; do i = 0, Nf(1)
            call vCross(x_xi(:,i,j),x_eta(:,i,j), self % normal(:,i,j))
            self % scal(i,j) = norm2(self % normal(:,i,j))
            self % normal(:,i,j) = self % normal(:,i,j) / self % scal(i,j)
         end do              ; end do

      else
!
!        ************************
!        Compute the x-coordinate
!        ************************
!
         crossTerm1 = x_xi(2,:,:) * self % x(3,:,:)    ! y_xi z
         crossTerm2 = x_eta(2,:,:) * self % x(3,:,:)   ! y_eta z
!
!        Compute the derivative of the cross terms
!        -----------------------------------------
         crossTerm1_eta = 0.0_RP
         crossTerm2_xi  = 0.0_RP
         do j = 0, Nf(2)   ; do i = 0, Nf(2)
            do l = 0, Nf(2)
               crossTerm1_eta(i,j) = crossTerm1_eta(i,j) + crossTerm1(i,l) * spA(2) % D(j,l)
            end do

            do l = 0, Nf(1)
               crossTerm2_xi(i,j) = crossTerm2_xi(i,j) + crossTerm2(l,j) * spA(1) % D(i,l)
            end do
         end do            ; end do       
!
!        Assign it to the normal
!        -----------------------
         self % normal(1,:,:) = crossTerm1_eta - crossTerm2_xi
!
!        ************************
!        Compute the y-coordinate
!        ************************
!
         crossTerm1 = x_xi(3,:,:) * self % x(1,:,:)    ! z_xi x
         crossTerm2 = x_eta(3,:,:) * self % x(1,:,:)   ! z_eta x
!
!        Compute the derivative of the cross terms
!        -----------------------------------------
         crossTerm1_eta = 0.0_RP
         crossTerm2_xi  = 0.0_RP
         do j = 0, Nf(2)   ; do i = 0, Nf(2)
            do l = 0, Nf(2)
               crossTerm1_eta(i,j) = crossTerm1_eta(i,j) + crossTerm1(i,l) * spA(2) % D(j,l)
            end do

            do l = 0, Nf(1)
               crossTerm2_xi(i,j) = crossTerm2_xi(i,j) + crossTerm2(l,j) * spA(1) % D(i,l)
            end do
         end do            ; end do       
!
!        Assign it to the normal
!        -----------------------
         self % normal(2,:,:) = crossTerm1_eta - crossTerm2_xi
!
!        ************************
!        Compute the z-coordinate
!        ************************
!
         crossTerm1 = x_xi(1,:,:) * self % x(2,:,:)    ! x_xi y
         crossTerm2 = x_eta(1,:,:) * self % x(2,:,:)   ! x_eta y
!
!        Compute the derivative of the cross terms
!        -----------------------------------------
         crossTerm1_eta = 0.0_RP
         crossTerm2_xi  = 0.0_RP
         do j = 0, Nf(2)   ; do i = 0, Nf(2)
            do l = 0, Nf(2)
               crossTerm1_eta(i,j) = crossTerm1_eta(i,j) + crossTerm1(i,l) * spA(2) % D(j,l)
            end do

            do l = 0, Nf(1)
               crossTerm2_xi(i,j) = crossTerm2_xi(i,j) + crossTerm2(l,j) * spA(1) % D(i,l)
            end do
         end do            ; end do       
!
!        Assign it to the normal
!        -----------------------
         self % normal(3,:,:) = crossTerm1_eta - crossTerm2_xi
!
!        ***********************
!        Compute scal and normal
!        ***********************
!
         do j = 0, Nf(2)   ; do i = 0, Nf(1)
            self % scal(i,j) = norm2(self % normal(:,i,j))
            self % normal(:,i,j) = self % normal(:,i,j) / self % scal(i,j)
         end do            ; end do
       
      end if

      if ( (side .eq. ELEFT) .or. (side .eq. EBACK) .or. (side .eq. EBOTTOM)) self % normal = -self % normal

   end subroutine ConstructMappedGeometryFace
!
!//////////////////////////////////////////////////////////////////////// 
!
      subroutine DestructMappedGeometryFace(self)
         implicit none
         !-------------------------------------------------------------------
         class(MappedGeometryFace), intent(inout) :: self
         !-------------------------------------------------------------------
         
         deallocate (self % x    )
         deallocate (self % scal  )
         deallocate (self % normal)
         
      end subroutine DestructMappedGeometryFace
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine computeMetricTermsConservativeForm(self, spAxi, spAeta, spAzeta, mapper)
         use PhysicsStorage
         implicit none
         type(MappedGeometry),    intent(inout) :: self
         type(NodalStorage),      intent(in)    :: spAxi, spAeta, spAzeta
         type(TransfiniteHexMap), intent(in)    :: mapper
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l
         real(kind=RP)  :: grad_x(NDIM,NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: auxgrad(NDIM,NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: coordsProduct(NDIM,0:self % Nx,0:self % Ny,0:self % Nz)
         real(kind=RP)  :: Jai(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
!
!        Compute the mapping gradient
!        ----------------------------
         grad_x = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               grad_x(:,1,i,j,k) = grad_x(:,1,i,j,k) + self % x(:,l,j,k) * spAxi % D(i,l)
            end do
      
            do l = 0, self % Ny
               grad_x(:,2,i,j,k) = grad_x(:,2,i,j,k) + self % x(:,i,l,k) * spAeta % D(j,l)
            end do

            do l = 0, self % Nz
               grad_x(:,3,i,j,k) = grad_x(:,3,i,j,k) + self % x(:,i,j,l) * spAzeta % D(k,l)
            end do
         end do         ; end do          ; end do
!
!        *****************************************
!        Compute the x-coordinates of the mappings
!        *****************************************
!
!        Compute coordinates combination
!        -------------------------------
         do k = 0, self % Nz    ; do j = 0, self % Ny  ; do i = 0, self % Nx
            coordsProduct(:,i,j,k) = self % x(3,i,j,k) * grad_x(2,:,i,j,k)
         end do            ; end do          ; end do
!
!        Compute its gradient
!        --------------------
         auxgrad = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               auxgrad(:,1,i,j,k) = auxgrad(:,1,i,j,k) + coordsProduct(:,l,j,k) * spAxi % D(i,l)
            end do
      
            do l = 0, self % Ny
               auxgrad(:,2,i,j,k) = auxgrad(:,2,i,j,k) + coordsProduct(:,i,l,k) * spAeta % D(j,l)
            end do

            do l = 0, self % Nz
               auxgrad(:,3,i,j,k) = auxgrad(:,3,i,j,k) + coordsProduct(:,i,j,l) * spAzeta % D(k,l)
            end do
         end do         ; end do          ; end do
!
!        Compute the curl
!        ----------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            Jai(1,i,j,k) = auxgrad(3,2,i,j,k) - auxgrad(2,3,i,j,k)
            Jai(2,i,j,k) = auxgrad(1,3,i,j,k) - auxgrad(3,1,i,j,k)
            Jai(3,i,j,k) = auxgrad(2,1,i,j,k) - auxgrad(1,2,i,j,k)
         end do         ; end do          ; end do
!
!        Assign to the first coordinate of each metrics
!        ----------------------------------------------
         self % jGradXi(1,:,:,:)   = -Jai(1,:,:,:)
         self % jGradEta(1,:,:,:)  = -Jai(2,:,:,:)
         self % jGradZeta(1,:,:,:) = -Jai(3,:,:,:)
!
!        *****************************************
!        Compute the y-coordinates of the mappings
!        *****************************************
!
!        Compute coordinates combination
!        -------------------------------
         do k = 0, self % Nz    ; do j = 0, self % Ny  ; do i = 0, self % Nx
            coordsProduct(:,i,j,k) = self % x(1,i,j,k) * grad_x(3,:,i,j,k)
         end do            ; end do          ; end do
!
!        Compute its gradient
!        --------------------
         auxgrad = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               auxgrad(:,1,i,j,k) = auxgrad(:,1,i,j,k) + coordsProduct(:,l,j,k) * spAxi % D(i,l)
            end do
      
            do l = 0, self % Ny
               auxgrad(:,2,i,j,k) = auxgrad(:,2,i,j,k) + coordsProduct(:,i,l,k) * spAeta % D(j,l)
            end do

            do l = 0, self % Nz
               auxgrad(:,3,i,j,k) = auxgrad(:,3,i,j,k) + coordsProduct(:,i,j,l) * spAzeta % D(k,l)
            end do
         end do         ; end do          ; end do
!
!        Compute the curl
!        ----------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            Jai(1,i,j,k) = auxgrad(3,2,i,j,k) - auxgrad(2,3,i,j,k)
            Jai(2,i,j,k) = auxgrad(1,3,i,j,k) - auxgrad(3,1,i,j,k)
            Jai(3,i,j,k) = auxgrad(2,1,i,j,k) - auxgrad(1,2,i,j,k)
         end do         ; end do          ; end do
!
!        Assign to the second coordinate of each metrics
!        -----------------------------------------------
         self % jGradXi(2,:,:,:)   = -Jai(1,:,:,:)
         self % jGradEta(2,:,:,:)  = -Jai(2,:,:,:)
         self % jGradZeta(2,:,:,:) = -Jai(3,:,:,:)
!
!        *****************************************
!        Compute the z-coordinates of the mappings
!        *****************************************
!
!        Compute coordinates combination
!        -------------------------------
         do k = 0, self % Nz    ; do j = 0, self % Ny  ; do i = 0, self % Nx
            coordsProduct(:,i,j,k) = self % x(2,i,j,k) * grad_x(1,:,i,j,k)
         end do            ; end do          ; end do
!
!        Compute its gradient
!        --------------------
         auxgrad = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               auxgrad(:,1,i,j,k) = auxgrad(:,1,i,j,k) + coordsProduct(:,l,j,k) * spAxi % D(i,l)
            end do
      
            do l = 0, self % Ny
               auxgrad(:,2,i,j,k) = auxgrad(:,2,i,j,k) + coordsProduct(:,i,l,k) * spAeta % D(j,l)
            end do

            do l = 0, self % Nz
               auxgrad(:,3,i,j,k) = auxgrad(:,3,i,j,k) + coordsProduct(:,i,j,l) * spAzeta % D(k,l)
            end do
         end do         ; end do          ; end do
!
!        Compute the curl
!        ----------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            Jai(1,i,j,k) = auxgrad(3,2,i,j,k) - auxgrad(2,3,i,j,k)
            Jai(2,i,j,k) = auxgrad(1,3,i,j,k) - auxgrad(3,1,i,j,k)
            Jai(3,i,j,k) = auxgrad(2,1,i,j,k) - auxgrad(1,2,i,j,k)
         end do         ; end do          ; end do
!
!        Assign to the third coordinate of each metrics
!        ----------------------------------------------
         self % jGradXi(3,:,:,:)   = -Jai(1,:,:,:)
         self % jGradEta(3,:,:,:)  = -Jai(2,:,:,:)
         self % jGradZeta(3,:,:,:) = -Jai(3,:,:,:)
!
!        ********************
!        Compute the Jacobian
!        ********************
!
         do k = 0, self % Nz  ; do j = 0, self % Ny   ; do i = 0, self % Nx
            self % jacobian(i,j,k) = jacobian3D(a1 = grad_x(:,1,i,j,k), &
                                                a2 = grad_x(:,2,i,j,k), &
                                                a3 = grad_x(:,3,i,j,k)   )
         end do               ; end do                ; end do

      end subroutine computeMetricTermsConservativeForm
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeMetricTermsCrossProductForm(self, spAxi, spAeta, spAzeta, mapper)       
!
!     -----------------------------------------------
!     Compute the metric terms in cross product form 
!     -----------------------------------------------
!
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry)   , intent(inout) :: self
         TYPE(NodalStorage)     , intent(in)    :: spAxi
         TYPE(NodalStorage)     , intent(in)    :: spAeta
         TYPE(NodalStorage)     , intent(in)    :: spAzeta
         TYPE(TransfiniteHexMap), intent(in)    :: mapper
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: i,j,k,l
         INTEGER       :: Nx, Ny, Nz
         REAL(KIND=RP) :: grad_x(3,3)         

         Nx = spAxi % N
         Ny = spAeta % N
         Nz = spAzeta % N
         
         DO k = 0, Nz
            DO j = 0,Ny
               DO i = 0,Nx
                  grad_x = 0.0_RP
                  do l = 0, Nx
                     grad_x(:,1) = grad_x(:,1) + self % X(:,l,j,k) * spAxi % D(i,l)
                  end do
                  do l = 0, Ny
                     grad_x(:,2) = grad_x(:,2) + self % X(:,i,l,k) * spAeta % D(j,l)
                  end do
                  do l = 0, Nz
                     grad_x(:,3) = grad_x(:,3) + self % X(:,i,j,l) * spAzeta % D(k,l)
                  end do
                 
                  CALL vCross( grad_x(:,2), grad_x(:,3), self % jGradXi  (:,i,j,k))
                  CALL vCross( grad_x(:,3), grad_x(:,1), self % jGradEta (:,i,j,k))
                  CALL vCross( grad_x(:,1), grad_x(:,2), self % jGradZeta(:,i,j,k))

                  self % jacobian(i,j,k) = jacobian3D(a1 = grad_x(:,1),a2 = grad_x(:,2),a3 = grad_x(:,3))
               END DO   
            END DO   
         END DO  

      END SUBROUTINE computeMetricTermsCrossProductForm
!
!///////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!     Returns the jacobian of the transformation computed from
!!     the three co-variant coordinate vectors.
!-------------------------------------------------------------------------------
!
      FUNCTION jacobian3D(a1,a2,a3)
!
      USE SMConstants
      IMPLICIT NONE

      REAL(KIND=RP)               :: jacobian3D
      REAL(KIND=RP), DIMENSION(3) :: a1,a2,a3,v
!
      CALL vCross(a2,a3,v)
      jacobian3D = vDot(a1,v)

      END FUNCTION jacobian3D
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns in result the cross product u x v
!-------------------------------------------------------------------------------
!
      SUBROUTINE vCross(u,v,result)
!
      IMPLICIT NONE
      
      REAL(KIND=RP), DIMENSION(3) :: u,v,result

      result(1) = u(2)*v(3) - v(2)*u(3)
      result(2) = u(3)*v(1) - v(3)*u(1)
      result(3) = u(1)*v(2) - v(1)*u(2)

      END SUBROUTINE vCross
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns the dot product u.v
!-------------------------------------------------------------------------------
!
      FUNCTION vDot(u,v)
!
      IMPLICIT NONE
      
      REAL(KIND=RP)               :: vDot
      REAL(KIND=RP), DIMENSION(3) :: u,v

      vDot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      END FUNCTION vDot
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns the 2-norm of u
!-------------------------------------------------------------------------------
!
      FUNCTION vNorm(u)
!
      IMPLICIT NONE
      
      REAL(KIND=RP)               :: vNorm
      REAL(KIND=RP), DIMENSION(3) :: u

      vNorm = SQRT(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))

      END FUNCTION vNorm
!
!///////////////////////////////////////////////////////////////////////////////
!
      
END Module MappedGeometryClass
