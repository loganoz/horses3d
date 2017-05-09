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
      INTEGER, PARAMETER :: LEFT   = 1, RIGHT  = 2, TOP  = 2, BOTTOM  = 1
      INTEGER, PARAMETER :: FRONT  = 1, BACK   = 2
      INTEGER, PARAMETER :: ELEFT  = 6, ERIGHT = 4, ETOP = 5, EBOTTOM = 3
      INTEGER, PARAMETER :: EFRONT = 1, EBACK  = 2 
!
!     -----
!     Class
!     -----
!
      TYPE MappedGeometry
            INTEGER                                         :: Nx, Ny, Nz                    ! Polynomial order
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: jGradXi, jGradEta, jGradZeta  ! 
            REAL(KIND=RP), DIMENSION(:,:,:)   , ALLOCATABLE :: x, y, z                       ! Position of points in absolute coordinates
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: xb                            ! 
            REAL(KIND=RP), DIMENSION(:,:,:)   , ALLOCATABLE :: jacobian, scal                ! 
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: normal                        ! Normal vectors on the nodes that are on the boundaries of the element (allocated with maximum order)
            
            CONTAINS
            
            PROCEDURE :: construct => ConstructMappedGeometry
            PROCEDURE :: destruct  => DestructMappedGeometry
      END TYPE MappedGeometry

!
!  ========
   CONTAINS 
!  ========
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ConstructMappedGeometry( self, spA, mapper )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
      CLASS(MappedGeometry)    :: self
      TYPE(TransfiniteHexMap)  :: mapper
      TYPE(NodalStorage)       :: spA
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: Nx, Ny, Nz, Nmax
      INTEGER       :: i, j, k
      REAL(KIND=RP) :: nrm
      REAL(KIND=RP) :: grad_x(3,3), jGrad(3)
      LOGICAL       :: useCrossProductMetrics = .TRUE. ! A switch for debugging purposes
                                                       ! So far, anisotropic representation is only implemented for cross product metrics
!
!     -----------
!     Allocations
!     -----------
!
      Nx        = spA % Nx
      Ny        = spA % Ny
      Nz        = spA % Nz
      Nmax      = MAX(Nx,Ny,Nz)
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
      
      ALLOCATE( self % JGradXi  (3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % JGradEta (3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % JGradZeta(3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % jacobian   (0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % x        (3,0:Nx,0:Ny,0:Nz)    )
      
      ALLOCATE( self % xb    (3,0:Nmax,0:Nmax,6)     )
      ALLOCATE( self % normal(3,0:Nmax,0:Nmax,6) )
      ALLOCATE( self % scal    (0:Nmax,0:Nmax,6)     )
!
!     --------------------------
!     Compute interior locations
!     --------------------------
!
      DO k = 0, Nz
         DO j= 0, Ny       
            DO i = 0,Nx 
               self % x(:,i,j,k) = mapper %  transfiniteMapAt([spA % xi(i), spA % eta(j), spa % zeta(k)])
            END DO
         END DO
      END DO
!
!     ----------------------
!     Compute face locations
!     ----------------------
!
      ! y-z planes
      DO j = 0, Nz
         DO i = 0, Ny
            self % xb(:,i,j,ELEFT)   = mapper % transfiniteMapAt([-1.0_RP    , spA % eta(i), spa % zeta(j)])
            self % xb(:,i,j,ERIGHT)  = mapper % transfiniteMapAt([ 1.0_RP    , spA % eta(i), spa % zeta(j)])
         END DO
      END DO 
      
      ! x-y planes
      DO j = 0, Ny
         DO i = 0, Nx
            self % xb(:,i,j,EBOTTOM) = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j),    -1.0_RP   ])
            self % xb(:,i,j,ETOP)    = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j),     1.0_RP   ])
         END DO
      END DO 
      
      ! x-z planes
      DO j = 0, Nz
         DO i = 0, Nx
            self % xb(:,i,j,EBACK)   = mapper % transfiniteMapAt([spA % xi(i),  1.0_RP     , spa % zeta(j)  ])
            self % xb(:,i,j,EFRONT)  = mapper % transfiniteMapAt([spA % xi(i), -1.0_RP     , spa % zeta(j)  ])
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
      IF (useCrossProductMetrics .OR. isHex8(mapper)) THEN 
      
         CALL computeMetricTermsCrossProductForm(self, spA, mapper)
!
!     ----------------
!     Boundary Normals - Must be evaluated at the boundaries!
!     ----------------
!
         ! y-z planes
         DO j = 0, Nz
            DO i = 0, Ny
!
!           ---------
!           Left face
!           ---------
!
               grad_x = mapper % metricDerivativesAt([-1.0_RP    , spA % eta(i), spa % zeta(j)])
               CALL vCross(grad_x(:,2), grad_x(:,3), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,ELEFT) = -jGrad/nrm
               self % scal(i,j,ELEFT)     = nrm
!
!           ----------
!           Right face
!           ----------
!
               grad_x = mapper % metricDerivativesAt([ 1.0_RP    , spA % eta(i), spa % zeta(j)])
               CALL vCross(grad_x(:,2), grad_x(:,3), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,ERIGHT) = jGrad/nrm
               self % scal(i,j,ERIGHT)     = nrm     
            END DO
         END DO
         
         ! x-y planes
         DO j = 0, Ny
            DO i = 0, Nx
!
!           -----------
!           bottom face
!           -----------
!
               grad_x = mapper % metricDerivativesAt([spA % xi(i), spA % eta(j),    -1.0_RP   ])
               CALL vCross(grad_x(:,1), grad_x(:,2), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,EBOTTOM) = -jGrad/nrm
               self % scal(i,j,EBOTTOM)     = nrm
!
!           --------
!           top face
!           --------
!
               grad_x = mapper % metricDerivativesAt([spA % xi(i), spA % eta(j),     1.0_RP   ])
               CALL vCross(grad_x(:,1), grad_x(:,2), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,ETOP) = jGrad/nrm
               self % scal(i,j,ETOP)     = nrm
            END DO
         END DO
         
         ! x-z planes
         DO j = 0, Nz
            DO i = 0, Nx
!
!           ----------
!           front face
!           ----------
!
               grad_x = mapper % metricDerivativesAt([spA % xi(i), -1.0_RP     , spa % zeta(j)  ])
               CALL vCross(grad_x(:,3), grad_x(:,1), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,EFRONT) = -jGrad/nrm
               self % scal(i,j,EFRONT)     = nrm
!
!           ---------
!           back face
!           ---------
!
               grad_x = mapper % metricDerivativesAt([spA % xi(i),  1.0_RP     , spa % zeta(j)  ])
               CALL vCross(grad_x(:,3), grad_x(:,1), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,EBACK) = jGrad/nrm
               self % scal(i,j,EBACK)     = nrm
              
            END DO
         END DO
         
      ELSE
       
         CALL computeMetricTermsConservativeForm(self, spA, mapper)
      
      ENDIF       

   END SUBROUTINE ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructMappedGeometry(self)
         IMPLICIT NONE 
         CLASS(MappedGeometry) :: self
         DEALLOCATE( self % jGradXi, self % jGradEta, self % jGradZeta, self % x, self % jacobian )
         DEALLOCATE( self % xb )
         DEALLOCATE( self % normal, self % scal )
      END SUBROUTINE DestructMappedGeometry
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeMetricTermsConservativeForm(self, spA, mapper)  
!
!     -----------------------------------------------
!     Compute the metric terms in conservative form 
!     -----------------------------------------------
!
      IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry)    :: self
         TYPE(NodalStorage)      :: spA
         TYPE(TransfiniteHexMap) :: mapper      
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: xiArray(0:spA % N, 3), w(0:spA % N, 3)
      REAL(KIND=RP) :: xi(3)
      
      real(KIND=RP) :: mappedxi(0:spA % N, 3)
      
      REAL(KIND=RP) :: grad_x(3, 3, 0:spA % N, 0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: xGauss(3, 0:spA % N, 0:spA % N, 0:spA % N)
      
      REAL(KIND=RP) :: xiDermat  (0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: etaDerMat (0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: zetaDerMat(0:spA % N, 0:spA % N)
      
      REAL(KIND=RP) :: IdentityMatrix(0:spA % N, 0:spA % N)
      
      REAL(KIND=RP) :: tArray(0:spA % N, 0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: dArray(0:spA % N, 0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: vArray(0:spA % N, 0:spA % N, 0:spA % N)
      
      REAL(KIND=RP) :: jGradXi  (3, 0:spA % N, 0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: jGradEta (3, 0:spA % N, 0:spA % N, 0:spA % N)
      REAL(KIND=RP) :: jGradZeta(3, 0:spA % N, 0:spA % N, 0:spA % N)
      
      REAL(KIND=RP) :: jGrad (3)
      REAL(KIND=RP) :: nrm
      
      REAL(KIND=RP) :: piN
      
      INTEGER :: i, j, k, m, l, n, noGauss
      
      INTEGER :: polOrder(3)
      
      REAL(KIND=RP), ALLOCATABLE :: xiInterpMat  (:,:)
      REAL(KIND=RP), ALLOCATABLE :: etaInterpMat (:,:)
      REAL(KIND=RP), ALLOCATABLE :: zetaInterpMat(:,:)      
!
!     ---------------------------
!     A convenience mapping array
!     ---------------------------
!
      INTEGER, DIMENSION(0:4) :: iCycle = (/3,1,2,3,1/)

      polOrder(:) = spa % N
      noGauss     = spA % N          
!
!     -------------------------------
!     Mappedxi is Legendre Gauss grid
!     -------------------------------
!
      DO k = 1,3
         DO j = 0, noGauss
            mappedxi(j,1) = spA % xi(j)
            mappedxi(j,2) = spA % eta(j)
            mappedxi(j,3) = spA % zeta(j)
         END DO
      END DO      
!
!     -------------------------------------------
!     Compute the mesh on the Chebyshev Lobatto 
!     grid and compute the gradients on that mesh
!     -------------------------------------------
!
      piN = PI/(noGauss)
      DO k = 1,3
         DO n = 0, noGauss
             xiArray(n,k) = - COS((n)*piN)
         END DO
      END DO
       
      DO l = 0,noGauss
         xi(3) = xiArray(l,3)
         DO m = 0,noGauss
            xi(2) = xiArray(m,2)
            DO n = 0,noGauss
               xi(1) = xiArray(n,1)
               CALL GeneralHexGradAndMap( xi, xGauss(:,n,m,l), grad_x(:,:,n,m,l), mapper%corners, mapper % faces )
            END DO
         END DO
      END DO

      CALL PolynomialDerivativeMatrix( noGauss, xiArray(:,1), xiDerMat )
      CALL PolynomialDerivativeMatrix( noGauss, xiArray(:,2), etaDerMat )
      CALL PolynomialDerivativeMatrix( noGauss, xiArray(:,3), zetaDerMat )
!
!     -----------------------------------------------------
!     Now compute metric terms at each grid point
!     This computes the jGradXi terms in conservative form.
!     See the notes for the derivations.
!     -----------------------------------------------------
!
!     ----------
!     First term
!     ----------
!
      iLoop: DO i = 1,3
         jLoop : DO j = 1,3

            tArray = xGauss(iCycle(j-1),:,:,:)*grad_x(iCycle(j+1),iCycle(i+1),:,:,:)

            SELECT CASE (i)
               CASE (1)
                  DO l = 0, noGauss
                     DO m = 0, noGauss   
                        CALL MatrixMultiplyDeriv(tArray(m,l,:), dArray(m,l,:), zetaDerMat, noGauss, transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradXi(j,:,:,:) = dArray
               CASE (2)
                  DO l = 0, noGauss
                     DO m = 0, noGauss   
                        CALL MatrixMultiplyDeriv(tArray(:,m,l), dArray(:,m,l), xiDerMat, noGauss, transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradEta(j,:,:,:) = dArray
               CASE (3)
                  DO l = 0, noGauss
                     DO m = 0, noGauss   
                        CALL MatrixMultiplyDeriv(tArray(m,:,l), dArray(m,:,l), etaDerMat, noGauss, transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO            
                  jGradZeta(j,:,:,:) = dArray
            END SELECT

         END DO jLoop
      END DO iLoop
!
!     -----------
!     Second term
!     -----------
!
      iLoop2: DO i = 1,3
         jLoop2 : DO j = 1,3

            tArray = xGauss(iCycle(j-1),:,:,:)*grad_x(iCycle(j+1),iCycle(i-1),:,:,:)

            SELECT CASE (i)
               CASE (1)
                  DO l = 0, noGauss
                     DO m = 0, noGauss   
                        CALL MatrixMultiplyDeriv(tArray(m,:,l), dArray(m,:,l), etaDerMat, noGauss, transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradXi(j,:,:,:) = jGradXi (j,:,:,:) - dArray
               CASE (2)
                  DO l = 0, noGauss
                     DO m = 0, noGauss   
                        CALL MatrixMultiplyDeriv(tArray(m,l,:), dArray(m,l,:), zetaDerMat, noGauss, transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradEta(j,:,:,:) = jGradEta(j,:,:,:) - dArray
               CASE (3)
                  DO l = 0, noGauss
                     DO m = 0, noGauss   
                        CALL MatrixMultiplyDeriv(tArray(:,m,l), dArray(:,m,l), xiDerMat, noGauss, transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradZeta(j,:,:,:) = jGradZeta(j,:,:,:) - dArray
            END SELECT

         END DO jLoop2
      END DO iLoop2    
      

!
!     ------------------------------------
!     Interpolate back onto the Gauss grid
!     ------------------------------------
!
      ALLOCATE( xiInterpMat  (0:noGauss,0:noGauss) )
      ALLOCATE( etaInterpMat (0:noGauss,0:noGauss) )
      ALLOCATE( zetaInterpMat(0:noGauss,0:noGauss) )
      
      CALL BarycentricWeights( noGauss, xiArray(:,1), w(:,1))
      CALL PolynomialInterpolationMatrix( noGauss, noGauss, xiArray(:,1), w(:,1), mappedxi(:,1), xiInterpmat)
      CALL BarycentricWeights( noGauss, xiArray(:,2), w(:,2))
      CALL PolynomialInterpolationMatrix( noGauss, noGauss, xiArray(:,2), w(:,2), mappedxi(:,2), etaInterpmat)
      CALL BarycentricWeights( noGauss, xiArray(:,3), w(:,3))
      CALL PolynomialInterpolationMatrix( noGauss, noGauss, xiArray(:,3), w(:,3), mappedxi(:,3), zetaInterpmat)      
      
      DO k = 1,3      
         
         tArray(:,:,:) = jGradXi(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         jGradXi(k,:,:,:) = vArray(:,:,:)
         
         tArray(:,:,:) = jGradEta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         jGradEta(k,:,:,:) = vArray(:,:,:)
         
         tArray(:,:,:) = jGradZeta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         jGradZeta(k,:,:,:) = vArray(:,:,:)

      END DO

      self % jGradXi   = jGradXi
      self % jGradEta  = jGradEta
      self % jGradZeta = jGradZeta              

!
!     ----------------------------------
!     Compute the jacobian at each point
!     ----------------------------------
!
      DO l = 0,noGauss
         DO m = 0,noGauss
            DO n = 0,noGauss
               tArray(n,m,l) = jacobian3D( grad_x(:,1,n,m,l), grad_x(:,2,n,m,l), grad_x(:,3,n,m,l) )
            END DO
         END DO
      END DO
      CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )

      self % jacobian = vArray

!
!     ----------------------------
!     Boundary Normals Computation
!     ----------------------------
!
!
!     ----------------------------------------
!     Interpolate to the faces. This is done 
!     by performing a 3D interpolation with 
!     Lobatto nodes in the direction of 
!     interpolation and Gauss nodes in the
!     other two. 
!     ----------------------------------------

      IdentityMatrix = 0.0_RP
      do i = 0, noGauss 
         IdentityMatrix(i,i) = 1.0_RP         
      enddo 

      CALL BarycentricWeights( noGauss, mappedxi(:,1), w(:,1))
      CALL PolynomialInterpolationMatrix( noGauss, noGauss, mappedxi(:,1), w(:,1), xiArray(:,1), xiInterpmat)
      CALL BarycentricWeights( noGauss, mappedxi(:,2), w(:,2))
      CALL PolynomialInterpolationMatrix( noGauss, noGauss, mappedxi(:,2), w(:,2), xiArray(:,2), etaInterpmat)
      CALL BarycentricWeights( noGauss, mappedxi(:,3), w(:,3))
      CALL PolynomialInterpolationMatrix( noGauss, noGauss, mappedxi(:,3), w(:,3), xiArray(:,3), zetaInterpmat)   
      

      DO k = 1,3
         tArray(:,:,:) = jGradXi(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, IdentityMatrix, IdentityMatrix )
         jGradXi(k,:,:,:) = vArray(:,:,:)


         tArray(:,:,:) = jGradEta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, IdentityMatrix, etaInterpMat, IdentityMatrix )
         jGradEta(k,:,:,:) = vArray(:,:,:)

         tArray(:,:,:) = jGradZeta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, IdentityMatrix, IdentityMatrix, zetaInterpMat )
         jGradZeta(k,:,:,:) = vArray(:,:,:)
      END DO
  
!
!     ---------------------------
!     Evaluation of the normals
!     To Do: divide by the sign 
!     of the Jacobian (see notes)
!     ---------------------------
!
      DO j = 0, noGauss
         DO i = 0, noGauss
!
!           ---------
!           Left face
!           ---------
!
            jGrad(:) = jGradXi (:,0,i,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,ELEFT) = -jGrad/nrm
            self % scal(i,j,ELEFT)     = nrm 
!
!           ----------
!           Right face
!           ----------
!
            jGrad(:) = jGradXi(:,noGauss,i,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,ERIGHT) = jGrad/nrm
            self % scal(i,j,ERIGHT)     = nrm 
!
!           -----------
!           bottom face
!           -----------
!
            jGrad(:) = jGradZeta(:,i,j,0)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,EBOTTOM) = -jGrad/nrm
            self % scal(i,j,EBOTTOM)     = nrm 
!
!           --------
!           top face
!           --------
!
            jGrad(:) = jGradZeta(:,i,j,noGauss)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,ETOP) = jGrad/nrm
            self % scal(i,j,ETOP)     = nrm 
!
!           ----------
!           front face
!           ----------
!
            jGrad(:) = jGradEta(:,i,0,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,EFRONT) = -jGrad/nrm
            self % scal(i,j,EFRONT)     = nrm 
!
!           ---------
!           back face
!           ---------
!
            jGrad(:) = jGradEta(:,i,noGauss,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,EBACK) = jGrad/nrm
            self % scal(i,j,EBACK)     = nrm 
            
         ENDDO
      ENDDO 
      
      DEALLOCATE(xiInterpMat)
      DEALLOCATE(etaInterpMat)
      DEALLOCATE(zetaInterpMat)
      
      END SUBROUTINE computeMetricTermsConservativeForm
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeMetricTermsConservativeDirectly(self, spA, mapper)  
!
!        -------------------------------------------------------------
!        Compute the metric terms by explicitly writing out the terms:
!
!        Ja1 = [(I(Y_eta*Z))_zeta - (I(Y_zeta*Z))_eta] x^ +
!              [(I(Z_eta*X))_zeta - (I(Z_zeta*X))_eta] y^ +
!              [(I(X_eta*Y))_zeta - (I(X_zeta*Y))_eta] z^ 
!
!        Ja2 = [(I(Y_zeta*Z))_xi - (I(Y_xi*Z))_zeta] x^ +
!              [(I(Z_zeta*X))_xi - (I(Z_xi*X))_zeta] y^ +
!              [(I(X_zeta*Y))_xi - (I(X_xi*Y))_zeta] z^ 
!
!        Ja3 = [(I(Y_xi*Z))_eta - (I(Y_eta*Z))_xi] x^ +
!              [(I(Z_xi*X))_eta - (I(Z_eta*X))_xi] y^ +
!              [(I(X_xi*Y))_eta - (I(X_eta*Y))_xi] z^ 
!
!        -------------------------------------------------------------
!
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry)    :: self
         TYPE(NodalStorage)      :: spA
         TYPE(TransfiniteHexMap) :: mapper
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: i,j,k, N
         REAL(KIND=RP) :: grad_x(3,3)
         N = spA % N
!
!        --------------------
!        Compute the Jacobian
!        --------------------
!
         DO k = 0, N
            DO j = 0,N
               DO i = 0,N
                  grad_x                 = mapper % metricDerivativesAt([spA % xi(i), spA % eta(j), spA % zeta(k)])
                  self % jacobian(i,j,k) = jacobian3D(a1 = grad_x(:,1),a2 = grad_x(:,2),a3 = grad_x(:,3))
               END DO   
            END DO   
         END DO  

         
         
      END SUBROUTINE computeMetricTermsConservativeDirectly
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeMetricTermsCrossProductForm(self, spA, mapper)       
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
         TYPE(MappedGeometry)    :: self
         TYPE(NodalStorage)      :: spA
         TYPE(TransfiniteHexMap) :: mapper
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: i,j,k
         INTEGER       :: Nx, Ny, Nz
         REAL(KIND=RP) :: grad_x(3,3)         
         Nx = spA % Nx
         Ny = spA % Ny
         Nz = spA % Nz
         
         DO k = 0, Nz
            DO j = 0,Ny
               DO i = 0,Nx
                  grad_x = mapper % metricDerivativesAt([spA % xi(i), spA % eta(j), spA % zeta(k)])
                 
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
