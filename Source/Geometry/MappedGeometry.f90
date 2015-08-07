!
!////////////////////////////////////////////////////////////////////////
!
!      MappedGeometry.f95
!      Created: 2008-06-19 15:58:02 -0400 
!      By: David Kopriva  
!
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
            INTEGER                                         :: N
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: jGradXi, jGradEta, jGradZeta
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: x
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: xb
            REAL(KIND=RP), DIMENSION(:,:,:)   , ALLOCATABLE :: jacobian, scal
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: normal
            
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
      INTEGER       :: N
      INTEGER       :: i, j, k
      REAL(KIND=RP) :: nrm
      REAL(KIND=RP) :: grad_x(3,3), jGrad(3)
!
!     -----------
!     Allocations
!     -----------
!
      N        = spA % N
      self % N = N
      
      ALLOCATE( self % JGradXi  (3,0:N,0:N,0:N) )
      ALLOCATE( self % JGradEta (3,0:N,0:N,0:N) )
      ALLOCATE( self % JGradZeta(3,0:N,0:N,0:N) )
      ALLOCATE( self % jacobian(0:N,0:N,0:N) )
      
      ALLOCATE( self % x(3,0:N,0:N,0:N)    )
      ALLOCATE( self % xb(3,0:N,0:N,6)     )
      ALLOCATE( self % normal(3,0:N,0:N,6) )
      ALLOCATE( self % scal(0:N,0:N,6)     )
!
!     --------------------------
!     Compute interior locations
!     --------------------------
!
      DO k = 0, N
         DO j= 0, N       
            DO i = 0,N 
               self % x(:,i,j,k) = mapper %  transfiniteMapAt([spA % xi(i), spA % eta(j), spa % zeta(k)])
            END DO
         END DO
      END DO
!
!     ----------------------
!     Compute face locations
!     ----------------------
!
      DO j = 0, N
         DO i = 0, N
            self % xb(:,i,j,ELEFT)   = mapper % transfiniteMapAt([-1.0_RP    , spA % eta(i), spa % zeta(j)])
            self % xb(:,i,j,ERIGHT)  = mapper % transfiniteMapAt([ 1.0_RP    , spA % eta(i), spa % zeta(j)])
            self % xb(:,i,j,EBOTTOM) = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j),    -1.0_RP   ])
            self % xb(:,i,j,ETOP)    = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j),     1.0_RP   ])
            self % xb(:,i,j,EBACK)   = mapper % transfiniteMapAt([spA % xi(i),  1.0_RP     , spa % zeta(j)  ])
            self % xb(:,i,j,EFRONT)  = mapper % transfiniteMapAt([spA % xi(i), -1.0_RP     , spa % zeta(j)  ])
         END DO
      END DO 
!
!     ------------
!     Metric terms
!     ------------
!
      IF (isHex8(mapper)) THEN 
         CALL computeMetricTermsCrossProductForm(self, spA, mapper)
      ELSE 
         !CALL computeMetricTermsCrossProductForm(self, spA, mapper)
         !self % jacobian = 0.0_RP
         CALL computeMetricTermsConservativeForm(self, spA, mapper)
      ENDIF 
!
!     ----------------
!     Boundary Normals - Must be evaluated at the boundaries!
!     ----------------
!
      DO j = 0, N
         DO i = 0, N
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
      REAL(KIND=RP) :: xiArray( spA % N + 1,3)
      REAL(KIND=RP) :: xi(3)
      REAL(KIND=RP) :: corners(3,8)
      
      REAL(KIND=RP) :: grad_x(3,3,spA % N + 1,spA % N + 1,spA % N + 1)
      REAL(KIND=RP) :: xGauss(3,  spA % N + 1,spA % N + 1,spA % N + 1)
      
      REAL(KIND=RP) :: xiDermat  (spA % N + 1, spA % N + 1)
      REAL(KIND=RP) :: etaDerMat (spA % N + 1, spA % N + 1)
      REAL(KIND=RP) :: zetaDerMat(spA % N + 1, spA % N + 1)
      
      REAL(KIND=RP) :: tArray( spA % N + 1, spA % N + 1, spA % N + 1 )
      REAL(KIND=RP) :: dArray( spA % N + 1, spA % N + 1, spA % N + 1 )
      REAL(KIND=RP) :: vArray( spA % N + 1, spA % N + 1, spA % N + 1 )
      
      REAL(KIND=RP) :: jGradXi  (3,spA % N + 1,spA % N + 1,spA % N + 1)
      REAL(KIND=RP) :: jGradEta (3,spA % N + 1,spA % N + 1,spA % N + 1)
      REAL(KIND=RP) :: jGradZeta(3,spA % N + 1,spA % N + 1,spA % N + 1)
      
      REAL(KIND=RP) :: piN
      
      INTEGER :: i,j,k,n,m,l
      
      INTEGER :: polOrder(3)
      real(KIND=RP) :: jacInv( spA % N + 1, spA % N + 1, spA % N + 1 )
      
      REAL(KIND=RP), ALLOCATABLE :: xiInterpMat  (:,:)
      REAL(KIND=RP), ALLOCATABLE :: etaInterpMat (:,:)
      REAL(KIND=RP), ALLOCATABLE :: zetaInterpMat(:,:)      
!
!     ---------------------------
!     A convenience mapping array
!     ---------------------------
!
      INTEGER, DIMENSION(0:4) :: iCycle = (/3,1,2,3,1/)

      polOrder(:) = spa % N + 1
      
      spA % N = spA % N + 1 
!
!     --------------------------------------
!     Compute the mesh on the Lobatto grid
!     and compute the gradients on that mesh
!     --------------------------------------
!
      DO k = 1,3
         piN = PI/(spA % N-1)
         DO n = 1, spA % N
            xiArray(n,k) = 0.5_RP*(1.0_RP - COS((n-1)*piN))
         END DO
      END DO
      DO j = 1,8
         DO i = 1,3
            corners(i,j) = mapper%corners(i,j)
         END DO
      END DO
      DO l = 1,spA % N
         xi(3) = xiArray(l,3)
         DO m = 1,spA % N
            xi(2) = xiArray(m,2)
            DO n = 1,spA % N
               xi(1) = xiArray(n,1)
               CALL GeneralHexGradAndMap( xi, xGauss(:,n,m,l), grad_x(:,:,n,m,l), corners, mapper % faces )
            END DO
         END DO
      END DO
      
      
      CALL SetStandardDerivativeMatrix( spA % N, xiArray(:,1), xiDerMat)
      CALL SetStandardDerivativeMatrix( spA % N, xiArray(:,2), etaDerMat)
      CALL SetStandardDerivativeMatrix( spA % N, xiArray(:,3), zetaDerMat)
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

            DO l = 1,spA % N
               DO m = 1,spA % N
                  DO n = 1,spA % N
                     tArray(n,m,l) = xGauss(iCycle(j-1),n,m,l)*grad_x(iCycle(j+1),iCycle(i+1),n,m,l)
                  END DO
               END DO
            END DO

            SELECT CASE (i)
               CASE (1)
                  CALL MMMultiply3D3( zetaDerMat, polOrder, tArray, dArray )
                  DO l = 1,spA % N
                     DO m = 1,spA % N
                        DO n = 1,spA % N
                           jGradXi(j,n,m,l) = dArray(n,m,l)
                        END DO
                     END DO
                  END DO
               CASE (2)
                  CALL MMMultiply3D1( xiDerMat, polOrder, tArray, dArray )
                  DO l = 1,spA % N
                     DO m = 1,spA % N
                        DO n = 1,spA % N
                           jGradEta(j,n,m,l) = dArray(n,m,l)
                        END DO
                     END DO
                  END DO
               CASE (3)
                  CALL MMMultiply3D2( etaDerMat, polOrder, tArray, dArray )
                  DO l = 1,spA % N
                     DO m = 1,spA % N
                        DO n = 1,spA % N
                           jGradZeta(j,n,m,l) = dArray(n,m,l)
                        END DO
                     END DO
                  END DO
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

            DO l = 1,spA % N
               DO m = 1,spA % N
                  DO n = 1,spA % N
                     tArray(n,m,l) = xGauss(iCycle(j-1),n,m,l)*grad_x(iCycle(j+1),iCycle(i-1),n,m,l)
                  END DO
               END DO
            END DO

            SELECT CASE (i)
               CASE (1)
                  CALL MMMultiply3D2( etaDerMat, polOrder, tArray, dArray )
                  DO l = 1,spA % N
                     DO m = 1,spA % N
                        DO n = 1,spA % N
                           jGradXi(j,n,m,l) = jGradXi(j,n,m,l)-dArray(n,m,l)
                        END DO
                     END DO
                  END DO
               CASE (2)
                  CALL MMMultiply3D3( zetaDerMat, polOrder, tArray, dArray )
                  DO l = 1,spA % N
                     DO m = 1,spA % N
                        DO n = 1,spA % N
                           jGradEta(j,n,m,l) = jGradEta(j,n,m,l) - dArray(n,m,l)
                        END DO
                     END DO
                  END DO
               CASE (3)
                  CALL MMMultiply3D1( xiDerMat, polOrder, tArray, dArray )
                  DO l = 1,spA % N
                     DO m = 1,spA % N
                        DO n = 1,spA % N
                           jGradZeta(j,n,m,l) = jGradZeta(j,n,m,l) - dArray(n,m,l)
                        END DO
                     END DO
                  END DO
            END SELECT

         END DO jLoop2
      END DO iLoop2      
!
!     ------------------------------------
!     Interpolate back onto the gauss grid
!     ------------------------------------
!
      ALLOCATE( xiInterpMat  (spA % N,spA % N) )
      ALLOCATE( etaInterpMat (spA % N,spA % N) )
      ALLOCATE( zetaInterpMat(spA % N,spA % N) )
      CALL MakeInterpMatFromTo( xiInterpmat  , spA % N, xiArray(:,1), spA % N, spA%xi(:) )
      CALL MakeInterpMatFromTo( etaInterpmat , spA % N, xiArray(:,2), spA % N, spA%xi(:) )
      CALL MakeInterpMatFromTo( zetaInterpmat, spA % N, xiArray(:,3), spA % N, spA%xi(:) )
      
      DO k = 1,3
         DO l = 1,spA % N
            DO m = 1,spA % N
               DO n = 1,spA % N
                  tArray(n,m,l) = jGradXi(k,n,m,l)
               END DO
            END DO
         END DO
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         DO l = 1,spA % N
            DO m = 1,spA % N
               DO n = 1,spA % N
                  jGradXi(k,n,m,l) = vArray(n,m,l)
               END DO
            END DO
         END DO
         DO l = 1,spA % N
            DO m = 1,spA % N
               DO n = 1,spA % N
                  tArray(n,m,l) = jGradEta(k,n,m,l)
               END DO
            END DO
         END DO
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         DO l = 1,spA % N
            DO m = 1,spA % N
               DO n = 1,spA % N
                  jGradEta(k,n,m,l) = vArray(n,m,l)
               END DO
            END DO
         END DO
         DO l = 1,spA % N
            DO m = 1,spA % N
               DO n = 1,spA % N
                  tArray(n,m,l) = jGradZeta(k,n,m,l)
               END DO
            END DO
         END DO
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         DO l = 1,spA % N
            DO m = 1,spA % N
               DO n = 1,spA % N
                  jGradZeta(k,n,m,l) = vArray(n,m,l)
               END DO
            END DO
         END DO
      END DO
!
!     ----------------------------------
!     Compute the jacobian at each point
!     ----------------------------------
!
      DO l = 1,spA % N
         DO m = 1,spA % N
            DO n = 1,spA % N
               tArray(n,m,l) = jacobian3D( grad_x(:,1,n,m,l), grad_x(:,2,n,m,l), grad_x(:,3,n,m,l) )
            END DO
         END DO
      END DO
      CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
      DO l = 1,spA % N
         DO m = 1,spA % N
            DO n = 1,spA % N
               jacInv(n,m,l) = 1.0_RP/vArray(n,m,l)
            END DO
         END DO
      END DO
      
      do k = 1,3
         do l = 1,spA % N
            do m = 1,spA % N
               do n = 1,spA % N
                  self % jGradXi(k, n-1, m-1, l-1 ) = jGradXi(k,n,m,l)
                  self % jGradEta(k, n-1, m-1, l-1 ) = jGradEta(k,n,m,l)
                  self % jGradZeta(k, n-1, m-1, l-1 ) = jGradZeta(k,n,m,l)
               enddo 
            enddo
         enddo
      enddo 

      do l = 1,spA % N
         do m = 1,spA % N
            do n = 1,spA % N
                self % jacobian(n-1,m-1,l-1)     = vArray(n,m,l)
            enddo 
         enddo
      enddo      

      spA % N = spA % N - 1

      END SUBROUTINE computeMetricTermsConservativeForm
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
         INTEGER       :: i,j,k, N
         REAL(KIND=RP) :: grad_x(3,3)         
         N = spA % N
         
         DO k = 0, N
            DO j = 0,N
               DO i = 0,N
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
!-------------------------------------------------------------------------------
!!    Copied and pasted subroutines from DSEM for the conservative metric terms
!!    computation
!-------------------------------------------------------------------------------
!
!                                                                       
!///////////////////////////////////////////////////////////////////////
!
!-----------------------------------------------------------
!!    Set an interpolation matrix from old to new points
!-----------------------------------------------------------
!
      SUBROUTINE SetInterpolationMatrix( nOld, xOld, nNew, xNew, b )
!

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                                , INTENT(IN)  :: nOld,nNew
      REAL(KIND=RP), DIMENSION(nOld)     , INTENT(IN)  :: xOld
      REAL(KIND=RP), DIMENSION(nNew)     , INTENT(IN)  :: xNew
      REAL(KIND=RP), DIMENSION(nNew,nOld), INTENT(OUT) :: b
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j,k
!                                                                       
      DO k = 1,nNew 
         DO j = 1,nOld 
            b(k,j) = EvaluateLagrangePoly(j,xNew(k),nOld,xOld)
         END DO
      END DO
!                                                                       
      END SUBROUTINE SetInterpolationMatrix
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!     Set the derivative matrix for point to point differentiation
!-------------------------------------------------------------------------------
!
      SUBROUTINE SetStandardDerivativeMatrix(n,x,b) 
!

!
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                        , INTENT(IN)  :: n
      REAL(KIND=RP), DIMENSION(n), INTENT(IN)  :: x
      REAL(KIND=RP)              , INTENT(OUT) :: b(n,n)
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: eps, sum
      INTEGER           :: k,j
!      
      eps = EPSILON(eps)
!
      DO k = 1,n 
         sum = 0._RP 
         DO j = 1,n 
            IF ( j == k ) CYCLE 
            b(k,j) = EvaluateLagrangePolyDeriv( j, x(k), n, x )
            IF ( dabs( b(k,j) ) < eps )      b(k,j) = 0._RP 
            sum = sum - b(k,j) 
         END DO 
         b(k,k) = sum 
      END DO 
!
      RETURN 
   END SUBROUTINE SetStandardDerivativeMatrix
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute a square matrix times a vector where the vector is
!!    stored as the first dimension of a three-dimensional array.
!---------------------------------------------------------------------
!
      SUBROUTINE MMMultiply3D1(A,n,x,y)
!

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER          , DIMENSION(3)             , INTENT(IN)  :: n
      REAL(KIND=RP), DIMENSION(n(1),n(1))     , INTENT(IN)  :: A
      REAL(KIND=RP), DIMENSION(n(1),n(2),n(3)), INTENT(IN)  :: x
      
      REAL(KIND=RP), DIMENSION(n(1),n(2),n(3)), INTENT(OUT) :: y
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: i,j,k,s
!
      y = 0.0_RP
      
      DO k = 1,n(3)
         DO j = 1,n(2)
            DO s = 1,n(1)
               DO i = 1,n(1)
                  y(i,j,k) = y(i,j,k) + A(i,s)*x(s,j,k)
               END DO
            END DO
         END DO
      END DO
      
      END SUBROUTINE MMMultiply3D1
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute a square matrix times a vector where the vector is
!!    stored as the second dimension of a three-dimensional array.
!---------------------------------------------------------------------
!
      SUBROUTINE MMMultiply3D2(A,n,x,y)
!

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER          , DIMENSION(3)             , INTENT(IN)  :: n
      REAL(KIND=RP), DIMENSION(n(2),n(2))     , INTENT(IN)  :: A
      REAL(KIND=RP), DIMENSION(n(1),n(2),n(3)), INTENT(IN)  :: x
      
      REAL(KIND=RP), DIMENSION(n(1),n(2),n(3)), INTENT(OUT) :: y
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: i,j,k,s
!
      y = 0.0_RP
      
      DO k = 1,n(3)
         DO s = 1,n(2)
            DO j = 1,n(2)
               DO i = 1,n(1)
                  y(i,j,k) = y(i,j,k) + A(j,s)*x(i,s,k)
               END DO
            END DO
         END DO
      END DO
      
      END SUBROUTINE MMMultiply3D2
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute a square matrix times a vector where the vector is
!!    stored as the third dimension of a three-dimensional array.
!---------------------------------------------------------------------
!
      SUBROUTINE MMMultiply3D3(A,n,x,y)
!

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER          , DIMENSION(3)             , INTENT(IN)  :: n
      REAL(KIND=RP), DIMENSION(n(3),n(3))     , INTENT(IN)  :: A
      REAL(KIND=RP), DIMENSION(n(1),n(2),n(3)), INTENT(IN)  :: x
      
      REAL(KIND=RP), DIMENSION(n(1),n(2),n(3)), INTENT(OUT) :: y
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: i,j,k,s
!
      y = 0.0_RP
      
      DO s = 1,n(3)
         DO k = 1,n(3)
            DO j = 1,n(2)
               DO i = 1,n(1)
                  y(i,j,k) = y(i,j,k) + A(k,s)*x(i,j,s)
               END DO
            END DO
         END DO
      END DO
      
      END SUBROUTINE MMMultiply3D3
!                                                                       
!///////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute the matrix that will interpolate from nOld points at xOld
!!    to nNew points at xNew.
!---------------------------------------------------------------------
!
      SUBROUTINE MakeInterpMatFromTo(interpMat,nOld,xOld,nNew,xNew) 
!
!     -------------------
!     Date: June 28, 2002
!     -------------------
!

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                           , INTENT(IN) :: nOld,nNew
      REAL(KIND=RP), DIMENSION(nOld), INTENT(IN) :: xOld
      REAL(KIND=RP), DIMENSION(nNew), INTENT(IN) :: xNew
      REAL(KIND=RP), DIMENSION(nNew,nOld)        :: interpMat

      INTEGER           :: j,k
!                                                                       
      DO k = 1,nNew 
         DO j = 1,nOld 
            interpMat(k,j) = EvaluateLagrangePoly(j,xNew(k),nOld,xOld)
         END DO
      END DO
!                                                                        
      END SUBROUTINE MakeInterpMatFromTo
      !
      !
      !
      SUBROUTINE Interp3DArray( inDim, inArray, outDim, outArray, interpXi, interpEta, interpZeta )
      

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER          , DIMENSION(3)                              :: inDim
      INTEGER          , DIMENSION(3)                              :: outDim
      REAL(KIND=RP), DIMENSION(inDim(1) ,inDim(2) , inDim(3) ) :: inArray
      REAL(KIND=RP), DIMENSION(outDim(1),outDim(2), outDim(3)) :: outArray
      REAL(KIND=RP), DIMENSION(outDim(1), inDim(1))            :: interpXi
      REAL(KIND=RP), DIMENSION(outDim(2), inDim(2))            :: interpEta
      REAL(KIND=RP), DIMENSION(outDim(3), inDim(3))            :: interpZeta
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(:)    , POINTER :: tempIn,tempOut
      REAL(KIND=RP), DIMENSION(:,:,:), POINTER :: tempArray
      
      INTEGER          , DIMENSION(3)              :: maxDim
      INTEGER                                      :: i,j,k
!
!     -------------------------------------------------
!     Allocate temporary storage for the interpolations
!     -------------------------------------------------
!
      ALLOCATE( tempIn ( MAXVAL(inDim)  ) )
      ALLOCATE( tempOut( MAXVAL(outDim) ) )
      maxDim = MAX(inDim,OutDim)
      ALLOCATE( tempArray(maxDim(1), maxDim(2), maxDim(3)) )
!
!     -----------------
!     Interpolate in Xi
!     -----------------
!
      DO k = 1, inDim(3)
         DO j = 1, inDim(2)
            DO i = 1, inDim(1)
               tempIn(i) = inArray(i,j,k)
            END DO 
            CALL InterpVectToVect( inDim(1), tempIn, outDim(1), tempOut, interpXi )
            DO i = 1, outDim(1)
               tempArray(i,j,k) = tempOut(i)
            END DO 
         END DO
      END DO 
!
!     ------------------
!     Interpolate in Eta
!     ------------------
!
      DO k = 1, inDim(3)
         DO i = 1, outDim(1)
            DO j = 1, inDim(2)
               tempIn(j) = tempArray(i,j,k)
            END DO 
            CALL InterpVectToVect( inDim(2), tempIn , outDim(2), tempOut, interpEta )
            DO j = 1, outDim(2)
               tempArray(i,j,k) = tempOut(j)
            END DO 
         END DO 
      END DO 
!
!     -------------------
!     Interpolate in Zeta
!     -------------------
!
      DO j = 1, outDim(2)
         DO i = 1, outDim(1)
            DO k = 1, inDim(3)
               tempIn(k) = tempArray(i,j,k)
            END DO 
            CALL InterpVectToVect( inDim(3), tempIn, outDim(3), tempOut, interpZeta )
            DO k = 1, outDim(3)
               outArray(i,j,k) = tempOut(k)
            END DO 
         END DO 
      END DO 
!
!     ---------------
!     Clean up memory
!     ---------------
!
      DEALLOCATE (tempIn)
      DEALLOCATE (tempOut)
      DEALLOCATE (tempArray)
      
      END SUBROUTINE Interp3DArray
      
!
!---------------------------------------------------------------------
!!    Compute at x the lagrange polynomial L_k of degree n-1    
!!    whose zeros are given by the z(i) 
!---------------------------------------------------------------------
!
      FUNCTION EvaluateLagrangePoly(k,x,n,z)
!
      
      IMPLICIT NONE
      REAL(KIND=RP) :: EvaluateLagrangePoly
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                         :: k !! Which poly
      INTEGER, INTENT(IN)                         :: n !! Order+1
      REAL(KIND=RP), DIMENSION(n), INTENT(IN) :: z !! Nodes
      REAL(KIND=RP)              , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: j
      REAL(KIND=RP) :: eLP
      
!                                                                       
      IF(k == 1)     THEN 
         eLP = (x - z(2))/(z(k) - z(2)) 
         DO  j = 3,n 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
      ELSE 
         eLP = (x - z(1))/(z(k) - z(1)) 
         DO  j = 2,k-1 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
         DO j = k+1,n 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
      END IF
      EvaluateLagrangePoly= eLp
!                                                                       
      END FUNCTION EvaluateLagrangePoly
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute at x the derivative of the lagrange polynomial L_k of
!!    degree n-1 whose zeros are given by the z(i)
!---------------------------------------------------------------------
!
      FUNCTION EvaluateLagrangePolyDeriv(k,x,n,z)
!
      
      IMPLICIT NONE
      REAL(KIND=RP) :: EvaluateLagrangePolyDeriv
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                         :: k !! Which poly
      INTEGER, INTENT(IN)                         :: n !! Order+1
      REAL(KIND=RP), DIMENSION(n), INTENT(IN) :: z !! Nodes
      REAL(KIND=RP)              , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: l,m
      REAL(KIND=RP) :: hp,poly
!                                                                       
      hp = 0.0_RP
      DO l = 1,n 
         if(l == k)     CYCLE
         poly = 1.0_RP
         DO m = 1,n 
            if(m == l)     CYCLE
            if(m == k)     CYCLE 
            poly = poly*(x - z(m))/(z(k) - z(m))
         END DO
         hp = hp + poly/(z(k) - z(l)) 
      END DO
      EvaluateLagrangePolyDeriv = hp 
!                                                                       
      END FUNCTION EvaluateLagrangePolyDeriv
!                                                                       
!///////////////////////////////////////////////////////////////////////
!
!                                                                       
!///////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Evaluate an interpolant at an array of points by matrix
!!    multiplication
!---------------------------------------------------------------------
!
      SUBROUTINE InterpVectToVect(nOld,oldVals,nNew,results,interpMat)
!
!     -------------------
!     Date: June 28, 2002
!     -------------------
!
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                           , INTENT(IN)  :: nOld,nNew
      REAL(KIND=RP), DIMENSION(nOld), INTENT(IN)  :: OldVals
      REAL(KIND=RP), DIMENSION(nNew), INTENT(OUT) :: results
      REAL(KIND=RP), DIMENSION(nNew,nOld)         :: interpMat
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER           :: i,k
      REAL(KIND=RP) :: sum
!                                                                       
      DO i = 1,nNew 
         sum = 0.0_RP 
         DO k = 1,nOld
            sum = sum + interpMat(i,k)*oldVals(k) 
         END DO
         results(i) = sum
      END DO
!                                                                        
      END SUBROUTINE InterpVectToVect
!                                                                       
!///////////////////////////////////////////////////////////////////////
!
      
END Module MappedGeometryClass
