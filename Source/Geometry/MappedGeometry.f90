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
      CLASS(MappedGeometry)     :: self
      TYPE(TransfiniteHexMap)  :: mapper
      TYPE(NodalStorage)       :: spA
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: N
      INTEGER       :: i, j, k
      REAL(KIND=RP) :: x(3), nrm
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
      CALL computeMetricTerms(self, spA, mapper)
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
            self % scal(i,j,ETOP)     = nrm
           
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
      SUBROUTINE computeMetricTerms(self, spA, mapper)  
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

      END SUBROUTINE computeMetricTerms
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

END Module MappedGeometryClass
