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
   USE Nodal2DStorageClass
   IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2, TOP = 2, BOTTOM = 1
      INTEGER, PARAMETER :: ELEFT = 4, ERIGHT = 2, ETOP = 3, EBOTTOM = 1
!
!     -----
!     Class
!     -----
!
      TYPE MappedGeometry
            INTEGER                                       :: N, M
            REAL(KIND=RP), DIMENSION(:,:)   , ALLOCATABLE :: X_xi, X_eta, Y_xi, Y_eta
            REAL(KIND=RP), DIMENSION(:,:)   , ALLOCATABLE :: x, y
            REAL(KIND=RP), DIMENSION(:,:,:) , ALLOCATABLE :: xb
            REAL(KIND=RP), DIMENSION(:,:)   , ALLOCATABLE :: jacobian, scal
            REAL(KIND=RP), DIMENSION(:,:,:) , ALLOCATABLE :: normal
      END TYPE MappedGeometry
   
      INTERFACE Construct
         MODULE PROCEDURE ConstructMappedGeometry
      END INTERFACE Construct

!
!  ========
   CONTAINS 
!  ========
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ConstructMappedGeometry( this, spA, N, M, mapper )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
      INTEGER, INTENT(IN)       :: N, M
      TYPE(TransfiniteQuadMap)  :: mapper
      TYPE(MappedGeometry)      :: this
      TYPE(Nodal2DStorage)      :: spA
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: i, j
      REAL(KIND=RP) :: x(2)
      REAL(KIND=RP) :: metricMatrix(2,2), nrm, s, jac
!
!     -----------
!     Allocations
!     -----------
!
      ALLOCATE( this%X_xi(0:N,0:N), this%X_eta(0:N,0:N), this%Y_xi(0:N,0:N), this%Y_eta(0:N,0:N))
      ALLOCATE( this%x(0:N,0:N), this%y(0:N,0:N), this%jacobian(0:N,0:N) )
      ALLOCATE( this%xb(2,0:MAX(N,M),4) )
      ALLOCATE( this%normal(0:MAX(N,M),2,4), this%scal(0:MAX(N,M),4) )
      
      this%N = N
      this%M = M
!
!     ----------------
!     Set the geometry
!     ----------------
!
      DO j= 0,M       
         DO i = 0,N 
            CALL EvaluateTransfiniteMapAt( mapper, spA%xi(i), spA%eta(j), x )
            this%x(i,j) = x(1)
            this%y(i,j) = x(2)
            CALL EvaluateMetricDerivatives( mapper, spA%xi(i), spA%eta(j), metricMatrix )
            this%X_xi(i,j)  = metricMatrix(1,1)
            this%X_eta(i,j) = metricMatrix(1,2)
            this%Y_xi(i,j)  = metricMatrix(2,1)
            this%Y_eta(i,j) = metricMatrix(2,2)
            this%jacobian(i,j) = metricMatrix(1,1)*metricMatrix(2,2) - &
                                 metricMatrix(1,2)*metricMatrix(2,1)
         END DO
      END DO
      DO i = 0, N 
         CALL EvaluateTransfiniteMapAt( mapper, spA%xi(i), -1.0_RP, x )
         this%xb(:,i,1) = x
         CALL EvaluateTransfiniteMapAt( mapper, spA%xi(i),  1.0_RP, x )
         this%xb(:,i,3) = x
      END DO
      DO j = 0, M 
         CALL EvaluateTransfiniteMapAt( mapper, -1.0_RP, spA%eta(j), x )
         this%xb(:,j,4) = x
         CALL EvaluateTransfiniteMapAt( mapper,  1.0_RP, spA%eta(j), x )
         this%xb(:,j,2) = x
      END DO
!
!     ----------------
!     Boundary Normals - Must be evaluated at the boundaries!
!     ----------------
!
      DO j = 0,M !\hat n^1
      i = 0
         CALL EvaluateMetricDerivatives( mapper, -1.0_RP, spA%eta(j), metricMatrix )
         jac = metricMatrix(1,1)*metricMatrix(2,2) - metricMatrix(1,2)*metricMatrix(2,1)
         nrm = SQRT( metricMatrix(2,2)**2 + metricMatrix(1,2)**2 )
         s   = -SIGN( 1.0_RP, jac ) !LEFT
         this%normal(j,1,4) =  s*metricMatrix(2,2)/nrm
         this%normal(j,2,4) = -s*metricMatrix(1,2)/nrm
         this%scal(j,4) = nrm
      i = N
         CALL EvaluateMetricDerivatives( mapper, 1.0_RP, spA%eta(j), metricMatrix )
         jac = metricMatrix(1,1)*metricMatrix(2,2) - metricMatrix(1,2)*metricMatrix(2,1)
         nrm = SQRT( metricMatrix(2,2)**2 + metricMatrix(1,2)**2 )
         s   = SIGN( 1.0_RP, jac ) !RIGHT
         this%normal(j,1,2) =  s*metricMatrix(2,2)/nrm
         this%normal(j,2,2) = -s*metricMatrix(1,2)/nrm
         this%scal(j,2) = nrm
      END DO
      
      DO i = 0,N !\hat n^2
      j = 0
         CALL EvaluateMetricDerivatives( mapper, spA%xi(i), -1.0_RP, metricMatrix )
         jac = metricMatrix(1,1)*metricMatrix(2,2) - metricMatrix(1,2)*metricMatrix(2,1)
         nrm = SQRT( metricMatrix(2,1)**2 + metricMatrix(1,1)**2 )
         s   = -SIGN( 1.0_RP, jac ) !DOWN
         this%normal(i,1,1) = -s*metricMatrix(2,1)/nrm
         this%normal(i,2,1) =  s*metricMatrix(1,1)/nrm
         this%scal(i,1) = nrm
      j = M
         CALL EvaluateMetricDerivatives( mapper, spA%xi(i), 1.0_RP, metricMatrix )
         jac = metricMatrix(1,1)*metricMatrix(2,2) - metricMatrix(1,2)*metricMatrix(2,1)
         nrm = SQRT( metricMatrix(2,1)**2 + metricMatrix(1,1)**2 )
         s   = SIGN( 1.0_RP, jac ) !UP
         this%normal(i,1,3) = -s*metricMatrix(2,1)/nrm
         this%normal(i,2,3) =  s*metricMatrix(1,1)/nrm
         this%scal(i,3) = nrm
      END DO
   END SUBROUTINE ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE DestructMappedGeometry(this)
      IMPLICIT NONE 
      TYPE(MappedGeometry) :: this
      DEALLOCATE( this%X_xi, this%X_eta, this%Y_xi, this%Y_eta, this%x, this%y, this%jacobian )
      DEALLOCATE( this%xb )
      DEALLOCATE( this%normal, this%scal )
   END SUBROUTINE DestructMappedGeometry
   
END Module MappedGeometryClass
