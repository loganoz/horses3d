!
!////////////////////////////////////////////////////////////////////////
!
!      TransfiniteMapsTests.f90
!      Created: May 20, 2015 at 1:57 PM 
!      By: David Kopriva  
!
!      Define tests for the TransfiniteMapClass transformation procedures
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE testCube
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP) :: corners(3,8)
         REAL(KIND=RP) :: du, dv, dw, u, v, w
         REAL(KIND=RP) :: x(3), p(3), e, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER       :: i, j, k, N, M, L
!
!        --------------------------------------------
!        Set up a cube with edges with length [1,2,3]
!        --------------------------------------------
!
         corners(:,1) = [0.0_RP,0.0_RP,0.0_RP]
         corners(:,2) = [1.0_RP,0.0_RP,0.0_RP]
         corners(:,3) = [1.0_RP,2.0_RP,0.0_RP]
         corners(:,4) = [0.0_RP,2.0_RP,0.0_RP]
         corners(:,5) = [0.0_RP,0.0_RP,3.0_RP]
         corners(:,6) = [1.0_RP,0.0_RP,3.0_RP]
         corners(:,7) = [1.0_RP,2.0_RP,3.0_RP]
         corners(:,8) = [0.0_RP,2.0_RP,3.0_RP]
         
         N = 4
         M = 4
         L = 4
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*i
                  p    = [u,v,w]
                  p    = 0.5*(p + 1.0_RP)
                  p(2) = 2*p(2)
                  p(3) = 3*p(3)
                  CALL Hex8TransfiniteMap([u,v,w],x,corners)
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "Transfinite map evaluation failed")
!
!        --------------
!        Gradients Test
!        --------------
!
         grad_x_Exact(1,:) = [0.5_RP, 0.0_RP, 0.0_RP]
         grad_x_Exact(2,:) = [0.0_RP, 1.0_RP, 0.0_RP]
         grad_x_Exact(3,:) = [0.0_RP, 0.0_RP, 1.5_RP]
        
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u = -1.0_RP + du*i
                  CALL GradHex8TransfiniteMap([u,v,w],grad_x,corners)
                  e = MAX(e,MAXVAL(ABS(grad_x - grad_x_Exact)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "GradHex8TransfiniteMap on cube failure")
         
      END SUBROUTINE testCube
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE testGenHexAsCube
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP)   , DIMENSION(3,8)   :: corners
         TYPE(FacePatch) , DIMENSION(6)     :: faceData
         REAL(KIND=RP)   , DIMENSION(2)     :: uKnots = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)   , DIMENSION(2)     :: vKnots = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)   , DIMENSION(3,2,2) :: points
         REAL(KIND=RP)                      :: du, dv, dw, u, v, w
         REAL(KIND=RP)                      :: x(3), p(3), e, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER                            :: i, j, k, N, M, L
         INTEGER                            :: iMax, jMax, kMax
         INTEGER                            :: iUnit
!
!        -----------------
!        Construct corners
!        -----------------
!
         corners(:,1) = [0.0_RP,0.0_RP,0.0_RP]
         corners(:,2) = [1.0_RP,0.0_RP,0.0_RP]
         corners(:,3) = [1.0_RP,2.0_RP,0.0_RP]
         corners(:,4) = [0.0_RP,2.0_RP,0.0_RP]
         corners(:,5) = [0.0_RP,0.0_RP,3.0_RP]
         corners(:,6) = [1.0_RP,0.0_RP,3.0_RP]
         corners(:,7) = [1.0_RP,2.0_RP,3.0_RP]
         corners(:,8) = [0.0_RP,2.0_RP,3.0_RP]
!
!        ---------------
!        Construct faces
!        ---------------
!
         points(:,1,1) = corners(:,1)
         points(:,2,1) = corners(:,2)
         points(:,2,2) = corners(:,6)
         points(:,1,2) = corners(:,5)
         
         CALL ConstructFacePatch(self = faceData(1),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = corners(:,4)
         points(:,2,1) = corners(:,3)
         points(:,2,2) = corners(:,7)
         points(:,1,2) = corners(:,8)
         
         CALL ConstructFacePatch(self = faceData(2),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = corners(:,1)
         points(:,2,1) = corners(:,2)
         points(:,2,2) = corners(:,3)
         points(:,1,2) = corners(:,4)
         
         CALL ConstructFacePatch(self = faceData(3),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = corners(:,5)
         points(:,2,1) = corners(:,6)
         points(:,2,2) = corners(:,7)
         points(:,1,2) = corners(:,8)
         
         CALL ConstructFacePatch(self = faceData(5),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = corners(:,1)
         points(:,2,1) = corners(:,4)
         points(:,2,2) = corners(:,8)
         points(:,1,2) = corners(:,5)
         
         CALL ConstructFacePatch(self = faceData(6),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = corners(:,2)
         points(:,2,1) = corners(:,3)
         points(:,2,2) = corners(:,7)
         points(:,1,2) = corners(:,6)
         
         CALL ConstructFacePatch(self = faceData(4),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         N = 8
         M = 8
         L = 8
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*i
                  p    = [u,v,w]
                  p    = 0.5*(p + 1.0_RP)
                  p(2) = 2*p(2)
                  p(3) = 3*p(3)
                  CALL GeneralHexTransfiniteMap([u,v,w],x,corners,faceData)
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,        &
                            tol           = 100*EPSILON(1.0_RP),           &
                            msg = "Test GeneralHexTransfiniteMap on cube")
!
!        --------------
!        Gradients Test
!        --------------
!
         grad_x_Exact(1,:) = [0.5_RP, 0.0_RP, 0.0_RP]
         grad_x_Exact(2,:) = [0.0_RP, 1.0_RP, 0.0_RP]
         grad_x_Exact(3,:) = [0.0_RP, 0.0_RP, 1.5_RP]
         
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u = -1.0_RP + du*i
                  CALL GradGeneralHexTransfiniteMap([u,v,w], grad_x, corners, faceData)
                  e = MAX(e,MAXVAL(ABS(grad_x - grad_x_Exact)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "Test GradGeneralHexTransfiniteMap on cube ")
         
!
!        -----------
!        Open a file
!        -----------
!
         iUnit = UnusedUnit()
         OPEN( UNIT = iUnit, FILE = "GeneralHexTransfiniteMap.tec" )
!
!        -------------
!        Set up header
!        -------------
!
         iMax = 8
         jMax = 8
         kMax = 8
         du = 2.0_RP/iMax
         dv = 2.0_RP/jMax
         dw = 2.0_RP/kMax
         
         WRITE(iUnit,*) 'VARIABLES = "X", "Y", "Z"'
         WRITE(iUnit,*) 'ZONE F=POINT, I=',iMax+1,', J=',jMax+1,', k=',kMax+1
         DO k = 0, kMax
            w = -1.0_RP + dw*k
            DO j = 0, jMax
               v = -1.0_RP + dv*j
               DO i = 0, iMax
                  u    = -1.0_RP + du*i
                  p    = [u,v,w]
                  p    = 0.5*(p + 1.0_RP)
                  p(2) = 2*p(2)
                  p(3) = 3*p(3)
                  CALL GeneralHexTransfiniteMap([u,v,w],x,corners,faceData)
                  WRITE(iUnit,*) x
               END DO   
            END DO
         END DO  
      END SUBROUTINE testGenHexAsCube
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE testFlaredHex8
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP) :: corners(3,8)
         REAL(KIND=RP) :: du, dv, dw, u, v, w
         REAL(KIND=RP) :: x(3), p(3), e, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER       :: i, j, k, N, M, L
         
         
         corners(:,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         corners(:,2) = [ 1.0_RP,-1.0_RP,0.0_RP]
         corners(:,3) = [ 1.0_RP, 1.0_RP,0.0_RP]
         corners(:,4) = [-1.0_RP, 1.0_RP,0.0_RP]
         corners(:,5) = [-2.0_RP,-2.0_RP,3.0_RP]
         corners(:,6) = [ 2.0_RP,-2.0_RP,3.0_RP]
         corners(:,7) = [ 2.0_RP, 2.0_RP,3.0_RP]
         corners(:,8) = [-2.0_RP, 2.0_RP,3.0_RP]
         
         N = 4
         M = 4
         L = 4
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*i
                  
                  p(1) = 0.5_RP*(1.0_RP-w)*u   + 0.5_RP*(1.0_RP+w)*2.0_RP*u
                  p(2) = 0.5_RP*(1.0_RP-w)*v   + 0.5_RP*(1.0_RP+w)*2.0_RP*v
                  p(3) = 0.5_RP*(1.0_RP-w)*0.0 + 0.5_RP*(1.0_RP+w)*3.0_RP

                  CALL Hex8TransfiniteMap([u,v,w],x,corners)
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "Transfinite map evaluation failed")
!
!        --------------
!        Gradients Test
!        --------------
!
         grad_x_Exact(1,:) = [1.0_RP, 0.0_RP, 0.0_RP]
         grad_x_Exact(2,:) = [0.0_RP, 2.0_RP, 0.0_RP]
         grad_x_Exact(3,:) = [0.0_RP, 0.0_RP, 3.0_RP]
        
         
         
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u = -1.0_RP + du*j
                  
                  grad_x_Exact(1,1) =  0.5_RP*(1.0_RP-w)   + 0.5_RP*(1.0_RP+w)*2.0_RP
                  grad_x_Exact(1,2) =  0.0_RP
                  grad_x_Exact(1,3) = -0.5_RP*u + 0.5_RP*2.0_RP*u
                  
                  grad_x_Exact(2,1) =  0.0_RP
                  grad_x_Exact(2,2) =  0.5_RP*(1.0_RP-w)   + 0.5_RP*(1.0_RP+w)*2.0_RP
                  grad_x_Exact(2,3) = -0.5_RP*v + 0.5_RP*2.0_RP*v
                  
                  grad_x_Exact(3,1) =  0.0_RP
                  grad_x_Exact(3,2) =  0.0_RP
                  grad_x_Exact(3,3) =  0.5_RP*3.0_RP
                  
                  CALL GradHex8TransfiniteMap([u,v,w],grad_x,corners)
                  e = MAX(e,MAXVAL(ABS(grad_x - grad_x_Exact)))
               END DO
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "GradHex8TransfiniteMap failure")
         
      END SUBROUTINE testFlaredHex8
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE testFlaredHex8AsGeneral
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP)   , DIMENSION(3,8)   :: corners
         TYPE(FacePatch) , DIMENSION(6)     :: faceData
         REAL(KIND=RP)   , DIMENSION(2)     :: uKnots = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)   , DIMENSION(2)     :: vKnots = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)   , DIMENSION(3,2,2) :: points
         REAL(KIND=RP)                      :: du, dv, dw, u, v, w
         REAL(KIND=RP)                      :: x(3), p(3), e, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER                            :: i, j, k, N, M, L
         INTEGER                            :: iUnit
         INTEGER                            :: iMax, jMax, kMax
!
!        -----------------
!        Construct corners
!        -----------------
!
         corners(:,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         corners(:,2) = [ 1.0_RP,-1.0_RP,0.0_RP]
         corners(:,3) = [ 1.0_RP, 1.0_RP,0.0_RP]
         corners(:,4) = [-1.0_RP, 1.0_RP,0.0_RP]
         corners(:,5) = [-2.0_RP,-2.0_RP,3.0_RP]
         corners(:,6) = [ 2.0_RP,-2.0_RP,3.0_RP]
         corners(:,7) = [ 2.0_RP, 2.0_RP,3.0_RP]
         corners(:,8) = [-2.0_RP, 2.0_RP,3.0_RP]
!
!        ---------------
!        Construct faces
!        ---------------
!
         points(:,1,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,2) = [ 2.0_RP,-2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP,-2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(1),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,2) = [ 2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP, 2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(2),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,2) = [ 1.0_RP, 1.0_RP,0.0_RP]
         points(:,1,2) = [-1.0_RP, 1.0_RP,0.0_RP]
         
         CALL ConstructFacePatch(self = faceData(3),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-2.0_RP,-2.0_RP,3.0_RP]
         points(:,2,1) = [ 2.0_RP,-2.0_RP,3.0_RP]
         points(:,2,2) = [ 2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP, 2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(5),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [-1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,2) = [-2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP,-2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(6),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [ 1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,2) = [ 2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [ 2.0_RP,-2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(4),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         N = 4
         M = 4
         L = 4
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*i
                  
                  p(1) = 0.5_RP*(1.0_RP-w)*u   + 0.5_RP*(1.0_RP+w)*2.0_RP*u
                  p(2) = 0.5_RP*(1.0_RP-w)*v   + 0.5_RP*(1.0_RP+w)*2.0_RP*v
                  p(3) = 0.5_RP*(1.0_RP-w)*0.0 + 0.5_RP*(1.0_RP+w)*3.0_RP

                  CALL GeneralHexTransfiniteMap([u,v,w],x,corners,faceData)
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,        &
                            tol           = 100*EPSILON(1.0_RP),           &
                            msg = "GeneralHexTransfiniteMap on flared hex8")
!
!        --------------
!        Gradients Test
!        --------------
!
         grad_x_Exact(1,:) = [1.0_RP, 0.0_RP, 0.0_RP]
         grad_x_Exact(2,:) = [0.0_RP, 2.0_RP, 0.0_RP]
         grad_x_Exact(3,:) = [0.0_RP, 0.0_RP, 3.0_RP]
         
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u = -1.0_RP + du*i
                  
                  grad_x_Exact(1,1) =  0.5_RP*(1.0_RP-w)   + 0.5_RP*(1.0_RP+w)*2.0_RP
                  grad_x_Exact(1,2) =  0.0_RP
                  grad_x_Exact(1,3) = -0.5_RP*u + 0.5_RP*2.0_RP*u
                  
                  grad_x_Exact(2,1) =  0.0_RP
                  grad_x_Exact(2,2) =  0.5_RP*(1.0_RP-w)   + 0.5_RP*(1.0_RP+w)*2.0_RP
                  grad_x_Exact(2,3) = -0.5_RP*v + 0.5_RP*2.0_RP*v
                  
                  grad_x_Exact(3,1) =  0.0_RP
                  grad_x_Exact(3,2) =  0.0_RP
                  grad_x_Exact(3,3) =  0.5_RP*3.0_RP
                  
                  CALL GradGeneralHexTransfiniteMap([u,v,w], grad_x, corners, faceData)
                  e = MAX(e,MAXVAL(ABS(grad_x - grad_x_Exact)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "GradGeneralHexTransfiniteMap on flared hex8 ")
!
!        -----------
!        Open a file
!        -----------
!
         iUnit = UnusedUnit()
         OPEN( UNIT = iUnit, FILE = "testFlaredHex8AsGeneral.tec" )
!
!        -------------
!        Set up header
!        -------------
!
         iMax = 8
         jMax = 8
         kMax = 8
         du = 2.0_RP/iMax
         dv = 2.0_RP/jMax
         dw = 2.0_RP/kMax
         
         WRITE(iUnit,*) 'VARIABLES = "X", "Y", "Z"'
         WRITE(iUnit,*) 'ZONE F=POINT, I=',iMax+1,', J=',jMax+1,', k=',kMax+1
         DO k = 0, kMax
            w = -1.0_RP + dw*k
            DO j = 0, jMax
               v = -1.0_RP + dv*j
               DO i = 0, iMax
                  u    = -1.0_RP + du*i
                  p    = [u,v,w]
                  p    = 0.5*(p + 1.0_RP)
                  p(2) = 2*p(2)
                  p(3) = 3*p(3)
                  CALL Hex8TransfiniteMap([u,v,w],x,corners)
                  WRITE(iUnit,*) x
               END DO   
            END DO
         END DO  
         
      END SUBROUTINE testFlaredHex8AsGeneral
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE testCubeWithHexMapper
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP)           :: corners(3,8)
         REAL(KIND=RP)           :: du, dv, dw, u, v, w
         REAL(KIND=RP)           :: x(3), p(3), e, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER                 :: i, j, k, N, M, L
         TYPE(TransfiniteHexMap) :: mapper
!
!        --------------------------------------------
!        Set up a cube with edges with length [1,2,3]
!        --------------------------------------------
!
         corners(:,1) = [0.0_RP,0.0_RP,0.0_RP]
         corners(:,2) = [1.0_RP,0.0_RP,0.0_RP]
         corners(:,3) = [1.0_RP,2.0_RP,0.0_RP]
         corners(:,4) = [0.0_RP,2.0_RP,0.0_RP]
         corners(:,5) = [0.0_RP,0.0_RP,3.0_RP]
         corners(:,6) = [1.0_RP,0.0_RP,3.0_RP]
         corners(:,7) = [1.0_RP,2.0_RP,3.0_RP]
         corners(:,8) = [0.0_RP,2.0_RP,3.0_RP]
         
         CALL mapper % constructWithCorners(corners)
         
         N = 4
         M = 4
         L = 4
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*j
                  p    = [u,v,w]
                  p    = 0.5*(p + 1.0_RP)
                  p(2) = 2*p(2)
                  p(3) = 3*p(3)
                  x = mapper % transfiniteMapAt([u,v,w])
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "Transfinite map evaluation failed")
!
!        --------------
!        Gradients Test
!        --------------
!
         grad_x_Exact(1,:) = [0.5_RP, 0.0_RP, 0.0_RP]
         grad_x_Exact(2,:) = [0.0_RP, 1.0_RP, 0.0_RP]
         grad_x_Exact(3,:) = [0.0_RP, 0.0_RP, 1.5_RP]
        
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u = -1.0_RP + du*i
                  grad_x = mapper % metricDerivativesAt([u,v,w])
                  e = MAX(e,MAXVAL(ABS(grad_x - grad_x_Exact)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "GradHex8TransfiniteMap on cube failure")
         
      END SUBROUTINE testCubeWithHexMapper
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE testFlaredHex8AsGeneralWithMapper
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         TYPE(FacePatch) , DIMENSION(6)     :: faceData
         REAL(KIND=RP)   , DIMENSION(2)     :: uKnots = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)   , DIMENSION(2)     :: vKnots = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)   , DIMENSION(3,2,2) :: points
         REAL(KIND=RP)                      :: du, dv, dw, u, v, w
         REAL(KIND=RP)                      :: x(3), p(3), e, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER                            :: i, j, k, N, M, L
         TYPE(TransfiniteHexMap)            :: mapper
!
!        ---------------
!        Construct faces
!        ---------------
!
         points(:,1,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,2) = [ 2.0_RP,-2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP,-2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(1),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,2) = [ 2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP, 2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(2),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,2) = [ 1.0_RP, 1.0_RP,0.0_RP]
         points(:,1,2) = [-1.0_RP, 1.0_RP,0.0_RP]
         
         CALL ConstructFacePatch(self = faceData(3),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-2.0_RP,-2.0_RP,3.0_RP]
         points(:,2,1) = [ 2.0_RP,-2.0_RP,3.0_RP]
         points(:,2,2) = [ 2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP, 2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(5),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [-1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [-1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,2) = [-2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [-2.0_RP,-2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(6),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         points(:,1,1) = [ 1.0_RP,-1.0_RP,0.0_RP]
         points(:,2,1) = [ 1.0_RP, 1.0_RP,0.0_RP]
         points(:,2,2) = [ 2.0_RP, 2.0_RP,3.0_RP]
         points(:,1,2) = [ 2.0_RP,-2.0_RP,3.0_RP]
         
         CALL ConstructFacePatch(self = faceData(4),points = points, &
                                                    uKnots = uKnots, &
                                                    vKnots = vKnots)
         CALL mapper % constructWithFaces(faceData)
         
         N = 4
         M = 4
         L = 4
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*i
                  
                  p(1) = 0.5_RP*(1.0_RP-w)*u   + 0.5_RP*(1.0_RP+w)*2.0_RP*u
                  p(2) = 0.5_RP*(1.0_RP-w)*v   + 0.5_RP*(1.0_RP+w)*2.0_RP*v
                  p(3) = 0.5_RP*(1.0_RP-w)*0.0 + 0.5_RP*(1.0_RP+w)*3.0_RP
                  
                  x = mapper % transfiniteMapAt([u,v,w])
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,        &
                            tol           = 100*EPSILON(1.0_RP),           &
                            msg = "GeneralHexTransfiniteMap on flared hex8")
!
!        --------------
!        Gradients Test
!        --------------
!
         grad_x_Exact(1,:) = [1.0_RP, 0.0_RP, 0.0_RP]
         grad_x_Exact(2,:) = [0.0_RP, 2.0_RP, 0.0_RP]
         grad_x_Exact(3,:) = [0.0_RP, 0.0_RP, 3.0_RP]
         
         e = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u = -1.0_RP + du*i
                  
                  grad_x_Exact(1,1) =  0.5_RP*(1.0_RP-w)   + 0.5_RP*(1.0_RP+w)*2.0_RP
                  grad_x_Exact(1,2) =  0.0_RP
                  grad_x_Exact(1,3) = -0.5_RP*u + 0.5_RP*2.0_RP*u
                  
                  grad_x_Exact(2,1) =  0.0_RP
                  grad_x_Exact(2,2) =  0.5_RP*(1.0_RP-w)   + 0.5_RP*(1.0_RP+w)*2.0_RP
                  grad_x_Exact(2,3) = -0.5_RP*v + 0.5_RP*2.0_RP*v
                  
                  grad_x_Exact(3,1) =  0.0_RP
                  grad_x_Exact(3,2) =  0.0_RP
                  grad_x_Exact(3,3) =  0.5_RP*3.0_RP
                  
                  grad_x = mapper % metricDerivativesAt([u,v,w])
                  e = MAX(e,MAXVAL(ABS(grad_x - grad_x_Exact)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,        &
                            tol = 100*EPSILON(1.0_RP),           &
                            msg = "GradGeneralHexTransfiniteMap on flared hex8 ")
         
      END SUBROUTINE testFlaredHex8AsGeneralWithMapper
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE testCylindricalMappingWithMapper  
         USE FTAssertions  
         USE SMConstants
         USE TransfiniteMapClass
         use FacePatchClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                            :: nUknots, nVKnots
         REAL(KIND=RP)   , ALLOCATABLE      :: uKnots(:)
         REAL(KIND=RP)   , ALLOCATABLE      :: vKnots(:)
         REAL(KIND=RP)   , ALLOCATABLE      :: points(:,:,:)
         TYPE(FacePatch) , DIMENSION(6)     :: faceData
         
         REAL(KIND=RP)                      :: du, dv, dw, u, v, w
         REAL(KIND=RP)                      :: x(3), p(3), e, e_grad, grad_x(3,3), grad_x_Exact(3,3)
         INTEGER                            :: i, j, k, N, M, L
         TYPE(TransfiniteHexMap)            :: mapper
         
         nUknots = 10
         nVKnots = 10
         ALLOCATE(uKnots(nUKnots))
         ALLOCATE(vKnots(nVKnots))
         ALLOCATE(points(3,nUKnots,nVKnots))
         
         DO i = 1, nUKnots
            uKnots(i) = -COS((i-1)*PI/(nUKnots-1)) 
         END DO  
         
         DO j = 1, nVKnots
            vKnots(j) = -COS((j-1)*PI/(nVKnots-1)) 
         END DO  
!
!        ------------------------------
!        Construct the six face patches
!        ------------------------------
!
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([uKnots(i), -1.0_RP, vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(1), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([uKnots(i), 1.0_RP, vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(2), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([uKnots(i), vKnots(j),-1.0_RP], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(3), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([uKnots(i), vKnots(j),1.0_RP], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(5), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([-1.0_RP, uKnots(i), vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(4), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([1.0_RP, uKnots(i), vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(6), uKnots, vKnots, points )
!
!        ---------------------
!        Construct the mapper 
!        ---------------------
!
         CALL mapper % constructWithFaces(faceData)
!
!        ------------------------------------------------------------------
!        Test the integrity of the points and derivatives.
!        They will only be good to within spectral accuracy.
!        ------------------------------------------------------------------
!
         N = 6
         M = 6
         L = 6
         du = 2.0_RP/N
         dv = 2.0_RP/M
         dw = 2.0_RP/L
!
!        ------------
!        Mapping test
!        ------------
!
         e      = 0.0_RP
         e_grad = 0.0_RP
         DO k = 0, L
            w = -1.0_RP + dw*k
            DO j = 0, M
               v = -1.0_RP + dv*j
               DO i = 0, N
                  u    = -1.0_RP + du*i
                  
                  x = mapper % transfiniteMapAt([u,v,w])
                  CALL cylindricalGeometry([u,v,w],p)
                  e = MAX(e,MAXVAL(ABS(p-x)))
                  
                  grad_x = mapper % metricDerivativesAt([u,v,w])
                  CALL cylindricalGeometryDerivatives([u,v,w],grad_x_Exact) 
                  e_grad = MAX(e_grad,MAXVAL(ABS(grad_x-grad_x_Exact)))
               END DO   
            END DO
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,      &
                            tol           = 1.0d-9, &
                            msg = "Test cylindrical geometry with mapper")
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = e,      &
                            tol           = 1.0d-9, &
                            msg = "Test cylindrical geometry derivatives with mapper")
                            
         CALL plotWithMapper(mapper, "CylindricalGeometry.tec")
         
      END SUBROUTINE testCylindricalMappingWithMapper
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cylindricalGeometry(u,x)
         USE SMConstants
         IMPLICIT NONE  
         REAL(KIND=RP) :: u(3) ! IN [-1,1]^3
         REAL(KIND=RP) :: x(3) ! physical space locations
         REAL(KIND=RP) :: r0 = 1.0_RP, rMax = 2.0_RP, theta0 = 0.0_RP, &
                          thetaMax = PI/2, z0 = 0.0_RP, zMax = 3.0_RP
         REAL(KIND=RP) :: r, theta, z
         
         r     = r0 + 0.5_RP*(u(1)+1.0_RP)*(rMax - r0)
         theta = theta0 + 0.5_RP*(u(2)+1.0_RP)*(thetaMax - theta0)
         z     = z0 + 0.5_RP*(u(3)+1.0_RP)*(zMax - z0)
         
         x(1) = r*COS(theta)
         x(2) = r*SIN(theta)
         x(3) = z
         
      END SUBROUTINE cylindricalGeometry
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cylindricalGeometryDerivatives(u,grad_x) 
         USE SMConstants 
         IMPLICIT NONE  
         REAL(KIND=RP) :: u(3), grad_x(3,3)
         
         REAL(KIND=RP) :: r0 = 1.0_RP, rMax = 2.0_RP, theta0 = 0.0_RP, &
                          thetaMax = PI/2, z0 = 0.0_RP, zMax = 3.0_RP
         REAL(KIND=RP) :: r, theta, z, drDxi, dthetaDEta, dZdZeta
         
         r     = r0 + 0.5_RP*(u(1)+1.0_RP)*(rMax - r0)
         theta = theta0 + 0.5_RP*(u(2)+1.0_RP)*(thetaMax - theta0)
         z     = z0 + 0.5_RP*(u(3)+1.0_RP)*(zMax - z0)
         
         drDxi      = 0.5_RP*(rMax - r0)
         dthetaDEta = 0.5_RP*(thetaMax - theta0)
         dZdZeta    = 0.5_RP*(zMax - z0)
         
         grad_x(:,1) = [drDxi*COS(theta), drDxi*SIN(theta),0.0_RP];
         grad_x(:,2) = [-r*SIN(theta)*dthetaDEta, r*COS(theta)*dthetaDEta,0.0_RP];
         grad_x(:,3) = [0.0_RP, 0.0_RP, dZdZeta]
         
      END SUBROUTINE cylindricalGeometryDerivatives
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE plotWithMapper(mapper, fName)
         USE SMConstants
         USE TransfiniteMapClass
         use Utilities, only:UnusedUnit
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(TransfiniteHexMap) :: mapper
         CHARACTER(LEN=*)        :: fName
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                 :: i, j, k
         INTEGER                 :: iUnit
         INTEGER                 :: iMax, jMax, kMax
         REAL(KIND=RP)           :: u, v, w, du, dv, dw, x(3)
!
!        -----------
!        Open a file
!        -----------
!
         iUnit = UnusedUnit()
         OPEN( UNIT = iUnit, FILE = fName )
!
!        -------------
!        Set up header
!        -------------
!
         iMax = 8
         jMax = 8
         kMax = 8
         du = 2.0_RP/iMax
         dv = 2.0_RP/jMax
         dw = 2.0_RP/kMax
         
         WRITE(iUnit,*) 'VARIABLES = "X", "Y", "Z"'
         WRITE(iUnit,*) 'ZONE F=POINT, I=',iMax+1,', J=',jMax+1,', k=',kMax+1
         
         DO k = 0, kMax
            w = -1.0_RP + dw*k
            DO j = 0, jMax
               v = -1.0_RP + dv*j
               DO i = 0, iMax
                  u    = -1.0_RP + du*i
                  x = mapper % transfiniteMapAt([u, v, w])
                  WRITE(iUnit,*) x
               END DO   
            END DO
         END DO  

      END SUBROUTINE plotWithMapper