!
!////////////////////////////////////////////////////////////////////////
!
!      FacePatchTests.f90
!      Created: May 20, 2015 at 12:34 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE quadraticFaceTest  
         IMPLICIT NONE  
         EXTERNAL :: quadraticFaceSurface, quadraticFaceSurfaceGradient
!
         CALL surfaceTest(3, 3, quadraticFaceSurface, quadraticFaceSurfaceGradient, 1.d-13)
         
      END SUBROUTINE quadraticFaceTest
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE flatFaceTest  
         IMPLICIT NONE  
         EXTERNAL :: flatFaceSurface     , flatFaceSurfaceGradient
!
         CALL surfaceTest(2, 2, flatFaceSurface, flatFaceSurfaceGradient, 1.d-13)
         
      END SUBROUTINE flatFaceTest
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE cubicFaceTest  
         IMPLICIT NONE  
         EXTERNAL :: cubicFaceSurface     , cubicFaceSurfaceGradient
!
         CALL surfaceTest(4, 4, cubicFaceSurface, cubicFaceSurfaceGradient, 1.d-13)
         
      END SUBROUTINE cubicFaceTest
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE sphericalFaceTest  
         IMPLICIT NONE  
         EXTERNAL :: sphericalFaceSurface, sphericalFaceSurfaceGradient
!
         CALL surfaceTest(10, 10, sphericalFaceSurface, sphericalFaceSurfaceGradient, 1.d-11)
         
      END SUBROUTINE sphericalFaceTest
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE surfaceTest(nUknots, nVKnots, faceSurfaceSubroutine, faceGradientSubroutine, tol)
         USE FTAssertions
         USE SMConstants
         USE FacePatchClass
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         INTEGER                                      :: nUknots, nVKnots
         REAL(KIND=RP)                                :: tol
         EXTERNAL                                     :: faceSurfaceSubroutine, faceGradientSubroutine
!
!        ---------------
!        Local variables
!        ---------------
!
         TYPE(FacePatch)                              :: face
         REAL(KIND=RP), DIMENSION(:)    , ALLOCATABLE :: uKnots, vKnots
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: points
         INTEGER                                      :: i, j
         INTEGER                                      :: N, M
         REAL(KIND=RP)                                :: eU, eV, eP, p(3), pInt(3)
         REAL(KIND=RP)                                :: du, dv, u, v
         REAL(KIND=RP)                                :: gradExact(3,2), grad(3,2)
         
!
!        --------------------------------
!        Construct a quadratic Face patch
!        --------------------------------
!
         ALLOCATE(uKnots(nUKnots))
         ALLOCATE(vKnots(nVKnots))
         ALLOCATE(points(3,nUKnots,nVKnots))
         
         DO i = 1, nUKnots
            uKnots(i) = -COS((i-1)*PI/(nUKnots-1)) 
         END DO  
         
         DO j = 1, nVKnots
            vKnots(j) = -COS((j-1)*PI/(nVKnots-1)) 
         END DO  
         
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL faceSurfaceSubroutine(uKnots(i), vKnots(j), points(:,i,j))
            END DO   
         END DO
         
         CALL ConstructFacePatch(self = face,uKnots = uKnots,vKnots = vKnots,points = points)
!
!        -----------------------------------
!        Make sure info was read in properly
!        -----------------------------------
!
         eU = MAXVAL(ABS(uKnots-face%uKnots))
         eV = MAXVAL(ABS(vKnots-face%vKnots))
         eP = MAXVAL(ABS(points-face%points))
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = eP,     &
                            tol           = tol,    &
                            msg = "Points not successfully stored")
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = eU,     &
                            tol           = tol,    &
                            msg = "uKnots not successfully stored")
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = eV,     &
                            tol           = tol,    &
                            msg = "vKnots not successfully stored")
!
!        -------------------------------
!        Test interpolation at the knots
!        -------------------------------
!
         eP = 0.0_RP
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL ComputeFacePoint(self = face,u = [uKnots(i),vKnots(j)],p = p)
               eP = MAXVAL(ABS(p - points(:,i,j)))
            END DO
         END DO
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = eP,     &
                            tol           = tol,    &
                            msg = "Test Interpolation at knots")
!
!        ----------------------------------
!        Test interpolation at other points
!        ----------------------------------
!
         N = 4*nUknots
         M = 4*nVKnots
         du = 2.0_RP/N
         dv = 2.0_RP/M
         eP = 0.0_RP
         DO j = 0, M
            v = j*dv - 1.0_RP
            DO i = 0, N
               u = i*du - 1.0_RP
               CALL faceSurfaceSubroutine(u,v,p)
               CALL ComputeFacePoint(self = face,u = [u,v],p = pInt)
               eP = MAXVAL(ABS(p-pInt))
            END DO   
         END DO  
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = eP,     &
                            tol           = tol,    &
                            msg = "Test interpolation at uniform points")
!
!        ------------------------------
!        Test derivatives along surface
!        ------------------------------
!
         DO j = 0, M
            v = j*dv - 1.0_RP
            DO i = 0, N
               u = i*du - 1.0_RP
               CALL faceGradientSubroutine(u,v,gradExact)
               CALL ComputeFaceDerivative(self = face,u = [u,v],grad = grad)
               eP = MAXVAL(ABS(gradExact-grad))
            END DO   
         END DO  
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue   = eP,     &
                            tol           = tol,    &
                            msg = "Test interpolation derivatives at uniform points")
!
!        ----
!        Done
!        ----
!
         CALL DestructFacePatch(face)
      END SUBROUTINE surfaceTest
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE quadraticFaceSurface(x,y,p)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: x, y
         REAL(KIND=RP) :: p(3)
         
         p(1) = x
         p(2) = y
         p(3) = x**2 + y**2
         
      END SUBROUTINE quadraticFaceSurface
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE quadraticFaceSurfaceGradient(x,y,grad)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: x, y
         REAL(KIND=RP) :: grad(3,2)
!
!        -------------
!        u derivatives
!        -------------
!
         grad(1,1) = 1.0_RP
         grad(2,1) = 0.0_RP
         grad(3,1) = 2*x
!
!        ----------------------
!        v derivatives         
!        ----------------------
!
         grad(1,2) = 0.0_RP
         grad(2,2) = 1.0_RP
         grad(3,2) = 2*y

      END SUBROUTINE quadraticFaceSurfaceGradient
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cubicFaceSurface(x,y,p)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: x, y
         REAL(KIND=RP) :: p(3)
         
         p(1) = x
         p(2) = y
         p(3) = x**3 + y**3
         
      END SUBROUTINE cubicFaceSurface
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cubicFaceSurfaceGradient(x,y,grad)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: x, y
         REAL(KIND=RP) :: grad(3,2)
!
!        -------------
!        u derivatives
!        -------------
!
         grad(1,1) = 1.0_RP
         grad(2,1) = 0.0_RP
         grad(3,1) = 3*x**2
!
!        ----------------------
!        v derivatives         
!        ----------------------
!
         grad(1,2) = 0.0_RP
         grad(2,2) = 1.0_RP
         grad(3,2) = 3*y**2

      END SUBROUTINE cubicFaceSurfaceGradient
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE flatFaceSurface(u,v,p)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: u, v
         REAL(KIND=RP) :: p(3)
         
         p(1) = u
         p(2) = v
         p(3) = u+v
         
      END SUBROUTINE flatFaceSurface
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE flatFaceSurfaceGradient(u,v,grad)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: u, v
         REAL(KIND=RP) :: grad(3,2)
!
!        -------------
!        u derivatives
!        -------------
!
         grad(1,1) = 1.0_RP
         grad(2,1) = 0.0_RP
         grad(3,1) = 1
!
!        ----------------------
!        v derivatives         
!        ----------------------
!
         grad(1,2) = 0.0_RP
         grad(2,2) = 1.0_RP
         grad(3,2) = 1

      END SUBROUTINE flatFaceSurfaceGradient

!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE sphericalFaceSurface(u,v,p)  
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: u, v
         REAL(KIND=RP) :: p(3)
         REAL(KIND=RP) :: theta, phi
         
         theta = 0.5_RP*(u+1.0_RP)*PI/4.0_RP
         phi   = 0.5_RP*(v+1.0_RP)*PI/4.0_RP
         
         p(1) = COS(theta)*SIN(phi)
         p(2) = SIN(theta)*SIN(phi)
         p(3) = COS(phi)
         
      END SUBROUTINE sphericalFaceSurface
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE sphericalFaceSurfaceGradient(u,v,grad)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP) :: u,v
         REAL(KIND=RP) :: grad(3,2)
         REAL(KIND=RP) :: theta, phi, p8
         
         theta = 0.5_RP*(u+1.0_RP)*PI/4.0_RP
         phi   = 0.5_RP*(v+1.0_RP)*PI/4.0_RP
         p8    = PI/8.0_RP
!
!        -------------
!        u derivatives
!        -------------
!
         grad(1,1) = -p8*SIN(theta)*SIN(phi)
         grad(2,1) =  p8*COS(theta)*SIN(phi)
         grad(3,1) =  0.0_RP
!
!        ----------------------
!        v derivatives         
!        ----------------------
!
         grad(1,2) =  p8*COS(theta)*COS(phi)
         grad(2,2) =  p8*SIN(theta)*COS(phi)
         grad(3,2) = -p8*SIN(phi)

      END SUBROUTINE sphericalFaceSurfaceGradient