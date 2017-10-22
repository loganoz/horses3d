!
!////////////////////////////////////////////////////////////////////////
!
!      MappedGeometryTests.f90
!      Created: May 22, 2015 at 12:46 PM 
!      By: David Kopriva  
!
!      Define tests for the MappedGeometry class
!
!////////////////////////////////////////////////////////////////////////
! 
      SUBROUTINE cubeTest
         USE NodalStorageClass
         USE MappedGeometryClass
         USE SMConstants
         USE FTAssertions
         
         IMPLICIT NONE  
         
         TYPE(NodalStorage)      :: spA
         TYPE(MappedGeometry)    :: geom
         TYPE(TransfiniteHexMap) :: mapper
         REAL(KIND=RP)           :: corners(3,8)
         REAL(KIND=RP)           :: e, x(3), xMap(3)
         INTEGER                 :: i,j,k, N
         
         N = 5
         CALL spA % construct(GAUSS,N,N,N)
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
!
!        ---------------------
!        Generate the geometry
!        ---------------------
!
         CALL geom % construct(spA, mapper)
!
!        -----------------
!        Test the geometry
!        -----------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  xMap =  mapper % transfiniteMapAt([spA % xi(i), spA % eta(j), spA % zeta(k)])
                  x    = geom % x(:,i,j,k)
                  e    = MAX(MAXVAL(ABS(x - xMap)),e)
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of nodal locations")
!
!        -----------
!        Face values
!        -----------
!
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([-1.0_RP, spA % eta(i), spA % zeta(j)])
               x    = geom % xb(:,i,j,ELEFT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Left face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([1.0_RP, spA % eta(i), spA % zeta(j)])
               x    = geom % xb(:,i,j,ERIGHT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Right face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j), -1.0_RP ])
               x    = geom % xb(:,i,j,EBOTTOM)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Bottom face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j), 1.0_RP ])
               x    = geom % xb(:,i,j,ETOP)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Top face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i), -1.0_RP, spA % zeta(j)])
               x    = geom % xb(:,i,j,EFRONT)
               e    = MAXVAL(ABS(x - xMap))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Front face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i),  1.0_RP, spA % zeta(j)])
               x    = geom % xb(:,i,j,EBACK)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Back face locations")
!
!        -----------------
!        Check the normals
!        -----------------
!
         e     = 0.0_RP
         xMap = [-1.0_RP,0.0_RP,0.0_RP]
         DO j = 0, N
            DO i = 0, N
               x    = geom % normal(:,i,j,ELEFT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Normal at left face")
         e = 0.0_RP
         xMap = [1.0_RP,0.0_RP,0.0_RP]
         DO j = 0, N
            DO i = 0, N
               x    = geom % normal(:,i,j,ERIGHT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Normal at right face")
         e = 0.0_RP
         xMap = [0.0_RP,0.0_RP,-1.0_RP]
         DO j = 0, N
            DO i = 0, N
               x    = geom % normal(:,i,j,EBOTTOM)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Normal at bottom face")
         e = 0.0_RP
         xMap = [0.0_RP,0.0_RP,1.0_RP]
         DO j = 0, N
            DO i = 0, N
               x    = geom % normal(:,i,j,ETOP)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Normal at top face")
         e = 0.0_RP
         xMap = [0.0_RP,-1.0_RP,0.0_RP]
         DO j = 0, N
            DO i = 0, N
               x    = geom % normal(:,i,j,EFRONT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Normal at front face")
         e = 0.0_RP
         xMap = [0.0_RP,1.0_RP,0.0_RP]
         DO j = 0, N
            DO i = 0, N
               x    = geom % normal(:,i,j,EBACK)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Normal at back face")
!
!        ----------------
!        Check the volume
!        ----------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
               e = MAX(e, ABS(geom % jacobian(i,j,k) - 0.75_RP))
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of Jacobian")
!
!        -----------------------
!        Metric terms -- jGradXi
!        -----------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  x    = geom % jGradXi(:,i,j,k)
                  e    = MAX(MAXVAL(ABS(x - [1.5_RP,0.0_RP,0.0_RP])),e)
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of JGradXi")
!
!        ------------------------
!        Metric terms -- jGradEta
!        ------------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  x    = geom % jGradEta(:,i,j,k)
                  e    = MAX(MAXVAL(ABS(x - [0.0_RP,0.75_RP,0.0_RP])),e)
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of JGradEta")
!
!        ------------------------
!        Metric terms -- Jacobian
!        ------------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  x    = geom % jGradZeta(:,i,j,k)
                  e    = MAX(MAXVAL(ABS(x - [0.0_RP,0.0_RP,0.5_RP])),e)
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of jGradZeta")
      END SUBROUTINE cubeTest
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cylinderTestGeometry  
         USE NodalStorageClass
         USE MappedGeometryClass
         USE TransfiniteMapClass
         USE SMConstants
         USE FTAssertions
         
         IMPLICIT NONE  
         
         TYPE(NodalStorage)      :: spA
         TYPE(MappedGeometry)    :: geom
         TYPE(TransfiniteHexMap) :: mapper
         REAL(KIND=RP)           :: corners(3,8), Jai(3,3), xMap(3)
         REAL(KIND=RP)           :: e, x(3), eScal, p(3), jac, exactJac
         REAL(KIND=RP)           :: u, v, w
         INTEGER                 :: i, j, k, N
         REAL(KIND=RP), EXTERNAL :: cylindricalJacobian
         EXTERNAL                :: cylindricalGeometry
         
         N = 10
         CALL spA % construct(GAUSS,N,N,N)
!
!        --------------------
!        Construct the mapper
!        --------------------
!
         CALL constructHexMapper(N, mapper, cylindricalGeometry)
!
!        ---------------------
!        Generate the geometry
!        ---------------------
!
         CALL geom % construct(spA, mapper)
!
!        ------------
!        Mapping test
!        ------------
!
         e = 0.0_RP
         DO k = 0, N
            w = spA % zeta(k)
            DO j = 0, N
               v = spA % eta(j)
               DO i = 0, N
                  u = spA % xi(i)
                  
                  CALL cylindricalGeometry([u,v,w],p)
                  
                  x = geom % x(:,i,j,k)
                  e = MAX(e,MAXVAL(ABS(p-x)))
               END DO   
            END DO
         END DO  
         
         CALL FTAssertEqual(expectedValue = 0.0_RP, &
                            actualValue = e,        &
                            tol = 1.d-9,           &
                            msg = "Mapped Geometry nodal locations max Error")
!
!        -----------
!        Face values
!        -----------
!
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([-1.0_RP, spA % eta(i), spA % zeta(j)])
               x    = geom % xb(:,i,j,ELEFT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Left face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([1.0_RP, spA % eta(i), spA % zeta(j)])
               x    = geom % xb(:,i,j,ERIGHT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Right face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j), -1.0_RP ])
               x    = geom % xb(:,i,j,EBOTTOM)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Bottom face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i), spA % eta(j), 1.0_RP ])
               x    = geom % xb(:,i,j,ETOP)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Top face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i), -1.0_RP, spA % zeta(j)])
               x    = geom % xb(:,i,j,EFRONT)
               e    = MAXVAL(ABS(x - xMap))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Front face locations")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               xMap = mapper % transfiniteMapAt([spA % xi(i),  1.0_RP, spA % zeta(j)])
               x    = geom % xb(:,i,j,EBACK)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Back face locations")
!
!        --------------------------------
!        Face values - should match faces
!        --------------------------------
!
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL ComputeFacePoint(self = mapper % faces(6),u = [spA % eta(i), spA % zeta(j)],p = xMap)
               x    = geom % xb(:,i,j,ELEFT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Left face locations matching input face")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL ComputeFacePoint(self = mapper % faces(4),u = [spA % eta(i), spA % zeta(j)],p = xMap)
               x    = geom % xb(:,i,j,ERIGHT)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Right face locations matching input face")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL ComputeFacePoint(self = mapper % faces(3),u = [spA % eta(i), spA % zeta(j)],p = xMap)
               x    = geom % xb(:,i,j,EBOTTOM)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Bottom face locations matching input face")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL ComputeFacePoint(self = mapper % faces(5),u = [spA % eta(i), spA % zeta(j)],p = xMap)
               x    = geom % xb(:,i,j,ETOP)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Top face locations matching input face")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL ComputeFacePoint(self = mapper % faces(1),u = [spA % eta(i), spA % zeta(j)],p = xMap)
               x    = geom % xb(:,i,j,EFRONT)
               e    = MAXVAL(ABS(x - xMap))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Front face locations matching input face")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL ComputeFacePoint(self = mapper % faces(2),u = [spA % eta(i), spA % zeta(j)],p = xMap)
               x    = geom % xb(:,i,j,EBACK)
               e    = MAX(MAXVAL(ABS(x - xMap)),e)
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Back face locations matching input face")
!
!        -------
!        Normals
!        -------
!
         e     = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL cylindricalFaceNormals([ spA % xi(i), -1.0_RP,  spA % zeta(j)], p, 1)
               x    = geom % normal(:,i,j,EFRONT)
               e    = MAX(e,MAXVAL(ABS(x - p)))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Error in normal at face 1")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL cylindricalFaceNormals([ spA % xi(i), 1.0_RP,  spA % zeta(j)], p, 2)
               x    = geom % normal(:,i,j,EBACK)
               e    = MAX(e,MAXVAL(ABS(x - p)))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Error in normal at face 2")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL cylindricalFaceNormals([ spA % xi(i), spA % zeta(j), -1.0_RP], p, 3)
               x    = geom % normal(:,i,j,EBOTTOM)
               e    = MAX(e,MAXVAL(ABS(x - p)))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Error in normal at face 3")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL cylindricalFaceNormals([ spA % xi(i), spA % zeta(j), 1.0_RP], p, 5)
               x    = geom % normal(:,i,j,ETOP)
               e    = MAX(e,MAXVAL(ABS(x - p)))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Error in normal at face 5")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
               CALL cylindricalFaceNormals([ -1.0_RP, spA % xi(i), spA % zeta(j)], p, 6)
               x    = geom % normal(:,i,j,ELEFT)
               e    = MAX(e,MAXVAL(ABS(x - p)))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Error in normal at face 6")
         e = 0.0_RP
         DO j = 0, N
            DO i = 0, N
                CALL cylindricalFaceNormals([ 1.0_RP, spA % xi(i), spA % zeta(j)], p, 4)
              x    = geom % normal(:,i,j,ERIGHT)
               e    = MAX(e,MAXVAL(ABS(x - p)))
            END DO   
         END DO   
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Error in normal at face 4")
         
!
!        -----------------------
!        Metric terms -- jGradXi
!        -----------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  x    = geom % jGradXi(:,i,j,k)
                  CALL cylindricalMetricTerms([ spA % xi(i), spA % eta(j),  spA % zeta(k)],Jai)
                  e    = MAX(e,MAXVAL(ABS(x - jAi(:,1))))
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 2.d-9,              &
                            msg           = "Computation of JGradXi")
!
!        ------------------------
!        Metric terms -- jGradEta
!        ------------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  x    = geom % jGradEta(:,i,j,k)
                  CALL cylindricalMetricTerms([ spA % xi(i), spA % eta(j),  spA % zeta(k)],Jai)
                  e    = MAX(e,MAXVAL(ABS(x - jAi(:,2))))
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of JGradEta")
!
!        ------------------------
!        Metric terms -- jGradZeta
!        ------------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  x    = geom % jGradZeta(:,i,j,k)
                  CALL cylindricalMetricTerms([ spA % xi(i), spA % eta(j),  spA % zeta(k)],Jai)
                  e    = MAX(e,MAXVAL(ABS(x(:) - jAi(:,3))))
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 5.d-9,              &
                            msg           = "Computation of jGradZeta")
!
!        ------------------------
!        Metric terms -- jacobian
!        ------------------------
!
         e = 0.0_RP
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
                  jac      = geom % jacobian(i,j,k)
                  exactJac = cylindricalJacobian([ spA % xi(i), spA % eta(j),  spA % zeta(k)])
                  e    = MAX(e,ABS(jac - exactJac))
               END DO   
            END DO   
         END DO
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-9,              &
                            msg           = "Computation of jacobian")
!
!        ------------
!        Plot results
!        ------------
!
         CALL plotGeometry(geom,"CylinderGeometry.tec")
         
      END SUBROUTINE cylinderTestGeometry
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE constructHexMapper(NKnots, mapper, mapSubroutine)  
         USE MappedGeometryClass
         IMPLICIT NONE  
         INTEGER                 :: NKnots
         TYPE(TransfiniteHexMap) :: mapper
         EXTERNAL                :: mapSubroutine
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
         REAL(KIND=RP)                      :: x(3)
         
         INTEGER                            :: i, j
         
         nUknots = NKnots
         nVKnots = NKnots
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
               CALL mapSubroutine([uKnots(i), -1.0_RP, vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(1), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL mapSubroutine([uKnots(i), 1.0_RP, vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(2), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL mapSubroutine([uKnots(i), vKnots(j),-1.0_RP], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(3), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL mapSubroutine([uKnots(i), vKnots(j),1.0_RP], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(5), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL mapSubroutine([-1.0_RP, uKnots(i), vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(6), uKnots, vKnots, points )

         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL mapSubroutine([1.0_RP, uKnots(i), vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(4), uKnots, vKnots, points )
!
!        ---------------------
!        Construct the mapper 
!        ---------------------
!
         CALL mapper % constructWithFaces(faceData)
         
      END SUBROUTINE constructHexMapper
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
      SUBROUTINE cylindricalFaceNormals(u, nHat, face)  
         USE SMConstants
         IMPLICIT NONE  
         INTEGER       :: face
         REAL(KIND=RP) :: u(3)    ! IN [-1,1]^3
         REAL(KIND=RP) :: nHat(3) ! physical space locations
         
         REAL(KIND=RP) :: r0 = 1.0_RP, rMax = 2.0_RP, theta0 = 0.0_RP, &
                          thetaMax = PI/2, z0 = 0.0_RP, zMax = 3.0_RP
         REAL(KIND=RP) :: r, theta, z

         r     = r0 + 0.5_RP*(u(1)+1.0_RP)*(rMax - r0)
         theta = theta0 + 0.5_RP*(u(2)+1.0_RP)*(thetaMax - theta0)
         z     = z0 + 0.5_RP*(u(3)+1.0_RP)*(zMax - z0)
         
         SELECT CASE ( face )
            CASE (6) 
               nHat  = -[ COS(theta), SIN(theta), 0.0_RP];
            CASE (4)
                nHat =  [ COS(theta), SIN(theta), 0.0_RP];
            CASE(3)
                nHat = -[0.0_RP, 0.0_RP, 1.0_RP];
            CASE (5)
                nHat =  [0.0_RP, 0.0_RP, 1.0_RP];
            CASE (1)
                nHat = [SIN(theta), -COS(theta), 0.0_RP];
            CASE(2)
                nHat =  [-SIN(theta), COS(theta), 0.0_RP];
            CASE DEFAULT 
         END SELECT 
         
      END SUBROUTINE cylindricalFaceNormals
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cylindricalMetricTerms(u,Jai)  
         USE SMConstants
         IMPLICIT NONE 
         REAL(KIND=RP) :: u(3)
         REAL(KIND=RP) :: Jai(3,3)
         
         REAL(KIND=RP) :: r0 = 1.0_RP, rMax = 2.0_RP, theta0 = 0.0_RP, &
                          thetaMax = PI/2, z0 = 0.0_RP, zMax = 3.0_RP
         REAL(KIND=RP) :: r, theta, z, drDxi, dthetaDEta, dZdZeta
         
         r     = r0     + 0.5_RP*(u(1)+1.0_RP)*(rMax - r0)
         theta = theta0 + 0.5_RP*(u(2)+1.0_RP)*(thetaMax - theta0)
         z     = z0     + 0.5_RP*(u(3)+1.0_RP)*(zMax - z0)
         
         drDxi      = 0.5_RP*(rMax - r0)
         dthetaDEta = 0.5_RP*(thetaMax - theta0)
         dZdZeta    = 0.5_RP*(zMax - z0)
!
         Jai(:,1) = [r*COS(theta)*dthetaDEta*dZdZeta, r*SIN(theta)*dthetaDEta*dZdZeta, 0.0_RP]
         Jai(:,2) = [-drDxi*SIN(theta)*dZdZeta, drDxi*COS(theta)*dZdZeta, 0.0_RP]
         Jai(:,3) = [0.0_RP, 0.0_RP, r*drDxi*dthetaDEta]
      END SUBROUTINE cylindricalMetricTerms
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION cylindricalJacobian(u)  RESULT(jac)
         USE SMConstants
         IMPLICIT NONE 
         REAL(KIND=RP) :: u(3)
         REAL(KIND=RP) :: jac
         
                  REAL(KIND=RP) :: r0 = 1.0_RP, rMax = 2.0_RP, theta0 = 0.0_RP, &
                          thetaMax = PI/2, z0 = 0.0_RP, zMax = 3.0_RP
         REAL(KIND=RP) :: r, theta, z, drDxi, dthetaDEta, dZdZeta
         
         r     = r0     + 0.5_RP*(u(1)+1.0_RP)*(rMax - r0)
         theta = theta0 + 0.5_RP*(u(2)+1.0_RP)*(thetaMax - theta0)
         z     = z0     + 0.5_RP*(u(3)+1.0_RP)*(zMax - z0)
         
         drDxi      = 0.5_RP*(rMax - r0)
         dthetaDEta = 0.5_RP*(thetaMax - theta0)
         dZdZeta    = 0.5_RP*(zMax - z0)
         
         jac = r*drDxi*dthetaDEta*dZdZeta

      END FUNCTION cylindricalJacobian
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE plotGeometry(geom, fName)
         USE MappedGeometryClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry) :: geom
         CHARACTER(LEN=*)     :: fName
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                 :: i, j, k, N, M, L
         INTEGER                 :: iUnit
         INTEGER                 :: iMax, jMax, kMax
         REAL(KIND=RP)           :: u, v, w, du, dv, dw, x(3)
         INTEGER, EXTERNAL       :: UnusedUnit
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
         iMax = geom % Nx
         jMax = geom % Ny
         kMax = geom % Nz
         
         WRITE(iUnit,*) 'VARIABLES = "X", "Y", "Z"'
         WRITE(iUnit,*) 'ZONE F=POINT, I=',iMax+1,', J=',jMax+1,', k=',kMax+1
         
         DO k = 0, kMax
            DO j = 0, jMax
               DO i = 0, iMax
                  WRITE(iUnit,*) geom % x(:,i,j,k)
               END DO   
            END DO
         END DO  
         
         CLOSE (iUnit)

      END SUBROUTINE plotGeometry
