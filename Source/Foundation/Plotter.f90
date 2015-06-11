!
!////////////////////////////////////////////////////////////////////////
!
!      Plotter.f90
!      Created: 2011-08-17 09:19:37 -0400 
!      By: David Kopriva  
!
!      The Plotter class separates the formatting of the data for
!      the plotting program from the data itself. It uses a delegate
!      in the form of the dataSource pointer to actually determine 
!      what is plotted beyond the mesh positions.
!
!      Two types of
!      plotters are available: Interpolating or noninterpolating, 
!      though for DGSEM approximation with Gauss points the only
!      useful one is interpolating.
!
!      To construct a plotter:
!         ConstructPlotter( plotter, fUnit, dataSource, newM )
!
!       Where,
!         fUnit      = file unit of plotFile
!         dataSource = an instance of PlotterDatasource or a subclass thereof
!         newN       = number of interpolation points in each direction (optional)
!
!      To use a plotter, call
!         ExportTocTecplot( plotter, elements )
!
!      where elements is the array of elements from the mesh.
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE DGSEMPlotterClass
      USE NodalStorageClass
      USE ElementClass
      USE PlotterDataSourceClass
      IMPLICIT NONE
      
      TYPE DGSEMPlotter
         INTEGER                              :: fUnit
         LOGICAL                              :: interpolate
         INTEGER                              :: newN, oldN
         CLASS(PlotterDatasource), POINTER    :: dataSource
      
         REAL(KIND=RP), ALLOCATABLE, PRIVATE  :: interpMatrix(:,:)
         REAL(KIND=RP), ALLOCATABLE, PRIVATE  :: tmpOldVector(:)
         REAL(KIND=RP), ALLOCATABLE, PRIVATE  :: tmpNewVector(:)
         REAL(KIND=RP), ALLOCATABLE, PRIVATE  :: tmp3Darray(:,:,:)
!
!        ========         
         CONTAINS
!        ========
!      
         PROCEDURE :: Construct => ConstructPlotter
         PROCEDURE :: Destruct => DestructPlotter
         PROCEDURE :: ExportToTecplot
      END TYPE DGSEMPlotter
      
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: oldXYZ(:,:,:,:)
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: newXYZ(:,:,:,:)
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: new3DState(:,:,:,:)
      
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: array3DNew(:,:,:)
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: array3DOld(:,:,:)
!
!     ========   
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructPlotter( self, spA, fUnit, dataSource, newN )
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(DGSEMPlotter)               :: self
         CLASS(NodalStorage)               :: spA
         CLASS(PlotterDataSource), POINTER :: dataSource
         INTEGER                           :: fUnit
         INTEGER, OPTIONAL                 :: newN
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP), DIMENSION(0:spA%N) :: w
         REAL(KIND=RP), ALLOCATABLE        :: x(:)
         REAL(KIND=RP)                     :: dx
         INTEGER                           :: j
         
         self % interpolate = .false.
         self % fUnit       = fUnit
         self % newN        = 0
         self % oldN        = spA % N
         self % dataSource  => dataSource ! Weak ownership
!
!        -------------------------------------
!        Construct interpolations, if necesary
!        -------------------------------------
!
         IF ( PRESENT(newN) )     THEN
            self % interpolate = .TRUE.
            self % newN        = newN
!
!           -------------------------
!           Allocate temporary arrays
!           -------------------------
!
            ALLOCATE( self % tmp3Darray(0:newN,0:newN,0:newN) )
            ALLOCATE( self % tmpNewVector(0:newN))
            ALLOCATE( self % tmpOldVector(0:spA % N))
            
            ALLOCATE(oldXYZ(0:spA % N,0:spA % N,0:spA % N,3))
            ALLOCATE(newXYZ(0:newN,0:newN,0:newN,3))
            ALLOCATE(new3DState(0:newN,0:newN,0:newN,numberOfOutputVariables()))
            ALLOCATE(array3DNew(0:newN,0:newN,0:newN))
            ALLOCATE(array3dOld(0:spA % N,0:spA % N,0:spA % N))
!
!           -----------------------------
!           Generate interpolation arrays
!           -----------------------------
!
            ALLOCATE(x(0:newN))
            dx = 2.0_RP/newN
            DO j = 0, newN
               x(j) = -1.0_RP + dx*j 
            END DO  
            
            ALLOCATE(self % interpMatrix(0:newN, 0:spA % N))
            
            CALL BarycentricWeights( N = spA % N, x = spA % xi, w = w)
            CALL PolynomialInterpolationMatrix(N        = spA % N,           &
                                               M        = newN,              &
                                               oldNodes = spA % xi,          &
                                               weights  = w,                 &
                                               newNodes = x,                 &
                                               T        = self % interpMatrix)
         END IF 
         
      END SUBROUTINE ConstructPlotter
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE DestructPlotter(self)  
         IMPLICIT NONE  
         CLASS(DGSEMPlotter) :: self
      END SUBROUTINE DestructPlotter
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplot(self, elements) 
      IMPLICIT NONE 
         
         CLASS(DGSEMPlotter)           :: self
         TYPE(Element)           :: elements(:)
         
         IF(self % interpolate)     THEN
            CALL ExportToTecplotI( self, elements )
         ELSE
            CALL ExportToTecplotNoI( self, elements )
         END IF 
         
      END SUBROUTINE ExportToTecplot
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplotNoI( self, elements ) 
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(DGSEMPlotter) :: self
         TYPE(Element)       :: elements(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                                  :: i, j, k, id, N, nPltVars, l
         CHARACTER(LEN=32)                        :: fmtString
         REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: outputVector
!
         nPltVars  = self % dataSource % NumberOfOutputVariables()
         fmtString = FormatString()

         ALLOCATE( outputVector(nPltVars) )
         
         WRITE(self % fUnit,*) self % dataSource % Title()
         WRITE(self % fUnit,*) self % dataSource % OutputVariableNames()
         
         N       = UBOUND(elements(1) % Q,1)
         
         DO id = 1, SIZE(elements) 
            WRITE(self % fUnit,*) "ZONE I=", N+1, ",J=",N+1, ",K=",N+1,", F=POINT"
            DO k = 0, N
               DO j= 0, N 
                  DO i = 0, N
                     CALL self % dataSource % OutputVectorFromStateVector( outputVector, elements(id) % Q(i,j,k,:) )
                     
                     WRITE(self % fUnit,fmtString) elements(id) % geom % x(1,i,j,k), &
                                                   elements(id) % geom % x(2,i,j,k), &
                                                   elements(id) % geom % x(3,i,j,k), &
                                                   (outputVector(l), l = 1, nPltVars)
                  END DO
               END DO
            END DO 
         END DO
         
         DEALLOCATE( outputVector )

      END SUBROUTINE ExportToTecplotNoI
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplotI( self, elements ) 
         IMPLICIT NONE
      
         
         CLASS(DGSEMPlotter) :: self
         TYPE(Element)       :: elements(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                                  :: i, j, k, id, N, nPltVars, l
         CHARACTER(LEN=32)                        :: fmtString
         REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: outputVector
!
         N = self % newN
         nPltVars  = self % dataSource % NumberOfOutputVariables()
         ALLOCATE( outputVector(nPltVars) )
         
         fmtString = FormatString()
!
!        ------------------------
!        Output plot file  header
!        ------------------------
!
         WRITE(self % fUnit,*) self % dataSource % Title()
         WRITE(self % fUnit,*) self % dataSource % OutputVariableNames()
         
         DO id = 1, SIZE(elements)
            WRITE(self % fUnit,*) "ZONE I=", N+1, ",J=",N+1, ",K=",N+1,", F=POINT"
!
!           -------------------------
!           Interpolate to new points
!           -------------------------
!
            CALL interpolateElementToFineMesh(self      = self, &
                                              e         = elements(id),&
                                              newPoints = newXYZ, &
                                              newState  = new3DState)
!
!           --------------------
!           Write out new points
!           --------------------
!
            DO k = 0, N
               DO j= 0, N 
                  DO i = 0, N
                     WRITE(self % fUnit,fmtString) newxyz(i,j,k,1), &
                                                   newxyz(i,j,k,2), &
                                                   newxyz(i,j,k,3), &
                                                   (new3DState(i,j,k,l), l = 1, nPltVars)
                  END DO
               END DO
            END DO 
         END DO
!
      END SUBROUTINE ExportToTecplotI
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE interpolateElementToFineMesh( self, e, newPoints, newState)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSemPlotter) :: self 
         TYPE( Element )    :: e
         REAL(KIND=RP)      :: newPoints(0:,0:,0:,1:)
         REAL(KIND=RP)      :: newState(0:,0:,0:,1:)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                     :: i, j, k, m
         INTEGER                     :: oldN, newN
!         
         oldN = self % oldN
         newN = self % newN
!
!        ------------------
!        Mesh interpolation
!        ------------------
!
         DO m = 1, 3
            DO k = 0, oldN
               DO j = 0, oldN    
                  DO i = 0, oldN
                     array3DOld(i,j,k) = e % geom % x(m,i,j,k) 
                  END DO  
               END DO   
            END DO   
            CALL coarseToFineInterpolation3D(self          = self,     &
                                             old3DArrayArg = array3DOld, &
                                             new3DArrayArg = array3DNew) 
            newPoints(:,:,:,m) = array3DNew
         END DO
!
!        -------------------
!        State interpolation
!        -------------------
!
         DO m = 1, numberOfOutputVariables()
            DO k = 0, oldN
               DO j = 0, oldN    
                  DO i = 0, oldN
                     array3DOld(i,j,k) = outputStateFromStateVector(indx = m,stateVector = e % Q(i,j,k,:)) 
                  END DO  
               END DO   
            END DO   

            CALL coarseToFineInterpolation3D(self          = self,     &
                                             old3DArrayArg = array3DOld, &
                                             new3DArrayArg = array3DNew) 
            newState(:,:,:,m) = array3DNew
         END DO
         
      END SUBROUTINE interpolateElementToFineMesh
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE coarseToFineInterpolation3D( self, old3DArrayArg, new3DArrayArg )  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSEMPlotter)         :: self
         REAL(KIND=RP), INTENT(IN)  :: old3DArrayArg(0:,0:,0:)
         REAL(KIND=RP), INTENT(OUT) :: new3DArrayArg(0:,0:,0:)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k
         INTEGER :: oldN, newN
         
         oldN = self % oldN
         newN = self % newN
!
!        ------------------
!        Interpolate in xi 
!        ------------------
!
         DO k = 0, oldN
            DO j = 0, oldN
               self % tmpOldVector = old3DArrayArg(:,j,k)
               CALL InterpolateToNewPoints( N       = oldN, &
                                            M       = newN, &
                                            T       = self % interpMatrix,&
                                            f       = self % tmpOldVector,&
                                            fInterp = self % tmpNewVector)
                DO i = 0, newN
                   self % tmp3Darray(i,j,k) = self % tmpNewVector(i) 
                END DO  
            END DO
         END DO  
!
!        ------------------
!        Interpolate in eta
!        ------------------
!
         DO k = 0, oldN
            DO i = 0, newN
               DO j = 0, oldN
                  self % tmpOldVector(j) = self % tmp3Darray(i,j,k) 
               END DO
               CALL InterpolateToNewPoints( N       = oldN, &
                                            M       = newN, &
                                            T       = self % interpMatrix,&
                                            f       = self % tmpOldVector,&
                                            fInterp = self % tmpNewVector)
                DO j = 0, newN
                   self % tmp3Darray(i,j,k) = self % tmpNewVector(j) 
                END DO  
            END DO   
         END DO
!
!        -------------------
!        Interpolate in zeta
!        -------------------
!
         DO j = 0, newN
            DO i = 0, newN
               DO k = 0, oldN
                  self % tmpOldVector(k) = self % tmp3Darray(i,j,k) 
               END DO
               CALL InterpolateToNewPoints( N       = oldN, &
                                            M       = newN, &
                                            T       = self % interpMatrix,&
                                            f       = self % tmpOldVector,&
                                            fInterp = self % tmpNewVector)
                DO k = 0, newN
                   new3DarrayArg(i,j,k) = self % tmpNewVector(k) 
                END DO  
            END DO   
         END DO  
         
      END SUBROUTINE coarseToFineInterpolation3D
      
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=32) FUNCTION FormatString() RESULT(valuesFMT)
         IMPLICIT NONE
         INTEGER           :: nPltVars
         CHARACTER(LEN=32) :: valuesStr
         
         nPltVars = NumberOfOutputVariables()
         WRITE(valuesStr,*) nPltVars+3
         valuesFMT = "(" // TRIM(valuesStr) // "E13.5)"
      END FUNCTION FormatString
      
   END MODULE DGSEMPlotterClass