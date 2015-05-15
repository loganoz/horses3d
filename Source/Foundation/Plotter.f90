!
!////////////////////////////////////////////////////////////////////////
!
!      Plotter.f90
!      Created: 2011-08-17 09:19:37 -0400 
!      By: David Kopriva  
!
!      The Plotter class separates the formatting of the data for
!      the plotting program from the data itself. Two types of
!      plotters are available: Interpolating or noninterpolating, 
!      though for DGSEM approximation with Gauss points the only
!      useful one is interpolating.
!
!      To construct a plotter:
!         ConstructPlotter( plotter, fUnit )
!         ConstructInterpolatingPlotter( plotter, fUnit, newN, newM )
!
!       Where,
!         fUnit = file unit
!         newN, newM = number of interpolation points
!
!      To use a plotter, call
!         ExportTocTecplot( plotter, sem )
!
!      where sem is the spectral element approximation of type DGSem.
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE PlotterClass
      USE DGSEMClass
      USE PlotterDataSource
      IMPLICIT NONE
      
      TYPE Plotter
         INTEGER :: fUnit
         LOGICAL :: interpolate
         INTEGER :: newN, newM
      END TYPE Plotter
!
!     ========   
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructPlotter( this, fUnit )
         IMPLICIT NONE 
         TYPE(Plotter) :: this
         INTEGER       :: fUnit
         this % interpolate = .false.
         this % fUnit       = fUnit
         this % newN        = 0
         this % newM        = 0
      END SUBROUTINE ConstructPlotter
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructInterpolatingPlotter( this, fUnit, newN, newM )
         IMPLICIT NONE 
         TYPE(Plotter) :: this
         INTEGER       :: fUnit
         INTEGER       :: newN, newM
         this % interpolate = .true.
         this % fUnit       = fUnit
         this % newN        = newN
         this % newM        = newM
      END SUBROUTINE ConstructInterpolatingPlotter
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplot(this, sem) 
      IMPLICIT NONE 
         
         TYPE(Plotter)           :: this
         TYPE(DGSem)             :: sem
         INTEGER                 :: N
         
         IF(this % interpolate)     THEN
            N     = SIZE(sem % mesh % elements(1) % Q,1) - 1
            CALL ExportToTecplotI( this, N, sem )
         ELSE
            CALL ExportToTecplotNoI( this, sem )
         END IF 
         
      END SUBROUTINE ExportToTecplot
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplotNoI( this, sem ) 
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(Plotter) :: this
         TYPE(DGSem)   :: sem
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                                  :: i, j, k, N, nPltVars, l
         CHARACTER(LEN=32)                        :: fmtString
         REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: outputVector
!
         nPltVars = NumberOfOutputVariables()
         ALLOCATE( outputVector(nPltVars) )
         fmtString = FormatString()
         
         WRITE(this % fUnit,*) Title()
         WRITE(this % fUnit,*) OutputVariableNames()
         
         N       = SIZE(sem % mesh % elements(1) % Q,1)
         
         DO k = 1, SIZE(sem % mesh % elements) 
            WRITE(this % fUnit,*) "ZONE I=", N+1, ",J=",N+1,", F=POINT"
            DO j= 0, N 
               DO i = 0, N
                  CALL OutputVectorFromStateVector( outputVector, sem % mesh % elements(k) % Q(i,j,:) )
                  
                  WRITE(this % fUnit,fmtString) sem % mesh % elements(k) % geom % x(i,j), &
                                              sem % mesh % elements(k) % geom % y(i,j), &
                                              (outputVector(l), l = 1, nPltVars)
               END DO
            END DO
         END DO
         
         DEALLOCATE( outputVector )

      END SUBROUTINE ExportToTecplotNoI
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplotI( this, N, sem ) 
      IMPLICIT NONE
      
         
      TYPE(Plotter) :: this
      TYPE(DGSem)   :: sem
      INTEGER       :: N
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER        :: fUnit
      INTEGER        :: newN, newM, nPltVars
      INTEGER        :: i, j, k, l
      REAL(KIND=RP), DIMENSION(0:this % newN, 0:this % newM ,2) :: fineMesh
      REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE          :: fineMeshSolution
      CHARACTER(LEN=32)                                     :: valuesFmt
!
!     ------------------------
!     Interpolation quantities
!     ------------------------
!
      REAL(KIND=RP), DIMENSION(0:this % newN,0:N) :: interpMatX
      REAL(KIND=RP), DIMENSION(0:this % newM,0:N) :: interpMatY
      REAL(KIND=RP), DIMENSION(0:N)             :: wx, wy
      REAL(KIND=RP), DIMENSION(0:this % newN)     :: newNodesX
      REAL(KIND=RP), DIMENSION(0:this % newM)     :: newNodesY
      REAL(KIND=RP)                             :: dx, dy
!
      newN     = this % newN
      newM     = this % newM
      fUnit    = this % fUnit
      nPltVars = NumberOfOutputVariables()
      ALLOCATE( fineMeshSolution(0:this % newN, 0:this % newM, nPltVars))
!
!     --------------------
!     Set up interpolation
!     --------------------
!
      dx = 2.0_RP/newN
      DO i = 0, newN 
         newNodesX(i) = -1.0_RP + i*dx
      END DO
      CALL BarycentricWeights( N, sem % spA % xi, wx )
      CALL PolynomialInterpolationMatrix( N, newN, sem % spA % xi, wx, newNodesX, interpMatX)
!
      dy = 2.0_RP/newM
      DO j = 0, newM 
         newNodesY(j) = -1.0_RP + j*dy
      END DO
      CALL BarycentricWeights( N, sem % spA % xi, wy )
      CALL PolynomialInterpolationMatrix( N, newM, sem % spA % xi, wy, newNodesY, interpMatY)
!
!     ----------------------------------
!     Interpolate and output the results
!     ----------------------------------
!
      WRITE(fUnit,*) TRIM( Title() )
      WRITE(fUnit,*) TRIM( OutputVariableNames() )
      
      valuesFMT = FormatString()
      
      DO k = 1, SIZE(sem % mesh % elements) 
      
         WRITE(fUnit,*) "ZONE I=", newN+1, ",J=",newM+1,", F=POINT"
!
!        ----------------------------------------------------------------------------------------
!        Rem: To avoid any jumps in the geometry at the element boundaries, don't interpolate
!        the locations. Instead, save the quadMap as part of the element and evaluate it here for
!        the NewX and newY points.
!        ----------------------------------------------------------------------------------------
!
         CALL InterpolateToFineMesh( sem, k, fineMeshSolution, fineMesh, N, newN, interpMatX, newM, interpMatY )
         DO j= 0, newM
            DO i = 0, newN
               WRITE(fUnit,valuesFmt) fineMesh(i,j,1), fineMesh(i,j,2), ( fineMeshSolution(i,j,l), l = 1, nPltVars )
            END DO
         END DO
         
      END DO
      
      DEALLOCATE( fineMeshSolution )

      END SUBROUTINE ExportToTecplotI
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=32) FUNCTION FormatString() RESULT(valuesFMT)
         IMPLICIT NONE
         INTEGER           :: nPltVars
         CHARACTER(LEN=32) :: valuesStr
         
         nPltVars = NumberOfOutputVariables()
         WRITE(valuesStr,*) nPltVars+2
         valuesFMT = "(" // TRIM(valuesStr) // "E13.5)"
      END FUNCTION FormatString
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE InterpolateToFineMesh( this, ID, interpSoln, newMesh, oldN, N, interpMatX, M, InterpMatY )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                              :: N, M, ID, oldN
      TYPE(DGSem)                          :: this
      REAL(KIND=RP), DIMENSION(0:,0:,1:)   :: interpSoln
      REAL(KIND=RP), DIMENSION(0:N,0:M,2)  :: newMesh
      REAL(KIND=RP), DIMENSION(0:N,0:oldN) :: interpMatX
      REAL(KIND=RP), DIMENSION(0:M,0:oldN) :: interpMatY
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(0:N,0:oldN)         :: uTmp
      REAL(KIND=RP), DIMENSION(0:N,0:oldN)         :: xTmp, yTmp
      REAL(KIND=RP), DIMENSION(0:M)                :: vNew
      REAL(KIND=RP), DIMENSION(0:oldN)             :: vOld
      INTEGER                                      :: i, j, k, nPlotValues
      REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: valuesToInterpolate
!
!     ------------------
!     Mesh interpolation
!     ------------------
!
      DO j = 0, oldN
         CALL InterpolateToNewPoints( oldN, N, interpMatX, this % mesh % elements(ID) % geom % x(:,j), xTmp(:,j) )
         CALL InterpolateToNewPoints( oldN, N, interpMatX, this % mesh % elements(ID) % geom % y(:,j), yTmp(:,j) )
      END DO
      
      DO i = 0, N 
         vOld = yTmp(i,:)
         CALL InterpolateToNewPoints( oldN, M, interpMatY, vOld , vNew )
         newMesh(i,:,2) = vNew
         vOld = xTmp(i,:)
         CALL InterpolateToNewPoints( oldN, M, interpMatY, vOld , vNew )
         newMesh(i,:,1) = vNew
      END DO
!
!     ------------------------------------
!     Compute variables to be interpolated
!     ------------------------------------
!
      nPlotValues = NumberOfOutputVariables()
      ALLOCATE( valuesToInterpolate(0:oldN,0:oldN, nPlotValues))
      DO j = 0, oldN
         DO i = 0, oldN
            CALL OutputVectorFromStateVector( valuesToInterpolate(i,j,:), this % mesh % elements(id) % Q(i,j,:) )
         END DO
      END DO
!
!     ----------------------
!     Solution interpolation
!     ----------------------
!
      DO k = 1, nPlotValues
         DO j = 0, oldN
            CALL InterpolateToNewPoints( oldN, N, interpMatX, valuesToInterpolate(:,j,k), uTmp(:,j) )
         END DO
         
         DO i = 0, N
            vOld = uTmp(i,:)
            CALL InterpolateToNewPoints( oldN, M, interpMatY, vOld , vNew )
            interpSoln(i,:,k) = vNew
         END DO
      END DO
      
      DEALLOCATE(valuesToInterpolate)

      END SUBROUTINE InterpolateToFineMesh
      
   END MODULE PlotterClass