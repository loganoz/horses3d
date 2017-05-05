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
!         ExportToTecplot( plotter, elements )
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
         LOGICAL                              :: interpolate =.FALSE.
         INTEGER                              :: newN
         CLASS(PlotterDatasource), POINTER    :: dataSource
         
         TYPE(Interpolator_t), ALLOCATABLE    :: interpMats(:,:,:)
!
!        ========         
         CONTAINS
!        ========
!      
         PROCEDURE :: Construct => ConstructPlotter
         PROCEDURE :: Destruct => DestructPlotter
         PROCEDURE :: ExportToTecplot
      END TYPE DGSEMPlotter
      
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: newXYZ(:,:,:,:)         ! Coordinates to plot
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: new3DState(:,:,:,:)     ! Value to ploit
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: x(:)                    ! Points in local element of plotted mesh
      REAL(KIND=RP), ALLOCATABLE, PRIVATE   :: array3DOld(:,:,:)       ! Temporary variable for storing the solution inside an element
!
!     ========   
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructPlotter( self, fUnit, dataSource, newN, spA )
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(DGSEMPlotter)               :: self
         CLASS(PlotterDataSource), POINTER :: dataSource
         INTEGER                           :: fUnit
         INTEGER            , OPTIONAL     :: newN
         CLASS(NodalStorage), OPTIONAL     :: spA(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP)                     :: dx      ! Distance between plotted points
         INTEGER                           :: j       ! Counter
         INTEGER                           :: NxMax   
         
         self % fUnit       = fUnit
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
            ALLOCATE(newXYZ    (3,                         0:newN,0:newN,0:newN))
            ALLOCATE(new3DState(numberOfOutputVariables(), 0:newN,0:newN,0:newN))
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
            
            NxMax = 0
            DO j=0, UBOUND(spA,1)
               IF (.NOT. spA(j) % Constructed) CYCLE
               NxMax = MAX(NxMax,spA(j) % N)    ! TODO: Change when anisotropic polynomials are implemented
            END DO
            
            ALLOCATE(self % interpMats(0:NxMax, 0:NxMax, 0:NxMax))
            
         END IF 
         
      END SUBROUTINE ConstructPlotter
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE DestructPlotter(self)  
         IMPLICIT NONE  
         CLASS(DGSEMPlotter) :: self
         
         IF (ALLOCATED(self % interpMats)) DEALLOCATE (self % interpMats)
      END SUBROUTINE DestructPlotter
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExportToTecplot(self, elements, spA) 
      IMPLICIT NONE 
         
         CLASS(DGSEMPlotter)           :: self
         TYPE(Element)                 :: elements(:)
         CLASS(NodalStorage)           :: spA(:)
         
         IF(self % interpolate)     THEN
            CALL ExportToTecplotI( self, elements , spA)
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
         
         DO id = 1, SIZE(elements) 
            N = elements(id) % N
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
      SUBROUTINE ExportToTecplotI( self, elements, spA ) 
         IMPLICIT NONE
      
         
         CLASS(DGSEMPlotter), TARGET :: self
         TYPE(Element)               :: elements(:)
         CLASS(NodalStorage)         :: spA(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                                  :: N      ! Polynomial order in destination mesh
         INTEGER                                  :: Nx     ! Polynomial orders in solution mesh 
         INTEGER                                  :: i, j, k, id, nPltVars, l, m
         CHARACTER(LEN=32)                        :: fmtString
         REAL(KIND=RP), POINTER                   :: Interp(:,:)
!
         N = self % newN
         nPltVars  = self % dataSource % NumberOfOutputVariables()
         
         fmtString = FormatString()
!
!        ------------------------
!        Output plot file  header
!        ------------------------
!
         WRITE(self % fUnit,*) self % dataSource % Title()
         WRITE(self % fUnit,*) self % dataSource % OutputVariableNames()
         
         DO id = 1, SIZE(elements)
            Nx = elements(id) % N
            IF (.NOT. self % interpMats(Nx,Nx,Nx) % Created) THEN
               CALL Create3DInterpolationMatrix (self % interpMats(Nx,Nx,Nx) % Mat,           & !> Interpolation matrix
                                                 Nx, Nx, Nx,                                  & !< Origin orders
                                                 N , N , N ,                                  & !< Destination orders
                                                 spA(Nx) % xi, spA(Nx) % eta, spA(Nx) % zeta, & !< Origin nodes
                                                 x , x , x )                                    !< Destination nodes
               self % interpMats(Nx,Nx,Nx) % Created = .TRUE.
            END IF
            Interp => self % interpMats(Nx,Nx,Nx) % Mat
            
            WRITE(self % fUnit,*) "ZONE I=", N+1, ",J=",N+1, ",K=",N+1,", F=POINT"
!
!           -------------------------
!           Interpolate to new points
!           -------------------------
!
            ! Coordinates:
            DO m = 1, 3
               CALL Interpolate3D(elements(id) % geom % x(m,:,:,:), newxyz(m,:,:,:), Interp, &
                                  Nx, Nx, Nx, N, N, N)
            END DO
            ! State:
            ALLOCATE(array3dOld(0:Nx,0:Nx,0:Nx))
            DO m = 1, numberOfOutputVariables()
               DO k = 0, Nx
                  DO j = 0, Nx    
                     DO i = 0, Nx
                        array3DOld(i,j,k) = outputStateFromStateVector(indx = m,stateVector = elements(id) % Q(i,j,k,:)) 
                     END DO  
                  END DO   
               END DO
               CALL Interpolate3D(array3DOld(:,:,:),new3DState(m,:,:,:), Interp, &
                                  Nx, Nx, Nx, N, N, N)
            END DO
            DEALLOCATE(array3dOld)
!
!           --------------------
!           Write out new points
!           --------------------
!
            DO k = 0, N
               DO j= 0, N 
                  DO i = 0, N
                     WRITE(self % fUnit,fmtString) newxyz(1,i,j,k), &
                                                   newxyz(2,i,j,k), &
                                                   newxyz(3,i,j,k), &
                                                   (new3DState(l,i,j,k), l = 1, nPltVars)
                  END DO
               END DO
            END DO 
         END DO
         
      END SUBROUTINE ExportToTecplotI
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
