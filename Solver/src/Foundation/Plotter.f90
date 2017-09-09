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
         
         TYPE(Interpolator_t), ALLOCATABLE    :: interpMats(:,:,:)     ! All possible interpolation matrices to output state (only needed ones are created)
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
         CLASS(NodalStorage), OPTIONAL     :: spA(0:,0:,0:)
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP)                     :: dx       ! Distance between plotted points
         INTEGER                           :: i,j,k    ! Counters
         INTEGER                           :: NxMax, & ! Maximum polynomial order in the x direction
                                              NyMax, & ! Maximum polynomial order in the y direction
                                              NzMax    ! Maximum polynomial order in the z direction
         
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
!           ---------------------------------------
!           Allocate spece for interpolation arrays
!           ---------------------------------------
!
            ALLOCATE(x(0:newN))
            dx = 2.0_RP/newN
            DO j = 0, newN
               x(j) = -1.0_RP + dx*j 
            END DO  
            
            NxMax = 0
            NyMax = 0
            NzMax = 0
            DO k=0, UBOUND(spA,3)
               DO j=0, UBOUND(spA,2)
                  DO i=0, UBOUND(spA,1)
                     IF (.NOT. spA(i,j,k) % Constructed) CYCLE
                     NxMax = MAX(NxMax,spA(i,j,k) % Nx)
                     NyMax = MAX(NyMax,spA(i,j,k) % Ny)
                     NzMax = MAX(NzMax,spA(i,j,k) % Nz)
                  END DO
               END DO
            END DO
            
            ALLOCATE(self % interpMats(0:NxMax, 0:NyMax, 0:NzMax))
            
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
         CLASS(NodalStorage)           :: spA(0:,0:,0:)
         
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
         INTEGER                                  :: i, j, k, id, Nxyz(3), nPltVars, l
         CHARACTER(LEN=32)                        :: fmtString
         REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: outputVector
!
         nPltVars  = self % dataSource % NumberOfOutputVariables()
         fmtString = FormatString()

         ALLOCATE( outputVector(nPltVars) )
         
         WRITE(self % fUnit,*) self % dataSource % Title()
         WRITE(self % fUnit,*) self % dataSource % OutputVariableNames()
         
         DO id = 1, SIZE(elements) 
            Nxyz = elements(id) % Nxyz
            WRITE(self % fUnit,*) "ZONE I=", Nxyz(1)+1, ",J=",Nxyz(2)+1, ",K=",Nxyz(3)+1,", F=POINT"
            DO k = 0, Nxyz(3)
               DO j = 0, Nxyz(2)
                  DO i = 0, Nxyz(1)
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
         CLASS(NodalStorage)         :: spA(0:,0:,0:)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                                  :: N        ! Polynomial order in destination mesh
         INTEGER                                  :: Nx,Ny,Nz ! Polynomial orders in solution mesh 
         INTEGER                                  :: i, j, k, id, nPltVars, l, m
         CHARACTER(LEN=32)                        :: fmtString
         REAL(KIND=RP), POINTER                   :: Interp(:,:)
         REAL(KIND=RP), ALLOCATABLE               :: array3DOld(:,:,:,:)       ! Temporary variable for storing the solution inside an element
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
            Nx = elements(id) % Nxyz(1)
            Ny = elements(id) % Nxyz(2)
            Nz = elements(id) % Nxyz(3)
            IF (.NOT. self % interpMats(Nx,Ny,Nz) % Created) THEN
               CALL Create3DInterpolationMatrix (self % interpMats(Nx,Ny,Nz) % Mat,           &                   !> Interpolation matrix
                                                 Nx, Ny, Nz,                                  &                   !< Origin orders
                                                 N , N , N ,                                  &                   !< Destination orders
                                                 spA(Nx,Ny,Nz) % xi, spA(Nx,Ny,Nz) % eta, spA(Nx,Ny,Nz) % zeta, & !< Origin nodes
                                                 x , x , x )                                                      !< Destination nodes
               self % interpMats(Nx,Ny,Nz) % Created = .TRUE.
            END IF
            Interp => self % interpMats(Nx,Ny,Nz) % Mat
            
            WRITE(self % fUnit,*) "ZONE I=", N+1, ",J=",N+1, ",K=",N+1,", F=POINT"
!
!           -------------------------
!           Interpolate to new points
!           -------------------------
!
            ! Coordinates:
            DO m = 1, 3
               CALL Interpolate3D(elements(id) % geom % x(m,:,:,:), newxyz(m,:,:,:), Interp, &
                                  Nx, Ny, Nz, N, N, N)
            END DO
            ! State:
            ALLOCATE(array3dOld(numberOfOutputVariables(),0:Nx,0:Ny,0:Nz))
            DO k = 0, Nz
               DO j = 0, Ny
                  DO i = 0, Nx
                     CALL self % dataSource % outputVectorFromStateVector(array3DOld(:,i,j,k), elements(id) % Q(i,j,k,:))
                  END DO  
               END DO   
            END DO
            
            DO m = 1, numberOfOutputVariables()
               CALL Interpolate3D(array3DOld(m,:,:,:),new3DState(m,:,:,:), Interp, &
                                  Nx, Ny, Nz, N, N, N)
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
