!
!////////////////////////////////////////////////////////////////////////
!
!      The PlotterDataSourceClass supplies the data needed for the
!      plotter. This can be generalized once bound procedures become
!      more common in fortran compilers.
!
!////////////////////////////////////////////////////////////////////////
!
      Module PlotterDataSourceClass 
      IMPLICIT NONE
      
      TYPE PlotterDataSource
         CONTAINS
         PROCEDURE, NOPASS :: numberOfOutputVariables
         PROCEDURE, NOPASS :: title
         PROCEDURE, NOPASS :: outputVariableNames
         PROCEDURE, NOPASS :: outputVectorFromStateVector
      END TYPE PlotterDataSource
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      INTEGER FUNCTION numberOfOutputVariables()
         IMPLICIT NONE 
         numberOfOutputVariables = 14
      END FUNCTION numberOfOutputVariables
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION title() 
      IMPLICIT NONE 
      title = ' TITLE = "Free-stream preservation Euler equations" '
      END FUNCTION title
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION outputVariableNames() 
         IMPLICIT NONE
         
         outputVariableNames = &
         ' VARIABLES = "x","y","z","rho","rhou","rhov","rhow","rhoe","u","v","w","p","S_rho","S_rhou","S_rhov","S_rhow","S_rhoe"'
      END FUNCTION outputVariableNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE outputVectorFromStateVector( outputVector, stateVector, x)
         USE SMConstants
         USE BoundaryConditionFunctions
         USE Physics
         IMPLICIT NONE
         REAL(KIND=RP), DIMENSION(:) :: outputVector, stateVector
         REAL(KIND=RP) :: x(3)
         REAL(KIND=RP) :: Source(N_EQN)
         
         REAL(KIND=RP)               :: rho, u, v, w, e, P
         
         rho = stateVector(1)
         P   = Pressure(stateVector)
         Source=0._RP
         
         IF (flowIsNavierStokes) THEN
            CALL ManufacturedSolutionSourceNS ( x, 0._RP, Source  )
         ELSE
            CALL ManufacturedSolutionSourceEuler( x, 0._RP, Source  )
         END IF
         
         outputVector(1) = rho
         outputVector(2) = stateVector(2) 
         outputVector(3) = stateVector(3) 
         outputVector(4) = stateVector(4) 
         outputVector(5) = stateVector(5) 
         outputVector(6) = stateVector(2) / rho
         outputVector(7) = stateVector(3) / rho
         outputVector(8) = stateVector(4) / rho
         outputVector(9) = P 
         outputVector(10) = Source(1)
         outputVector(11) = Source(2)
         outputVector(12) = Source(3)
         outputVector(13) = Source(4)
         outputVector(14) = Source(5)
         
      END SUBROUTINE outputVectorFromStateVector
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION outputStateFromStateVector( indx, stateVector ) RESULT(state)
         USE SMConstants
         IMPLICIT NONE
         INTEGER                     :: indx
         REAL(KIND=RP), DIMENSION(:) :: stateVector
         REAL(KIND=RP)               :: state
         
         state = stateVector(indx)
          
      END FUNCTION outputStateFromStateVector
      
      END Module PlotterDataSourceClass
