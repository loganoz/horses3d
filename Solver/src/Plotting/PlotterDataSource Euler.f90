!
!////////////////////////////////////////////////////////////////////////
!
!      PlotterDataSourceClass.f90
!      Created: 2011-08-17 12:34:31 -0400 
!      By: David Kopriva  
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
         PROCEDURE :: numberOfOutputVariables
         PROCEDURE :: title
         PROCEDURE :: outputVariableNames
         PROCEDURE :: outputVectorFromStateVector
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
         numberOfOutputVariables = 3
      END FUNCTION numberOfOutputVariables
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION title() 
      IMPLICIT NONE 
      title = ' TITLE = "Euler Flow" '
      END FUNCTION title
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION outputVariableNames() 
         IMPLICIT NONE
         
         outputVariableNames = ' VARIABLES = "x","y","p","Mach","Entropy"'
      END FUNCTION outputVariableNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE outputVectorFromStateVector( outputVector, stateVector)
         USE SMConstants
         USE PDEModule, ONLY: Temperature, Pressure, gamma, gammaM2
         IMPLICIT NONE
         REAL(KIND=RP), DIMENSION(:) :: outputVector, stateVector
         REAL(KIND=RP)               :: q2, a2
         
         outputVector(1) = Pressure(stateVector)
         
         q2 = (stateVector(2)**2 + stateVector(3)**2)/stateVector(1)**2
         a2 = gamma*outputVector(1)/stateVector(1)
         
         outputVector(2) = SQRT(q2/a2)
         outputVector(3) = outputVector(1)/stateVector(1)**gamma - 1.0_RP/gammaM2
         
      END SUBROUTINE outputVectorFromStateVector
      
      END Module PlotterDataSourceClass
