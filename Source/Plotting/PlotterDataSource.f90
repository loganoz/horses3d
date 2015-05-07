!
!////////////////////////////////////////////////////////////////////////
!
!      PlotterDataSource.f90
!      Created: 2011-08-17 12:34:31 -0400 
!      By: David Kopriva  
!
!      The PlotterDataSource supplies the data needed for the
!      plotter. This can be generalized once bound procedures become
!      more common in fortran compilers.
!
!////////////////////////////////////////////////////////////////////////
!
      Module PlotterDataSource 
      IMPLICIT NONE
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      INTEGER FUNCTION NumberOfOutputVariables()
         IMPLICIT NONE 
         NumberOfOutputVariables = 3
      END FUNCTION NumberOfOutputVariables
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION Title() 
      IMPLICIT NONE 
      title = ' TITLE = "Euler Flow" '
      END FUNCTION Title
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION OutputVariableNames() 
         IMPLICIT NONE
         
         OutputVariableNames = ' VARIABLES = "x","y","p","Mach","Entropy"'
      END FUNCTION OutputVariableNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE OutputVectorFromStateVector( outputVector, stateVector)
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
         
      END SUBROUTINE OutputVectorFromStateVector
      
      END Module PlotterDataSource
