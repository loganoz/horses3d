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
         numberOfOutputVariables = 5
      END FUNCTION numberOfOutputVariables
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION title() 
      IMPLICIT NONE 
      title = ' TITLE = "Generic solution" '
      END FUNCTION title
!
!////////////////////////////////////////////////////////////////////////
!
      CHARACTER(LEN=132) FUNCTION outputVariableNames() 
         IMPLICIT NONE
         
         outputVariableNames = ' VARIABLES = "x","y","z","rho","rhou","rhov","rhow","rhoe"'
      END FUNCTION outputVariableNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE outputVectorFromStateVector( outputVector, stateVector)
         USE SMConstants
         IMPLICIT NONE
         REAL(KIND=RP), DIMENSION(:) :: outputVector, stateVector
         
         outputVector = stateVector(1:numberOfOutputVariables())
         
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
