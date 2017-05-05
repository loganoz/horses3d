!
!////////////////////////////////////////////////////////////////////////
!
!      Setup.f90
!      Created: Jan 11, 2017 at 09:10 AM 
!      By: Juan Manzanero (extracted from David Kopriva's NSLite3DMain.f90)
!
!////////////////////////////////////////////////////////////////////////
!
      Module mainKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: cflKey                  = "cfl"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: meshFileNameKey         = "mesh file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartKey              = "restart"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartFileNameKey      = "restart file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveFileNameKey         = "save file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfTimeStepsKey    = "number of time steps"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: outputIntervalKey       = "output interval"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: convergenceToleranceKey = "convergence tolerance"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfBoundariesKey   = "number of boundaries"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: plotFileNameKey         = "plot file name"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(9) :: mainKeywords =  [cflKey,                  &
                                                                          meshFileNameKey,         &
                                                                          restartKey,              &
                                                                          restartFileNameKey,      &
                                                                          saveFileNameKey,         &
                                                                          numberOfTimeStepsKey,    &
                                                                          outputIntervalKey,       &
                                                                          convergenceToleranceKey, &
                                                                          plotFileNameKey]
      END MODULE mainKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
