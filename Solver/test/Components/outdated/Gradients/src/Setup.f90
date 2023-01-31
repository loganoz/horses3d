!
!////////////////////////////////////////////////////////////////////////
!
!      Setup.f90
!      Created: June 19, 2015 at 12:52 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      Module mainKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: discretizationNodesKey  = "discretization nodes"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveGradientsToSolutionKey  = "save gradients with solution"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: inviscidDiscretizationKey  = "inviscid discretization"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: viscousDiscretizationKey  = "viscous discretization"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: splitFormKey               = "split form"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: machNumberKey           = "mach number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: reynoldsNumberKey       = "reynolds number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: aoaThetaKey             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: aoaPhiKey               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: flowIsNavierStokesKey   = "flowisnavierstokes"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: polynomialOrderKey      = "polynomial order"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: cflKey                  = "cfl"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: meshFileNameKey         = "mesh file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartKey              = "restart"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartFileNameKey      = "restart file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfTimeStepsKey    = "number of time steps"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: outputIntervalKey       = "output interval"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: convergenceToleranceKey = "convergence tolerance"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfPlotPointsKey   = "number of plot points"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfBoundariesKey   = "number of boundaries"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: solutionFileNameKey     = "solution file name"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(20) :: mainKeywords =  [machNumberKey,           &
                                                                          discretizationNodesKey,  &
                                                                          inviscidDiscretizationKey, &
                                                                          viscousDiscretizationKey, &
                                                                          splitFormKey, &
                                                                          saveGradientsToSolutionKey, &
                                                                          reynoldsNumberKey,       &
                                                                          aoaThetaKey,             &
                                                                          aoaPhiKey,               &
                                                                          flowIsNavierStokesKey,   &
                                                                          polynomialOrderKey,      &
                                                                          cflKey,                  &
                                                                          meshFileNameKey,         &
                                                                          restartKey,              &
                                                                          restartFileNameKey,      &
                                                                          numberOfTimeStepsKey,    &
                                                                          outputIntervalKey,       &
                                                                          convergenceToleranceKey, &
                                                                          numberOfPlotPointsKey,   &
                                                                          solutionFileNameKey]
      END MODULE mainKeywordsModule

