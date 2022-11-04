!
! ////////////////////////////////////////////////////////////////////
!HORSES3D to Foam Result
!     main.f90
!
!!
!!     Modification History:
!!       version 0.0
!
!    Main Program horses2foam
!
!////////////////////////////////////////////////////////////////////////////////////////
!
    PROGRAM horses2foam
        USE SMConstants
        USE SharedSpectralBasis
        USE readHorses2Foam
        USE meshToFoam
        IMPLICIT NONE
        character(LEN=LINE_LENGTH),  parameter  :: CONFIG_FILE ="horses2foam.convert"
        character(LEN=LINE_LENGTH)              :: meshFile, boundaryFile, outputVariables
		character(LEN=LINE_LENGTH)           :: resultFile
        INTEGER                                 :: polyOrder(3)
		LOGICAL                                 :: generateMesh

!
!  Write Header Log
!  ------------------------	
        write(STD_OUT,'(A,A)') "/*-------------------------------------- HORSES3D - FOAM FILE --------------------------------------*\"
        write(STD_OUT,'(A,A)') "####################################################################################################"
        write(STD_OUT,'(A,A)') "#                                                                                                  #"
        write(STD_OUT,'(A,A)') "#            HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver             #"
        write(STD_OUT,'(A,A)') "#                              Foam File for Paraview Post-Processing                              #"
        write(STD_OUT,'(A,A)') "#                                                                                                  #"
        write(STD_OUT,'(A,A)') "####################################################################################################"
        write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
!
!  Read File Input "horses2foam.convert"
!  ------------------------
        CALL readHorses2FoamConvert(CONFIG_FILE, meshFile, boundaryFile, resultFile, outputVariables, polyOrder, generateMesh)

!
!  Construct Spectral basis
!  ------------------------
        CALL ConstructSpectralBasis
!
!  Generate Foam Mesh and Result
!  ------------------------
        CALL generateFoamMesh (meshFile, boundaryFile, resultFile, outputVariables, polyOrder, generateMesh)

!
!  Destruct Spectral basis
!  -----------------------
        call DestructSpectralBasis
		
        write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
        write(STD_OUT,'(A,A)') "/*------------------------------------------- END PROGRAM ------------------------------------------*\"
        write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
    END PROGRAM horses2foam
