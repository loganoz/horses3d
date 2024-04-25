!
! ////////////////////////////////////////////////////////////////////
! HORSES3D Converter
!     Main program horsesConverter
!	This add-ons has 3 capabilities:
!      i  . convertSolutionToDiffMesh : Convert .hsol from NS simulation to another mesh (.hmesh) with identical geometry
!      ii . convertHorses2OFMesh      : Convert .hmesh file to OpenFOAM polyMesh with element can be discretised into p-order GL nodes
!      iii. convertOFVTK2Horses       : Convert vtk result from OpenFOAM at which the vtk has points data to .hsol 
!                                       (cell to point data interpolation are performed in OpenFOAM and identical nodes)
!
!////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
    PROGRAM horsesConverter
        USE SMConstants
        USE SharedSpectralBasis
        USE getTaskConverter
        USE convertSolution
		USE MPI_Process_Info
		USE convertMesh2OF
		USE convertVTK2Horses
        IMPLICIT NONE
        character(LEN=LINE_LENGTH),  parameter  :: CONFIG_FILE ="horsesConverter.convert"
        character(LEN=LINE_LENGTH)              :: meshFile1, boundaryFile1, resultFile1
		character(LEN=LINE_LENGTH)              :: meshFile2, boundaryFile2, VTKfile
        INTEGER                                 :: polyOrder(3), tasktype
		REAL(KIND=RP)           				:: Ref(4)
		
		CALL MPI_Process % Init
!
!  Write Header Log
!  ------------------------	
        write(STD_OUT,'(A,A)') "/*------------------------------------- HORSES3D - CONVERTER -------------------------------------*\"
        write(STD_OUT,'(A,A)') "####################################################################################################"
        write(STD_OUT,'(A,A)') "#                                                                                                  #"
        write(STD_OUT,'(A,A)') "#            HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver             #"
        write(STD_OUT,'(A,A)') "#        Convert .hsol to different Mesh, Mesh to OpenFOAM polyMesh, OpenFOAM VTK to .hsol         #"
        write(STD_OUT,'(A,A)') "#                            input control file : horsesConverter.convert                          #"
        write(STD_OUT,'(A,A)') "#                                                                                                  #"
        write(STD_OUT,'(A,A)') "####################################################################################################"
        write(STD_OUT,'(A,A)') "\*------------------------------------------------------------------------------------------------*/"
!
!  Read File Input "horsesConverter.convert"
!  ------------------------
        CALL readGetTaskInput(CONFIG_FILE, meshFile1, boundaryFile1, resultFile1, meshFile2, boundaryFile2, polyOrder, tasktype, VTKfile, Ref)
!
!  Construct Spectral basis
!  ------------------------
        CALL ConstructSpectralBasis
		
		IF (tasktype==1) then
!
!  			Perform .hsol convertion to different Mesh 
!  			------------------------------------------
			CALL convertSolutionToDiffMesh (meshFile1, boundaryFile1, resultFile1, meshFile2, boundaryFile2, polyOrder)
			
		ELSE IF (tasktype==2) then
!
!  			Convert Horses Mesh to OpenFOAM Mesh 
!  			------------------------------------
			CALL convertHorses2OFMesh (meshFile1, boundaryFile1,polyOrder)
		ELSE IF (tasktype==3) then
!
!  			Convert Openfoam VTK result to Horses 
!  			------------------------------------
			CALL convertOFVTK2Horses (meshFile1, boundaryFile1, polyOrder, VTKfile, Ref)	
	    ELSE	
			write(STD_OUT,'(A,A)') "ERROR--Task Keyword is not recognized"	
			write(STD_OUT,'(A,A)') "Available Task Keyword: meshInterpolation, horsesMesh2OF, or OF2Horses"	
			CALL EXIT(0)
	    END IF 

!
!  Destruct Spectral basis
!  -----------------------
        call DestructSpectralBasis
		
        write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
        write(STD_OUT,'(A,A)') "/*------------------------------------------- END PROGRAM ------------------------------------------*\"
        write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
    END PROGRAM horsesConverter
!
!////////////////////////////////////////////////////////////////////////
!   
