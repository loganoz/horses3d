!
! ////////////////////////////////////////////////////////////////////
! HORSES3D Converter
!     Main program horsesConverter
!
!      This module convert Horses3D mesh storage into OpenFOAM polyMesh files
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE convertMesh2OF
    USE SMConstants
    USE InterpolationMatrices
    USE SharedSpectralBasis
    USE foamCreateMeshFileConverter
    USE Storage
    USE NodalStorageClass
	use SolutionFile
	use PhysicsStorage
	use convertSolution

!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE convertHorses2OFMesh (meshName, boundaryFile, Nout)
            IMPLICIT NONE
            CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: meshName
			CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: boundaryFile
			INTEGER, INTENT(IN)     				   :: Nout(NDIM)

!
!        ---------------
!        Local variables
!        ---------------
!
            type(Mesh_t)                               :: mesh
            integer                                    :: eID
            real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
            integer                                    :: i,fid, iSol
			integer                       			   :: pos, pos2
			character(len=LINE_LENGTH) 				   :: dir, time
			
!
!  Write Header Log
!  ------------------------	
        write(STD_OUT,'(A,A)') "/*-------------------------------------- HORSES3D - FOAM FILE --------------------------------------*\"
        write(STD_OUT,'(A,A)') "####################################################################################################"
        write(STD_OUT,'(A,A)') "#                                                                                                  #"
        write(STD_OUT,'(A,A)') "#            HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver             #"
        write(STD_OUT,'(A,A)') "#                            Convert HORSES3D mesh to OpenFOAM polyMesh                            #"
        write(STD_OUT,'(A,A)') "#                                                                                                  #"
        write(STD_OUT,'(A,A)') "####################################################################################################"
        write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
		write(STD_OUT,'(A,A)') "\*Require Input: Mesh Filename 1, Boundary Filename 1, Polynomial Order                             */"
		write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
!
!       Describe Input
!       --------------
		write(STD_OUT,'(/)')
		write(STD_OUT,'(10X,A,A)') "Input Control File:"
		write(STD_OUT,'(10X,A,A)') "-------------------"
		write(STD_OUT,'(30X,A,A25,A30)') "->","Task: ", "horsesMesh2OF"
		write(STD_OUT,'(30X,A,A25,A30)') "->","Mesh Filename 1: ", trim(meshName)
		write(STD_OUT,'(30X,A,A25,A30)') "->","Boundary Filename 1: ", trim(boundaryFile)
		write(STD_OUT,'(30X,A,A25,I5,I5,I5)') "->","Polynomial Mesh: ", Nout(1), Nout(2), Nout(3)	
		write(STD_OUT,'(30X,A,A25,A30)') "->","Discretization nodes: ", "Gauss-Lobatto"		
		boundaryFileName=boundaryFile
		hasBoundaries=.true.		
!
!        Read the mesh
!        -------------
		 call mesh % ReadMesh(meshName)
!
!        Create GaussLobatto Nodes-Nout order
!        ------------------------------------		 
		 call addNewSpectralBasis(spA, Nout, GAUSSLOBATTO)
		 xi   = spA(Nout(1))% x
		 eta  = spA(Nout(2)) % x
		 zeta = spA(Nout(3)) % x
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
			
            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)     
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)       
!
!           Perform interpolation
!           ---------------------

            call ProjectStoragePoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T)
            end associate
         end do
		 
!        Create the result directory
!        ---------------------------
		 CALL getcwd(dir)
		 CALL system('mkdir foamFiles')
		 CALL chdir('foamFiles')
		 
!        Generate Mesh
!        -------------
		 CALL system('mkdir constant')
	     CALL chdir('constant')
	     CALL system('mkdir polyMesh')
		 CALL chdir('polyMesh')
		 CALL createFoamMesh (mesh)
		 CALL chdir(trim(dir))

!        Create pointer for paraview
!        ---------------------------
		 CALL chdir('foamFiles')
		 open(fid, file="paraView.foam", status="unknown", action="write")
		 close(fid)
		 
		 write(STD_OUT,'(/)')
		 write(STD_OUT,'(10X,A,A)') "Finish - horsesMesh2OF"
		 write(STD_OUT,'(10X,A,A)') "----------------------"

        END SUBROUTINE convertHorses2OFMesh
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine convert horses3d mesh type variables into foamMesh using foamMeshStorage Module
!      Data are then written into files
!
        SUBROUTINE createFoamMesh (mesh)
            USE foamMeshStorageConverter
            IMPLICIT NONE
            TYPE(Mesh_t)                 ,INTENT(INOUT)              :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
            type(foamMesh_t)              			   :: foamMesh  
            INTEGER                                    :: fid
            INTEGER                                    :: nFace

!        Construct foamMesh Variable
            CALL foamMesh % Construct(mesh)
			
			write(STD_OUT,'(10X,A,A)') "Creating Output Mesh File:"
            write(STD_OUT,'(10X,A,A)') "-------------------------"
!
!------------------
!        File: points
!------------------
!           Create the file
            CALL createFilePointsHeader (fid)
!           Write elements
            CALL writeFoamPoints (fid,mesh,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "points"
!
!------------------
!        File: faces
!------------------
!           Create the file
            CALL createFileFacesHeader (fid)
!           Write faces
            CALL writeFoamFaces (fid,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "faces"
!
!------------------
!        File: neighbour
!------------------
!           Create the file
            CALL createFileNeighbourHeader (fid)
!           Write faces´s neighbour cell
            CALL writeFoamNeighbour (fid,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "neighbour"
!
!------------------
!        File: owner
!------------------
!           Create the file
            CALL createFileOwnerHeader (fid)
!           Write faces´s owner cell
            CALL writeFoamOwner (fid,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "owner"
!
!------------------
!        File: boundary
!------------------
!           Create the file
            CALL createFileBoundaryHeader (fid)
!           Write boundary
            CALL writeFoamBoundary (fid,foamMesh, mesh)
!           Close the file
            CLOSE(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "boundary"
!
!------------------
!        Foam Mesh Statistics
!------------------
            nFace=foamMesh % nFacesShared+foamMesh % nFacesUnshared+foamMesh % nFacesBoundaries - foamMesh % nMultipleFaces
			write(STD_OUT,*)
            write(STD_OUT,'(10X,A,A)') "Output Mesh Statistics:"
            write(STD_OUT,'(10X,A,A)') "----------------------"
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Points:", foamMesh % nPointsUnique
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Faces:", nFace
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Cells:", foamMesh % nCells
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Boundaries:", size(foamMesh % boundaries)
			write(STD_OUT,*)
!
!        Deallocate

        END SUBROUTINE createFoamMesh
!
!////////////////////////////////////////////////////////////////////////
!
!
      character(len=LINE_LENGTH) function getFormat(nData)
         implicit none
         INTEGER            ,INTENT(IN) :: nData

         INTEGER        :: i

         getFormat = '(A,'

         DO i=1,nData
            write(getFormat,'(A,A)')trim(getFormat),'E18.10,'
         END DO
         write(getFormat,'(A,A)')trim(getFormat),'A)'

      end function getFormat
END MODULE convertMesh2OF
