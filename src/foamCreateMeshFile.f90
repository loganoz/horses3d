!
! ////////////////////////////////////////////////////////////////////
!   HORSES3D to foam - foamCreateMeshFile Module
!
!      This module converts foamMeshStorage data into foamMesh files
!       foamMesh files consist of 'points', 'faces', 'owner', 'neighbour' and 'boundary files
!
!////////////////////////////////////////////////////////////////////////////////////////
!
MODULE foamCreateMeshFile
    USE SMConstants
    USE Storage
    IMPLICIT NONE

!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFilePointsHeader (fid)
            IMPLICIT NONE
            INTEGER                   , INTENT(OUT)     :: fid
!        ---------------
!        Local variables
!        ---------------
!
            CHARACTER(LEN=LINE_LENGTH)     :: fileName, fileClass, fileNote, fileLocation, fileObject
!
            fileName    = 'points'
            fileClass   ='vectorField;'
            fileNote    =''
            fileLocation='"constant/polyMesh";'
            fileObject  ='points;'
!

            CALL createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
        END SUBROUTINE createFilePointsHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFileFacesHeader (fid)
            IMPLICIT NONE
            INTEGER                   , INTENT(OUT)     :: fid
!        ---------------
!        Local variables
!        ---------------
!
            CHARACTER(LEN=LINE_LENGTH)     :: fileName, fileClass, fileNote, fileLocation, fileObject
!
            fileName    = 'faces'
            fileClass   ='faceList;'
            fileNote    =''
            fileLocation='"constant/polyMesh";'
            fileObject  ='faces;'
!

            CALL createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
        END SUBROUTINE createFileFacesHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFileNeighbourHeader (fid)
            IMPLICIT NONE
            INTEGER                   , INTENT(OUT)     :: fid
!        ---------------
!        Local variables
!        ---------------
!
            CHARACTER(LEN=LINE_LENGTH)     :: fileName, fileClass, fileNote, fileLocation, fileObject
!
            fileName    = 'neighbour'
            fileClass   ='labelList;'
            fileNote    =''
            fileLocation='"constant/polyMesh";'
            fileObject  ='neighbour;'
!

            CALL createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
        END SUBROUTINE createFileNeighbourHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFileOwnerHeader (fid)
            IMPLICIT NONE
            INTEGER                   , INTENT(OUT)     :: fid
!        ---------------
!        Local variables
!        ---------------
!
            CHARACTER(LEN=LINE_LENGTH)     :: fileName, fileClass, fileNote, fileLocation, fileObject
!
            fileName    = 'owner'
            fileClass   ='labelList;'
            fileNote    =''
            fileLocation='"constant/polyMesh";'
            fileObject  ='owner;'
!
            CALL createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
        END SUBROUTINE createFileOwnerHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFileBoundaryHeader (fid)
            IMPLICIT NONE
            INTEGER                   , INTENT(OUT)     :: fid
!        ---------------
!        Local variables
!        ---------------
!
            CHARACTER(LEN=LINE_LENGTH)     :: fileName, fileClass, fileNote, fileLocation, fileObject
!
            fileName    = 'boundary'
            fileClass   ='polyBoundaryMesh;'
            fileNote    =''
            fileLocation='"constant/polyMesh";'
            fileObject  ='boundary;'
!
            CALL createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
        END SUBROUTINE createFileBoundaryHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFileResultHeader (name, nData, fid)
            IMPLICIT NONE
            CHARACTER(LEN=LINE_LENGTH), INTENT(IN)      :: name
            INTEGER                   , INTENT(IN)      :: nData
            INTEGER                   , INTENT(OUT)     :: fid
!        ---------------
!        Local variables
!        ---------------
!
            CHARACTER(LEN=LINE_LENGTH)     :: fileName, fileClass, fileNote, fileLocation, fileObject
!
            fileName    = trim(name)
            IF (nData.EQ.1) THEN
                fileClass   ='pointScalarField;'
            ELSE IF (nData.EQ.3) THEN
                fileClass   ='pointVectorField;'
            ELSE
                fileClass   ='pointTensorField;'
            END IF

            fileNote    =''
            fileLocation='"0";'
            write(fileObject,'(A,A)')trim(name),';'
!

            CALL createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
        END SUBROUTINE createFileResultHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE createFileHeader (fileName, fileClass, fileNote, fileLocation, fileObject, fid)
            IMPLICIT NONE
            CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: fileName, fileClass, fileNote, fileLocation, fileObject
            INTEGER                   , INTENT(OUT)    :: fid

!        ---------------
!        Local variables
!        ---------------
!
!
            open(newunit = fid, file =trim(fileName), action = "write", status = "unknown")

            write(fid,'(A)')"/*-------------------------------------- HORSES3D - FOAM FILE --------------------------------------*\"
            write(fid,'(A)')"####################################################################################################"
            write(fid,'(A)')"#                                                                                                  #"
            write(fid,'(A)')"#            HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver             #"
            write(fid,'(A)')"#                              Foam File for Paraview Post-Processing                              #"
            write(fid,'(A)')"#                                                                                                  #"
            write(fid,'(A)')"####################################################################################################"
            write(fid,'(A)')"\*--------------------------------------------------------------------------------------------------*/"
            write(fid,'(A)')"FoamFile"
            write(fid,'(A)')"{"
            write(fid,'(A)')'    format      ascii;'
            write(fid,'(A,A)')'    class       ',trim(fileClass)
            IF (LEN(trim(fileNote)).GT.1) THEN
                write(fid,'(A,A)')'    note        ',trim(fileNote)
            END IF
            write(fid,'(A,A)')'    location    ',trim(fileLocation)
            write(fid,'(A,A)')'    object      ',trim(fileObject)
            write(fid,'(A)')'}'
            write(fid,'(A)')"//--------------------------------------------------------------------------------------------------//"
            write(fid,*)
            write(fid,*)

        END SUBROUTINE createFileHeader
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE writeFoamPoints (fid,mesh,foamMesh)
            USE foamMeshStorage
            IMPLICIT NONE
            INTEGER                     , INTENT(IN)        :: fid
            TYPE(Mesh_t )               , INTENT(IN)        :: mesh
            TYPE(foamMesh_t )           , INTENT(IN)        :: foamMesh
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k, eID

            write(fid,'(I0)') foamMesh % nPoints
            write(fid,'(A)')'('

!        Write the points
!        ----------------
            do eID = 1, size(mesh % elements)
                associate ( e => mesh % elements(eID) )
                do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                    write(fid,'(A,E18.10,E18.10,E18.10,A)')'( ', e % xOut(:,i,j,k), ' )'
                end do                ; end do                ; end do
                end associate
            end do
         write(fid,'(A)')')'

        END SUBROUTINE writeFoamPoints
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE writeFoamFaces (fid,foamMesh)
            USE foamMeshStorage
            IMPLICIT NONE
            INTEGER                     , INTENT(IN)        :: fid
            TYPE(foamMesh_t )           , INTENT(IN)        :: foamMesh
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j

            write (fid,'(I0)')foamMesh % nFacesShared+foamMesh % nFacesUnshared+foamMesh % nFacesBoundaries
            write(fid,'(A)')'('
            DO i=1,foamMesh % nFacesShared
                write(fid,'(A,I10,I10,I10,I10,A)')'4 ( ', foamMesh % faceShared (i) % facePoints(:), ' )'
            END DO
            DO i=1,foamMesh % nFacesUnshared
                write(fid,'(A,I10,I10,I10,I10,A)')'4 ( ', foamMesh % faceUnshared (i) % facePoints(:), ' )'
            END DO
            DO i=1,size(foamMesh % boundaries)
                ASSOCIATE (e=> foamMesh % boundaries(i))
                DO j=1, size (e % faceBoundaries)
                    write(fid,'(A,I10,I10,I10,I10,A)')'4 ( ',e % faceBoundaries(j) % facePoints(:), ' )'
                END DO
                END ASSOCIATE
            END DO
            write(fid,'(A)')')'
        END SUBROUTINE writeFoamFaces
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE writeFoamOwner (fid,foamMesh)
            USE foamMeshStorage
            IMPLICIT NONE
            INTEGER                     , INTENT(IN)        :: fid
            TYPE(foamMesh_t )           , INTENT(IN)        :: foamMesh
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j

            write (fid,'(I0)')foamMesh % nFacesShared+foamMesh % nFacesUnshared+foamMesh % nFacesBoundaries
            write(fid,'(A)')'('
            DO i=1,foamMesh % nFacesShared
                write(fid,'(I0)') foamMesh % faceShared (i)     % faceOwner
            END DO
            DO i=1,foamMesh % nFacesUnshared
                write(fid,'(I0)') foamMesh % faceUnshared (i)   % faceOwner
            END DO
            DO i=1,size(foamMesh % boundaries)
                ASSOCIATE (e=> foamMesh % boundaries(i))
                DO j=1, size (e % faceBoundaries)
                    write(fid,'(I0)') e % faceBoundaries(j) % faceOwner
                END DO
                END ASSOCIATE
            END DO
            write(fid,'(A)')')'
        END SUBROUTINE writeFoamOwner
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE writeFoamNeighbour (fid,foamMesh)
            USE foamMeshStorage
            IMPLICIT NONE
            INTEGER                     , INTENT(IN)        :: fid
            TYPE(foamMesh_t )           , INTENT(IN)        :: foamMesh
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i
            write (fid,'(I0)')foamMesh % nFacesShared
            write(fid,'(A)')'('
            DO i=1,foamMesh % nFacesShared
                write(fid,'(I0)') foamMesh % faceShared (i) % faceNeighbour
            END DO
            write(fid,'(A)')')'
        END SUBROUTINE writeFoamNeighbour
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE writeFoamBoundary (fid,foamMesh, mesh)
            USE foamMeshStorage
            IMPLICIT NONE
            INTEGER                     , INTENT(IN)        :: fid
            TYPE(foamMesh_t )           , INTENT(IN)        :: foamMesh
            TYPE(Mesh_t )               , INTENT(IN)        :: mesh
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i

            write (fid,'(I0)')size(foamMesh % boundaries)
            write(fid,'(A)')'('
            DO i=1,size(foamMesh % boundaries)
                write(fid,'(A,A)')'    ',mesh % boundaries (i) % Name
                write(fid,'(A)') '    {'
                write(fid,'(A)') '        type            patch;'
                write(fid,'(A,I0,A)') '        nFaces          ', foamMesh % boundaries(i) % nFace,';'
                write(fid,'(A,I0,A)') '        startFace       ', foamMesh % boundaries(i) % faceStart,';'
                write(fid,'(A)') '    }'
            END DO
            write(fid,'(A)')')'
            write(fid,'(A)')''
            write(fid,'(A)')"//------------------------------------------ END OF FILE -------------------------------------------//"
        END SUBROUTINE writeFoamBoundary

END MODULE foamCreateMeshFile
