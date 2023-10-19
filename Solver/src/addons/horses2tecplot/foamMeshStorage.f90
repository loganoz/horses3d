!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!   HORSES3D to Foam Result - foamMeshStorage Module
!
!      This module convert Horses3D mesh storage to foamMesh Storage
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
MODULE foamMeshStorage
    USE SMConstants
    USE Storage
    IMPLICIT NONE

    TYPE foamFaces_t
        INTEGER                     :: facePoints    (4)
        INTEGER                     :: faceOwner
        INTEGER                     :: faceNeighbour
        INTEGER                     :: faceBoundary
    END TYPE foamFaces_t
    !
    TYPE foamBoundaries_t
        INTEGER                     :: faceStart
        INTEGER                     :: nFace
        TYPE(foamFaces_t)           ,ALLOCATABLE    :: faceBoundaries(:)
    END TYPE foamBoundaries_t
    !
    TYPE foamMesh_t
        INTEGER                     :: nCells
        INTEGER                     :: Nout(3)
        INTEGER                     :: nPoints
        INTEGER                     :: nFacesShared
        INTEGER                     :: nFacesUnshared
        INTEGER                     :: nFacesBoundaries
        TYPE(foamFaces_t)           ,ALLOCATABLE    :: faceUnshared(:)
        TYPE(foamFaces_t)           ,ALLOCATABLE    :: faceShared(:)
        TYPE(foamBoundaries_t)      ,ALLOCATABLE    :: boundaries(:)
        CONTAINS
            PROCEDURE   :: Construct     => constructMeshFaces
    END TYPE foamMesh_t
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine construct faces for foam Mesh format with derived type: foamMesh_t with input type Mesh_t
!
        SUBROUTINE constructMeshFaces(self, mesh)
            IMPLICIT NONE
            CLASS(foamMesh_t)                       :: self
            TYPE(Mesh_t),       INTENT(INOUT)       :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: eID              ! Spectral (not finite) element counter
            INTEGER         :: i,j,k
            INTEGER         :: pointID          ! Finite Points counter
            INTEGER         :: N(3)
            INTEGER         :: faceUnsharedCount, faceSharedCount, faceBoundaryCount
            INTEGER         :: nFacesBoundary, nFacesBoundaries
            INTEGER , ALLOCATABLE, DIMENSION(:,:,:)    :: CellPoints
            INTEGER , ALLOCATABLE, DIMENSION(:,:)      :: boundarySet
            INTEGER , ALLOCATABLE, DIMENSION(:)        :: boundaryFace

            self % Nout=mesh % elements(1) % Nout
            N = self % Nout
            self % nCells   = mesh % no_of_elements * ((N(1))*(N(2))*(N(3)))
            self % nPoints  = mesh % no_of_elements * ((N(1)+1)*(N(2)+1)*(N(3)+1))
            self % nFacesShared = mesh % no_of_elements * ((N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))

            faceSharedCount  =0
            faceUnsharedCount=0
            faceBoundaryCount=0
            nFacesBoundaries =0
!
!        Calculate the number of FE Faces on each Boundary -- Tested and Reliable
!        ------------------------------------------------------------------
           ASSOCIATE ( e => mesh % boundaries)
           ALLOCATE( self % boundaries(size(e)))
           DO i=1, size(e)
                nFacesBoundary=0
                DO j=1, e(i) % no_of_faces
                    SELECT CASE(e(i) % elementSides(j))
                        CASE(1,2)
                            nFacesBoundary = nFacesBoundary + N(1)*N(3)
                        CASE(3,5)
                            nFacesBoundary = nFacesBoundary + N(1)*N(2)
                        CASE(4,6)
                            nFacesBoundary = nFacesBoundary + N(2)*N(3)
                    END SELECT
                 END DO
                 nFacesBoundaries=nFacesBoundaries+nFacesBoundary
                 ALLOCATE(self % boundaries(i) % faceBoundaries (nFacesBoundary))
            END DO
            END ASSOCIATE
            self % nFacesBoundaries = nFacesBoundaries
            self % boundaries(:) % nFace = 0
!
!        Store the number of Faces Unshared without boundary - Tested
!        --------------------------------------------
             self % nFacesUnshared   = (mesh % no_of_elements *2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1)))-nFacesBoundaries
!
!        Allocate point, faceUnshared and face shared
!        --------------------------------------------
            ALLOCATE (self % faceUnshared (1:self % nFacesUnshared),self % faceShared (1:self % nFacesShared))
            ALLOCATE(CellPoints(0:N(1),0:N(2),0:N(3)))
!
!        Construct Spectral Boundaries Array -- Tested and Reliable
!        --------------------------------------------
            CALL constructSpectralBoundariesArray(mesh, boundarySet)
!
!        Assign FE faces for shared and unshared without boundary  to foamMesh -- Tested
!        --------------------------------------------
            DO eID = 1 , size(mesh % elements)

                associate ( e => mesh % elements(eID) )
                pointID=((N(1)+1)*(N(2)+1)*(N(3)+1))*(eID-1)           ! FE Point ID after eID
!
!               Assign points at Spectral Element eID for FE faces extraction and store Point into foamMesh  -- Tested
!               -------------------------------------------------------------------------------------------
                DO k = 0, N(3) ; DO j = 0, N(2) ; DO i = 0, N(1)
                    CellPoints(i,j,k)=pointID
                    pointID=pointID+1
                END DO                ; END DO                ; END DO
!
!               Assign Boundary Faces from single Spectral Element eID into Finite Elements -- Tested
!               -------------------------------------------------------------------------------------------
                CALL finiteElementBoundariesFace (N, eID, size(boundarySet,2), boundarySet, boundaryFace)
!
!               Assign Faces from single Spectral Element eID into Finite Elements -- Tested
!               -------------------------------------------------------------------------------------------
                CALL finiteElementFacesSet(CellPoints, N, eID, boundaryFace, faceSharedCount, faceUnsharedCount, self)
                END ASSOCIATE
            END DO
!
!           Assign boundary faces of all FE cell into foamMesh based on its boundary type -- Tested
!           --------------------------------------------
            DO i=1,size(self % boundaries)
                k=0
                IF (i.EQ.1) THEN
                    self % boundaries(i) % faceStart = size(self % faceShared) + size(self % faceUnshared)
                ELSE
                    self % boundaries(i) % faceStart = self % boundaries(i-1) % faceStart +self % boundaries(i-1) % nFace
                END IF
            END DO

        END SUBROUTINE constructMeshFaces
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine handle one spectral element´s nodes and construct the foamMesh faces
!
!
        SUBROUTINE finiteElementFacesSet(CellPoints,N,eID,boundaryFace, faceSharedCount, faceUnsharedCount, self)
            IMPLICIT NONE
            INTEGER                                         ,INTENT(IN)   :: N(3)
            INTEGER     ,DIMENSION(0:N(1),0:N(2),0:N(3))    ,INTENT(IN)   :: CellPoints
            INTEGER                                         ,INTENT(IN)   :: eID
            INTEGER                                         ,INTENT(INOUT):: faceUnsharedCount
            INTEGER                                         ,INTENT(INOUT):: faceSharedCount
            INTEGER     ,DIMENSION(1:6*N(1)*N(2)*N(3))      ,INTENT(IN)   :: boundaryFace
            TYPE(foamMesh_t)                                ,INTENT(INOUT):: self
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k
            INTEGER         :: feID,feIDstart,faceSharedStart, faceUnsharedStart, nUnsharedCount
            INTEGER         :: nStart, nEnd
            INTEGER         ,DIMENSION(1:2,1:2,1:2)                 ::  PointSet
            INTEGER         ,DIMENSION(1:7,1:6*N(1)*N(2)*N(3))      ::  FaceSet
            INTEGER         ,DIMENSION(1:7,1:(N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))      &
                                          ::  FaceSetShared
            INTEGER         ,DIMENSION(1:7,1:(2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1))))  &
                                          :: FaceSetUnshared


            feID              =(eID-1)*N(1)*N(2)*N(3)
            feIDstart         = feID
            faceSharedStart   = ((N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))*(eID-1)
            faceUnsharedStart = 2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1))*(eID-1)

            DO k=0,N(3)-1
                DO j=0,N(2)-1
                    DO i=0,N(1)-1
                        feID    =feID+1
                        nStart  =(feID-1)*6+1-feIDstart*6
                        nEnd    =(feID-1)*6+6-feIDstart*6
                        PointSet(:,:,:)=CellPoints(i:i+1,j:j+1,k:k+1)
                        FaceSet(1:4,nStart:nEnd)=facePointsOrder(PointSet)
                        FaceSet(5,nStart:nEnd)=feID-1
                        FaceSet(6,nStart:nEnd)=0
                        FaceSet(7,nStart:nEnd)=0
                    END DO
                END DO
            END DO

!
!        Assign Boundary ID to FE faces
!        --------------------------------------------
            FaceSet(7,:)=boundaryFace(:)
!
!        Sort array from min to max at 1st point
!        --------------------------------------------
            CALL sortItegerArrayMinMax(SHAPE(FaceSet), FaceSet, 1, 0)
!
!        Search for double faces, assign neighbour, and construct the FaceSetShared
!        --------------------------------------------
            CALL assignNeighbour(N,SHAPE(FaceSet) , FaceSet, FaceSetShared, FaceSetUnshared)
!
!        Sort array from min to max at boundary ID
!        --------------------------------------------
            CALL sortItegerArrayMinMax(SHAPE(FaceSetUnshared), FaceSetUnshared, 7, 0)
!
!
!        --------------------------------------------
            nUnsharedCount=0
            DO i=1, 2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1))
                IF (FaceSetUnshared(7,i).EQ.0) THEN
                    nUnsharedCount=nUnsharedCount+1
                ELSE IF (FaceSetUnshared(7,i).GT.0) THEN
                    ASSOCIATE(e => self % boundaries(FaceSetUnshared(7,i)))
                        e % nFace = e % nFace + 1
                        e % faceBoundaries (e % nFace) % facePoints     =  FaceSetUnshared(1:4,i)
                        e % faceBoundaries (e % nFace) % faceOwner      =  FaceSetUnshared(5,i)
                        e % faceBoundaries (e % nFace) % faceNeighbour  =  FaceSetUnshared(6,i)
                        e % faceBoundaries (e % nFace) % faceBoundary   =  FaceSetUnshared(7,i)
                    END ASSOCIATE
                END IF
            END DO
!
!        Assign Shared Faces
!        --------------------------------------------
            DO i=1,((N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))
                faceSharedCount=faceSharedCount+1
                self % faceShared(faceSharedCount) % facePoints     = FaceSetShared(1:4,i)
                self % faceShared(faceSharedCount) % faceOwner      = FaceSetShared(5,i)
                self % faceShared(faceSharedCount) % faceNeighbour  = FaceSetShared(6,i)
                self % faceShared(faceSharedCount) % faceBoundary   = FaceSetShared(7,i)
            END DO
!
!        Assign Unshared Faces
!        --------------------------------------------
            DO i=1,nUnsharedCount
                faceUnsharedCount=faceUnsharedCount+1
                self % faceUnshared(faceUnsharedCount) % facePoints     = FaceSetUnshared(1:4,i)
                self % faceUnshared(faceUnsharedCount) % faceOwner      = FaceSetUnshared(5,i)
                self % faceUnshared(faceUnsharedCount) % faceNeighbour  = FaceSetUnshared(6,i)
                self % faceUnshared(faceUnsharedCount) % faceBoundary   = FaceSetUnshared(7,i)
            END DO

        END SUBROUTINE finiteElementFacesSet
!
!////////////////////////////////////////////////////////////////////////
!
! Given a set of points ID of a Hexa element, this function return the 6 set of faces construct from Point ID
! Point order pointing to the element center
!
        FUNCTION facePointsOrder(CellPoints) RESULT (FaceSet)
            IMPLICIT NONE
            INTEGER   ,DIMENSION(1:2,1:2,1:2)   ,INTENT(IN)     ::  CellPoints
            INTEGER   ,DIMENSION(1:4,1:6)                       ::  FaceSet
!
!        ---------------
!        Local variables
!        ---------------
!
!           Assign Face : Order follows horses .bmesh file
            FaceSet(1,1)=CellPoints(1,1,1)
            FaceSet(2,1)=CellPoints(1,1,2)
            FaceSet(3,1)=CellPoints(2,1,2)
            FaceSet(4,1)=CellPoints(2,1,1)

            FaceSet(1,2)=CellPoints(1,2,1)
            FaceSet(2,2)=CellPoints(2,2,1)
            FaceSet(3,2)=CellPoints(2,2,2)
            FaceSet(4,2)=CellPoints(1,2,2)

            FaceSet(1,3)=CellPoints(1,1,1)
            FaceSet(2,3)=CellPoints(2,1,1)
            FaceSet(3,3)=CellPoints(2,2,1)
            FaceSet(4,3)=CellPoints(1,2,1)

            FaceSet(1,4)=CellPoints(2,1,1)
            FaceSet(2,4)=CellPoints(2,1,2)
            FaceSet(3,4)=CellPoints(2,2,2)
            FaceSet(4,4)=CellPoints(2,2,1)

            FaceSet(1,5)=CellPoints(1,1,2)
            FaceSet(2,5)=CellPoints(1,2,2)
            FaceSet(3,5)=CellPoints(2,2,2)
            FaceSet(4,5)=CellPoints(2,1,2)

            FaceSet(1,6)=CellPoints(1,1,1)
            FaceSet(2,6)=CellPoints(1,2,1)
            FaceSet(3,6)=CellPoints(1,2,2)
            FaceSet(4,6)=CellPoints(1,1,2)

        END FUNCTION facePointsOrder
!
!////////////////////////////////////////////////////////////////////////
!
! This Recursive Subroutine sort Integer ArraySet at nColumn from its Minimum to Maximum Value
! When call for the first time k must be set to 0 -- WARNING DO NOT USE FOR HUGE ARRAY
!
        RECURSIVE SUBROUTINE sortItegerArrayMinMax(sizeArray, ArraySet, nColumn, k)
            IMPLICIT NONE
            INTEGER         ,DIMENSION(2)                               ,INTENT(IN)       :: sizeArray
            INTEGER         ,DIMENSION(1:sizeArray(1),1:sizeArray(2))   ,INTENT(INOUT)    :: ArraySet
            INTEGER                                                     ,INTENT(IN)       :: nColumn
            INTEGER                                                     ,INTENT(IN)       :: k
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         ,DIMENSION(1:sizeArray(1))  :: BufferSet
            INTEGER                                     :: i
!
!        ---------------
!        Perform Sorting Value
!        ---------------
!
            SELECT CASE(k)
            CASE(0)
                DO i=1,sizeArray(2)-1
                    IF (ArraySet(nColumn,i).GT.ArraySet(nColumn,i+1)) THEN
                        BufferSet(:)    = ArraySet(:,i)
                        ArraySet(:,i)   = ArraySet(:,i+1)
                        ArraySet(:,i+1) = BufferSet(:)
                        IF(i.GT.1) THEN
                            CALL sortItegerArrayMinMax(sizeArray, ArraySet, nColumn, i)
                        END IF
                    END IF
                END DO
            CASE DEFAULT ! For Recursive
                IF (ArraySet(nColumn,k-1).GT.ArraySet(nColumn,k)) THEN
                    BufferSet(:)    = ArraySet(:,k-1)
                    ArraySet(:,k-1) = ArraySet(:,k)
                    ArraySet(:,k)   = BufferSet(:)
                    IF(k.GT.1) THEN
                        CALL sortItegerArrayMinMax(sizeArray, ArraySet, nColumn, k-1)
                    END IF
                END IF
            END SELECT

        END SUBROUTINE sortItegerArrayMinMax
!
!////////////////////////////////////////////////////////////////////////
!
! This Subroutine assign Neighbour cell to the 1st face which are used by 2 cell
!
        SUBROUTINE assignNeighbour(N,sizeArray, FaceSet, FaceSetShared, FaceSetUnshared)
            IMPLICIT NONE
            INTEGER                                                     ,INTENT(IN)     :: N(3)
            INTEGER         ,DIMENSION(2)                               ,INTENT(IN)     :: sizeArray
            INTEGER         ,DIMENSION(1:sizeArray(1),1:sizeArray(2))   ,INTENT(INOUT)  :: FaceSet
            INTEGER         ,DIMENSION(1:sizeArray(1),1:(N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))  &
                    ,INTENT(OUT)    :: FaceSetShared
            INTEGER         ,DIMENSION(1:sizeArray(1),1:2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1)))  &
                    ,INTENT(OUT)    :: FaceSetUnshared
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k
            k=0
!
!        ---------------
!        Search for double face (1st and 3rd Point Identical) -> Assign neighbour cell to both
!        ---------------
!
            DO i=1,sizeArray(2)-5
                DO j=1,5
                    IF ((FaceSet(1,i).EQ.FaceSet(1,i+j)).AND.(FaceSet(3,i).EQ.FaceSet(3,i+j))) THEN
                        k=k+1
                        FaceSet(6,i)=FaceSet(5,i+j)             ! Assign Neighbour both
                        FaceSet(6,i+j)=1                        ! Assign 1 to prevent error on unshared face filter
                        FaceSetShared(:,k)=FaceSet(:,i)         ! Store the first
                        EXIT
                    END IF
                END DO
            END DO

            k=0
            DO i=1,sizeArray(2)
                IF (FaceSet(6,i).EQ.0) THEN
                    k=k+1
                    FaceSetUnshared(:,k)=FaceSet(:,i)
                END IF
            END DO
        END SUBROUTINE assignNeighbour
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
! This Subroutine rearrange Boundary Set from spectra mesh boundary into array [CellID, elementSide, boundaryID]
! The array is sorted from low to high Cell ID
!
        SUBROUTINE constructSpectralBoundariesArray(mesh, boundarySet)
            IMPLICIT NONE
            TYPE(Mesh_t)                                        ,INTENT(IN)     :: mesh
            INTEGER   , ALLOCATABLE, DIMENSION(:,:)  ,INTENT(OUT)    :: boundarySet
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k=0
            INTEGER         :: nFacesTotal
!
!        ---------------
!        Construct Array for Spectral Boundary and sort the Cell from low to high value
!        ---------------
!
            nFacesTotal=SUM(mesh % boundaries (:) % no_of_faces)
            ALLOCATE(boundarySet(1:3,1:nFacesTotal))

            DO i=1, size(mesh % boundaries)
                DO j=1, mesh % boundaries(i) % no_of_faces
                    k=k+1
                    boundarySet(1,k)= (mesh % boundaries(i) % elements(j)) - 1  ! foam Cell start from 0
                    boundarySet(2,k)= mesh % boundaries(i) % elementSides(j)    ! face Side
                    boundarySet(3,k)= i                                         ! boundaryID
                END DO
            END DO

!        Sorting Cell Value
!        ---------------
            CALL sortItegerArrayMinMax(SHAPE(boundarySet), boundarySet, 1, 0)

        END SUBROUTINE constructSpectralBoundariesArray
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
! This Subroutine rearrange Boundary Set from spectra mesh boundary into array [CellID, elementSide, boundaryID]
! The array is sorted from low to high Cell ID
!
        SUBROUTINE finiteElementBoundariesFace (N, eID, sizeBoundary, boundarySet, boundaryFace)
            IMPLICIT NONE

            INTEGER                                  ,INTENT(IN)    :: N(3)
            INTEGER                                  ,INTENT(IN)    :: eID
            INTEGER                                  ,INTENT(IN)    :: sizeBoundary
            INTEGER   ,DIMENSION(3,sizeBoundary)     ,INTENT(IN)    :: boundarySet
            INTEGER   ,ALLOCATABLE  , DIMENSION(:)   ,INTENT(OUT)   :: boundaryFace
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k
            INTEGER         :: nStart, nEnd

            k=0
            nStart=0
            nEnd=0
!
!           Allocate boundaryFace for each faces in FE (including double faces)
!           -----------------------------------
            ALLOCATE (boundaryFace(1:6*N(1)*N(2)*N(3)))
!
!        ---------------
!        Count Number of Spectra Faces for element eID
!        ---------------
!
            DO i=1,size(boundarySet,2)
                IF((boundarySet(1,i).EQ.(eID-1)).AND.(k.EQ.0)) THEN
                    k=1
                    nStart=i
                    nEnd  =i
                ELSE IF((boundarySet(1,i).EQ.(eID-1)).AND.(k.NE.0)) THEN
                    nEnd  =i
                ELSE IF (boundarySet(1,i).GT.(eID-1)) THEN
                    EXIT
                END IF
            END DO
            k=0

!
!        ---------------
!        Count Number of Spectra Faces for element eID
!        ---------------
!
            boundaryFace(:)=0               ! Fill 0
            DO i=1, nEnd-nStart+1
                SELECT CASE(boundarySet(2,nStart+i-1))
                    CASE(1)
                        DO j=1, N(3)
                            DO k=1,N(1)
                                boundaryFace((((j-1)*N(1)*N(2)+k)-1)*6+1)=boundarySet(3,nStart+i-1)
                            END DO
                        END DO
                    CASE(2)
                        DO j=1, N(3)
                            DO k=1,N(1)
                                boundaryFace((((j)*N(1)*N(2)-N(1)+k)-1)*6+2)=boundarySet(3,nStart+i-1)
                            END DO
                        END DO
                    CASE(3)
                        DO j=1, N(1)*N(2)
                            boundaryFace((j-1)*6+3)=boundarySet(3,nStart+i-1)
                        END DO
                    CASE(4)
                        DO j=N(1), (N(1)*N(2)*N(3)), N(1)
                            boundaryFace((j-1)*6+4)=boundarySet(3,nStart+i-1)
                        END DO
                    CASE(5)
                        DO j=(N(1)*N(2)*N(3))-(N(1)*N(2))+1, (N(1)*N(2)*N(3))
                            boundaryFace((j-1)*6+5)=boundarySet(3,nStart+i-1)
                        END DO
                    CASE(6)
                        DO j=1, (N(1)*N(2)*N(3)), N(1)
                            boundaryFace((j-1)*6+6)=boundarySet(3,nStart+i-1)
                        END DO
                 END SELECT
            END DO

        END SUBROUTINE finiteElementBoundariesFace
!
END MODULE foamMeshStorage
!
!////////////////////////////////////////////// END OF FILE //////////////////////////////////////////////
!
