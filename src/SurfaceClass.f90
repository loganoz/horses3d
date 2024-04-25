#include "Includes.h"
Module SurfaceClass 
    use SMConstants
    use HexMeshClass
    use ZoneClass
    use ElementConnectivityDefinitions, only: NODES_PER_FACE, FACES_PER_ELEMENT
    Implicit None

    private
    public Surface_t, SurfaceEdge, SurfaceFace_t, SurfaceElement_t

    type SurfaceEdge

        real(kind=RP), dimension(3,2)                       :: corners

        contains

            procedure :: construct          => EdgeConstruct
            procedure :: isEqual            => EdgeIsEqualToEdge
            procedure :: isConnected        => EdgeIsConnectedToEdge
            procedure :: getBCPostion       => EdgeGetBCPosition

    end type SurfaceEdge

    type SurfaceFace_t

        class(SurfaceEdge), dimension(:), allocatable       :: edges
        integer                                             :: eID
        integer                                             :: globaleID
        integer                                             :: fID
        integer                                             :: numberOfConnections

        contains

            procedure :: construct          => FaceConstruct
            procedure :: destruct           => FaceDestruct 
            procedure :: setNoConnections   => FaceSetNoOfConnections
            procedure :: isConnected        => FaceIsConnected
            procedure :: isFullConnected    => FaceIsFullConnected
            procedure :: isTwiceEdConnected => FaceIsConnectedByEdgeTwice
            procedure :: getBCPostion       => FaceGetBCPosition
            procedure :: shareCorner        => FaceShareCorner 
            procedure :: reconstructPeriod  => FaceReconstructPeriodic 
            procedure :: updateEdgesPeriod  => FaceUpdateEdgesOnPeriodic

    end type SurfaceFace_t

    type SurfaceElement_t

        class(SurfaceFace_t), dimension(:), allocatable     :: faces
        integer                                             :: eID
        integer                                             :: globaleID
        integer                                             :: fID
        integer, dimension(:), allocatable                  :: extrafIDs
        logical                                             :: needSecondFace
        logical                                             :: isInBCZone

        contains

            procedure :: construct          => ElementConstruct
            procedure :: destruct           => ElementDestruct 
            procedure :: updateIsInZone     => ElementUpdateIsInZone
            procedure :: setNeedSecond      => ElementSetNeedSecond
            procedure :: setNeedNotSecond   => ElementSetNotNeedSecond
            procedure :: getNotConnectedN   => ElementGetNotConnectedN
            procedure :: reconstructPeriod  => ElementReconstructPeriodic

    end type SurfaceElement_t

    type Surface_t

        integer                                              :: totalNumberOfPoints
        integer                                              :: totalNumberOfFaces
        integer, dimension(:), allocatable                   :: globaleIDs
        integer, dimension(:), allocatable                   :: fIDs
        class(Zone_t), allocatable                           :: surfaceZone
        character(len=LINE_LENGTH)                           :: Name
        character(len=LINE_LENGTH)                           :: fileName

        contains

            procedure :: construct          => SurfaceConstruct
            procedure :: destruct           => SurfaceDestruct
            procedure :: writeToTecplot     => SurfaceWriteSingleZoneToTecplot
            procedure :: saveToFile         => SurfaceSaveToFile

    end type Surface_t

    integer,          parameter                         :: NO_OF_OUTPUTVARIABLES = 2
    character(len=8), parameter                         :: PRECISION_FORMAT = "(E18.10)"

    contains
!
!////////////////////////////////////////////////////////////////////////
! surface class methods
!////////////////////////////////////////////////////////////////////////
!
    Subroutine SurfaceConstruct(self, mesh, fileName, eIDs, gIDs, fIDs, N)

        use FileReadingUtilities, only: getFileName, RemovePath
        use MPI_Utilities
        implicit none

        class(Surface_t)                                :: self
        type(HexMesh), intent(in)                       :: mesh
        character(len=LINE_LENGTH), intent(in)          :: fileName
        integer, dimension(N), intent(in)               :: eIDs, gIDs, fIDs
        integer, intent(in)                             :: N

        ! local variables
        integer                                         :: totalN

        totalN = N

        self % totalNumberOfFaces = totalN
        self % totalNumberOfPoints = NODES_PER_FACE * totalN
        self % fileName = fileName
        self % name = RemovePath(getFileName(fileName))

        allocate(self % surfaceZone)
        call self % surfaceZone % CreateFictitious(-1, self % name, totalN, fIDs)

        allocate( self % globaleIDs(totalN), self % fIDs(totalN) )
        self % globaleIDs = gIDs
        self % fIDs = fIDs

        print *, "surf constructed"

    End Subroutine SurfaceConstruct
! 
!////////////////////////////////////////////////////////////////////////
!
    Subroutine SurfaceDestruct(self)

        implicit none
        class(Surface_t), intent(inout)                   :: self

        safedeallocate(self % globaleIDs)
        safedeallocate(self % fIDs)
        safedeallocate(self % surfaceZone)

    End Subroutine SurfaceDestruct
! 
!////////////////////////////////////////////////////////////////////////
!
    Subroutine SurfaceWriteSingleZoneToTecplot(self, mesh, meshName)

        use FileReadingUtilities, only: getFileName
        implicit none

        class(Surface_t)                                :: self
        type(HexMesh), intent(in)                       :: mesh
        character(len=LINE_LENGTH), intent(in)          :: meshName

        ! local variables
        integer                                         :: fd, fID, i, row
        character(len=LINE_LENGTH)                      :: formatout, title, tecplotFileName
        integer                                         :: Nx,Ny
        integer, dimension(NODES_PER_FACE)              :: corners
        logical                                         :: facesHasErrors
        integer                                         :: totalNumberOfPoints, totalNumberOfFaces, skipped

        totalNumberOfPoints = self % totalNumberOfPoints
        totalNumberOfFaces = self % totalNumberOfFaces
        tecplotFileName = trim(getFileName(self % fileName)) // ".tec"
        open ( newunit = fD , file = trim(tecplotFileName) , status = "unknown" , action = "write" ) 

!       Add the title
!       -------------
        write(title,'(A,A,A)') '"Generated from ',trim(meshName),'"'
        write(fd,'(A,A)') "TITLE = ", trim(title)

!       Add the variables
!       -----------------
        write(fd,'(A,A)') 'VARIABLES = "x","y","z"', '"eID","fID"'

!       Add the zone info
!       -----------------
        write(fd,'(A,I0,A,I0,A,A,A)') "ZONE N=", totalNumberOfPoints,", E=", totalNumberOfFaces, &
                                      ',ET=QUADRILATERAL, F=FEPOINT, T="surface_', trim(self % Name), '"'
                      
!       Write the points
!       ----------------
        write(formatout,'(A,I0,A,A)') "(",3+no_of_outputVariables,PRECISION_FORMAT,")"
        do fID = 1 , self % totalNumberOfFaces
            ! if (self % fIDs(fID) .eq. 0) cycle
            associate ( f => mesh % faces(self % fIDs(fID)) )
                Nx = f % Nf(1)
                Ny = f % Nf(2)
                associate (x => f % geom % x)
                    write(fd,trim(formatout)) x(:,0,0), real(self % globaleIDs(fID),KIND=RP), real(self % fIDs(fID),KIND=RP)
                    write(fd,trim(formatout)) x(:,Nx,0), real(self % globaleIDs(fID),KIND=RP), real(self % fIDs(fID),KIND=RP)
                    write(fd,trim(formatout)) x(:,Nx,Ny), real(self % globaleIDs(fID),KIND=RP), real(self % fIDs(fID),KIND=RP)
                    write(fd,trim(formatout)) x(:,0,Ny), real(self % globaleIDs(fID),KIND=RP), real(self % fIDs(fID),KIND=RP)
                end associate
            end associate
        end do

!       Write the connections
!       ---------------------
        row = 0
        do fID = 1 , self % totalNumberOfFaces
            ! if (self % fIDs(fID) .eq. 0) cycle
            row = row + 1
            corners = [(i,i=1,NODES_PER_FACE)] + NODES_PER_FACE*(row-1)

            write(fd,*) corners
        end do

        close(fd)

    End Subroutine SurfaceWriteSingleZoneToTecplot
! 
!////////////////////////////////////////////////////////////////////////
!
    Subroutine SurfaceSaveToFile(self, mesh, saveNodes)

        implicit none

        class(Surface_t)                                :: self
        type(HexMesh), intent(in)                       :: mesh
        logical, intent(in)                             :: saveNodes

        ! local variables
        integer                                         :: fd, fID

        open ( newunit = fD , file = trim(self % fileName) , status = "unknown" , action = "write" ) 

!       Add total numberOfFaces
!       -------------
        write(fd,*) self % totalNumberOfFaces

        write(fd,*) saveNodes
                      
!       Write the points
!       ----------------
        if (saveNodes) then
            do fID = 1 , self % totalNumberOfFaces
                write(fd,*) self % globaleIDs(fID), self % fIDs(fID), mesh % faces(self % fIDs(fID)) % nodeIDs
            end do
        else
            do fID = 1 , self % totalNumberOfFaces
                write(fd,*) self % globaleIDs(fID), self % fIDs(fID)
            end do
        end if 

        close(fd)

    END SUBROUTINE SurfaceSaveToFile
!
!////////////////////////////////////////////////////////////////////////
! edge methods
!////////////////////////////////////////////////////////////////////////
!
    Subroutine EdgeConstruct(self, f, side)
!        *****************************************************************************
!           This subroutine creates and edge based on the desired side
!           being counterclock numbered
!               ----(3)----
!              |           |
!              |           |
!             (4)         (2)
!              |           |
!              |           |
!               ----(1)----
!        *****************************************************************************

        use FaceClass
        implicit none

        class(SurfaceEdge)                              :: self
        type(Face), intent(in)                          :: f
        integer, intent(in)                             :: side

        !local variables
        integer                                         :: Nx, Ny

        Nx = f % Nf(1)
        Ny = f % Nf(2)

        associate(x => f % geom % x)
            select case (side)
                case (1)
                    self % corners(:,1) = x(:,0,0)
                    self % corners(:,2) = x(:,Nx,0)
                case (2)
                    self % corners(:,1) = x(:,Nx,0)
                    self % corners(:,2) = x(:,Nx,Ny)
                case (3)
                    self % corners(:,1) = x(:,Nx,Ny)
                    self % corners(:,2) = x(:,0,Ny)
                case (4)
                    self % corners(:,1) = x(:,0,Ny)
                    self % corners(:,2) = x(:,0,0)
            end select
        end associate

    End Subroutine EdgeConstruct
!
!////////////////////////////////////////////////////////////////////////
!
    Function EdgeIsEqualToEdge(self, otherEdge) result(isEqual)

        use Utilities, only: AlmostEqual
        implicit none

        class(SurfaceEdge)                              :: self
        class(SurfaceEdge), intent(in)                  :: otherEdge
        logical                                         :: isEqual

        !local variables
        integer                                         :: i, j, equalCorners, usedCorner
        isEqual = .false.
        equalCorners = 0
        usedCorner = 0
        do i = 1, 2
            ext_corner_loop: do j = 1, 2
                if (j .eq. usedCorner) cycle ext_corner_loop
                if (testEqualCorners(self % corners(:,i), otherEdge % corners(:,j))) then
                    equalCorners = equalCorners + 1
                    usedCorner = j
                    exit ext_corner_loop
                end if 
            end do ext_corner_loop
        end do
        if (equalCorners .eq. 2) isEqual = .true.

    End Function EdgeIsEqualToEdge
!
!////////////////////////////////////////////////////////////////////////
!
    Function EdgeIsConnectedToEdge(self, otherEdge) result(isConnected)

        use Utilities, only: AlmostEqual
        implicit none

        class(SurfaceEdge)                              :: self
        class(SurfaceEdge), intent(in)                  :: otherEdge
        logical                                         :: isConnected

        !local variables
        integer                                         :: i, j, equalCorners, usedCorner
        isConnected = .false.
        equalCorners = 0
        usedCorner = 0
        do i = 1, 2
            ext_corner_loop: do j = 1, 2
                if (j .eq. usedCorner) cycle ext_corner_loop
                if (testEqualCorners(self % corners(:,i), otherEdge % corners(:,j))) then
                    equalCorners = equalCorners + 1
                    usedCorner = j
                    exit ext_corner_loop
                end if 
            end do ext_corner_loop
        end do
        if (equalCorners .eq. 1) isConnected = .true.

    End Function EdgeIsConnectedToEdge
!
!////////////////////////////////////////////////////////////////////////
!
    Function EdgeGetBCPosition(self, BCDimension) result(pos)

        class(SurfaceEdge)                              :: self
        integer, intent(in)                             :: BCDimension
        real(kind=RP), dimension(2)                     :: pos

        pos(1) = self % corners(BCDimension,1)
        pos(2) = self % corners(BCDimension,2)

    End Function EdgeGetBCPosition
!
!////////////////////////////////////////////////////////////////////////
! surface face methods
!////////////////////////////////////////////////////////////////////////
!
    Subroutine FaceConstruct(self, mesh, eID, geID, fID)

        implicit none

        class(SurfaceFace_t)                            :: self
        type(HexMesh), intent(in)                       :: mesh
        integer, intent(in)                             :: eID, geID, fID

        !local variables
        integer                                         :: i

        self % eID = eID
        self % globaleID = geID
        self % fID = fID

        allocate(self % edges(NODES_PER_FACE))
        do i = 1, NODES_PER_FACE
            call self % edges(i) % construct(mesh % faces(fID), i)
        end do

        ! This is change only for faces at boundaries, for non closed surfaces
        self % numberOfConnections = NODES_PER_FACE

    End Subroutine FaceConstruct
! 
!////////////////////////////////////////////////////////////////////////
!
    Subroutine FaceDestruct(self)

        implicit none
        class(SurfaceFace_t), intent(inout)             :: self

        safedeallocate(self % edges)

    End Subroutine FaceDestruct
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine FaceSetNoOfConnections(self, mode)

        use Utilities, only: AlmostEqual
        implicit none

        class(SurfaceFace_t)                            :: self
        integer, intent(in)                             :: mode

        select case (mode)
        case (1)
            self % numberOfConnections = self % numberOfConnections - 1
        case (2)
            self % numberOfConnections = self % numberOfConnections - 2
        case default
            self % numberOfConnections = self % numberOfConnections - 1
    end select

    End Subroutine FaceSetNoOfConnections
!
!////////////////////////////////////////////////////////////////////////
!
    Function FaceIsConnected(self,otherFace) result(isCon)

        implicit none

        class(SurfaceFace_t)                            :: self
        class(SurfaceFace_t), intent(in)                :: otherFace
        logical                                         :: isCon

        !local variables
        integer                                         :: k, edgeID

        isCon = .false.
        do edgeID = 1, NODES_PER_FACE
            do k = 1, NODES_PER_FACE
                if (self % edges(edgeID) % isEqual(otherFace % edges(k))) then
                    isCon = .true.
                    return
                end if 
            end do
        end do 

    End Function FaceIsConnected
!
!////////////////////////////////////////////////////////////////////////
!
    Function FaceIsConnectedByEdgeTwice(self, surfaceFaces, N) result(isTwice)

        implicit none

        class(SurfaceFace_t)                            :: self
        type(SurfaceFace_t), dimension(N), intent(in)   :: surfaceFaces
        integer, intent(in)                             :: N
        logical                                         :: isTwice

        !local variables
        integer                                         :: i, k, equalEdges, edgeID
        integer, dimension(NODES_PER_FACE)              :: usedEdges

        equalEdges = 0
        usedEdges = 0

        face_loop: do i = 1, N
            if ((self % fID .eq. surfaceFaces(i) % fID) .and. (self % globaleID .eq. surfaceFaces(i) % globaleID)) cycle face_loop
            if (equalEdges .ge. NODES_PER_FACE) exit face_loop
            this_edge_loop: do edgeID = 1, NODES_PER_FACE
                ext_edge_loop: do k = 1, NODES_PER_FACE
                    if (self % edges(edgeID) % isEqual(surfaceFaces(i) % edges(k))) then
                        equalEdges = equalEdges + 1
                        usedEdges(equalEdges) = edgeID
                        cycle face_loop
                    end if 
                end do ext_edge_loop 
            end do this_edge_loop 
        end do face_loop

        isTwice = .false.
        do i = 1, equalEdges
            do k = i+1, equalEdges
                if (usedEdges(i) .eq. usedEdges(k)) then
                    isTwice = .true.
                    return
                end if
            end do
        end do

    End Function FaceIsConnectedByEdgeTwice
!
!////////////////////////////////////////////////////////////////////////
!
    Function FaceIsFullConnected(self, surfaceFaces, N) result(isCon)

        implicit none

        class(SurfaceFace_t)                            :: self
        class(SurfaceFace_t), dimension(N), intent(in)  :: surfaceFaces
        integer, intent(in)                             :: N
        logical                                         :: isCon

        !local variables
        integer                                         :: i, k, equalEdges, edgeID
        ! integer, dimension(NODES_PER_FACE)              :: usedEdges

        equalEdges = 0
        ! usedEdges = 0

        face_loop: do i = 1, N
            if ((self % fID .eq. surfaceFaces(i) % fID) .and. (self % globaleID .eq. surfaceFaces(i) % globaleID)) cycle face_loop
            this_edge_loop: do edgeID = 1, NODES_PER_FACE
                ! if (any(usedEdges .eq. edgeID)) cycle this_edge_loop
                ext_edge_loop: do k = 1, NODES_PER_FACE
                    if (self % edges(edgeID) % isEqual(surfaceFaces(i) % edges(k))) then
                        equalEdges = equalEdges + 1
                        ! usedEdges(equalEdges) = edgeID
                        cycle face_loop
                    end if 
                end do ext_edge_loop 
            end do this_edge_loop 

        end do face_loop
        isCon = .false.
        if (equalEdges .eq. self % numberOfConnections) isCon = .true.

    End Function FaceIsFullConnected
!
!////////////////////////////////////////////////////////////////////////
!
    Function FaceGetBCPosition(self, BCDimension, BCSide) result(pos)

        implicit none

        class(SurfaceFace_t)                            :: self
        integer, intent(in)                             :: BCDimension, BCSide
        real(kind=RP)                                   :: pos

        !local variables
        real(kind=RP), dimension(NODES_PER_FACE,2)      :: cornersPosition
        integer                                         :: i

        do i = 1, NODES_PER_FACE
            cornersPosition(i,:) = self % edges(i) % getBCPostion(BCDimension)
        end do

        select case (BCSide)
        case (1)
            pos = minval(cornersPosition)
        case (2)
            pos = maxval(cornersPosition)
        end select

    End Function FaceGetBCPosition
!
!////////////////////////////////////////////////////////////////////////
!
    Function FaceShareCorner(self, otherFace) result(share)

        implicit none

        class(SurfaceFace_t)                            :: self
        class(SurfaceFace_t), intent(in)                :: otherFace
        logical                                         :: share

        !local variables
        integer                                         :: k, edgeID

        share = .false.
         do edgeID = 1, NODES_PER_FACE
            do k = 1, NODES_PER_FACE
                if (self % edges(edgeID) % isConnected(otherFace % edges(k))) then
                    share = .true.
                    return
                end if 
            end do
        end do 

    End Function FaceShareCorner
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine FaceReconstructPeriodic(self, old, new, oldNewMap)

        class(SurfaceFace_t)                                        :: self
        real(kind=RP), dimension(NDIM,NODES_PER_FACE), intent(in)   :: old, new
        integer, dimension(NODES_PER_FACE), intent(in)              :: oldNewMap

        !local variables
        integer                                                     :: i, j, k

        do i = 1, NODES_PER_FACE
            do j = 1, 2
                old_loop: do k = 1, NODES_PER_FACE
                    if (testEqualCorners(self % edges(i) % corners(:,j), old(:,k))) then
                        self % edges(i) % corners(:,j) = new(:,oldNewMap(k))
                        exit old_loop
                    end if 
                end do old_loop
            end do
        end do

    End Subroutine FaceReconstructPeriodic
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine FaceUpdateEdgesOnPeriodic(self, surfaceElements, N)

        implicit none

        class(SurfaceFace_t)                                :: self
        class(SurfaceElement_t), dimension(N), intent(in)   :: surfaceElements
        integer, intent(in)                                 :: N

        !local variables
        integer                                             :: i, j
        integer                                             :: eIndex, faceIndex

        ! get the index of the element and the face
        do i = 1, N
            if (self % globaleID .eq. surfaceElements(i) % globaleID) then
                eIndex = i
                exit
            end if 
        end do

        do i = 1, FACES_PER_ELEMENT
            if (self % fID .eq. surfaceElements(eIndex) % faces(i) % fID) then
                faceIndex = i
                exit
            end if 
        end do

        !replace corners of edges with element ones that have been modified by reconstructPeriod
        do i = 1, NODES_PER_FACE
            self % edges(i) % corners = surfaceElements(eIndex) % faces(faceIndex) % edges(i) % corners
        end do


    End Subroutine FaceUpdateEdgesOnPeriodic
!
!////////////////////////////////////////////////////////////////////////
! surface face methods
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementConstruct(self, mesh, eID, geID, fID, zoneMarkers, M, notcheckBC)

        implicit none

        class(SurfaceElement_t)                         :: self
        type(HexMesh), intent(in)                       :: mesh
        integer, intent(in)                             :: eID, geID, fID, M
        integer, dimension(M), intent(in)               :: zoneMarkers
        logical, intent(in)                             :: notcheckBC

        !local variables
        integer                                         :: i, localFaceID

        self % eID = eID
        self % globaleID = geID
        self % fID = fID

        allocate(self % faces(NUM_OF_NEIGHBORS))
        do i = 1, NUM_OF_NEIGHBORS
            localFaceID = mesh % elements(eID) % faceIDs(i)
            call self % faces(i) % construct(mesh, eID, geID, localFaceID)
        end do

        self % isInBCZone = .false.
        if (.not. notcheckBC) call self % updateIsInZone(mesh, zoneMarkers, M)

        ! is false by default, will be change if necessary in future
        self % needSecondFace = .false.

    End Subroutine ElementConstruct
! 
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementDestruct(self)

        implicit none
        class(SurfaceElement_t), intent(inout)          :: self

        !local variables
        integer                                         :: i

        safedeallocate(self % extrafIDs)
        do i = 1, NUM_OF_NEIGHBORS
            call self % faces(i) % destruct()
        end do
        safedeallocate(self % faces)

    End Subroutine ElementDestruct
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementUpdateIsInZone(self, mesh, zoneMarkers, M)

        implicit none

        class(SurfaceElement_t)                         :: self
        type(HexMesh), intent(in)                       :: mesh
        integer, dimension(M), intent(in)               :: zoneMarkers
        integer, intent(in)                             :: M

        !local variables
        integer                                         :: zoneID, fID, i, j
        integer, dimension(2)                           :: eIDs
        real(kind=RP)                                   :: limits

        do j = 1, M
            zoneID = zoneMarkers(j)
                do i = 1, mesh % zones(zoneID) % no_of_faces
                    fID = mesh % zones(zoneID) % faces(i)
                    eIDs = mesh % faces(fID) % elementIDs
                    if (any(eIDs .eq. self % eID)) then
                        self % isInBCZone = .true.
                        return
                    end if 
                end do
        end do 

    End Subroutine ElementUpdateIsInZone
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementSetNeedSecond(self, N)

        implicit none

        class(SurfaceElement_t)                         :: self
        integer, intent(in)                             :: N

        if (N .eq. 0) return
        allocate(self % extrafIDs(N))
        self % needSecondFace = .true.
        self % extrafIDs(1) = 0

    End Subroutine ElementSetNeedSecond
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementSetNotNeedSecond(self)

        implicit none

        class(SurfaceElement_t)                         :: self

        if (.not. self % needSecondFace) return
        if (self % extrafIDs(1) .ne. 0) return
        safedeallocate(self % extrafIDs)
        self % needSecondFace = .false.

    End Subroutine ElementSetNotNeedSecond
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementGetNotConnectedN(self, mesh, surfaceElements, N, newFaceID)

        use ElementConnectivityDefinitions, only: normalAxis
        implicit none

        class(SurfaceElement_t), target                 :: self
        type(HexMesh), intent(in)                       :: mesh
        type(SurfaceElement_t), dimension(N)            :: surfaceElements
        integer, intent(in)                             :: N
        integer, intent(out)                            :: newFaceID

        !local variables
        integer                                         :: i, j, k, numE, eID, targetEID, edgeIndex, targertEfID, faceIndex, normalFace
        integer, dimension(2)                           :: thisFaceIndexes
        class(SurfaceFace_t), pointer                   :: surfaceFace
        integer, dimension(:), allocatable              :: allElements, connectedElements, possibleElements

        if ( .not. self % needSecondFace ) return
        newFaceID = 0
        targetEID = 0

        do i = 1, NUM_OF_NEIGHBORS
            if (self % faces(i) % fID .eq. self % fID) surfaceFace => self % faces(i)
        end do

    ! first try with the opposite element to the actual face
        ! first get the opposite element
        do i = 1, NUM_OF_NEIGHBORS
            if (mesh % elements(self % eID) % faceIDs(i) .ne. self % fID) cycle
            thisFaceIndexes(1) = i
            exit
        end do
        normalFace = normalAxis(thisFaceIndexes(1)) * (-1)
        ! thisFaceIndexes(2) = findloc(normalAxis, normalFace, dim=1)
        thisFaceIndexes(2) = maxloc(merge(1.0, 0.0, normalAxis == normalFace), dim=1)

        do eID = 1, N
            if (mesh % elements(self % eID) % Connection(thisFaceIndexes(2)) % globID .eq. surfaceElements(eID) % globaleID) then
                targetEID = eID
                exit
                end if 
        end do

        if (targetEID .ne. 0) then

            ! get face of self which is connected to the element face and is connected to the old face
            do j = 1, NUM_OF_NEIGHBORS
                if (surfaceElements(targetEID) % faces(j) % fID .eq. surfaceElements(targetEID) % fID) then
                   targertEfID = j
                   exit
                end if 
            end do
            do i = 1, NUM_OF_NEIGHBORS
                if (self % faces(i) % isConnected(surfaceElements(targetEID) % faces(targertEfID))) then
                    if (.not. surfaceFace % isConnected(self % faces(i))) cycle
                    newFaceID = self % faces(i) % fID
                    self % extrafIDs(1) = newFaceID
                    exit
                end if 
            end do
        end if 
        if (newFaceID .ne. 0) return

! ------------------
        ! if not, get all connected elements

        allocate(allElements(NUM_OF_NEIGHBORS+3))
        numE = 0
        elems_loop: do eID = 1, N
            if (numE .ge. NUM_OF_NEIGHBORS+3) exit elems_loop
            if ( self % globaleID .eq. surfaceElements(eID)  % globaleID ) cycle elems_loop
            this_elem_loop: do i = 1, NUM_OF_NEIGHBORS
                if (associated(surfaceFace, target=self % faces(i))) cycle this_elem_loop
                ext_elem_loop: do j = 1, NUM_OF_NEIGHBORS
                    if (self % faces(i) % isConnected(surfaceElements(eID) % faces(j))) then
                        numE = numE + 1
                        allElements(numE) = eID
                        cycle elems_loop
                    end if 
                end do ext_elem_loop 
            end do this_elem_loop
        end do elems_loop 

        allocate(connectedElements(numE))
        connectedElements = allElements(1:numE)
        deallocate(allElements)

        ! get the element that is not connected to the already found surface
        ! first get the elements with non shared edges
        allocate(allElements(numE))
        targetEID = 0
        i = 0
        con_elems_loop: do k = 1, numE
            eID = connectedElements(k)
             do j = 1, NUM_OF_NEIGHBORS
                if (surfaceFace % isConnected(surfaceElements(eID) % faces(j))) cycle con_elems_loop
            end do
            i = i + 1
            allElements(i) = eID
        end do con_elems_loop 

        ! now reduce to non connected by corner if necessary
        numE = i
        if (numE .eq. 1) then
            targetEID = allElements(1)
            deallocate(allElements)
        else
            allocate(possibleElements(numE))
            possibleElements = allElements(1:numE)
            deallocate(allElements)
            i = 0
            pos_elems_loop: do k = 1, numE
                eID = possibleElements(k)
                do j = 1, NUM_OF_NEIGHBORS
                    if (surfaceFace % shareCorner(surfaceElements(eID) % faces(j))) cycle pos_elems_loop
                end do
                targetEID = eID
                exit
                i = i + 1
            end do pos_elems_loop
            ! if (i .gt. 1) print *, "Warning: more than one element is possible to be neighbour"
        end if 

        if (targetEID .eq. 0) then 
            ! print *, "Warning: for element ", self%globaleID, "a neighbour was not found"
            return
        end if

        ! get face of self which is connected to the element face and is connected to the old face
        do j = 1, NUM_OF_NEIGHBORS
            if (surfaceElements(targetEID) % faces(j) % fID .eq. surfaceElements(targetEID) % fID) then
               targertEfID = j
               exit
            end if 
        end do
        do i = 1, NUM_OF_NEIGHBORS
            if (self % faces(i) % isConnected(surfaceElements(targetEID) % faces(targertEfID))) then
                if (.not. surfaceFace % isConnected(self % faces(i))) cycle
                newFaceID = self % faces(i) % fID
                self % extrafIDs(1) = newFaceID
                exit
            end if 
        end do

        if (newFaceID .ne. 0) return

        ! if not get face of self which is connected to some other element face and is connected to the old face
        this_face_loop: do i = 1, NUM_OF_NEIGHBORS
            do j = 1, NUM_OF_NEIGHBORS
                if (self % faces(i) % isConnected(surfaceElements(targetEID) % faces(j))) then
                    if (.not. surfaceFace % isConnected(self % faces(i))) cycle this_face_loop
                    newFaceID = self % faces(i) % fID
                    self % extrafIDs(1) = newFaceID
                    exit this_face_loop
                end if 
            end do
        end do this_face_loop

    End Subroutine ElementGetNotConnectedN
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine ElementReconstructPeriodic(self)

        use ElementConnectivityDefinitions, only: NODES_PER_ELEMENT
        implicit none

        class(SurfaceElement_t), target                     :: self

        !local variables
        class(SurfaceFace_t), pointer                       :: faceNew=>null(), faceOpposite=>null()
        integer                                             :: faceConnections
        real(kind=RP), dimension(:,:), allocatable          :: allCorners
        real(kind=RP), dimension(NDIM,FACES_PER_ELEMENT)    :: newCorners, notChangeCorners, oldCorners
        integer                                             :: i,j,k
        integer, dimension(NODES_PER_FACE)                  :: oldNewCornersMap
        real(kind=RP), dimension(NODES_PER_FACE)            :: residual

        ! get the new face and the opposite to it: the one without any connections and the one with FACES_PER_ELEMENT-2 connections
        ! --------------------------------------
        faces_loop: do i = 1, FACES_PER_ELEMENT
            faceConnections = 0
            do j = 1, FACES_PER_ELEMENT
                if (i .eq. j) cycle
                if (self % faces(i) % isConnected(self % faces(j))) faceConnections = faceConnections + 1
            end do
            if (faceConnections .eq. 0) then
                faceNew => self % faces(i)
                cycle faces_loop
            elseif (faceConnections .eq. FACES_PER_ELEMENT-2) then
                faceOpposite => self % faces(i)
                cycle faces_loop
            end if 
        end do faces_loop

        if (.not. associated(faceNew)) then
            return
        end if

        ! get corners arrays
        ! -----------------

        ! first get all the cornes, new and original
        allocate( allCorners(NDIM,NODES_PER_FACE*FACES_PER_ELEMENT) )
        do i = 1, FACES_PER_ELEMENT
            do j = 1, NODES_PER_FACE
                allCorners(:,(i-1)*NODES_PER_FACE+j) = self % faces(i) % edges(j) % corners(:,1)
            end do
        end do 
        ! get new and opposite corners
        do j = 1, NODES_PER_FACE
            newCorners(:,j) = faceNew % edges(j) % corners(:,1)
            notChangeCorners(:,j) = faceOpposite % edges(j) % corners(:,1)
        end do
        ! get old corners, the ones that are not new nor opposite corners, remove duplicated
        k = 0
        outer: do i = 1, size(allCorners(1,:))
            if (k .ge. NODES_PER_FACE) exit outer
            do j = 1, NODES_PER_FACE
                if (testEqualCorners(newCorners(:,j), allCorners(:,i))) cycle outer
                if (testEqualCorners(notChangeCorners(:,j), allCorners(:,i))) cycle outer
            end do
            do j = 1, k
                if (testEqualCorners(oldCorners(:,j), allCorners(:,i))) cycle outer
            end do
            k = k + 1
            oldCorners(:,k) = allCorners(:,i)
        end do outer

        ! map old vs new
        ! --------------
        do i = 1, NODES_PER_FACE
            residual = huge(newCorners(1,1))
            do j = 1, NODES_PER_FACE
                residual(j) = POW2(oldCorners(1,i) - newCorners(1,j)) + POW2(oldCorners(2,i) - newCorners(2,j))
            end do
            oldNewCornersMap(i) = minloc(residual,dim=1)
        end do

        do i = 1, FACES_PER_ELEMENT
            if (associated(faceNew, target=self % faces(i))) cycle
            if (associated(faceOpposite, target=self % faces(i))) cycle
            call self % faces(i) % reconstructPeriod(oldCorners, newCorners, oldNewCornersMap)
        end do

        nullify(faceNew, faceOpposite)

    End Subroutine ElementReconstructPeriodic
!
!////////////////////////////////////////////////////////////////////////
! helper procedures
!////////////////////////////////////////////////////////////////////////
!
    Function testEqualCorners(cornerA, cornerB) result(isEqual)

        use Utilities, only: AlmostEqual
        implicit none

        real(kind=RP), dimension(3), intent(in)         :: cornerA, cornerB
        logical                                         :: isEqual

        isEqual = AlmostEqual(cornerA(1), cornerB(1)) .and. &
                  AlmostEqual(cornerA(2), cornerB(2)) .and. &
                  AlmostEqual(cornerA(3), cornerB(3))

    End Function testEqualCorners
!
!////////////////////////////////////////////////////////////////////////
!

End Module SurfaceClass