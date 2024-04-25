#include "Includes.h"
Module FWHPreSurface  !
    use SMConstants
    use Headers
    use HexMeshClass
    use FTValueDictionaryClass
    use SurfaceClass
    Implicit None

    private
    public extractSurface

    interface
        Function elementInGeo(r, ratio, centerPosition, x, y, N, lengthAspect, filter) result(isIn)
            use SMConstants
            implicit none
            real(kind=RP)               , intent(in)    :: r, ratio, lengthAspect
            real(kind=RP), dimension(N) , intent(in)    :: x, y
            real(kind=RP), dimension(2) , intent(in)    :: centerPosition
            integer                     , intent(in)    :: N
            logical, optional           , intent(in)    :: filter
            logical                                     :: isIn
        End Function elementInGeo

        Function distaceToGeo(r, ratio, centerPosition, x, y) result(rdiff)

            use SMConstants
            implicit none

            real(kind=RP), intent(in)                   :: r, ratio, x, y
            real(kind=RP), dimension(2), intent(in)     :: centerPosition
            real(kind=RP)                               :: rdiff
        End Function distaceToGeo

    end interface

!==============
    contains
!==============
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine extractSurface(mesh, controlVariables)

        use mainKeywordsModule
        use FileReadingUtilities, only: getFileName, getRealArrayFromString, RemovePath, getCharArrayFromString, getIntArrayFromString
        use Utilities,            only: toLower

        implicit none

        type(HexMesh)                                       :: mesh
        type(FTValueDictionary)                             :: controlVariables

        !local variables
        type(Surface_t), allocatable                       :: surface
        real(kind=RP)                                       :: radii, ratio, lengthAspect
        real(kind=RP), dimension(:), allocatable            :: centerPosition
        integer, dimension(:), allocatable                  :: eIDs, globaleIDs, numberOfFacesperElement, firstFIDs
        integer                                             :: numberOfElements, numberOfFaces, numberOfSecondFaces, numberOfVerifiedFaces
        character(len=LINE_LENGTH)                          :: geometry, fileName, meshName, center, boundaryNames, newElems, discElems
        character(len=LINE_LENGTH), allocatable             :: boundaryNamesArr(:)
        procedure(elementInGeo), pointer                    :: isInGeometry => null()
        procedure(distaceToGeo), pointer                    :: rToGeometry => null()
        type(SurfaceFace_t), dimension(:), allocatable     :: firstFaces, allCreatedFaces
        type(SurfaceFace_t), dimension(:), allocatable     :: secondFaces, oldFaces, oldTemp
        type(SurfaceFace_t), allocatable     :: tempSecFaces
        type(SurfaceElement_t), dimension(:), allocatable  :: elementsOfFaces
        integer                                             :: fID, eID, eIndex
        integer                                             :: i, k, numberOfBCZones, j, numberOfNewFaces
        logical                                             :: connectedAtBoundary, discardElems, includeElems, useFilter, saveForMPI
        integer, dimension(:), allocatable                  :: zoneMarkers

        integer, dimension(:), allocatable                  :: allE, allgE, allF
        integer, dimension(:), allocatable                  :: nE, newElemID, allNewElemID, allExcludedElemID, discElemID, secIdsNot
        integer, dimension(:,:), allocatable                :: oldIds, secIds, realSecIds, oldIdsTemp, secIdsTemp, bcFaces
        integer                                             :: nNe, NelemDir, numberOfBCElemNew, numberOfBCElemDisc, numberOfBCFaces

!       Get the solution and mesh file names
!       ------------------------------------
        fileName = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = LINE_LENGTH )
        fileName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".surfaceZone"
        meshName = trim(controlVariables % stringValueForKey( meshFileNameKey, requestedLength = LINE_LENGTH ))


!       Get the geometry info
!       ---------------------
        geometry = controlVariables % stringValueForKey("geometry", LINE_LENGTH)
        call toLower(geometry)
        radii = controlVariables % doublePrecisionValueForKey("radii")
        ratio = controlVariables % getValueOrDefault("ratio", 1.0_RP)
        connectedAtBoundary = controlVariables % getValueOrDefault("include boundaries", .false.)
        center = trim(controlVariables % stringValueForKey("center position",LINE_LENGTH))
        centerPosition = getRealArrayFromString(center)
        lengthAspect = controlVariables % getValueOrDefault("surface tolerance", 0.7_RP)

!       Get the BC info
!       ----------------
        if (controlVariables%containsKey("boundary names")) then
            boundaryNames = trim(controlVariables % stringValueForKey("boundary names",LINE_LENGTH))
        else
            boundaryNames = ""
        end if
        call toLower(boundaryNames)
        call getCharArrayFromString(boundaryNames, LINE_LENGTH, boundaryNamesArr)
        numberOfBCZones = size(boundaryNamesArr)
        allocate(zoneMarkers(numberOfBCZones))
        zoneMarkers = 0
        bc_loop: do k = 1, numberOfBCZones
            do i = 1, size(mesh % zones)
                if (trim(mesh % zones(i) % Name) .eq. trim(boundaryNamesArr(k))) then
                    zoneMarkers(k) = mesh % zones(i) % marker
                    cycle bc_loop
                end if 
            end do
        end do bc_loop

!       Get the misc info
!       -----------------
        NelemDir = controlVariables % getValueOrDefault("elements in direction", 0)
        useFilter = controlVariables % getValueOrDefault("use filtering", .false.)
        saveForMPI = controlVariables % getValueOrDefault("save for mpi", .false.)

        includeElems = controlVariables % containsKey("elements to include")
        newElems = trim(controlVariables % stringValueForKey("elements to include",LINE_LENGTH))
        if (newElems .eq. "") newElems = "[0]"
        newElemID = getIntArrayFromString(newElems)
        numberOfBCElemNew = size(newElemID)

        discardElems = controlVariables % containsKey("elements to exclude")
        discElems = trim(controlVariables % stringValueForKey("elements to exclude",LINE_LENGTH))
        if (discElems .eq. "") discElems = "[0]"
        discElemID = getIntArrayFromString(discElems)
        numberOfBCElemDisc = size(discElemID)

!       get the apropiated functions
!       ----------------------------

        select case (trim(geometry))
        case ("circle")
            isInGeometry => isInCircle
            rToGeometry  => distaceToCircle
        case ("ellipse")
            isInGeometry => isInEllipse
            rToGeometry  => distanceToEllipse
        case default
            print *, "Warning, no geometry defined, using circle as default"
            isInGeometry => isInCircle
            rToGeometry  => distaceToCircle
        end select

!       get the forced elements to include
!       ----------------------------
        allocate(allNewElemID(numberOfBCElemNew*NelemDir))
        allNewElemID = 0
        if (includeElems) then
            k = 1
            do i = 1, numberOfBCElemNew
                call getAdditionalElements(mesh, newElemID(i), NelemDir, zoneMarkers, numberOfBCZones, nE, nNe)
                allNewElemID(k:k+nNe-1) = nE
                k = k + nNe
                safedeallocate(nE)
            end do
        end if 

!       get the forced elements to exclude
!       ----------------------------
        allocate(allExcludedElemID(numberOfBCElemDisc*NelemDir))
        allExcludedElemID = 0
        if (discardElems) then
            k = 1
            do i = 1, numberOfBCElemDisc
                call getAdditionalElements(mesh, discElemID(i), NelemDir, zoneMarkers, numberOfBCZones, nE, nNe)
                allExcludedElemID(k:k+nNe-1) = nE
                k = k + nNe
                safedeallocate(nE)
            end do
        end if 

!       get the elements ids
!       --------------------
        call getElements(mesh, isInGeometry, radii, ratio, lengthAspect, centerPosition, allNewElemID, allExcludedElemID, eIDs, &
                         globaleIDs, useFilter)
        numberOfElements = size(eIDs)

!       get the first faces ids and construct the faces
!       ----------------------------------------------
        allocate(firstFIDs(numberOfElements))
        allocate(firstFaces(numberOfElements))
        do i = 1, numberOfElements
            call getFirstFace(mesh, rToGeometry, eIDs(i), radii, ratio, centerPosition, firstFIDs(i))
            call firstFaces(i) % construct(mesh, eIDs(i), globaleIDs(i), firstFIDs(i))
        end do     

!       construct elements
!       ------------------
        allocate(elementsOfFaces(numberOfElements))
        do i = 1, numberOfElements
            call elementsOfFaces(i) % construct(mesh, eIDs(i), globaleIDs(i), firstFIDs(i), zoneMarkers, numberOfBCZones, .false.)
            call elementsOfFaces(i) % reconstructPeriod()
        end do     

!       update faces connections
!       ------------------------
        ! if (.not. connectedAtBoundary) then
            do i = 1, numberOfElements
            call firstFaces(i)%updateEdgesPeriod(elementsOfFaces, numberOfElements)
                if (.not. elementsOfFaces(i) % isInBCZone) cycle
                call firstFaces(i) % setNoConnections(mode=1)
            end do
        ! end if 

!       get the elements that need second faces
!       ---------------------------------------
        allocate(numberOfFacesperElement(numberOfElements))
        numberOfFacesperElement = 2
        do i = 1, numberOfElements
            if ( firstFaces(i) % isFullConnected(firstFaces, numberOfElements) ) then
                numberOfFacesperElement(i) = 1
            else
                call elementsOfFaces(i) % setNeedSecond(numberOfFacesperElement(i) - 1)
            end if
        end do     

!       get first total numbers
!       -----------------------
        numberOfFaces = sum(numberOfFacesperElement)
        numberOfSecondFaces = numberOfFaces - numberOfElements
        numberOfVerifiedFaces = numberOfElements

        allocate(oldIds(numberOfElements,3))
        oldIds(:,1) = eIDs
        oldIds(:,2) = globaleIDs
        oldIds(:,3) = firstFIDs
        allocate(oldFaces(numberOfElements))
        oldFaces = firstFaces

!       start loop to obtain extra faces
!       --------------------------------

        do j = 1, 4
        ! print *, "j: ", j

!       get the second faces ids
!       ------------------------
        allocate(secIds(numberOfSecondFaces,3))

        k = 1
        do i = 1, numberOfElements
            if ( (.not. elementsOfFaces(i) % needSecondFace) .or. (elementsOfFaces(i) % extrafIDs(1) /= 0) ) cycle
            call elementsOfFaces(i) % getNotConnectedN(mesh, elementsOfFaces, numberOfElements, secIds(k,3))
            secIds(k,1) = eIDs(i)
            secIds(k,2) = globaleIDs(i)
            if ( k .ge. numberOfSecondFaces ) exit
            k = k + 1
        end do     
        if (k .le. 1) exit
        numberOfNewFaces = k

!       create connected secondFaces and discard the non well positioned
!       ----------------------------------------------------------------
        if (j .le. 1) then
            allocate( allCreatedFaces(numberOfElements+1) )
            allCreatedFaces(1:numberOfElements) = oldFaces
            allocate( secIdsNot(numberOfNewFaces) )

            k = 0
            do i = 1, numberOfNewFaces
                if (secIds(i,3) .eq. 0) then
                    ! eIndex = findloc(globaleIDs, secIds(i,1), dim=1)
                    eIndex = maxloc(merge(1.0, 0.0, globaleIDs == secIds(i,1)), dim=1)
                    call elementsOfFaces(eIndex) % setNeedNotSecond()
                    if (any(secIdsNot .eq. secIds(i,2))) cycle
                    k = k + 1
                    secIdsNot(k) = secIds(i,2)
                    cycle
                end if
                if (any(secIdsNot .eq. secIds(i,2))) cycle
                allocate(tempSecFaces)
                call tempSecFaces % construct(mesh, secIds(i,1), secIds(i,2), secIds(i,3))
                allCreatedFaces(numberOfVerifiedFaces+1) = tempSecFaces
                if (tempSecFaces % isTwiceEdConnected(allCreatedFaces, numberOfElements+1)) then
                ! if (tempSecFaces % fee(allCreatedFaces, numberOfElements+1)) then
                    if (any(secIdsNot .eq. secIds(i,2))) cycle
                    k = k + 1
                    secIdsNot(k) = secIds(i,2)
                    ! eIndex = findloc(globaleIDs, secIds(i,1), dim=1)
                    eIndex = maxloc(merge(1.0, 0.0, globaleIDs == secIds(i,1)), dim=1)
                    call elementsOfFaces(eIndex) % setNeedNotSecond()
                else
                    do eID = 1, numberOfElements
                        if (.not. elementsOfFaces(eID) % needSecondFace) cycle
                        if (elementsOfFaces(eID) % globaleID .eq. secIds(i,2)) cycle
                        if (any(secIdsNot .eq. elementsOfFaces(eID) % globaleID)) cycle
                        if (.not. any(secIds(:,2) .eq. elementsOfFaces(eID) % globaleID)) cycle

                        ! see if with the new face the original face of the element if connected, to remove second connection
                        if (allCreatedFaces(eID) % isFullConnected(allCreatedFaces, numberOfElements+1)) then
                            k = k + 1
                            call elementsOfFaces(eID) % setNeedNotSecond()
                            secIdsNot(k) = elementsOfFaces(eID) % globaleID
                            exit
                        end if
                    end do
                end if
                call tempSecFaces % destruct
                deallocate(tempSecFaces)
            end do
            deallocate(allCreatedFaces)
            if (k .gt. 0) then
                numberOfNewFaces = numberOfNewFaces - k
                numberOfFaces = numberOfElements + numberOfNewFaces
                numberOfSecondFaces = numberOfFaces - numberOfElements
                allocate(secIdsTemp(numberOfNewFaces,3))
                k = 1
                do i = 1, size(secIds, dim=1)
                    if (k .gt. numberOfNewFaces) exit
                    if (any(secIds(i,2) .eq. secIdsNot)) cycle
                    secIdsTemp(k,:) = secIds(i,:)
                    k = k+1
                end do
                deallocate(secIds)
                allocate(secIds(numberOfNewFaces,3))
                secIds = secIdsTemp(1:numberOfNewFaces,:)
                deallocate(secIdsTemp)
            end if 
            deallocate(secIdsNot)
        end if

!       create all secondFaces
!       ----------------------
        allocate( secondFaces(numberOfNewFaces) )
        do i = 1, numberOfNewFaces
            call secondFaces(i) % construct(mesh, secIds(i,1), secIds(i,2), secIds(i,3))
            call secondFaces(i)%updateEdgesPeriod(elementsOfFaces, numberOfElements)
        end do
        deallocate(secIds)

!       update faces connections
!       ------------------------
        ! if (.not. connectedAtBoundary) then
            do i = 1, numberOfSecondFaces
                eID = secondFaces(i) % globaleID
                ! eIndex = findloc(globaleIDs, eID, dim=1)
                eIndex = maxloc(merge(1.0, 0.0, globaleIDs == eID), dim=1)
                if (.not. elementsOfFaces(eIndex) % isInBCZone) cycle
                if (.not. elementsOfFaces(eIndex) % needSecondFace) cycle
                call secondFaces(i) % setNoConnections(mode=1)
            end do
        ! end if 

        ! set num of con for sec faces at bc
        allocate( allCreatedFaces(numberOfFaces) )
        allCreatedFaces(1:numberOfVerifiedFaces) = oldFaces
        allCreatedFaces(numberOfVerifiedFaces+1:numberOfFaces) = secondFaces

!       get all possible secondFaces not full connected
!       -----------------------------------------------
        allocate( realSecIds(numberOfFaces,3) )
        k = 0
        do i = numberOfVerifiedFaces + 1, numberOfFaces
            if (allCreatedFaces(i) % isFullConnected(allCreatedFaces, numberOfFaces)) then
                k = k + 1
                realSecIds(k,1) = allCreatedFaces(i) % eID
                realSecIds(k,2) = allCreatedFaces(i) % globaleID
                realSecIds(k,3) = allCreatedFaces(i) % fID
            else
                eID = allCreatedFaces(i) % globaleID
                ! eIndex = findloc(globaleIDs, eID, dim=1)
                eIndex = maxloc(merge(1.0, 0.0, globaleIDs == eID), dim=1)
                if (.not. elementsOfFaces(eIndex) % needSecondFace) cycle
                elementsOfFaces(eIndex) % extrafIDs(1) = 0
            end if 
        end do

!       destroy previous faces
!       ----------------------
        do i = 1, numberOfSecondFaces
            call secondFaces(i) % destruct()
        end do
        deallocate(secondFaces)
        deallocate(allCreatedFaces)

        numberOfSecondFaces = k

!       exit loop if not new faces are created
!       --------------------------------------
        if (numberOfSecondFaces .eq. 0) exit

!       create real secondFaces
!       -----------------------
        allocate( secondFaces(numberOfSecondFaces) )
        do i = 1, numberOfSecondFaces
            call secondFaces(i) % construct(mesh, realSecIds(i,1), realSecIds(i,2), realSecIds(i,3))
            call secondFaces(i) % updateEdgesPeriod(elementsOfFaces, numberOfElements)
        end do

        k = numberOfVerifiedFaces - numberOfElements

!       update old ids with the real new ones
!       -------------------------------------
        allocate( oldIdsTemp(numberOfVerifiedFaces,3) )
        oldIdsTemp = oldIds
        deallocate(oldIds)
        allocate( oldIds(numberOfVerifiedFaces+numberOfSecondFaces,3) )
        oldIds(1:numberOfVerifiedFaces,:) = oldIdsTemp
        oldIds(numberOfVerifiedFaces+1 : numberOfVerifiedFaces+numberOfSecondFaces,:) = realSecIds
        deallocate(oldIdsTemp)
        deallocate(realSecIds)

!       update old Faces with the real new ones
!       ---------------------------------------
        allocate( oldTemp(numberOfVerifiedFaces) )
        oldTemp = oldFaces
        deallocate(oldFaces)
        allocate( oldFaces( numberOfVerifiedFaces+numberOfSecondFaces ) )
        oldFaces(1:numberOfVerifiedFaces) = oldTemp
        oldFaces(numberOfVerifiedFaces+1 : numberOfVerifiedFaces+numberOfSecondFaces) = secondFaces
        deallocate(oldTemp)
        deallocate(secondFaces)

        numberOfVerifiedFaces = numberOfVerifiedFaces + numberOfSecondFaces

!       update elements that need secondFaces based on all real secondFaces
!       -------------------------------------------------------------------
        do i = 1, numberOfVerifiedFaces
            eID = oldFaces(i) % globaleID
            ! eIndex = findloc(globaleIDs, eID, dim=1)
            eIndex = maxloc(merge(1.0, 0.0, globaleIDs == eID), dim=1)
            if ( (.not. elementsOfFaces(eIndex) % needSecondFace) .or. (elementsOfFaces(eIndex) % extrafIDs(1) .ne. 0) ) cycle
            if (oldFaces(i) % isFullConnected(oldFaces, numberOfVerifiedFaces)) then
                call elementsOfFaces(eIndex) % setNeedNotSecond()
                numberOfFaces = numberOfFaces - 1
            end if 
        end do

        numberOfSecondFaces = numberOfFaces - numberOfVerifiedFaces

        ! exit loop if no more faces are needed
        if (numberOfFaces .eq. numberOfVerifiedFaces) exit
        end do
        if (numberOfFaces .ne. numberOfVerifiedFaces) print *, "Warning, not all possible faces have been created. Check results first"

!       get boundary faces if needed
!       ----------------------------
        if (connectedAtBoundary) then
            call getBCFaces(mesh, isInGeometry, radii, ratio, centerPosition, lengthAspect, allExcludedElemID, allNewElemID, zoneMarkers, numberOfBCZones, bcFaces, numberOfBCFaces)
            allocate(oldIdsTemp(numberOfVerifiedFaces+numberOfBCFaces,3))
            oldIdsTemp(1:numberOfVerifiedFaces,:) = oldIds
            oldIdsTemp(numberOfVerifiedFaces+1:,:) = bcFaces
            call move_alloc(oldIdsTemp,oldIds)
            numberOfVerifiedFaces = numberOfVerifiedFaces + numberOfBCFaces
        end if 

!       create final surface and create output files
!       --------------------------------------------
        allocate(surface)
        call surface % construct(mesh, fileName, oldIds(:,1), oldIds(:,2), oldIds(:,3), numberOfVerifiedFaces)
        call surface % writeToTecplot(mesh, meshName)
        call surface % saveToFile(mesh, saveForMPI)

        call surface % destruct()

    End Subroutine extractSurface
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine getElements(mesh, isInGeometry, radii, ratio, lengthAspect, centerPosition, toInclude, toExclude, eIDArray, geIDArray, useFilter)

        use ElementConnectivityDefinitions, only: NODES_PER_ELEMENT, FACES_PER_ELEMENT, normalAxis
        use MeshTypes,                      only: emptyBCName
        implicit none

        type(HexMesh)                     , intent(in)          :: mesh
        procedure(elementInGeo)                                 :: isInGeometry
        real(kind=RP)                     , intent(in)          :: radii, ratio, lengthAspect
        real(kind=RP), dimension(2)       , intent(in)          :: centerPosition
        integer, dimension(:), intent(in)                       :: toInclude, toExclude
        integer, dimension(:), allocatable, intent(out)         :: geIDArray,eIDArray
        logical, intent(in)                                     :: useFilter

        !local variables
        integer                                                 :: eID, i, j, k
        integer, dimension(:,:), allocatable                    :: allElements, elementsTemp
        integer                                                 :: nElements, extrudedDirectionIndex, numberOfNeighbours, totalNeighbours
        integer, dimension(:), allocatable                      :: neighboursIndex

        call getAllInteriorElements(mesh, isInGeometry, radii, ratio, centerPosition, lengthAspect, useFilter, toExclude, toInclude, allElements, nElements)

        if (.not. useFilter) then
            allocate(elementsTemp(nElements,2))
            elementsTemp = allElements(1:nElements,:)
            deallocate(allElements)
            allocate(allElements(nElements,2))
            k = 0
            do eID = 1, nElements
                ! works if there's only one element neighbour for each face
                totalNeighbours = FACES_PER_ELEMENT
                numberOfNeighbours = 0
                associate ( e => mesh % elements(elementsTemp(eID,1)) )
                    do j = 1, FACES_PER_ELEMENT
                        if ( any(elementsTemp(:,2) .eq. e % Connection(j) % globID) ) then
                            numberOfNeighbours = numberOfNeighbours + 1
                        end if 
                        ! remove one constrain for each boundary of the element
                        if (trim(e % boundaryName(j)) .ne. emptyBCName) then
                            totalNeighbours = totalNeighbours - 1
                        end if 
                    end do
                end associate
                ! is at the surface boundary if there isn't an element inside the surface in the "wall" outgoing surface direction
                if (numberOfNeighbours .lt. totalNeighbours) then
                    k = k + 1
                    allElements(k,:) = elementsTemp(eID,:)
                end if 
            end do 
            nElements = k
            deallocate(elementsTemp)
        end if 

        allocate(geIDArray(nElements), eIDArray(nElements))
        eIDArray = allElements(1:nElements,1)
        geIDArray = allElements(1:nElements,2)
        deallocate(allElements)

    End Subroutine getElements
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine getAllInteriorElements(mesh, isInGeometry, radii, ratio, centerPosition, lengthAspect, useFilter, toExclude, toInclude, allElements, nElements)


        use ElementConnectivityDefinitions, only: NODES_PER_ELEMENT
        implicit none

        type(HexMesh)                     , intent(in)          :: mesh
        procedure(elementInGeo)                                 :: isInGeometry
        real(kind=RP)                     , intent(in)          :: radii, ratio, lengthAspect
        real(kind=RP), dimension(2)       , intent(in)          :: centerPosition
        logical, intent(in)                                     :: useFilter
        integer, dimension(:), intent(in)                       :: toInclude, toExclude
        integer, dimension(:,:), allocatable, intent(out)       :: allElements
        integer, intent(out)                                    :: nElements

        !local variables
        real(kind=RP), dimension(:), pointer                    :: x, y, z
        real(kind=RP), dimension(3,NODES_PER_ELEMENT), target   :: xx
        integer                                                 :: eID, k
        integer                                                 :: Nx, Ny, Nz

        allocate(allElements(mesh % no_of_elements,2))
        k = 1

        do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )

                find_or_not: if (any(toExclude .eq. e % globID)) then
                    cycle
                elseif (any(toInclude .eq. e % globID)) then

                    allElements(k,1) = e % eID
                    allElements(k,2) = e % globID
                    k = k + 1

                else find_or_not

                    Nx = e % Nxyz(1)
                    Ny = e % Nxyz(2)
                    Nz = e % Nxyz(3)

                    ! conrners position
                    xx(:,1) = e % geom % x(:,0,0,0)
                    xx(:,2) = e % geom % x(:,Nx,0,0)
                    xx(:,3) = e % geom % x(:,0,Ny,0)
                    xx(:,4) = e % geom % x(:,Nx,Ny,0)
                    xx(:,5) = e % geom % x(:,0,0,Nz)
                    xx(:,6) = e % geom % x(:,Nx,0,Nz)
                    xx(:,7) = e % geom % x(:,0,Ny,Nz)
                    xx(:,8) = e % geom % x(:,Nx,Ny,Nz)

                    x => xx(1,:)
                    y => xx(2,:)
                    z => xx(3,:)

                    ! todo: configure to call x,z or y,z too
                    if (isInGeometry(radii, ratio, centerPosition, x, y, NODES_PER_ELEMENT, lengthAspect, filter=useFilter)) then
                        allElements(k,1) = e % eID
                        allElements(k,2) = e % globID
                        k = k + 1
                    end if
                    nullify(x,y,z)
                end if find_or_not 

            end associate
        end do

        nElements = k - 1

    End Subroutine getAllInteriorElements
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine getBCFaces(mesh, isInGeometry, radii, ratio, centerPosition, lengthAspect, toExclude, toInclude, zoneMarkers, nBC, BCFaces, nBCFaces)

        use ElementConnectivityDefinitions, only: FACES_PER_ELEMENT
        implicit none

        type(HexMesh)                     , intent(in)          :: mesh
        procedure(elementInGeo)                                 :: isInGeometry
        real(kind=RP)                     , intent(in)          :: radii, ratio, lengthAspect
        real(kind=RP), dimension(2)       , intent(in)          :: centerPosition
        integer, dimension(:), intent(in)                       :: toInclude, toExclude
        integer, dimension(nBC), intent(in)                     :: zoneMarkers
        integer, intent(in)                                     :: nBC
        integer, dimension(:,:), allocatable, intent(out)       :: BCFaces
        integer, intent(out)                                    :: nBCFaces

        !local variables
        integer                                                 :: eID, i, k, fID
        integer, dimension(:,:), allocatable                    :: allElements, elementsTemp
        integer                                                 :: nElements
        character(len=LINE_LENGTH), dimension(nBC)              :: boundaryNames

        do i = 1, nBC
            boundaryNames(i) = trim(mesh % zones(zoneMarkers(i)) % Name)
        end do

        !get all BC elements IDs
        call getAllInteriorElements(mesh, isInGeometry, radii, ratio, centerPosition, lengthAspect, .false., toExclude, toInclude, allElements, nElements)
        allocate(elementsTemp(nElements,2))
        elementsTemp = allElements(1:nElements,:)
        deallocate(allElements)
        allocate(allElements(nElements,2))
        k = 0
        elems_loop:do eID = 1, nElements
            associate ( e => mesh % elements(elementsTemp(eID,1)) )
                do i = 1, FACES_PER_ELEMENT
                    if (any(boundaryNames .eq. e % boundaryName(i))) then
                        k = k + 1
                        allElements(k,:) = elementsTemp(eID,:)
                        cycle elems_loop
                    end if 
                end do
            end associate
        end do elems_loop

        nBCFaces = k
        deallocate(elementsTemp)
        allocate(elementsTemp(nBCFaces,2))
        elementsTemp = allElements(1:nBCFaces,:)
        call move_alloc(elementsTemp, allElements)
        
        ! get ID of the face that is on the corresponded BC
        allocate(BCFaces(nBCFaces,3))
        elem_loop:do eID = 1, nBCFaces
            do i = 1, FACES_PER_ELEMENT
                fID = mesh % elements(allElements(eID,1)) % faceIDs(i)
                if (any(zoneMarkers .eq. mesh % faces(fID) % zone)) then
                    BCFaces(eID,1:2) = allElements(eID,:)
                    BCFaces(eID,3) = fID
                    cycle elem_loop
                end if 
            end do
            print *, "Warning, the element ", allElements(eID,1), "didn't found a BC face"
        end do elem_loop

    End Subroutine getBCFaces
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine getFirstFace(mesh, distanceToGeometry, eID, radii, ratio, centerPosition, fID)

        type(HexMesh), intent(in)                   :: mesh
        procedure(distaceToGeo)                     :: distanceToGeometry
        integer, intent(in)                         :: eID
        real(kind=RP), intent(in)                   :: radii, ratio
        real(kind=RP), dimension(2), intent(in)     :: centerPosition
        integer, intent(out)                        :: fID

        ! local variables
        integer                                     :: i, faceIndex
        integer                                     :: Nx,Ny
        integer, dimension(6)                       :: fIDs
        real(kind=RP), dimension(2)                 :: faceCenter
        real(kind=RP), dimension(6)                 :: faceDistance

        fIDs = mesh % elements(eID) % faceIDs
        do i = 1, 6
            associate ( f => mesh % faces(fIDs(i)) )
                Nx = f % Nf(1)
                Ny = f % Nf(2)
                associate (x => f % geom % x)
                    if (x(3,0,0) .eq. x(3,Nx,Ny)) then
                        faceDistance(i) = huge(x(1,0,0))
                        cycle
                    end if
                    faceCenter(1) = (x(1,0,0) + x(1,Nx,Ny)) / 2.0_RP
                    faceCenter(2) = (x(2,0,0) + x(2,Nx,Ny)) / 2.0_RP
                end associate
            end associate
            faceDistance(i) = abs(distanceToGeometry(radii, ratio, centerposition, faceCenter(1), faceCenter(2)))
        end do

        faceIndex = minloc(faceDistance, dim=1)
        fID = fIDs(faceIndex)

    End Subroutine getFirstFace
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine getAdditionalElements(mesh, BCeID, numberOfElements, zoneMarkers, zN, newEIDs, realNumberOfElements)

        use ElementConnectivityDefinitions, only: normalAxis, FACES_PER_ELEMENT
        use ElementClass
        implicit none

        type(HexMesh), target, intent(in)               :: mesh
        integer, intent(in)                             :: BCeID, numberOfElements, zN
        integer, dimension(zN), intent(in)              :: zoneMarkers
        integer, dimension(:), allocatable, intent(out) :: newEIDs
        integer, intent(out)                            :: realNumberOfElements

        ! local variables
        integer                                         :: i, j, neID, zoneID
        class(Element), pointer                         :: e, BCe
        integer                                         :: thisMarker, fBCindex, conIndex, conNormal
        integer, dimension(:), allocatable              :: nEs

        realNumberOfElements = 0
        BCe => mesh%elements(BCeID)

        zone_loop: do i = 1, zN
            zoneID = zoneMarkers(i)
            face_loop: do j = 1, FACES_PER_ELEMENT
                if ( trim(mesh % zones(zoneID) % Name) /= trim(BCe % boundaryName(j)) ) cycle face_loop
                thisMarker = i
                fBCindex = j
                ! thisMarker = zoneID
                exit zone_loop
            end do face_loop
        end do zone_loop

        e => BCe
        allocate( nEs(numberOfElements) )

        i=1
        do i = 1, numberOfElements
            nEs(i) = e%eID
            realNumberOfElements = realNumberOfElements + 1
            conNormal = normalAxis(fBCindex) * (-1)
            ! conIndex = findloc(normalAxis, conNormal, dim=1)
            conIndex = maxloc(merge(1.0, 0.0, normalAxis == conNormal), dim=1)
            neID = e % Connection(conIndex) % globID
            if (neID .eq. 0) exit
            e => mesh % elements(neID)
        end do

        !todo check that new bc have been reached

        allocate(newEIDs(realNumberOfElements))
        newEIDs = nEs(1:realNumberOfElements)
        deallocate(nEs)


    End Subroutine getAdditionalElements
!
!////////////////////////////////////////////////////////////////////////
! isInGeometry functions
!////////////////////////////////////////////////////////////////////////
!
    Function isInEllipse(mayorAxis, ratio, centerPosition, x, y, N, lengthAspect, filter) result(isIn)

        real(kind=RP)               , intent(in)    :: mayorAxis, ratio, lengthAspect
        real(kind=RP), dimension(N) , intent(in)    :: x, y
        real(kind=RP), dimension(2) , intent(in)    :: centerPosition
        integer                     , intent(in)    :: N
        logical, optional           , intent(in)    :: filter
        logical                                     :: isIn

        ! local variables
        integer                                     :: i, j, k
        real(kind=RP)                               :: rEllipse, rdiff, eSize
        real(kind=RP)                               :: r, xx, yy, xs(2), rs(2), xEllipse, yEllipse
        real(kind=RP)                               :: xo, yo
        real(kind=RP), dimension(:,:), allocatable  :: rCalcsqr
        logical                                     :: onlyExternal

        if (present(filter)) then
            onlyExternal = filter
        else
            onlyExternal = .false.
        end if 

        isIn = .false.
        xo = centerPosition(1)
        yo = centerPosition(2)
        allocate(rCalcsqr(N,3))
        k = 1
        do i = 1, N
            rEllipse = POW2(x(i) - xo) + POW2( (y(i) - yo) / ratio )
            if ( POW2(mayorAxis) .lt. rEllipse ) return
            rCalcsqr(k,1) = x(i)
            rCalcsqr(k,2) = y(i)
            rCalcsqr(k,3) = POW2(x(i) - xo) + POW2(y(i) - yo)
            k = k + 1
        end do
        if (onlyExternal) then
            k = maxloc(rCalcsqr(:,3),dim=1)
            xx = rCalcsqr(k,1)
            yy = rCalcsqr(k,2)
            rdiff = distanceToEllipse(mayoraxis, ratio, centerposition, xx, yy)
            eSize = getElementSize(x, y, centerPosition, N)
            if ( abs(rdiff) .gt. eSize*lengthAspect ) return
        end if
        isIn = .true.

    End Function isInEllipse
!
!////////////////////////////////////////////////////////////////////////
!
    Function isInCircle(r, ratio, centerPosition, x, y, N, lengthAspect, filter) result(isIn)

        real(kind=RP)               , intent(in)    :: r, ratio, lengthAspect
        real(kind=RP), dimension(N) , intent(in)    :: x, y
        real(kind=RP), dimension(2) , intent(in)    :: centerPosition
        integer                     , intent(in)    :: N
        logical, optional           , intent(in)    :: filter
        logical                                     :: isIn

        ! local variables
        integer                                     :: i, j, k
        real(kind=RP)                               :: rdiff, eSize
        real(kind=RP)                               :: xo, yo
        real(kind=RP), dimension(:), allocatable    :: rCalcsqr
        logical                                     :: onlyExternal

        if (present(filter)) then
            onlyExternal = filter
        else
            onlyExternal = .false.
        end if 

        isIn = .false.
        xo = centerPosition(1)
        yo = centerPosition(2)
        allocate(rCalcsqr(N))
        k = 1
        do i = 1, N
            rCalcsqr(k) = POW2(x(i) - xo) + POW2(y(i) - yo)
            if ( POW2(r) .lt. rCalcsqr(k) ) return
            k = k + 1
        end do
        if (onlyExternal) then
            rdiff = r - sqrt(maxval(rCalcsqr))
            eSize = sqrt(POW2(x(N) - x(1)) + POW2(y(N) - y(1)))
            if ( abs(rdiff) .gt. eSize*lengthAspect ) return
        end if
        isIn = .true.

    End Function isInCircle

!
!////////////////////////////////////////////////////////////////////////
! distance to geometry functions
!////////////////////////////////////////////////////////////////////////
!
    Function distanceToEllipse(mayorAxis, ratio, centerPosition, x, y) result(rdiff)

        use Utilities, only: AlmostEqual
        implicit none

        real(kind=RP), intent(in)                   :: mayorAxis, ratio, x, y
        real(kind=RP), dimension(2), intent(in)     :: centerPosition
        real(kind=RP)                               :: rdiff

        ! local variables
        integer                                     :: k
        real(kind=RP)                               :: a, b, c, m, yint ! line equation and quadratic coefficients
        real(kind=RP)                               :: xo, yo
        real(kind=RP)                               :: r, xs(2), rs(2), xEllipse, yEllipse, rPoint

        xo = centerPosition(1)
        yo = centerPosition(2)

        vertical_case: if (AlmostEqual(x, xo)) then
            a = 1
            b = -2*yo
            c = POW2(yo) - POW2(mayorAxis * ratio)
            xs(1) = (-b + sqrt(POW2(b) - 4*a*c))/(2.0_RP*a)
            xs(2) = (-b - sqrt(POW2(b) - 4*a*c))/(2.0_RP*a)
            rs(1) = abs(y - xs(1))
            rs(2) = abs(y - xs(2))
            k = minloc(rs,dim=1)
            yEllipse = xs(k)
            r = abs(yEllipse - yo)
            rPoint = abs(y - yo)
        else vertical_case
            ! get line equation from center to corner of element
            m = (y - yo) / (x - xo)
            yint = yo - m*xo
            ! solve quadratic equation of intersection between ellipse and line
            a = 1 + POW2(m/ratio)
            b = -2*xo + 2*m*(yint-yo)/POW2(ratio)
            c = POW2((yint-yo)/ratio) + POW2(xo) - POW2(mayorAxis)
            xs(1) = (-b + sqrt(POW2(b) - 4*a*c))/(2.0_RP*a)
            xs(2) = (-b - sqrt(POW2(b) - 4*a*c))/(2.0_RP*a)
            rs(1) = abs(x - xs(1))
            rs(2) = abs(x - xs(2))
            k = minloc(rs,dim=1)
            xEllipse = xs(k)
            yEllipse = m*xEllipse + yint
            ! get the distance from center to intersection
            r= sqrt(POW2(xEllipse - xo) + POW2(yEllipse - yo))
            rPoint = sqrt(POW2(x - xo) + POW2(y - yo))
        end if vertical_case 

        rdiff = rPoint - r

    End Function distanceToEllipse
!
!////////////////////////////////////////////////////////////////////////
!
    Function distaceToCircle(r, ratio, centerPosition, x, y) result(rdiff)

        use Utilities, only: AlmostEqual
        implicit none

        real(kind=RP), intent(in)                   :: r, ratio, x, y
        real(kind=RP), dimension(2), intent(in)     :: centerPosition
        real(kind=RP)                               :: rdiff

        ! local variables
        real(kind=RP)                               :: xo, yo, rPoint

        xo = centerPosition(1)
        yo = centerPosition(2)

        rPoint = sqrt(POW2(x - xo) + POW2(y - yo))
        rdiff = rPoint - r

    End Function distaceToCircle
!
!////////////////////////////////////////////////////////////////////////
! extra helper functions
!////////////////////////////////////////////////////////////////////////
!
    Function getElementSize(x, y, centerPosition, N) result(eS)

        real(kind=RP), dimension(N), intent(in)     :: x, y
        real(kind=RP), dimension(2), intent(in)     :: centerPosition
        integer , intent(in)                        :: N
        real(kind=RP)                               :: eS

        ! local variables
        integer                                     :: i, j, k
        real(kind=RP)                               :: r, dx, dy, eps
        real(kind=RP)                               :: xo, yo

        ! eps = 0.15_RP
        ! eps = 0.25_RP
        eps = 0.3_RP
        xo = centerPosition(1)
        yo = centerPosition(2)
        r = POW2(x(1) - xo) + POW2(y(1) - yo)
        dx = (abs(x(1) - xo)) / r
        dy = (abs(y(1) - yo)) / r

        if (dy .le. eps) then
            eS = maxval(x) - minval(x)
        elseif (dx .le. eps) then
            eS = maxval(y) - minval(y)
        else
            eS = sqrt( POW2(maxval(x) - minval(x)) + POW2(maxval(y) - minval(y)) )
        end if 

    End Function getElementSize
!
!////////////////////////////////////////////////////////////////////////
!
End Module FWHPreSurface