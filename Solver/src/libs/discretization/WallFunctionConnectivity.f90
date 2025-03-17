!
!//////////////////////////////////////////////////////
!
! This module stores the connection of each face of the wall that will be used in the 
! Wall Function with the first normal neighbour
! element and other needed variables
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
Module WallFunctionConnectivity  !

    use SMConstants
    use HexMeshClass
    Implicit None
!   
!  *****************************
!  Default everything to private
!  *****************************
!
    private
!
!  ****************
!  Public variables
!  ****************
!
    public useWallFunc, wallFunBCs
!
!  ******************
!  Public definitions
!  ******************
!
    public Initialize_WallConnection, WallFunctionGatherFlowVariables, WallUpdateMeanV, WallStartMeanV
!

    logical                                                          :: useWallFunc
    character(len=BC_STRING_LENGTH), dimension(:), allocatable       :: wallFunBCs
    integer, dimension(:), allocatable                               :: wallElemIds, wallNormalDirection, wallNormalIndex, wallFaceID
    real(kind=RP), dimension(:,:,:,:), allocatable                   :: meanVelocity
    real(kind=RP)                                                    :: timeCont

    contains 
!   
!------------------------------------------------------------------------------------------------------------------------
!
    Subroutine Initialize_WallConnection(controlVariables, mesh)
!     *******************************************************************
!        This subroutine initializes and store the arrays index that are
!        used to represent the connection of the face to the neighbour element.
!        Also calls the definitions of the wall function to update the 
!        parameter of the models based on the controlVariables.
!     *******************************************************************
!
        use FTValueDictionaryClass
        use Utilities,            only: toLower
        use FileReadingUtilities, only: getCharArrayFromString
        use ElementConnectivityDefinitions, only: normalAxis, FACES_PER_ELEMENT
        use MPI_Process_Info
        use WallFunctionDefinitions, only: Initialize_Wall_Function, wallFuncIndex, STD_WALL, ABL_WALL, u_tau0, useAverageV
        use Headers
#ifdef _HAS_MPI_
       use mpi
#endif
        implicit none
        class(FTValueDictionary),  intent(in)  :: controlVariables
        class(HexMesh), intent(inout)          :: mesh

!       ---------------
!       Local variables
!       ---------------
!
        character(len=LINE_LENGTH)             :: wallBC_str
        integer                                :: numberFacesWall, numberBC
        integer, dimension(:), allocatable     :: zonesWall
        integer                                :: i, j, nz, k, fID, efID, ff
        integer                                :: actualElementID, linkedElementID, normalIndex, oppositeIndex, oppositeNormalIndex
        integer                                :: allFaces, ierr

        call Initialize_Wall_Function(controlVariables, useWallFunc)
        if (.not. useWallFunc) then
            return
        end if

        ! get BC where the Wall Function will be applied
        wallBC_str = controlVariables % stringValueForKey("wall function boundaries", LINE_LENGTH)
        call toLower(wallBC_str)
        call getCharArrayFromString(wallBC_str, BC_STRING_LENGTH, wallFunBCs)

        ! get zones and number of faces for wall function
        numberBC = size(wallFunBCs)
        numberFacesWall = 0
        allocate(zonesWall(numberBC))

        do i = 1, numberBC
            do nz = 1, size(mesh % zones)
                if (trim(mesh % zones(nz) % Name) .eq. trim(wallFunBCs(i))) then
                    zonesWall(i) = nz
                    numberFacesWall = numberFacesWall + mesh % zones(nz) % no_of_faces
                    exit
                end if
            end do
        end do

        if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
            call mpi_allreduce(numberFacesWall, allFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            allFaces = numberFacesWall
        end if
        if (allFaces .eq. 0) then
            useWallFunc = .false.
            if (MPI_Process % isRoot) write(STD_OUT,'(A)') "No wall BC found, the wall function will be deactivated"
            return
        end if

        allocate( wallElemIds(numberFacesWall), wallNormalIndex(numberFacesWall), wallNormalDirection(numberFacesWall), &
                 wallFaceID(numberFacesWall) )

        !get for each face of the wall, the linked element, normalDirection and index
        k = 0
        do j = 1, numberBC
            nz = zonesWall(j)
            do i = 1, mesh % zones(nz) % no_of_faces
                k = k + 1
                fID = mesh % zones(nz) % faces(i)
                actualElementID = mesh % faces(fID) % ElementIDs(1)
                associate ( e => mesh % elements(actualElementID) )
                    elem_loop:do efID = 1, FACES_PER_ELEMENT
                        if ( trim(wallFunBCs(j)) .eq. trim(e % boundaryName(efID)) ) then
                            normalIndex = normalAxis(efID)
                            exit elem_loop
                        end if
                    end do elem_loop
                    oppositeIndex = -1 * normalIndex
                    ! use the maxloc line if the compiler doesn't support findloc
                    ! oppositeIndex = findloc(normalAxis,oppositeIndex,dim=1)
                    oppositeIndex = maxloc(merge(1.0, 0.0, normalAxis == oppositeIndex),dim=1)
                    linkedElementID = e % Connection(oppositeIndex) % globID
                end associate
                linkedElementID = getElementID(mesh, linkedElementID)
                wallFaceID(k) = fID
                wallElemIds(k) = linkedElementID
                !get the normalIndex of the linked element instead of actual element, needed for rotated meshes
                do ff = 1, FACES_PER_ELEMENT
                    if (mesh % elements(linkedElementID) % Connection(ff) % globID .eq. mesh % elements(actualElementID) % globID) then
                        oppositeNormalIndex = normalAxis(ff)
                        exit
                    end if 
                end do
                wallNormalDirection(k) = abs(oppositeNormalIndex)
                wallNormalIndex(k) = getNormalIndex(mesh % faces(fID), mesh % elements(linkedElementID), wallNormalDirection(k))
            end do
        end do

!       Initialize u_tau storage for the wall faces
!       -------------------------------------------
        do i = 1, numberFacesWall
            fID = wallFaceID(i)
            associate( f => mesh%faces(fID) )
                f % storage(1) % u_tau_NS = u_tau0
            end associate
        end do 

!       Allocate average velocity if requested, only for same p in all faces
!       -------------------------------------------
        if (useAverageV) then
            fID = wallFaceID(1)
            associate( f => mesh%faces(fID) )
                allocate( meanVelocity(NDIM,numberFacesWall,0:f % Nf(1),0:f % nf(2)) )
            end associate
            ! meanVelocity = 0.0_RP
            call WallStartMeanV(mesh)
            timeCont = 0.0_RP
        end if

!       Describe the Wall function
!       --------------------------
        if ( .not. MPI_Process % isRoot ) return
        write(STD_OUT,'(/)')
        call Subsection_Header("Wall function")

        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of faces: ", allFaces
        select case (wallFuncIndex)
            case (STD_WALL)
                write(STD_OUT,'(30X,A,A28,A)') "->", "Wall Function Law: ", "Reichardt"
            case (ABL_WALL)
                write(STD_OUT,'(30X,A,A28,A10)') "->", "Wall Function Law: ", "ABL"
        end select
        if (useAverageV) write(STD_OUT,'(30X,A,A28)') "->", "Wall Function using time averaged values"

    End Subroutine Initialize_WallConnection

    integer Function getNormalIndex(f, e, normalDirection)
        use FaceClass
        use ElementClass
        implicit none

        class(Element), intent(in)                                       :: e
        class(Face), intent(in)                                          :: f
        integer, intent(in)                                              :: normalDirection

!       ---------------
!       Local variables
!       ---------------
!
        real(kind=RP), dimension(NDIM)                                   :: x0, xN, xf
        real(kind=RP)                                                    :: dx(2)
        integer                                                          :: N, indexArray(2), minIndex

        xf = f % geom % x(:,0,0)
        x0 = e % geom % x(:,0,0,0)
        select case (normalDirection)
        case (1)
            N = e % Nxyz(1)
            xN = e % geom % x(:,N,0,0)
        case (2)
            N = e % Nxyz(2)
            xN = e % geom % x(:,0,N,0)
        case (3)
            N = e % Nxyz(3)
            xN = e % geom % x(:,0,0,N)
        case default
           write(STD_OUT,'(A)') "Error: wallNormalDirection not found in axisMap"
           errorMessage(STD_OUT)
           error stop 
        end select

        indexArray = [0,N]
        dx(1) = norm2(xf-x0)
        dx(2) = norm2(xf-xN)
        minIndex = minloc(dx,dim=1)

        getNormalIndex = indexArray(minIndex)

    End Function getNormalIndex

    Subroutine WallFunctionGatherFlowVariables(mesh, f, V, rho, mu, dWall, Vavg)

!     *******************************************************************
!        This subroutine get the flow variables of the neighbour element
!        of the face.
!     *******************************************************************
!
        use PhysicsStorage
        use FaceClass
        use VariableConversion, only: get_laminar_mu_kappa
        use WallFunctionDefinitions, only: useAverageV
        implicit none

        class(HexMesh), intent(in)                                       :: mesh
        class(Face), intent(in)                                          :: f
        real(kind=RP), dimension(NDIM,0:f%Nf(1),0:f%Nf(2)), intent(out)  :: V
        real(kind=RP), dimension(0:f%Nf(1),0:f%Nf(2)), intent(out)       :: rho
        real(kind=RP), dimension(0:f%Nf(1),0:f%Nf(2)), intent(out)       :: mu
        real(kind=RP), dimension(0:f%Nf(1),0:f%Nf(2)), intent(out)       :: dWall
        real(kind=RP), dimension(NDIM,0:f%Nf(1),0:f%Nf(2)), intent(out)  :: Vavg

!       ---------------
!       Local variables
!       ---------------
!
        ! real(kind=RP), dimension(NCONS,0:f%Nf(1),0:f%Nf(2))             :: Q
        real(kind=RP), dimension(:,:,:), allocatable                     :: Q, x
        ! real(kind=RP), dimension(NDIM,0:f%Nf(1),0:f%Nf(2))              :: x
        real(kind=RP), dimension(NDIM)                                  :: dWallVector
        real(kind=RP)                                                   :: kappa, invRho
        integer                                                         :: i, j, fInd

        call WallGetFaceConnectedQ(mesh, f, Q, x, fInd)
        do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
                rho(i,j) = Q(IRHO,i,j)
                invRho = 1.0_RP / rho(i,j)
                V(:,i,j) = Q(IRHOU:IRHOW,i,j) * invRho
                call get_laminar_mu_kappa(Q(:,i,j),mu(i,j),kappa)
                dWallVector(:) = x(:,i,j) - f % geom % x(:,i,j)
                dWall(i,j) = norm2(dWallVector)
            end do
        end do
        
        if (useAverageV) then
            Vavg(:,:,:) = meanVelocity(:,fInd,:,:)
            ! Vavg = reshape( meanVelocity(:,fInd,0:f%Nf(1),0:f%Nf(2)), /NDIM,0:f%Nf(1),0:f%Nf(2)/ )
        else
            Vavg = 0.0_RP
        end if 

    End Subroutine WallFunctionGatherFlowVariables

    Subroutine WallGetFaceConnectedQ(mesh,f,Q,x,faceIndex)
!     *******************************************************************
!        This subroutine get the flow solution of the neighbour element
!        of the face.
!     *******************************************************************
        use PhysicsStorage
        use FaceClass
        type(HexMesh), intent(in)                                        :: mesh
        class(Face), intent(in)                                          :: f
        real(kind=RP), dimension(:,:,:), allocatable, intent(out)        :: Q
        real(kind=RP), dimension(:,:,:), allocatable, intent(out)        :: x
        integer, intent(out)                                             :: faceIndex
!       Local variables
        integer                                                         :: eID, solIndex
        ! integer                                                         :: faceIndex, eID, solIndex

        allocate( Q(NCONS,0:f % Nf(1),0:f % Nf(2)), x(NDIM,0:f % Nf(1),0:f % Nf(2)) )

        ! use the maxloc line if the compiler doesn't support findloc
        ! faceIndex = findloc(wallFaceID, f % ID, dim=1)
        faceIndex = maxloc(merge(1.0, 0.0, wallFaceID == f % ID),dim=1)
        eID = wallElemIds(faceIndex)
        solIndex = wallNormalIndex(faceIndex)

        ! select the local coordinate directions of the face base on the definition of axisMap in HexElementConnectivityDefinitions
        associate ( e => mesh % elements(eID) )
            select case (wallNormalDirection(faceIndex))
            case (1)
                Q(:,:,:) = e % storage % Q(:,solIndex,:,:)
                x(:,:,:) = e % geom % x(:,solIndex,:,:)
            case (2)
                Q(:,:,:) = e % storage % Q(:,:,solIndex,:)
                x(:,:,:) = e % geom % x(:,:,solIndex,:)
            case (3)
                Q(:,:,:) = e % storage % Q(:,:,:,solIndex)
                x(:,:,:) = e % geom % x(:,:,:,solIndex)
            case default
               write(STD_OUT,'(A)') "Error: wallNormalDirection not found in axisMap"
               errorMessage(STD_OUT)
               error stop 
            end select
        end associate

    End Subroutine WallGetFaceConnectedQ

    Subroutine WallUpdateMeanV(mesh, dt)
        use PhysicsStorage
        use WallFunctionDefinitions, only: useAverageV
        implicit none
        type(HexMesh), intent(in)                       :: mesh
        real(kind=RP),  intent(in)                      :: dt
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                                         :: fIndex, fID, i, j, fInd
        real(kind=RP), dimension(:,:,:), allocatable    :: Q, x
        real(kind=RP)                                   :: invRho
        real(kind=RP), dimension(NDIM)                  :: localV
        logical                                         :: saveLocal

        ! create separate to set initial conditions
        if (.not. useAverageV) return

        do fIndex = 1, size(wallFaceID)
            fID = wallFaceID(fIndex)
            associate( f => mesh%faces(fID) )
                call WallGetFaceConnectedQ(mesh, f, Q, x, fInd)
                do j = 0, f % Nf(2)
                    do i = 0, f % Nf(1)
                        invRho = 1.0_RP / Q(IRHO,i,j)
                        localV(:) = Q(IRHOU:IRHOW,i,j) * invRho
                        meanVelocity(:,fIndex,i,j) = ( (meanVelocity(:,fIndex,i,j) * timeCont) + localV(:) * dt ) / (timeCont+dt)
                    end do
                end do
            end associate
        end do 

        timeCont = timeCont + dt

    End Subroutine WallUpdateMeanV

    Subroutine WallStartMeanV(mesh)
        use PhysicsStorage
        use WallFunctionDefinitions, only: useAverageV
        implicit none
        type(HexMesh), intent(in)                       :: mesh
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                                         :: fIndex, fID, i, j, fInd
        real(kind=RP), dimension(:,:,:), allocatable    :: Q, x
        real(kind=RP)                                   :: invRho
        real(kind=RP), dimension(NDIM)                  :: localV
        logical                                         :: saveLocal

        ! create separate to set initial conditions
        if (.not. useAverageV) return

        do fIndex = 1, size(wallFaceID)
            fID = wallFaceID(fIndex)
            associate( f => mesh%faces(fID) )
                call WallGetFaceConnectedQ(mesh, f, Q, x, fInd)
                do j = 0, f % Nf(2)
                    do i = 0, f % Nf(1)
                        invRho = 1.0_RP / Q(IRHO,i,j)
                        localV(:) = Q(IRHOU:IRHOW,i,j) * invRho
                        meanVelocity(:,fIndex,i,j) = localV(:)
                    end do
                end do
            end associate
        end do 

    End Subroutine WallStartMeanV

    integer Function getElementID(mesh, globeID)

!     *******************************************************************
!        This subroutine get the element ID based on the global element ID
!        needed specially for MPI, to be sure that the element is in the mesh
!        partition.
!     *******************************************************************
!
        implicit none
        class(HexMesh), intent(in)             :: mesh
        integer, intent(in)                    :: globeID
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                                :: eID, solIndex

        do eID = 1, mesh % no_of_elements
            if (mesh % elements(eID) % globID .eq. globeID) then
                getElementID = eID
                return
            end if
        end do

        ! if the code reach this point the element does not exists
        write(STD_OUT,'(A,I0)') "Error: The element of the wall function does not exists in the mesh or mesh partition, global ID: ", globeID
        errorMessage(STD_OUT)
        error stop 

    End Function getElementID

End Module WallFunctionConnectivity
#endif