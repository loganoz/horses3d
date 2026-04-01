#include "Includes.h"
module LambVectorMeshInterpolation
    use SMConstants
    use omp_lib

    implicit none

    private 
    public InterpolateLambVector

    type LambElementStorage_t
        real(rp), allocatable :: Lamb(:,:,:,:) ! Lamb vector
    end type LambElementStorage_t

    logical :: saveGradients = .false.

    contains

    subroutine InterpolateLambVector(controlVariables)
        use FTValueDictionaryClass
        use HexMeshClass
        use Physics_CAAKeywordsModule
        use FileReadingUtilities
        implicit None
        TYPE(FTValueDictionary)                :: controlVariables

        ! Local variables
        type(HexMesh)                           :: M_in, M_out
        integer                             :: Nmax_in, Nmax_out, Nmax

        character(len=LINE_LENGTH) :: dir, basename, auxstr, meshFilename
        character(len=LINE_LENGTH), allocatable :: filenames(:)
        character(len=LINE_LENGTH), allocatable :: itera_i(:)
        integer :: i, pos, nExistingMeshFiles
        logical :: meshFileExists
        logical :: onlyOneMesh

        nExistingMeshFiles = 0
        
        !
        ! Find all the different files
        !
        ! Get all the Lamb vector files
        dir = controlVariables % stringValueForKey(LambReadDirectoryKey, LINE_LENGTH)
        basename = controlVariables % stringValueForKey(LambVectorFileNameKey, LINE_LENGTH)
        call getFilenamesAndTimes(dir, basename, filenames)
        ! Get iterations
        allocate(itera_i(size(filenames)))
        do i = 1, size(filenames)
            auxstr = getFileName(getFileName(RemovePath(filenames(i))))
            pos = index(auxstr,'_',BACK=.true.)
            itera_i(i) = trim(auxstr(pos+1:))
            ! Check if the mesh for such iteration exists
            auxstr = controlVariables % stringValueForKey("stats mesh file name", requestedLength = LINE_LENGTH)
            meshFilename = trim(RemovePath(getFileName(auxstr))) // "_" // trim(itera_i(i)) // ".hmesh"
            inquire(file="MESH/" // trim(meshFilename), exist=meshFileExists)
            if (.not. meshFileExists) then
                print *, "File ", trim(meshFilename), " does not exist."
            else
                nExistingMeshFiles = nExistingMeshFiles + 1
            end if
        end do


        if (nExistingMeshFiles == 1) then
            onlyOneMesh = .true.
            print *, "WARNING!!"
            print *, "Only a single mesh file was found in the MESH directory."
            print *, "All the Lamb vector data are interpolated with respect to this mesh."
        end if
        if (nExistingMeshFiles < size(filenames)) then
            print *, "ERROR: There are more Lamb vector files than mesh files."
            error stop
        end if        
        
        !
        ! Interpolate fields for all the time instants
        !
        if (onlyOneMesh) then
            auxstr = controlVariables % stringValueForKey("stats mesh file name", requestedLength = LINE_LENGTH)
            meshFilename = "MESH/" // trim(RemovePath(getFileName(auxstr))) // "_" // trim(itera_i(1)) // ".hmesh"
            call readInputOutputMesh(controlVariables, meshFilename, M_in, M_out, Nmax_in, Nmax_out)
        end if
        do i = 1, size(filenames)
            if (.not. onlyOneMesh) then
                auxstr = controlVariables % stringValueForKey("stats mesh file name", requestedLength = LINE_LENGTH)
                meshFilename = "MESH/" // trim(RemovePath(getFileName(auxstr))) // "_" // trim(itera_i(i)) // ".hmesh"
                call readInputOutputMesh(controlVariables, meshFilename, M_in, M_out, Nmax_in, Nmax_out)
            end if
            call InterpolateLamb_Mesh(controlVariables, M_in, M_out, Nmax_in, Nmax_out, filenames(i))
            if (.not. onlyOneMesh) then
                ! Destruct M_in
                safedeallocate(M_in % elements)
                safedeallocate(M_in % nodes)
                safedeallocate(M_in % Nx)
                safedeallocate(M_in % Ny)
                safedeallocate(M_in % Nz)
                ! Destruct M_out
                safedeallocate(M_out % elements)
                safedeallocate(M_out % nodes)
                safedeallocate(M_out % Nx)
                safedeallocate(M_out % Ny)
                safedeallocate(M_out % Nz)
            end if
        end do
    end subroutine InterpolateLambVector

    subroutine readInputOutputMesh(controlVariables, meshInhmesh_name, M_in, M_out, Nmax_in, Nmax_out)
        use FTValueDictionaryClass
        use HexMeshClass
        use pAdaptationClass          , only: GetMeshPolynomialOrders
        implicit None
        TYPE(FTValueDictionary)                :: controlVariables
        character(len=*), intent(in)            :: meshInhmesh_name
        class(HexMesh)                          :: M_in, M_out
        integer                             :: Nmax_in, Nmax_out


        !
        ! Local variables
        !
        type(HexMesh)                           :: M_in_aux
        CHARACTER(LEN=LINE_LENGTH)              :: fileName
        integer :: nodeType_in, nodeType_out
        integer, allocatable                :: Nx(:), Ny(:), Nz(:)
        integer                             :: eID

        !
        ! Read mesh M_in
        !
        ! First, read the hmesh file to get the nodes for each element
        ! fileName = controlVariables % stringValueForKey("stats mesh file name", requestedLength = LINE_LENGTH)
        print *, "Reading mesh M_in from hmesh: ", trim(meshInhmesh_name)
        call readMesh_FromHmeshFile(meshInhmesh_name, M_in_aux, nodeType_in, Nmax_in)
        print *, "DONE"
        ! Second, put the polynomial degree of each element in Nx, Ny, Nz array
        allocate(Nx(M_in_aux % no_of_elements))
        allocate(Ny(M_in_aux % no_of_elements))
        allocate(Nz(M_in_aux % no_of_elements))
        do eID = 1, M_in_aux % no_of_elements
            Nx(eID) = M_in_aux % elements(eID) % Nxyz(1)
            Ny(eID) = M_in_aux % elements(eID) % Nxyz(2)
            Nz(eID) = M_in_aux % elements(eID) % Nxyz(3)
        end do
        ! Second, read the msh file to construct all the geometry related data
        fileName = controlVariables % stringValueForKey("stats geo mesh file name", requestedLength = LINE_LENGTH)
        print *, "Reading mesh M_in from geo: ", trim(fileName)
        call ConstructSimpleMesh(M_in, fileName, AllNx=Nx, AllNy=Ny, AllNz=Nz)
        print *, "DONE"
        ! Assign data from hmesh to M_in
        do eID = 1, M_in % no_of_elements
            M_in % elements(eID) % Nxyz = M_in_aux % elements(eID) % Nxyz
            M_in % elements(eID) % geom % x = M_in_aux % elements(eID) % geom % x
        end do
        M_in % nodeType = nodeType_in
        call M_in % PrepareForIO()
        deallocate(Nx)
        deallocate(Ny)
        deallocate(Nz)
        
        !
        ! Read mesh M_out
        !
        print *, "Reading polynomial order for M_out"
        call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax_out)
        print *, "DONE"
        fileName = controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH)
        print *, "Reading mesh M_out: ", trim(fileName)
        call ConstructSimpleMesh(M_out, fileName, AllNx=Nx, AllNy=Ny, AllNz=Nz)
        print *, "DONE"
        deallocate(Nx)
        deallocate(Ny)
        deallocate(Nz)
        ! Get the node type of the mesh M_out
        call getNodeType(controlVariables, nodetype_out)
        M_out % nodeType = nodeType_out

    end subroutine readInputOutputMesh

    subroutine InterpolateLamb_Mesh(controlVariables, M_in, M_out, Nmax_in, Nmax_out, lambFilename)
        ! Interpolate fields defined in mesh M_in onto the mesh M_out.
        use FTValueDictionaryClass
        use pAdaptationClass          , only: GetMeshPolynomialOrders
        use HexMeshClass
        use StorageClass
        use NodalStorageClass
        use ProbeClass
        use TransfiniteMapClass
        use SolutionFile, only: NO_OF_SAVED_REFS
        use FileReadingUtilities
        implicit None
        TYPE(FTValueDictionary)                 :: controlVariables
        class(HexMesh)                          :: M_in, M_out
        integer, intent(in)                     :: Nmax_in, Nmax_out
        character(len=*)                        :: lambFilename

        !
        ! Local variables
        !
        CHARACTER(LEN=LINE_LENGTH)              :: fileName
        integer :: nodeType_in, nodeType_out
        integer                             :: Nmax
        integer                             :: eID, i, j, k, ii, jj, kk

        type(NodalStorage_t), target, allocatable :: NodalStorage_in(:)
        type(NodalStorage_t), pointer :: NodalStorage_out(:)

        type(LambElementStorage_t), allocatable :: storage_e_in(:), storage_e_out(:) ! To save the stats
        integer :: no_of_stats
        logical :: allocateSoundVelocity
        type(TransfiniteHexMap), pointer :: hexMap, hex8Map, genHexMap
        type(Probe_t) :: probe
        real(kind=RP) :: refs(NO_OF_SAVED_REFS)
        integer :: iter
        real(rp) :: time
        
        ! Initialize the nodal storage
        if (M_in % nodeType .ne. M_out % nodeType) then
            allocate(NodalStorage_in(0:Nmax_in))
            allocate(NodalStorage_out(0:Nmax_out))
        else
            Nmax = max(Nmax_in, Nmax_out)
            allocate(NodalStorage_in(0:Nmax))
            NodalStorage_out => NodalStorage_in
        end if
        ! Construct nodal storage data for the input mesh
        do eID = 1, M_in % no_of_elements
            do k = 0, M_in % elements(eID) % Nxyz(3)
                call NodalStorage_in(k) % Construct( M_in % nodeType, k )
            end do
            do j = 0, M_in % elements(eID) % Nxyz(1)
                call NodalStorage_in(j) % Construct( M_in % nodeType, j )
            end do
            do i = 0, M_in % elements(eID) % Nxyz(1)
                call NodalStorage_in(i) % Construct( M_in % nodeType, i )
            end do
        end do
        ! Construct nodal storage data for the output mesh
        do eID = 1, M_out % no_of_elements
            do k = 0, M_out % Nz(eID)
                call NodalStorage_out(k) % Construct( M_out % nodeType, k )
            end do
            do j = 0, M_out % Ny(eID)
                call NodalStorage_out(j) % Construct( M_out % nodeType, j )
            end do
            do i = 0, M_out % Nx(eID)
                call NodalStorage_out(i) % Construct( M_out % nodeType, i )
            end do
        end do


        !
        ! Read stats in M_in
        !
        ! Initialize storage for mesh M_in
        allocate(storage_e_in(M_in % no_of_elements))
        fileName = controlVariables % stringValueForKey("stats qbase file name", requestedLength = LINE_LENGTH)
        call readLambVector(M_in, storage_e_in, lambFilename)


        ! Initialize storage for mesh M_out
        allocate(storage_e_out(M_out % no_of_elements))
        call allocateStorage(M_out, storage_e_out)

        !$omp parallel default(shared) private(hex8Map,genHexMap,hexMap)
        ! Allocate the transfinite mappings
        allocate(hex8Map)
        allocate(genHexMap)
        
        !$omp do schedule(runtime) private(eID,i,j,k,ii,jj,kk,probe)
        do eID = 1, M_out % no_of_elements

            ! Create transfinite mapping
            if (M_out % elements(eID) % SurfInfo % IsHex8) then
               call hex8Map % setCorners(M_out % elements(eID) % SurfInfo % corners)
               hexMap => hex8Map
            else
               CALL genHexMap % destruct()
               CALL genHexMap % constructWithFaces(M_out % elements(eID) % SurfInfo % facePatches)
               hexMap => genHexMap
            end if

            do k = 0, M_out % elements(eID) % Nxyz(3)
                do j = 0, M_out % elements(eID) % Nxyz(2)
                    do i = 0, M_out % elements(eID) % Nxyz(1)
                        !
                        ! Create the node as a probe
                        !
                        call assignProbeCoordinates(probe, NodalStorage_out, M_out % elements(eID) % Nxyz, hexMap, eID, i, j, k)

                        !
                        ! Find element of M_in containing the probe
                        !
                        ! This is needed because the global variable NodalStorage is used inside the function
                        NodalStorage => NodalStorage_in
                        call findProbeInMesh(M_in, probe)

                        !
                        ! Interpolation
                        !
                        ! Allocate the polynomials for interpolation
                        safedeallocate(probe % lxi  ) ; allocate( probe % lxi(0 : M_in % elements(probe % eID) % Nxyz(1)) )
                        safedeallocate(probe % leta ) ; allocate( probe % leta(0 : M_in % elements(probe % eID) % Nxyz(2)) )
                        safedeallocate(probe % lzeta) ; allocate( probe % lzeta(0 : M_in % elements(probe % eID) % Nxyz(3)) )
                        probe % lxi = NodalStorage_in(M_in % elements(probe % eID) % Nxyz(1)) % lj(probe % xi(1))
                        probe % leta = NodalStorage_in(M_in % elements(probe % eID) % Nxyz(2)) % lj(probe % xi(2))
                        probe % lzeta = NodalStorage_in(M_in % elements(probe % eID) % Nxyz(3)) % lj(probe % xi(3))
                        ! Compute high-order interpolation
                        do kk = 0, M_in % elements(probe % eID) % Nxyz(3)
                            do jj = 0, M_in % elements(probe % eID) % Nxyz(2)
                                do ii = 0, M_in % elements(probe % eID) % Nxyz(1)
                                    storage_e_out(eID) % Lamb(:,i,j,k) = storage_e_out(eID) % Lamb(:,i,j,k) + storage_e_in(probe % eID) % Lamb(:,ii,jj,kk) * probe % lxi(ii) * probe % leta(jj) * probe % lzeta(kk)
                                end do
                            end do
                        end do   

                    end do
                end do
            end do
        end do
        !$omp end do

        deallocate(hex8Map)
        deallocate(genHexMap)

        !$omp end parallel


        ! Export field values at probes
        call M_out % PrepareForIO()
        fileName = "./RESULTS/" // trim(removePath(getFileName(getFileName(lambFilename)))) // '.InterpolatedLamb.hsol'
        call saveLambVector(M_out, storage_e_out, refs, iter, time, trim(fileName))

    end subroutine InterpolateLamb_Mesh

    subroutine ConstructSimpleMesh(mesh, meshFileName, AllNx, AllNy, AllNz)
        use SMConstants
        use HexMeshClass
        use readHDF5
        use readSpecM
        use Read_GMSH, only: CheckGMSHversion
        use readGMSH
        use FileReadingUtilities, only: getFileExtension
        implicit none
        !-arguments----------------------------------------------
        type(HexMesh)               :: mesh
        character(len=*)            :: meshFileName
        integer, intent(in)         :: AllNx(:), AllNy(:), AllNz(:)
        !-local-variables----------------------------------------
        integer                    :: gmsh_version
        character(len=LINE_LENGTH) :: ext

        ext = getFileExtension(trim(meshFileName))
        if (trim(ext)=='h5') then
            call ConstructSimpleMesh_FromHDF5File_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
        elseif (trim(ext)=='mesh') then
            call ConstructSimpleMesh_FromSpecFile_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
        elseif (trim(ext)=='msh') then
            call CheckGMSHversion (meshFileName, gmsh_version)
            select case (gmsh_version)
            case (4)
                call ConstructSimpleMesh_FromGMSHFile_v4_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
            case (2)
                call ConstructSimpleMesh_FromGMSHFile_v2_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
            case default
                error stop "ReadMeshFile :: Unrecognized GMSH version."
            end select
        else
            error stop 'Mesh file extension not recognized.'
        end if
        mesh % no_of_allElements = mesh % no_of_elements

    end subroutine ConstructSimpleMesh

    subroutine readMesh_FromHmeshFile(meshName, mesh, nodeType, Nmax)
        ! Read mesh coordinates from .hmesh file. Purely sequential code.
        use HexMeshClass
        use SolutionFile
        implicit none
        character(len=*), intent(in)     :: meshName
        class(HexMesh) :: mesh
        integer, intent(out) :: nodeType
        integer, intent(out) :: Nmax

        ! Local variables
        integer, dimension(:), allocatable     :: arrayDimensions
        integer                                :: dimensionsSize
        integer                                :: fid, eID
        integer                                :: no_of_elements, npa, ndof
        integer                                :: lb, ub

        ! Get mesh node type and type of file
        nodeType = getSolutionFileNodeType(meshName)

        ! Get number of elements
        no_of_elements = getSolutionFileNoOfElements(meshName)
        mesh % no_of_elements = no_of_elements

        ! Allocate elements
        allocate(mesh % elements(no_of_elements))
        dimensionsSize = 4
        allocate(arrayDimensions(dimensionsSize))

        ! Read coordinates
        fid = putSolutionFileInReadDataMode(meshName)
        npa = 0
        do eID = 1, no_of_elements

            call getSolutionFileArrayDimensions(fid,arrayDimensions)

            ! Allocate memory for the coordinates
            ! arrayDimensions(2:dimensionsSize) = arrayDimensions(2:dimensionsSize) - 1
            mesh % elements(eID) % Nxyz(1:dimensionsSize-1) = arrayDimensions(2:dimensionsSize) - 1
            allocate( mesh % elements(eID) % geom % x(NDIM,0:mesh % elements(eID) % Nxyz(1),0:mesh % elements(eID) % Nxyz(2),0:mesh % elements(eID) % Nxyz(3)) )

            ! Read data
            read(fid) mesh % elements(eID) % geom % x

            ndof = (mesh % elements(eID) % Nxyz(1) + 1) * (mesh % elements(eID) % Nxyz(2) + 1) * (mesh % elements(eID) % Nxyz(3) + 1)
            npa = npa + ndof

            Nmax = max(Nmax, maxval(mesh % elements(eID) % Nxyz))

        end do

        ! Close file
        close(fid)

    end subroutine readMesh_FromHmeshFile

    subroutine getNodeType(controlVariables, nodetype)
        use mainKeywordsModule, only: discretizationNodesKey
        use NodalStorageClass, only: GAUSS, GAUSSLOBATTO
        use FTValueDictionaryClass
        implicit none
        TYPE(FTValueDictionary)                :: controlVariables
        integer, intent(out) :: nodeType

        if (.not. controlVariables % containsKey(discretizationNodesKey)) then
            call controlVariables % addValueForKey("Gauss",discretizationNodesKey)
        end if

        select case ( trim(controlVariables % stringValueForKey(trim(discretizationNodesKey), requestedLength = LINE_LENGTH)) )
            case("Gauss")
                nodeType = GAUSS
            case("Gauss-Lobatto")
                nodeType = GAUSSLOBATTO
            case default
                print*, "Unknown discretization nodes."
                print*, "Options available are:"
                print*, "   * Gauss"
                print*, "   * Gauss-Lobatto"
                errorMessage(STD_OUT)
                error stop
        end select
    end subroutine getNodeType

    subroutine write_points_vtk(filename, npa, xd_dp_pa, f_pa)
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        integer, intent(in) :: npa
        real(rp), intent(in) :: xd_dp_pa(3,npa)
        real(rp), optional, intent(in) :: f_pa(:,:)

        ! local
        integer :: i
        integer :: nf
        integer :: unit

        if (present(f_pa)) nf = size(f_pa,1)

        open(newunit=unit,file=filename,status='replace',action='write',form='formatted')

        ! header
        write(unit,'(A)') '# vtk DataFile Version 3.0'
        write(unit,'(A)') 'point data'
        write(unit,'(A)') 'ASCII'
        write(unit,'(A)') 'DATASET POLYDATA'

        ! points
        write(unit,'(A,I0,A)') 'POINTS ', npa, ' double'
        do i=1,npa
            write(unit,'(3ES24.16)') xd_dp_pa(1,i),xd_dp_pa(2,i),xd_dp_pa(3,i)
        end do

        ! vertices
        write(unit,'(A,I0,1X,I0)') 'VERTICES ', npa, 2*npa
        do i=0,npa-1
            write(unit,'(2I10)') 1,i
        end do

        if (present(f_pa)) then
            ! point data
            write(unit,'(A,I0)') 'POINT_DATA ', npa

            if (nf == 1) then

                write(unit,'(A)') 'SCALARS field double'
                write(unit,'(A)') 'LOOKUP_TABLE default'

                do i=1,npa
                    write(unit,'(ES24.16)') f_pa(1,i)
                end do

            elseif (nf == 3) then

                write(unit,'(A)') 'VECTORS field double'

                do i=1,npa
                    write(unit,'(3ES24.16)') f_pa(1,i),f_pa(2,i),f_pa(3,i)
                end do

            else

                write(*,*) 'Error: f_pa must have size 1 or 3 in first dimension'
                stop

            endif
        end if

        close(unit)

    end subroutine write_points_vtk

    subroutine allocateStorage(M_out, storage_e_out)
        use HexMeshClass
        implicit none
        type(HexMesh), intent(in)          :: M_out
        type(LambElementStorage_t)          :: storage_e_out(:)
        
        integer :: eID

        do eID = 1, M_out % no_of_elements
            ! Allocate Lamb vector
            allocate(storage_e_out(eID) % Lamb(1:NDIM, 0:M_out % elements(eID) % Nxyz(1), 0:M_out % elements(eID) % Nxyz(2), 0:M_out % elements(eID) % Nxyz(3)) )
            storage_e_out(eID) % Lamb = 0.0_rp
        end do
    end subroutine allocateStorage

    subroutine assignProbeCoordinates(probe, NodalStorage_out, Nxyz, hexMap, eID, i, j, k)
        use ProbeClass
        use NodalStorageClass
        use TransfiniteMapClass
        implicit none
        type(Probe_t) :: probe
        type(NodalStorage_t), pointer :: NodalStorage_out(:)
        integer, intent(in) :: Nxyz(NDIM)
        type(TransfiniteHexMap), intent(in) :: hexMap
        integer, intent(in) :: eID, i, j, k

        ! Local variables
        TYPE(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
        real(rp) :: x(NDIM)

        probe % ID = eID
        spAxi => NodalStorage_out(Nxyz(1))
        spAeta => NodalStorage_out(Nxyz(2))
        spAzeta => NodalStorage_out(Nxyz(3))
        ! Coordinates of the point in the computational space
        x = [spAxi % x(i), spAeta % x(j), spAzeta % x(k)]
        ! Coordinates of the point in the physical space
        probe % x = hexMap % transfiniteMapAt(x)

    end subroutine assignProbeCoordinates

    subroutine findProbeInMesh(M_in, probe)
        ! Find the requested point in the mesh
        use ProbeClass
        use NodalStorageClass
        use HexMeshClass
        implicit none
        type(HexMesh) :: M_in
        type(Probe_t) :: probe

        probe % active = M_in % FindPointWithCoords(probe % x, probe % eID, probe % xi)
        ! Check whether the probe is located in another partition
        call probe % LookInOtherPartitions
        ! Disable the probe if the point is not found
        if ( .not. probe % active ) then
            ! if ( MPI_Process % isRoot ) then
                write(STD_OUT,'(A,I0,A)') "Probe ", probe % ID, " was not successfully initialized."
                print*, "Probe is set to inactive."
            ! end if
            ! return
        end if
    end subroutine findProbeInMesh

    subroutine readLambVector(mesh, storage_e, filename)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit None
        CLASS(HexMesh)                  :: mesh
        type(LambElementStorage_t) :: storage_e(mesh % no_of_elements)
        character(len=*)                :: fileName

        !
        ! Local variables
        !
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos


        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        if ( (fileType .ne. SOLUTION_FILE) ) then
            print *, "The file ", fileName, " is not a solution file."
            errorMessage(STD_OUT)
            error stop
        end if
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
                "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % Lamb(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize sound velocity
                read(fID, pos=pos) storage_e(eID) % Lamb(:,:,:,:)

            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readLambVector

    subroutine saveLambVector(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(LambElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        ! Local variables
        integer                          :: fid, eID
        integer                          :: no_stat_s
        integer(kind=AddrInt)            :: pos 
        real(kind=RP), allocatable       :: Q(:,:,:,:)

        ! Create new file
        call CreateNewSolutionFile(trim(name),SOLUTION_FILE, mesh % nodeType, mesh % no_of_Allelements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        do eID = 1, mesh % no_of_elements
            pos = POS_INIT_DATA + (mesh % elements(eID) % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*mesh % elements(eID) % offsetIO*SIZEOF_RP
            call writeArray(fid, storage_e(eID) % Lamb, position=pos)
        end do
        close(fID)

        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveLambVector

    subroutine getFilenamesAndTimes(dir, basename, filenames)
        use FileReadingUtilities, only: readFilesByGlob
        use FTValueDictionaryClass
        implicit none
        character(len=*), intent(in) :: dir
        character(len=*), intent(in) :: basename
        character(len=*), allocatable :: filenames(:)

        ! Read the filenames matching pattern
        call readFilesByGlob(dir, basename, filenames)

    end subroutine

end module LambVectorMeshInterpolation