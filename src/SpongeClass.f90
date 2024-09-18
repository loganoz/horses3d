!//////////////////////////////////////////////////////
!
!This class represents the numerical source term as a sponge to complement BC

#include "Includes.h"
Module SpongeClass  !
#if defined(NAVIERSTOKES)
    use SMConstants
    use HexMeshClass
    use SolutionFile
    use PhysicsStorage_NS, only: NCONS
    Implicit None

    private
    public sponge

    !definition of sponge class 
    type sponge_t
        integer                                                  :: nElements        ! number of elements in sponge region
        integer                                                  :: nElementsAll     ! number of elements in sponge region in all partitions
        integer, dimension(:), allocatable                       :: elementIndexMap  ! map from eID of mesh to sponge arrays
        real(kind=RP), dimension(:,:,:,:), allocatable           :: intensity        ! Intensity of the sponge in the domain, includes the amplitude and the ramp, precomputed for all elements in sponge region
        real(kind=RP)                                            :: amplitude        ! amplitude of the source term
        real(kind=RP)                                            :: delta            ! temporal filter width
        real(kind=RP)                                            :: rampWidth        ! length of the ramp width
        real(kind=RP), dimension(NDIM)                           :: x0               ! position of start of the sponge (for cylindrical in the center)
        real(kind=RP)                                            :: radious          ! radious of the ramp zone in cylindrical/cirular sponges
        real(kind=RP), dimension(:,:,:,:,:), allocatable         :: Qbase            ! Base flow (moving average or constant), for every variable in each GAUSS or GL node
        real(kind=RP), dimension(NDIM)                           :: axis             ! axis vector of the sponge. In cylindrical axis of the cylinder, in cartesian, the aligned vector
        character(len=STRING_CONSTANT_LENGTH)                    :: shapeType        ! shape of the sponge, either cartesian (aligned with an axis) or cylindrical
        character(len=STRING_CONSTANT_LENGTH)                    :: initialFileName  ! file name to load the initial base flow
        character(len=STRING_CONSTANT_LENGTH)                    :: solutionFileName ! file name to write base flow
        logical                                                  :: readBaseFLowFlag     ! read base flow from file or use instantaneous Q to start
        logical                                                  :: useMovingAverage ! to use Qbase as a moving average
        logical                                                  :: isActive = .false.

        contains

        procedure :: construct      => spongeConstruct
        procedure :: destruct       => spongeDestruct
        procedure :: creatRamp
        procedure :: addSource
        procedure :: initializeBaseFlow
        procedure :: updateBaseFlow
        procedure :: readBaseFlow
        procedure :: writeBaseFlow

    end type sponge_t

  type(sponge_t)                                                :: sponge

  contains

!/////////////////////////////////////////////////////////////////////////
!           SPONGE PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

    Subroutine spongeConstruct(self, mesh, controlVariables)
        use FileReadingUtilities, only: getRealArrayFromString
        use FTValueDictionaryClass
        use MPI_Process_Info
        use Headers
        use mainKeywordsModule
        use FileReadingUtilities, only: getFileName
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(in)                               :: mesh
        type(FTValueDictionary), intent(in)                     :: controlVariables

        !local variables
        character(len=STRING_CONSTANT_LENGTH)                   :: coordinates, axis, fileNameControl, solution_file, restart_name, restart_baseflow_name
        logical                                                 :: fileExists

        self % isActive = .false.
        if (.not. controlVariables % logicalValueForKey("use sponge")) return

        self % amplitude = controlVariables % getValueOrDefault("sponge amplitude",1.0_RP)
        self % rampWidth = controlVariables % getValueOrDefault("sponge ramp width",1.0_RP)
        self % delta = controlVariables % getValueOrDefault("sponge temporal width",1.0_RP)
        coordinates = controlVariables % stringValueForKey("sponge start position", requestedLength = STRING_CONSTANT_LENGTH)
        self % x0 = getRealArrayFromString(trim(coordinates))
        self % useMovingAverage = controlVariables % logicalValueForKey("sponge use moving average")
        self % radious = controlVariables % getValueOrDefault("sponge radious",0.0_RP)

        axis = controlVariables % stringValueForKey("sponge axis", requestedLength = STRING_CONSTANT_LENGTH)
        self % axis = getRealArrayFromString(trim(axis))
        ! self % axis = [0,0,1]
        ! for now only cylindrical sponges are available, todo: cartesian
        self % shapeType = trim(controlVariables % stringValueForKey("sponge shape", requestedLength = STRING_CONSTANT_LENGTH))

!       Get the file name of the solution
!       ---------------------------------
        solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STRING_CONSTANT_LENGTH )
        solution_file = trim(getFileName(solution_file))
        write(self % solutionFileName,'(2A)')  TRIM(solution_file),'_baseflow'

!       Get the file name for the base flow
!       -----------------------------------
        self % readBaseFLowFlag = .false.
        restart_name = controlVariables % stringValueForKey( restartFileNameKey, requestedLength = STRING_CONSTANT_LENGTH )
        restart_name = trim(getFileName(restart_name))
        write(restart_baseflow_name,'(2A)')  TRIM(restart_name),'_baseflow.hsol'
        inquire(file=trim(restart_baseflow_name), exist=fileExists)
        if (fileExists) then
            self % initialFileName = trim(restart_baseflow_name)
            self % readBaseFLowFlag = .true.
        else
            if (controlVariables % containsKey("sponge base flow file")) then
                self % initialFileName = trim(controlVariables % stringValueForKey("sponge base flow file", requestedLength = STRING_CONSTANT_LENGTH))
                self % readBaseFLowFlag = .true.
            end if
        end if 

        ! create arrays and pre calculate values
        call self % creatRamp(mesh)
        call self % initializeBaseFlow(mesh)

         self % isActive = .true.
         if ( .not. MPI_Process % isRoot ) return
         call Subsection_Header("Sponge")
         write(STD_OUT,'(30X,A,A28,A)') "->", "Shape: ", self % shapeType
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of elements: ", self % nElementsAll
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Amplitude: ", self % amplitude
         write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "Ramp start: ", self % x0
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Ramp width: ", self % rampWidth
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Ramp radious: ", self % radious
         write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "Axis: ", self % axis
         write(STD_OUT,'(30X,A,A28,L1)') "->", "Use moving average: ", self % useMovingAverage
         if (self % readBaseFLowFlag) write(STD_OUT,'(30X,A,A28,A)') "->", "Initial base flow file: ", self % initialFileName

    End Subroutine spongeConstruct
!
    Subroutine spongeDestruct(self)
        Implicit None
        class(sponge_t), intent(inout)                   :: self

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return
        safedeallocate(self % elementIndexMap)
        safedeallocate(self % intensity)
        safedeallocate(self % Qbase)

    End Subroutine spongeDestruct
!
    Subroutine creatRamp(self, mesh)
        use MPI_Process_Info
#ifdef _HAS_MPI_
        use mpi
#endif
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(in)                               :: mesh

        !local variables
        real(kind=RP), dimension(:,:,:), allocatable            :: xStar, sigma
        real(kind=RP), dimension(NDIM)                          :: r_vector
        logical, dimension(:), allocatable                      :: hasSponge
        integer                                                 :: i, j, k, eID, counter, spongeEID, ierr
        integer, dimension(NDIM)                                :: Nxyz

        ! it wont work for p-refinement or p adaption
        Nxyz = mesh % elements(1) % Nxyz

        allocate( xStar(0:Nxyz(1),0:Nxyz(2),0:Nxyz(3)), sigma(0:Nxyz(1),0:Nxyz(2),0:Nxyz(3)), hasSponge(mesh % no_of_elements) )
        hasSponge = .false.
        ! self % elementIndexMap = 0
        counter = 0
        do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
                do k = 0, Nxyz(3) ; do j = 0, Nxyz(2) ; do i = 0, Nxyz(1)
                    r_vector(:) = e % geom % x(:,i,j,k) - self % x0(:)
                    ! in this case xStar is actually rdiff ^2
                    xStar(i,j,k) = sum(r_vector*r_vector) - POW2(self % radious)
                end do         ; end do          ; end do
                if (any(xStar .ge. 0.0_RP)) then
                    hasSponge(eID) = .true.
                    counter = counter + 1
                end if 
            end associate
        end do

        self % nElements = counter
        if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
      call mpi_allreduce(self % nElements, self % nElementsAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            self % nElementsAll = self % nElements
        end if

        ! create mapping array
        allocate( self % intensity(Nxyz(1),Nxyz(2),Nxyz(3),self%nElements), self  %  elementIndexMap(self % nElements) )
        counter = 0
        do eID = 1, mesh % no_of_elements
            if (hasSponge(eID)) then
                counter = counter + 1
                self % elementIndexMap(counter) = eID
            end if 
        end do

        ! now calculate the ramp amplitude
        do spongeEID = 1, self % nElements
            eID = self % elementIndexMap(spongeEID)
            sigma = 0
            associate(e => mesh % elements(eID))
                do k = 0, Nxyz(3) ; do j = 0, Nxyz(2) ; do i = 0, Nxyz(1)
                    r_vector(:) = e % geom % x(:,i,j,k) - self % x0(:)
                    ! remove components in the axis direction
                    r_vector(:) = r_vector(:) - sum((e % geom % x(:,i,j,k) - self % x0(:))*self % axis(:))*self % axis(:)
                    ! xStar is non-dimensionalized by the width of the ramp, since the difference is squared, the width is too
                    xStar(i,j,k) = (sum(r_vector*r_vector) - POW2(self % radious)) / POW2(self % rampWidth)
                end do         ; end do          ; end do
            end associate
            ! limit xStar to [0,1] since after 1 should be constant at the amplitude value
            xStar = max(0.0_RP,xStar)
            xStar = min(1.0_RP,xStar)
            ! Sponge Ramping Function, taken from Beck, A., and Munz, C.-D., Direct Aeroacoustic Simulations Based on High Order Discontinuous Galerkin Schemes, Springer, Cham, 2018, Vol. 579
            sigma  = 6.0_RP*xStar**5.0_RP - 15.0_RP*xStar**4.0_RP + 10.0_RP*xStar**3.0_RP
            ! limit sigms <=1 since after 1 should be constant at the amplitude value
            sigma  = MIN(1.0_RP,sigma)
            ! pre computed intensity, including amplitude and ramp damping
            self % intensity(:,:,:,spongeEID) = self % amplitude * sigma(:,:,:)
        end do 

        deallocate(xStar, sigma, hasSponge)

    End Subroutine creatRamp
!
    Subroutine addSource(self,mesh)
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(inout)                            :: mesh
    
        !local variables
        integer                                                 :: i, j, k, eID, spongeEID
        integer, dimension(NDIM)                                :: Nxyz

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

        Nxyz = mesh % elements(1) % Nxyz

!$omp do schedule(runtime) private(i,j,k)
        do spongeEID = 1, self % nElements
            eID = self % elementIndexMap(spongeEID)
            associate(e => mesh % elements(eID))
                do k = 0, Nxyz(3) ; do j = 0, Nxyz(2) ; do i = 0, Nxyz(1)
                    e % storage % S_NS(:,i,j,k) = - self % intensity(i,j,k,spongeEID) * &
                                                  (e % storage % Q(:,i,j,k) - self % Qbase(:,i,j,k,eID))
                end do         ; end do          ; end do
            end associate
        end do
!$omp end do

    End Subroutine addSource
!
    Subroutine initializeBaseFlow(self,mesh)
        use ElementClass
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(in)                               :: mesh

        !local variables
        integer                                                 :: eID
        integer, dimension(NDIM)                                :: Nxyz

        Nxyz = mesh % elements(1) % Nxyz
        allocate(self % Qbase(NCONS, 0: Nxyz(1), 0: Nxyz(2), 0: Nxyz(3), mesh % no_of_elements))
        if (self % readBaseFLowFlag) then
            call self % readBaseFlow(mesh)
        else
!$omp do schedule(runtime)
            do eID = 1, mesh % no_of_elements
                associate(e => mesh % elements(eID))
                    self % Qbase(:,:,:,:,eID) = e % storage % Q(:,:,:,:) 
                end associate
            end do
!$omp end do
        end if 
!
    End Subroutine initializeBaseFlow
!
    ! advance in time base flow as single euler step, taken from Beck 2018
    Subroutine updateBaseFlow(self, mesh, dt)
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(inout)                            :: mesh
        real(kind=RP), intent(in)                               :: dt

        !local variables
        integer                                                 :: eID, spongeEID

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!       Only update for moving average
!       ------------------------
        if (.not. self % useMovingAverage) return

!$omp do schedule(runtime)
        do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
                self % Qbase(:,:,:,:,eID) = self % Qbase(:,:,:,:,eID) + (e % storage % Q(:,:,:,:) - self % Qbase(:,:,:,:,eID)) * (dt / self%delta)
            end associate
        end do
!$omp end do

    End Subroutine updateBaseFlow
!
    Subroutine  readBaseFlow(self,mesh)
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(in)                               :: mesh

        !local variables
        INTEGER                        :: fID, eID, fileType, no_of_elements, flag 
        integer                        :: padding, pos
        integer                        :: Nxp1, Nyp1, Nzp1, no_of_eqs, array_rank
        real(kind=RP), allocatable     :: Q(:,:,:,:)
        integer, dimension(NDIM)                                :: Nxyz

!       Get the file type
!       -----------------
        fileType = getSolutionFileType(trim(self % initialFileName))

        select case (fileType)
        case(SOLUTION_FILE)
           padding = 1*NCONS
        case default
           print*, "Only solution file format is accepted"
           errorMessage(STD_OUT)
           stop
        end select
!
!       Read the number of elements
!       ---------------------------
        no_of_elements = getSolutionFileNoOfElements(trim(self % initialFileName))

        if ( no_of_elements .ne. mesh % no_of_allElements ) then
           write(STD_OUT,'(A,A)') "The number of elements stored in the restart file ", &
                                  "do not match that of the mesh file"
           errorMessage(STD_OUT)
           stop
        end if
!
!       Read the terminator indicator
!       -----------------------------
        flag = getSolutionFileDataInitFlag(trim(self % initialFileName))

        if ( flag .ne. BEGINNING_DATA ) then
           print*, "Beginning data flag was not found in the file."
           errorMessage(STD_OUT)
           stop
        end if
!
!       Read elements data
!       ------------------
        fID = putSolutionFileInReadDataMode(trim(self % initialFileName))
        Nxyz = mesh % elements(1) % Nxyz
        allocate(Q(NCONS, 0:Nxyz(1), 0:Nxyz(2), 0:Nxyz(3)))
        do eID = 1, size(mesh % elements)
           associate( e => mesh % elements(eID) )
               pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + padding*e % offsetIO*SIZEOF_RP
               read(fID, pos=pos) array_rank
               read(fID) no_of_eqs, Nxp1, Nyp1, Nzp1
               if (      ((Nxp1-1) .ne. e % Nxyz(1)) &
                    .or. ((Nyp1-1) .ne. e % Nxyz(2)) &
                    .or. ((Nzp1-1) .ne. e % Nxyz(3)) &
                    .or. (no_of_eqs .ne. NCONS )       ) then
                  write(STD_OUT,'(A,I0,A)') "Error reading restart file: wrong dimension for element "&
                                              ,eID,"."

                  write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Element dimensions: ", e % Nxyz(1), &
                                                                        " ,", e % Nxyz(2), &
                                                                        " ,", e % Nxyz(3), &
                                                                        "."

                  write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Restart dimensions: ", Nxp1-1, &
                                                                        " ,", Nyp1-1, &
                                                                        " ,", Nzp1-1, &
                                                                        "."

                  errorMessage(STD_OUT)
                  stop
               end if

               read(fID) Q
               self % Qbase(:,:,:,:,eID) = Q(:,:,:,:)
          end associate
        end do
        deallocate(Q)
!
!       Close the file
!       --------------
        close(fID)

    End Subroutine  readBaseFlow
!
    Subroutine  writeBaseFlow(self,mesh,iter,time,last)
        use FluidData, only: thermodynamics, refValues, dimensionless
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(in)                               :: mesh
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        logical, intent(in), optional          :: last
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                          :: fid, eID
        integer(kind=AddrInt)            :: pos, padding
        real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
        real(kind=RP), allocatable       :: Q(:,:,:,:)
        integer, dimension(NDIM)                                :: Nxyz
        character(len=LINE_LENGTH)                              :: name
        logical                                                 :: notAddIter

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!       Only update for moving average
!       ------------------------
        if (.not. self % useMovingAverage) return

        if (present(last)) then
            notAddIter = last
        else
            notAddIter = .false.
        end if 

        if (notAddIter) then
            write(name,'(2A)')  trim(self % solutionFileName),'.hsol'
        else
            write(name,'(2A,I10.10,A)')  trim(self % solutionFileName),'_',iter,'.hsol'
        end if 
!
!       Gather reference quantities
!       ---------------------------
        refs(GAMMA_REF) = thermodynamics % gamma
        refs(RGAS_REF)  = thermodynamics % R
        refs(RHO_REF)   = refValues      % rho
        refs(V_REF)     = refValues      % V
        refs(T_REF)     = refValues      % T
        refs(MACH_REF)  = dimensionless  % Mach
        refs(RE_REF)    = dimensionless  % Re

!       Create new file
!       ---------------
        call CreateNewSolutionFile(trim(name),SOLUTION_FILE, mesh % nodeType, mesh % no_of_allElements, iter, time, refs)
        
        padding = NCONS
!
!       Write arrays
!       ------------
        fID = putSolutionFileInWriteDataMode(trim(name))
        ! do eID = 1, self % no_of_elements
        Nxyz = mesh % elements(1) % Nxyz
        allocate(Q(NCONS, 0:Nxyz(1), 0:Nxyz(2), 0:Nxyz(3)))
        do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + padding*e % offsetIO*SIZEOF_RP
                Q(1:NCONS,:,:,:)  = self % Qbase(:,:,:,:,eID)
                call writeArray(fid, Q, position=pos)
            end associate
        end do
        deallocate(Q)

        close(fid)
!
!       Close the file
!       --------------
        call SealSolutionFile(trim(name))

    End Subroutine  writeBaseFlow
!
#endif
End Module SpongeClass
