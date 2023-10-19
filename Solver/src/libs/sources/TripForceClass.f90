!//////////////////////////////////////////////////////
!
!This class represents the numerical force use to induce turbulence, or tripping
!following JFM by Schlatter & Orlu 2012: 
!Turbulent boundary layers at moderate Reynolds numbers: infow length and tripping effects

#include "Includes.h"
Module TripForceClass
  use SMConstants
  use HexMeshClass
  use NodeClass
  use Utilities
  use FTValueDictionaryClass
  use PhysicsStorage
#ifdef HAS_MKL
  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R
#endif
  Implicit None

  private
  public randomTrip

  !definition of gtripClass 
  type gtripClass
    real(kind=RP), dimension(:), allocatable                 :: g, z
    real(kind=RP), dimension(:,:), allocatable               :: h ! the first two dimensions are the old and new values of the random signal respectively, the 3 one is the steady state
    real(kind=RP)                                            :: At, As, ts, zs
    integer                                                  :: Ncutz, N, tsubstep ! N is the number of points to be evaluated in the spanwise direction
    logical                                                  :: isSimetric

  contains

    procedure :: construct      => gtripConstruct
    procedure :: destruct       => gtripDestruct
    procedure :: updateInTime   => gTripForce
    procedure :: randSignal
    procedure :: forceAtz       => getgTripForceAtZ

  end type gtripClass

  !definition of the complete tripClass 
  type ftripClass
      integer                                               :: numberOfTrips
      real(kind=RP), dimension(:,:), allocatable            :: centerPositions
      real(kind=RP), dimension(2)                           :: spatialGaussianAtenuation
      real(kind=RP), dimension(:,:), allocatable            :: normals
      type(gtripClass)                                      :: gTrip
      logical                                               :: active = .false.

      contains

    procedure :: construct      => ftripConstruct
    procedure :: destruct       => ftripDestruct
    procedure :: getTripForce
    procedure :: getTripSource

  end type ftripClass

  type(ftripClass)                                          :: randomTrip

  contains
!/////////////////////////////////////////////////////////////////////////
!           FTRIP PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

  Subroutine ftripConstruct(self, mesh, controlVariables)
    use FileReadingUtilities, only: getRealArrayFromString, getCharArrayFromString
    use Utilities,            only: toLower
    use MPI_Process_Info
    use Headers
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

    class(ftripClass)                                         :: self
    type(HexMesh), intent(in)                                 :: mesh
    type(FTValueDictionary), intent(in)                       :: controlVariables

    !local variables
    real(kind=RP), parameter                                  :: min_normal = 1.0e-15_RP
    integer                                                   :: N, fID, i, partitionRank, ierr
    real(kind=RP)                                             :: x1, x2, y1, y2
    character(len=LINE_LENGTH)                                :: gaussAten_str, trip_zones_str
    character(len=LINE_LENGTH), dimension(:), allocatable     :: trip_zones
    integer, dimension(2)                                     :: faceIndex
    integer, dimension(:), allocatable                        :: trip_zones_index
    real(kind=RP), dimension(NDIM)                            :: temp_normals
    real(kind=RP), dimension(2)                               :: temp_center

    if (controlVariables % containsKey("trip center")) then
        x1 = controlVariables % doublePrecisionValueForKey("trip center")
    else 
        error stop "Trip center must be defined"
    end if

    if (controlVariables % containsKey("trip attenuation")) then
        gaussAten_str = trim(controlVariables % stringValueForKey("trip attenuation", LINE_LENGTH))
        self % spatialGaussianAtenuation = getRealArrayFromString(gaussAten_str)
    else 
        error stop "Trip attenuation must be defined"
    end if
    
    if (controlVariables % containsKey("trip zone")) then
        trip_zones_str = trim(controlVariables % stringValueForKey("trip zone", LINE_LENGTH))
        call toLower(trip_zones_str)
        call getCharArrayFromString(trip_zones_str, LINE_LENGTH, trip_zones)
    else 
        error stop "Trip zone must be defined"
    end if
    
    N = 1
    if (controlVariables % containsKey("trip center 2")) then
        x2 = controlVariables % doublePrecisionValueForKey("trip center 2")
        N = 2
    end if
    
    if (size(trip_zones) .ne. N) error stop "The length of the trip zones is not the same as the trip center length"

    allocate(self % centerPositions(N,2), self % normals(N,NDIM), trip_zones_index(N))
    ! only x of centerPositions is read, y is calculated
    ! if two x are set, the first is search at y>=0 and the second at y<0
    self % numberOfTrips = N

    trip_zones_index = 0
    do i = 1, size(mesh % zones)
      if (trim(mesh % zones(i) % Name) .eq. trim(trip_zones(1))) trip_zones_index(1) = i
      if (N .eq. 2 .and. trim(mesh % zones(i) % Name) .eq. trim(trip_zones(2))) trip_zones_index(2) = i
    end do

    if (trip_zones_index(1) .eq. 0) error stop "Trip zone not found"

    call getFaceTripIndex(mesh % zones(trip_zones_index(1)), mesh, x1, 1, fID, faceIndex, partitionRank, yNegative = .false.)
    if (partitionRank .eq. MPI_Process % rank) then
      temp_center = mesh % faces(fID) % geom % x(1:2,faceIndex(1),faceIndex(2))
      temp_normals = mesh % faces(fID) % geom % normal(:,faceIndex(1),faceIndex(2))
      do i = 1, NDIM
        if (abs(temp_normals(i)) .le. min_normal) temp_normals(i) = 0.0_RP
      end do
    end if
    if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
          call mpi_Bcast(temp_center, 2, MPI_DOUBLE, partitionRank, MPI_COMM_WORLD, ierr)
          call mpi_Bcast(temp_normals, NDIM, MPI_DOUBLE, partitionRank, MPI_COMM_WORLD, ierr)
#endif
    end if 
    self % centerPositions(1,:) = temp_center
    self % normals(1,:) = temp_normals

    if (N .eq. 2) then
      call getFaceTripIndex(mesh % zones(trip_zones_index(2)), mesh, x2, 1, fID, faceIndex, partitionRank, yNegative = .true.)
      if (partitionRank .eq. MPI_Process % rank) then
        temp_center = mesh % faces(fID) % geom % x(1:2,faceIndex(1),faceIndex(2))
        temp_normals = mesh % faces(fID) % geom % normal(:,faceIndex(1),faceIndex(2))
        do i = 1, NDIM
          if (abs(temp_normals(i)) .le. min_normal) temp_normals(i) = 0
        end do
      end if
      if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
        call mpi_Bcast(temp_center, 2, MPI_DOUBLE, partitionRank, MPI_COMM_WORLD, ierr)
        call mpi_Bcast(temp_normals, NDIM, MPI_DOUBLE, partitionRank, MPI_COMM_WORLD, ierr)
#endif
      end if 
      self % centerPositions(2,:) = temp_center
      self % normals(2,:) = temp_normals
    end if

    self % active = .true.

!   Describe the Trip
!   --------------------------
    if ( MPI_Process % isRoot ) then
        write(STD_OUT,'(/)')
        call Subsection_Header("Trip Source Term")
        ! write(STD_OUT,'(30X,A,A28,A)') "->", "Trip Type: ", "Random"
        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of trips: ", N
        if (N .eq. 2) write(STD_OUT,'(30X,A,A28,A,ES10.3,A,ES10.3,A)') "->", "Second Trip center: ", "[",self % centerPositions(2,1),",",self % centerPositions(2,2),"]"
        write(STD_OUT,'(30X,A,A28,A,ES10.3,A,ES10.3,A)') "->", "First Trip center: ", "[",self % centerPositions(1,1),",",self % centerPositions(1,2),"]"
    end if

    call self % gTrip % construct(mesh, controlVariables)

  End Subroutine ftripConstruct
!
  Subroutine ftripDestruct(self)

    class(ftripClass), intent(inout)                          :: self

    safedeallocate(self % centerPositions)
    safedeallocate(self % normals)
    call self % gTrip % destruct()

  End Subroutine ftripDestruct
!
  Subroutine getTripForce(self, x, F2, tripIndex)
      use MPI_Process_Info
    implicit none

    class(ftripClass)                                         :: self
    real(kind=RP), dimension(NDIM), intent(in)                :: x
    real(kind=RP), intent(out)                                :: F2
    integer, intent(out)                                      :: tripIndex

    !local variables
    real(kind=RP)                                             :: g_z
    real(kind=RP), dimension(2,2)                             :: lims, gaussLimits
    integer                                                   :: i

    F2 = 0.0_RP
    tripIndex = 0
    gaussLimits(:,1) = 0.0_RP - self % spatialGaussianAtenuation(:) * 9.0_RP/4.0_RP ! aprox 3 std deviation of gaussian function
    gaussLimits(:,2) = 0.0_RP + self % spatialGaussianAtenuation(:) * 9.0_RP/4.0_RP
    do i = 1, self % numberOfTrips
        lims(1,:) = gaussLimits(1,:) + self % centerPositions(i,1)
        lims(2,:) = gaussLimits(2,:) + self % centerPositions(i,2)
        if( (x(1)>lims(1,1) .and. x(1)<lims(1,2)) .and. (x(2)>lims(2,1) .and. x(2)<lims(2,2)) ) then
            g_z = self % gTrip % forceAtz(x(3))
            F2 = exp( -((x(1) - self % centerPositions(i,1)) / self % spatialGaussianAtenuation(1)) ** 2.0_RP &
                      -((x(2) - self % centerPositions(i,2)) / self % spatialGaussianAtenuation(2)) ** 2.0_RP ) * g_z
            tripIndex = i
            exit
        end if
    end do

  End Subroutine getTripForce
!
  Subroutine getTripSource(self, x, S)
    implicit none

    class(ftripClass)                                         :: self
    real(kind=RP), dimension(NDIM), intent(in)                :: x
    real(kind=RP), dimension(NCONS), intent(inout)            :: S

    !local variables
    real(kind=RP)                                             :: F
    integer                                                   :: tripIndex

    if (.not. self % active) return

    call self % getTripForce(x, F, tripIndex)
    S(IRHOU:IRHOW) = S(IRHOU:IRHOW) - F * self % normals(tripIndex,:)

  End Subroutine getTripSource
!
!/////////////////////////////////////////////////////////////////////////
!           GTRIP PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////
!
  Subroutine  gtripConstruct(self, mesh, controlVariables)

    use Headers
    use MPI_Process_Info
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

    class(gtripClass)                                         :: self
    type(HexMesh), intent(in)                                 :: mesh
    type(FTValueDictionary), intent(in)                       :: controlVariables

    !local variables
    integer, parameter                                        :: seedSizeOriginal = 2
    integer                                                   :: seedSize, seedSizeNeed
    integer, dimension(:), allocatable                        :: seed, seedTemp
    real(kind=RP)                                             :: dz, zMax, zMin, allzMax, allzMin
    integer                                                   :: eID, k, N, Nz, i, j
    integer                                                   :: ierr

#ifdef HAS_MKL
!
!   ----------------------
!   Read Control variables
!   ----------------------
!
    if (controlVariables % containsKey("trip time scale")) then
        self % ts = controlVariables % doublePrecisionValueForKey("trip time scale")
    else
        error stop "trip time scale must be defined"
    end if
    if (controlVariables % containsKey("trip number of modes")) then
        self % Ncutz = controlVariables % IntegerValueForKey("trip number of modes")
    else
        error stop "trip number of modes must be defined"
    end if
    if (controlVariables % containsKey("trip z points")) then
        N = controlVariables % IntegerValueForKey("trip z points")
    else
        error stop "trip z points must be defined"
    end if

    self % N = N
    if (self % Ncutz .ge. self % N) error stop 'The trip cuttoff number of modes is bigger than the number of points in z'

    self % At = controlVariables % getValueOrDefault("trip amplitude",1.0)
    self % As = controlVariables % getValueOrDefault("trip amplitude steady",0.0)
    self % isSimetric = .FALSE.

    !this needs to be updated for different simulations, but remain constant for the same one if is restarted
    ! the numbers are obtained from a call to RANDOM_SEED with a get; and are hardcoded to the control file or blank and use default
    allocate(seedTemp(seedSizeOriginal))
    seedTemp(1) = controlVariables % getValueOrDefault("random seed 1",930187532)
    seedTemp(2) = controlVariables % getValueOrDefault("random seed 2",597734650)
    if ( MPI_Process % isRoot ) then
        call random_seed(size=seedSizeNeed)
        if (seedSizeNeed .gt. seedSizeOriginal) then
            seedSize = seedSizeNeed
            allocate(seed(seedSize))
            seed(1:seedSizeOriginal) = seedTemp(:)
            do i = seedSizeOriginal+1, seedSizeNeed
                j = mod(i+1, 2) + 1
                seed(i) = seedTemp(j)
            end do
        else
            seedSize = seedSizeOriginal
            allocate(seed(seedSize))
            seed = seedTemp
        end if
        deallocate(seedTemp)
        call random_seed (PUT=seed(1:seedSize))
    end if

    !search for the maximum and minimum values of z in the whole mesh
!!$omp parallel do private(eID,k) reduction(MIN:zMin) reduction(MAX:zMax) schedule(runtime)
    zMin = mesh % nodes(mesh % elements(1) % nodeIDs(1)) % x(3)
    zMax = mesh % nodes(mesh % elements(1) % nodeIDs(1)) % x(3)
    do eID = 1, mesh % no_of_elements
       associate ( e => mesh % elements(eID) )
        do k = 1, 6
           associate ( xx => mesh % nodes(e % nodeIDs(k))  %  x(3))
             if (xx < zMin) zMin = xx
             if (xx > zMax) zMax = xx
         end associate
         end do
       end associate
    end do
!!$omp end parallel do
    if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
      call mpi_allreduce(zMin, allzMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
      call mpi_allreduce(zMax, allzMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    else
      allzMin = zMin
      allzMax = zMax
    end if

    allocate(self % g(0:N-1), self % h(0:N-1,3), self % z(0:N-1))
    ! todo: check why this is used
    ! self % Ncutz = int((xz(2) - xz(1))/ self % zs)
    self % zs = (allzMax - allzMin)/ real(self % Ncutz,RP)

    ! create evenly spaced z array
    dz = (allzMax - allzMin) / (N-1)
    self % z(0:N-1) = [(allzMin + k*dz, k=0, N-1)] ! implicit loop, mini cycle

    ! start of simulation, if is a restart it will be updated
    self % tsubstep = -1

    !get the steady state random signal
    call randSignal(self, self % h(:,3))
    !get the first random signal, for i=0
    call randSignal(self, self % h(:,2))

!   Describe the Trip
!   --------------------------
    if ( .not. MPI_Process % isRoot ) return
    write(STD_OUT,'(30X,A,A28,A,I0,A,I0,A)') "->", "Random seeds: ", "[",seed(1),",",seed(2),"]"
    write(STD_OUT,'(30X,A,A28,F7.4)') "->", "Trip Amplitude: ", self % At
    write(STD_OUT,'(30X,A,A28,ES10.3)') "->", "Trip Time scale: ", self % ts
    write(STD_OUT,'(30X,A,A28,ES10.3)') "->", "Trip Space scale: ", self % zs

#else
      error stop 'MKL not linked correctly'
#endif

  End Subroutine gtripConstruct 
!
  Subroutine  gtripDestruct(self)

    class(gtripClass), intent(inout)                          :: self
 
    if (ALLOCATED(self%g)) DEALLOCATE (self%g)
    if (ALLOCATED(self%h)) DEALLOCATE (self%h)
    if (ALLOCATED(self%z)) DEALLOCATE (self%z)
    self % tsubstep = -1

  End Subroutine gtripDestruct 
!
  Subroutine randSignal(self, signal)
  ! creates a random signal (array) with the first Ncutz modes random and .le. 1.0 and 0 for modes greater than Ncutz.
  ! The array is created in an evenly distributed z (spanwise direction), it needs to be interpolated at the mesh points later.

    use MPI_Process_Info
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

    class(gtripClass)                                         :: self
    real(kind=RP), dimension(0:self%N-1), intent(out)         :: signal

    !local variables
    integer                                                   :: i, ierr
    integer                                                   :: nh ! half of cuttoff number
    complex(kind=CP)                                          :: fourierArg
    real(kind=RP), dimension(self%Ncutz/2+1)                  :: ranNum
    real(kind=RP), dimension(0:self%N+1)                      :: fourierCoefficients, preSignal

    if ( MPI_Process % isRoot ) then
        nh = self%Ncutz/2
        ! the signal coefficients are initialized at zero, it'll be changed for frequencies <= Ncutz
        fourierCoefficients = 0_RP

        ! the mean value, n=0 is random and ampliated by sqrt 2 for the real part and 0 for the img
        ! the coef array is packed in CCS format
        call random_number(ranNum)
        fourierArg = exp(ImgI*2.0_RP*PI*ranNum(1))
        fourierCoefficients(0) = real(fourierArg)*sqrt(2.0_RP)
        fourierCoefficients(1) = 0

        ! for the other values of n less than the cut off, the random amplitude is ampliated by sqrt 2 for the real part and 0 for the
        ! img for the symmetric case
        ! for the non symmetric case, the amplitude is set to 1.0 and the phase is random (both real and part are saved)
        if (self%isSimetric) then
          do i = 1, nh
            fourierArg = exp(ImgI*2.0_RP*PI*ranNum(i+1))
            fourierCoefficients(2*i) = real(fourierArg)*sqrt(2.0_RP)
            fourierCoefficients(2*i+1) = 0
          end do  
        else
          do i = 1, nh
            fourierArg = exp(ImgI*2.0_RP*PI*ranNum(i+1))
            fourierCoefficients(2*i) = real(fourierArg)
            fourierCoefficients(2*i+1) = aimag(fourierArg)
          end do
        end if

        call backfft(preSignal, self % N, fourierCoefficients)
        ! scale and give only the part of the array of interest
    end if
    if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
          call mpi_Bcast(preSignal, self % N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
#endif
    end if
    signal = 1.0/real(self % Ncutz)*preSignal(0:self % N-1)

  End Subroutine randSignal
!
  Subroutine  gTripForce(self, t)
    !compute the g term of the trip force equation

    implicit none

    class(gtripClass)                                         :: self
    real(kind=RP), intent(in)                                 :: t

    !local variables
    integer                                                   :: i !integer quotient for local time and time cuttoff scale (ts)
    integer                                                   :: j !time iterator
    real(kind=RP)                                             :: b, p ! terms of gTripForce equation, related to time interpolant

    i = int(t/self % ts)
    p = (t - real(i)*self % ts)/self % ts
    b = p*p*(3.0_RP - 2.0_RP*p)

    do j = self % tsubstep + 1, i
      self % h(:,1) = self % h(:,2)
      call randSignal(self,self%h(:,2) )
    end do

    self % tsubstep = i
    self % g = self % At * ( (1.0_RP-b) * self % h(:,1) + b * self % h(:,2) ) + self % As * self % h(:,3)

  End Subroutine gTripForce
!
  Subroutine forfft(coef,x,N,coeffCCS)

    integer, intent(in)                                       :: N
    real(kind=RP), dimension(0:N-1), intent(in)               :: x
    real(kind=RP), dimension(0:N+1), intent(out)               :: coeffCCS
    complex(kind=CP), dimension(0:int(N/2)), intent(out)           :: coef

    ! local variables
    real(kind=RP), dimension(0:N+1)                           :: y
#ifdef HAS_MKL
    type(dfti_descriptor), pointer                            :: hand
#endif
    integer                                                   :: status, i, nh

#ifdef HAS_MKL

    y(0:N-1) = x

    status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 1, N)
    status = DftiCommitDescriptor(hand)
    status = DftiComputeForward( hand, y )
    status = DftiFreeDescriptor(hand)

    nh = N/2
    do i= 0,nh
      coef(i) = cmplx(y(2*i), y(2*i+1))
    end do
    ! if (mod(N,2).ne.0) nh = nh +1
    ! do i= 1,nh-1
    !   coef(nh+i) = conjg(coef(nh-i))
    ! end do
    coeffCCS = y
#else
      error stop 'MKL not linked correctly'
#endif

  End Subroutine forfft 
!
  Subroutine backfft(x,N,coeffCCS)

    integer, intent(in)                                       :: N
    real(kind=RP), dimension(0:N+1), intent(in)               :: coeffCCS
    real(kind=RP), dimension(0:N-1), intent(out)              :: x

    ! local variables
    real(kind=RP), dimension(0:N+1)                           :: y
#ifdef HAS_MKL
    type(dfti_descriptor), pointer                            :: hand
#endif
    integer                                                   :: status, i

#ifdef HAS_MKL

    y = coeffCCS

    status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 1, N)
    status = DftiCommitDescriptor(hand)
    status = DftiComputeForward( hand, y )
    status = DftiFreeDescriptor(hand)

    x = y(0:N-1)

#else
      error stop 'MKL not linked correctly'
#endif

  End Subroutine backfft 
!
  Function getgTripForceAtZ(self, specificZ) result(gz)
    ! returns the value of the g term at a given z position, using linear interpolation
    ! it assumes that the array g corresponds position by position to the array z

    implicit none

    class(gtripClass)                                         :: self
    real(kind=RP), intent(in)                                 :: specificZ
    real(kind=RP)                                             :: gz   !g term at a z position

    !local variables
    integer                                                   :: i

    !starts in 1 in case there is a value .eq. to the 0 position
    do i = 1, self % N - 1
      if (specificZ .ge. self % z(i)) exit
    end do  

    gz = self % g(i-1) + ( self % g(i) - self % g(i-1) ) / ( self % z(i) - self % z(i-1) ) * (specificZ - self%z(i-1))

  End Function getgTripForceAtZ 
!
!/////////////////////////////////////////////////////////////////////////
!           HELPER PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////
!
  Subroutine getFaceTripIndex(zone, mesh, x, x_dim, faceID, faceIndex, partitionRank, yNegative)
    use ZoneClass
    use MPI_Process_Info
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

    type(Zone_t), intent(in)                                  :: zone
    type(HexMesh), intent(in)                                 :: mesh
    real(kind=RP), intent(in)                                 :: x
    integer, intent(in)                                       :: x_dim
    logical, intent(in), optional                             :: yNegative
    integer, intent(out)                                      :: faceID
    integer, dimension(2), intent(out)                        :: faceIndex
    integer, intent(out)                                      :: partitionRank

    !local variables
    integer                                                   :: i, j, faceZoneID, fID, yIndex, allfID
    integer                                                   :: nx, ny, ierr
    real(kind=RP)                                             :: x_diff, minx_diff, lmin, lmax
    logical                                                   :: onlyYNegative
    integer                                                   :: tempRank

    if (present(yNegative)) then
        onlyYNegative = yNegative
    else
        onlyYNegative = .false.
    end if 

    yIndex = x_dim + 1
    if (yIndex .gt. NDIM) yIndex = yIndex - NDIM

    faceID = 0
    do faceZoneID = 1, zone % no_of_faces
      fID = zone % faces(faceZoneID)
      associate( f => mesh % faces(fID) )
        nx = f % Nf(1)
        ny = f % Nf(2)
        lmin = minval([f % geom % x(x_dim,0,0), f % geom % x(x_dim,nx,ny)])
        lmax = maxval([f % geom % x(x_dim,0,0), f % geom % x(x_dim,nx,ny)])
        if ( (x .ge. lmin) .and. (x .le. lmax) ) then
          if (onlyYNegative .eqv. (f % geom % x(yIndex,0,0) .lt. 0.0_RP)) then
            faceID = fID
            exit
          end if
        end if 
      end associate
    end do

    ! if partition does not found the face, just store the info of the one that did it, and which partition number is
    if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_allreduce(faceID, allfID, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

      if (allfID .eq. 0) error stop "trip center not found in any face of the zone"

      ! start value for allreduce
      tempRank = 0
      if (allfID .eq. faceID) tempRank = MPI_Process % rank
#ifdef _HAS_MPI_
      call mpi_allreduce(tempRank, partitionRank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    else
      if (faceID .eq. 0) error stop "trip center not found in any face of the zone"
      partitionRank = 0
    end if

    !only calculate faceIndex if its in the partition
    faceIndex = 0
    if (partitionRank .eq. MPI_Process % rank) then
        minx_diff = abs(mesh % faces(faceID) % geom % x(x_dim,0,0) - x)
        do j = 0, ny
          do i = 0, nx
            x_diff = abs(mesh % faces(faceID) % geom % x(x_dim,i,j) - x)
            if (x_diff .le. minx_diff) then
                minx_diff = x_diff
                faceIndex = [i,j]
            end if
          end do
        end do
    end if

  End Subroutine getFaceTripIndex
!
End Module TripForceClass