#include "Includes.h"
Module LambVectorInterpolation  !
#if defined(ACOUSTIC) 
    use SMConstants
    Implicit None

    private
    public LambInterpolation

    ! Definition of class 
    type :: LambVectorInterpolation_t
        logical :: isActive
        integer :: interpolationType
        character(len=100), allocatable :: filenames(:)
        real(rp), allocatable :: times(:)
        integer :: currTimeIndex
        contains
        procedure :: construct      => LVI_construct
        procedure :: interpolate    => LVI_interpolate
    end type LambVectorInterpolation_t


    enum, bind(C)
        enumerator :: TYPE_CONSTANT = 1, TYPE_LINEAR
    end enum

    type(LambVectorInterpolation_t) :: LambInterpolation

    contains

    subroutine LVI_construct(self, controlVariables)
        use Physics_CAAKeywordsModule
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit none
        class(LambVectorInterpolation_t) :: self
        type(FTValueDictionary), intent(in) :: controlVariables

        ! Local variables
        character(len=KEYWORD_LENGTH) :: keyword

        self % isActive = controlVariables % logicalValueForKey(LAMB_VECTOR_KEY)

        if (controlVariables % containsKey(LambInterpolationTypeKey)) then
            keyword = controlVariables % stringValueForKey(LambInterpolationTypeKey, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case (INTERPOLATION_CONSTANT_NAME)
               self % interpolationType = TYPE_CONSTANT
            case (INTERPOLATION_LINEAR_NAME)
                self % interpolationType = TYPE_LINEAR
            end select
        else
            self % interpolationType = TYPE_CONSTANT
        end if

        call getFilenamesAndTimes(self, controlVariables)

        self % currTimeIndex = -1

    end subroutine LVI_construct

    subroutine LVI_interpolate(self, mesh, time)
        use HexMeshClass
        implicit none
        class(LambVectorInterpolation_t) :: self
        class(HexMesh) :: mesh
        real(rp), intent(in) :: time

        if (.not. self % isActive) return

        select case (self % interpolationType)
        case (TYPE_CONSTANT)
            call LVI_interpolateConstant(self, mesh, time)
        case (TYPE_LINEAR)
            call LVI_interpolateLinear(self, mesh, time)
        end select

    end subroutine LVI_interpolate

    subroutine LVI_interpolateConstant(self, mesh, time)
        use HexMeshClass
        implicit none
        class(LambVectorInterpolation_t) :: self
        class(HexMesh) :: mesh
        real(rp), intent(in) :: time
        
        ! Local variables
        character(len=100) :: filename
        integer :: eID

        ! Read data from the only file that exists
        self % currTimeIndex = 1
        filename = self % filenames(self % currTimeIndex)
        call mesh % LoadLambVector(filename, indLamb_=1)
        do eID = 1, mesh % no_of_elements
            mesh % elements(eID) % storage % Lamb(:,:,:,:) = mesh % elements(eID) % storage % Lambread(:,:,:,:,1)
        end do

        ! Deactivate to avoid file reading in next iterations
        self % isActive = .false.

    end subroutine LVI_interpolateConstant

    subroutine LVI_interpolateLinear(self, mesh, time)
        use HexMeshClass
        implicit none
        class(LambVectorInterpolation_t) :: self
        class(HexMesh) :: mesh
        real(rp), intent(in) :: time
        
        ! Local variables
        character(len=100) :: filename
        integer :: newTimeIndex
        integer :: eID
        real(rp) :: t0, t1, alpha
        
        newTimeIndex = findPreviousNext(self, time)

        ! Check if the time interval has changed
        if (self % currTimeIndex .ne. newTimeIndex) then
            self % currTimeIndex = newTimeIndex
            ! Read previous time instant data
            filename = self % filenames(self % currTimeIndex)
            call mesh % LoadLambVector(filename, indLamb_=1)

            ! Read next time instant data
            filename = self % filenames(self % currTimeIndex + 1)
            call mesh % LoadLambVector(filename, indLamb_=2)
        end if

        ! Interpolate
        t0 = self % times(self % currTimeIndex)
        t1 = self % times(self % currTimeIndex+1)
        alpha = (time - t0) / (t1 - t0)
        do eID = 1, mesh % no_of_elements
            mesh % elements(eID) % storage % Lamb(:,:,:,:) = mesh % elements(eID) % storage % Lambread(:,:,:,:,1) * (1.0_rp - alpha) + alpha *  mesh % elements(eID) % storage % Lambread(:,:,:,:,2)
        end do

    end subroutine LVI_interpolateLinear

    subroutine getFilenamesAndTimes(self, controlVariables)
        use FileReadingUtilities, only: readFilesByGlob
        use Utilities, only: BubblesortWithFriend
        use SolutionFile, only: getSolutionFileTimeAndIteration
        use FTValueDictionaryClass
        use Physics_CAAKeywordsModule
        implicit none
        type(LambVectorInterpolation_t) :: self
        type(FTValueDictionary), intent(in) :: controlVariables


        integer, allocatable :: perm(:)
        character(len=100) :: dir, basename, filename
        integer :: i, iter

        ! Read the filenames matching pattern
        dir = controlVariables % stringValueForKey(LambReadDirectoryKey, LINE_LENGTH)
        basename = controlVariables % stringValueForKey(LambVectorFileNameKey, LINE_LENGTH)
        call readFilesByGlob(dir, basename, self % filenames)

        ! Read times from file and sort them
        allocate(self % times(size(self % filenames)))
        allocate(perm(size(self % filenames)))
        do i = 1, size(self % filenames)
            call getSolutionFileTimeAndIteration(trim(self % filenames(i)), iter, self % times(i))
            perm(i) = i
        end do
        call BubblesortWithFriend(self % times, perm)
        self % filenames = self % filenames(perm)

    end subroutine

    function findPreviousNext(self, time) result(i)
        implicit none
        class(LambVectorInterpolation_t), intent(in) :: self
        real(rp), intent(in) :: time

        integer :: i

        i = 1
        do
            if ((self % times(i) <= time) .and. (time < self % times(i+1))) then
                return
            end if
            i = i+1
        end do
    end function findPreviousNext

#endif
End Module LambVectorInterpolation
