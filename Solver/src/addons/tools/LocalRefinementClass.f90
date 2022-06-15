!This class represents the local p-refinement of a mesh

Module  LocalRefinement  !

    use FileReadingUtilities, only: getRealArrayFromString
    use FTValueDictionaryClass
    use SMConstants
    use pAdaptationClass          , only: GetMeshPolynomialOrders

    Implicit None
!     ---------------------
!     Main local refinement type
!     ---------------------

    type LocalRef_t

        real(kind=RP), dimension(:,:), allocatable          :: xlim   ! limits for refinement in x direction
        real(kind=RP), dimension(:,:), allocatable          :: ylim
        real(kind=RP), dimension(:,:), allocatable          :: zlim
        integer, dimension(:), allocatable         :: Nx           ! polynomial order for each selection whithin the limits in x diredirection
        integer, dimension(:), allocatable         :: Ny
        integer, dimension(:), allocatable         :: Nz
        integer, dimension(3)                      :: globalNxyz   ! default polynomial order, outside the limits
        integer                                    :: lenRegions   ! number of regions, it has to be the same for each diredirection

        contains
           procedure   :: Construct               => LocalRef_Construct
           procedure   :: Destruct                => LocalRef_Destruct
           procedure   :: getOrderOfPosition

    end type LocalRef_t
!
!////////////////////////////////////////////////////////////////////////
    contains
!////////////////////////////////////////////////////////////////////////
!
    Subroutine  LocalRef_Construct(self, controlVariables)

        IMPLICIT NONE

        class(LocalRef_t)                         :: self
        type(FTValueDictionary), intent(in)     :: controlVariables

        !local variables
        character(len=LINE_LENGTH)              :: NxName, NyName, NzName, limsName
        integer                                 :: lenRegions, lenRegionz, lenRegiony
        real(kind=RP),dimension(:),allocatable  :: tempArray
        integer, allocatable                    :: Nx(:), Ny(:), Nz(:) !temporal variables for getting the default pol orders
        integer                                 :: Nmax

        !get polynomial orders for all regions in all diredirections
        NxName = trim(controlVariables%stringValueForKey("x regions orders",LINE_LENGTH))
        NyName = trim(controlVariables%stringValueForKey("y regions orders",LINE_LENGTH))
        NzName = trim(controlVariables%stringValueForKey("z regions orders",LINE_LENGTH))

        !get the size of each array, that is remove count of [] and all "," (there are one less than the size)
        lenRegions = (len_trim(NxName) - 1) / 2
        lenRegiony = (len_trim(NyName) - 1) / 2
        lenRegionz = (len_trim(NzName) - 1) / 2

        Allocate (self % xlim(2,lenRegions), self % ylim(2,lenRegions), self % zlim(2,lenRegions) )
        Allocate ( self % Nx(lenRegions), self % Ny(lenRegions), self % Nz(lenRegions))

        self % lenRegions              = lenRegions

        !get global pol orders
        call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
        self % globalNxyz(1) = Nx(1)
        self % globalNxyz(2) = Ny(1)
        self % globalNxyz(3) = Nz(1)

        self % Nx = getRealArrayFromString(NxName)
        self % Ny = getRealArrayFromString(NyName)

        !get the limits of the regions for all directions
        allocate (tempArray(2*lenRegions))

        limsName = trim(controlVariables%stringValueForKey("x regions limits",LINE_LENGTH))
        tempArray = getRealArrayFromString(limsName)
        ! convert to a two dimensional array
        self % xlim = reshape(source= tempArray, shape=[2,lenRegions])

        limsName = trim(controlVariables%stringValueForKey("y regions limits",LINE_LENGTH))
        tempArray = getRealArrayFromString(limsName)
        self % ylim = reshape(source= tempArray, shape=[2,lenRegions])

        !check if z arrays are set in control file
        if (lenRegionz .gt. 0) then 
            self % Nz = getRealArrayFromString(NzName)
            limsName = trim(controlVariables%stringValueForKey("z regions limits",LINE_LENGTH))
            tempArray = getRealArrayFromString(limsName)
            self % zlim = reshape(source= tempArray, shape=[2,lenRegions])
        else
            ! set default values
            self % Nz = Nz(1)
            self % zlim(1,:) = -10
            self % zlim(2,:) = 10
        end if

        deallocate(tempArray)

!
    End Subroutine LocalRef_Construct 
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine  LocalRef_Destruct(self)

        class(LocalRef_t), intent(inout) :: self

        if (ALLOCATED(self%xlim)) DEALLOCATE (self%xlim)
        if (ALLOCATED(self%ylim)) DEALLOCATE (self%ylim)
        if (ALLOCATED(self%zlim)) DEALLOCATE (self%zlim)
        if (ALLOCATED(self%Nx)) DEALLOCATE (self%Nx)
        if (ALLOCATED(self%Ny)) DEALLOCATE (self%Ny)
        if (ALLOCATED(self%Nz)) DEALLOCATE (self%Nz)
        self % globalNxyz = 1
        self % lenRegions = 0

    End Subroutine LocalRef_Destruct 
!
!////////////////////////////////////////////////////////////////////////
!
    Subroutine  getOrderOfPosition(self, coordinates, Nx, Ny, Nz)

        implicit none

        class (LocalRef_t)                       :: self
        real(kind=RP), dimension(3,8), intent(in)       :: coordinates
        integer,intent(out)                    :: Nx, Ny, Nz

        !local variables
        integer                                :: i

        ! set default value if the coordinates are not in the limits
        Nx = self % globalNxyz(1)
        Ny = self % globalNxyz(2)
        Nz = self % globalNxyz(3)

        ! search for the coordinates whithin the limits and set the order of the corresponded one
         associate (cx => coordinates(1,:)); associate (cy => coordinates(2,:)); associate (cz => coordinates(3,:))
            do i = 1, self%lenRegions
                if (all(cx .ge. self%xlim(1,i)) .and. all(cx .le. self%xlim(2,i))) then
                    if (all(cy .ge. self%ylim(1,i)) .and. all(cy .le. self%ylim(2,i))) then
                        if (all(cz .ge. self%zlim(1,i)) .and. all(cz .le. self%zlim(2,i))) then
                            Nx = self % Nx(i)
                            Ny = self % Ny(i)
                            Nz = self % Nz(i)
                            return
                        end if
                    end if
                end if  
            end do  
        end associate; end associate; end associate

    End Subroutine getOrderOfPosition 

End Module  LocalRefinement 
