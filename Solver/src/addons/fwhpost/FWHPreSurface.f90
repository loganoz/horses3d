#include "Includes.h"
Module FWHPreSurface  !
    use SMConstants
    use Headers
    use HexMeshClass
    use FTValueDictionaryClass
    Implicit None

    contains

    Subroutine getElements(mesh, controlVariables)

        type(HexMesh)                               :: mesh
        type(FTValueDictionary)                     :: controlVariables

        !local variables
        real(kind=RP)                               :: r
        real(kind=RP), dimension(:), allocatable    :: x, y, z
        integer                                     :: eID, i, j, k

        r = 0.4_RP

                        print *, "start"
        do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
                ! do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                !     x = e % geom % x(:,i,j,k)
                ! end do                  ; end do                ; end do
                allocate( x(e % Nxyz(1)) )
                allocate( y(e % Nxyz(2)) )
                allocate( z(e % Nxyz(3)) )
                x = e % geom % x(1,:,0,0)
                y = e % geom % x(2,0,:,0)
                z = e % geom % x(3,0,0,:)
                
                ! print *, "z: ", z, "z1", z(1)
                if (z(1) .gt. 0.8_RP) then
                        ! print *, "gz"
                    if (isInCircle(r, x, y, e%Nxyz(1), e%Nxyz(2))) then
                        print *, "eID: ", e%eID
                        if (e%eID .eq. 304) then
                            print *, "x: ", x, "y: ", y
                            print *, "r2", POW2(x(i) - 0.5_RP) + POW2(y(j) - 0.5_RP)
                        end if
                    end if
                end if 
                deallocate(x,y,z)
            end associate
        end do

    End Subroutine getElements

    Function isInCircle(r, x, y, Nx, Ny) result(isIn)

        real(kind=RP)                               :: r
        real(kind=RP), dimension(Nx)                :: x
        real(kind=RP), dimension(Ny)                :: y
        integer                                     :: Nx, Ny
        logical                                     :: isIn

        ! local variables
        integer                                     :: i, j

        ! isIn = .true.
        ! isIn = .false.
        do j = 0, Ny ; do i = 0, Nx
            ! isIn = isIn .and. ( POW2(r) .gt. (POW2(x(i) - 0.5_RP) + POW2(y(j) - 0.5_RP)) )
            if ( POW2(r) .lt. (POW2(x(i) - 0.5_RP) + POW2(y(j) - 0.5_RP)) ) then
                isIn = .false.
                return
            end if
        end do       ; end do
        isIn = .true.

    End Function isInCircle

End Module FWHPreSurface
