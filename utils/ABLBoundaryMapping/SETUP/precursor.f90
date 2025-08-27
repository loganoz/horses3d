Module precursor  
    use smconstants
    Implicit None

    private

    public read_wind_precursor_table, filename, num_dirs, file_names, times, directories, &
           read_vtk_mesh, read_vtk_vectorfield, read_vtk_scalarfield, time_interp_weights, &
           i_tl_prev, i_tu_prev, dirname, iformat, &
           points, centroids, polygons, num_points, num_polygons, &
           field_data_U_l, field_data_U_u, field_data_p_u, field_data_p_l, field_data_T_u, field_data_T_l, &
           ind, tnbr, til, ntri, &
           normalize, z_order, &
           hashtable, simplices, outside, conn_mat, Tmat, weights, &
           find_closest_point, compute_affine_transform, compute_barycentric_weights, &
           max_x, min_x, max_y, min_y

    character(len=256) :: filename, dirname, iformat 
    integer :: num_dirs, i_tl_prev, i_tu_prev
    character(len=256), dimension(3) :: file_names
    real(kind=RP), dimension(:), allocatable :: times
    character(len=256), dimension(:), allocatable :: directories
    real(kind=RP), allocatable :: points(:,:), centroids(:,:)
    integer, allocatable :: polygons(:,:)
    integer :: num_points, num_polygons
    real(kind=RP), allocatable :: field_data_U_l(:, :), field_data_U_u(:, :), field_data_p_u(:), field_data_p_l(:), field_data_T_u(:), field_data_T_l(:)
    INTEGER(kind=4), ALLOCATABLE :: ind(:), tnbr(:, :), til(:, :)
    integer, allocatable :: hashtable(:), simplices(:), outside(:)
    integer(kind=4) :: ntri
    integer(kind=4), allocatable :: conn_mat(:,:)
    REAL(kind=RP), ALLOCATABLE :: Tmat(:, :, :), weights(:, :)
    real(kind=RP) :: max_x, min_x, max_y, min_y
    
    contains

    subroutine read_wind_precursor_table()
        implicit none
        character(len=256) :: line, trimmed_line
        integer :: i, j, ios, unit
        real(kind=RP) :: temp_time
        character(len=256) :: temp_dir

        filename = 'PRECURSOR/precursor.dat'
      
        ! Open the input file
        open(unit=10, file=filename, status='old', action='read')
      
        ! Skip comment lines
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) then
                print *, 'Error reading file or end of file reached'
                stop
            end if
            trimmed_line = trim(adjustl(line))
            if (trimmed_line(1:1) /= '#') exit
        end do
      
        ! Read the number of directories
        read(trimmed_line, '(I10)', iostat=ios) num_dirs
        if (ios /= 0) then
            print *, 'Error reading number of directories'
            stop
        end if
      
        ! Skip comment lines
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) then
                print *, 'Error reading file or end of file reached'
                stop
            end if
            trimmed_line = trim(adjustl(line))
            if (trimmed_line(1:1) /= '#') exit
        end do
      
        ! Read the file names
        read(trimmed_line, '(A)', iostat=ios) trimmed_line
        if (ios /= 0) then
            print *, 'Error reading file names'
            stop
        end if
        read(trimmed_line, *, iostat=ios) file_names(1), file_names(2), file_names(3)
        if (ios /= 0) then
            print *, 'Error parsing file names'
            stop
        end if
      
        ! Skip comment lines
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) then
                print *, 'Error reading file or end of file reached'
                stop
            end if
            trimmed_line = trim(adjustl(line))
            if (trimmed_line(1:1) /= '#') exit
        end do

        ! Read file format
        read(trimmed_line, '(A)', iostat=ios) iformat
        if (ios /= 0) then
            print *, 'Error reading file format'
            stop
        end if
        if (iformat /= 'ascii' .and. iformat /= 'binary') then
            print *, 'Error in file format. Only "ascii" or "binary" allowed'
            stop
        end if

        ! Skip comment lines
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) then
                print *, 'Error reading file or end of file reached'
                stop
            end if
            trimmed_line = trim(adjustl(line))
            if (trimmed_line(1:1) /= '#') exit
        end do
      
        ! Allocate arrays for times and directories
        allocate(times(num_dirs))
        allocate(directories(num_dirs))
      
        ! Read the times and directories
        do i = 1, num_dirs
            read(trimmed_line, *, iostat=ios) times(i), directories(i)
            if (ios /= 0) then
                stop
            end if
            if (i < num_dirs) then
                read(10, '(A)', iostat=ios) line
                if (ios /= 0) then
                    stop
                end if
                trimmed_line = trim(adjustl(line))
            end if
        end do
      
        ! Close the input file
        close(10)
      
        ! Sort the times and corresponding directories in ascending order
        do i = 1, num_dirs-1
            do j = i+1, num_dirs
                if (times(i) > times(j)) then
                    ! Swap times
                    temp_time = times(i)
                    times(i) = times(j)
                    times(j) = temp_time
      
                    ! Swap directories
                    temp_dir = directories(i)
                    directories(i) = directories(j)
                    directories(j) = temp_dir
                end if
            end do
        end do

        ! Initialise indexing for precursor files
        i_tl_prev = -1 ! -1 guarantees that in first time iteration source files are read
        i_tu_prev = -1
      
    end subroutine read_wind_precursor_table

    subroutine read_vtk_mesh(filename, num_points, points, num_polygons, polygons)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: num_points, num_polygons
        real(kind=8), allocatable, intent(out) :: points(:,:)
        integer, allocatable, intent(out) :: polygons(:,:)
        integer :: i, ios, unit
        character(len=256) :: line, dummy
      
        if (iformat == 'ascii') then

            ! Open the file
            open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening file: ', filename
                stop
            end if
        
            ! Read the header
            do
                read(unit, '(A)', iostat=ios) line
                if (index(line, 'POINTS') > 0) exit
            end do
        
            ! Read the number of polygons
            read(line, '(A)', iostat=ios) line
            read(line, *) dummy, num_points
        
            ! Allocate the points array
            IF (ALLOCATED(points)) THEN
                DEALLOCATE(points)
            END IF
            allocate(points(3, num_points))
        
            ! Read the points
            do i = 1, num_points
                read(unit, *, iostat=ios) points(:, i)
                if (ios /= 0) then
                    print *, 'Error reading point ', i
                    stop
                end if
            end do
        
            ! Read until the polygons section
            do
            read(unit, '(A)', iostat=ios) line
            if (index(line, 'POLYGONS') > 0) exit
            end do
        
            ! Read the number of polygons
            read(line, '(A)', iostat=ios) line
            read(line, *) dummy, num_polygons
        
            ! Allocate the polygons array
            IF (ALLOCATED(polygons)) THEN
                DEALLOCATE(polygons)
            END IF
            allocate(polygons(5, num_polygons))
        
            ! Read the polygons
            do i = 1, num_polygons
                read(unit, *, iostat=ios) polygons(:, i)
                if (ios /= 0) then
                    print *, 'Error reading polygon ', i
                    stop
                end if
            end do
        
            ! Close the file
            close(unit)
        
        else if (iformat == 'binary') then

            ! Open the binary input file
            open(newunit=unit, file=filename, status='old', form='unformatted', access='stream', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening input file'
                stop
            end if

            ! Read number of points in mesh
            read(unit) num_points

            ! Allocate the points array
            IF (ALLOCATED(points)) THEN
                DEALLOCATE(points)
            END IF
            allocate(points(3, num_points))

            ! Read the points
            read(unit) points

            ! Read number of polygons
            read(unit) num_polygons

            ! Allocate the polygons array
            IF (ALLOCATED(polygons)) THEN
                DEALLOCATE(polygons)
            END IF
            allocate(polygons(5, num_polygons))

            ! Read the polygons
            read(unit) polygons

            ! Close the files
            close(unit)

        else
            print *, 'Unrecognised input file format: ', iformat
            stop
        end if

    end subroutine read_vtk_mesh
      
    subroutine read_vtk_vectorfield(filename, field_data)
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=RP), allocatable, intent(out) :: field_data(:,:)
        integer :: i, ios, unit
        character(len=256) :: line, field_name, field_type
        integer :: num_cells, num_components, dummy_int, dummy_int_aux
        real(kind=RP) :: dummy_real
      
        if (iformat == 'ascii') then

            ! Open the file
            open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening file: ', filename
                stop
            end if
        
            ! Read until the field data section
            do
                read(unit, '(A)', iostat=ios) line
                if (index(line, 'FIELD') > 0) exit
            end do
            read(unit, '(A)', iostat=ios) line
            read(line, *) field_name, num_components, num_cells, field_type
        
            ! Allocate the field data array
            IF (ALLOCATED(field_data)) THEN
                DEALLOCATE(field_data)
            END IF
            allocate(field_data(num_components, num_cells))
        
            ! Read the field data
            do i = 1, num_cells
                read(unit, *, iostat=ios) field_data(:, i)
                if (ios /= 0) then
                    print *, 'Error reading field data at cell ', i
                    stop
                end if
            end do
        
            ! Close the file
            close(unit)

        else if (iformat == 'binary') then

            ! Open the binary input file
            open(newunit=unit, file=filename, status='old', form='unformatted', access='stream', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening input file'
                stop
            end if

            ! Skip points
            read(unit) dummy_int
            do i = 1, dummy_int * 3
                read(unit) dummy_real
            end do

            ! Skip polygons
            read(unit) dummy_int
            do i = 1, dummy_int * 5
                read(unit) dummy_int_aux
            end do

            ! Read number of cells
            read(unit) num_cells
            read(unit) num_components

            ! Allocate the field data array with the exact number of elements
            IF (ALLOCATED(field_data)) THEN
                DEALLOCATE(field_data)
            END IF
            allocate(field_data(num_components, num_cells))

            ! Read cell data
            read(unit) field_data

            ! Close the file
            close(unit)

        else
            print *, 'Unrecognised input file format: ', iformat
            stop
        end if

    end subroutine read_vtk_vectorfield

    subroutine read_vtk_scalarfield(filename, field_data)
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=RP), allocatable, intent(out) :: field_data(:)
        integer :: i, ios, unit, num_cells, num_components, num_read
        character(len=256) :: line, field_name, field_type
        integer :: current_index
        integer :: dummy_int, dummy_int_aux
        real(kind=RP) :: dummy_real

        if (iformat == 'ascii') then
    
            ! Open the file
            open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening file: ', filename
                stop
            end if
        
            ! Read until the field data section
            do
                read(unit, '(A)', iostat=ios) line
                if (index(line, 'FIELD') > 0) exit
            end do
            read(unit, '(A)', iostat=ios) line
            read(line, *) field_name, num_components, num_cells, field_type
        
            ! Allocate the field data array with the exact number of elements
            IF (ALLOCATED(field_data)) THEN
                DEALLOCATE(field_data)
            END IF
            allocate(field_data(num_cells))
    
        
            ! Read the field data
            current_index = 1
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                read(line, *, iostat=ios) (field_data(current_index + i - 1), i = 1, 10)
                if (ios == 0) then
                    current_index = current_index + count(field_data(current_index:current_index + 9) /= 0.0)
                else
                    exit
                end if
            end do
        
            ! Close the file
            close(unit)

        else if (iformat == 'binary') then

            ! Open the binary input file
            open(newunit=unit, file=filename, status='old', form='unformatted', access='stream', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening input file'
                stop
            end if

            ! Skip points
            read(unit) dummy_int
            do i = 1, dummy_int * 3
                read(unit) dummy_real
            end do

            ! Skip polygons
            read(unit) dummy_int
            do i = 1, dummy_int * 5
                read(unit) dummy_int_aux
            end do

            ! Read number of cells
            read(unit) num_cells
            read(unit) num_components

            ! Allocate the field data array with the exact number of elements
            IF (ALLOCATED(field_data)) THEN
                DEALLOCATE(field_data)
            END IF
            allocate(field_data(num_cells))

            ! Read cell data
            read(unit) field_data

            ! Close the file
            close(unit)

        else
            print *, 'Unrecognised input file format: ', iformat
            stop
        end if

    end subroutine read_vtk_scalarfield

    function time_interp_weights(t1, t2, t) result(alpha)
        implicit none
        real(kind=RP), intent(in) :: t1, t2, t
        real(kind=RP) :: alpha
        real(kind=RP) :: t_lower, t_upper
      
        ! Determine the lower and upper bounds
        if (t1 < t2) then
            t_lower = t1
            t_upper = t2
        else
            t_lower = t2
            t_upper = t1
        end if
      
        ! Calculate the interpolation factor alpha
        alpha = (t - t_lower) / (t_upper - t_lower)
    end function time_interp_weights

    function normalize(value, min_val, max_val, new_min, new_max) result(norm_value)
        implicit none
        real(kind=RP), intent(in) :: value
        real(kind=RP), intent(in) :: min_val, max_val, new_min, new_max
        integer :: norm_value
        norm_value = int((value - min_val) / (max_val - min_val) * (new_max - new_min) + new_min)
    end function normalize
    
    function z_order(x, y, k) result(z)
        implicit none
        integer, intent(in) :: x, y, k
        integer :: z
        integer :: i, mask_x, mask_y
    
        z = 0
        do i = 0, k-1
            mask_x = ibits(x, i, 1)
            mask_y = ibits(y, i, 1)
            z = z + ishft(mask_x, 2*i) + ishft(mask_y, 2*i + 1)
        end do
    end function z_order

    subroutine find_closest_point(n, src, tar, srci)

        integer, intent(in) :: n
        real(kind=RP), intent(in) :: src(2, n)
        real(kind=RP), intent(in) :: tar(2)
        integer(kind=4), intent(inout) :: srci
        integer, parameter :: dp = kind(1.0d0)
        integer :: i
        real(kind=RP) :: dist
        real(dp) :: min_dist 
        real(kind=RP) :: dx, dy
      
        ! Initialize minimum distance to a large value
        min_dist = huge(1.0_dp)
        srci = -1
      
        ! Find the closest point
        do i = 1, n
          dx = src(1, i) - tar(1)
          dy = src(2, i) - tar(2)
          dist = dx*dx + dy*dy
          if (dist < min_dist) then
            min_dist = dist
            srci = i
          end if
        end do
      
    end subroutine find_closest_point

    subroutine compute_affine_transform(v1, v2, v3, transform)
        real(kind=RP), intent(in) :: v1(2), v2(2), v3(2)  ! Vertices of the simplex
        real(kind=RP), intent(out) :: transform(3, 2)     ! Transform matrix
        
        real(kind=RP) :: T(2, 2), T_inv(2, 2)
        real(kind=RP) :: r(2)
        integer(kind=4) :: i
        
        ! Compute the matrix T
        do i = 1, 2
            T(i, 1) = v2(i) - v1(i)
            T(i, 2) = v3(i) - v1(i)
        end do
        
        ! Compute the inverse of T
        call invert_2x2(T, T_inv)
        
        ! Compute the vector r
        r = v1
        
        ! Fill the transform matrix
        transform(1:2, 1:2) = T_inv
        transform(3, 1:2) = r
    end subroutine compute_affine_transform

    subroutine invert_2x2(A, A_inv)
        real(kind=RP), intent(in) :: A(2, 2)
        real(kind=RP), intent(out) :: A_inv(2, 2)
        real(kind=RP) :: det
      
        ! Compute the determinant
        det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
      
        ! Check for singular matrix
        if (det == 0.0) then
            print *, "Matrix is singular and cannot be inverted"
            stop
        end if
      
        ! Compute the inverse
        A_inv(1, 1) = A(2, 2) / det
        A_inv(1, 2) = -A(1, 2) / det
        A_inv(2, 1) = -A(2, 1) / det
        A_inv(2, 2) = A(1, 1) / det
    end subroutine invert_2x2

    subroutine compute_barycentric_weights(nptar, transform, xy_target, weights)
        implicit none
        integer, parameter :: ndim = 2
        integer(kind=4), intent(in) :: nptar
        real(kind=RP), intent(in) :: transform(nptar, ndim+1, ndim)
        real(kind=RP), intent(in) :: xy_target(ndim, nptar)
        real(kind=RP), intent(out) :: weights(nptar, ndim+1)
        real(kind=RP) :: bary_coords(nptar, ndim)
      
        integer :: i, j, k
        real(kind=RP) :: diff(ndim)
      
        ! Initialize bary_coords and weights
        bary_coords = 0.0
        weights = 0.0
      
        ! Loop over each target point
        do i = 1, nptar
          bary_coords(i, :) = matmul(transform(i, :ndim, :), xy_target(:, i)-transform(i, ndim+1, :))   
          ! Compute weights
          do j = 1, ndim
            weights(i, j) = bary_coords(i, j)
          end do
          weights(i, ndim+1) = 1.0 - sum(bary_coords(i, :))
        end do
    end subroutine compute_barycentric_weights

End Module precursor