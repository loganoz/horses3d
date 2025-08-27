program ascii_to_binary_vtk
    implicit none
    character(len=256) :: arg, file_path, output_path
    integer :: ios, argc
    character(len=256) :: line, dummy
    integer :: num_points, num_polygons, total_entries, num_cell_data, num_components, i
    real(kind=8), allocatable :: points(:,:), cell_data(:,:)
    integer, allocatable :: polygons(:,:)
    real(kind=8) :: x, y, z
    integer :: poly(5), current_index

    ! Get the number of command-line arguments
    argc = 0
    do
        call get_command_argument(argc, arg)
        if (len_trim(arg) == 0) exit

        write (*, *) trim(arg)
        argc = argc + 1
    end do

    ! Check if the correct number of arguments is provided
    if (argc /= 3) then
        print *, 'Usage: ascii_to_binary_vtk <file_path> <output_path>'
        stop
    end if

    ! Get the command-line arguments
    call get_command_argument(1, file_path)
    call get_command_argument(2, output_path)

    ! Open the ASCII input file
    open(unit=10, file=file_path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening input file'
        stop
    end if

    ! Open the binary output file
    open(unit=20, file=output_path, status='replace', form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening output file'
        stop
    end if

    ! Read file info until points header
    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (trim(line) == 'DATASET POLYDATA') exit
    end do

    ! Read points header
    read(10, '(A)', iostat=ios) line
    read(line, *) dummy, num_points, dummy

    ! Allocate array for points
    allocate(points(3, num_points))

    ! Read points
    do i = 1, num_points
        read(10, *) x, y, z
        points(:, i) = (/ x, y, z /)
    end do

    ! Read polygons header
    read(10, '(A)', iostat=ios) line  ! Empty line
    read(10, '(A)', iostat=ios) line
    read(line, *) dummy, num_polygons, total_entries

    ! Allocate array for polygons
    allocate(polygons(5, num_polygons))

    ! Read polygons
    do i = 1, num_polygons
        read(10, *) poly
        polygons(:, i) = poly
    end do

    ! Read cell data headers
    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (trim(line) == 'FIELD attributes 1') exit
    end do

    ! Read cell data field header
    read(10, '(A)', iostat=ios) line
    read(line, *) dummy, num_components, num_cell_data, dummy

    ! Allocate array for cell data
    allocate(cell_data(num_components, num_cell_data))

    ! Read cell data
    if (num_components > 1) then
        do i = 1, num_cell_data
            read(10, *) cell_data(:, i)
        end do
    else if (num_components == 1) then
        current_index = 1
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
            read(line, *, iostat=ios) (cell_data(1, current_index + i - 1), i = 1, 10)
            if (ios == 0) then
                current_index = current_index + count(cell_data(1, current_index:current_index + 9) /= 0.0)
            else
                exit
            end if
        end do
    else 
        print *, 'Unexpected number of components in cell data'
        stop
    end if

    ! Write points to binary file
    write(20) num_points
    write(20) points

    ! Write polygons to binary file
    write(20) num_polygons
    write(20) polygons

    ! Write cell data to binary file
    write(20) num_cell_data
    write(20) num_components
    write(20) cell_data

    ! Close the files
    close(10)
    close(20)

    print *, 'Conversion to binary complete.'


end program ascii_to_binary_vtk