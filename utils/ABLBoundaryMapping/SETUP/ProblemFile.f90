!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva
!  
!      Modified: November 21, 2024 at 11:33 AM
!      By: Eduardo Jané
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
!
!
!//////////////////////////////////////////////////////////////////////// 
!
!
!//////////////////////////////////////////////////////////////////////// 
!
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
            use SMConstants
            use precursor
            IMPLICIT NONE

            integer :: i, j, ierror
            real :: centroid(3), point(3)

            ! Read precursor input file
            call read_wind_precursor_table()
            ! print *, 'Number of directories: ', num_dirs
            ! print *, 'File names: ', trim(file_names(1)),', ', trim(file_names(2)),', ',trim(file_names(3))
            ! do i = 1, num_dirs
            !     print *, 'Time: ', times(i), ', Directory: ', trim(directories(i))
            ! end do

            ! Load mesh from first VTK file in list (assuming that all files have the same discretisation!)
            dirname = 'PRECURSOR'
            filename = trim(dirname) // '/' // trim(directories(1)) // '/' // trim(file_names(1))
            call read_vtk_mesh(filename, num_points, points, num_polygons, polygons)

            ! Allocate memory
            allocate(centroids(2, num_polygons))   
            allocate(ind(num_polygons))
            allocate(tnbr(3, num_polygons*2))
            allocate(til(3, num_polygons*2))

            ! Compute centroid points for each polygon
            ! PRINT *, 'Computing polygon centroids'
            DO i = 1, num_polygons
            centroid(:) = (/ 0.0, 0.0, 0.0 /)
            DO j = 2, 5
                point = points(:, polygons(j, i) + 1) ! VTK indexes points from 0
                centroid(1) = centroid(1) + point(1)
                centroid(2) = centroid(2) + point(2)
                centroid(3) = centroid(3) + point(3)
            END DO
            centroid(1) = centroid(1) / 4.0
            centroid(2) = centroid(2) / 4.0
            centroid(3) = centroid(3) / 4.0
            centroids(:, i) = centroid(2:) ! in this case, X is constant
            ind(i) = i
            END DO

            ! Call the dtris2 subroutine to compute the Delaunay triangulation
            ! PRINT *, 'Performing 2D Delaunay tesselation'
            CALL dtris2(num_polygons, centroids, ind, ntri, til, tnbr, ierror)

         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use precursor
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
#if defined(NAVIERSTOKES)
            integer :: i, j, k, l, n_bp, b_id, f_id, ix, iy, ierror, itri, iedg, swap
            real(kind=RP), dimension(:,:,:), allocatable :: xyz
            real(kind=RP), allocatable :: bpoints(:, :)
            real(kind=RP) :: transform(3, 2)
!
!           Build inlet points hash table, connectivity and interpolation weights matrices
!           ------------------------------------------------------------------------------
            b_id = 0
            do i = 1, size(mesh % zones)
                if (trim(mesh % zones(i) % name) == 'inlet') then ! TO DO: Name of the BC should be an input in the precursor driver file
                    b_id = i
                    exit
                end if
            end do

            n_bp = 0
            do i = 1, mesh % zones(b_id) % no_of_faces
                f_id = mesh % zones(b_id) % faces(i) 
                xyz = mesh % faces(f_id) % geom % x 
                do j = 0, size(xyz, 2) - 1
                    do k = 0, size(xyz, 3) - 1
                        n_bp = n_bp + 1
                    end do
                end do
            end do
            ! print *, 'Number of points in boundary ', trim(mesh % zones(b_id) % name), ': ', n_bp

            if (n_bp > 0) then

                allocate(bpoints(2, n_bp))
                l = 1
                do i = 1, mesh % zones(b_id) % no_of_faces
                    f_id = mesh % zones(b_id) % faces(i)
                    xyz = mesh % faces(f_id) % geom % x 
                    ! print *, '1st dim: ', size(xyz, 1)
                    ! print *, '2nd dim: ', size(xyz, 2)
                    ! print *, '3rd dim: ', size(xyz, 3)
                    do j = 0, size(xyz, 2) - 1
                        do k = 0, size(xyz, 3) - 1
                            bpoints(1, l) = xyz(2, j, k) ! hard coded for yz plane
                            bpoints(2, l) = xyz(3, j, k)
                            l = l + 1
                            ! print *, 'Boundary point (x, y, z) for (j, k): ', xyz(:, j, k), j, k
                        end do
                    end do
                end do

                max_x = maxval(bpoints(1,:)) ! Hard coded for yz plane -> vals needed for normalization in hash table
                min_x = minval(bpoints(1,:))
                max_y = maxval(bpoints(2,:))
                min_y = minval(bpoints(2,:))

                allocate(hashtable(n_bp))
                allocate(simplices(n_bp))
                allocate(outside(n_bp))

                itri = 1
                iedg = 1
                do i = 1, n_bp
                    ix = normalize(bpoints(1, i), min_x, max_x, 0.0_RP, real(2**RP - 1, kind=RP))
                    iy = normalize(bpoints(2, i), min_y, max_y, 0.0_RP, real(2**RP - 1, kind=RP))
                    hashtable(i) = z_order(ix, iy, RP)
                    ! print *, 'Preprocessed (x,y) -> hash: ', bpoints(:, i), ' -> ', hashtable(i)
                    call walkt2(bpoints(1, i), bpoints(2, i), ntri, centroids, til, tnbr, itri, iedg, ierror)
                    simplices(i) = itri
                    if (iedg < 0) then
                        call find_closest_point(num_polygons, centroids, bpoints(:, i), outside(i))
                    else
                        outside(i) = 0
                    end if
                end do

                ! Compute Affine transform matrix (from cartesian to local barycentric coordinates) for each tesselation triangle
                ! PRINT *, 'Computing cartesian-to-barycentric affine transforms for tessellation'
                
                ! Dynamically allocate the array
                allocate(Tmat(ntri, 3, 2))

                DO i = 1, ntri
                    CALL compute_affine_transform(centroids(:, til(1, i)), centroids(:, til(2, i)), centroids(:, til(3, i)), transform)
                    Tmat(i, :, :) = transform   
                    ! Re-arrange til of current triangle after transform (first vertex becomes last to match transform definition)
                    swap = til(1, i)
                    til(1, i) = til(2, i)
                    til(2, i) = til(3, i)
                    til(3, i) = swap
                END DO

                ! Compute the interpolation weights for each target point
                ! PRINT *, 'Computing interpolation weights matrix'
                allocate(weights(n_bp, 3))
                CALL compute_barycentric_weights(n_bp, Tmat(simplices, :, :), bpoints, weights)

                ! Loop over each target point
                ! PRINT *, 'Building connectivity matrix'
                allocate(conn_mat(n_bp, 3))
                do i = 1, n_bp
                    ! Get the connectivity matrix for the current simplex
                    conn_mat(i, :) = til(:, simplices(i))
                end do

            end if

#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            use precursor
#if defined(_HAS_MPI_)
            use mpi
#endif
#if defined(NAVIERSTOKES)
#endif
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
integer        :: eid, i, j, k
real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
real(kind=RP)  :: Q(NCONS), phi, theta
#endif
#if defined(_HAS_MPI_)
integer        :: irank, ierr
#endif

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
    associate ( gammaM2 => dimensionless_ % gammaM2, &
                gamma => thermodynamics_ % gamma )
    theta = refvalues_ % AOAtheta*(pi/180.0_RP)
    phi   = refvalues_ % AOAphi*(pi/180.0_RP)

    do eID = 1, mesh % no_of_elements
    associate( Nx => mesh % elements(eID) % Nxyz(1), &
                ny => mesh % elemeNts(eID) % nxyz(2), &
                Nz => mesh % elements(eID) % Nxyz(3) )
    do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
        qq = 1.0_RP
        u  = qq*cos(theta)*cos(phi)
        v  = qq*sin(theta)*cos(phi)
        w  = qq*sin(phi)

        q(1) = 1.0_RP
        p    = 1.0_RP/(gammaM2)
        q(2) = q(1)*u
        q(3) = q(1)*v
        q(4) = q(1)*w
        q(5) = p/(gamma - 1._RP) + 0.5_RP*q(1)*(u**2 + v**2 + w**2)

        mesh % elements(eID) % storage % q(:,i,j,k) = q 
    end do;        end do;        end do
    end associate
    end do

    end associate

#endif  

#if defined(_HAS_MPI_)
    call mpi_comm_rank(MPI_COMM_WORLD, irank, ierr)

    if (irank == 0) then
        write(*, *)
        write(*, *)
        write(STD_OUT,'(15X,A)') "Boundary precursor"
        write(STD_OUT,'(15X,A)') "------------------"
        flush(6)
        write(STD_OUT,'(30X,A,A30,A)') "->", "Boundary name: ", "inlet"
        write(STD_OUT,'(30X,A,A30,I0)') "->", "Precursor snapshots: ", num_dirs
        write(STD_OUT,'(30X,A,A30,A,F7.3,A,F7.3,A)') "->", "Snapshots time range: ", "[", times(1), ", ", times(num_dirs), "]"
        write(STD_OUT,'(30X,A,A30,I0)') "->", "Source points: ", num_polygons
        write(STD_OUT,'(30X,A,A30,I0)') "->", "2D Delaunay elements: ", ntri
    end if
#endif

         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            use precursor
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_

            real(kind=RP) :: alpha, u, v, w, temp, p, rho, t_dim, L_ref
            real(kind=RP) :: rho_0, p_0, temp_ref, patm
            real(kind=RP) :: u_nd, v_nd, w_nd, p_nd, rho_nd
            integer :: i, j, i_tl, i_tu, ix_norm, iy_norm, hash
            integer, dimension(:), allocatable :: indices

            associate( gamma   => thermodynamics_ % gamma , &
                       R       => thermodynamics_ % R , &
                       gammaM2 => dimensionless_ % gammaM2, &
                       uref    => refValues_ % V, &
                       pref    => refValues_ % p, & ! pref here is rhorhef*uref**2
                       Tref    => refValues_ % T, &
                       rhoref  => refValues_ % rho)

                ! Compute time with dimensions
                L_ref = 1.0_RP ! Default reference length is 1 m
                t_dim = t / uref * L_ref

                ! Find the interval of source files for time-interpolation
                i_tl = 0
                i_tu = 0
                do i = 1, num_dirs-1
                    if ((times(i) <= t_dim) .and. (t_dim <= times(i+1))) then
                        i_tl = i
                        i_tu = i+1
                        exit
                    end if
                end do

                ! Stay in final interval if exceeded
                if (t_dim > times(num_dirs)) then
                    i_tl = num_dirs-1
                    i_tu = num_dirs 
                    alpha = 1.0_RP
                else
                    ! Get time-interpolation weights
                    alpha = time_interp_weights(times(i_tl), times(i_tu), t_dim)
                end if

                ! ! Get time-interpolation weights
                ! alpha = time_interp_weights(times(i_tl), times(i_tu), t_dim)

                if((i_tl_prev /= i_tl) .and. (i_tu_prev /= i_tu)) then
                    ! Read source data from VTK files
                    ! PRINT *, 'Reading input *.vtk polydata'
                    dirname = 'PRECURSOR'
                    filename = trim(dirname) // '/' // trim(directories(i_tl)) // '/' // trim(file_names(1))
                    call read_vtk_vectorfield(filename, field_data_U_l)
                    ! PRINT *, 'Read lower limit U file ', trim(filename)
                    filename = trim(dirname) // '/' // trim(directories(i_tl)) // '/' // trim(file_names(2))
                    call read_vtk_scalarfield(filename, field_data_T_l)
                    ! PRINT *, 'Read lower limit T file ', trim(filename)
                    filename = trim(dirname) // '/' // trim(directories(i_tl)) // '/' // trim(file_names(3))
                    call read_vtk_scalarfield(filename, field_data_p_l)
                    ! PRINT *, 'Read lower limit p file ', trim(filename)
        
                    filename = trim(dirname) // '/' // trim(directories(i_tu)) // '/' // trim(file_names(1))
                    call read_vtk_vectorfield(filename, field_data_U_u)
                    ! PRINT *, 'Read upper limit U file ', trim(filename)
                    filename = trim(dirname) // '/' // trim(directories(i_tu)) // '/' // trim(file_names(2))
                    call read_vtk_scalarfield(filename, field_data_T_u)
                    ! PRINT *, 'Read upper limit T file ', trim(filename)
                    filename = trim(dirname) // '/' // trim(directories(i_tu)) // '/' // trim(file_names(3))
                    call read_vtk_scalarfield(filename, field_data_p_u)
                    ! PRINT *, 'Read upper limit p file ', trim(filename)
                end if
        
                i_tl_prev = i_tl
                i_tu_prev = i_tu

                ! Interpolate U from source to target points
                ix_norm = normalize(x(2), min_x, max_x, 0.0_RP, real(2**RP - 1, kind=RP))
                iy_norm = normalize(x(3), min_y, max_y, 0.0_RP, real(2**RP - 1, kind=RP))
                hash = z_order(ix_norm, iy_norm, RP)
                indices = findloc(hashtable, hash)
                ! print *, 'hashtable, hash:', hashtable, ',', hash
                ! print *, 'indices: ', indices
                i = indices(1)
                ! print *, 'Boundary (x,y) -> hash -> index in hashtable: ', x(2), x(3), ' -> ', hash, ' -> ', i

                if (outside(i) == 0) then
                    ! print *, 'Point with hash ', hash, 'and index in hashtable ', i, ' is not outside source plane' 
                    u = 0.0
                    v = 0.0
                    w = 0.0
                    temp = 0.0
                    p = 0.0
                    do j = 1, 3 !(ndim+1)
                        ! print *, 'Adding value of simplex node ', j
                        ! print *, 'Weight ', weights(i, j)
                        ! print *, 'Node in source data ', conn_mat(i, j)
                        u = u + weights(i, j) * ((1.0 - alpha)*field_data_U_l(1, conn_mat(i, j)) + alpha*field_data_U_u(1, conn_mat(i, j)))
                        v = v + weights(i, j) * ((1.0 - alpha)*field_data_U_l(2, conn_mat(i, j)) + alpha*field_data_U_u(2, conn_mat(i, j)))
                        w = w + weights(i, j) * ((1.0 - alpha)*field_data_U_l(3, conn_mat(i, j)) + alpha*field_data_U_u(3, conn_mat(i, j)))
                        temp = temp + weights(i, j) * ((1.0 - alpha)*field_data_T_l(conn_mat(i, j)) + alpha*field_data_T_u(conn_mat(i, j)))
                        p = p + weights(i, j) * ((1.0 - alpha)*field_data_p_l(conn_mat(i, j)) + alpha*field_data_p_u(conn_mat(i, j)))
                    end do
                else
                    ! print *, 'Point with hash ', hash, 'and index in hashtable ', i, ' is outside source plane'
                    u = (1.0 - alpha)*field_data_U_l(1, outside(i)) + alpha*field_data_U_u(1, outside(i))
                    v = (1.0 - alpha)*field_data_U_l(2, outside(i)) + alpha*field_data_U_u(2, outside(i))
                    w = (1.0 - alpha)*field_data_U_l(3, outside(i)) + alpha*field_data_U_u(3, outside(i))
                    temp = (1.0 - alpha)*field_data_T_l(outside(i)) + alpha*field_data_T_u(outside(i))
                    p = (1.0 - alpha)*field_data_p_l(outside(i)) + alpha*field_data_p_u(outside(i))
                end if

                ! print *, 'Interpolated u at (y,z) point ', x(2), ',', x(3), ': ', u
                ! print *, 'Interpolated v at (y,z) point ', x(2), ',', x(3), ': ', v
                ! print *, 'Interpolated w at (y,z) point ', x(2), ',', x(3), ': ', w

                temp_ref = 300.0_RP ! This should be a constant specified in the input file
                rho_0 = 1.177_RP ! density of air at 300K
                rho = rho_0 * (1.0 - ( (temp - temp_ref)/temp_ref )) ! Boussinesq buoyancy density

                ! Compute dimensionsless variables
                patm = R * rhoref * Tref ! As defined in PhysicsStorage_NS.f90 for reference
                rho_nd = rho / rhoref
                ! rho_nd = 1.0_RP
                u_nd = u / uref
                v_nd = v / uref
                w_nd = w / uref
                !p_nd = 1.0_RP / (gammaM2) ! This is 1atm/pref with pref = rho*uref**2
                p_nd = (p + patm) / pref ! Add atmospheric pressure and divide by reference value. In this case, source 'p' is not kinematic, so there is no need to multiply by rho. 

                !print *, 'rho, rhoref, p, pref: ', rho, rhoref, p, pref

                Q(1) = rho_nd
                Q(2) = rho_nd*u_nd
                Q(3) = rho_nd*v_nd
                Q(4) = rho_nd*w_nd
                Q(5) = p_nd/(gamma - 1.0_RP) + 0.5_RP*rho_nd*(u_nd*u_nd + v_nd*v_nd + w_nd*w_nd)

            end associate

         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           -------------------------------------------------------------------------
!           Called before every time-step to allow periodic operations to be performed
!           Here:
!              * Reading mean momentum/velocity in the volume
!              * Obtaining fTurbulentChannel according to it
!           -------------------------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
            use precursor
#endif
            IMPLICIT NONE
            !-arguments---------------------------------------------------
            CLASS(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
#if defined(NAVIERSTOKES)

#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_) ! , Q, dt
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            !-arguments--------------------------------------------------
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(out) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
            
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            real(kind=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime


         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
