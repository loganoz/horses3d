#include "Includes.h"

module ActuatorLine
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
    use SMConstants
    use MPI_Process_Info
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

private
public farm, ConstructFarm, DestructFarm, UpdateFarm, ForcesFarm, WriteFarmForces
!
!  ******************************
!  DEFINE TURBINE, BLADE, AIRFOIL
!  ******************************
! 
!
    type airfoil_t
    integer                         :: num_aoa
    integer                         :: num_Re 
    real(KIND=RP), allocatable      :: aoa(:)  ! in rad
    real(KIND=RP), allocatable      :: Re(:)  ! reynolds number
    real(KIND=RP), allocatable      :: cl(:,:)
    real(KIND=RP), allocatable      :: cd(:,:)
    end type

    type blade_t
    real(KIND=RP), allocatable      :: r_R(:)  ! in m
    real(KIND=RP), allocatable      :: chord(:)   ! in m
    real(KIND=RP), allocatable      :: twist(:)    ! in rad
    real(KIND=RP)                   :: azimuth_angle   ! in rad
    integer, allocatable            :: num_airfoils(:)    ! for different Re
    CHARACTER(LEN=30), allocatable  :: airfoil_files(:,:)  ! file names for Cl-Cd
    type(airfoil_t), allocatable    :: airfoil_t(:)  ! airfoil data AoA-Cl-Cd
    real(KIND=RP), allocatable      :: local_velocity(:)  ! local flow speed at blade section, im m/s
    real(KIND=RP), allocatable      :: local_angle(:)  ! local flow angle at blade section, in rad
    real(KIND=RP), allocatable      :: point_xyz_loc(:,:) ! x,y location of blade points
    real(KIND=RP), allocatable      :: local_torque(:)  ! N.m
    real(KIND=RP), allocatable      :: local_thrust(:)  ! N
    real(KIND=RP), allocatable      :: local_thrust_temp(:)  ! N
    real(KIND=RP), allocatable      :: local_rotor_force(:)  ! N
    real(KIND=RP), allocatable      :: local_rotor_force_temp(:)  ! N
    real(KIND=RP), allocatable      :: local_velocity_temp(:)
    real(KIND=RP), allocatable      :: local_angle_temp(:)
    real(KIND=RP), allocatable      :: local_root_bending(:)  ! N.m
    real(KIND=RP), allocatable      :: gauss_epsil_delta(:)  ! HO element size, for calculate force Gaussian
    real(KIND=RP)                   :: tip_c1,tip_c2  ! tip force corrections
    real(KIND=RP), allocatable      :: local_gaussian_sum(:) ! necessary for Gaussian weighted average
    real(KIND=RP), allocatable      :: local_Re(:) ! local Re based on local conditions and the chord of the airfoil at the blade section
    real(KIND=RP), allocatable      :: local_Re_temp(:)
    end type

    type turbine_t
    integer                        :: num_blades=3  ! number of blades -> hardcoded 3 blades
    real(KIND=RP)                  :: radius ! turb radius in mm
    real(KIND=RP)                  :: blade_pitch ! turb radius in rad
    real(KIND=RP)                  :: rot_speed ! rad/s
    real(KIND=RP)                  :: hub_cood_x,hub_cood_y,hub_cood_z ! hub height in mm
    real(KIND=RP)                  :: normal_x, normal_y, normal_z ! rotor normal pointing backward
    type(blade_t)                  :: blade_t(3) ! hardcoded 3 blades
    integer                        :: num_blade_sections  ! number of 2D section for Cl-Cd data
    real(KIND=RP)                  :: blade_torque(3) ! N.m
    real(KIND=RP)                  :: blade_thrust(3) ! N
    real(KIND=RP)                  :: blade_root_bending(3) !N.m
    real(KIND=RP)                  :: Cp ! turbine power coef.
    real(KIND=RP)                  :: Ct ! turbine thrust coef.
    real(KIND=RP), allocatable     :: average_conditions(:,:) ! time and blade average local variables at the blade section
    end type
                                   
    type Farm_t
    integer                        :: num_turbines
    type(turbine_t), allocatable   :: turbine_t(:)
    real(KIND=RP)                  :: gauss_epsil  ! force Gaussian shape
    real(KIND=RP)                  :: time
    real(KIND=RP)                  :: tolerance_factor
    integer                        :: epsilon_type
    logical                        :: calculate_with_projection
    logical                        :: active = .false.
    logical                        :: save_average = .false.
    logical                        :: save_instant = .false.
    logical                        :: verbose = .false.
    logical                        :: averageSubElement = .true.
    character(len=LINE_LENGTH)     :: file_name
    integer                        :: number_iterations
    integer                        :: save_iterations

    end type

    type(Farm_t)                   :: farm

    ! max 10 airfoils file names per section
    integer, parameter             :: MAX_AIRFOIL_FILES = 10
    integer, dimension(:), allocatable         :: numElementsPerTurbine, elementsActuated, turbineOfElement

!  ========
contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructFarm(self, controlVariables, t0, mesh)
       use FTValueDictionaryClass
       use mainKeywordsModule, only: solutionFileNameKey, restartFileNameKey
       use FileReadingUtilities      , only: getFileName
       use PhysicsStorage
       use HexMeshClass
       use fluiddata
       use MPI_Process_Info
       implicit none
       class(farm_t) , intent(inout)                :: self
       TYPE(FTValueDictionary), intent(in)          :: controlVariables
       real(kind=RP), intent(in)                    :: t0
       type(HexMesh), intent(in)                    :: mesh
!        ---------------
!        Local variables
!        ---------------
!
         integer     ::  i, j, k, ii, fid, n_aoa, n_airfoil
         CHARACTER(LEN=LINE_LENGTH) :: arg, char1
         CHARACTER(LEN=LINE_LENGTH) :: solution_file
         CHARACTER(LEN=5)           :: file_id
         real(kind=RP), dimension(:), allocatable   :: initial_azimutal
         character(len=STRING_CONSTANT_LENGTH)  :: restart_name, restart_operations_name
         logical                    :: fileExists
         integer        :: nelem, eID, eIndex
         real(kind=RP)  :: tolerance, r_square
         real(kind=RP)  :: delta, delta_temp
         integer        :: delta_count, delta_paritions, ierr

    if (.not. controlVariables % logicalValueForKey("use actuatorline")) return

    self % time = t0

    self % epsilon_type = controlVariables % getValueOrDefault("actuator epsilon type", 0)
    self % calculate_with_projection = controlVariables % getValueOrDefault("actuator calculate with projection", .false.)
    self % save_average = controlVariables % getValueOrDefault("actuator save average", .false.)
    self % save_instant = controlVariables % getValueOrDefault("actuator save instant", .false.)
    self % save_iterations = controlVariables % getValueOrDefault("actuator save iteration", 1)
    self % verbose = controlVariables % getValueOrDefault("actuator verbose", .false.)
    self % averageSubElement = controlVariables % getValueOrDefault("actuator average subelement", .true.)
    self % tolerance_factor = controlVariables % getValueOrDefault("actuator tolerance", 0.2_RP)

    restart_name = controlVariables % stringValueForKey( restartFileNameKey, requestedLength = STRING_CONSTANT_LENGTH )
    restart_name = trim(getFileName(restart_name))
    write(restart_operations_name,'(2A)')  TRIM(restart_name),'_Actuator_Line_operations.dat'
    inquire(file=trim(restart_operations_name), exist=fileExists)

    arg='./ActuatorDef/Act_ActuatorDef.dat'
    OPEN( newunit = fid,file=trim(arg),status="old",action="read")

    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,*) self%num_turbines

    if (self % verbose .and. MPI_Process % isRoot) then
        print *,'-------------------------'
        print *,achar(27)//'[34m READING FARM DEFINITION'
        write(*,*) "Number of turbines in farm:", self%num_turbines
    endif

    READ(fid,'(A132)') char1

    allocate(self%turbine_t(self%num_turbines))

    do k = 1, self%num_turbines
       READ(fid,*) self%turbine_t(k)%hub_cood_x, self%turbine_t(k)%hub_cood_y, self%turbine_t(k)%hub_cood_z
    ENDDO

    READ(fid,'(A132)') char1

    do k = 1, self%num_turbines
       READ(fid,*) self%turbine_t(k)%radius
    ENDDO

    READ(fid,'(A132)') char1

    do k = 1, self%num_turbines
       READ(fid,*) self%turbine_t(k)%normal_x, self%turbine_t(k)%normal_y, self%turbine_t(k)%normal_z
    ENDDO
    
    READ(fid,'(A132)') char1

    do k = 1, self%num_turbines
       READ(fid,*) self%turbine_t(k)%rot_speed
    ENDDO

    READ(fid,'(A132)') char1

    do k = 1, self%num_turbines
       READ(fid,*) self%turbine_t(k)%blade_pitch
    ENDDO
    
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1

    ! Read blade info, we assume all 3 blades are the same for one turbine
    do k = 1, self%num_turbines
        READ(fid,*) self%turbine_t(k)%num_blade_sections
    enddo

    do k=1, self%num_turbines
     associate (num_blade_sections => self%turbine_t(k)%num_blade_sections)

      do j=1, self%turbine_t(k)%num_blades
         allocate( self%turbine_t(k)%blade_t(j)%r_R(0:num_blade_sections),self%turbine_t(k)%blade_t(j)%chord(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%twist(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%num_airfoils(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%airfoil_files(num_blade_sections,MAX_AIRFOIL_FILES),self%turbine_t(k)%blade_t(j)%airfoil_t(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%local_velocity(num_blade_sections), self%turbine_t(k)%blade_t(j)%local_angle(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%point_xyz_loc(num_blade_sections,3),self%turbine_t(k)%blade_t(j)%local_torque(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%local_thrust(num_blade_sections),self%turbine_t(k)%blade_t(j)%local_root_bending(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%local_rotor_force(num_blade_sections),self%turbine_t(k)%blade_t(j)%local_gaussian_sum(num_blade_sections), &
         self%turbine_t(k)%blade_t(j)%local_Re(num_blade_sections), self%turbine_t(k)%blade_t(j)%gauss_epsil_delta(num_blade_sections) )

         do i=1, num_blade_sections
            self%turbine_t(k)%blade_t(j)%airfoil_files(i,:)=' '
         enddo
      enddo
    endassociate
   enddo

    if (self%calculate_with_projection) then
        do k=1, self%num_turbines
          associate (num_blade_sections => self%turbine_t(k)%num_blade_sections)
          do j=1, self%turbine_t(k)%num_blades
             allocate( self%turbine_t(k)%blade_t(j)%local_thrust_temp(num_blade_sections), &
                       self%turbine_t(k)%blade_t(j)%local_rotor_force_temp(num_blade_sections) , &
                       self%turbine_t(k)%blade_t(j)%local_velocity_temp(num_blade_sections) , &
                       self%turbine_t(k)%blade_t(j)%local_angle_temp(num_blade_sections) , &
                       self%turbine_t(k)%blade_t(j)%local_Re_temp(num_blade_sections))
          enddo
          endassociate
       enddo
   end if

   ! now read each blade definition, 1 per turbine
   do k=1, self%num_turbines
       self%turbine_t(k)%blade_t(1)%r_R(0) = 0.0_RP
          READ(fid,'(A132)') char1 ! leave one comment line per turbine for easy reading of the dat file
       do i = 1, self%turbine_t(k)%num_blade_sections
          READ(fid,*) self%turbine_t(k)%blade_t(1)%r_R(i), self%turbine_t(k)%blade_t(1)%chord(i), &
                      self%turbine_t(k)%blade_t(1)%twist(i), self%turbine_t(k)%blade_t(1)%num_airfoils(i)
          
          ! one file per Re for each airfoil
          self%turbine_t(k)%blade_t(1)%airfoil_t(i)%num_Re = self%turbine_t(k)%blade_t(1)%num_airfoils(i)
                        
          do n_airfoil = 1, self%turbine_t(k)%blade_t(1)%num_airfoils(i)
              READ(fid,*) self%turbine_t(k)%blade_t(1)%airfoil_files(i,n_airfoil)
          enddo
       enddo
   enddo

! read numerical parameters
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1

     READ(fid,*) self%gauss_epsil

     READ(fid,'(A132)') char1
     READ(fid,*) self%turbine_t(1)%blade_t(1)%tip_c1,self%turbine_t(1)%blade_t(1)%tip_c2

     do k=1, self%num_turbines
       self%turbine_t(k)%blade_t(1)%tip_c1 = self%turbine_t(1)%blade_t(1)%tip_c1
       self%turbine_t(k)%blade_t(1)%tip_c2 = self%turbine_t(1)%blade_t(1)%tip_c2
    end do

    close(fid)

    if (self % verbose .and. MPI_Process % isRoot) then
        print *,achar(27)//'[34m END OF READING FARM DEFINITION'
    end if

    if (MPI_Process % isRoot) then

        call Subsection_Header("Actuator Line")
        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of turbines: ", self % num_turbines
        do k = 1, self%num_turbines
            write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of blade sections: ", self%turbine_t(k)%num_blade_sections
        end do

        select case (self % epsilon_type)
        case (0)
            write(STD_OUT,'(30X,A,A28,ES10.3)') "->", 'Fixed Epsilon value: ',self%gauss_epsil
        case (1)
            write(STD_OUT,'(30X,A,A)') "->", 'Epsilon calculated based on drag value'
        case (2)
            write(STD_OUT,'(30X,A)') 'Epsilon calculated based on element size and polynomial order'
            write(STD_OUT,'(30X,A,A28,F10.3)') "->", 'Constant for Epsilon: ',self%gauss_epsil
            if (self%calculate_with_projection) write(STD_OUT,'(30X,A)') 'Warining, epsilon calculated using properties of element 1'
        case default
            write(STD_OUT,'(30X,A,A28,ES10.3)') "->", 'Fixed Epsilon value: ',self%gauss_epsil
        end select
        write(STD_OUT,'(30X,A,A28,F10.3,F10.3)') "->", 'Tip correction constants: ', self%turbine_t(1)%blade_t(1)%tip_c1, self%turbine_t(1)%blade_t(1)%tip_c2
        write(STD_OUT,'(30X,A,A28,L1)') "->", "Projection formulation: ", self % calculate_with_projection
        if (.not. self%calculate_with_projection) write(STD_OUT,'(30X,A,A28,L1)') "->", "Average sub-Element: ", self % averageSubElement
        write(STD_OUT,'(30X,A,A28,L1)') "->", "Save blade average values: ", self % save_average
        if (fileExists)  write(STD_OUT,'(30X,A)') 'Using restaring operations of turbines'
    end if

   ! now read each airfoil, only for blade 1 of each turbine
    do k = 1, self%num_turbines
      do i = 1, self%turbine_t(k)%num_blade_sections
          arg=trim('./ActuatorDef/'//trim(self%turbine_t(k)%blade_t(1)%airfoil_files(i,1)))
          OPEN( newunit = fid,file=trim(arg),status="old",action="read")
          READ(fid,'(A132)') char1
          READ(fid,*) self%turbine_t(k)%blade_t(1)%airfoil_t(i)%num_aoa 
          close(fid)

          associate (num_re => self%turbine_t(k)%blade_t(1)%airfoil_t(i)%num_Re, num_aoa => self%turbine_t(k)%blade_t(1)%airfoil_t(i)%num_aoa)

             do j=1, self%turbine_t(k)%num_blades
               allocate( self%turbine_t(k)%blade_t(j)%airfoil_t(i)%aoa(num_aoa), &
                         self%turbine_t(k)%blade_t(j)%airfoil_t(i)%cl(num_aoa,num_re), &
                         self%turbine_t(k)%blade_t(j)%airfoil_t(i)%cd(num_aoa,num_re), &
                         self%turbine_t(k)%blade_t(j)%airfoil_t(i)%Re(num_re) )
             enddo

             do n_airfoil = 1, num_re

                arg=trim('./ActuatorDef/'//trim(self%turbine_t(k)%blade_t(1)%airfoil_files(i,n_airfoil)))
                OPEN( newunit = fid,file=trim(arg),status="old",action="read")
                READ(fid,'(A132)') char1

                READ(fid,*) n_aoa
                if ( n_aoa .ne. num_aoa ) then
                    print *, "Error: not same number of AoA in all files for same blade section, file: ", trim(arg)
                    call exit(99)
                end if

                READ(fid,'(A132)') char1
                READ(fid,*) self%turbine_t(k)%blade_t(1)%airfoil_t(i)%Re(n_airfoil)

                if (self % verbose .and. MPI_Process % isRoot) then
                    print *,'-------------------------'
                    print *,achar(27)//'[34m READING FARM AIRFOIL DATA (Cl-Cd)'
                    print*, 'reading: ', trim(arg)
                    write(*,*) 'The number of AoA in the file is: ', num_aoa,' '//achar(27)//'[0m '    
                end if 
    
                READ(fid,'(A132)') char1

                do ii = 1,  num_aoa
                     READ(fid,*) self%turbine_t(k)%blade_t(1)%airfoil_t(i)%aoa(ii), self%turbine_t(k)%blade_t(1)%airfoil_t(i)%cl(ii,n_airfoil), &
                                 self%turbine_t(k)%blade_t(1)%airfoil_t(i)%cd(ii,n_airfoil)
                     ! file is in deg, convert to rad
                     self%turbine_t(k)%blade_t(1)%airfoil_t(i)%aoa(ii) = self%turbine_t(k)%blade_t(1)%airfoil_t(i)%aoa(ii) * PI / 180.0_RP
                enddo

                close(fid)
             end do ! number of airfoil files

          endassociate

      enddo ! number of blade sections
    enddo ! number of turbines

    do k = 1, self%num_turbines
      do i = 1, self%turbine_t(k)%num_blade_sections
          self%turbine_t(k)%blade_t(1)%point_xyz_loc(i,1) =   self%turbine_t(k)%hub_cood_x
      enddo
    enddo

    !all blades of each turbine are the same
    do k=1, self%num_turbines
       do j=1, self%turbine_t(k)%num_blades 
          self%turbine_t(k)%blade_t(j)=self%turbine_t(k)%blade_t(1)
       enddo
    enddo

    allocate(initial_azimutal(self%num_turbines))
    if (fileExists) then
        initial_azimutal = 0.0_RP
        ! call self % readOperationFile
    else
        do k = 1, self%num_turbines
            initial_azimutal(k) = self%turbine_t(k)%rot_speed*t0 * Lref / refValues%V
        end do
    end if
    ! azimuthal angle for the 3 blades
    ! azimuth_angle angle of blades is the angle to respect to +y axis, the angular velocity vector will point to +x
    do k = 1, self%num_turbines
        self%turbine_t(k)%blade_t(1)%azimuth_angle = initial_azimutal(k)
        self%turbine_t(k)%blade_t(2)%azimuth_angle = initial_azimutal(k) + PI*2.0_RP/3.0_RP
        self%turbine_t(k)%blade_t(3)%azimuth_angle = initial_azimutal(k) + PI*4.0_RP/3.0_RP
    end do

    ! average_conditions are thrust and rotor force
    if (MPI_Process % isRoot) then
      if (self % save_average) then
          do k = 1, self%num_turbines
              allocate ( self%turbine_t(k)%average_conditions(self%turbine_t(k)%num_blade_sections,5) )
          end do
          do k = 1, self%num_turbines
              self % turbine_t(k) % average_conditions(:,:) = 0.0_RP
          end do
          self % number_iterations = 0
      end if
    end if
!
!   Get the solution file name
!   --------------------------
    solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = LINE_LENGTH )
    solution_file = trim(getFileName(solution_file))
    self % file_name = trim(solution_file)
!
!   Get the elements that are relevant
!   ----------------------------------
    allocate(numElementsPerTurbine(self%num_turbines))
    numElementsPerTurbine = 0
    element_loop:do eID = 1, mesh%no_of_elements
        do k=1, self%num_turbines
            tolerance = self%tolerance_factor*self%turbine_t(k)%radius
            r_square = POW2(minval(mesh%elements(eID)%geom%x(2,:,:,:))-self%turbine_t(k)%hub_cood_y) + POW2(minval(mesh%elements(eID)%geom%x(3,:,:,:))-self%turbine_t(k)%hub_cood_z)
            if( r_square <= POW2(self%turbine_t(k)%radius+tolerance) &
                .and. minval(mesh%elements(eID)%geom%x(1,:,:,:)) < self%turbine_t(k)%hub_cood_x+tolerance &
                .and. maxval(mesh%elements(eID)%geom%x(1,:,:,:)) >self%turbine_t(k)%hub_cood_x-tolerance) then
                numElementsPerTurbine(k) = numElementsPerTurbine(k) + 1
                cycle element_loop
            end if
        end do
    end do element_loop
    nelem = sum(numElementsPerTurbine)
    allocate(elementsActuated(nelem),turbineOfElement(nelem))
    eIndex = 0
    element_loop2:do eID = 1, mesh%no_of_elements
        do k=1, self%num_turbines
            tolerance = self%tolerance_factor*self%turbine_t(k)%radius
            r_square = POW2(minval(mesh%elements(eID)%geom%x(2,:,:,:))-self%turbine_t(k)%hub_cood_y) + POW2(minval(mesh%elements(eID)%geom%x(3,:,:,:))-self%turbine_t(k)%hub_cood_z)
            if( r_square <= POW2(self%turbine_t(k)%radius+tolerance) &
                .and. minval(mesh%elements(eID)%geom%x(1,:,:,:)) < self%turbine_t(k)%hub_cood_x+tolerance &
                .and. maxval(mesh%elements(eID)%geom%x(1,:,:,:)) >self%turbine_t(k)%hub_cood_x-tolerance) then
                eIndex = eIndex + 1
                elementsActuated(eIndex) = eID
                turbineOfElement(eIndex) = k
                cycle element_loop2
            end if
        end do
    end do element_loop2

    ! precalculate epsilon for projection mode
    if (self % calculate_with_projection) then
        if (MPI_Process % doMPIAction) then
            if (nelem .gt. 0) then
              delta_temp = (mesh % elements(elementsActuated(1)) % geom % Volume / product(mesh % elements(elementsActuated(1)) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
              delta_count = 1
            else
              delta_temp = 0.0_RP
              delta_count = 0
            end if 
#ifdef _HAS_MPI_
            call mpi_allreduce(delta_temp, delta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(delta_count, delta_paritions, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
            delta = delta / real(delta_paritions,kind=RP)
        else    
            delta = (mesh % elements(elementsActuated(1)) % geom % Volume / product(mesh % elements(elementsActuated(1)) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
        end if
!$omp do     schedule(runtime)private(i,j,k)
        do k = 1, self%num_turbines
          do j = 1, self%turbine_t(k)%num_blades
             do i = 1, self%turbine_t(k)%num_blade_sections
               self % turbine_t(k) % blade_t(j) % gauss_epsil_delta(i) = delta
             end do
          enddo
        enddo
!$omp end do
    end if

!
!   Create output files
!   -------------------
    if (MPI_Process % isRoot) then
      do k=1, self%num_turbines
        write(file_id, '(I3.3)') k

        write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Line_Forces_turb_", trim(file_id) , ".dat"
        open ( newunit = fID , file = trim(arg) , status = "unknown" ,    action = "write" ) 
        write(fid,'(10(2X,A24))') "time", "thrust_1", "blade_torque_1", "blade_root_bending_1", "thrust_2", "blade_torque_2", "blade_root_bending_2", "thrust_3", "blade_torque_3", "blade_root_bending_3"
        close(fid)
!
        write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Line_CP_CT_turb_", trim(file_id) , ".dat"
        open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" ) 
        write(fid,'(3(2X,A24))') "time", "Cp(power_coef)", "Ct(thust_coef)"
        close(fid)
!
        if (self % save_average) then
            write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Line_average_turb_", trim(file_id) , ".dat"
            open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" ) 
            write(fid,'(6(2X,A24))') "R", "U", "AoA", "Re", "Tangential_Force", "Axial_Force"
            close(fid)
        end if
      end do
    end if
!
       self % active = .true.
!
   end subroutine ConstructFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
   subroutine DestructFarm(self)
   implicit none
   class(Farm_t), intent(inout)       :: self

   ! integer                            :: i, j, k

   ! do i=1, self%num_turbines
   !  do j=1, self%turbine_t(j)%num_blades
   !      do k=1, self%turbine_t(i)%num_blade_sections
   !          deallocate ( self%turbine_t(i)%blade_t(j)%airfoil_t(k)%aoa, &
   !                  self%turbine_t(i)%blade_t(j)%airfoil_t(k)%Re, &
   !                  self%turbine_t(i)%blade_t(j)%airfoil_t(k)%cl, &
   !                  self%turbine_t(i)%blade_t(j)%airfoil_t(k)%cd )
   !      end do
   !      deallocate( self%turbine_t(i)%blade_t(j)%r_R, &
   !                  self%turbine_t(i)%blade_t(j)%chord, &
   !                  self%turbine_t(i)%blade_t(j)%twist, &
   !                  self%turbine_t(i)%blade_t(j)%local_velocity, &
   !                  self%turbine_t(i)%blade_t(j)%local_angle, &
   !                  self%turbine_t(i)%blade_t(j)%point_xyz_loc, &
   !                  self%turbine_t(i)%blade_t(j)%local_torque, &
   !                  self%turbine_t(i)%blade_t(j)%local_thrust, &
   !                  self%turbine_t(i)%blade_t(j)%local_rotor_force, &
   !                  self%turbine_t(i)%blade_t(j)%local_root_bending, &
   !                  self%turbine_t(i)%blade_t(j)%local_Re, &
   !                  self%turbine_t(i)%blade_t(j)%num_airfoils, &
   !                  self%turbine_t(i)%blade_t(j)%airfoil_files, &
   !                  self%turbine_t(i)%blade_t(j)%airfoil_t, &
   !                  self%turbine_t(i)%blade_t(j)%local_gaussian_sum, &
   !                  self%turbine_t(i)%blade_t(j)%local_thrust_temp )

   !  end do
   !  safedeallocate(self%turbine_t(i)%average_conditions)
   ! end do

   deallocate(self%turbine_t)

   end subroutine DestructFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!

   subroutine UpdateFarm(self,time, mesh)
   use fluiddata
   use HexMeshClass
   use PhysicsStorage
   use MPI_Process_Info
   implicit none

   class(Farm_t), intent(inout)      :: self
   real(kind=RP), intent(in)         :: time
   type(HexMesh), intent(in)         :: mesh

   !local variables
   integer                           :: ii, jj, i, j, k, kk
   real(kind=RP)                     :: dt, interp, delta_temp
   logical                           :: found
   integer                           :: eID, ierr
   real(kind=RP), dimension(NDIM)    :: x, xi
   real(kind=RP), dimension(NCONS)   :: Q, Qtemp

   if (.not. self % active) return

   dt = time - self % time
   self % time = time

   dt = dt * Lref / refValues%V
   ! only for constant rot_speed
   ! theta = self%turbine_t(1)%rot_speed * t
   interp = 1.0_RP

   projection_cond:if (self%calculate_with_projection) then
!
!    ----------------------------------------------------------------------------------
!    calculate for all mesh points its contribution based on the gaussian interpolation
!    ----------------------------------------------------------------------------------
!
!$omp do schedule(runtime)private(ii,jj,kk)
   do kk = 1, self%num_turbines
      do jj = 1, self%turbine_t(kk)%num_blades

         self%turbine_t(kk)%blade_t(jj)%azimuth_angle = self%turbine_t(kk)%blade_t(jj)%azimuth_angle + self%turbine_t(kk)%rot_speed*dt

         self%turbine_t(kk)%blade_t(jj)%local_rotor_force(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_thrust(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_torque(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_root_bending(:) = 0.0_RP

         self%turbine_t(kk)%blade_t(jj)%local_rotor_force_temp(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_thrust_temp(:)=0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_velocity_temp(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_angle_temp(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_Re_temp(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(:)= 0.0_RP
      
         do ii = 1, self%turbine_t(kk)%num_blade_sections
           ! y,z coordinate of every acutator line point
           self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,2) = self%turbine_t(kk)%hub_cood_y + self%turbine_t(kk)%blade_t(jj)%r_R(ii) * cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle)
           self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,3) = self%turbine_t(kk)%hub_cood_z + self%turbine_t(kk)%blade_t(jj)%r_R(ii) * sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle)
      
         end do
      enddo
    enddo
!$omp end do
!

! no projection
   else projection_cond
!
!    ----------------------------------------------------------------
!    use the local Q based on the position of the actuator line point
!    ----------------------------------------------------------------
!
!$omp do schedule(runtime)private(ii,jj,kk,eID,Q,Qtemp,delta_temp,xi,found)
    do kk = 1, self%num_turbines
      do jj = 1, self%turbine_t(kk)%num_blades

         self%turbine_t(kk)%blade_t(jj)%azimuth_angle = self%turbine_t(kk)%blade_t(jj)%azimuth_angle + self%turbine_t(kk)%rot_speed*dt

         self%turbine_t(kk)%blade_t(jj)%local_rotor_force(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_thrust(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_torque(:) = 0.0_RP
         self%turbine_t(kk)%blade_t(jj)%local_root_bending(:) = 0.0_RP
!
         do ii = 1, self%turbine_t(kk)%num_blade_sections
           ! y,z coordinate of every acutator line point
           self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,2) = self%turbine_t(kk)%hub_cood_y + self%turbine_t(kk)%blade_t(jj)%r_R(ii) * cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle)
           self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,3) = self%turbine_t(kk)%hub_cood_z + self%turbine_t(kk)%blade_t(jj)%r_R(ii) * sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle)
!
!          -----------------------------------
!          get the elements of each line point
!          -----------------------------------
!
           x = [self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,1),self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,2),self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,3)]
           call FindActuatorPointElement(mesh, x, eID, xi, found)
           if (found) then
             ! averaged state values of the cell
             Qtemp = element_averageQ(mesh,eID, xi, self % averageSubElement)
             delta_temp = (mesh % elements(eID) % geom % Volume / product(mesh % elements(eID) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
           else
             Qtemp = 0.0_RP
             delta_temp = 0.0_RP
           end if
           if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
             call mpi_allreduce(Qtemp, Q, NCONS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
             call mpi_allreduce(delta_temp, self%turbine_t(kk)%blade_t(jj)%gauss_epsil_delta(ii), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
           else
               Q = Qtemp
               self % turbine_t(kk) % blade_t(jj) % gauss_epsil_delta(ii) = delta_temp

           end if
           if (all(Q .eq. 0.0_RP)) then
             print*, "Actuator line point not found in mesh, x: ", x
             call exit(99)
           end if
           call FarmGetLocalForces(self, ii, jj, kk, Q, interp, self%turbine_t(kk)%blade_t(jj)%local_angle(ii), self%turbine_t(kk)%blade_t(jj)%local_velocity(ii), & 
                                   self%turbine_t(kk)%blade_t(jj)%local_Re(ii), self%turbine_t(kk)%blade_t(jj)%local_thrust(ii), self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii))
         end do
      enddo
    enddo

!$omp end do

   end if projection_cond
!
   end subroutine UpdateFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ForcesFarm(self, mesh, time)
   use PhysicsStorage
   use HexMeshClass
   use fluiddata
   implicit none

   class(Farm_t) , intent(inout)     :: self
   type(HexMesh), intent(in)         :: mesh
   real(kind=RP),intent(in)          :: time

! local vars
   real(kind=RP)                     :: Non_dimensional, t, interp
   integer                           :: ii,jj, kk
   integer                           :: i,j, k
   integer                           :: eID, eIndex
   real(kind=RP), dimension(NDIM)    :: actuator_source
   real(kind=RP)                     :: local_angle
   real(kind=RP)                     :: local_velocity
   real(kind=RP)                     :: local_Re
   real(kind=RP)                     :: local_thrust
   real(kind=RP)                     :: local_rotor_force
#if defined(MULTIPHASE)
   real(kind=RP)                     :: invSqrtRho
#endif

    if (.not. self % active) return

    Non_dimensional = POW2(refValues % V) * refValues % rho / Lref
    t = time * Lref / refValues % V

    if (self % calculate_with_projection) then

        do eIndex = 1, size(elementsActuated)
            eID = elementsActuated(eIndex)
            kk = turbineOfElement(eIndex)

            do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                actuator_source(:) = 0.0_RP

                do jj = 1, self % turbine_t(kk) % num_blades
                    do ii = 1, self % turbine_t(kk) % num_blade_sections
                        interp = GaussianInterpolation(self%epsilon_type, mesh % elements(eID) % geom % x(:,i,j,k), self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,:), &
                                                       self%turbine_t(kk)%blade_t(jj)%chord(ii), self%gauss_epsil,self%turbine_t(kk)%blade_t(jj)%gauss_epsil_delta(ii))
                        call FarmGetLocalForces(self, ii, jj, kk, mesh%elements(eID)%storage%Q(:,i,j,k), interp, local_angle, local_velocity, local_Re, local_thrust, local_rotor_force)

                        ! minus account action-reaction effect, is the force on the fliud
                        actuator_source(1) = actuator_source(1) - local_thrust
                        actuator_source(2) = actuator_source(2) - (-local_rotor_force*sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) )
                        actuator_source(3) = actuator_source(3) - local_rotor_force*cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) 

                        !acumulate in temporal variables, for each time step as the non temp are recalculated for each element
                        self%turbine_t(kk)%blade_t(jj)%local_thrust_temp(ii)=self%turbine_t(kk)%blade_t(jj)%local_thrust_temp(ii)+local_thrust
                        self%turbine_t(kk)%blade_t(jj)%local_rotor_force_temp(ii)=self%turbine_t(kk)%blade_t(jj)%local_rotor_force_temp(ii)+local_rotor_force
                        self%turbine_t(kk)%blade_t(jj)%local_velocity_temp(ii)=self%turbine_t(kk)%blade_t(jj)%local_velocity_temp(ii)+local_velocity*interp
                        self%turbine_t(kk)%blade_t(jj)%local_angle_temp(ii)=self%turbine_t(kk)%blade_t(jj)%local_angle_temp(ii)+local_angle*interp
                        self%turbine_t(kk)%blade_t(jj)%local_Re_temp(ii)=self%turbine_t(kk)%blade_t(jj)%local_Re_temp(ii)+local_Re*interp

                        self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)=self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)+interp

                    enddo
                enddo
                ! NS(IRHOU:IRHOW) = NS(IRHOU:IRHOW) + actuator_source(:) / Non_dimensional
                actuator_source = actuator_source/ Non_dimensional

#if defined(NAVIERSTOKES)
                mesh % elements(eID) % storage % S_NS(IRHOU,i,j,k) = actuator_source(1)
                mesh % elements(eID) % storage % S_NS(IRHOV,i,j,k) = actuator_source(2)
                mesh % elements(eID) % storage % S_NS(IRHOW,i,j,k) = actuator_source(3)
#elif defined(INCNS)
                mesh % elements(eID) % storage % S_NS(INSRHOU,i,j,k) = actuator_source(1)
                mesh % elements(eID) % storage % S_NS(INSRHOV,i,j,k) = actuator_source(2)
                mesh % elements(eID) % storage % S_NS(INSRHOW,i,j,k) = actuator_source(3)
#elif defined(MULTIPHASE)
                invSqrtRho = 1.0_RP / sqrt(mesh % elements(eID) % storage % rho(i,j,k))
                mesh % elements(eID) % storage % S_NS(IMSQRHOU,i,j,k) = actuator_source(1)*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOV,i,j,k) = actuator_source(2)*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOW,i,j,k) = actuator_source(3)*invSqrtRho
#endif
            end do                  ; end do                ; end do
        end do
    
    else ! no projection

!$omp do schedule(runtime) private(i,j,k,ii,jj,kk,actuator_source,eID,interp)
        do eIndex = 1, size(elementsActuated)
            eID = elementsActuated(eIndex)
            ! only one turbine is associated for one element
            kk = turbineOfElement(eIndex)

            do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                actuator_source(:) = 0.0_RP

                do jj = 1, self % turbine_t(kk) % num_blades
                    do ii = 1, self % turbine_t(kk) % num_blade_sections
                        interp = GaussianInterpolation(self%epsilon_type, mesh % elements(eID) % geom % x(:,i,j,k), self%turbine_t(kk)%blade_t(jj)%point_xyz_loc(ii,:), &
                                                       self%turbine_t(kk)%blade_t(jj)%chord(ii), self%gauss_epsil,self%turbine_t(kk)%blade_t(jj)%gauss_epsil_delta(ii))
        
                        ! minus account action-reaction effect, is the force on the fluid
                        actuator_source(1) = actuator_source(1) - self%turbine_t(kk)%blade_t(jj)%local_thrust(ii) * interp
                        actuator_source(2) = actuator_source(2) - (-self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii)*sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle)  * interp)
                        actuator_source(3) = actuator_source(3) - self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii)*cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) * interp 

                    enddo
                enddo
                ! NS(IRHOU:IRHOW) = NS(IRHOU:IRHOW) + actuator_source(:) / Non_dimensional
                actuator_source = actuator_source/ Non_dimensional
#if defined(NAVIERSTOKES)
                mesh % elements(eID) % storage % S_NS(IRHOU,i,j,k) = actuator_source(1)
                mesh % elements(eID) % storage % S_NS(IRHOV,i,j,k) = actuator_source(2)
                mesh % elements(eID) % storage % S_NS(IRHOW,i,j,k) = actuator_source(3)
#elif defined(INCNS)
                mesh % elements(eID) % storage % S_NS(INSRHOU,i,j,k) = actuator_source(1)
                mesh % elements(eID) % storage % S_NS(INSRHOV,i,j,k) = actuator_source(2)
                mesh % elements(eID) % storage % S_NS(INSRHOW,i,j,k) = actuator_source(3)
#elif defined(MULTIPHASE)
                invSqrtRho = 1.0_RP / sqrt(mesh % elements(eID) % storage % rho(i,j,k))
                mesh % elements(eID) % storage % S_NS(IMSQRHOU,i,j,k) = actuator_source(1)*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOV,i,j,k) = actuator_source(2)*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOW,i,j,k) = actuator_source(3)*invSqrtRho
#endif
            end do                  ; end do                ; end do
        end do
!$omp end do

    endif

   end subroutine  ForcesFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine WriteFarmForces(self,time,iter,last)
   use fluiddata
   use PhysicsStorage
   use MPI_Process_Info
   implicit none

   class(Farm_t), intent(inout)  :: self
   real(kind=RP),intent(in)      :: time
   integer, intent(in)           :: iter
   logical, optional             :: last
   integer                       :: fid
   CHARACTER(LEN=LINE_LENGTH)    :: arg
   real(kind=RP)                 :: t
   integer                       :: ii, jj, kk
   logical                       :: isLast
   logical                       :: save_instant
   integer                       :: ierr
   real(kind=RP), dimension(:), allocatable :: local_thrust_temp, local_rotor_force_temp, local_gaussian_sum, local_velocity_temp, local_angle_temp, local_Re_temp
   CHARACTER(LEN=5)           :: file_id

   if (.not. self % active) return

   if (present(last)) then
       isLast = last
   else
       isLast = .false.
   end if

   ! for last iteration save average values only, all calculations and comunications were done before
   if (isLast) then
      if ( .not. self % save_average ) return
      if ( .not. MPI_Process % isRoot ) return
      do kk=1, self%num_turbines
          write(file_id, '(I3.3)') kk
          write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Line_average_turb_", trim(file_id) , ".dat"
          open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
          do ii = 1, self % turbine_t(kk) % num_blade_sections
            write(fid,"(6(2X,ES24.16))") self%turbine_t(kk)%blade_t(1)%r_R(ii), self%turbine_t(kk)%average_conditions(ii,:)
          end do
          close(fid)
      end do
      return
    end if

   save_instant = self%save_instant .and. ( mod(iter,self % save_iterations) .eq. 0 )
   t = time * Lref / refValues%V

   if (self%calculate_with_projection) then
     ! this is necessary for Gaussian weighted sum
     
    do kk = 1, self%num_turbines
         do jj = 1, self%turbine_t(kk)%num_blades
             self%turbine_t(kk)%blade_t(jj)%local_thrust(:) = 0.0_RP
             self%turbine_t(kk)%blade_t(jj)%local_rotor_force(:) = 0.0_RP
             self%turbine_t(kk)%blade_t(jj)%local_velocity(:) = 0.0_RP
             self%turbine_t(kk)%blade_t(jj)%local_angle(:) = 0.0_RP
             self%turbine_t(kk)%blade_t(jj)%local_Re(:) = 0.0_RP
         end do
    end do
  
     if ( (MPI_Process % doMPIAction) ) then
       do kk = 1, self%num_turbines
         associate (num_blade_sections => self%turbine_t(kk)%num_blade_sections)
           allocate( local_thrust_temp(num_blade_sections), local_rotor_force_temp(num_blade_sections), local_gaussian_sum(num_blade_sections),&
                     local_velocity_temp(num_blade_sections), local_angle_temp(num_blade_sections), local_Re_temp(num_blade_sections) )

           do jj = 1, self%turbine_t(kk)%num_blades
             local_thrust_temp = self%turbine_t(kk)%blade_t(jj)%local_thrust_temp
             local_rotor_force_temp = self%turbine_t(kk)%blade_t(jj)%local_rotor_force_temp
             local_velocity_temp = self%turbine_t(kk)%blade_t(jj)%local_velocity_temp
             local_angle_temp = self%turbine_t(kk)%blade_t(jj)%local_angle_temp
             local_Re_temp = self%turbine_t(kk)%blade_t(jj)%local_Re_temp
             local_gaussian_sum = self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum
  
#ifdef _HAS_MPI_
                call mpi_allreduce(local_thrust_temp, self%turbine_t(kk)%blade_t(jj)%local_thrust_temp, num_blade_sections, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call mpi_allreduce(local_rotor_force_temp, self%turbine_t(kk)%blade_t(jj)%local_rotor_force_temp, num_blade_sections, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call mpi_allreduce(local_velocity_temp, self%turbine_t(kk)%blade_t(jj)%local_velocity_temp, num_blade_sections, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call mpi_allreduce(local_angle_temp, self%turbine_t(kk)%blade_t(jj)%local_angle_temp, num_blade_sections, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call mpi_allreduce(local_Re_temp, self%turbine_t(kk)%blade_t(jj)%local_Re_temp, num_blade_sections, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call mpi_allreduce(local_gaussian_sum, self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum, num_blade_sections, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  
           end do
           deallocate(local_thrust_temp,local_rotor_force_temp,local_gaussian_sum,local_velocity_temp,local_angle_temp,local_Re_temp)
         endassociate
       end do
     end if
         
!$omp do schedule(runtime)private(ii,jj,kk)
     do kk = 1, self%num_turbines
        do jj = 1, self%turbine_t(kk)%num_blades
             do ii = 1, self%turbine_t(kk)%num_blade_sections
     
                 self%turbine_t(kk)%blade_t(jj)%local_thrust(ii)=self%turbine_t(kk)%blade_t(jj)%local_thrust(ii)+self%turbine_t(kk)%blade_t(jj)%local_thrust_temp(ii)/self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)
                 self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii)=self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii)+self%turbine_t(kk)%blade_t(jj)%local_rotor_force_temp(ii)/self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)
                 self%turbine_t(kk)%blade_t(jj)%local_angle(ii)=self%turbine_t(kk)%blade_t(jj)%local_angle(ii)+self%turbine_t(kk)%blade_t(jj)%local_angle_temp(ii)/self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)
                 self%turbine_t(kk)%blade_t(jj)%local_velocity(ii)=self%turbine_t(kk)%blade_t(jj)%local_velocity(ii)+self%turbine_t(kk)%blade_t(jj)%local_velocity_temp(ii)/self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)
                 self%turbine_t(kk)%blade_t(jj)%local_Re(ii)=self%turbine_t(kk)%blade_t(jj)%local_Re(ii)+self%turbine_t(kk)%blade_t(jj)%local_Re_temp(ii)/self%turbine_t(kk)%blade_t(jj)%local_gaussian_sum(ii)
     
             enddo
         enddo
       enddo
!$omp end do
     
   end if

   if ( .not. MPI_Process % isRoot ) return

   ! save in memory the time step forces for each element blade and the whole blades
   call FarmUpdateBladeForces(self)

!write output torque thrust to file
      do kk=1, self%num_turbines
        write(file_id, '(I3.3)') kk
        write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Line_Forces_turb_", trim(file_id) , ".dat"

        open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
        write(fid,"(10(2X,ES24.16))") t, &
        self%turbine_t(kk)%blade_thrust(1),self%turbine_t(kk)%blade_torque(1),self%turbine_t(kk)%blade_root_bending(1), &
        self%turbine_t(kk)%blade_thrust(2),self%turbine_t(kk)%blade_torque(2),self%turbine_t(kk)%blade_root_bending(2), &
        self%turbine_t(kk)%blade_thrust(3),self%turbine_t(kk)%blade_torque(3),self%turbine_t(kk)%blade_root_bending(3)
        close(fid)

        write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Line_CP_CT_turb_", trim(file_id) , ".dat"
        open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
        write(fid,"(10(2X,ES24.16))") t, self%turbine_t(kk)%Cp, self%turbine_t(kk)%Ct
        close(fid)

        if (save_instant) then
          do jj = 1, self%turbine_t(kk)%num_blades
            write(arg , '(2A,I3.3,A,I10.10,3A)') trim(self%file_name), "_Actuator_Line_instant_",jj ,"_" ,iter, "_turb_", trim(file_id), ".dat"
            open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" ) 
            write(fid,'(6(2X,A24))') "R", "U", "AoA", "Re", "Tangential_Force", "Axial_Force"
            do ii = 1, self % turbine_t(kk) % num_blade_sections
              write(fid,"(6(2X,ES24.16))") self%turbine_t(kk)%blade_t(1)%r_R(ii), &
                   self%turbine_t(kk)%blade_t(jj)%local_velocity(ii), &
                   ( self%turbine_t(kk)%blade_t(jj)%local_angle(ii) - (self%turbine_t(kk)%blade_t(jj)%twist(ii) + self%turbine_t(kk)%blade_pitch) ) * 180.0_RP / PI, &
                   self%turbine_t(kk)%blade_t(jj)%local_Re(ii), &
                   self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii), &
                   self%turbine_t(kk)%blade_t(jj)%local_thrust(ii)
            end do
            close(fid)
          end do
        end if 
    end do
  !endif

end subroutine WriteFarmForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Subroutine FarmGetLocalForces(self, ii, jj, kk, Q, interp, local_angle, local_velocity, local_Re, local_thrust, local_rotor_force)
        use PhysicsStorage
        use fluiddata
#if defined(NAVIERSTOKES)
        use VariableConversion, only: Temperature, SutherlandsLaw
#endif
        implicit none
        class(Farm_t)                                 :: self
        integer, intent(in)                           :: ii, jj, kk
        real(kind=RP), dimension(NCONS), intent(in)   :: Q
        real(kind=RP), intent(in)                     :: interp
        real(kind=RP), intent(out)                    :: local_angle
        real(kind=RP), intent(out)                    :: local_velocity
        real(kind=RP), intent(out)                    :: local_Re
        real(kind=RP), intent(out)                    :: local_thrust
        real(kind=RP), intent(out)                    :: local_rotor_force

        !local variables
        real(kind=RP)                                 :: density, Cl, Cd, aoa, g1_func, tip_correct
        ! real(kind=RP)                                 :: angle_temp
        real(kind=RP)                                 :: wind_speed_axial, wind_speed_rot
        real(kind=RP)                                 :: lift_force, drag_force
        real(kind=RP)                                 :: T, muL
#if defined(MULTIPHASE)
        real(kind=RP)                                 :: rho, invSqrtRho
#endif
!
!       -----------------------------
!       get airfoil related variables
!       -----------------------------
!
#if defined(NAVIERSTOKES)
        density = Q(IRHO) * refValues % rho

        ! project [v.w] in the rotational direction (theta in cylindrical coordinates)
        wind_speed_axial = (Q(IRHOU)/Q(IRHO)) * refValues % V ! our x is the z in cylindrical
        wind_speed_rot = ( -Q(IRHOV)*sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) + Q(IRHOW)*cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) ) / Q(IRHO) * refValues % V

#elif defined(INCNS)
        density = Q(INSRHO) * refValues % rho

        ! project [v.w] in the rotational direction (theta in cylindrical coordinates)
        wind_speed_axial = (Q(INSRHOU)/Q(INSRHO)) * refValues % V ! our x is the z in cylindrical
        wind_speed_rot = ( -Q(INSRHOV)*sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) + Q(INSRHOW)*cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) ) / Q(INSRHO) * refValues % V
#elif defined(MULTIPHASE)
        rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2)) * Q(IMC)
        rho = min(max(rho, dimensionless % rho_min),dimensionless % rho_max)
        invSqrtRho = 1.0_RP/sqrt(rho)
        density = rho * refValues % rho

        ! project [v.w] in the rotational direction (theta in cylindrical coordinates)
        wind_speed_axial = (Q(IMSQRHOU)*invSqrtRho) * refValues % V ! our x is the z in cylindrical
        wind_speed_rot = ( -Q(IMSQRHOV)*sin(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) + Q(IMSQRHOW)*cos(self%turbine_t(kk)%blade_t(jj)%azimuth_angle) )*invSqrtRho * refValues % V

#endif

        tip_correct = 1.0_RP
        aoa = 0.0_RP

#if defined(NAVIERSTOKES)
         T     = Temperature(Q)
         muL = SutherlandsLaw(T) * refValues % mu
#elif defined(INCNS)
         muL = refValues % mu
#elif defined(MULTIPHASE)
         muL = refValues % mu
#endif

        local_velocity = sqrt( POW2(self % turbine_t(kk) % rot_speed*self % turbine_t(kk) % blade_t(jj) % r_R(ii) - wind_speed_rot) + POW2(wind_speed_axial) )
        local_angle = atan( wind_speed_axial / (self % turbine_t(kk) % rot_speed*self % turbine_t(kk) % blade_t(jj) % r_R(ii) - wind_speed_rot) ) 

        ! alpha = phi - gamma, gamma = blade pitch + airfoil local twist
        aoa = local_angle - (self % turbine_t(kk) % blade_t(jj) % twist(ii) + self % turbine_t(kk) % blade_pitch) 
        local_Re = local_velocity * self % turbine_t(kk) % blade_t(jj) % chord(ii) * density / muL
        call Get_Cl_Cl_from_airfoil_data(self % turbine_t(kk) % blade_t(jj) % airfoil_t(ii), aoa, local_Re, Cl, Cd)
!
!       ---------------
!       tip correction
!       ---------------
!
        g1_func=exp(-self % turbine_t(kk) % blade_t(jj) % tip_c1*(self % turbine_t(kk) % num_blades*self % turbine_t(kk) % rot_speed*self % turbine_t(kk) % radius/refValues % V-self % turbine_t(kk) % blade_t(jj) % tip_c2))+0.1
        ! in this it should be the local_angle at the tip

        ! without perturbation
        ! angle_temp = atan(refValues%V/(self%turbine_t(1)%rot_speed*self%turbine_t(1)%radius))
        ! only axial wind speed
        ! angle_temp = atan(wind_speed_axial/(self%turbine_t(1)%rot_speed*self%turbine_t(1)%radius))

        ! 1 possibility: use the global radius and angle (i.e. the tip values)
        ! tip_correct = 2.0_RP/PI*(acos( exp(-g1_func*self%turbine_t(1)%num_blades*(self%turbine_t(1)%radius-self%turbine_t(1)%blade_t(jj)%r_R(ii)) / &
                      ! (2.0_RP*self%turbine_t(1)%radius*sin(angle_temp))) ))

        ! 2 possibility: use the local radius and angle
        tip_correct = 2.0_RP/PI*(acos( exp(-g1_func*self%turbine_t(kk)%num_blades*(self%turbine_t(kk)%radius-self%turbine_t(kk)%blade_t(jj)%r_R(ii)) / abs(2.0_RP*self%turbine_t(kk)%blade_t(jj)%r_R(ii)*sin(local_angle))) ))
!
!       --------------------------------
!       Save forces on the blade segment
!       --------------------------------
!
         ! lift=Cl*1/2.rho*v_local^2*Surface ; and surface=section area of the blade S=length*chord
         ! lift and drag are multiply by the gaussian interp for the case of each mesh node contributing to the force (for local only is 1)

         lift_force = 0.5_RP * density * Cl * tip_correct * POW2(local_velocity) * self%turbine_t(kk)%blade_t(jj)%chord(ii) &
                      * (self%turbine_t(kk)%blade_t(jj)%r_R(ii) - self%turbine_t(kk)%blade_t(jj)%r_R(ii-1)) * interp

        drag_force = 0.5_RP * density * Cd * tip_correct * POW2(local_velocity) * self%turbine_t(kk)%blade_t(jj)%chord(ii) &
                      * (self%turbine_t(kk)%blade_t(jj)%r_R(ii) - self%turbine_t(kk)%blade_t(jj)%r_R(ii-1)) * interp

        local_rotor_force = lift_force * sin(local_angle) - drag_force * cos(local_angle)
                                  
        local_thrust = lift_force * cos(local_angle) + drag_force * sin(local_angle)
                                                     !

    End Subroutine FarmGetLocalForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Subroutine FarmUpdateBladeForces(self)
        use fluiddata
        Implicit None

        class(Farm_t), intent(inout)      :: self
        !local variables
        integer                           :: ii, jj, kk
        real(kind=RP), dimension(:), allocatable  :: aoa
!
!       ------------------------------
!       Save forces on the whole blade
!       ------------------------------
!
!$omp do schedule(runtime)private(ii,jj,kk)
    do kk = 1, self%num_turbines
      do jj = 1, self%turbine_t(kk)%num_blades
          self%turbine_t(kk)%blade_thrust(jj) = 0.0_RP
          self%turbine_t(kk)%blade_torque(jj) = 0.0_RP
          self%turbine_t(kk)%blade_root_bending(jj) = 0.0_RP

          do ii = 1, self%turbine_t(kk)%num_blade_sections

              self%turbine_t(kk)%blade_t(jj)%local_torque(ii) = self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii)*self%turbine_t(kk)%blade_t(jj)%r_R(ii)

              self%turbine_t(kk)%blade_t(jj)%local_root_bending(ii) = sqrt(POW2(self%turbine_t(kk)%blade_t(jj)%local_thrust(ii)) + &
                  POW2(self%turbine_t(kk)%blade_t(jj)%local_rotor_force(ii))) * self%turbine_t(kk)%blade_t(jj)%r_R(ii)

              self%turbine_t(kk)%blade_thrust(jj)=self%turbine_t(kk)%blade_thrust(jj)+self%turbine_t(kk)%blade_t(jj)%local_thrust(ii)
              self%turbine_t(kk)%blade_torque(jj)=self%turbine_t(kk)%blade_torque(jj)+self%turbine_t(kk)%blade_t(jj)%local_torque(ii)
              self%turbine_t(kk)%blade_root_bending(jj)=self%turbine_t(kk)%blade_root_bending(jj)+self%turbine_t(kk)%blade_t(jj)%local_root_bending(ii)
         enddo
      enddo
    enddo
!$omp end do

    do kk = 1, self%num_turbines

      self%turbine_t(kk)%Cp = 2.0_RP * (self%turbine_t(kk)%blade_torque(1)+self%turbine_t(kk)%blade_torque(2)+self%turbine_t(kk)%blade_torque(3)) * self%turbine_t(kk)%rot_speed / &
                              (refValues%rho * POW3(refValues%V) * pi * POW2(self%turbine_t(kk)%radius))

      self%turbine_t(kk)%Ct = 2.0_RP * (self%turbine_t(kk)%blade_thrust(1)+self%turbine_t(kk)%blade_thrust(2)+self%turbine_t(kk)%blade_thrust(3)) / &
                              (refValues%rho * POW2(refValues%V) * pi * POW2(self%turbine_t(kk)%radius))

    enddo
!
!        ----------------------
!        Save average variables
!        ----------------------
!
    if (self % save_average) then
      do kk = 1, self%num_turbines
         self % turbine_t(kk) % average_conditions = self % turbine_t(kk) % average_conditions * real(self % number_iterations,RP)

         !saving only for blade 1
         jj =1
         allocate(aoa(self%turbine_t(kk)%num_blade_sections))
         aoa = self%turbine_t(kk)%blade_t(jj)%local_angle(:) - (self%turbine_t(kk)%blade_t(jj)%twist(:) + self%turbine_t(kk)%blade_pitch) 
         self % turbine_t(kk) % average_conditions(:,1) = self % turbine_t(kk) % average_conditions(:,1) + self%turbine_t(kk)%blade_t(jj)%local_velocity
         self % turbine_t(kk) % average_conditions(:,2) = self % turbine_t(kk) % average_conditions(:,2) + aoa * 180.0_RP / PI
         self % turbine_t(kk) % average_conditions(:,3) = self % turbine_t(kk) % average_conditions(:,3) + self%turbine_t(kk)%blade_t(jj)%local_Re
         self % turbine_t(kk) % average_conditions(:,4) = self % turbine_t(kk) % average_conditions(:,4) + self%turbine_t(kk)%blade_t(jj)%local_rotor_force
         self % turbine_t(kk) % average_conditions(:,5) = self % turbine_t(kk) % average_conditions(:,5) + self%turbine_t(kk)%blade_t(jj)%local_thrust
         deallocate(aoa)
      end do

         self % number_iterations = self % number_iterations + 1
      do kk = 1, self%num_turbines
         self % turbine_t(kk) % average_conditions = self % turbine_t(kk) % average_conditions / real(self % number_iterations,RP)
      end do
    end if 
!
    End Subroutine FarmUpdateBladeForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Function GaussianInterpolation(epsilon_type, x, x_point, chord, gauss_epsil, delta, Cd)
        implicit none
        integer, intent(in)                     :: epsilon_type
        real(kind=RP), intent(in)               :: x(NDIM)
        real(kind=RP), intent(in)               :: x_point(NDIM)
        real(kind=RP), intent(in)               :: chord
        real(kind=RP), intent(in)               :: gauss_epsil
        real(kind=RP), intent(in)               :: delta
        real(kind=RP), intent(in), optional     :: Cd
        real(kind=RP)                           :: GaussianInterpolation

        !local variables
        real(kind=RP)                           :: epsil

        select case (epsilon_type)
        case (0)
! EPSILON - option 1 (from file)
            epsil = gauss_epsil
        case (1)
! EPSILON - option 2
            if (present(Cd)) then
                epsil = max(chord/4.0_RP,chord*Cd/2.0_RP)
            else
                epsil = gauss_epsil
            end if
        case (2)
! EPSILON - option 3 (k is from file)
! eps = k*delta; k is in gauss_epsil, gauss_epsil_delta is obtained in UpdateFarm
            epsil = gauss_epsil * delta
        case default
            epsil = gauss_epsil
        end select

        GaussianInterpolation = exp( -(POW2(x(1) - x_point(1)) + POW2(x(2) - x_point(2)) + POW2(x(3) - x_point(3))) / POW2(epsil) ) / ( POW3(epsil) * pi**(3.0_RP/2.0_RP) )

    End Function GaussianInterpolation
!
!///////////////////////////////////////////////////////////////////////////////////////
!
! based on HexMesh_FindPointWithCoords, without curvature only in the precalculated list
    Subroutine FindActuatorPointElement(mesh, x, eID, xi, success)
       use HexMeshClass
       Implicit None

       type(HexMesh), intent(in)                     :: mesh
       real(kind=RP), dimension(NDIM), intent(in)    :: x       ! physical space
       integer, intent(out)                          :: eID 
       real(kind=RP), dimension(NDIM), intent(out)   :: xi      ! computational space
       logical, intent(out)                          :: success
       !
       logical                                       :: found
       integer                                       :: eIndex

       success = .false.
       found = .false.

!      Search in linear (not curved) mesh (faster and safer)
!      For AL the mesh is expected to be linear
!      -----------------------------------------------------

       do eIndex = 1, size(elementsActuated)
          eID = elementsActuated(eIndex)
          found = mesh % elements(eID) % FindPointInLinElement(x, mesh % nodes)
          if ( found ) exit
       end do

!      If found in linear mesh, use FindPointWithCoords in that element
       if (found) then
           success = mesh % elements(eID) % FindPointWithCoords(x, mesh % dir2D_ctrl, xi)
       end if
!
    End Subroutine FindActuatorPointElement
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    subroutine Get_Cl_Cl_from_airfoil_data(airfoil, aoa, Re, Cl_out, Cd_out)
         implicit none
      
         type (airfoil_t), intent(in)   :: airfoil
         real(KIND=RP), intent(in)      :: aoa, Re
         real(KIND=RP), intent(out)     :: Cl_out, Cd_out
         integer                        :: i,k
         real(kind=RP), dimension(2)    :: Cl_inter, Cd_inter
               
         Cl_out=0.0_RP
         Cd_out=0.0_RP

         if (airfoil%num_Re .eq. 1) then
             do i=1, airfoil % num_aoa-1
                 if (airfoil%aoa(i+1)>=aoa .and. airfoil%aoa(i)<=aoa ) then
                    Cl_out=InterpolateAirfoilData(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cl(i,1),airfoil%cl(i+1,1),aoa)
                    Cd_out=InterpolateAirfoilData(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cd(i,1),airfoil%cd(i+1,1),aoa)
                    exit
                 endif
             end do
        
         else

             do k=1, airfoil % num_Re-1
                 if (airfoil%Re(k+1)>=Re .and. airfoil%Re(k)<=Re ) then
                     do i=1, airfoil % num_aoa-1
                         if (airfoil%aoa(i+1)>=aoa .and. airfoil%aoa(i)<=aoa ) then

                            Cl_inter(1) = InterpolateAirfoilData(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cl(i,k),airfoil%cl(i+1,k),aoa)
                            Cl_inter(2) = InterpolateAirfoilData(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cl(i,k+1),airfoil%cl(i+1,k+1),aoa)
                            Cl_out=InterpolateAirfoilData(airfoil%Re(k),airfoil%Re(k+1),Cl_inter(1),Cl_inter(2),Re)

                            Cd_inter(1) = InterpolateAirfoilData(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cd(i,k),airfoil%cd(i+1,k),aoa)
                            Cd_inter(2) = InterpolateAirfoilData(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cd(i,k+1),airfoil%cd(i+1,k+1),aoa)
                            Cd_out=InterpolateAirfoilData(airfoil%Re(k),airfoil%Re(k+1),Cd_inter(1),Cd_inter(2),Re)

                            exit
                         endif
                     end do
                 endif
             end do

         end if

    end subroutine Get_Cl_Cl_from_airfoil_data

! linear interpolation given two points; returns y for new_x following line coefs (a,b) with y=ax+b
function InterpolateAirfoilData(x1,x2,y1,y2,new_x)
   implicit none
    
   real(KIND=RP), intent(in)    :: x1, x2, y1, y2, new_x
   real(KIND=RP)                :: a, b, InterpolateAirfoilData

    if(abs(x1-x2)<1.0e-6_RP) then 
      a=100.0_RP
   else
      a=(y1- y2)/(x1- x2)
   endif
    b= y1-a*x1;
    InterpolateAirfoilData=a*new_x+b
end function

function full_element_averageQ(mesh,eID)
   use HexMeshClass
   use PhysicsStorage
   use NodalStorageClass
   implicit none

   type(HexMesh), intent(in)    :: mesh
   integer, intent(in)          :: eID 
   integer                      :: k, j, i

   integer                      :: total_points
   real(kind=RP), dimension(NCONS)   :: full_element_averageQ, Qsum

 
   Qsum(:) = 0.0_RP
   total_points = 0
   do k = 0, mesh%elements(eID) % Nxyz(3)   ; do j = 0, mesh%elements(eID) % Nxyz(2) ; do i = 0, mesh%elements(eID) % Nxyz(1)
       Qsum(:)=Qsum(:)+mesh%elements(eID) % Storage % Q(:,i,j,k)
       total_points=total_points + 1
   end do                  ; end do                ; end do

   full_element_averageQ(:) = Qsum(:) / real(total_points,RP)

end function full_element_averageQ

Function semi_element_averageQ(mesh,eID,xi) result(Qe)
   use HexMeshClass
   use PhysicsStorage
   use NodalStorageClass

   Implicit None
   type(HexMesh), intent(in)    :: mesh
   integer, intent(in)          :: eID 
   real(kind=RP), dimension(NDIM), intent(in) :: xi
   real(kind=RP), dimension(NCONS)   :: Qe

   real(kind=RP), dimension(NCONS)   :: Qsum

   integer                      :: k, j, i, direction, N, ind
   integer, dimension(NDIM)     :: firstNodeIndex
   integer                      :: total_points
   type(NodalStorage_t), pointer :: spAxi

   ! fist get the sub element nodes index
   do direction = 1, NDIM

     N = mesh % elements(eID) % Nxyz(direction)
     spAxi   => NodalStorage(N)

     do ind = 0, N
         firstNodeIndex(direction) = ind-1
         if (xi(direction) .le. spAxi%x(ind)) exit
     end do

     if (firstNodeIndex(direction) .eq. -1) firstNodeIndex(direction) = 0

   end do

   nullify(spAxi)

   ! now average on the sub element
   Qsum(:) = 0.0_RP
   total_points = 0
   do k = firstNodeIndex(IZ), firstNodeIndex(IZ)+1   ; do j = firstNodeIndex(IY), firstNodeIndex(IY)+1 ; do i = firstNodeIndex(IX),firstNodeIndex(IX)+1
       Qsum(:) = Qsum(:) + mesh % elements(eID) % Storage % Q(:,i,j,k)
       total_points = total_points + 1
   end do                  ; end do                ; end do

   Qe(:) = Qsum(:) / real(total_points,RP)

End Function semi_element_averageQ

Function element_averageQ(mesh,eID,xi,averageSubElement) result(Qe)
   use HexMeshClass
   use PhysicsStorage
   Implicit None
   type(HexMesh), intent(in)    :: mesh
   integer, intent(in)          :: eID 
   logical, intent(in)          :: averageSubElement
   real(kind=RP), dimension(NDIM), intent(in) :: xi
   real(kind=RP), dimension(NCONS)   :: Qe

    if (averageSubElement) then
        Qe = semi_element_averageQ(mesh, eid, xi)
    else
        Qe = full_element_averageQ(mesh, eid)
    end if 
End Function element_averageQ

#endif
end module 

