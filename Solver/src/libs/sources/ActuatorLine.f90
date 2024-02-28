#include "Includes.h"

module ActuatorLine
#if defined(NAVIERSTOKES)
    use SMConstants
    use MPI_Process_Info
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

private
public farm
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
    real(KIND=RP), allocatable      :: local_lift(:)  ! blade sectional Lift
    real(KIND=RP), allocatable      :: local_drag(:)  ! blade sectional Lift
    real(KIND=RP), allocatable      :: point_xyz_loc(:,:) ! x,y location of blade points
    real(KIND=RP), allocatable      :: local_torque(:)  ! N.m
    real(KIND=RP), allocatable      :: local_thrust(:)  ! N
    real(KIND=RP), allocatable      :: local_thrust_temp(:)  ! N
    real(KIND=RP), allocatable      :: local_rotor_force(:)  ! N
    real(KIND=RP), allocatable      :: local_rotor_force_temp(:)  ! N
    real(KIND=RP), allocatable      :: local_root_bending(:)  ! N.m
    real(KIND=RP)                   :: gauss_epsil  ! force Gaussian shape
    real(KIND=RP)                   :: tip_c1,tip_c2  ! tip force corrections
    real(KIND=RP), allocatable      :: local_gaussian_sum(:) ! necessary for Gaussian weighted average
    real(KIND=RP), allocatable      :: local_Re(:) ! local Re based on local conditions and the chord of the airfoil at the blade section
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
    integer                        :: epsilon_type
    logical                        :: calculate_with_projection
    logical                        :: active = .false.
    logical                        :: save_average = .false.
    logical                        :: save_instant = .false.
    character(len=LINE_LENGTH)     :: file_name
    integer                        :: number_iterations
    integer                        :: save_iterations

   contains
        procedure   :: ConstructFarm
        procedure   :: DestructFarm
        procedure   :: UpdateFarm
        procedure   :: ForcesFarm             
        procedure   :: WriteFarmForces
        procedure   :: GaussianInterpolation
        procedure   :: FarmUpdateLocalForces
        procedure   :: FarmUpdateBladeForces
    end type

    type(Farm_t)                  :: farm

     ! max 10 airfoils file names per section
     integer, parameter           :: MAX_AIRFOIL_FILES = 10

!  ========
contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructFarm(self, controlVariables, t0)
       use FTValueDictionaryClass
       use mainKeywordsModule, only: solutionFileNameKey
       use FileReadingUtilities      , only: getFileName
       use PhysicsStorage
       use fluiddata
       use MPI_Process_Info
       implicit none
       class(farm_t) , intent(inout)                :: self
       TYPE(FTValueDictionary), intent(in)          :: controlVariables
       real(kind=RP), intent(in)                    :: t0
!        ---------------
!        Local variables
!        ---------------
!
         integer     ::  i, j, k, ii, fid, io, n_aoa
         CHARACTER(LEN=LINE_LENGTH) :: arg, char1
         CHARACTER(LEN=LINE_LENGTH) :: solution_file
         real(kind=RP)              :: initial_azimutal

    if (.not. controlVariables % logicalValueForKey("use actuatorline")) return

    self % epsilon_type = controlVariables % getValueOrDefault("actuator epsilon type", 0)
    self % calculate_with_projection = controlVariables % getValueOrDefault("actuator calculate with projection", .false.)
    self % save_average = controlVariables % getValueOrDefault("actuator save average", .false.)
    self % save_instant = controlVariables % getValueOrDefault("actuator save instant", .false.)
    self % save_iterations = controlVariables % getValueOrDefault("actuator save iteration", 1)

    arg='./ActuatorDef/Act_ActuatorDef.dat'
    OPEN( newunit = fid,file=trim(arg),status="old",action="read")

    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,*) self%num_turbines

    print *,'-------------------------'
    print *,achar(27)//'[34m READING FARM DEFINITION'
    write(*,*) "Number of turbines in farm:", self%num_turbines

    READ(fid,'(A132)') char1

    allocate(self%turbine_t(self%num_turbines))

    do i = 1, self%num_turbines
       READ(fid,*) self%turbine_t(i)%hub_cood_x, self%turbine_t(i)%hub_cood_y, self%turbine_t(i)%hub_cood_z
    ENDDO

    READ(fid,'(A132)') char1

    do i = 1, self%num_turbines
       READ(fid,*) self%turbine_t(i)%radius
    ENDDO

    READ(fid,'(A132)') char1

    do i = 1, self%num_turbines
       READ(fid,*) self%turbine_t(i)%normal_x, self%turbine_t(i)%normal_y, self%turbine_t(i)%normal_z
    ENDDO
    
    READ(fid,'(A132)') char1

    do i = 1, self%num_turbines
       READ(fid,*) self%turbine_t(i)%rot_speed
    ENDDO

    READ(fid,'(A132)') char1

    do i = 1, self%num_turbines
       READ(fid,*) self%turbine_t(i)%blade_pitch
    ENDDO
    

    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1

    ! Read blade info, we assume all 3 blades are the same for one turbine
    READ(fid,*) self%turbine_t(1)%num_blade_sections

    write(*,*) "Number of blade sections:", self%turbine_t(1)%num_blade_sections

     associate (num_blade_sections => self%turbine_t(1)%num_blade_sections)

     READ(fid,'(A132)') char1

     do i=1, self%num_turbines
      do j=1, self%turbine_t(i)%num_blades
     allocate( self%turbine_t(i)%blade_t(j)%r_R(0:num_blade_sections),self%turbine_t(i)%blade_t(j)%chord(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%twist(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%num_airfoils(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%airfoil_files(num_blade_sections,MAX_AIRFOIL_FILES),self%turbine_t(i)%blade_t(j)%airfoil_t(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%local_velocity(num_blade_sections), self%turbine_t(i)%blade_t(j)%local_angle(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%local_lift(num_blade_sections), self%turbine_t(i)%blade_t(j)%local_drag(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%point_xyz_loc(num_blade_sections,3),self%turbine_t(i)%blade_t(j)%local_torque(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%local_thrust(num_blade_sections),self%turbine_t(i)%blade_t(j)%local_root_bending(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%local_rotor_force(num_blade_sections),self%turbine_t(i)%blade_t(j)%local_gaussian_sum(num_blade_sections), &
     self%turbine_t(i)%blade_t(j)%local_Re(num_blade_sections) )

         do k=1, num_blade_sections
            self%turbine_t(i)%blade_t(j)%airfoil_files(k,:)=' '
         enddo
      ENDDO
   enddo

    if (self%calculate_with_projection) then
        do i=1, self%num_turbines
          do j=1, self%turbine_t(i)%num_blades
         allocate( self%turbine_t(i)%blade_t(j)%local_thrust_temp(num_blade_sections), &
                   self%turbine_t(i)%blade_t(j)%local_rotor_force_temp(num_blade_sections) )
          enddo
       enddo
   end if

    endassociate

    self%turbine_t(1)%blade_t(1)%r_R(0) = 0.0_RP
   do i = 1, self%turbine_t(1)%num_blade_sections
      READ(fid,*) self%turbine_t(1)%blade_t(1)%r_R(i), self%turbine_t(1)%blade_t(1)%chord(i), &
                  self%turbine_t(1)%blade_t(1)%twist(i), self%turbine_t(1)%blade_t(1)%num_airfoils(i)
      
      ! one file per Re for each airfoil
      self%turbine_t(1)%blade_t(1)%airfoil_t(i)%num_Re = self%turbine_t(1)%blade_t(1)%num_airfoils(i)
                    
      do j = 1, self%turbine_t(1)%blade_t(1)%num_airfoils(i)   
            READ(fid,*) self%turbine_t(1)%blade_t(1)%airfoil_files(i,j)  
      enddo
   ENDDO


     ! all turbines have the same blades
     !do i=1, self%num_turbines
     !    do j=1, self%turbine_t(i)%num_blades
     !       self%turbine_t(i)%blade_t(j)=self%turbine_t(1)%blade_t(1)
     !    ENDDO
     ! enddo
     ! write(*,*) "All turbines have the same blades"
   !  write(*,*) self%turbine_t(1)%blade_t(1)%airfoil_files(2,2)

! read numerical parameters
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1

     READ(fid,*) self%turbine_t(1)%blade_t(1)%gauss_epsil    

     READ(fid,'(A132)') char1
     READ(fid,*) self%turbine_t(1)%blade_t(1)%tip_c1,self%turbine_t(1)%blade_t(1)%tip_c2


    print*,'Gaussian value for actuator line',self%turbine_t(1)%blade_t(1)%gauss_epsil
    print*,'Tip correction constants', self%turbine_t(1)%blade_t(1)%tip_c1,self%turbine_t(1)%blade_t(1)%tip_c2

     print *, "Use initial azimuthal angle", controlVariables%logicalValueForKey("actuator initial azimuthal")
     print *, "Use projection formulation: ", self % calculate_with_projection
     print *, "Save blade average values: ", self % save_average

    close(fid)

  do i = 1, self%turbine_t(1)%num_blade_sections
      arg=trim('./ActuatorDef/'//trim(self%turbine_t(1)%blade_t(1)%airfoil_files(i,1)))
      OPEN( newunit = fid,file=trim(arg),status="old",action="read")
      READ(fid,'(A132)') char1
      READ(fid,*) self%turbine_t(1)%blade_t(1)%airfoil_t(i)%num_aoa 
      close(fid)

      associate (num_re => self%turbine_t(1)%blade_t(1)%airfoil_t(i)%num_Re, num_aoa => self%turbine_t(1)%blade_t(1)%airfoil_t(i)%num_aoa)

         do ii=1, self%num_turbines
            do j=1, self%turbine_t(ii)%num_blades
                  allocate( self%turbine_t(ii)%blade_t(j)%airfoil_t(i)%aoa(num_aoa), &
                        self%turbine_t(ii)%blade_t(j)%airfoil_t(i)%cl(num_aoa,num_re), &
                        self%turbine_t(ii)%blade_t(j)%airfoil_t(i)%cd(num_aoa,num_re), &
                        self%turbine_t(ii)%blade_t(j)%airfoil_t(i)%Re(num_re) )
            enddo
         enddo

         do k = 1, num_re

            arg=trim('./ActuatorDef/'//trim(self%turbine_t(1)%blade_t(1)%airfoil_files(i,k)))
            OPEN( newunit = fid,file=trim(arg),status="old",action="read")
            READ(fid,'(A132)') char1

            READ(fid,*) n_aoa
            if ( n_aoa .ne. num_aoa ) then
                print *, "Error: not same number of AoA in all files for same blade section, file: ", trim(arg)
                call exit(99)
            end if

            READ(fid,'(A132)') char1
            READ(fid,*) self%turbine_t(1)%blade_t(1)%airfoil_t(i)%Re(k)

            print *,'-------------------------'
            print *,achar(27)//'[34m READING FARM AIRFOIL DATA (Cl-Cd)'
            print*, 'reading: ', trim(arg)
            write(*,*) 'The number of AoA in the file is: ', num_aoa,' '//achar(27)//'[0m '    
    
            READ(fid,'(A132)') char1

           do ii = 1,  num_aoa
                READ(fid,*) self%turbine_t(1)%blade_t(1)%airfoil_t(i)%aoa(ii), self%turbine_t(1)%blade_t(1)%airfoil_t(i)%cl(ii,k), &
                            self%turbine_t(1)%blade_t(1)%airfoil_t(i)%cd(ii,k)
                ! file is in deg, convert to rad
                self%turbine_t(1)%blade_t(1)%airfoil_t(i)%aoa(ii) = self%turbine_t(1)%blade_t(1)%airfoil_t(i)%aoa(ii) * PI / 180.0_RP
           enddo

           close(fid)
         end do ! number of airfoil files

      endassociate

    !all airfoils of all blades of all turbines are the same
    !do ii=1, self%num_turbines
    !  do j=1, self%turbine_t(i)%num_blades
    !     do k=1, self%turbine_t(1)%blade_t(1)%num_airfoils(j) 
    !     self%turbine_t(ii)%blade_t(j)%airfoil_t(k)=self%turbine_t(1)%blade_t(1)%airfoil_t(1)
    !     enddo
    !  ENDDO
    !enddo
      
 enddo ! number of blade sections

!$omp do schedule(runtime)private(i)
   do i = 1, self%turbine_t(1)%num_blade_sections
      self%turbine_t(1)%blade_t(1)%point_xyz_loc(i,1) =   self%turbine_t(1)%hub_cood_x
   enddo
!$omp end do

   !all airfoils of all blades of all turbines are the same
   do ii=1, self%num_turbines
      do j=1, self%turbine_t(ii)%num_blades 
         self%turbine_t(ii)%blade_t(j)=self%turbine_t(1)%blade_t(1)
         enddo
   enddo

    ! azimuthal angle for the 3 blades
    ! initial azimuthal angle valid for restaring a simulation with same rotational speed and refValues
    if (controlVariables % logicalValueForKey("actuator initial azimuthal")) then
        initial_azimutal = self%turbine_t(1)%rot_speed*(t0*Lref/refValues%V)
    else
        initial_azimutal = 0.0_RP
    end if
    ! azimuth_angle angle of blades is the angle to respect to +y axis, the angular velocity vector will point to +x
    self%turbine_t(:)%blade_t(1)%azimuth_angle = initial_azimutal
    self%turbine_t(:)%blade_t(2)%azimuth_angle = initial_azimutal + PI*2.0_RP/3.0_RP
    self%turbine_t(:)%blade_t(3)%azimuth_angle = initial_azimutal + PI*4.0_RP/3.0_RP

    ! average_conditions are thrust and rotor force
    if (MPI_Process % isRoot) then
      if (self % save_average) then
          if (self % calculate_with_projection) then
              allocate ( self%turbine_t(1)%average_conditions(self%turbine_t(1)%num_blade_sections,2) )
          else
              ! additional average_conditions are velocity, AoA and Re, writen before forces
              allocate ( self%turbine_t(1)%average_conditions(self%turbine_t(1)%num_blade_sections,5) )
          end if
          self % turbine_t(1) % average_conditions(:,:) = 0.0_RP
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
!   Create output files
!   -------------------
    if (MPI_Process % isRoot) then
      write(arg , '(A,A)') trim(self%file_name) , "_Actuator_Line_Forces.dat"
      open ( newunit = fID , file = trim(arg) , status = "unknown" ,    action = "write" ) 
      write(fid,*) 'time, thrust_1, blade_torque_1, blade_root_bending_1,thrust_2, blade_torque_12 blade_root_bending_2,thrust_3, blade_torque_3, blade_root_bending_3'
      close(fid)
!
      write(arg , '(A,A)') trim(self%file_name) , "_Actuator_Line_CP_CT.dat"
      open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" ) 
      write(fid,*) 'time, Cp (power coef.), Ct (thust coef.)'
      close(fid)
!
      if (self % save_average) then
          write(arg , '(A,A)') trim(self%file_name) , "_Actuator_Line_average.dat"
          open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" ) 
          if (self % calculate_with_projection) then
              write(fid,*) 'R, Tangential_Force, Axial_Force'
          else
              write(fid,*) 'R, U, AoA, Re, Tangential_Force, Axial_Force'
          end if
          close(fid)
      end if
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
   !                  self%turbine_t(i)%blade_t(j)%local_lift, &
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
   integer                           :: ii, jj, i, j, k
   real(kind=RP)                     :: theta,t, interp, tolerance
   logical                           :: found, allfound
   integer                           :: eID, ierr
   real(kind=RP), dimension(NDIM)    :: x, xe
   real(kind=RP), dimension(NCONS)   :: Q, Qtemp
   real(kind=RP), dimension(:), allocatable  :: aoa

    if (.not. self % active) return

   t = time * Lref / refValues%V
   theta = self%turbine_t(1)%rot_speed * t
   interp = 1.0_RP
   tolerance=0.2_RP*self%turbine_t(1)%radius

   projection_cond:if (self%calculate_with_projection) then
!
!    ----------------------------------------------------------------------------------
!    calculate for all mesh points its contribution based on the gaussian interpolation
!    ----------------------------------------------------------------------------------
!
      if ( (MPI_Process % doMPIAction) ) then
            print*, "MPI not implemented yet for AL projection mode"
            call exit(99)
      end if
!$omp do schedule(runtime)private(ii,jj)
      do jj = 1, self%turbine_t(1)%num_blades

         self%turbine_t(1)%blade_t(jj)%local_lift(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_drag(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_rotor_force(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_rotor_force_temp(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_thrust(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_thrust_temp(:)=0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_torque(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_root_bending(:) = 0.0_RP
         self%turbine_t(1)%blade_t(jj)%local_gaussian_sum(:)= 0.0_RP
      
             do ii = 1, self%turbine_t(1)%num_blade_sections
                ! y,z coordinate of every acutator line point
                self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,2) = self%turbine_t(1)%hub_cood_y + self%turbine_t(1)%blade_t(jj)%r_R(ii) * cos(theta+self%turbine_t(1)%blade_t(jj)%azimuth_angle)
                self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,3) = self%turbine_t(1)%hub_cood_z + self%turbine_t(1)%blade_t(jj)%r_R(ii) * sin(theta+self%turbine_t(1)%blade_t(jj)%azimuth_angle)
      
              end do
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
!$omp do schedule(runtime)private(ii,jj,eID,Q,Qtemp,x,xe,found)
     do jj = 1, self%turbine_t(1)%num_blades
       self%turbine_t(1)%blade_t(jj)%local_lift(:) = 0.0_RP
       self%turbine_t(1)%blade_t(jj)%local_drag(:) = 0.0_RP
       self%turbine_t(1)%blade_t(jj)%local_rotor_force(:) = 0.0_RP
       self%turbine_t(1)%blade_t(jj)%local_thrust(:) = 0.0_RP
       self%turbine_t(1)%blade_t(jj)%local_torque(:) = 0.0_RP
       self%turbine_t(1)%blade_t(jj)%local_root_bending(:) = 0.0_RP
!
       do ii = 1, self%turbine_t(1)%num_blade_sections
          ! y,z coordinate of every acutator line point
          self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,2) = self%turbine_t(1)%hub_cood_y + self%turbine_t(1)%blade_t(jj)%r_R(ii) * cos(theta+self%turbine_t(1)%blade_t(jj)%azimuth_angle)
          self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,3) = self%turbine_t(1)%hub_cood_z + self%turbine_t(1)%blade_t(jj)%r_R(ii) * sin(theta+self%turbine_t(1)%blade_t(jj)%azimuth_angle)
!
!         -----------------------------------
!         get the elements of each line point
!         -----------------------------------
!
          x = [self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,1),self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,2),self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,3)]
          found = mesh % FindPointWithCoords(x, eID, xe)
          if (found) then
            ! averaged state values of the cell
            Qtemp = element_averageQ(mesh,eID)
          else
            Qtemp = 0.0_RP
          end if
          if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
            call mpi_allreduce(Qtemp, Q, NCONS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
          else
              Q = Qtemp
          end if
          if (all(Q .eq. 0.0_RP)) then
            print*, "Actuator line point not found in mesh, x: ", x
            call exit(99)
          end if
          call FarmUpdateLocalForces(self, ii, jj, Q, theta, interp)
        end do
     enddo
!$omp end do

   end if projection_cond
!
   end subroutine UpdateFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ForcesFarm(self, x, Q, NS, time)
   use PhysicsStorage
   use fluiddata
   implicit none

   class(Farm_t) , intent(inout)     :: self
   real(kind=RP),intent(in)          :: x(NDIM)
   real(kind=RP),intent(in)          :: Q(NCONS)
   real(kind=RP),intent(inout)       :: NS(NCONS)
   real(kind=RP),intent(in)          :: time

! local vars
   real(kind=RP)                     :: tolerance, Non_dimensional, t, theta, interp, local_gaussian
   integer                           :: ii,jj, LAST_SECTION
   real(kind=RP), dimension(NDIM)    :: actuator_source

    if (.not. self % active) return

    ! 20% of the radius for max of the rotor thickness
    tolerance=0.2_RP*self%turbine_t(1)%radius

        ! turbine is pointing backwards as x positive
if( POW2(x(2)-self%turbine_t(1)%hub_cood_y)+POW2(x(3)-self%turbine_t(1)%hub_cood_z) <= POW2(self%turbine_t(1)%radius+tolerance) &
    .and. (x(1) < self%turbine_t(1)%hub_cood_x+tolerance .and. x(1)>self%turbine_t(1)%hub_cood_x-tolerance)) then

      Non_dimensional = POW2(refValues % V) * refValues % rho / Lref
      t = time * Lref / refValues % V

      theta = self%turbine_t(1)%rot_speed * t

      actuator_source(:) = 0.0_RP
      local_gaussian=0.0_RP

 if (self%calculate_with_projection) then

 
   do jj = 1, self%turbine_t(1)%num_blades
      do ii = 1, self%turbine_t(1)%num_blade_sections
          interp = GaussianInterpolation(self, ii, jj, x)
          call FarmUpdateLocalForces(self, ii, jj,  Q, theta, interp)

          ! minus account action-reaction effect, is the force on the fliud
          actuator_source(1) = actuator_source(1) - self%turbine_t(1)%blade_t(jj)%local_thrust(ii) 
          actuator_source(2) = actuator_source(2) - (-self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)*sin(self%turbine_t(1)%rot_speed*t + self%turbine_t(1)%blade_t(jj)%azimuth_angle) )
          actuator_source(3) = actuator_source(3) - self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)*cos(self%turbine_t(1)%rot_speed*t + self%turbine_t(1)%blade_t(jj)%azimuth_angle) 

          !acumulate in temporal variables, for each time step as the non temp are recalculated for each element
          self%turbine_t(1)%blade_t(jj)%local_thrust_temp(ii)=self%turbine_t(1)%blade_t(jj)%local_thrust_temp(ii)+self%turbine_t(1)%blade_t(jj)%local_thrust(ii)
          self%turbine_t(1)%blade_t(jj)%local_rotor_force_temp(ii)=self%turbine_t(1)%blade_t(jj)%local_rotor_force_temp(ii)+self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)

         self%turbine_t(1)%blade_t(jj)%local_gaussian_sum(ii)=self%turbine_t(1)%blade_t(jj)%local_gaussian_sum(ii)+interp
     
         local_gaussian=local_gaussian+interp

      enddo
  enddo
  
     ! actuator_source(:)=actuator_source(:)/local_gaussian
    
    NS(IRHOU:IRHOW) = NS(IRHOU:IRHOW) + actuator_source(:) / Non_dimensional



else ! no projection

        ! LAST_SECTION=self%turbine_t(1)%num_blade_sections

       do jj = 1, self%turbine_t(1)%num_blades
          
          do ii = 1, self%turbine_t(1)%num_blade_sections

            interp = GaussianInterpolation(self, ii, jj, x)
    
            ! minus account action-reaction effect, is the force on the fliud
            actuator_source(1) = actuator_source(1) - self%turbine_t(1)%blade_t(jj)%local_thrust(ii) * interp
            actuator_source(2) = actuator_source(2) - (-self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)*sin(self%turbine_t(1)%rot_speed*t + self%turbine_t(1)%blade_t(jj)%azimuth_angle) )
            actuator_source(3) = actuator_source(3) - self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)*cos(self%turbine_t(1)%rot_speed*t + self%turbine_t(1)%blade_t(jj)%azimuth_angle) 

            !local_gaussian=local_gaussian+interp

         enddo
     enddo
     
        !actuator_source(:)=actuator_source(:)/local_gaussian

       NS(IRHOU:IRHOW) = NS(IRHOU:IRHOW) + actuator_source(:) / Non_dimensional

    endif
    
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
   integer                       :: fid, io
   CHARACTER(LEN=LINE_LENGTH)    :: arg
   real(kind=RP)                 :: t
   integer                       :: ii, jj
   logical                       :: saveAverage
   logical                       :: save_instant

    if (.not. self % active) return
    if ( .not. MPI_Process % isRoot ) return

    if (present(last)) then
        saveAverage = last
    else
        saveAverage = .false.
    end if

   save_instant = self%save_instant .and. ( mod(iter,self % save_iterations) .eq. 0 )
   t = time * Lref / refValues%V

   
 if (self%calculate_with_projection) then
   ! this is necessary for Gaussian weighted sum
   
   do jj = 1, self%turbine_t(1)%num_blades
       self%turbine_t(1)%blade_t(jj)%local_thrust(:) = 0.0_RP
       self%turbine_t(1)%blade_t(jj)%local_rotor_force(:) = 0.0_RP
   end do

   !$omp do schedule(runtime)private(ii,jj)
      do jj = 1, self%turbine_t(1)%num_blades
     
           do ii = 1, self%turbine_t(1)%num_blade_sections
   
               self%turbine_t(1)%blade_t(jj)%local_thrust(ii)=self%turbine_t(1)%blade_t(jj)%local_thrust(ii)+self%turbine_t(1)%blade_t(jj)%local_thrust_temp(ii)/self%turbine_t(1)%blade_t(jj)%local_gaussian_sum(ii)
               self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)=self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)+self%turbine_t(1)%blade_t(jj)%local_rotor_force_temp(ii)/self%turbine_t(1)%blade_t(jj)%local_gaussian_sum(ii)
   
           enddo
       enddo
   !$omp end do
   
   end if

   ! save in memory the time step forces for each element blade and the whole blades
   call self % FarmUpdateBladeForces()

!write output torque thrust to file
      write(arg , '(A,A)') trim(self%file_name) , "_Actuator_Line_Forces.dat"

      open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
      write(fid,"(10(2X,ES24.16))") t, &
      self%turbine_t(1)%blade_thrust(1),self%turbine_t(1)%blade_torque(1),self%turbine_t(1)%blade_root_bending(1), &
      self%turbine_t(1)%blade_thrust(2),self%turbine_t(1)%blade_torque(2),self%turbine_t(1)%blade_root_bending(2), &
      self%turbine_t(1)%blade_thrust(3),self%turbine_t(1)%blade_torque(3),self%turbine_t(1)%blade_root_bending(3)
      close(fid)

      write(arg , '(A,A)') trim(self%file_name) , "_Actuator_Line_CP_CT.dat"
      open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
      write(fid,"(10(2X,ES24.16))") t, self%turbine_t(1)%Cp, self%turbine_t(1)%Ct
      close(fid)

    if (self % save_average .and. saveAverage) then
      write(arg , '(A,A)') trim(self%file_name) , "_Actuator_Line_average.dat"
      open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
      if (self%calculate_with_projection) then
          do ii = 1, self % turbine_t(1) % num_blade_sections
              write(fid,"(3(2X,ES24.16))") self%turbine_t(1)%blade_t(1)%r_R(ii), self%turbine_t(1)%average_conditions(ii,:)
          end do
      else
          do ii = 1, self % turbine_t(1) % num_blade_sections
              write(fid,"(6(2X,ES24.16))") self%turbine_t(1)%blade_t(1)%r_R(ii), self%turbine_t(1)%average_conditions(ii,:)
          end do
      end if 
      close(fid)
    end if

     if (save_instant) then
         do jj = 1, self%turbine_t(1)%num_blades
          write(arg , '(2A,I3.3,A,I10.10,A)') trim(self%file_name) , "_Actuator_Line_instant_",jj ,"_" ,iter, ".dat"
          open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" ) 
          write(fid,*) 'R, U, AoA, Re'
          do ii = 1, self % turbine_t(1) % num_blade_sections
              write(fid,"(4(2X,ES24.16))") self%turbine_t(1)%blade_t(1)%r_R(ii), &
                 self%turbine_t(1)%blade_t(jj)%local_velocity(ii), &
                 ( self%turbine_t(1)%blade_t(jj)%local_angle(ii) - (self%turbine_t(1)%blade_t(jj)%twist(ii) + self%turbine_t(1)%blade_pitch) ) * 180.0_RP / PI, &
                 self%turbine_t(1)%blade_t(jj)%local_Re(ii)
          end do
          close(fid)
         end do
     end if 
  !endif

end subroutine WriteFarmForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Subroutine FarmUpdateLocalForces(self, ii, jj, Q, theta, interp)
        use PhysicsStorage
        use fluiddata
        use VariableConversion, only: Temperature, SutherlandsLaw
        implicit none
        class(Farm_t)                                 :: self
        integer, intent(in)                           :: ii, jj
        real(kind=RP), dimension(NCONS), intent(in)   :: Q
        real(kind=RP), intent(in)                     :: theta
        real(kind=RP), intent(in)                     :: interp

        !local variables
        real(kind=RP)                                 :: density, Cl, Cd, aoa, g1_func, tip_correct, angle_temp
        real(kind=RP)                                 :: wind_speed_axial, wind_speed_rot
        real(kind=RP)                                 :: lift_force, drag_force
        real(kind=RP)                                 :: T, muL
!
!       -----------------------------
!       get airfoil related variables
!       -----------------------------
!
        ! option 1, not recommended, ignore LES velocity directions, just use U0 in x-direction
        ! wind_speed_axial =  refValues % V ! wind goes in the x-direction
        ! wind_speed_axial = sqrt( POW2(Q(IRHOU)) + POW2(Q(IRHOV)) ) / Q(IRHO) * refValues_%V ! wind goes in the x-direction
        ! wind_speed_rot = 0.0_RP

        ! option 2, project [v.w] in the rotational direction (theta in cylindrical coordinates)
        wind_speed_axial = (Q(IRHOU)/Q(IRHO)) * refValues % V ! our x is the z in cylindrical
        wind_speed_rot = ( -Q(IRHOV)*sin(theta+self%turbine_t(1)%blade_t(jj)%azimuth_angle) + Q(IRHOW)*cos(theta+self%turbine_t(1)%blade_t(jj)%azimuth_angle) ) / Q(IRHO) * refValues % V

        density = Q(IRHO) * refValues % rho

        tip_correct = 1.0_RP
        aoa = 0.0_RP

         T     = Temperature(Q)
         muL = SutherlandsLaw(T) * refValues % mu

        self%turbine_t(1)%blade_t(jj)%local_velocity(ii) = sqrt( POW2(self%turbine_t(1)%rot_speed*self%turbine_t(1)%blade_t(jj)%r_R(ii) - wind_speed_rot) + &
                                                                 POW2(wind_speed_axial) )
        self%turbine_t(1)%blade_t(jj)%local_angle(ii) = atan( wind_speed_axial / (self%turbine_t(1)%rot_speed*self%turbine_t(1)%blade_t(jj)%r_R(ii) - wind_speed_rot) ) 

        ! alpha = phi - gamma, gamma = blade pitch + airfoil local twist
        aoa = self%turbine_t(1)%blade_t(jj)%local_angle(ii) - (self%turbine_t(1)%blade_t(jj)%twist(ii) + self%turbine_t(1)%blade_pitch) 
        self%turbine_t(1)%blade_t(jj)%local_Re(ii) = self%turbine_t(1)%blade_t(jj)%local_velocity(ii) * self%turbine_t(1)%blade_t(jj)%chord(ii) * density / muL
        call Get_Cl_Cl_from_airfoil_data(self%turbine_t(1)%blade_t(jj)%airfoil_t(ii), aoa, self%turbine_t(1)%blade_t(jj)%local_Re(ii), Cl, Cd)
!
!       ---------------
!       tip correction
!       ---------------
!
        g1_func=exp(-self%turbine_t(1)%blade_t(1)%tip_c1*(self%turbine_t(1)%num_blades*self%turbine_t(1)%rot_speed*self%turbine_t(1)%radius/refValues%V-self%turbine_t(1)%blade_t(1)%tip_c2))+0.1
        ! in this it should be the local_angle at the tip

        ! without perturbation
        ! angle_temp = atan(refValues%V/(self%turbine_t(1)%rot_speed*self%turbine_t(1)%radius))
        ! only axial wind speed
        ! angle_temp = atan(wind_speed_axial/(self%turbine_t(1)%rot_speed*self%turbine_t(1)%radius))

        ! 1 possibility: use the global radius and angle (i.e. the tip values)
        ! tip_correct = 2.0_RP/PI*(acos( exp(-g1_func*self%turbine_t(1)%num_blades*(self%turbine_t(1)%radius-self%turbine_t(1)%blade_t(jj)%r_R(ii)) / &
                      ! (2.0_RP*self%turbine_t(1)%radius*sin(angle_temp))) ))

        ! 2 possibility: use the local radius and angle
        tip_correct = 2.0_RP/PI*(acos( exp(-g1_func*self%turbine_t(1)%num_blades*(self%turbine_t(1)%radius-self%turbine_t(1)%blade_t(jj)%r_R(ii)) / abs(2.0_RP*self%turbine_t(1)%blade_t(jj)%r_R(ii)*sin(self%turbine_t(1)%blade_t(jj)%local_angle(ii)))) ))
!
!       --------------------------------
!       Save forces on the blade segment
!       --------------------------------
!
         ! lift=Cl*1/2.rho*v_local^2*Surface ; and surface=section area of the blade S=length*chord
         ! lift and drag are multiply by the gaussian interp for the case of each mesh node contributing to the force (for local only is 1)

         lift_force = 0.5_RP * density * Cl * tip_correct * POW2(self%turbine_t(1)%blade_t(jj)%local_velocity(ii)) &
                      * self%turbine_t(1)%blade_t(jj)%chord(ii) * (self%turbine_t(1)%blade_t(jj)%r_R(ii) - self%turbine_t(1)%blade_t(jj)%r_R(ii-1)) * interp

        drag_force = 0.5_RP * density * Cd * tip_correct * POW2(self%turbine_t(1)%blade_t(jj)%local_velocity(ii)) &
                      * self%turbine_t(1)%blade_t(jj)%chord(ii) * (self%turbine_t(1)%blade_t(jj)%r_R(ii) - self%turbine_t(1)%blade_t(jj)%r_R(ii-1)) * interp

        self%turbine_t(1)%blade_t(jj)%local_lift(ii) =  lift_force
        self%turbine_t(1)%blade_t(jj)%local_drag(ii) =  drag_force

        self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii) = lift_force * sin(self%turbine_t(1)%blade_t(jj)%local_angle(ii)) &
                                                              - drag_force * cos(self%turbine_t(1)%blade_t(jj)%local_angle(ii))
                                  
        self%turbine_t(1)%blade_t(jj)%local_thrust(ii) = lift_force * cos(self%turbine_t(1)%blade_t(jj)%local_angle(ii)) & 
                                                         + drag_force * sin(self%turbine_t(1)%blade_t(jj)%local_angle(ii))

!
    End Subroutine FarmUpdateLocalForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Subroutine FarmUpdateBladeForces(self)
        use fluiddata
        Implicit None

        class(Farm_t), intent(inout)      :: self
        !local variables
        integer                           :: ii, jj, i, j, k
        real(kind=RP), dimension(:), allocatable  :: aoa
!
!       ------------------------------
!       Save forces on the whole blade
!       ------------------------------
!
!$omp do schedule(runtime)private(ii,jj)
     do jj = 1, self%turbine_t(1)%num_blades
          self%turbine_t(1)%blade_thrust(jj) = 0.0_RP
          self%turbine_t(1)%blade_torque(jj) = 0.0_RP
          self%turbine_t(1)%blade_root_bending(jj) = 0.0_RP

          do ii = 1, self%turbine_t(1)%num_blade_sections

              self%turbine_t(1)%blade_t(jj)%local_torque(ii) = self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii)*self%turbine_t(1)%blade_t(jj)%r_R(ii)

              self%turbine_t(1)%blade_t(jj)%local_root_bending(ii) = sqrt(POW2(self%turbine_t(1)%blade_t(jj)%local_thrust(ii)) + &
                  POW2(self%turbine_t(1)%blade_t(jj)%local_rotor_force(ii))) * self%turbine_t(1)%blade_t(jj)%r_R(ii)

              self%turbine_t(1)%blade_thrust(jj)=self%turbine_t(1)%blade_thrust(jj)+self%turbine_t(1)%blade_t(jj)%local_thrust(ii)
              self%turbine_t(1)%blade_torque(jj)=self%turbine_t(1)%blade_torque(jj)+self%turbine_t(1)%blade_t(jj)%local_torque(ii)
              self%turbine_t(1)%blade_root_bending(jj)=self%turbine_t(1)%blade_root_bending(jj)+self%turbine_t(1)%blade_t(jj)%local_root_bending(ii)
         enddo
     enddo
!$omp end do

     self%turbine_t(1)%Cp = 2.0_RP * (self%turbine_t(1)%blade_torque(1)+self%turbine_t(1)%blade_torque(2)+self%turbine_t(1)%blade_torque(3)) * self%turbine_t(1)%rot_speed / &
                           (refValues%rho * POW3(refValues%V) * pi * POW2(self%turbine_t(1)%radius))

     self%turbine_t(1)%Ct = 2.0_RP * (self%turbine_t(1)%blade_thrust(1)+self%turbine_t(1)%blade_thrust(2)+self%turbine_t(1)%blade_thrust(3)) / &
                           (refValues%rho * POW2(refValues%V) * pi * POW2(self%turbine_t(1)%radius))
!
!        ----------------------
!        Save average variables
!        ----------------------
!
     if (self % save_average) then
         self % turbine_t(1) % average_conditions = self % turbine_t(1) % average_conditions * real(self % number_iterations,RP)

         !saving only for blade 1
         jj =1
         if (self%calculate_with_projection) then
             self % turbine_t(1) % average_conditions(:,1) = self % turbine_t(1) % average_conditions(:,1) + self%turbine_t(1)%blade_t(jj)%local_rotor_force
             self % turbine_t(1) % average_conditions(:,2) = self % turbine_t(1) % average_conditions(:,2) + self%turbine_t(1)%blade_t(jj)%local_thrust
         else
             allocate(aoa(self%turbine_t(1)%num_blade_sections))
             aoa = self%turbine_t(1)%blade_t(jj)%local_angle(:) - (self%turbine_t(1)%blade_t(jj)%twist(:) + self%turbine_t(1)%blade_pitch) 
             self % turbine_t(1) % average_conditions(:,1) = self % turbine_t(1) % average_conditions(:,1) + self%turbine_t(1)%blade_t(jj)%local_velocity
             self % turbine_t(1) % average_conditions(:,2) = self % turbine_t(1) % average_conditions(:,2) + aoa * 180.0_RP / PI
             self % turbine_t(1) % average_conditions(:,3) = self % turbine_t(1) % average_conditions(:,3) + self%turbine_t(1)%blade_t(jj)%local_Re
             self % turbine_t(1) % average_conditions(:,4) = self % turbine_t(1) % average_conditions(:,4) + self%turbine_t(1)%blade_t(jj)%local_rotor_force
             self % turbine_t(1) % average_conditions(:,5) = self % turbine_t(1) % average_conditions(:,5) + self%turbine_t(1)%blade_t(jj)%local_thrust
             deallocate(aoa)
         end if

         self % number_iterations = self % number_iterations + 1
         self % turbine_t(1) % average_conditions = self % turbine_t(1) % average_conditions / real(self % number_iterations,RP)
     end if 
!
    End Subroutine FarmUpdateBladeForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Function GaussianInterpolation(self, ii, jj, x, Cd)
        implicit none
        class(Farm_t), intent(in)               :: self
        integer, intent(in)                     :: ii, jj
        real(kind=RP), intent(in)               :: x(NDIM)
        real(kind=RP), intent(in), optional     :: Cd
        real(kind=RP)                           :: GaussianInterpolation

        !local variables
        real(kind=RP)                           :: epsil

        select case (self%epsilon_type)
        case (0)
! EPSILON - opcion 1 (se lee del fichero)
            epsil = self%turbine_t(1)%blade_t(1)%gauss_epsil
        case (1)
! EPSILON - opcion 2
            if (present(Cd)) then
                epsil = max(self%turbine_t(1)%blade_t(jj)%chord(ii)/4.0_RP,self%turbine_t(1)%blade_t(jj)%chord(ii)*Cd/2.0_RP)
            else
                epsil = self%turbine_t(1)%blade_t(1)%gauss_epsil
            end if
        case default
            epsil = self%turbine_t(1)%blade_t(1)%gauss_epsil
        end select

        GaussianInterpolation = exp( -(POW2(x(1) - self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,1)) + &
                  POW2(x(2) - self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,2)) + POW2(x(3) - self%turbine_t(1)%blade_t(jj)%point_xyz_loc(ii,3))) / POW2(epsil) ) / ( POW3(epsil) * pi**(3.0_RP/2.0_RP) )

    End Function GaussianInterpolation
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

function element_averageQ(mesh,eID)
   use HexMeshClass
   use PhysicsStorage
   implicit none

   type(HexMesh), intent(in)    :: mesh
   integer, intent(in)          :: eID 
   integer                      :: k, j, i

   integer                      :: total_points
   real(kind=RP), dimension(NCONS)   :: element_averageQ, Qsum

 
   Qsum(:) = 0.0_RP
   total_points = 0
   do k = 0, mesh%elements(eID) % Nxyz(3)   ; do j = 0, mesh%elements(eID) % Nxyz(2) ; do i = 0, mesh%elements(eID) % Nxyz(1)
       Qsum(:)=Qsum(:)+mesh%elements(eID) % Storage % Q(:,k,j,i)
       total_points=total_points + 1
   end do                  ; end do                ; end do

   element_averageQ(:) = Qsum(:) / real(total_points,RP)

end function element_averageQ

#endif
end module 

