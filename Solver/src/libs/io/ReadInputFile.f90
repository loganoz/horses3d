!
!////////////////////////////////////////////////////////////////////////
!
!      ReadInputFile.f90
!      Created: June 10, 2015 at 3:09 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
   module FileReaders
      use FileReadingUtilities, only: GetKeyword, GetValueAsString
      implicit none
      
      private
      public ReadControlFile, ReadOrderFile
      
   contains
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------
!     Control file reader
!     -------------------
      SUBROUTINE ReadControlFile (controlVariables)
         USE SMConstants
         USE FTValueDictionaryClass
         USE SharedBCModule
         USE mainKeywordsModule
         use MPI_Process_Info 
         use Utilities, only: toLower
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         CHARACTER(LEN=LINE_LENGTH) :: inputLine
         CHARACTER(LEN=LINE_LENGTH) :: keyword, keywordValue
         CHARACTER(LEN=LINE_LENGTH) :: boundaryName
         CHARACTER(LEN=LINE_LENGTH) :: boundaryType
         CHARACTER(LEN=LINE_LENGTH) :: boundaryValue
         CHARACTER(LEN=LINE_LENGTH) :: arg
         character(len=LINE_LENGTH) :: boundaryNameControlVariable
         INTEGER                    :: numberOfBCs, k
         INTEGER                    :: ist
         logical                                 :: isInsideHagstagZone
!
!        ---------------------------------------
!        External functions from FileReading.f90
!        ---------------------------------------
!
         integer                                 :: fid, io
         interface
            subroutine PreprocessInputLine(line)
               implicit none
               character(len=*), intent(inout) :: line
            end subroutine PreprocessInputLine
            subroutine WriteDefaultControlFile
            end subroutine WriteDefaultControlfile
         end interface
!
!        -----------------------------------------------
!        Read the input file.
!
!        we use dictionaries to store the input file 
!        parameters.
!        -----------------------------------------------
!

         if ( command_argument_count() .eq. 0 ) then
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/,/,A)') "*** ERROR: Missing input file"
               write(STD_OUT,'(A)') "*** ERROR: Syntax is: "
               write(STD_OUT,'(A)') "*** ERROR:             >> HORSES3D ControlFile.control"
               write(STD_OUT,'(A)') "*** ERROR: Default control file written to ./Sample.control"
               call WriteDefaultControlFile
               write(STD_OUT,'(A,/,/)') "*** ERROR: Stopping."
            end if
            stop
         end if
!
!        -----------------
!        Open control file
!        -----------------
!
         CALL get_command_argument(1, arg)
         OPEN(newunit=fid,file=trim(arg),status="old",action="read", iostat=io)

         if ( io .gt. 0 ) then
            write(STD_OUT,'(/,/,A,A,A)') '*** ERROR: Wrong control file "', trim(arg),'".'
            write(STD_OUT,'(A)')         '*** ERROR: Stopping.'
            stop
         end if

         isInsideHagstagZone = .false.

         DO
            READ(fid,'(A132)', IOSTAT = ist) inputLine

            IF(ist /= 0 ) EXIT !.OR. inputLine(1:1) == '/'

            call PreprocessInputLine(inputLine)
            if ( len_trim(inputLine) .le. 0 ) cycle

            if ( index(inputLine,'#define') .ne. 0 ) then
               isInsideHagstagZone = .true.

            elseif ( (index(inputLine,'#end') .ne. 0) .and. (isInsideHagstagZone) ) then
               isInsideHagstagZone = .false.

            end if

            if ( isInsideHagstagZone ) cycle
            
            keyword      = ADJUSTL(GetKeyword(inputLine))
            keywordValue = ADJUSTL(GetValueAsString(inputLine))
            CALL toLower(keyword)
            CALL controlVariables % addValueForKey(keywordValue,TRIM(keyword))
            
            IF(keyword == numberOfBoundariesKey) THEN 
!
!              ---------------------------------------------------------------------------
!              We will store the type and values of the boundaries in dictionaries so that
!              we can associate a name of a boundary curve found in the mesh file with a
!              particular value and type of boundary conditions.
!              ---------------------------------------------------------------------------
!
               numberOfBCs = controlVariables%integerValueForKey(numberOfBoundariesKey)
               
               DO k = 1, numberOfBCs 
                  READ(fid,*) boundaryName, boundaryValue, boundaryType
                  CALL toLower(boundaryName)
                  CALL toLower(boundaryType)
                  CALL bcTypeDictionary % addValueForKey(boundaryType, boundaryName)
                  CALL bcValueDictionary % addValueForKey(boundaryValue, boundaryName)
                  write(boundaryNameControlVariable,'(A,I0)') "BoundaryName",k
                  call controlVariables % addValueForKey(boundaryName,trim(boundaryNameControlVariable))
               END DO
            END IF
            
            IF(keyword == "adaptation conforming boundaries") THEN
!
!              ---------------------------------------------------------------------------
!              If there's p-adaptation the user may want to keep a conforming representation 
!              on certain boundaries.
!              TODO: Move this elsewhere after implementing the adaptator definition
!                    as a block
!              ---------------------------------------------------------------------------   
!
               numberOfBCs = controlVariables%integerValueForKey("adaptation conforming boundaries")
               
               DO k = 1, numberOfBCs 
                  READ(fid,*) boundaryName
                  
                  CALL toLower(boundaryName)
                  CALL conformingBoundariesDic % addValueForKey(boundaryName, boundaryName)
                  
               END DO
            END IF
            
         END DO

         CLOSE(UNIT=fid)

      END SUBROUTINE ReadControlFile
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------------
!     Subroutine that reads input file containing polynomial orders for mesh
!     ----------------------------------------------------------------------
      subroutine ReadOrderFile(filename, Nx, Ny, Nz)
         implicit none
         !-arguments----------------------------------------
         character(len=*), intent(in) :: filename          !<  Name of file containing polynomial orders to initialize
         integer, allocatable         :: Nx(:),Ny(:),Nz(:) !>  Polynomial orders for each element
         !-local-variables----------------------------------
         integer                      :: fd       ! File unit
         integer                      :: nelem    ! Number of elements
         integer                      :: i        ! counter
         !--------------------------------------------------
         
         open(newunit = fd, FILE = filename )   
            READ(fd,*) nelem
            
            allocate(Nx(nelem),Ny(nelem),Nz(nelem))
            
            do i = 1, nelem
               READ(fd,*) Nx(i), Ny(i), Nz(i)
            ENDDO
         close(UNIT=fd)
         
      end subroutine ReadOrderFile
   end module FileReaders
      
      subroutine PreprocessInputLine(line)
!
!        ******************************************************************
!        This function eliminates all text at the RHS of a comment (! or /)
!        ******************************************************************
!
         implicit none
         character(len=*), intent(inout)  :: line
!
!        ---------------
!        Local variables
!        ---------------
!
         character, parameter :: comments(1) = (/"!"/)
         integer              :: pos, com

         do com = 1, size(comments)
            pos = index(line, comments(com))

            if ( pos .gt. 0 ) then
               line = line(1:pos-1)
            end if
         end do

         IF ( line(1:1) == '/') line = ""

      end subroutine PreprocessInputLine

      subroutine WriteDefaultControlFile
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fid

         open(newunit=fid, file="./Sample.control", action="write", status="unknown")

         write(fid,'(A)') "!"
         write(fid,'(A)') "!       *******************"
         write(fid,'(A)') "!       Sample control file"
         write(fid,'(A)') "!       *******************"
         write(fid,'(A)') "!"
         write(fid,'(A)') "!-------------------------- Configuration:-"

         write(fid,'(A40,A,A)') "Mesh file name", " = ", "MESH_FILE_NAME"
         write(fid,'(A40,A,A)') "Solution file name", " = ", "SOLUTION_FILE_NAME.hsol"
         write(fid,'(A40,A,A)') "Save gradients to solution", " = ", ".false."
         write(fid,'(A40,A,A)') "2D mesh offset direction", " = ", "3D"
         write(fid,'(A40,A,A)') "Restart", " = ", ".false."
         write(fid,'(A40,A,A)') "Restart file name", " = ", "RESTART_FILE_NAME.hsol"


         write(fid,'(/,A)') "!-------------------- Physical parameters:-"
         write(fid,'(A)') "!                        ** Navier-Stokes"
         write(fid,'(A40,A,A)') "Flow equations", " = ", "NS  ! Euler"
         write(fid,'(A40,A,A)') "Mach number", " = ", "0.1"
         write(fid,'(A40,A,A)') "Reynolds number", " = ", "1000.0"
         write(fid,'(A40,A,A)') "Froude number", " = ", "1.0d+300"
         write(fid,'(A40,A,A)') "Prandtl number", " = ", "0.72"
         write(fid,'(A40,A,A)') "AoA Theta", " = ", "0.0"
         write(fid,'(A40,A,A)') "AoA Phi", " = ", "0.0"
         write(fid,'(A40,A,A)') "Gravity direction", " = ", "[x,y,z]"
         write(fid,'(A40,A,A)') "Compute gradients", " = ", ".true."

         write(fid,'(/,A)') "!                        ** Cahn-Hilliard"
         write(fid,'(A40,A,A)') "Peclet number", " = ", "1.0"
         write(fid,'(A40,A,A)') "Capilar number", " = ", "1.0"
         write(fid,'(A40,A,A)') "Interface width (dimensionless)", " = ", "1.0"
         write(fid,'(A40,A,A)') "Density ratio (rho2/rho1)", " = ", "1.0"
         write(fid,'(A40,A,A)') "Viscosity ratio (mu2/mu1)", " = ", "1.0"
         write(fid,'(A40,A,A)') "Wall contact angle", " = ", "0.0"

         write(fid,'(/,A)') "!------------------------- Discretization:-"
         write(fid,'(A40,A,A)') "Polynomial order", " = ", "5"
         write(fid,'(A40,A,A)') "Discretization nodes", " = ", "Gauss  ! Gauss-Lobatto"
         write(fid,'(A40,A,A)') "Riemann solver", " = ", "Roe  ! Standard Roe/Low dissipation Roe"
         write(fid,'(A40,A,A)') "Lambda stabilization", " = ", "1.0"
         write(fid,'(A40,A,A)') "Inviscid discretization", " = ", "Standard  ! Split-form"
         write(fid,'(A40,A,A)') "Split form", " = ", "Ducros   ! Pirozzoli/Morinishi/Kennedy-Gruber/Entropy conserving/Entropy and energy conserving"
         write(fid,'(A40,A,A)') "Viscous discretization", " = ", "BR1   ! IP/BR2"
         write(fid,'(A40,A,A)') "Cahn-Hilliard discretization", " = ", "IP   ! BR1/BR2"
         write(fid,'(A40,A,A)') "Penalty parameter", " = ", "1.0"
         write(fid,'(A40,A,A)') "Interior penalty variant", " = ", "SIPG   ! IIPG/NIPG"


         write(fid,'(/,A)') "!----------------------- Time integration:-"
         write(fid,'(A40,A,A)') "Time integration", " = ", "Explicit   ! IMEX/Implicit"
         write(fid,'(A40,A,A)') "CFL", " = ", "0.4"
         write(fid,'(A40,A,A)') "dCFL", " = ", "0.4"
         write(fid,'(A40,A,A)') "dt", " = ", "0.01"
         write(fid,'(A40,A,A)') "Number of time steps", " = ", "100"
         write(fid,'(A40,A,A)') "Output interval", " = ", "1"
         write(fid,'(A40,A,A)') "Convergence tolerance", " = ", "1.0e-10"
         write(fid,'(A40,A,A)') "Simulation type", " = ", "steady-state   ! time-accurate"
         write(fid,'(A40,A,A)') "Final time", " = ", "1.0"
         write(fid,'(A40,A,A)') "Autosave mode", " = ", "Iteration  ! Time"
         write(fid,'(A40,A,A)') "Autosave interval", " = ", "100   ! 1.0"


         write(fid,'(/,A)') "!-------------------- Boundary conditions:-"
         write(fid,'(A40,A,A)') "Number of boundaries", " = ", "8"
         write(fid,'(A,1X,A,1X,A)') "bname1", "0.0", "freeslipwall"
         write(fid,'(A,1X,A,1X,A)') "bname2", "0.0", "noslipadiabaticwall"
         write(fid,'(A,1X,A,1X,A)') "bname3", "0.0", "noslipisothermalwall"
         write(fid,'(A,1X,A,1X,A)') "bname4", "0.0", "inflow"
         write(fid,'(A,1X,A,1X,A)') "bname5", "0.0", "outflow"
         write(fid,'(A,1X,A,1X,A)') "bname6", "0.0", "periodic+"
         write(fid,'(A,1X,A,1X,A)') "bname7", "0.0", "periodic-"
         write(fid,'(A,1X,A,1X,A)') "bname8", "0.0", "user-defined"

         write(fid,'(/,/,/,A,A)') "!", "#define probe 1"
         write(fid,'(A,A20,A,A)') "!", "Name", " = ", "Probe example"
         write(fid,'(A,A20,A,A)') "!", "Position", " = ", "[x,y,z]"
         write(fid,'(A,A20,A,A)') "!", "Variable", " = ", "pressure/velocity/u/v/w/mach/K"  
         write(fid,'(A,A)') "!", "#end"

         write(fid,'(/,A,A)') "!", "#define surface monitor 1"
         write(fid,'(A,A20,A,A)') "!", "Name", " = ", "Surface monitor example"
         write(fid,'(A,A20,A,A)') "!", "Marker", " = ", "bnameX"
         write(fid,'(A,A20,A,A)') "!", "Variable", " = ", "mass-flow/flow/pressure-force/viscous-force/force/lift/drag/pressure-average"
         write(fid,'(A,A20,A,A)') "!", "Direction", " = ", "[x,y,z]"
         write(fid,'(A,A20,A,A)') "!", "Reference surface", " = ", "1.0"
         write(fid,'(A,A)') "!", "#end"

         write(fid,'(/,A,A)') "!", "#define volume monitor 1"
         write(fid,'(A,A20,A,A)') "!", "Name", " = ", "Surface monitor example"
         write(fid,'(A,A20,A,A)') "!", "Variable", " = ", "kinetic energy/kinetic energy rate/enstrophy/entropy/entropy rate/mean velocity"
         write(fid,'(A,A)') "!", "#end"


         write(fid,'(/,A,A)') "!", "#define statistics"
         write(fid,'(A,A20,A,A)') "!", "Variable", " = ", "kinetic energy/kinetic energy rate/enstrophy/entropy/entropy rate/mean velocity"
         write(fid,'(A,A20,A,A)') "!", "Sampling interval", " = ", "10"
         write(fid,'(A,A20,A,A)') "!", "Starting iteration", " = ", "0"
         write(fid,'(A,A20,A,A)') "!", "Starting time", " = ", "0.0"
         write(fid,'(A,A)') "!", "!Real time control commands:"
         write(fid,'(A,A)') "!", "@start"
         write(fid,'(A,A)') "!", "@stop"
         write(fid,'(A,A)') "!", "@pause"
         write(fid,'(A,A)') "!", "@reset"
         write(fid,'(A,A)') "!", "@dump"
         write(fid,'(A,A)') "!", "#end"

      end subroutine WriteDefaultControlFile
