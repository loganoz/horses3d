module VolumeMonitorClass
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use MPI_Process_Info
#include "Includes.h"

   private 
   public VOLUME
#if defined(NAVIERSTOKES)
   public KINETIC_ENERGY, KINETIC_ENERGY_RATE, ENSTROPHY
   public ENTROPY, ENTROPY_RATE
#elif defined(CAHNHILLIARD)
   public FREE_ENERGY
#endif
   public VolumeMonitor_t

!
!  *******************************
!  Volume monitor class definition
!  *******************************
!
   type VolumeMonitor_t
      logical                         :: active
      integer                         :: ID
      real(kind=RP), allocatable      :: values(:)
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: variable
      contains
         procedure   :: Initialization => VolumeMonitor_Initialization
         procedure   :: Update         => VolumeMonitor_Update
         procedure   :: WriteLabel     => VolumeMonitor_WriteLabel
         procedure   :: WriteValues    => VolumeMonitor_WriteValue
         procedure   :: WriteToFile    => VolumeMonitor_WriteToFile
   end type VolumeMonitor_t
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////
!
!           VOLUME MONITOR PROCEDURES
!           -------------------------
!///////////////////////////////////////////////////////////////////////////
!
      subroutine VolumeMonitor_Initialization( self , mesh , ID, solution_file , FirstCall)
!
!        *****************************************************************************
!              This subroutine initializes the volume monitor. The following
!           data is obtained from the case file:
!              -> Name: The monitor name (10 characters maximum)
!              -> Variable: The variable to be monitorized.
!        *****************************************************************************
!  
         use ParamfileRegions
         use Utilities, only: toLower
         implicit none
         class(VolumeMonitor_t) :: self
         class(HexMesh)         :: mesh
         integer                :: ID
         character(len=*)       :: solution_file
         logical, intent(in)    :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         integer                          :: fID
         integer                          :: pos
!
!        Get monitor ID
!        --------------
         self % ID = ID
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define volume monitor " , self % ID
         
         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "Name"              , self % monitorName      , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "Variable"          , self % variable         , in_label , "# end" ) 
!
!        Enable the monitor
!        ------------------
         self % active = .true.
         allocate ( self % values(BUFFER_SIZE) )
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
         call toLower(self % variable)

#if defined(NAVIERSTOKES)
         select case ( trim ( self % variable ) )
         case ("kinetic energy")
         case ("kinetic energy rate")
         case ("enstrophy")
         case ("entropy")
         case ("entropy rate")
         case ("mean velocity")
         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Kinetic energy"
               print*, "   * Kinetic energy rate"
               print*, "   * Enstrophy"
               print*, "   * Entropy"
               print*, "   * Entropy rate"
               print*, "   * Mean velocity"
               stop "Stopped."

            end if
         end select

#elif defined(CAHNHILLIARD)
         select case ( trim ( self % variable ) )
         case ("free energy")
         case default
            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Free energy"
               stop "Stopped."

            end if
         end select
#endif
!
!        Prepare the file in which the monitor is exported
!        -------------------------------------------------
         write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % monitorName) , ".volume"  
!
!        Create file
!        -----------
         if (FirstCall) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
!
!        Write the file headers
!        ----------------------
            write( fID , '(A20,A  )') "Monitor name:      ", trim(self % monitorName)
            write( fID , '(A20,A  )') "Selected variable: " , trim(self % variable)

            write( fID , * )
            write( fID , '(A10,2X,A24,2X,A24)' ) "Iteration" , "Time" , trim(self % variable)

            close ( fID )
         end if

      end subroutine VolumeMonitor_Initialization

      subroutine VolumeMonitor_Update ( self , mesh , bufferPosition )
!
!        *******************************************************************
!           This subroutine updates the monitor value computing it from
!           the mesh. It is stored in the "bufferPosition" position of the 
!           buffer.
!        *******************************************************************
!
         use VolumeIntegrals
         implicit none
         class   (  VolumeMonitor_t ) :: self
         class   (  HexMesh       )   :: mesh
         integer                      :: bufferPosition
!
!        Compute the volume integral
!        ---------------------------
         select case ( trim(self % variable) )
#if defined(NAVIERSTOKES)
         case ("kinetic energy")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, KINETIC_ENERGY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("kinetic energy rate")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, KINETIC_ENERGY_RATE) / ScalarVolumeIntegral(mesh, VOLUME)
   
         case ("enstrophy")
            self % values(bufferPosition) = 0.5_RP * ScalarVolumeIntegral(mesh, ENSTROPHY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("entropy")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("entropy rate")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY_RATE) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("mean velocity")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, VELOCITY) / ScalarVolumeIntegral(mesh, VOLUME)

#elif defined(CAHNHILLIARD)
         case ("free energy")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, FREE_ENERGY) / ScalarVolumeIntegral(mesh, VOLUME)
            
#endif
         end select


      end subroutine VolumeMonitor_Update

      subroutine VolumeMonitor_WriteLabel ( self )
!
!        *************************************************************
!              This subroutine writes the label for the volume
!           monitor, when invoked from the time integrator Display
!           procedure.
!        *************************************************************
!
         implicit none
         class(VolumeMonitor_t)             :: self

         write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))

      end subroutine VolumeMonitor_WriteLabel
   
      subroutine VolumeMonitor_WriteValue ( self , bufferLine ) 
!
!        *************************************************************
!              This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
!
         implicit none
         class(VolumeMonitor_t) :: self
         integer                 :: bufferLine

         write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values ( bufferLine ) 

      end subroutine VolumeMonitor_WriteValue 

      subroutine VolumeMonitor_WriteToFile ( self , iter , t , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(VolumeMonitor_t) :: self
         integer                 :: iter(:)
         real(kind=RP)           :: t(:)
         integer                 :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i
         integer                    :: fID

         if ( MPI_Process % isRoot ) then
            open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
         
            do i = 1 , no_of_lines
               write( fID , '(I10,2X,ES24.16,2X,ES24.16)' ) iter(i) , t(i) , self % values(i)

            end do
        
            close ( fID )
         end if

         if ( no_of_lines .ne. 0 ) self % values(1) = self % values(no_of_lines)
      
      end subroutine VolumeMonitor_WriteToFile
end module VolumeMonitorClass
