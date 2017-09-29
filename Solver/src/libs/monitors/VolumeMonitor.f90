module VolumeMonitor
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
#include "Includes.h"


   private 
   public VOLUME, KINETIC_ENERGY, KINETIC_ENERGY_RATE, ENSTROPHY
   public VolumeMonitor_t

!
!  *******************************
!  Volume monitor class definition
!  *******************************
!
   type VolumeMonitor_t
      logical                         :: active
      integer                         :: ID
      real(kind=RP)                   :: values(BUFFER_SIZE)
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
      subroutine VolumeMonitor_Initialization( self , mesh , ID, solution_file )
!
!        *****************************************************************************
!              This subroutine initializes the volume monitor. The following
!           data is obtained from the case file:
!              -> Name: The monitor name (10 characters maximum)
!              -> Variable: The variable to be monitorized.
!        *****************************************************************************
!  
         use ParamfileRegions
         implicit none
         class(VolumeMonitor_t) :: self
         class(HexMesh)         :: mesh
         integer                :: ID
         character(len=*)       :: solution_file
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
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
         select case ( trim ( self % variable ) )
         
         case ("Kinetic energy")

         case ("Kinetic energy rate")

         case ("Enstrophy")

         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Kinetic energy"
               print*, "   * Kinetic energy rate"
               print*, "   * Enstrophy"
               stop "Stopped."

            end if
         end select

!
!        Prepare the file in which the monitor is exported
!        -------------------------------------------------
         write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % monitorName) , ".volume"  
!
!        Create file
!        -----------
         open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
!
!        Write the file headers
!        ----------------------
         write( fID , '(A20,A  )') "Monitor name:      ", trim(self % monitorName)
         write( fID , '(A20,A  )') "Selected variable: " , trim(self % variable)

         write( fID , * )
         write( fID , '(A10,2X,A24,2X,A24)' ) "Iteration" , "Time" , trim(self % variable)

         close ( fID ) 

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
         class   (  VolumeMonitor_t )  :: self
         class   (  HexMesh       ) :: mesh
         integer                       :: bufferPosition
!
!        Compute the volume integral
!        ---------------------------
         select case ( trim(self % variable) )
         case ("kinetic energy")
            self % values(bufferPosition) = ScalarVolumeIntegral(mesh, spA, KINETIC_ENERGY)

         case ("kinetic energy rate")
   
         case ("enstrophy")

         end select
         self % values(bufferPosition) = 1.0_RP

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

         open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
         
         do i = 1 , no_of_lines
            write( fID , '(I10,2X,ES24.16,2X,ES24.16)' ) iter(i) , t(i) , self % values(i)

         end do
        
         close ( fID )

         self % values(1) = self % values(no_of_lines)
      
      end subroutine VolumeMonitor_WriteToFile
end module VolumeMonitor
