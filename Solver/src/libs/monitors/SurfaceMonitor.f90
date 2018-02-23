#if defined(NAVIERSTOKES)
module SurfaceMonitorClass
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use MPI_Process_Info
#include "Includes.h"
   
   private
   public   SurfaceMonitor_t


!
!  ********************************
!  Surface monitor class definition
!  ********************************
!
   type SurfaceMonitor_t
      logical                         :: active
      logical                         :: isDimensionless
      integer                         :: ID
      real(kind=RP)                   :: direction(NDIM)
      integer                         :: marker
      real(kind=RP), allocatable      :: referenceSurface
      real(kind=RP), allocatable      :: values(:)
      real(kind=RP)                   :: dynamicPressure
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: variable
      contains
         procedure   :: Initialization => SurfaceMonitor_Initialization
         procedure   :: Update         => SurfaceMonitor_Update
         procedure   :: WriteLabel     => SurfaceMonitor_WriteLabel
         procedure   :: WriteValues    => SurfaceMonitor_WriteValue
         procedure   :: WriteToFile    => SurfaceMonitor_WriteToFile
   end type SurfaceMonitor_t

   contains
!
!/////////////////////////////////////////////////////////////////////////
!
!           SURFACE MONITOR PROCEDURES
!           --------------------------
!/////////////////////////////////////////////////////////////////////////
!
      subroutine SurfaceMonitor_Initialization( self , mesh , ID, solution_file , FirstCall)
!
!        *****************************************************************************
!              This subroutine initializes the surface monitor. The following
!           data is obtained from the case file:
!              -> Name: The monitor name (10 characters maximum)
!              -> Marker: The surface marker in which the monitor will be computed.
!              -> Variable: The variable to be monitorized.
!              -> Reference surface (optional): Reference surface for lift/drag coefficients
!              -> Direction (optional): Direction in which the forces are computed
!        *****************************************************************************
!  
         use ParamfileRegions
         implicit none
         class(SurfaceMonitor_t) :: self
         class(HexMesh)          :: mesh
         integer                 :: ID
         character(len=*)        :: solution_file
         logical, intent(in)     :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         character(len=STR_LEN_MONITORS)  :: directionName
         integer, allocatable             :: marker
         character(len=STR_LEN_MONITORS)  :: markerName
         integer                          :: pos
         integer                          :: fID
         integer                          :: zoneID
         real(kind=RP)                    :: directionValue(NDIM)
         interface
            function getArrayFromString( line ) result ( array )
               use SMConstants
               implicit none
               character(len=*),    intent(in)  :: line
               real(kind=RP), allocatable       :: array(:)
            end function getArrayFromString
         end interface

!
!        Get monitor ID
!        --------------
         self % ID = ID
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define surface monitor " , self % ID
         
         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "Name"              , self % monitorName      , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "Marker"            , markerName              , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "Variable"          , self % variable         , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "Reference surface" , self % referenceSurface , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "Direction"         , directionName        , in_label , "# end" ) 
!
!        Enable the monitor
!        ------------------
         self % active = .true.
         allocate ( self % values(BUFFER_SIZE) )
!
!        Get the surface marker
!        ----------------------
         self % marker = -1
         do zoneID = 1, size(mesh % zones)
            if ( trim(mesh % zones(zoneID) % name) .eq. trim(markerName) ) then
               self % marker = zoneID
               exit
            end if
         end do

         if ( self % marker .eq. -1 ) then
            self % active = .false.
            write(*,'(A,I0)') "Warning: Marker not specified for surface monitor ", self % ID
            write(*,'(A,I0,A)') "     Surface monitor ", self % ID, " disabled."
         end if
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
!
!        ****************************************
         select case ( trim ( self % variable ) )
!        ****************************************
!
            case ("mass-flow")
               self % isDimensionless = .false.

            case ("flow")
               self % isDimensionless = .false.

            case ("pressure-force")
               self % isDimensionless = .false.
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for pressure-force in surface monitor " , self % ID , "."
                  stop "Stopped"

               else
                  directionValue = getArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

            case ("viscous-force")
               self % isDimensionless = .false.
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for pressure-force in surface monitor " , self % ID , "."
                  stop "Stopped"

               else
                  directionValue = getArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

            case ("force")
               self % isDimensionless = .false.

               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for pressure-force in surface monitor " , self % ID , "."
                  stop "Stopped"

               else
                  directionValue = getArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if


            case ("lift")
               self % isDimensionless = .true.

               if ( .not. allocated ( self % referenceSurface ) ) then
                  print*, "Reference surface not specified for lift surface monitor " , self % ID , "."
                  stop "Stopped"
               end if
               
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for lift in surface monitor " , self % ID , "."
                  print*, "    ...  Using [0,1,0] as default."
                  self % direction = [0._RP,1._RP,0._RP]

               else
                  directionValue = getArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

               self % dynamicPressure = 0.5_RP * refValues % rho * POW2(refValues % V)* self % referenceSurface

            case ("drag")
               self % isDimensionless = .true.

               if ( .not. allocated ( self % referenceSurface ) ) then
                  print*, "Reference surface not specified for drag surface monitor " , self % ID , "."
                  stop "Stopped"
               end if
               
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for drag in surface monitor " , self % ID , "."
                  print*, "    ...  Using [1,0,0] as default."
                  self % direction = [1._RP,0._RP,0._RP]

               else
                  directionValue = getArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

               self % dynamicPressure = 0.5_RP * refValues % rho * refValues % V * refValues % V * self % referenceSurface

            case ("pressure-average")
               self % isDimensionless = .false.

            case default

               if ( len_trim (self % variable) .eq. 0 ) then
                  print*, "Variable was not specified for surface monitor " , self % ID , "."
               else
                  print*, 'Variable "',trim(self % variable),'" surface monitor ', self % ID, ' not implemented yet.'
                  print*, "Options available are:"
                  print*, "   * mass-flow"
                  print*, "   * flow"
                  print*, "   * pressure-force"
                  print*, "   * viscous-force"
                  print*, "   * force"
                  print*, "   * lift"
                  print*, "   * drag"
                  print*, "   * pressure-average"
                  stop "Stopped."

               end if
!
!        **********
         end select
!        **********
!
!        Prepare the file in which the monitor is exported
!        -------------------------------------------------
         write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % monitorName) , ".surface"  
!
!        Create file
!        -----------
         if (FirstCall) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
!
!        Write the file headers
!        ----------------------
            write( fID , '(A20,A  )') "Monitor name:      ", trim(self % monitorName)
            write( fID , '(A20,I0 )') "Surface marker:    ", self % marker
            write( fID , '(A20,A  )') "Selected variable: " , trim(self % variable)

            if ( self % isDimensionless ) then
               write(fID , '(A20,ES24.10)') "Dynamic pressure: " , self % dynamicPressure
            end if

            write( fID , * )
            write( fID , '(A10,2X,A24,2X,A24)' ) "Iteration" , "Time" , trim(self % variable)

            close ( fID )
         end if
      end subroutine SurfaceMonitor_Initialization

      subroutine SurfaceMonitor_Update ( self, mesh, bufferPosition )
!
!        *******************************************************************
!           This subroutine updates the monitor value computing it from
!           the mesh. It is stored in the "bufferPosition" position of the 
!           buffer.
!        *******************************************************************
!
         use SurfaceIntegrals
         implicit none
         class   (  SurfaceMonitor_t )   :: self
         class   (  HexMesh       )      :: mesh
         integer                         :: bufferPosition
         real(kind=RP)                   :: F(NDIM)
         real(kind=RP)                   :: direction(NDIM)

         select case ( trim ( self % variable ) )

         case ("mass-flow")
            self % values(bufferPosition) = ScalarSurfaceIntegral(mesh, self % marker, MASS_FLOW)

         case ("flow")
            self % values(bufferPosition) = ScalarSurfaceIntegral(mesh, self % marker, FLOW)

         case ("pressure-force")
            F = VectorSurfaceIntegral(mesh, self % marker, PRESSURE_FORCE)
            F = refValues % rho * POW2(refValues % V) * POW2(refValues % L) * F
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("viscous-force")
            F = VectorSurfaceIntegral(mesh, self % marker, VISCOUS_FORCE)
            F = refValues % rho * POW2(refValues % V) * POW2(refValues % L) * F
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("force")
            F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE)
            F = refValues % rho * POW2(refValues % V) * POW2(refValues % L) * F
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("lift")
            F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE)
            F = 2.0_RP * POW2(refValues % L) * F / self % referenceSurface
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("drag")

            F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE)
            F = 2.0_RP * POW2(refValues % L) * F / self % referenceSurface
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("pressure-average")
            self % values(bufferPosition) = ScalarSurfaceIntegral(mesh, self % marker, PRESSURE_FORCE) / ScalarSurfaceIntegral(mesh, self % marker, SURFACE)
  
         end select
         

      end subroutine SurfaceMonitor_Update

      subroutine SurfaceMonitor_WriteLabel ( self )
!
!        *************************************************************
!              This subroutine writes the label for the surface
!           monitor, when invoked from the time integrator Display
!           procedure.
!        *************************************************************
!
         implicit none
         class(SurfaceMonitor_t), intent(in)   :: self

         write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))

      end subroutine SurfaceMonitor_WriteLabel
   
      subroutine SurfaceMonitor_WriteValue ( self , bufferLine ) 
!
!        *************************************************************
!              This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
!
         implicit none
         class(SurfaceMonitor_t) :: self
         integer                 :: bufferLine

         write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values ( bufferLine ) 

      end subroutine SurfaceMonitor_WriteValue 

      subroutine SurfaceMonitor_WriteToFile ( self , iter , t , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(SurfaceMonitor_t) :: self
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
      
      end subroutine SurfaceMonitor_WriteToFile
!
end module SurfaceMonitorClass
#endif
