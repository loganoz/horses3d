!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
module MonitorsClass
   use SMConstants
   use NodalStorageClass
   use HexMeshClass
   use MonitorDefinitions
   use ResidualsMonitorClass
   use VolumeMonitorClass
#if defined(NAVIERSTOKES)
   use StatisticsMonitor
   use ProbeClass
   use SurfaceMonitorClass
#endif
   implicit none
!
#include "Includes.h"

   private
   public      Monitor_t , ConstructMonitors
!
!  *****************************
!  Main monitor class definition
!  *****************************
!  
   type Monitor_t
      character(len=LINE_LENGTH)           :: solution_file
      integer                              :: no_of_probes
      integer                              :: no_of_surfaceMonitors
      integer                              :: no_of_volumeMonitors
      integer                              :: bufferLine
      integer                , allocatable :: iter(:)
      integer                              :: dt_restriction
      logical                              :: write_dt_restriction
      real(kind=RP)          , allocatable :: t(:)
      real(kind=RP)          , allocatable :: SimuTime(:)
      type(Residuals_t)                    :: residuals
      class(VolumeMonitor_t),  allocatable :: volumeMonitors(:)
#if defined(NAVIERSTOKES)
      class(Probe_t),          allocatable :: probes(:)
      class(SurfaceMonitor_t), allocatable :: surfaceMonitors(:)
      type(StatisticsMonitor_t)            :: stats
#endif
      contains
         procedure   :: WriteLabel      => Monitor_WriteLabel
         procedure   :: WriteUnderlines => Monitor_WriteUnderlines
         procedure   :: WriteValues     => Monitor_WriteValues
         procedure   :: UpdateValues    => Monitor_UpdateValues
         procedure   :: WriteToFile     => Monitor_WriteToFile
   end type Monitor_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      function ConstructMonitors( mesh, controlVariables ) result(Monitors)
         use FTValueDictionaryClass
         use mainKeywordsModule
         implicit none
         class(HexMesh), intent(in)           :: mesh
         class(FTValueDictionary), intent(in) :: controlVariables
         type(Monitor_t)                      :: Monitors
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: fID , io
         integer                         :: i
         character(len=STR_LEN_MONITORS) :: line
         character(len=STR_LEN_MONITORS) :: solution_file                                            
         logical, save                   :: FirstCall = .TRUE.
         interface
            character(len=LINE_LENGTH) function getFileName(inputLine)
               use SMConstants
               character(len=*)   :: inputLine
            end function getFileName
         end interface
!
!        Setup the buffer
!        ----------------
!         if (controlVariables % containsKey("monitors flush interval") ) then
!            BUFFER_SIZE = controlVariables % integerValueForKey("monitors flush interval")
!         else
!            BUFFER_SIZE = BUFFER_SIZE_DEFAULT
!         end if
         allocate ( Monitors % SimuTime(BUFFER_SIZE), Monitors % t(BUFFER_SIZE), Monitors % iter(BUFFER_SIZE) )
!
!        Get the solution file name
!        --------------------------
         solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_MONITORS )
!
!        Remove the *.hsol termination
!        -----------------------------
         solution_file = trim(getFileName(solution_file))
         Monitors % solution_file = trim(solution_file)
!
!        Search in case file for probes, surface monitors, and volume monitors
!        ---------------------------------------------------------------------
         call getNoOfMonitors( Monitors % no_of_probes, Monitors % no_of_surfaceMonitors, Monitors % no_of_volumeMonitors )
!
!        Initialize
!        ----------
         call Monitors % residuals % Initialization( solution_file , FirstCall )

         allocate ( Monitors % volumeMonitors ( Monitors % no_of_volumeMonitors )  )
         do i = 1 , Monitors % no_of_volumeMonitors
            call Monitors % volumeMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall  )
         end do

#if defined(NAVIERSTOKES)
         call Monitors % stats     % Construct(mesh)

         allocate ( Monitors % probes ( Monitors % no_of_probes )  )
         do i = 1 , Monitors % no_of_probes
            call Monitors % probes(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do

         allocate ( Monitors % surfaceMonitors ( Monitors % no_of_surfaceMonitors )  )
         do i = 1 , Monitors % no_of_surfaceMonitors
            call Monitors % surfaceMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
#endif

         Monitors % write_dt_restriction = controlVariables % logicalValueForKey( "write dt restriction" )
         
         Monitors % bufferLine = 0
         
         FirstCall = .FALSE.
      end function ConstructMonitors

      subroutine Monitor_WriteLabel ( self )
!
!        ***************************************************
!           This subroutine prints the labels for the time
!         integrator Display procedure.
!        ***************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)              :: self
         integer                       :: i 
      
         if ( .not. MPI_Process % isRoot ) return
!
!        Write "Iteration" and "Time"
!        ----------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) "Iteration"
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "Time"
!
!        Write residuals labels
!        ----------------------
         call self % residuals % WriteLabel
!
!        Write volume monitors labels
!        -----------------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteLabel
         end do

#if defined(NAVIERSTOKES)
!
!        Write probes labels
!        -------------------
         do i = 1 , self % no_of_probes
            call self % probes(i) % WriteLabel
         end do
!
!        Write surface monitors labels
!        -----------------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteLabel
         end do

         call self % stats % WriteLabel
#endif
!
!        Write label for dt restriction
!        ------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "dt restr."

         write(STD_OUT , *) 

      end subroutine Monitor_WriteLabel

      subroutine Monitor_WriteUnderlines( self ) 
!
!        ********************************************************
!              This subroutine displays the underlines for the
!           time integrator Display procedure.
!        ********************************************************
!
         use PhysicsStorage
         use MPI_Process_Info
         implicit none
         class(Monitor_t)                         :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: i
         character(len=MONITOR_LENGTH), parameter :: dashes = "----------"

         if ( .not. MPI_Process % isRoot ) return
!
!        Print dashes for "Iteration" and "Time"
!        ---------------------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) trim ( dashes ) 
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
!
!        Print dashes for residuals
!        --------------------------
         do i = 1 , NTOTALVARS
            write(STD_OUT , '(3X,A10)' , advance = "no" ) trim(dashes)
         end do
!
!        Print dashes for volume monitors
!        --------------------------------
         do i = 1 , self % no_of_volumeMonitors
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % volumeMonitors(i) % monitorName ) + 2 ) )
         end do
#if defined(NAVIERSTOKES)
!
!        Print dashes for probes
!        -----------------------
         do i = 1 , self % no_of_probes
            if ( self % probes(i) % active ) then
               write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % probes(i) % monitorName ) + 2 ) )
            end if
         end do
!
!        Print dashes for surface monitors
!        ---------------------------------
         do i = 1 , self % no_of_surfaceMonitors
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % surfaceMonitors(i) % monitorName ) + 2 ) )
         end do

         if ( self % stats % state .ne. 0 ) write(STD_OUT,'(3X,A10)',advance="no") trim(dashes)
#endif

         
!
!        Print dashes for dt restriction
!        -------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
         
         write(STD_OUT , *) 

      end subroutine Monitor_WriteUnderlines

      subroutine Monitor_WriteValues ( self )
!
!        *******************************************************
!              This subroutine prints the values for the time
!           integrator Display procedure.
!        *******************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)           :: self
         integer                    :: i

         if ( .not. MPI_Process % isRoot ) return
!
!        Print iteration and time
!        ------------------------
         write ( STD_OUT , ' ( I10            ) ' , advance = "no" ) self % iter    ( self % bufferLine ) 
         write ( STD_OUT , ' ( 1X,A,1X,ES10.3 ) ' , advance = "no" ) "|" , self % t ( self % bufferLine ) 
!
!        Print residuals
!        ---------------
         call self % residuals % WriteValues( self % bufferLine )
!
!        Print volume monitors
!        ---------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteValues ( self % bufferLine )
         end do

#if defined(NAVIERSTOKES)
!
!        Print probes
!        ------------
         do i = 1 , self % no_of_probes
            call self % probes(i) % WriteValues ( self % bufferLine )
         end do
!
!        Print surface monitors
!        ----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteValues ( self % bufferLine )
         end do

         call self % stats % WriteValue
#endif
!
!        Print dt restriction
!        --------------------
         if (self % write_dt_restriction) then
            select case (self % dt_restriction)
               case (DT_FIXED) ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Fixed'
               case (DT_DIFF)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Diffusive'
               case (DT_CONV)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Convective'
            end select
         end if

         write(STD_OUT , *) 

      end subroutine Monitor_WriteValues

      subroutine Monitor_UpdateValues ( self, mesh, t , iter, maxResiduals )
!
!        ***************************************************************
!              This subroutine updates the values for the residuals,
!           for the probes, surface and volume monitors.
!        ***************************************************************
!        
         use PhysicsStorage
         use StopwatchClass
         implicit none
         class(Monitor_t)    :: self
         class(HexMesh)      :: mesh
         real(kind=RP)       :: t
         integer             :: iter
         real(kind=RP)       :: maxResiduals(NTOTALVARS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i 
!
!        Move to next buffer line
!        ------------------------
         self % bufferLine = self % bufferLine + 1
!
!        Save time, iteration and CPU-time
!        -----------------------
         self % t       ( self % bufferLine )  = t
         self % iter    ( self % bufferLine )  = iter
         self % SimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("Solver")
!
!        Compute current residuals
!        -------------------------
         call self % residuals % Update( mesh, maxResiduals, self % bufferLine )
!
!        Update volume monitors
!        ----------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % Update( mesh , self % bufferLine )
         end do

#if defined(NAVIERSTOKES)
!
!        Update probes
!        -------------
         do i = 1 , self % no_of_probes
            call self % probes(i) % Update( mesh , self % bufferLine )
         end do
!
!        Update surface monitors
!        -----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % Update( mesh , self % bufferLine )
         end do
!
!        Update statistics
!        -----------------
         call self % stats % Update(mesh, iter, t, trim(self % solution_file) )
#endif

!
!        Update dt restriction
!        ---------------------
         if (self % write_dt_restriction) self % dt_restriction = mesh % dt_restriction 
         
      end subroutine Monitor_UpdateValues

      subroutine Monitor_WriteToFile ( self , mesh, force) 
!
!        ******************************************************************
!              This routine has a double behaviour:
!           force = .true.  -> Writes to file and resets buffers
!           force = .false. -> Just writes to file if the buffer is full
!        ******************************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)        :: self
         class(HexMesh)          :: mesh
         logical, optional       :: force
!        ------------------------------------------------
         integer                 :: i 
         logical                 :: forceVal

         if ( present ( force ) ) then
            forceVal = force

         else
            forceVal = .false.

         end if

         if ( forceVal ) then 
!
!           In this case the monitors are exported to their files and the buffer is reseted
!           -------------------------------------------------------------------------------
            call self % residuals % WriteToFile ( self % iter , self % t, self % SimuTime , self % bufferLine )
   
            do i = 1 , self % no_of_volumeMonitors
               call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

#if defined(NAVIERSTOKES)   
            do i = 1 , self % no_of_probes
               call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
   
            do i = 1 , self % no_of_surfaceMonitors
               call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
!
!              Write statistics
!              ----------------
            if ( self % bufferLine .eq. 0 ) then
               i = 1
            else
               i = self % bufferLine
            end if
            call self % stats % WriteFile(mesh, self % iter(i), self % t(i), self % solution_file)
#endif
!
!           Reset buffer
!           ------------
            self % bufferLine = 0

         else
!
!           The monitors are exported just if the buffer is full
!           ----------------------------------------------------
            if ( self % bufferLine .eq. BUFFER_SIZE ) then

               call self % residuals % WriteToFile ( self % iter , self % t, self % SimuTime , BUFFER_SIZE )

               do i = 1 , self % no_of_volumeMonitors
                  call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

#if defined(NAVIERSTOKES)
               do i = 1 , self % no_of_probes
                  call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine ) 
               end do

               do i = 1 , self % no_of_surfaceMonitors
                  call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do
#endif
!
!              Reset buffer
!              ------------
               self % bufferLine = 0

            end if
         end if

      end subroutine Monitor_WriteToFile
!
!//////////////////////////////////////////////////////////////////////////////
!
!        Auxiliars
!
!//////////////////////////////////////////////////////////////////////////////
!
   subroutine getNoOfMonitors(no_of_probes, no_of_surfaceMonitors, no_of_volumeMonitors)
      use ParamfileRegions
      implicit none
      integer, intent(out)    :: no_of_probes
      integer, intent(out)    :: no_of_surfaceMonitors
      integer, intent(out)    :: no_of_volumeMonitors
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
!
!     Initialize
!     ----------
      no_of_probes = 0
      no_of_surfaceMonitors = 0
      no_of_volumeMonitors = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find monitors
!     ------------------------------------
readloop:do 
         read ( fID , '(A)' , iostat = io ) line

         if ( io .lt. 0 ) then
!
!           End of file
!           -----------
            line = ""
            exit readloop

         elseif ( io .gt. 0 ) then
!
!           Error
!           -----
            errorMessage(STD_OUT)
            stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineprobe' ) .gt. 0 ) then
               no_of_probes = no_of_probes + 1

            elseif ( index ( line , '#definesurfacemonitor' ) .gt. 0 ) then
               no_of_surfaceMonitors = no_of_surfaceMonitors + 1 

            elseif ( index ( line , '#definevolumemonitor' ) .gt. 0 ) then
               no_of_volumeMonitors = no_of_volumeMonitors + 1 

            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

end subroutine getNoOfMonitors

end module MonitorsClass
!
!///////////////////////////////////////////////////////////////////////////////////
!
