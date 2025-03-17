!
!///////////////////////////////////////////////////////////////////////////
!
!     This class is a stopwatch to measure computation times. It handles 
!  multiple simultaneous events, and allows for pausing and resetting each
!  one individually. It measures both the elapsed time, and the total
!  CPU time (different when computing in parallel).
!
!  The user just interacts with the Stopwatch class. The process to measure
!  times is:
!
!     1/ Create an event. Each event is identified by a name (character)
!           call Stopwatch % CreateNewEvent("EventName")
!
!     2/ Start the event:
!           call Stopwatch % Start("EventName")
!
!     3/ Pause the event:
!           call Stopwatch % Pause("EventName")
!
!     4/ Get the times:
!           call Stopwatch % ElapsedTime("EventName")
!           call Stopwatch % CPUTime("EventName")
!
!     5/ Reset the event (if desired):
!           call Stopwatch % Reset("EventName")
!
!  NOTE: The time can be measured while the event is running. It just takes
!  a measure of the actual time and compares it to the previous register.
!
!///////////////////////////////////////////////////////////////////////////
!
module StopwatchClass
   use SMConstants
   implicit none

   private
   public   Stopwatch

   integer, parameter   :: IYEAR = 1
   integer, parameter   :: IMONTH = 2
   integer, parameter   :: IDAY = 3
   integer, parameter   :: IHOUR = 5
   integer, parameter   :: IMINUTES = 6
   integer, parameter   :: ISECONDS = 7
   integer, parameter   :: IMILLISECONDS = 8

   integer, parameter      :: STR_LEN_STOPWATCH = 128
   integer, parameter      :: reference_date(8) = [2017,6,1,0,0,0,0,0]

   type Event_t
      character(len=STR_LEN_STOPWATCH)    :: name
      logical                             :: running
      real(kind=RP)                       :: tstart
      real(kind=RP)                       :: tstartCPU
      real(kind=RP)                       :: elapsedTime
      real(kind=RP)                       :: lastElapsedTime
      real(kind=RP)                       :: CPUTime
      class(Event_t), pointer, private    :: next  => NULL()
      contains
         procedure      :: Start                => Event_Start
         procedure      :: Pause                => Event_Pause
         procedure      :: Reset                => Event_Reset
         procedure      :: GetElapsedTime       => Event_GetElapsedTime
         procedure      :: GetLastElapsedTime   => Event_GetLastElapsedTime
         procedure      :: GetCPUTime           => Event_GetCPUTime
         procedure      :: Destruct             => Event_Destruct
   end type Event_t

   type Stopwatch_t
      class(Event_t),   pointer, private  :: head => NULL()
      integer, private                    :: no_of_events = 0
      contains
         procedure   :: CreateNewEvent    => Stopwatch_CreateNewEvent
         procedure   :: Start             => Stopwatch_Start
         procedure   :: Pause             => Stopwatch_Pause
         procedure   :: Reset             => Stopwatch_Reset
         procedure   :: ElapsedTime       => Stopwatch_ElapsedTime
         procedure   :: lastElapsedTime   => Stopwatch_lastElapsedTime
         procedure   :: CPUTime           => Stopwatch_CPUTime
         procedure   :: Destruct          => Stopwatch_Destruct
         procedure   :: WriteSummaryFile  => Stopwatch_WriteSummaryFile
   end type Stopwatch_t

   type(Stopwatch_t)    :: Stopwatch 

   contains
!
!/////////////////////////////////////////////////////////////////////
!
!           STOPWATCH PROCEDURES
!           --------------------
!/////////////////////////////////////////////////////////////////////
!
   function NewStopwatch()
      implicit none
      type(Stopwatch_t)       :: NewStopwatch

      NewStopwatch % head => NULL()
      NewStopwatch % no_of_events = 0

   end function NewStopwatch

   subroutine Stopwatch_CreateNewEvent( self , Name )
      implicit none
      class(Stopwatch_t)      :: self
      character(len=*)    :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                 :: i
      class(Event_t), pointer :: event

      if ( self % no_of_events .eq. 0 ) then
         self % no_of_events = 1 
         allocate(self % head)

         select type (ev => self % head)
            type is (Event_t)
               ev = NewEvent(trim(Name))
         end select

      else

         event => self % head
         do i = 1 , self % no_of_events-1
            if ( trim(event % name) .eq. trim(name) ) then
               print*, "Warning: The event '", trim(Name), "' already exists."
               return
            end if
    
            event => event % next

         end do

         if ( trim(event % name) .eq. trim(name) ) then
            print*, "Warning: The event already exists."
         end if
        
         self % no_of_events = self % no_of_events + 1 
         allocate(event % next)
         event => event % next
   
         select type (ev => event)
            type is (Event_t)
               ev = NewEvent(trim(Name))
         end select

      end if
      
   end subroutine Stopwatch_CreateNewEvent

   subroutine Stopwatch_Start(self , Name)
      implicit none
      class(Stopwatch_t)               :: self
      character(len=*) :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 

      event => self % head
      do i = 1 , self % no_of_events

         if ( trim(event % Name) .eq. trim(Name) ) then
            call event % Start()
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

   end subroutine Stopwatch_Start

   subroutine Stopwatch_Pause(self , Name)
      implicit none
      class(Stopwatch_t)               :: self
      character(len=*) :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd

      event => self % head
      do i = 1 , self % no_of_events

         if ( trim(event % Name) .eq. trim(Name) ) then
            call event % Pause
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

   end subroutine Stopwatch_Pause

   function Stopwatch_ElapsedTime(self,Name,units) result(elapsedTime)
      implicit none
      class(Stopwatch_t)                     :: self
      character(len=*)                       :: Name
      character(len=*), intent(in), optional :: units
      real(kind=RP)                          :: elapsedTime
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd
      character(len=STR_LEN_STOPWATCH)       :: chosenUnits

      if ( present(units) ) then
         chosenUnits = units

      else
         chosenUnits = "sec"

      end if

      event => self % head
      do i = 1 , self % no_of_events 

         if ( trim(event % Name) .eq. trim(Name) ) then
            elapsedTime = event % GetElapsedTime(chosenUnits)
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

      elapsedTime = -1.0_RP

   end function Stopwatch_ElapsedTime
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------
!  Stopwatch_lastElapsedTime:
!  Gets last elapsed time (time between last start and pause -or now)
!  ------------------------------------------------------------------
   function Stopwatch_lastElapsedTime(self,Name,units) result(lastElapsedTime)
      implicit none
      !-arguments---------------------------------------------
      class(Stopwatch_t)                     :: self
      character(len=*)                       :: Name
      character(len=*), intent(in), optional :: units
      real(kind=RP)                          :: lastElapsedTime
      !-local-variables---------------------------------------
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd
      character(len=STR_LEN_STOPWATCH)       :: chosenUnits
      !-------------------------------------------------------
      
      if ( present(units) ) then
         chosenUnits = units

      else
         chosenUnits = "sec"

      end if

      event => self % head
      do i = 1 , self % no_of_events 

         if ( trim(event % Name) .eq. trim(Name) ) then
            lastElapsedTime = event % GetLastElapsedTime(chosenUnits)
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

      lastElapsedTime = -1.0_RP

   end function Stopwatch_lastElapsedTime
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Stopwatch_CPUTime(self,Name,units) result(CPUTime)
      implicit none
      class(Stopwatch_t)                     :: self
      character(len=*)                       :: Name
      character(len=*), intent(in), optional :: units
      real(kind=RP)                          :: CPUTime
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd
      character(len=STR_LEN_STOPWATCH)       :: chosenUnits

      if ( present(units) ) then
         chosenUnits = units

      else
         chosenUnits = "sec"

      end if

      event => self % head
      do i = 1 , self % no_of_events 

         if ( trim(event % Name) .eq. trim(Name) ) then
            CPUTime = event % GetCPUTime(chosenUnits)
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

      CPUTime = -1.0_RP

   end function Stopwatch_CPUTime

   subroutine Stopwatch_Reset(self,Name) 
      implicit none
      class(Stopwatch_t)               :: self
      character(len=*) :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd

      event => self % head
      do i = 1 , self % no_of_events

         if ( trim(event % Name) .eq. trim(Name) ) then
            call event % Reset
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

   end subroutine Stopwatch_Reset

   subroutine Stopwatch_Destruct(self)
      implicit none
      class(Stopwatch_t)      :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i
      class(Event_t), pointer    :: event , nextEvent

      event => self % head
      do i = 1 , self % no_of_events 
         nextEvent => event % next
         call event % Destruct  
         deallocate(event)
         event => nextEvent

      end do

      self % no_of_events = 0
      self % head => NULL()

   end subroutine Stopwatch_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------
!  Stopwatch_lastElapsedTime:
!  Gets last elapsed time (time between last start and pause -or now)
!  ------------------------------------------------------------------
   subroutine Stopwatch_WriteSummaryFile (this, fileRoot,units)
      implicit none
      !-arguments---------------------------------------------
      class(Stopwatch_t)         , intent(inout) :: this
      character(len=*)           , intent(in)    :: fileRoot
      character(len=*), optional , intent(in)    :: units
      !-local-variables---------------------------------------
      class(Event_t), pointer          :: event
      character(len=LINE_LENGTH)       :: fileName
      character(len=STR_LEN_STOPWATCH) :: chosenUnits
      integer                          :: fd, i
      !-------------------------------------------------------
      
      if ( present(units) ) then
         chosenUnits = units
      else
         chosenUnits = "sec"
      end if
      
      fileName = trim(fileRoot) // ".Stopwatch.info"
      
      open (newunit=fd, file = fileName)
      
      write(fd,'(A)') '#Stopwatch information file'
      write(fd,'(A7,43x,2(x,A15))') '# Event', 'Elapsed time(s)', 'CPU-Time (s)'
      
      event => this % head
      do i = 1 , this % no_of_events 
         write(fd,'(A50,2(x,ES15.5))') event % name, event % GetElapsedTime(chosenUnits), event % GetCPUTime(chosenUnits)
         event => event % next
      end do
      
      close(fd)
      
   end subroutine Stopwatch_WriteSummaryFile
!
!/////////////////////////////////////////////////////////////////////
!
!           EVENT PROCEDURES
!           ----------------
!/////////////////////////////////////////////////////////////////////
!
   function NewEvent(Name)
      implicit none
      character(len=*), intent(in) :: Name
      type(Event_t)                :: NewEvent
      
      NewEvent % Name        = trim(Name)
      NewEvent % running     = .false.
      NewEvent % tstart      = 0.0_RP
      NewEvent % tstartCPU   = 0.0_RP
      NewEvent % elapsedTime = 0.0_RP
      NewEvent % CPUTime     = 0.0_RP
      NewEvent % next        => NULL()

   end function NewEvent

   subroutine Event_Start(self)
      implicit none
      class(Event_t)    :: self

      self % running = .true.
      self % tstart = GetCurrentSeconds()
      call CPU_TIME( self % tstartCPU )

   end subroutine Event_Start

   subroutine Event_Pause(self)
      implicit none
      class(Event_t)    :: self
      real(kind=RP)     :: tEnd 
      real(kind=RP)     :: tCPUEnd

      if ( .not. self % running ) return

      tEnd = GetCurrentSeconds()
      call CPU_TIME(tCPUEnd)

      self % elapsedTime     = tEnd - self % tStart + self % elapsedTime
      self % lastElapsedTime = tEnd - self % tStart
      self % cpuTime         = tCPUEnd - self % tStartCPU + self % CPUTime

      self % tStart = 0.0_RP
      self % tStartCPU = 0.0_RP
      self % running = .false.

   end subroutine Event_Pause

   function Event_GetElapsedTime(self,units) result (time)
      implicit none
      class(Event_t)               :: self
      character(len=*), intent(in) :: units
      real(kind=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP)     :: tActual

      if ( self % running ) then
         tActual = GetCurrentSeconds()
         time = tActual - self % tStart + self % elapsedTime

      else
         time = self % elapsedTime

      end if

      if ( (trim(units) .eq. "minutes") .or. (trim(units) .eq. "Minutes") .or. (trim(units) .eq. "min") ) then
         time = time / 60.0_RP
   
      else if ( (trim(units) .eq. "hours") .or. (trim(units) .eq. "Hours") .or. (trim(units) .eq. "h") ) then
         time = time / 3600.0_RP

      else if ( (trim(units) .eq. "days") .or. (trim(units) .eq. "Days") .or. (trim(units) .eq. "d") ) then
         time = time / (24.0_RP * 3600.0_RP)
         
      end if

   end function Event_GetElapsedTime
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------
!  Event_GetLastElapsedTime:
!  Gets last elapsed time (time between last start and pause -or now)
!  ------------------------------------------------------------------
   function Event_GetLastElapsedTime(self,units) result (time)
      implicit none
      !-arguments---------------------------------------------
      class(Event_t)               :: self
      character(len=*), intent(in) :: units
      real(kind=RP)                :: time
      !-local-variables---------------------------------------
      real(kind=RP)     :: tActual
      !-------------------------------------------------------

      if ( self % running ) then
         tActual = GetCurrentSeconds()
         time = tActual - self % tStart

      else
         time = self % lastElapsedTime

      end if

      if ( (trim(units) .eq. "minutes") .or. (trim(units) .eq. "Minutes") .or. (trim(units) .eq. "min") ) then
         time = time / 60.0_RP
   
      else if ( (trim(units) .eq. "hours") .or. (trim(units) .eq. "Hours") .or. (trim(units) .eq. "h") ) then
         time = time / 3600.0_RP

      else if ( (trim(units) .eq. "days") .or. (trim(units) .eq. "Days") .or. (trim(units) .eq. "d") ) then
         time = time / (24.0_RP * 3600.0_RP)
         
      end if

   end function Event_GetLastElapsedTime
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Event_GetCPUTime(self,units) result (time)
      implicit none
      class(Event_t)               :: self
      character(len=*), intent(in) :: units
      real(kind=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP)     :: tActual

      if ( self % running ) then
         call CPU_TIME(tActual)
         time = tActual - self % tStartCPU + self % CPUTime

      else
         time = self % CPUTime

      end if

      if ( (trim(units) .eq. "minutes") .or. (trim(units) .eq. "Minutes") .or. (trim(units) .eq. "min") ) then
         time = time / 60.0_RP
   
      else if ( (trim(units) .eq. "hours") .or. (trim(units) .eq. "Hours") .or. (trim(units) .eq. "h") ) then
         time = time / 3600.0_RP

      else if ( (trim(units) .eq. "days") .or. (trim(units) .eq. "Days") .or. (trim(units) .eq. "d") ) then
         time = time / (24.0_RP * 3600.0_RP)
         
      end if

   end function Event_GetCPUTime

   subroutine Event_Reset(self)
      implicit none
      class(Event_t)          :: self
      
      self % tStart = 0.0_RP
      self % tStartCPU = 0.0_RP
      self % elapsedTime = 0.0_RP
      self % CPUTime    = 0.0_RP
      self % running = .false.

   end subroutine Event_Reset

   subroutine Event_Destruct(self)
      implicit none
      class(Event_t)       :: self

      self % tStart = 0.0_RP
      self % tStartCPU = 0.0_RP
      self % elapsedTime = 0.0_RP
      self % CPUTime    = 0.0_RP
      self % running = .false.
      self % next => NULL()

   end subroutine Event_Destruct
!
!/////////////////////////////////////////////////////////////////////
!
!           AUXILIARY PROCEDURES
!           -------------------
!////////////////////////////////////////////////////////////////////
!
   function GetCurrentSeconds() result(sec)
      implicit none
      real(kind=RP)         :: sec
!
!     ---------------
!     Local variables
!     ---------------
!
      integer, dimension(8) :: currentDate
      integer               :: no_of_days

      call date_and_time (values = currentDate)

      
      no_of_days = DaysBetweenTwoDates( reference_date(IDAY) , reference_date(IMONTH) , reference_date(IYEAR) , &
                                        currentDate(IDAY)    , currentDate(IMONTH)    , currentDate(IYEAR) )

      sec = no_of_days * 24.0_RP * 3600.0_RP + SecondsToMidnight(currentDate(IHOUR),&
                                                                 currentDate(IMINUTES),&
                                                                 currentDate(ISECONDS),&
                                                                 currentDate(IMILLISECONDS) )
   end function GetCurrentSeconds
   
   function SecondsToMidnight(hour,minute,seconds,millisecs) result (secs)
      implicit none
      integer, intent(in)     :: hour
      integer, intent(in)     :: minute
      integer, intent(in)     :: seconds
      integer, intent(in)     :: millisecs
      real(kind=RP)           :: secs

      secs = seconds + millisecs/1000.0_RP + minute * 60.0_RP + hour * 3600.0_RP

   end function SecondsToMidnight

   function DaysBetweenTwoDates(d1,m1,y1,d2,m2,y2) result (days)
      implicit none
      integer, intent(in)        :: d1
      integer, intent(in)        :: m1
      integer, intent(in)        :: y1
      integer, intent(in)        :: d2
      integer, intent(in)        :: m2
      integer, intent(in)        :: y2
      integer                    :: days

      days = DaysBetweenTwoYears(y1,y2) + DaysToNewYearsDay(d2,m2,y2) - DaysToNewYearsDay(d1,m1,y1)

   end function DaysBetweenTwoDates

   function DaysBetweenTwoYears(year1,year2) result(days)
      implicit none
      integer, intent(in)        :: year1 
      integer, intent(in)        :: year2
      integer                    :: days
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: no_of_leapYears
      integer        :: i 

      no_of_leapYears = 0

      do i = year1 , year2 - 1
         if ( checkIfLeapYear(i) ) then
            no_of_leapYears = no_of_leapYears + 1 
         end if
      end do

      days = 365 * (year2-year1) + no_of_leapYears

   end function DaysBetweenTwoYears

   function DaysToNewYearsDay(day,month,year) result(days)
      implicit none
      integer, intent(in)        :: day
      integer, intent(in)        :: month
      integer, intent(in)        :: year
      integer                    :: days
      integer, parameter         :: days_per_month(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: i

      days = day

      do i = 1 , month-1
         days = days + days_per_month(i)
      end do

      if ( checkIfLeapYear(year) .and. (month .gt. 2)) then
         days = days + 1 
      end if 

   end function DaysToNewYearsDay

   function checkIfLeapYear(year) result (leap)
      implicit none
      integer, intent(in)     :: year
      logical                 :: leap

      if ( mod(year,4) .ne. 0 ) then
         leap = .false.
         return
      end if

      if ( (mod(year,100) .eq. 0) .and. (mod(year,400) .ne. 0 )) then
         leap = .false.
         return

      else
         leap = .true.
         return
      end if

   end function checkIfLeapYear

end module StopwatchClass