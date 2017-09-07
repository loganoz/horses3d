!
! /////////////////////////////////////////////////////////////////////
!
!
!     TimerClass.f90
!
!!     Created on: June 6, 2002
!!
!!     @Author     David Kopriva
!!
!!     Modification History:
!!
!
!      MODULE TimerClass
!
!!        A timer has the ability to time the execution of a code. 
!!        It works like a stopwatch. To start a timer, call "StartTimer"
!!        and stop it with "StopTimer". To check the elapsed time, use
!!        "ElapsedTimeOnTimer" to return the elapsed time. This routine
!!        has an optional argument "resultUnits" that tells what units
!!        ( seconds, minutes or hours ) the returned values will have.
!!
!        TYPE Timer
!        USES PrecisionDefinitions
!
!        PUBLIC DATA:
!           None
!        PUBLIC METHODS:
!
!        Constructor/Destructor
!
!           SUBROUTINE ConstructTimer( this )
!           SUBROUTINE DestructTimer( this )
!
!        Running a timer
!
!           SUBROUTINE StartTimer( this )
!           SUBROUTINE RestartTimer(this)
!           SUBROUTINE StopTimer( this )
!
!        Enquiry functions
!
!        REAL( KIND=RP ) FUNCTION ElapsedTimeOnTimer( this, resultUnits )
!                 CHARACTER( LEN=7 ), OPTIONAL, INTENT( IN ) :: resultUnits
!        CHARACTER(LEN=10) FUNCTION StartTimeOnTimer( this )
!        CHARACTER(LEN=10) FUNCTION EndTimeOnTimer  ( this )
!        SUBROUTINE        PrintStartTime  ( this, unit )
!        SUBROUTINE        PrintEndTime    ( this, unit )
!
! /////////////////////////////////////////////////////////////////////
!
!  ******
   MODULE TimerClass
!  ******
!
     USE SMConstants         
     IMPLICIT NONE
!
!    ---------------
!    Type definition
!    ---------------
!         
!    ------------------------------------------------
!!   Provides an interface to the F90 timing routine
!    ------------------------------------------------
!
     TYPE Timer
        CHARACTER( LEN=8 )    :: date
        CHARACTER( LEN=10 )   :: start_time_str, end_time_str
        CHARACTER( LEN=5 )    :: zone
        INTEGER, DIMENSION(8) :: start_time_values, end_time_values
        REAL(KIND=RP)         :: elapsedTime
!
     END TYPE Timer
!
!    ========
     CONTAINS
!    ========
!
!     /////////////////////////////////////////////////////////////////
!
!     ---------------------------------
!!    Constuctor:Initial setup of timer
!     ---------------------------------
!
      SUBROUTINE ConstructTimer( this )
!
      TYPE( Timer ) :: this
!
      CALL DATE_AND_TIME( this%date, this%start_time_str, this%zone, &
                         this%start_time_values )
      this%end_time_str    = this%start_time_str
      this%end_time_values = this%start_time_values
      this%elapsedTime     = 0.0_RP
!
      END SUBROUTINE ConstructTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     --------------------------------------------
!!    Destructor:This destructor has nothing to do
!     --------------------------------------------
!
      SUBROUTINE DestructTimer( this )
!
      TYPE( Timer ) :: this
      
      END SUBROUTINE DestructTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------------------
!!    Start the stopwatch: Call this when timing is to start.
!     ---------------------------------------------------------------
!
      SUBROUTINE StartTimer( this )
!
      TYPE( Timer ) :: this

      CALL DATE_AND_TIME( this%date, this%start_time_str, this%zone, &
                         this%start_time_values )
      this%elapsedTime = 0.0_RP
      
      END SUBROUTINE StartTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------------------
!!    Start the stopwatch: Call this when timing is to start, but 
!!    accumulate the time.
!     ---------------------------------------------------------------
!
      SUBROUTINE RestartTimer( this )
!
      TYPE( Timer )     :: this
      REAL(KIND=RP) :: elapsed

      elapsed            = ElapsedTimeOnTimer( this )
      this%elapsedTime   = elapsed
      CALL DATE_AND_TIME( this%date, this%start_time_str, this%zone, &
                         this%start_time_values )
      
      END SUBROUTINE RestartTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!!    Stop the stopwatch: Call this to stop this timer.
!     ----------------------------------------------------------
!
      SUBROUTINE StopTimer( this )
!
      TYPE( Timer ) :: this

      CALL DATE_AND_TIME( this%date, this%end_time_str, this%zone, &
                         this%end_time_values )
      
      END SUBROUTINE StopTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!     Reads and returns the elapsed time on the
!!                   timer. This routine will give the wrong answer for 
!!                   leap year near march 1 or any computation run over
!!                   new year's eve. You shouldn't be working
!!                   at midnight on New Year's eve anyway.The input
!!                   variable resultUnits is a Character*7 optional
!!                   variable that is either "seconds", "minutes" or
!!                   "hours". The default is minutes.
!     ----------------------------------------------------------------
!
      REAL( KIND=RP ) FUNCTION ElapsedTimeOnTimer( this, resultUnits )
!
      TYPE( Timer ), INTENT( INOUT )             :: this
      CHARACTER( LEN=* ), OPTIONAL, INTENT( IN ) :: resultUnits
      
      INTEGER, DIMENSION(12) :: no_of_days = (/31, 28, 31, 30, 31, 30, 31, &
                                                  31, 30, 31, 30, 31/)
      INTEGER, DIMENSION(8)  :: elapsed
      REAL( KIND=RP )    :: elapsed_time
      INTEGER                :: month
!
      elapsed = this%end_time_values - this%start_time_values
!
      elapsed_time = 0.0_RP
      DO month = this%start_time_values(2), this%end_time_values(2)-1
         elapsed_time = no_of_days( month ) + elapsed_time
      END DO
      elapsed_time = elapsed_time*24*60

      elapsed_time = elapsed_time &
                   + ( elapsed(3)*24._RP + elapsed(5) )*60._RP &
                   + elapsed(6) + ( elapsed(7) &
                   + elapsed(8)/1000._RP )/60._RP
      ElapsedTimeOnTimer = elapsed_time + this%elapsedTime
      
      IF (PRESENT( resultUnits ) )     THEN
         SELECT CASE( resultUnits )
            CASE( "seconds" )
               ElapsedTimeOnTimer = ElapsedTimeOnTimer*60._RP
            CASE( "hours" )
               ElapsedTimeOnTimer = ElapsedTimeOnTimer/60._RP
         END SELECT
      END IF
      
      END FUNCTION ElapsedTimeOnTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------
!!    StartTimeOnTimer: Return the time the timer was started
!     -------------------------------------------------------
!
      CHARACTER(LEN=10) FUNCTION StartTimeOnTimer( this )
!
      TYPE( Timer ) :: this
!
      StartTimeOnTimer = this%start_time_str
!
      END FUNCTION StartTimeOnTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------
!!    EndTimeOnTimer: Return the time the timer was stopped
!     -------------------------------------------------------
!
      CHARACTER(LEN=10) FUNCTION EndTimeOnTimer( this )
!
      TYPE( Timer ) :: this
!
      EndTimeOnTimer = this%end_time_str
!
      END FUNCTION EndTimeOnTimer
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------
!!    PrintStartTime: Print the time the timer was started
!     ----------------------------------------------------
!
      SUBROUTINE PrintStartTime( this, unit )
!
      TYPE( Timer ) :: this
      INTEGER       :: unit
!
      WRITE( unit, * ) this%start_time_str
!
      END SUBROUTINE PrintStartTime
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------
!!    PrintStartTime: Print the time the timer was started
!     ----------------------------------------------------
!
      SUBROUTINE PrintEndTime( this, unit )
!
      TYPE( Timer ) :: this
      INTEGER       :: unit
!
      WRITE( unit, * ) this%end_time_str
!
      END SUBROUTINE PrintEndTime
!
!  **********     
   END MODULE TimerClass
!  **********     
