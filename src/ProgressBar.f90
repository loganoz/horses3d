!+++++++++++++++++++++++++++++++++++++++++
!> \brief Module containing a set of classes that can be used as progress-bars
!! and status counters. These counters only write to the standard output.
!!
!! This class is part of the online Fortran 2003 Tutorial by
!! Danny E.P. Vanpoucke (http://dannyvanpoucke.be)
!! and comes without warranty.
!!
!! This module makes use of:
!! - NOTHING
!<-----------------------------------------------------------------------
module progressbarsmodule
    use SMConstants
    implicit none
    private



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> \brief Displays a text progress bar in the command-line.
    !!
    !! By making use of the carriage return character it is possible to generate a
    !! progress bar, by continuously overwriting the same line.\n
    !! The progress bar appear as:\n
    !! \par
    !! \b example \n
    !!            Your Progress Message, 43% [========
    !!
    !<-----------------------------------------------------------------------
    type, public :: TProgressBar
        private
        logical :: init
        logical :: running
        logical :: done
        character(len=255) :: message
        character(len=30) :: progressString
        character(len=20) :: bar
        real :: progress
    contains
        private
        procedure,pass(this),public :: initialize
        procedure,pass(this),public :: reset
        procedure,pass(this),public :: run
        procedure,pass(this),private:: printbar
        procedure,pass(this),private:: updateBar
    end type TProgressBar


contains

    !+++++++++++++++++++++++++++++++++++
    !>\brief Initialize the TProgressBar object.
    !!
    !! @param[in] msg String containing message to be displayed with the progressbar.[\b OPTIONAL ,\b DEFAULT = none]
    !<----------------------------------
    subroutine initialize(this,msg)
        class(TProgressBar) :: this
        character(len=*),intent(in),optional :: msg

        call this%reset()
        if(present(msg)) this%message=msg
        this%init=.true.

    end subroutine initialize
    !+++++++++++++++++++++++++++++++++++
    !>\brief Reset the TProgressBar object.
    !<----------------------------------
    subroutine reset(this)
        class(TProgressBar) :: this

        this%init=.false.
        this%done=.false.
        this%running=.false.
        this%message=""
        this%progressString=""
        this%progress=0.0

    end subroutine reset
    !+++++++++++++++++++++++++++++++++++
    !>\brief Run the TProgressBar object.
    !!
    !! @param[in] pct Real value providing a the % of progress.
    !! @param[in] Ix Integer number providing the number of digits to be used in the %. [\b OPTIONAL ,\b DEFAULT = 2]
    !! @param[in] msg String containing message to be displayed with the progressbar.[\b OPTIONAL ,\b DEFAULT = none/ provided by the initialisation]
    !<----------------------------------
    subroutine run(this,pct,Ix,msg)
        class(TProgressBar) :: this
        real::pct
        integer, intent(in), optional :: Ix
        character(len=*),intent(in),optional :: msg
        if (.not. this%init) call this%initialize(msg)
        if (.not. this%done) then
            this%running=.true.
            this%progress=pct
            call this%updateBar(Ix)
            call this%printbar()
            if (abs(pct-100.0)<1.0E-6) then
                this%done=.true.
                write(*,'(A6)') "] done"
            end if
        end if

    end subroutine run
    !+++++++++++++++++++++++++++++++++++
    !>\brief Update the bar of the TProgressBar object.
    !!
    !! @param[in] Ix Integer number providing the number of digits to be used in the %. [\b OPTIONAL ,\b DEFAULT = 2]
    !<----------------------------------
    subroutine updateBar(this,Ix)
        class(TProgressBar) :: this
        integer, intent(in), optional :: Ix

        integer :: Ixx,np
        character(len=50)::fm, fb

        Ixx=2
        if (present(Ix)) then
            if (Ix>=0) Ixx=Ix
        end if
        fb="========================"

        write(fm,'(A2,I0,A1,I0,A1)')"(F",Ixx+5,".",Ixx,")"
        write(this%progressString,trim(fm)) this%progress
        ! every new bar section per 5% (20 pieces)
        np=NINT(this%progress/5)

        if(np/=len_trim(this%bar)) then! things changed=redo bar
            this%bar=fb(1:np)
        end if

    end subroutine updateBar
    !+++++++++++++++++++++++++++++++++++
    !>\brief Print the TProgressBar object.
    !!
    !<----------------------------------
    subroutine printbar(this)
        ! use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
        class(TProgressBar) :: this

        character(len=50)::fm

        !fmt='(A'//trim(adjustl(ls))//',2X,1I3,1A1,2X,1A1,256A1)'
        fm='(A1,A,X,2A,2X,A1,A)'
        !start with carriage return to stay on the same line.
        write(STD_OUT,trim(fm), advance='NO') achar(13),&
            &trim(this%message),trim(adjustl(this%progressString)),'%','[',trim(adjustl(this%bar))
        ! flush(OUTPUT_UNIT)
        flush(STD_OUT)
    end subroutine printbar


end module progressbarsmodule