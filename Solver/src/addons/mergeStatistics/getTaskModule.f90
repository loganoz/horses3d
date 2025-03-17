!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!     This module explores the command line flags introduced by the user. Software usage is:
!
!           $ horses.mergeStats *.stats.hsol --intitial-iteration=INTEGER
!
!  Initial iteration flag is optional.
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module getTaskModule
   use SMConstants
   use SolutionFile, only: getSolutionFileType, getSolutionFileTimeAndIteration, STATS_FILE
   use StatsAveragingModule, only: hasGradients
   implicit none
   
   private
   public   getTask
   
   integer, parameter   :: EXPORT_GAUSS = 0
   integer, parameter   :: EXPORT_HOMOGENEOUS = 1

   character(len=*), parameter   :: INITIAL_ITER_FLAG="--initial-iteration="
   character(len=*), parameter   :: FILE_NAME_FLAG="--file-name="
   character(len=*), parameter   :: GRADIENTS_FLAG="--has-gradients"

   contains

      subroutine getTask(num_of_statFiles, fileNames, weights, initial_iteration, fileName)
         implicit none
         !-arguments----------------------------------------------------
         integer,                                 intent(out) :: num_of_statFiles
         character(len=LINE_LENGTH), allocatable, intent(out) :: fileNames(:)
         real(kind=RP)             , allocatable, intent(out) :: weights(:)
         integer                                , intent(out) :: initial_iteration
         character(len=LINE_LENGTH)             , intent(out) :: fileName
         !-local-variables----------------------------------------------
         integer                    :: no_of_arguments
         integer                    :: i, sol, pos, pos2
         integer                    :: io, fid
         integer, allocatable       :: last_iter(:)
         integer, allocatable       :: iter_number(:)
         character(len=LINE_LENGTH) :: auxiliarName
         real(kind=RP)  :: time
         !--------------------------------------------------------------
!
!        Get number of command arguments         
!        -------------------------------
         no_of_arguments = command_argument_count()
!
!        Exit if no input arguments are specified
!        ----------------------------------------
         if ( no_of_arguments .eq. 0 ) then
            write(STD_OUT,'(A)') "No statistics file(s) specified"
            stop
         end if
!
!        Check if the solution file is present
!        -------------------------------------
!
!        Loop to get number of files
!        ---------------------------
         num_of_statFiles = 0
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            open(newunit=fid, file=trim(auxiliarName), action="read", form="unformatted", access="stream", iostat=io)
            close(fid)
            if ( io .ne. 0 ) cycle

            if ( getSolutionFileType(auxiliarName) .eq. STATS_FILE ) then
               num_of_statFiles = num_of_statFiles + 1 
            end if
         end do
!
!        Loop to get solution file names
!        -------------------------------
         if ( num_of_statFiles .ne. 0 ) then
            allocate( fileNames  (num_of_statFiles) )
            allocate( last_iter  (num_of_statFiles) )
            allocate( iter_number(num_of_statFiles) )
            allocate( weights    (num_of_statFiles) )
            
            sol = 0
            do i = 1, no_of_arguments
               call get_command_argument(i, auxiliarName)
               open(newunit=fid, file=trim(auxiliarName), action="read", form="unformatted", access="stream", iostat=io)
               close(fid)
               if ( io .ne. 0 ) cycle

               if ( getSolutionFileType(auxiliarName) .eq. STATS_FILE ) then   
                  sol = sol + 1 
                  fileNames(sol) = trim(auxiliarName)
                  call getSolutionFileTimeAndIteration(trim(auxiliarName), last_iter(sol), time)
               end if
            end do
         else
            write(STD_OUT,'(A)') "No statistics file(s) specified"
            stop
         end if
!
!        Get the initial iteration
!        -------------------------
         initial_iteration = 1
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName), INITIAL_ITER_FLAG)
            
            if ( pos .ne. 0 ) then
               read(auxiliarName(pos+len_trim(INITIAL_ITER_FLAG):len_trim(auxiliarName)),*) initial_iteration
               exit
            end if
         end do
         
!
!        Get the file name
!        -----------------
         fileName = ''
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName), FILE_NAME_FLAG)
            
            if ( pos .ne. 0 ) then
               fileName = auxiliarName(pos+len_trim(FILE_NAME_FLAG):len_trim(auxiliarName))
               exit
            end if
         end do
         if (fileName == '') fileName = 'Merged.stats.hsol'
         
         ! order files?

!        Get if has gradients
!        -----------------
         hasGradients = .false.
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName), GRADIENTS_FLAG)
            
            if ( pos .ne. 0 ) then
               hasGradients = .true.
               exit
            end if
         end do
         
!
!        Get the weights
!        ---------------
         do i = num_of_statFiles, 2, -1
            if (last_iter(i) .le. last_iter(i-1) ) print*, 'ERROR: files are not ordered'
            iter_number(i) = last_iter(i) - last_iter(i-1)
         end do
         iter_number(1) = last_iter(1) - (initial_iteration -   1)
         
         weights = iter_number / real(sum(iter_number),RP)
         
!        Describe
!        --------
         
         write(*,'(A60,A14,A14,A14,A14)') 'Statistics file', 'Initial iter', 'Final iter', 'Iter number', 'Weight'
         write(*,'(A60,A14,A14,A14,A14)') '---------------', '------------', '----------', '-----------', '------'
         
         ! First file
         write(*,'(A60,I14,I14,I14,F14.6)') trim( fileNames(1) ), initial_iteration, last_iter(1), iter_number(1), weights(1)
         
         do i = 2, num_of_statFiles
            write(*,'(A60,I14,I14,I14,F14.6)') trim( fileNames(i) ), last_iter(i-1)+1, last_iter(i), iter_number(i), weights(i)
         end do
         
      end subroutine getTask
end module getTaskModule
