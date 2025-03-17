!
!//////////////////////////////////////////////////////
!
!  Program to merge statistics files
!  Usage:
!
!  horses.mergeStats *.hsol --initial-iteration=INTEGER --file-name=CHARACTER
!
!  NOTES:
!  *  Only usable with statistics files that are obtained with the "reset interval" keyword
!     and/or with individual consecutive simulations 
!  *  Only constant time-stepping is supported
!  *  Dynamic p-adaptation is currently not supported
!
program MergeStatistics
   use SMConstants
   use getTaskModule       , only: getTask
   use Headers  !           , only: Main_Header, Section_Header
   use StatsAveragingModule, only: Initialize_StatsAveragingModule, Finalize_StatsAveragingModule, PerformAveraging
   use MPI_Process_Info    , only: MPI_Process
   implicit none
!
!  -----------------
!  Program variables
!  -----------------
!
   integer                                 :: num_of_statFiles, initial_iteration
   real(kind=RP), allocatable              :: weights(:)
   character(len=LINE_LENGTH), allocatable :: fileNames(:)
   character(len=LINE_LENGTH)              :: fileName
!
!  ------------------------
!  Initializations
!  ------------------------
!
   call MPI_Process % Init
   
   call Main_Header("HORSES Merge Statistics Utility",__DATE__,__TIME__)
   
   write(STD_OUT,'(/,/)')
   call Section_Header("Job description")
   write(STD_OUT,'(/,/)')
   
   call getTask(num_of_statFiles, fileNames, weights, initial_iteration, fileName)
   
   call Initialize_StatsAveragingModule(num_of_statFiles, fileNames)
   
!
!  ------------------------
!  Average statistics files
!  ------------------------
!
   call PerformAveraging(num_of_statFiles, fileNames, weights, fileName)
!
!  ---------
!  Finish up
!  ---------
!
   call Finalize_StatsAveragingModule
   
   call MPI_Process % Close
   
end program MergeStatistics