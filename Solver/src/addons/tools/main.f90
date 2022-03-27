!
!   @File:    main.f90
!   @Author:  Oscar Marino (oscar.marino@upm.es)
!   @Created: May 3 2020
!   @Last revision date: 28 Feb 2022
!   @Last revision author: 
!   @Last revision commit: 
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
Program main
    use SMConstants
    use FTValueDictionaryClass
    use FileReaders               , only: ReadControlFile 
    use Headers
    use StopwatchClass
    use SharedBCModule
    use MPI_Process_Info
    use LocalRefinementTool
    use FWHTools
#ifdef _HAS_MPI_
      use mpi
#endif
    Implicit None

    TYPE( FTValueDictionary)                :: controlVariables
    character(len=LINE_LENGTH)              :: toolType

!   -----------------------------------------
!   Start measuring the total simulation time
!   -----------------------------------------
!
    call Stopwatch % CreateNewEvent("TotalTime")
    call Stopwatch % Start("TotalTime")

!   ---------------
!   Initializations
!   ---------------
!
    call MPI_Process % Init
    if ( MPI_Process % doMPIAction ) then
        write(STD_OUT,'(A)') "HORSES additional tool cannot be run with MPI"
        call exit(99)
    end if
!
    call controlVariables % initWithSize(16)
    call ConstructSharedBCModule
!
    call ReadControlFile( controlVariables )
!
!   --------------
!   Main execution
!   --------------
!
    toolType = controlVariables % stringValueForKey("tool type", requestedLength = LINE_LENGTH) 
    select case (trim(toolType))
    case ("fwh post")
        call Main_Header("HORSES additional proccesing tools: FWH Post-Proccesing",__DATE__,__TIME__)
        call FWHTool(controlVariables, isSurface=.false.)
    case ("fwh surface")
        call Main_Header("HORSES additional proccesing tools: FWH surface Pre-Proccesing",__DATE__,__TIME__)
        call FWHTool(controlVariables, isSurface=.true.)
    case ("local p refinement")
        call Main_Header("HORSES additional proccesing tools: Local Refinement Pre-Proccesing",__DATE__,__TIME__)
        call LocalRef(controlVariables)
    case default
        call Main_Header("HORSES additional proccesing tools",__DATE__,__TIME__)
        write(STD_OUT,'(A)') "The requested tool type is not implemented"
        write(STD_OUT,'(A)') "Implemented types are:"
        write(STD_OUT,'(A)') "  * fwh post"
        write(STD_OUT,'(A)') "  * fwh surface"
        write(STD_OUT,'(A)') "  * local p refinement"
    end select
!
!   ---------
!   Finish up
!   ---------
!
    call Stopwatch % destruct
    call destructSharedBCModule

    call MPI_Process % Close
   
End Program main
