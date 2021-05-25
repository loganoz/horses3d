#include "Includes.h"
Module FWHPostProc  !

    use SMConstants
    use FWHGeneralClass
    use DGSEMClass
    use SolutionFile
    Implicit None

    contains

    Subroutine LoadAllFiles(controlVariables, sem)

        implicit none
        TYPE(DGSem)                                             :: sem
        TYPE( FTValueDictionary), intent(in)                    :: controlVariables

        !local variables
        integer                                                 :: i, numberOfFiles, reason, fID
        character(len=LINE_LENGTH), dimension(:), allocatable   :: fileNames
        character(LEN=LINE_LENGTH)                              :: pattern, fullExpression
        real(kind=RP)                                           :: r

        call sem % fwh % autosaveConfig(controlVariables, 0.0_RP)

        ! get files in temporal txt
        pattern = controlVariables % stringValueForKey("accoustic files pattern", LINE_LENGTH)
        write(fullExpression,'(A,A,A)') "ls ", trim(pattern), " > horses_temporal_file.txt"
        call system(trim(fullExpression))

        open ( newunit = fID , file = "horses_temporal_file.txt", status = "old" , action = "read" ) 
        ! get the number of files that match the pattern
        i = 0
        do
            read(fID,fmt='(a)',iostat=reason) r
            if (reason/=0) exit
            i = i+1
        end do
        numberOfFiles = i

        ! get files names in array
        allocate(fileNames(numberOfFiles))
        rewind(fID)
        do i = 1, numberOfFiles
            read(fID, '(A)') fileNames(i)
        end do 
        close(fID)

        do i = 1, numberOfFiles
            call LoadSingleFile(fileNames(i), sem % mesh, sem % fwh, i)
        end do 

        call sem % fwh % writeToFile(force=.TRUE.)

        sem % numberOftimeSteps = numberOfFiles

    End Subroutine LoadAllFiles

    Subroutine LoadSingleFile(fileName, mesh, fwh, actualIter)

        Implicit None
        character(len=*), intent(in)                            :: fileName
        class (HexMesh), intent(inout)                          :: mesh
        class(FWHClass), intent(inout)                          :: fwh
        integer, intent(in)                                     :: actualIter

        !local variables
        integer                                                 :: iter
        real(kind=RP)                                           :: time
        integer                                                 :: fID, eID, fileType, no_of_elements, flag, nodetype

!       Get the file type
!       -----------------
        fileType = getSolutionFileType(trim(fileName))

        select case (fileType)
            case(MESH_FILE)
               print*, "The selected restart file is a mesh file"
               errorMessage(STD_OUT)
               stop

            case(SOLUTION_FILE)
               print*, "The selected file is a variable solution file"
               errorMessage(STD_OUT)
               stop

            case(SOLUTION_AND_GRADIENTS_FILE)
               print*, "The selected file is a variable solution file"
               errorMessage(STD_OUT)
               stop

            case(STATS_FILE)
               print*, "The selected file is a statistics file"
               errorMessage(STD_OUT)
               stop

            case(ZONE_SOLUTION_FILE)

            case default
               print*, "Unknown restart file format"
               errorMessage(STD_OUT)
               stop
        end select
!
!       Get the node type
!       -----------------
        nodeType = getSolutionFileNodeType(trim(fileName))

        if ( nodeType .ne. mesh % nodeType ) then
           print*, "Solution file uses a different discretization nodes than the mesh."
           errorMessage(STD_OUT)
        end if
!
!       Read the number of elements
!       ---------------------------
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))

        if ( no_of_elements .ne. fwh % sourceZone % no_of_faces ) then
           write(STD_OUT,'(A,A)') "The number of faces stored in the file ", &
                                  "do not match that of the surface defined"
           errorMessage(STD_OUT)
           stop
        end if
!
!       Read the initial iteration and time
!       -----------------------------------
        call getSolutionFileTimeAndITeration(trim(fileName), iter, time)
!
!       Read the terminator indicator
!       -----------------------------
        flag = getSolutionFileDataInitFlag(trim(fileName))

        if ( flag .ne. BEGINNING_DATA ) then
           print*, "Beginning data flag was not found in the file."
           errorMessage(STD_OUT)
           stop
        end if

!
!       Load Solution to FWH Surface
!       -----------------------------
        call fwh % loadSourceSol(trim(fileName), mesh)
!
!       Update and Write if necessary
!       -----------------------------
        call fwh % updateValues(mesh, time, iter, isFromFile=.TRUE.)
        call fwh % writeToFile()
!
!       Display
!       -----------------------------
        print *, "iter: ", actualIter, "file used: ", trim(fileName), "with time: ", time

    End Subroutine LoadSingleFile

End Module FWHPostProc
