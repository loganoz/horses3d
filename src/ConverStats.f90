#include "Includes.h"
Module ConverStats  !
    use SMConstants
    use SolutionFile
    use FTValueDictionaryClass
    Implicit None

    contains

    Subroutine ConvertStatsForRestart(controlVariables)
        use FileReadingUtilities, only: getFileName
        implicit none
        TYPE( FTValueDictionary), intent(in)                    :: controlVariables

        !local variables
        character(LEN=LINE_LENGTH)                              :: fileName, restartFileName
        integer                                                 :: fid, eID, newFid
        real(kind=RP), allocatable, dimension(:,:,:,:)          :: Q,stats,Q_x
        integer                                                 :: no_of_elements
        integer                                                 :: NVARS, NSTAT
        integer                                                 :: padding, dimensionsSize, nodeType
        integer, dimension(NDIM)                                :: Nsol
        integer                                                 :: iter
        integer, dimension(:),  allocatable                     :: arrayDimensions
        real(kind=RP)                                           :: time
        logical                                                 :: hasGradients
        real(kind=RP)                                           :: refs(NO_OF_SAVED_REFS)

        fileName = controlVariables % stringValueForKey("stats file", LINE_LENGTH)
        hasGradients = controlVariables % logicalValueForKey("has gradients")
        !remove hsol and stats extensions
        restartFileName = trim(getFileName(fileName))
        restartFileName = trim(getFileName(restartFileName))
        write(restartFileName,'(A,A)') trim(restartFileName), '_restart.hsol'

         NVARS = 5
         NSTAT = 9

        padding = NVARS
        select case ( getSolutionFileType(trim(fileName)) )
         case (STATS_FILE)

         case default
           print*, "File expected to be a stats file"
           errorMessage(STD_OUT)
           error stop
        end select

         no_of_elements = getSolutionFileNoOfElements(fileName)
            dimensionsSize = 4
         allocate(arrayDimensions(dimensionsSize))
         call getSolutionFileTimeAndIteration(trim(fileName),iter,time)

         refs = getSolutionFileReferenceValues(trim(fileName))
         nodeType = getSolutionFileNodeType(fileName)

         fid = putSolutionFileInReadDataMode(fileName)

         ! NVARS = arrayDimensions(1)

!        Create new file
!        ---------------
        call CreateNewSolutionFile(trim(restartFileName),SOLUTION_FILE, nodeType, &
                                       no_of_elements, iter, time, refs)

        newFid = putSolutionFileInWriteDataMode(trim(restartFileName))

        do eID = 1, no_of_elements
            call getSolutionFileArrayDimensions(fid,arrayDimensions)

!   
            Nsol(1:3) = arrayDimensions(2:4) - 1
!
!           Allocate memory for the statistics
!           ----------------------------------
            allocate( stats(1:NSTAT,0:Nsol(1), 0:Nsol(2), 0:Nsol(3)) )
            read(fid) stats
            allocate( Q(1:NVARS,0:Nsol(1),0:Nsol(2),0:Nsol(3)) )
!
!           Read data
!           ---------
            read(fid) Q
!
            if ( hasGradients ) then
               allocate( Q_x(1:NVARS,0:Nsol(1),0:Nsol(2),0:Nsol(3)) )
!              Read data
!              ---------
               read(fid) Q_x
               read(fid) Q_x
               read(fid) Q_x
            end if

            ! pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + padding*e % offsetIO * SIZEOF_RP
            call writeArray(newFid, Q)
            safedeallocate(stats)
            safedeallocate(Q)
            safedeallocate(Q_x)

        end do
         close(fid)
         close(newFid)
         call SealSolutionFile(trim(restartFileName))

    End Subroutine ConvertStatsForRestart
End Module ConverStats