
!//////////////////////////////////////////////////////
!
!This class represents the general behaviour of the Ffowcs Williams and Hawckings aero acoustic analogy
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
Module FWHGeneralClass  !

    use SMConstants
    use FWHDefinitions, only: OB_BUFFER_SIZE_DEFAULT, OB_BUFFER_SIZE_DEFAULT, STR_LEN_OBSERVER
    use FWHObseverClass
    use HexMeshClass
    use FileReadingUtilities, only: getFileName
    use SurfaceMesh, only: surfacesMesh, FWH_POSITION, SURFACE_TYPE_FWH
    Implicit None

!
!   *****************************
!   Main FWH class definition
!   *****************************
    type FWHClass

        character(len=LINE_LENGTH)                                        :: solution_file
        integer                                                           :: numberOfObservers = 0
        integer                                                           :: bufferLine
        integer, dimension(:), allocatable                                :: iter
        real(kind=RP), dimension(:), allocatable                          :: t
        class(ObserverClass), dimension(:), allocatable                   :: observers
        integer                                                           :: totalNumberOfFaces
        logical                                                           :: isSolid
        logical                                                           :: isActive = .false.
        logical                                                           :: firstWrite
        logical                                                           :: interpolate

        contains

            procedure :: construct      => FWHConstruct
            procedure :: destruct       => FWHDestruct
            procedure :: updateValues   => FWHUpate
            procedure :: writeToFile    => FWHWriteToFile

    end type FWHClass

    contains

    Subroutine FWHConstruct(self, mesh, controlVariables)

        use FTValueDictionaryClass
        use mainKeywordsModule
        use FileReadingUtilities, only: getCharArrayFromString
        use FWHDefinitions,       only: getMeanStreamValues
        use Headers
        use Utilities,            only: toLower
#ifdef _HAS_MPI_
        use mpi
#endif
        implicit none

        class(FWHClass)                                     :: self
        class(HexMesh), intent(in)                          :: mesh
        class(FTValueDictionary), intent(in)                :: controlVariables

!       ---------------
!       Local variables
!       ---------------
!
        integer                                             :: fID , io
        integer                                             :: i
        character(len=STR_LEN_OBSERVER)                     :: line
        character(len=STR_LEN_OBSERVER)                     :: solution_file
        integer                                             :: no_of_zones, no_of_face_i, ierr, no_of_faces
        logical, save                                       :: FirstCall = .TRUE.

!        look if the acoustic analogy calculations are set to be computed
!        --------------------------------
        !TODO read acoustic analogy type and return if is not defined, check for FWH if is defined and not FWH stop and send error
        if (.not. controlVariables % containsKey("acoustic analogy")) then
            self % isActive = .FALSE.
            ! print *, "FWH not activated"
            return
        end if

!       check the that sourceZone is FWH
!       ----------------------------------
        if (surfacesMesh % surfaceTypes(FWH_POSITION) .ne. SURFACE_TYPE_FWH) then
            self % isActive = .FALSE.
            print *, "FWH surface not found, the FWH routines will not deactivated"
            return
        end if 

!       Setup the buffer
!       ----------------
        if (controlVariables % containsKey("observers flush interval") ) then
           OB_BUFFER_SIZE = controlVariables % integerValueForKey("observers flush interval")
        end if

        self % isActive = .TRUE.
        allocate( self % t(OB_BUFFER_SIZE), self % iter(OB_BUFFER_SIZE) )

!       Get the general configuration of control file
!       First get the surface as a zone
!       -------------------------------
        self % isSolid   = .not. controlVariables % logicalValueForKey("acoustic analogy permeable")

!       Get the solution file name
!       --------------------------
        solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_OBSERVER )
!
!       Remove the *.hsol termination
!       -----------------------------
        solution_file = trim(getFileName(solution_file))
        self % solution_file = trim(solution_file)

!       Search in case file for observers
!       ---------------------------------------------------------------------
        if (mesh % child) then ! Return doing nothing if this is a child mesh
           self % numberOfObservers = 0
        else
           self % numberOfObservers = getNoOfObservers()
        end if

!       Set interpolate attribute as TRUE by default
        ! todo: read from constrol variables
        self % interpolate = .TRUE.
        ! self % interpolate = .FALSE.

!       Initialize observers
!       --------------------
        call getMeanStreamValues()
        no_of_faces = surfacesMesh % totalFaces(SURFACE_TYPE_FWH)
        allocate( self % observers(self % numberOfObservers) )
        do i = 1, self % numberOfObservers
            call self % observers(i) % construct(surfacesMesh % zones(SURFACE_TYPE_FWH) , mesh, i, self % solution_file, FirstCall, &
                                                 self % interpolate, no_of_faces, surfacesMesh % elementSide(:,1))
        end do 

        self % bufferLine = 0
        self % firstWrite = .FALSE.
        
        FirstCall = .FALSE.

!        Describe the zones
!        ------------------
         if ( .not. MPI_Process % isRoot ) return
         call Subsection_Header("Fictitious FWH zone")
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of faces: ", no_of_faces
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of observers: ", self % numberOfObservers
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of integrals: ", self % numberOfObservers * no_of_faces
         write(STD_OUT,'(30X,A,A28,L1)') "->", "Save zone solution: ", controlVariables % containsKey("acoustic save timestep")

    End Subroutine FWHConstruct

    Subroutine FWHUpate(self, mesh, t, iter, isFromFile)

        implicit none

        class(FWHClass)                                     :: self
        class(HexMesh)                                      :: mesh
        real(kind=RP), intent(in)                           :: t
        integer, intent(in)                                 :: iter
        logical, intent(in), optional                       :: isFromFile

!       ---------------
!       Local variables
!       ---------------
!
        integer                                             :: i 
        logical                                             :: prolong

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!       Check if prolong is necessary
!       ------------------------
        if (present(isFromFile)) then
            prolong = .not. isFromFile
        else
            prolong = .TRUE.
        end if 
!
!       Move to next buffer line
!       ------------------------
        self % bufferLine = self % bufferLine + 1
!
!       Save time and iteration
!       -----------------------
        self % t       ( self % bufferLine )  = t
        self % iter    ( self % bufferLine )  = iter
!
!       Save Solution to elements faces of fwh surface
!       -----------------------
        if (prolong) call SourceProlongSolution(surfacesMesh % zones(SURFACE_TYPE_FWH), mesh, surfacesMesh % elementSide(:,1))

!       see if its regular or interpolated
!       -----------------------
        if (.not. self % firstWrite) then
            do i = 1, self % numberOfObservers
                call self % observers(i) % update(mesh, self % isSolid, self % bufferLine, self % interpolate)
            end do 
        else
            do i = 1, self % numberOfObservers
                call self % observers(i) % updateOneStep(mesh, self % bufferLine, self % isSolid, t)
            end do 
        end if

    End Subroutine FWHUpate

    Subroutine FWHWriteToFile(self, force)
!
!        ******************************************************************
!              This routine has a double behaviour:
!           force = .true.  -> Writes to file and resets buffers
!           force = .false. -> Just writes to file if the buffer is full
!        ******************************************************************
!
        use MPI_Process_Info
        implicit none

        class(FWHClass)                                     :: self
        logical, optional                                   :: force

!       ---------------
!       Local variables
!       ---------------
        integer                                             :: i 
        logical                                             :: forceVal

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

        if ( present ( force ) ) then
           forceVal = force
        else
           forceVal = .false.
        end if

        if ( forceVal ) then 
!
!           In this case the observers are exported to their files and the buffer is reset
!           ------------------------------------------------------------------------------
            if (.not. self % firstWrite .and. self % interpolate) then
                do i =1, self % numberOfObservers
                    call self % observers(i) % interpolateSol(self % t, self % bufferLine)
                end do
                self % firstWrite = .TRUE.
            end if
            do i =1, self%numberOfObservers
                call self % observers(i) % writeToFile(self % iter, self % t, self % bufferLine)
            end do
!               Reset buffer
!               ------------
                self % bufferLine = 0
        else
!               The observers are exported just if the buffer is full
!               ----------------------------------------------------
            if ( self % bufferLine .eq. OB_BUFFER_SIZE ) then
                if (.not. self % firstWrite .and. self % interpolate) then
                    do i =1, self % numberOfObservers
                        call self % observers(i) % interpolateSol(self % t, self % bufferLine)
                    end do
                    self % firstWrite = .TRUE.
                end if
                do i =1, self%numberOfObservers
                    call self % observers(i) % writeToFile(self % iter, self % t, self % bufferLine)
                end do
!               Reset buffer
!               ------------
                self % bufferLine = 0
            end if
        end if

    End Subroutine FWHWriteToFile

    Subroutine FWHDestruct(self)

        implicit none
        class(FWHClass), intent(inout)                   :: self

!       ---------------
!       Local variables
!       ---------------
        integer                                          :: i 

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

        safedeallocate(self % iter)
        safedeallocate(self % t)
        ! safedeallocate(self % sourceZone)
        do i = 1, self % numberOfObservers
            call self % observers(i) % destruct
        end do
        safedeallocate(self % observers)

    End Subroutine FWHDestruct
!
!//////////////////////////////////////////////////////////////////////////////
!
!        Auxiliars
!
!//////////////////////////////////////////////////////////////////////////////
!
    Function getNoOfObservers() result(no_of_observers)
      use ParamfileRegions
      implicit none
      integer                        :: no_of_observers
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
      no_of_observers = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find the observers
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
            error stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineacousticobserver' ) .gt. 0 ) then
               no_of_observers = no_of_observers + 1

            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

    End Function getNoOfObservers

End Module FWHGeneralClass