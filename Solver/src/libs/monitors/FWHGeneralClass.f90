!
!   @File:    ObserverClass.f90
!   @Author:  Oscar Marino (oscar.marino@upm.es)
!   @Created: Mar 25 2020
!   @Last revision date: 
!   @Last revision author: 
!   @Last revision commit: 
!
!//////////////////////////////////////////////////////
!
!This class represents the general behaviour of the Ffowcs Williams and Hawckings aero accoustic analogy

#include "Includes.h"
Module FWHGeneralClass  !

    use SMConstants
    use FWHDefinitions, only: OB_BUFFER_SIZE_DEFAULT, OB_BUFFER_SIZE_DEFAULT, STR_LEN_OBSERVER
    use FWHObseverClass
    use HexMeshClass
    use ZoneClass
    use FileReadingUtilities      , only: getFileName
    use AutosaveClass
    Implicit None

!
!   *****************************
!   Main FWH class definition
!   *****************************
    type FWHClass

        character(len=LINE_LENGTH)                                        :: solution_file
        integer                                                           :: numberOfObservers = 0
        integer                                                           :: bufferLine
        real(kind=RP)                                                     :: dt_update
        integer, dimension(:), allocatable                                :: iter
        real(kind=RP), dimension(:), allocatable                          :: t
        class(ObserverClass), dimension(:), allocatable                   :: observers
        class(Zone_t), allocatable                                        :: sourceZone
        integer                                                           :: totalNumberOfFaces
        integer, dimension(:), allocatable                                :: globalFid
        integer, dimension(:), allocatable                                :: faceOffset
        integer, dimension(:), allocatable                                :: elementSide
        logical                                                           :: isSolid
        logical                                                           :: isActive = .false.
        logical                                                           :: firstWrite
        logical                                                           :: interpolate
        logical                                                           :: saveSourceSolFile
        logical                                                           :: saveSourceMeshFile
        type(Autosave_t)                                                  :: autosave

        contains

            procedure :: construct      => FWHConstruct
            procedure :: destruct       => FWHDestruct
            procedure :: autosaveConfig => FWHSaveSolutionConfiguration
            procedure :: updateValues   => FWHUpate
            procedure :: writeToFile    => FWHWriteToFile
            procedure :: saveSourceSol  => FWHSaveSourceSolution
            procedure :: loadSourceSol  => FWHLoadSourceSolution
            procedure :: saveSourceMesh => FWHSaveSourceMesh

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
        ! integer, dimension(:), allocatable                  :: facesIDs, faces_per_zone, zonesIDs, eSides
        integer, dimension(:), allocatable                  :: facesIDs, faces_per_zone, zonesIDs
        logical, save                                       :: FirstCall = .TRUE.
        character(len=LINE_LENGTH)                          :: zones_str, zones_str2, surface_file
        character(len=LINE_LENGTH), allocatable             :: zones_names(:), zones_temp(:), zones_temp2(:)

!        look if the accoustic analogy calculations are set to be computed
!        --------------------------------
        !TODO read accoustic analogy type and return if is not defined, check for FWH if is defined and not FWH stop and send error
        if (.not. controlVariables % containsKey("accoustic analogy")) then
            self % isActive = .FALSE.
            ! print *, "FWH not activated"
            return
        end if

!        Setup the buffer
!        ----------------
         if (controlVariables % containsKey("observers flush interval") ) then
            OB_BUFFER_SIZE = controlVariables % integerValueForKey("observers flush interval")
         end if

        self % isActive = .TRUE.
        allocate( self % t(OB_BUFFER_SIZE), self % iter(OB_BUFFER_SIZE) )

!       Get the general configuration of control file
!       First get the surface as a zone
!       -------------------------------
        self % isSolid   = .not. controlVariables % logicalValueForKey("accoustic analogy permable")
        if (controlVariables % containsKey("accoustic solid surface")) then
            zones_str = controlVariables % stringValueForKey("accoustic solid surface", LINE_LENGTH)
            zones_str2 = controlVariables % stringValueForKey("accoustic solid surface cont", LINE_LENGTH)
            call toLower(zones_str)
            call toLower(zones_str2)
            call getCharArrayFromString(zones_str, LINE_LENGTH, zones_temp)
            if (zones_str2 .ne. "") then
                no_of_zones = size(zones_temp)
                call getCharArrayFromString(zones_str2, LINE_LENGTH, zones_temp2)
                no_of_zones = no_of_zones + size(zones_temp2)
                allocate(zones_names(no_of_zones))
                zones_names(1:size(zones_temp)) = zones_temp
                zones_names(size(zones_temp)+1:no_of_zones) = zones_temp2
                safedeallocate(zones_temp)
                safedeallocate(zones_temp2)
            else
                no_of_zones = size(zones_temp)
                allocate(zones_names(no_of_zones))
                zones_names = zones_temp
                safedeallocate(zones_temp)
            end if 

    !       Get the zones ids from the mesh and for each, the number of faces
    !       --------------------------
            allocate( faces_per_zone(no_of_zones), zonesIDs(no_of_zones) )
            do i = 1, no_of_zones
                zonesIDs(i) = getZoneID(zones_names(i), mesh)
                if (zonesIDs(i) .eq. -1) then
                    write(*,'(A,A,A)') "Warning: Accoustic surface ", trim(zones_names(i)), " not found in the mesh, will be ignored"
                    faces_per_zone(i) = 0
                else
                    faces_per_zone(i) = size(mesh % zones(zonesIDs(i)) % faces)
                end if
            end do 

    !       Get the faces Ids of all zones as a single array
    !       --------------------------
            allocate( facesIDs(SUM(faces_per_zone)) , self % elementSide(SUM(faces_per_zone)))
            no_of_face_i = 1
            do i = 1, no_of_zones
                if (zonesIDs(i) .eq. -1) cycle
                facesIDs(no_of_face_i:no_of_face_i+faces_per_zone(i)-1) = mesh % zones(zonesIDs(i)) % faces
                no_of_face_i = no_of_face_i + faces_per_zone(i) 
            end do 
            ! the side is always 1 since is at a face of a boundary
            ! eSides = 1
            self % elementSide = 1

            deallocate(zonesIDs, zones_names) 
        elseif (controlVariables % containsKey("accoustic surface file")) then
            allocate( faces_per_zone(1) )
            surface_file = controlVariables % stringValueForKey("accoustic surface file", LINE_LENGTH)
            call SourceLoadSurfaceFromFile(mesh, surface_file, facesIDs, faces_per_zone(1), self % elementSide)
        else
            stop "Accoustic surface for integration is not defined"
        end if

        ! create self sourceZone using facesIDs
!       --------------------------
        allocate( self % sourceZone )
        call self % sourceZone % CreateFicticious(-1, "FW_Surface", SUM(faces_per_zone), facesIDs)
        deallocate(facesIDs, faces_per_zone)

!       Gather the total number of faces
!       ------------------
        if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
            call mpi_allreduce(self % sourceZone % no_of_faces, no_of_faces, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            no_of_faces = self % sourceZone % no_of_faces
        end if
        self % totalNumberOfFaces = no_of_faces

!       Get arrays for I/O of surface solution file
!       --------------------------
        allocate(self % globalFid(self % sourceZone % no_of_faces))
        allocate(self % faceOffset(self % sourceZone % no_of_faces))
        call SourcePrepareForIO(self % sourceZone, mesh, no_of_faces, self % globalFid, self % faceOffset)

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

!       Set interpolate atribute as TRUE by default
        ! todo: read from constrol variables
        self % interpolate = .TRUE.
        ! self % interpolate = .FALSE.

!       Get whether the source surface solution and the mesh is saved in files
!       -----------------------------------------------------------------------
        self % saveSourceSolFile = controlVariables % logicalValueForKey("accoustic solution save")
        self % saveSourceMeshFile = controlVariables % logicalValueForKey("accoustic surface mesh save")

!       Initialize observers
!       --------------------
        call getMeanStreamValues()
        allocate( self % observers(self % numberOfObservers) )
        do i = 1, self % numberOfObservers
            call self % observers(i) % construct(self % sourceZone, mesh, i, self % solution_file, FirstCall, &
                                                 self % interpolate, no_of_faces, self % elementSide)
        end do 

        self % bufferLine = 0
        self % firstWrite = .FALSE.
        
        FirstCall = .FALSE.

!        Describe the zones
!        ------------------
         if ( .not. MPI_Process % isRoot ) return
         call Subsection_Header("Ficticious FWH zone")
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of faces: ", no_of_faces
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of observers: ", self % numberOfObservers
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of integrals: ", self % numberOfObservers * no_of_faces
         write(STD_OUT,'(30X,A,A28,L1)') "->", "Save zone solution: ", self % saveSourceSolFile

    End Subroutine FWHConstruct

    Subroutine FWHSaveSolutionConfiguration(self, controlVariables, t0)

        use FTValueDictionaryClass
        implicit none

        class(FWHClass)                                     :: self
        class(FTValueDictionary), intent(in)                :: controlVariables
        real(kind=RP), intent(in)                           :: t0

!       ---------------
!       Local variables
!       ---------------
        real(kind=RP)                                       :: dtSave

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) then
            self % autosave = Autosave_t(.FALSE.,.FALSE.,huge(1),huge(1.0_RP),nextAutosaveTime=huge(1.0_RP),mode=AUTOSAVE_UNDEFINED)
            return
        end if

!       configure save type, used for update, write and save file
!       --------------------------
        if (controlVariables % containsKey("accoustic save timestep")) then
            dtSave = controlVariables % doublePrecisionValueForKey("accoustic save timestep")
            self % autosave = Autosave_t(.FALSE.,.TRUE.,huge(1),dtSave,nextAutosaveTime=dtSave+t0,mode=AUTOSAVE_BY_TIME)
        else
            !Save each time step
            self % autosave = Autosave_t(.FALSE.,.TRUE.,1,huge(1.0_RP),nextAutosaveTime=huge(1.0_RP),mode=AUTOSAVE_BY_ITERATION)
        end if

    End Subroutine FWHSaveSolutionConfiguration

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
        if (prolong) call SourceProlongSolution(self % sourceZone, mesh, self % elementSide)

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
!           In this case the observers are exported to their files and the buffer is reseted
!           -------------------------------------------------------------------------------
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
        safedeallocate(self % sourceZone)
        do i = 1, self % numberOfObservers
            call self % observers(i) % destruct
        end do
        safedeallocate(self % observers)

    End Subroutine FWHDestruct

     SUBROUTINE FWHSaveSourceSolution(self, mesh, iter, time)
        IMPLICIT NONE
!
!       ------------------------------------
!       Save the results to the surface solution file
!       ------------------------------------
!
!       ----------------------------------------------
        class(FWHClass)                                     :: self
        class(HexMesh), intent(in)                          :: mesh
        integer, intent(in)                                 :: iter              !< Time step
        real(kind=RP), intent(in)                           :: time              !< Simu time

        !local variables
        character(len=LINE_LENGTH)                          :: FinalName      !  Final name for particular file
!       ----------------------------------------------
        
        if(.not. self % saveSourceSolFile .or. .not. self % isActive) return

        WRITE(FinalName,'(2A,I10.10,A)')  TRIM(self % solution_file),'_',iter,'.fwhs.hsol'
        call SourceSaveSolution(self % sourceZone, mesh, time, iter, FinalName, self % totalNumberOfFaces, self % globalFid, self % faceOffset)
     
     END SUBROUTINE FWHSaveSourceSolution

     Subroutine FWHLoadSourceSolution(self, fileName, mesh)

        implicit none
        class(FWHClass)                                      :: self
        character(len=*), intent(in)                         :: fileName
        class (HexMesh), intent(inout)                       :: mesh

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

        call SourceLoadSolution(self % sourceZone, mesh, fileName, self % globalFid, self % faceOffset, self % elementSide)

     End Subroutine FWHLoadSourceSolution

     Subroutine FWHSaveSourceMesh(self, mesh, initial_iteration)
        IMPLICIT NONE
!
!       ------------------------------------
!       Save the mesh to the surface mesh file
!       ------------------------------------
!
!       ----------------------------------------------
        class(FWHClass)                                     :: self
        class(HexMesh), intent(in)                          :: mesh
        integer, intent(in)                                 :: initial_iteration

        !local variables
        character(len=LINE_LENGTH)                          :: FinalName      !  Final name for particular file
!       ----------------------------------------------
        
        if(.not. self % saveSourceMeshFile .or. .not. self % isActive) return

        WRITE(FinalName,'(2A,I10.10)')  TRIM(self % solution_file),'_',initial_iteration
        call SourceSaveMesh(self % sourceZone, mesh, FinalName, self % totalNumberOfFaces, self % globalFid, self % faceOffset)

     End Subroutine FWHSaveSourceMesh
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
            stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineaccousticobserver' ) .gt. 0 ) then
               no_of_observers = no_of_observers + 1

            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

    End Function getNoOfObservers

    integer Function getZoneID(zone_name, mesh) result(n)

        character(len=*), intent(in)                        :: zone_name
        class(HexMesh), intent(in)                          :: mesh

        !local variables
        integer                                             :: zoneID

         n = -1
         do zoneID = 1, size(mesh % zones)
            if ( trim(mesh % zones(zoneID) % name) .eq. trim(zone_name) ) then
               n = zoneID
               exit
            end if
         end do

    End Function getZoneID

End Module FWHGeneralClass
