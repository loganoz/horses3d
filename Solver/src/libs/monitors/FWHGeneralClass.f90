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
!This class represents the general behaviour of the Fwoc Williams and Hawckings aero accoustic analogy

#include "Includes.h"
Module FWHGeneralClass  !

    use SMConstants
    ! use MonitorDefinitions
    use FWHDefinitions, only: OB_BUFFER_SIZE_DEFAULT, OB_BUFFER_SIZE_DEFAULT, STR_LEN_OBSERVER
    use FWHObseverClass
    use HexMeshClass
    use ZoneClass
    use FileReadingUtilities      , only: getFileName
    Implicit None

!
!   *****************************
!   Main FWH class definition
!   *****************************
    type FWHClass

        character(len=LINE_LENGTH)                                        :: solution_file
        integer                                                           :: numberOfObservers
        integer                                                           :: bufferLine
        real(kind=RP)                                                     :: dt_update
        integer, dimension(:), allocatable                                :: iter
        real(kind=RP), dimension(:), allocatable                          :: t
        class(ObserverClass), dimension(:), allocatable                   :: observers
        class(Zone_t), allocatable                                        :: sourceZone
        logical                                                           :: isSolid
        logical                                                           :: isActive
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
        use Monopole
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
        integer                                             :: no_of_zones, no_of_face_i
        integer, dimension(:), allocatable                  :: facesIDs, faces_per_zone, zonesIDs
        logical, save                                       :: FirstCall = .TRUE.
        character(len=LINE_LENGTH)                          :: zones_str
        character(len=LINE_LENGTH), allocatable             :: zones_names(:)

!        look if the accoustic analogy calculations are set to be computed
!        --------------------------------
        !TODO read accoustic analogy type and return if is not defined, check for FWH if is defined and not FWH stop and send error
        if (.not. controlVariables % containsKey("accoustic analogy")) then
            self % isActive = .FALSE.
            print *, "FWH not activated"
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
!       --------------------------
        self % isSolid   = .not. controlVariables % logicalValueForKey("accoustic analogy permable")
        ! if (self % isSolid) then
            if (controlVariables % containsKey("accoustic solid surface")) then
                zones_str = controlVariables%stringValueForKey("accoustic solid surface", LINE_LENGTH)
            else 
                stop "Accoustic surface for integration is not defined"
            end if
            call getCharArrayFromString(zones_str, LINE_LENGTH, zones_names)

    !       Get the zones ids fo the mesh and for each the number of faces
    !       --------------------------
            no_of_zones = size(zones_names)
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
            allocate( facesIDs(SUM(faces_per_zone)) )
            no_of_face_i = 1
            do i = 1, no_of_zones
                if (zonesIDs(i) .eq. -1) cycle
                facesIDs(no_of_face_i:no_of_face_i+faces_per_zone(i)-1) = mesh % zones(zonesIDs(i)) % faces
                no_of_face_i = no_of_face_i + faces_per_zone(i) 
            end do 
        ! else
        !     stop "Permeable surfaces not implemented yet"
        ! end if

        ! create self sourceZone using facesIDs
!       --------------------------
        allocate( self % sourceZone )
        call self % sourceZone % CreateFicticious(-1, "FW_Surface", SUM(faces_per_zone), facesIDs)

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
        self % interpolate = .TRUE.
        ! self % interpolate = .FALSE.
!       Initialize observers
!       ----------
        call getMeanStreamValues()
        allocate( self % observers(self % numberOfObservers) )
        do i = 1, self % numberOfObservers
            call self % observers(i) % construct(self % sourceZone, mesh, i, self % solution_file, FirstCall, self % interpolate)
        end do 

        self % bufferLine = 0
        self % firstWrite = .FALSE.
        
        FirstCall = .FALSE.

!        Describe the zones
!        ------------------
         if ( .not. MPI_Process % isRoot ) return
         call Subsection_Header("Ficticious FWH zone")
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of faces: ", self % sourceZone % no_of_faces
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of observers: ", self % numberOfObservers
         write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of integrals: ", self % numberOfObservers * self % sourceZone % no_of_faces

        call setVals()

    End Subroutine FWHConstruct

    Subroutine FWHUpate(self, mesh, t, iter, dt, isFirst)

        implicit none

        class(FWHClass)                                     :: self
        class(HexMesh)                                      :: mesh
        real(kind=RP), intent(in)                           :: t, dt
        integer, intent(in)                                 :: iter
        logical, intent(in), optional                        :: isFirst

!       ---------------
!       Local variables
!       ---------------
!
        integer                                             :: i 

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

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
        call SourceProlongSolution(self % sourceZone, mesh, t)
!       see if its regular or interpolated
!       -----------------------
        if (.not. self % firstWrite) then
            do i = 1, self % numberOfObservers
                ! call self % observers(i) % update(mesh, self % bufferLine, self % isSolid, dt, isFirst)
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
