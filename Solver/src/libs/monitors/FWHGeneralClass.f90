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
    use MonitorDefinitions
    use FWHObseverClass
    use HexMeshClass
    use ZoneClass
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
        class(Zone_t)                                                     :: sourceZone
        logical                                                           :: isSolid

        contains

            procedure :: construct      => FWHConstruct
            procedure :: destruct       => FWHDestruct
            procedure :: updateValues   => FWHUpate
            procedure :: writeToFile    => FWHWriteToFile

    end type FWHClass
           ! se debe construir desde la clase general de FW, esta debe hacer algo similar a la de monitores: crear update, escribir,
           ! crear archivo de escritura, allocar, leer de control file, etc...
           ! debe crear la clase zona que tenga todas las caras de source (aglomerar todas las zone de BC dadas)

    contains

    Subroutine FWHConstruct(self, mesh, controlVariables)
        use FTValueDictionaryClass
        use mainKeywordsModule
        use FileReadingUtilities, only: getCharArrayFromString
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
        character(len=STR_LEN_MONITORS)                     :: line
        character(len=STR_LEN_MONITORS)                     :: solution_file
        integer, dimension(:), allocatable                  :: facesIDs
        logical, save                                       :: FirstCall = .TRUE.
        character(len=LINE_LENGTH)                          :: zones_str, zones_names(:)

        allocate( self % t(BUFFER_SIZE), self % iter(BUFFER_SIZE) )

!       Get the general configuration of control file
!       --------------------------
        !TODO read accoustic analogy type and return if is not defined, check for FWH if is defined and not FWH stop and send error
        self % isSolid   = .not. controlVariables % logicalValueForKey("accoustic analogy permable")
        if self % isSolid then
            if (controlVariables % containsKey("accoustic solid surface")) then
                zones_str = controlVariables%stringValueForKey("accoustic solid surface", LINE_LENGTH)
            else 
                stop "Accoustic surface for integration is not defined"
            end if
            call getCharArrayFromString(zones_str, LINE_LENGTH, zones_names)
            print *, zones_names
        else
            stop "Permeable surfaces not implemented yet"
        end if
! ------------------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------------------
        ! todo create self sourceZone using zones_names --------------------------------------
! ------------------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------------------

!       Get the solution file name
!       --------------------------
        solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_MONITORS )
!
!       Remove the *.hsol termination
!       -----------------------------
        solution_file = trim(getFileName(solution_file))
        self % solution_file = trim(solution_file)

!       Search in case file for probes, surface monitors, and volume monitors
!       ---------------------------------------------------------------------
        if (mesh % child) then ! Return doing nothing if this is a child mesh
           self % numberOfObservers = 0
        else
           self % numberOfObservers = getNoOfObservers()
        end if

!       Initialize
!       ----------
        allocate( self % observers(self % numberOfObservers) )
        do i = 1, self%numberOfObservers
        print *, i
            ! call self % observers(i) % construct(self % sourceZone, mesh, i, self % solution_file, FirstCall, self % isSolid)
        end do 

        self % bufferLine = 0
        
        FirstCall = .FALSE.

    End Subroutine FWHConstruct

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
!     Read the whole file to find monitors
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

            if ( index ( line , '#defineobserver' ) .gt. 0 ) then
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
