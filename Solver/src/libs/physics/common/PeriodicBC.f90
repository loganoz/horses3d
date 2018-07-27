!
!//////////////////////////////////////////////////////
!
!   @File:    PeriodicBC.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jul 25 15:26:43 2018
!   @Last revision date: Thu Jul 26 22:00:46 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fb773e7c8706f4b4ef1f5bf9693a2b44f6c12dd2
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module PeriodicBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use Utilities, only: toLower, almostEqual
   use FluidData
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
!
!  ****************
!  Public variables
!  ****************
!
!
!  ******************
!  Public definitions
!  ******************
!
   public PeriodicBC_t
!
!  ****************************
!  Static variables definitions
!  ****************************
!
!
!  ****************
!  Class definition
!  ****************
!
   type, extends(GenericBC_t) ::  PeriodicBC_t
      character(len=LINE_LENGTH) :: associatedbname   
      contains
         procedure         :: Destruct          => PeriodicBC_Destruct
         procedure         :: GetPeriodicPair   => PeriodicBC_GetPeriodicPair
   end type PeriodicBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface PeriodicBC_t
      module procedure ConstructPeriodicBC
   end interface PeriodicBC_t
!
!  *******************
!  Function prototypes
!  *******************
!
   abstract interface
   end interface
!
!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      function ConstructPeriodicBC(bname)
!
!        ********************************************************************
!        · Definition of the periodic boundary condition in the control file:
!              #define boundary bname
!                 type             = periodic
!                 coupled boundary = bname
!              #end
!        ********************************************************************
!
         implicit none
         type(PeriodicBC_t)             :: ConstructPeriodicBC
         character(len=*), intent(in) :: bname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: fid, io
         character(len=LINE_LENGTH) :: boundaryHeader
         character(len=LINE_LENGTH) :: currentLine
         character(len=LINE_LENGTH) :: keyword, keyval
         logical                    :: inside
         type(FTValueDIctionary)    :: bcdict
         interface
            subroutine PreprocessInputLine(line)
               implicit none
               character(len=*), intent(inout) :: line
            end subroutine PreprocessInputLine
         end interface

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructPeriodicBC % BCType = "periodic"
         ConstructPeriodicBC % bname  = bname
         call ToLower(ConstructPeriodicBC % bname)

         write(boundaryHeader,'(A,A)') "#define boundary ",trim(bname)
         call toLower(boundaryHeader)
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do 
            read(fid, '(A)', iostat=io) currentLine

            IF(io .ne. 0 ) EXIT

            call PreprocessInputLine(currentLine)
            call toLower(currentLine)

            if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructPeriodicBC % bname)) 
            end if
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit


               keyword  = ADJUSTL(GetKeyword(currentLine))
               keyval   = ADJUSTL(GetValueAsString(currentLine))
               call ToLower(keyword)
      
               call bcdict % AddValueForKey(keyval, trim(keyword))

            end if

         end do
!
!        Analyze the gathered data
!        -------------------------
         if ( .not. bcdict % ContainsKey("coupled boundary") ) then
            print*, 'Select a "coupled boundary" for boundary',trim(bname)
            errorMessage(STD_OUT)
         else
            ConstructPeriodicBC % associatedbname = bcdict % StringValueForKey("coupled boundary", LINE_LENGTH)
         end if
         close(fid)
         call bcdict % Destruct
   
      end function ConstructPeriodicBC

      subroutine PeriodicBC_GetPeriodicPair(self, bname)
         implicit none
         class(PeriodicBC_t), intent(in)  :: self
         character(len=*), intent(out) :: bname

         bname = self % associatedbname
      end subroutine PeriodicBC_GetPeriodicPair
!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine PeriodicBC_Destruct(self)
         implicit none
         class(PeriodicBC_t)    :: self

      end subroutine PeriodicBC_Destruct
end module PeriodicBCClass
