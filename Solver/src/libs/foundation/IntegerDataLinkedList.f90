!
!//////////////////////////////////////////////////////
!
!   @File:    IntegerDataLinkedList.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 13:29:58 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module IntegerDataLinkedList
   implicit none

   private
   public   IntegerDataLinkedList_t

    type IntegerDataLinkedList_t
      logical                       :: allowRepetitions
      class(IntegerData_t), pointer :: head => NULL()
      integer                       :: no_of_entries = 0
      contains
         procedure   :: Add        => IntegerDataLinkedList_Add
         procedure   :: ExportToArray => IntegerDataLinkedList_ExportToArray
         procedure   :: Destruct      => IntegerDataLinkedList_Destruct
    end type IntegerDataLinkedList_t
   
    type IntegerData_t
      integer                          :: value
      class(IntegerData_t), pointer    :: next => NULL()
    end type IntegerData_t

    interface IntegerDataLinkedList_t
      module procedure  ConstructIntegerDataLinkedList
    end interface 
!
!   ========
    contains 
!   ========
!
      function ConstructIntegerDataLinkedList(allowRepetitions)
         implicit none
         logical, intent(in), optional     :: allowRepetitions
         type(IntegerDataLinkedList_t) :: ConstructIntegerDataLinkedList 

         ConstructIntegerDataLinkedList % head => NULL()
         ConstructIntegerDataLinkedList % no_of_entries = 0
         
         if ( present(allowRepetitions) ) then
            ConstructIntegerDataLinkedList % allowRepetitions = allowRepetitions
         else
            ConstructIntegerDataLinkedList % allowRepetitions = .false.
         end if

      end function ConstructIntegerDataLinkedList

      subroutine IntegerDataLinkedList_Add( self , value ) 
         implicit none
         class(IntegerDataLinkedList_t) :: self
         integer                        :: value
         class(IntegerData_t), pointer  :: current
         class(IntegerData_t), pointer  :: previous
         class(IntegerData_t), pointer  :: newdata
         logical                        :: insert


         if ( self % no_of_entries .eq. 0 ) then
            allocate( self % head ) 
            self % head % value = value
            self % no_of_entries = 1

         else
            current => self % head
            if ( (.not. self % allowRepetitions) .and. &
                 (current % value .eq. value)   ) return
         
            do while ( associated(current % next) ) 
               current => current % next
               if ( (.not. self % allowRepetitions) .and. &
                    (current % value .eq. value)   ) return

            end do

            allocate(current % next)
            current % next % value = value
            self % no_of_entries = self % no_of_entries + 1 
         end if

      end subroutine IntegerDataLinkedList_Add

      subroutine IntegerDataLinkedList_ExportToArray( self , array, sorted ) 
         use SortingModule
         implicit none
         class(IntegerDataLinkedList_t) :: self
         integer, allocatable           :: array(:)
         logical, optional              :: sorted
         class(IntegerData_t), pointer       :: current
         integer                          :: i

         allocate( array( self % no_of_entries ) ) 

         current => self % head

         do i = 1 , self % no_of_entries
            array(i) = current % value
            current => current % next
         end do

         if ( present(sorted) ) then
            if ( sorted ) call QSort(array)
         end if

      end subroutine IntegerDataLinkedList_ExportToArray

      subroutine IntegerDataLinkedList_Destruct(self)
         implicit none
         class(IntegerDataLinkedList_t)   :: self
         class(IntegerData_t), pointer    :: data, nextdata
         integer     :: i

         data => self % head
         do i = 1, self % no_of_entries
            nextdata => data % next

            deallocate(data)
            data => nextdata
         end do
         
         self % no_of_entries = 0

      end subroutine IntegerDataLinkedList_Destruct

end module IntegerDataLinkedList