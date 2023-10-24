module IntegerDataLinkedList
   implicit none

   private
   public   IntegerDataLinkedList_t

    type IntegerDataLinkedList_t
      logical                      :: allowRepetitions
      type(IntegerData_t), pointer :: head => NULL()
      integer                      :: no_of_entries = 0
      contains
         procedure   :: Add           => IntegerDataLinkedList_Add
         procedure   :: ExportToArray => IntegerDataLinkedList_ExportToArray
         procedure   :: check         => CheckInteger
         procedure   :: Destruct      => IntegerDataLinkedList_Destruct
    end type IntegerDataLinkedList_t
   
    type IntegerData_t
      integer                      :: value
      type(IntegerData_t), pointer :: next => NULL()
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
         logical, intent(in), optional :: allowRepetitions
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
         type(IntegerData_t), pointer   :: current
         type(IntegerData_t), pointer   :: previous
         type(IntegerData_t), pointer   :: newdata
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

         nullify(current)
         nullify(previous)
         nullify(newdata)

      end subroutine IntegerDataLinkedList_Add

      subroutine IntegerDataLinkedList_ExportToArray( self , array, sorted ) 
         use Utilities, only: QSort
         implicit none
         class(IntegerDataLinkedList_t) :: self
         integer, allocatable           :: array(:)
         logical, optional              :: sorted
         type(IntegerData_t), pointer   :: current
         integer                        :: i

         allocate( array( self % no_of_entries ) ) 

         current => self % head

         do i = 1 , self % no_of_entries
            array(i) = current % value
            current => current % next
         end do

         if ( present(sorted) ) then
            if ( sorted ) call QSort(array)
         end if

         nullify(current)

      end subroutine IntegerDataLinkedList_ExportToArray
      
      
     logical function CheckInteger( self, value ) result( found ) 
   
         implicit none
         !-arguments-------------------------------------------------
         class(IntegerDataLinkedList_t), intent(inout) :: self
         integer,                        intent(in)    :: value
         !-local-variables-------------------------------------------
         type(IntegerData_t), pointer :: current => null()
      
         found = .false.
      
         if( .not. associated(self% head) ) return
      
         current => self% head

         if( current% value .eq. value ) then
            found = .true.
            nullify(current)
            return
         end if
      
         do while( associated(current% next) )
            current => current% next
            if ( current% value .eq. value  ) then
               found = .true.
               nullify(current)
               return
            end if
         end do
       
         nullify(current)
      
      end function CheckInteger

      elemental subroutine IntegerDataLinkedList_Destruct(self)
         implicit none
         class(IntegerDataLinkedList_t), intent(inout) :: self
         type(IntegerData_t), pointer                  :: data, nextdata
         integer                                       :: i

         data => self % head
         do i = 1, self % no_of_entries
            nextdata => data % next

            deallocate(data)
            data => nextdata
         end do
         
         self % no_of_entries = 0

         nullify(data)
         nullify(nextdata)

      end subroutine IntegerDataLinkedList_Destruct

end module IntegerDataLinkedList