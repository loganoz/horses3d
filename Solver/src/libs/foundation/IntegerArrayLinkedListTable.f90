module IntegerArrayLinkedListTable
   use SMConstants, only: STD_OUT
   use Utilities,   only: Qsort
   implicit none

   private
   public   Table_t

   integer, parameter   :: DATA_SIZE = 3

   type Entry_t
      integer                 :: global_pos
      integer                 :: val(DATA_SIZE)
      class(Entry_t), pointer :: next => NULL()
   end type Entry_t

   type LinkedList_t
      integer                    :: no_of_entries = 0
      class(Entry_t), pointer    :: head => NULL()
      contains
         procedure :: AddEntry      => LinkedList_AddEntry
         procedure :: ContainsEntry => LinkedList_ContainsEntry
         procedure :: Describe      => LinkedList_Describe
         procedure :: Destruct      => LinkedList_Destruct
   end type LinkedList_t

   type Table_t
      integer                         :: no_of_lists
      integer                         :: no_of_entries = 0
      type(LinkedList_t), allocatable :: lists(:)
      contains 
         procedure :: AddEntry      => Table_AddEntry
         procedure :: ContainsEntry => Table_ContainsEntry
         procedure :: Destruct      => Table_Destruct
         procedure :: Describe      => Table_Describe
   end type Table_t

   interface Table_t
      module procedure ConstructTableWithSize
   end interface Table_t

   interface LinkedList_t
      module procedure ConstructLinkedList
   end interface LinkedList_t

   contains
!
!/////////////////////////////////////////////////////////////////
!
!        TABLE PROCEDURES
!        ----------------
!
!/////////////////////////////////////////////////////////////////
!
      function ConstructTableWithSize(t_size)
         implicit none
         integer, intent(in)  :: t_size
         type(Table_t)        :: ConstructTableWithSize
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: t_id
   
         ConstructTableWithSize % no_of_lists = t_size
!
!        Allocate the table
!        ------------------
         allocate(ConstructTableWithSize % lists(t_size))
!
!        Initialize the linked lists
!        ---------------------------
         do t_id = 1 , t_size
            ConstructTableWithSize % lists(t_id) = LinkedList_t()
         end do

      end function ConstructTableWithSize

      subroutine Table_AddEntry(self, val)
         implicit none
         class(Table_t), intent(inout)    :: self
         integer,        intent(in)       :: val(DATA_SIZE+1)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: val_ordered(DATA_SIZE+1)
         integer  :: global_position
         integer  :: t_id
!
!        Order value
!        -----------
         val_ordered = val
         call QSort(val_ordered)
         
         self % no_of_entries = self % no_of_entries + 1 
         global_position = self % no_of_entries

         call self % lists(val_ordered(1)) % AddEntry(val_ordered(2:DATA_SIZE+1), global_position)
   
      end subroutine Table_AddEntry

      integer function Table_ContainsEntry(self, val)
         implicit none
         class(Table_t), intent(inout)    :: self
         integer,        intent(in)       :: val(DATA_SIZE+1)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: val_ordered(DATA_SIZE+1)
!
!        Order value
!        -----------
         val_ordered = val
         call QSort(val_ordered)

         Table_ContainsEntry =  self % lists(val_ordered(1)) % ContainsEntry(val_ordered(2:DATA_SIZE+1))
   
      end function Table_ContainsEntry
   
      subroutine Table_Destruct(self)
         implicit none
         class(Table_t), intent(inout)    :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: t_pos

         do t_pos = 1, self % no_of_lists
            call self % lists(t_pos) % Destruct
         end do

         deallocate(self % lists)

      end subroutine Table_Destruct

      subroutine Table_Describe(self)
         implicit none
         class(Table_t), intent(in)  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: t_id

         do t_id = 1, self % no_of_lists
            write(STD_OUT, '(A,A20,I0)') "->", "Entries in list: ", t_id
            call self % lists(t_id) % Describe
         end do

      end subroutine Table_Describe

!
!//////////////////////////////////////////////////////////////////
!
!        LINKED LIST PROCEDURES
!        ----------------------
!
!//////////////////////////////////////////////////////////////////
!
      function ConstructLinkedList()
         implicit none
         type(LinkedList_t)   :: ConstructLinkedList
!
!        ---------------
!        Local variables
!        ---------------
!
         
         ConstructLinkedList % no_of_entries = 0
         ConstructLinkedList % head => NULL()

         allocate(ConstructLinkedList % head)

      end function ConstructLinkedList
   
      subroutine LinkedList_AddEntry(self, new_val, global_position)
         implicit none
         class(LinkedList_t), intent(inout) :: self
         integer, intent(in)                :: new_val(DATA_SIZE)
         integer, intent(in)                :: global_position
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: e_id
         class(Entry_t),   pointer :: current 

         current => self % head

         do e_id = 1, self % no_of_entries
!
!           Check if is equal to any of the entries
!           ---------------------------------------
            if (all(current % val .eq. new_val) ) then
               return
            endif

            current => current % next
         end do            
!
!        Store the content in current (already allocated)
!        ------------------------------------------------
         current % val = new_val
         current % global_pos = global_position
         allocate(current % next)

         self % no_of_entries = self % no_of_entries + 1 

      end subroutine LinkedList_AddEntry

      integer function LinkedList_ContainsEntry(self, new_val)
         implicit none
         class(LinkedList_t), intent(inout) :: self
         integer, intent(in)                :: new_val(DATA_SIZE)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: e_id
         class(Entry_t),   pointer :: current 

         current => self % head

         do e_id = 1, self % no_of_entries
!
!           Check if is equal to any of the entries
!           ---------------------------------------
            if (all(current % val .eq. new_val) ) then
               LinkedList_ContainsEntry = current % global_pos
               return
            endif

            current => current % next
         end do            

         LinkedList_ContainsEntry = 0

      end function LinkedList_ContainsEntry

      subroutine LinkedList_Destruct(self)
         implicit none
         class(LinkedList_t), intent(inout) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                 :: e_id
         class(Entry_t), pointer :: current, prev

         current => self % head

         do e_id = 1, self % no_of_entries
            prev    => current
            current => current % next

            deallocate(prev)
         end do
!
!        Deallocate the entry that is allocated but unused
!        -------------------------------------------------
         deallocate(current)
         
         self % no_of_entries = 0

      end subroutine LinkedList_Destruct

      subroutine LinkedList_Describe(self)
         implicit none
         class(LinkedList_t), intent(in)  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: e_id, d_id
         class(Entry_t), pointer :: current

         current => self % head

         do e_id = 1, self % no_of_entries
            write(STD_OUT,'(20X, I0,A,I0)', advance = "no") e_id, ": [", current % val(1)
            do d_id = 1, DATA_SIZE-1
               write(STD_OUT,'(A,I0)', advance = "no") ", ", current % val(d_id+1)
            end do
            write(STD_OUT,'(A)') "]"

            current => current % next
         end do

      end subroutine LinkedList_Describe

end module IntegerArrayLinkedListTable