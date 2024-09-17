module RealDataLinkedList
   use SMConstants
   implicit none

   private
   public   RealDataLinkedList_t

    type RealDataLinkedList_t
      type(RealData_t), pointer :: head => NULL()
      integer                   :: no_of_entries = 0
      contains
         procedure   :: Append   => RealDataLinkedList_Append
         procedure   :: Load     => RealDataLinkedList_Load
         procedure   :: check    => CheckReal
         procedure   :: Destruct => RealDataLinkedList_Destruct
    end type RealDataLinkedList_t
   
   type RealData_t
      real(kind=RP)             :: value
      type(RealData_t), pointer :: next => NULL()
    end type RealData_t
!
!   ========
    contains 
!   ========
!
      subroutine RealDataLinkedList_Append( self , value ) 
         implicit none
         class(RealDataLinkedList_t) :: self
         real(kind=RP)               :: value
         type(RealData_t), pointer   :: current

         if ( self % no_of_entries .eq. 0 ) then
            allocate( self % head ) 
            self % head % value = value
            self % no_of_entries = 1

         else
            current => self % head
            
            do while ( associated(current % next) ) 
               current => current % next
            end do

            allocate(current % next)
            current % next % value = value
            self % no_of_entries = self % no_of_entries + 1 

         end if

         nullify(current)

      end subroutine RealDataLinkedList_Append

      subroutine RealDataLinkedList_Load( self , array ) 
         implicit none
         class(RealDataLinkedList_t) :: self
         real(kind=RP), allocatable  :: array(:)
         type(RealData_t), pointer   :: current
         integer                     :: i

         allocate( array( self % no_of_entries ) ) 

         current => self % head

         do i = 1 , self % no_of_entries
            array(i) = current % value
            current => current % next
         end do

         nullify(current)

      end subroutine RealDataLinkedList_Load
      
      logical function CheckReal( self, value ) result( found ) 
         use Utilities
         implicit none
         !-arguments-------------------------------------------------
         class(RealDataLinkedList_t), intent(inout) :: self
         real(kind=RP),               intent(in)    :: value
         !-local-variables-------------------------------------------
         type(RealData_t), pointer :: current => null()
      
         found = .false.
      
         if( .not. associated(self% head) ) return
      
         current => self% head

         if( AlmostEqual(current% value, value) ) then
            found = .true.
            nullify(current)
            return
         end if
      
         do while( associated(current% next) )
            current => current% next
            if ( AlmostEqual(current% value, value)  ) then
               found = .true.
               nullify(current)
               return
            end if
         end do
      
         nullify(current)
      
      end function CheckReal
      
      elemental subroutine RealDataLinkedList_Destruct( self )
         implicit none
         !-arguments------------------------------------------------
         class(RealDataLinkedList_t), intent(inout) :: self
         type(RealData_t), pointer                  :: data, nextdata
         integer                                    :: i
         !----------------------------------------------------------
         
         data => self % head
         do i = 1, self % no_of_entries
            nextdata => data % next
            deallocate(data)
            data => nextdata
         end do

         self % no_of_entries = 0

         nullify(data)
         nullify(nextdata)
         
      end subroutine RealDataLinkedList_Destruct

end module RealDataLinkedList