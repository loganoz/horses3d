!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
module RealDataLinkedList
   use SMConstants
   implicit none

   private
   public   RealDataLinkedList_t

    type RealDataLinkedList_t
      class(RealData_t), pointer    :: head => NULL()
      integer                       :: no_of_entries = 0
      contains
         procedure   :: Append => RealDataLinkedList_Append
         procedure   :: Load   => RealDataLinkedList_Load
    end type RealDataLinkedList_t
   
    type RealData_t
      real(kind=RP)  :: value
      class(RealData_t), pointer    :: next => NULL()
    end type RealData_t
!
!   ========
    contains 
!   ========
!
      subroutine RealDataLinkedList_Append( self , value ) 
         implicit none
         class(RealDataLinkedList_t)      :: self
         real(kind=RP)                    :: value
         class(RealData_t), pointer       :: current

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

      end subroutine RealDataLinkedList_Append

      subroutine RealDataLinkedList_Load( self , array ) 
         implicit none
         class(RealDataLinkedList_t)      :: self
         real(kind=RP), allocatable       :: array(:)
         class(RealData_t), pointer       :: current
         integer                          :: i

         allocate( array( self % no_of_entries ) ) 

         current => self % head

         do i = 1 , self % no_of_entries
            array(i) = current % value

            current => current % next
         end do

      end subroutine RealDataLinkedList_Load

end module RealDataLinkedList
