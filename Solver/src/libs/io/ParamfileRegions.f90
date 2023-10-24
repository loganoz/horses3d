module ParamfileRegions
   use SMConstants
   use Utilities, only: toLower
   implicit none
!
!  ******************
!  File configuration
!  ******************
!
   private
   public  readValueInRegion , getSquashedLine
   character, parameter       :: comment = '!'
   character, parameter       :: equal(2) = ['=',':'] 
   integer,   parameter       :: STR_LEN_PARAM = 512

   interface readValueInRegion
      module procedure readCharacterValueInRegion , readIntegerValueInRegion , readRealValueInRegion , readLogicalValueInRegion
   end interface readValueInRegion
!
   
!
!  ========
   contains
!  ========
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine readCharacterValueInRegion(fileName , label , var , in_label , out_label)
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         character(len=*)             :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
!        -------------------------------------------------------------
         integer               :: fID
         character(len=STR_LEN_PARAM)            :: auxstr
         integer               :: io
         integer               :: position
         logical               :: inside
         
         inside = .false.
!
!        Open file
!        ---------
         open ( newunit = fID , file = trim(fileName) , status = 'old' , action = 'read' ) 

         do 
            read( fID , '(A)' , iostat = io ) auxstr
            call toLower(auxstr)
            if (io .lt. 0) then
               var = ""
               exit

            elseif (io .gt. 0) then
               print*, "IOSTAT returned a positive number"
               error stop "Stopped."
         
            else
               
!
!              Check if inside a zone
!              ----------------------
               if ( getSquashedLine(auxstr) .eq. getSquashedLine(in_label) ) then
                  inside = .true.
                  cycle
               elseif( getSquashedLine(auxstr) .eq. getSquashedLine(out_label) ) then
                  inside = .false.
                  cycle
               end if

!              Removed commented part of the string
!              ------------------------------------
               if (inside) then
                  position = index(trim(auxstr) , comment)
                  if ( position .gt. 0 ) then
                     auxstr = auxstr(1:position-1)
                  end if
!   
!                 Look for the label
!                 ------------------
                  position = index(trim(auxstr) , trim(label) )
                  if ( position .eq. 0) then
                     cycle
                  end if
!   
!                 The label is present. Check whether the equal is present
!                 --------------------------------------------------------
                  position = max(index(trim(auxstr) , equal(1)) , index(trim(auxstr) , equal(2)) )
                  if ( position .eq. 0 ) then
                     cycle
                  end if

                  if ( getSquashedLine(auxstr(1:position-1)) .eq. getSquashedLine(label) ) then
!
!                    The label matches the one present in the string
!                    -----------------------------------------------
                     auxstr = auxstr(position+1:)
                  else
                     cycle
                  end if
!   
!                 Get the value
!                 -------------
                  var = trim(adjustl(auxstr))

                  exit
               end if

            end if
         end do

!
!        Close file
!        ----------
         close ( fID ) 
      

      end subroutine readCharacterValueInRegion

      subroutine readIntegerValueInRegion(fileName , label , var , in_label , out_label )
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         integer, allocatable         :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
         character(len=STR_LEN_PARAM) :: auxstr
         integer                      :: io

         call readCharacterValueInRegion(fileName , label, auxstr , in_label , out_label)

         auxstr = trim(adjustl(auxstr))
         
         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if
 

      end subroutine readIntegerValueInRegion

      subroutine readRealValueInRegion(fileName , label , var , in_label , out_label)
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         real(kind=RP), allocatable   :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
         character(len=STR_LEN_PARAM) :: auxstr
         integer                      :: io
      
         call readCharacterValueInRegion(fileName , label, auxstr , in_label , out_label)

         auxstr = trim(adjustl(auxstr))
         
         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if
      
      end subroutine readRealValueInRegion

      subroutine readLogicalValueInRegion(fileName , label , var , in_label , out_label)
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         logical, allocatable         :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
         character(len=STR_LEN_PARAM) :: auxstr
         integer                      :: io
      
         call readCharacterValueInRegion(fileName , label, auxstr , in_label , out_label)

         auxstr = trim(adjustl(auxstr))

         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if

      end subroutine readLogicalValueInRegion

      function getSquashedLine(line) result (squashed)
         implicit none
         character(len=*), intent(in)     :: line
         character(len=STR_LEN_PARAM)     :: squashed
         integer                          :: pos
!
!        First remove comments
!        ---------------------
         pos = index(trim(line) , comment)
         if ( pos .gt. 0 ) then
            squashed = trim(adjustl(line(1:pos-1)))

         else
            squashed = trim(adjustl(line))
         end if

         pos = index(trim(adjustl(squashed)),' ')

         if (pos .eq. 0) then
            squashed = trim(adjustl(squashed))
         else
            squashed = squashed(1:pos-1) // squashed(pos+1:)
         end if

         do
            pos = index(trim(adjustl(squashed)),' ')

            if (pos .eq. 0) then
               squashed = trim(adjustl(squashed))
               return
            else
               squashed = squashed(1:pos-1) // squashed(pos+1:)
            end if
         end do
            
         write(squashed,'(A)')  trim(adjustl(squashed))

      end function getSquashedLine

end module ParamfileRegions