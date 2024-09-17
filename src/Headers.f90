!
!/////////////////////////////////////////////////////////////////////////////////////
!
!  Set of routines to uniformize the screen outputs.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
MODULE Headers
   use SMConstants
   implicit none

        INTEGER,        SAVE     :: iter
        INTEGER,        SAVE     :: loop_size
        INTEGER,        SAVE     :: no_of_points
        INTEGER,        SAVE     :: prev

        PRIVATE
        PUBLIC  Main_header, Ruler_Reset, Ruler_Update, Section_Header
        PUBLIC  SubSection_Header, Ruler_Header_Reset , Ruler_Header_Update
!
!       ========
        CONTAINS
!       ========
!       
        SUBROUTINE Main_header(title,date,time)
                use MPI_Process_Info
                IMPLICIT NONE
                CHARACTER(LEN = *)      :: title
                CHARACTER(LEN = *)      :: date
                CHARACTER(LEN = *)      :: time
                INTEGER, PARAMETER   :: siz = 100
                CHARACTER(LEN = siz)    :: ast1,ast2,astTitle
                INTEGER                 :: i

                if ( .not. MPI_Process % isRoot ) return

                ast1 = ''
                DO i = 1 , siz
                        ast1 = TRIM(ast1) // '#'
                END DO
                ast2 = ''
                ast2(1:1) = '#'
                ast2(siz:siz) = '#'

                astTitle = ''
                astTitle(1:1) = '#'
                astTitle(siz/2-LEN_TRIM(title)/2:siz/2 - LEN_TRIM(title)/2 + LEN_TRIM(title)-1) = TRIM(title)
                astTitle(siz:siz) = '#'

                WRITE(*,'(A)') ast1
                WRITE(*,'(A)') ast2
                WRITE(*,'(A)') ast2
                WRITE(*,'(A)') astTitle
                WRITE(*,'(A)') ast2
                WRITE(*,'(A)') ast2
                WRITE(*,'(A)') ast1
                WRITE(*,'(A,A,A,A)') "Compiled ",trim(date), ", ",trim(time)

        END SUBROUTINE Main_header
     
        SUBROUTINE Section_header(title)
           use MPI_Process_Info
           IMPLICIT NONE
           CHARACTER(LEN=*)     :: title
           INTEGER, PARAMETER   :: siz = 86
           CHARACTER(LEN=siz)      :: dotted_title
           INTEGER        :: i

           if ( .not. MPI_Process % isRoot ) return
     
           dotted_title = ''
           
           dotted_title(1:LEN_TRIM(title)) = TRIM(title)
     
           DO i = 1 , siz - LEN_TRIM(title)
              dotted_title = TRIM(dotted_title) // '.'
           END DO   
           WRITE(*,'(5X,A,A)') "\\\\ ",TRIM(dotted_title)
           FLUSH(6)
           
        END SUBROUTINE Section_header
      
        SUBROUTINE SubSection_header(title)
           use MPI_Process_Info
           IMPLICIT NONE
           CHARACTER(LEN=*)     :: title
           INTEGER        :: i
           character(len=len_trim(title))  :: auxstring
     
           if ( .not. MPI_Process % isRoot ) return

           do i = 1 , len_trim(title)
            auxstring(i:i) = "-"
           end do
           WRITE(STD_OUT,'(15X,A)') TRIM(title)
           write(STD_OUT,'(15X,A)') trim(auxstring)
           FLUSH(6)
           
        END SUBROUTINE SubSection_header
   
        SUBROUTINE Ruler_Header_Reset(title , loop_in)
                IMPLICIT NONE
                character(len=*)          :: title
                integer                   :: loop_in
                integer, parameter        :: siz = 86
      
                iter          = 0
                loop_size     = loop_in
                no_of_points  = siz - len_trim(title)
                prev          = 0
                
                call flush(STD_OUT)               
                write(* , '(/,5X,A,A)',advance = "no") "\\\\ ",trim(title)
                call flush(STD_OUT)               

        end subroutine Ruler_Header_Reset


        SUBROUTINE Ruler_Header_Update()
                IMPLICIT NONE
        
                iter = iter + 1     
                IF(FLOOR(1.0d0*iter*no_of_points/loop_size) .GT. prev) THEN
                        WRITE(*,'(A)',ADVANCE='NO') '.'
                        prev = FLOOR(1.0d0*iter*no_of_points/loop_size)
                        call flush(6)
                END IF
                
                
        END SUBROUTINE Ruler_Header_Update
                
 
        SUBROUTINE Ruler_Reset(title,n_in,loop_in)
                IMPLICIT NONE
                CHARACTER(LEN=*)  :: title
                INTEGER                 :: loop_in
                INTEGER                 :: n_in

                WRITE(* , '(20X,A)' , ADVANCE = 'NO' ) TRIM(title)
                FLUSH(6)
                iter = 0
                loop_size = loop_in
                no_of_points    = n_in - LEN_TRIM(title)
                prev    = 0
        END SUBROUTINE Ruler_Reset

        SUBROUTINE Ruler_Update()
                IMPLICIT NONE
        
                iter = iter + 1     
                IF(FLOOR(1.0d0*iter*no_of_points/loop_size) .GT. prev) THEN
                        WRITE(*,'(A)',ADVANCE='NO') '.'
                        prev = FLOOR(1.0d0*iter*no_of_points/loop_size)
                        FLUSH(6)
                END IF
                
                
        END SUBROUTINE Ruler_Update
                
                
END MODULE Headers

