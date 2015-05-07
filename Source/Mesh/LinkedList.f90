!
!////////////////////////////////////////////////////////////////////////
!
!      ListAndHash.f95
!      Created: 2008-06-24 09:22:03 -0400 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      Module LinkedListClass
      USE ListDataClass
      IMPLICIT NONE 
!
!-------------------------------------------------------------------
! A basic linked list class
!-------------------------------------------------------------------
!
      TYPE LLRecord
         TYPE(LLRecord), POINTER :: next
         TYPE(ListData)          :: d
      END TYPE LLRecord
      
      TYPE LinkedList
         TYPE(LLRecord), POINTER :: head, tail, current
      END TYPE LinkedList
      
!      PRIVATE :: LLRecord
      
      INTERFACE Construct
         MODULE PROCEDURE ConstructList
      END INTERFACE Construct
      INTERFACE Destruct
         MODULE PROCEDURE DestructList
      END INTERFACE Destruct
      INTERFACE PRINT
         MODULE PROCEDURE PrintList
      END INTERFACE PRINT
!
!     ========
      CONTAINS 
!     ========
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructList( this )
         IMPLICIT NONE 
         TYPE(LinkedList) :: this
         NULLIFY(this%head)! => NULL()
         NULLIFY(this%tail)! => NULL()
         NULLIFY(this%current) !=> NULL()
      END SUBROUTINE ConstructList
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructList( this )
         IMPLICIT NONE 
         TYPE(LinkedList)        :: this
         TYPE(LLRecord), POINTER :: p
         this%current => this%head
         DO WHILE (ASSOCIATED(this%current) )
            p => this%current%next
            CALL Destruct(this%current%d)
            DEALLOCATE(this%current)
            this%current => p
         END DO
      END SUBROUTINE DestructList
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Add( this, d )
         IMPLICIT NONE 
         TYPE(LinkedList) :: this
         TYPE(ListData)   :: d
         TYPE(LLREcord), POINTER   :: newRecord
         
         ALLOCATE( newRecord )
         IF ( .NOT.ASSOCIATED(this%tail) )     THEN
            this%head => newRecord
            this%tail => newRecord
         ELSE
            this%tail%next => newRecord
            this%tail => newRecord
         END IF
         this%current => this%tail
         NULLIFY(this%current%next)! => NULL()
         CALL Copy( d, this%current%d )
      END SUBROUTINE Add
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE GetCurrentData( this, d )
      IMPLICIT NONE 
         TYPE(LinkedList) :: this
         TYPE(ListData)   :: d
         CALL copy( this%current%d, d )
      END SUBROUTINE GetCurrentData
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetToStart( this )
      IMPLICIT NONE 
         TYPE(LinkedList) :: this
         this%current => this%head
      END SUBROUTINE SetToStart
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE MoveToNext( this )
         IMPLICIT NONE 
         TYPE(LinkedList) :: this
         this%current => this%current%next
      END SUBROUTINE MoveToNext
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintList( this )
         IMPLICIT NONE 
         TYPE(LinkedList)        :: this
         this%current => this%head
         DO WHILE (ASSOCIATED(this%current) )
            CALL PRINT(this%current%d)
            CALL MoveToNext( this )
         END DO
      END SUBROUTINE PrintList
      
      END Module LinkedListClass
