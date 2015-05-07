!
!////////////////////////////////////////////////////////////////////////
!
!      ListDataClass.f
!      Created: 2008-06-24 09:43:48 -0400 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      Module ListDataClass 
      IMPLICIT NONE
      
      TYPE ListData
         INTEGER             :: id, i, j, key
      END TYPE ListData
      
      INTERFACE Construct
         MODULE PROCEDURE ConstructListData
         MODULE PROCEDURE ConstructCornerData
      END INTERFACE Construct
      INTERFACE Destruct
         MODULE PROCEDURE DestructListData
      END INTERFACE Destruct
      INTERFACE PRINT
         MODULE PROCEDURE PrintListData
      END INTERFACE PRINT
!
!     ========
      CONTAINS 
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructListData( this ) 
      IMPLICIT NONE 
         TYPE(ListData) :: this
         this%id = 0
         this%i = 0
         this%j = 0
         this%key = 0
      END SUBROUTINE ConstructListData
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructCornerData( this, id, i, j) 
      IMPLICIT NONE 
         TYPE(ListData) :: this
         INTEGER        :: id, i, j
         this%id = id
         this%i = i
         this%j = j
         this%key = 0
      END SUBROUTINE ConstructCornerData
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructListData( this ) 
         IMPLICIT NONE 
         TYPE(ListData) :: this
      END SUBROUTINE DestructListData
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Copy( this, that ) 
         IMPLICIT NONE 
         TYPE(ListData) :: this, that
         that = this
      END SUBROUTINE Copy
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintListData( this ) 
         IMPLICIT NONE 
         TYPE(ListData) :: this
         PRINT *, this%id, this%i, this%j, this%key
      END SUBROUTINE PrintListData
      
      END Module ListDataClass
      