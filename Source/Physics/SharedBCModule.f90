!
!////////////////////////////////////////////////////////////////////////
!
!      SharedBCModule.f90
!      Created: 2011-07-05 10:15:35 -0400 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      Module SharedBCModule
         USE FTValueDictionaryClass
         IMPLICIT NONE
         TYPE(FTValueDictionary) :: bcValueDictionary, bcTypeDictionary
         
         CONTAINS
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE constructSharedBCModule  
            IMPLICIT NONE
            CALL bcValueDictionary % initWithSize(8)
            CALL bcTypeDictionary  % initWithSize(8)
         END SUBROUTINE constructSharedBCModule
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE destructSharedBCModule  
            IMPLICIT NONE
            CALL bcValueDictionary % release()
            CALL bcTypeDictionary  % release()
         END SUBROUTINE destructSharedBCModule
         
      END Module SharedBCModule
