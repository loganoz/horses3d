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

         private
         public bcValueDictionary, bcTypeDictionary, zoneNameDictionary, conformingBoundariesDic
         public constructSharedBCModule, destructSharedBCModule

         TYPE(FTValueDictionary) :: bcValueDictionary, bcTypeDictionary
         type(FTValueDictionary) :: zoneNameDictionary, conformingBoundariesDic
         
         CONTAINS
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE constructSharedBCModule  
            IMPLICIT NONE
            CALL bcValueDictionary % initWithSize(8)
            CALL bcTypeDictionary  % initWithSize(8)
            CALL zoneNameDictionary % initWithSize(8)
            call conformingBoundariesDic % initWithSize(8)
         END SUBROUTINE constructSharedBCModule
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE destructSharedBCModule  
            IMPLICIT NONE
            CALL bcValueDictionary % destruct()
            CALL bcTypeDictionary  % destruct()
            call zoneNameDictionary % destruct()
            call conformingBoundariesDic % destruct()
         END SUBROUTINE destructSharedBCModule
         
      END Module SharedBCModule
