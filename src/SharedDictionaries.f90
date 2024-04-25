      Module SharedBCModule
         USE FTValueDictionaryClass
         IMPLICIT NONE

         private
         public zoneNameDictionary
         public constructSharedBCModule, destructSharedBCModule

         type(FTValueDictionary) :: zoneNameDictionary
         
         CONTAINS
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE constructSharedBCModule  
            IMPLICIT NONE
            CALL zoneNameDictionary % initWithSize(8)
         END SUBROUTINE constructSharedBCModule
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE destructSharedBCModule  
            IMPLICIT NONE
            call zoneNameDictionary % destruct()
         END SUBROUTINE destructSharedBCModule
         
      END Module SharedBCModule