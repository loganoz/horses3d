        !COMPILER-GENERATED INTERFACE MODULE: Tue May 17 10:04:22 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE B3HS_HASH_KEY_JENKINS__genmod
          INTERFACE 
            FUNCTION B3HS_HASH_KEY_JENKINS(KEY,RANGE) RESULT(CODE)
              CHARACTER(*), INTENT(IN) :: KEY
              INTEGER(KIND=4), INTENT(IN) :: RANGE
              INTEGER(KIND=4) :: CODE
            END FUNCTION B3HS_HASH_KEY_JENKINS
          END INTERFACE 
        END MODULE B3HS_HASH_KEY_JENKINS__genmod
