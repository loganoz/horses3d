        !COMPILER-GENERATED INTERFACE MODULE: Tue May 24 12:21:51 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDPERIODICOPERATION__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDPERIODICOPERATION(MESH,TIME,DT,       &
     &MONITORS)
              USE RESIDUALSMONITORCLASS
              USE MONITORSCLASS
              USE MPI_IBMUTILITIES
              USE KDCLASS
              USE TESSELLATIONTYPES
              USE IBMCLASS
              USE MPI_FACE_CLASS
              USE TRANSFINITEMAPCLASS
              USE FACEPATCHCLASS
              USE CONNECTIVITYCLASS
              USE ELEMENTCLASS
              USE MAPPEDGEOMETRYCLASS
              USE FACECLASS
              USE NODECLASS
              USE STORAGECLASS
              USE HEXMESHCLASS
              CLASS (HEXMESH) :: MESH
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: DT
              TYPE (MONITOR_T), INTENT(IN) :: MONITORS
            END SUBROUTINE USERDEFINEDPERIODICOPERATION
          END INTERFACE 
        END MODULE USERDEFINEDPERIODICOPERATION__genmod
