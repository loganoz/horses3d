        !COMPILER-GENERATED INTERFACE MODULE: Tue May 24 12:21:51 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDFINALIZE__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDFINALIZE(MESH,TIME,ITER,MAXRESIDUAL,  &
     &MULTIPHASE_,MONITORS,ELAPSEDTIME,CPUTIME)
              USE RESIDUALSMONITORCLASS
              USE MONITORSCLASS
              USE FLUIDDATA
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
              INTEGER(KIND=4) :: ITER
              REAL(KIND=8) :: MAXRESIDUAL
              TYPE (MULTIPHASE_T), INTENT(IN) :: MULTIPHASE_
              TYPE (MONITOR_T), INTENT(IN) :: MONITORS
              REAL(KIND=8), INTENT(IN) :: ELAPSEDTIME
              REAL(KIND=8), INTENT(IN) :: CPUTIME
            END SUBROUTINE USERDEFINEDFINALIZE
          END INTERFACE 
        END MODULE USERDEFINEDFINALIZE__genmod
