        !COMPILER-GENERATED INTERFACE MODULE: Tue May 17 10:08:23 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DISPLAYSIMULATIONSTATISTICS__genmod
          INTERFACE 
            SUBROUTINE DISPLAYSIMULATIONSTATISTICS(ITER,MESH)
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
              INTEGER(KIND=4), INTENT(IN) :: ITER
              TYPE (HEXMESH), INTENT(IN) :: MESH
            END SUBROUTINE DISPLAYSIMULATIONSTATISTICS
          END INTERFACE 
        END MODULE DISPLAYSIMULATIONSTATISTICS__genmod
