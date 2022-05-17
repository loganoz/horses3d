        !COMPILER-GENERATED INTERFACE MODULE: Tue May 17 10:06:13 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDINITIALCONDITION__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDINITIALCONDITION(MESH,THERMODYNAMICS_,&
     &DIMENSIONLESS_,REFVALUES_)
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
              TYPE (THERMODYNAMICS_T), INTENT(IN) :: THERMODYNAMICS_
              TYPE (DIMENSIONLESS_T), INTENT(IN) :: DIMENSIONLESS_
              TYPE (REFVALUES_T), INTENT(IN) :: REFVALUES_
            END SUBROUTINE USERDEFINEDINITIALCONDITION
          END INTERFACE 
        END MODULE USERDEFINEDINITIALCONDITION__genmod
