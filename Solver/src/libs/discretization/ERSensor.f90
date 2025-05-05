#include "Includes.h"
Module ERSensor  !
    use DGSEMClass
    use HexMeshClass
    use SMConstants
    use FTValueDictionaryClass
    Implicit None

private

public initializeERFilters, filterSolutions, getER

    type(DGSem)                                  :: sem_filterA
    type(DGSem)                                  :: sem_filterB

    contains

    Subroutine initializeERFilters(sem,controlVariables)
        Implicit None
        type(DGSem)            , intent(in)                      :: sem                !<  Fine sem class
        type(FTValueDictionary), intent(in)                      :: controlVariables   !<  Input variables

        type(FTValueDictionary)                                  :: controlVariablesFilterA
        type(FTValueDictionary)                                  :: controlVariablesFilterB
        integer                                                  :: ntemp
        logical                                                  :: success
        integer, dimension(:), allocatable                       :: Nx_A, Ny_A, Nz_A
        integer, dimension(:), allocatable                       :: Nx_B, Ny_B, Nz_B

        print *, "initializing ER sensor"
        controlVariablesFilterA = controlVariables
        controlVariablesFilterB = controlVariables

        ntemp = controlVariablesFilterA % integerValueForKey("polorder filter a")
        print *, "polynomial filter A: ", ntemp
        call controlVariablesFilterA % removeObjectForKey("polynomial order")
        call controlVariablesFilterA % addValueForKey(ntemp,TRIM("polynomial order"))

        call getSimpleOrders(controlVariablesFilterA,Nx_A,Ny_A,Nz_A)
        call sem_filterA % construct (controlVariables = controlVariablesFilterA, Nx_ = Nx_A, Ny_ = Ny_A, Nz_ = Nz_A, success = success, ChildSem = .true. )

        ntemp = controlVariablesFilterB % integerValueForKey("polorder filter b")
        print *, "polynomial filter B: ", ntemp
        call controlVariablesFilterB % removeObjectForKey("polynomial order")
        call controlVariablesFilterB % addValueForKey(ntemp,TRIM("polynomial order"))

        call getSimpleOrders(controlVariablesFilterB,Nx_B,Ny_B,Nz_B)
        call sem_filterB % construct (controlVariables = controlVariablesFilterB, Nx_ = Nx_B, Ny_ = Ny_B, Nz_ = Nz_B, success = success, ChildSem = .true. )
        ! print *, "end initializing ER sensor"

    End Subroutine initializeERFilters

    Subroutine filterSolutions(mesh)
        Implicit None

        TYPE(HexMesh), intent(in)           :: mesh

        integer     :: eID
        integer     :: nodeType

        ! print *, "filtering solution"
        nodeType = mesh % nodeType

        do eID =1, mesh % no_of_elements
            call mesh%elements(eID)%storage%InterpolateSolution(sem_filterA%mesh%elements(eID)%storage, nodeType, with_gradients=.false.)
            call mesh%elements(eID)%storage%InterpolateSolution(sem_filterB%mesh%elements(eID)%storage, nodeType, with_gradients=.false.)
        end do
        ! print *, "end filtering solution"

    End Subroutine filterSolutions

    Subroutine getER(mesh)
        Implicit None
        TYPE(HexMesh), intent(inout)           :: mesh

        integer                              :: eID
        real(kind=RP), dimension(NDIM)       :: integralSolution, integralFilterA, integralFilterB
        real(kind=RP)                        :: ER

        ! print *, "getting ER sensor"

        do eID =1, mesh % no_of_elements

            call elementMomentumInteral(mesh % elements(eID), integralSolution)
            call elementMomentumInteral(sem_filterA % mesh % elements(eID), integralFilterA)
            call elementMomentumInteral(sem_filterB % mesh % elements(eID), integralFilterB)

            ER = sum(POW2(integralSolution - integralFilterA)) / (sum(POW2(integralSolution - integralFilterB)) + 1.0E-15_RP)

            mesh % elements(eID) % storage % sensor = ER

        end do
!
    End Subroutine getER

    Subroutine elementMomentumInteral(e,integral)
        use ElementClass
        use NodalStorageClass
        use PhysicsStorage
        Implicit None
        TYPE(Element), intent(in)           :: e
        real(kind=RP), intent(out)          :: integral(NDIM)

        integer                              :: i, j, k

        integral= 0.0_RP

        associate ( wx => NodalStorage(e % Nxyz(1)) % w, &
                    wy => NodalStorage(e % Nxyz(2)) % w, &
                    wz => NodalStorage(e % Nxyz(3)) % w    )
           do k = 0, e % Nxyz(3)  ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
              integral= integral + wx(i) * wy(j) * wz(k) * e % storage % Q(IRHOU:IRHOW,i,j,k) * e % geom % jacobian(i,j,k)
           end do            ; end do           ; end do
!
        end associate
    End Subroutine elementMomentumInteral

    Subroutine getSimpleOrders(controlVariables, Nx, Ny, Nz)
        use FileReaders                     , only: ReadOrderFile
        use ReadMeshFile                    , only: NumOfElemsFromMeshFile
        Implicit None
        integer, allocatable   , intent(inout) :: Nx(:), Ny(:), Nz(:)  
        type(FTValueDictionary), intent(in)                      :: controlVariables

        integer              :: nelem
!
!     *************************************
!     Read the simulation polynomial orders
!     *************************************
!
      if (controlVariables % containsKey("polynomial order file")) then
         call ReadOrderFile( controlVariables % stringValueForKey("polynomial order file", requestedLength = LINE_LENGTH), &
                             Nx, Ny, Nz )
      else
         nelem = NumOfElemsFromMeshFile( controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH) )
         allocate( Nx(nelem), Ny(nelem), Nz(nelem) )
         
         if (controlVariables % containsKey("polynomial order")) then
            Nx = controlVariables % integerValueForKey("polynomial order")
            Ny = Nx
            Nz = Nx
         else
            if (controlVariables % containsKey("polynomial order i") .AND. &
                controlVariables % containsKey("polynomial order j") .AND. &
                controlVariables % containsKey("polynomial order k") ) then
               Nx = controlVariables % integerValueForKey("polynomial order i")
               Ny = controlVariables % integerValueForKey("polynomial order j")
               Nz = controlVariables % integerValueForKey("polynomial order k")
            else
               error stop "The polynomial order(s) must be specified"
            end if
         end if
      end if

    End Subroutine getSimpleOrders

End Module ERSensor
