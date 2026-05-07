!//////////////////////////////////////////////////////
!
!This class represents the source term for Acoustic Perturbation Equation

#include "Includes.h"
Module APESourceClass  !
#if defined(ACOUSTIC) 
    use SMConstants
    use HexMeshClass
    use SolutionFile
    use PhysicsStorage, only: NCONS
    use PhysicsStorage_CAA
    use StorageClass, only: ElementStorage_t
    Implicit None

    private
    public apeSource, addAPESource, constructAPESource

    ! Acoustic Perturbation Equation number
    integer, parameter :: APEEQ_4 = 4

    ! Definition of apeSource class 
    type apeSource_t
        logical :: isActive = .false.
        integer :: apeeq
        logical :: useLambVector = .false.
        logical :: useCube = .false.
        integer, allocatable :: cube_elements(:)
    contains
    end type apeSource_t

  type(apeSource_t)                                                :: apeSource

  contains

!/////////////////////////////////////////////////////////////////////////
!           APE PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

    subroutine constructAPESource(self, mesh, controlVariables)
        use Physics_CAAKeywordsModule
        use FTValueDictionaryClass
        use HexMeshClass
        use FileReadingUtilities, only: GetRealArrayFromString
        implicit none
        type(apeSource_t) :: self
        type(HexMesh)    :: mesh
        type(FTValueDictionary), intent(in)                     :: controlVariables

        ! Local variables
        real(rp), allocatable :: cubeCoords(:)

        self % isActive = controlVariables % logicalValueForKey(SOURCE_TERM_KEY)
        if (.not. self % isActive) return

        ! Read APE equation from control file
        self % apeeq = controlVariables % getValueOrDefault(APE_NUMBER_KEY, APEEQ_4)
        if (self % apeeq .ne. APEEQ_4) then
            print *, "Only APE-4 equation is implemented."
            error stop
        end if

        self % useLambVector = controlVariables % logicalValueForKey(LAMB_VECTOR_KEY)

        if (self % useLambVector) then
            if (controlVariables % containsKey(LambVectorCubeKey)) then
                self % useCube = .true.
                ! Read cube coordinates
                cubeCoords = GetRealArrayFromString( controlVariables % StringValueForKey(LambVectorCubeKey,requestedLength = LINE_LENGTH))
                if (size(cubeCoords) .ne. 6) then
                    print *, "ERROR: size(cubeCoords) != 6"
                    print *, "Cube coordinates read: ", cubeCoords
                    print *, "[x0, x1, y0, y1, z0, z1] is expected."
                    error stop
                end if
                ! Get elements inside cube
                call mesh % findElementsInsideCube(cubeCoords, self % cube_elements)
            end if
        end if

    end subroutine constructAPESource
!
    Subroutine addAPESource(self, mesh)
        type(apeSource_t)                                    :: self
        type(HexMesh), intent(inout)                         :: mesh

        ! Check if is activated
        if (.not. self % isActive) return

        select case (self % apeeq)
        case (APEEQ_4)
            if (self % useLambVector) then
                if (self % useCube) then
                    call addAPE4Source_withoutLamb(self, mesh)
                    call addAPE4Source_onlyLamb(self, mesh)
                else
                    call addAPE4Source(self, mesh)
                end if
            else
                call addAPE4Source_withoutLamb(self, mesh)
            end if
        end select
    end subroutine addAPESource

    !/////////////////////////////////////////////////////////////////////////
    ! APE-4
    !/////////////////////////////////////////////////////////////////////////
    Subroutine addAPE4Source(self, mesh)
        type(apeSource_t)                                    :: self
        type(HexMesh), intent(inout)                         :: mesh
    
        ! Local variables
        integer                                                 :: i, j, k, eID
        real(rp)                                                :: qrho, qvel(NDIM), qp, aux(NDIM)

!$omp do schedule(runtime) private(i,j,k,eID,qrho,qvel,qp,aux)
            do eID = 1, mesh % no_of_elements
                qrho = 0.0
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                    qvel =  e % storage % Lambbase(:,i,j,k) - e % storage % Lamb(:,i,j,k)
                    aux = e % storage % Qbase(IBRHO,i,j,k) * e % storage % Q(ICAAU:ICAAW,i,j,k) + &
                          e % storage % Qbase(IBU:IBW,i,j,k) * e % storage % Q(ICAAP,i,j,k) / e % storage % Qbase(IBA2,i,j,k)
                    qp = dot_product(aux, e % storage % grada2base(:,i,j,k))
                    e % storage % S_NS(ICAARHO,i,j,k) = e % storage % S_NS(ICAARHO,i,j,k) + qrho
                    e % storage % S_NS(ICAAU:ICAAW,i,j,k) = e % storage % S_NS(ICAAU:ICAAW,i,j,k) + qvel
                    e % storage % S_NS(ICAAP,i,j,k) = e % storage % S_NS(ICAAP,i,j,k) + qp
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

    End Subroutine addAPE4Source

    Subroutine addAPE4Source_withoutLamb(self, mesh)
        type(apeSource_t)                                    :: self
        type(HexMesh), intent(inout)                         :: mesh
    
        ! Local variables
        integer                                                 :: i, j, k, eID
        real(rp)                                                :: qrho, qp, aux(NDIM)

!$omp do schedule(runtime) private(i,j,k,eID,qrho,qp,aux)
            do eID = 1, mesh % no_of_elements
                qrho = 0.0
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                    aux = e % storage % Qbase(IBRHO,i,j,k) * e % storage % Q(ICAAU:ICAAW,i,j,k) + &
                          e % storage % Qbase(IBU:IBW,i,j,k) * e % storage % Q(ICAAP,i,j,k) / e % storage % Qbase(IBA2,i,j,k)
                    qp = dot_product(aux, e % storage % grada2base(:,i,j,k))
                    e % storage % S_NS(ICAARHO,i,j,k) = e % storage % S_NS(ICAARHO,i,j,k) + qrho
                    e % storage % S_NS(ICAAP,i,j,k) = e % storage % S_NS(ICAAP,i,j,k) + qp
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

    End Subroutine addAPE4Source_withoutLamb
!
    Subroutine addAPE4Source_onlyLamb(self, mesh)
        type(apeSource_t)                                    :: self
        type(HexMesh), intent(inout)                         :: mesh
    
        ! Local variables
        integer                                                 :: i, j, k, eID
        real(rp)                                                :: qvel(NDIM)

!$omp do schedule(runtime) private(i,j,k,eID,qvel)
            do eID = 1, size(self % cube_elements)
               associate ( e => mesh % elements(self % cube_elements(eID)) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                    qvel =  e % storage % Lambbase(:,i,j,k) - e % storage % Lamb(:,i,j,k)
                    e % storage % S_NS(ICAAU:ICAAW,i,j,k) = e % storage % S_NS(ICAAU:ICAAW,i,j,k) + qvel
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

    End Subroutine addAPE4Source_onlyLamb
#endif
End Module APESourceClass
