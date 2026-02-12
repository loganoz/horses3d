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
        integer :: apeeq
    contains
    end type apeSource_t

  type(apeSource_t)                                                :: apeSource

  contains

!/////////////////////////////////////////////////////////////////////////
!           APE PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

    subroutine constructAPESource(self, controlVariables)
        use FTValueDictionaryClass
        implicit none
        type(apeSource_t) :: self
        type(FTValueDictionary), intent(in)                     :: controlVariables

        ! Read APE equation from control file
        self % apeeq = controlVariables % getValueOrDefault("APE equation", APEEQ_4)
        if (self % apeeq .ne. APEEQ_4) then
            print *, "Only APE-4 equation is implemented."
            error stop
        end if

    end subroutine constructAPESource
!
    Subroutine addAPESource(self, mesh)
        type(apeSource_t)                                    :: self
        type(HexMesh), intent(inout)                         :: mesh

        select case (self % apeeq)
        case (APEEQ_4)
            call addAPE4Source(self, mesh)
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

        qrho = 0.0

!$omp do schedule(runtime) private(i,j,k,eID)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                    qvel =  e % storage % Lambbase(:,i,j,k) - e % storage % Lamb_NS(:,i,j,k)
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
!
#endif
End Module APESourceClass
