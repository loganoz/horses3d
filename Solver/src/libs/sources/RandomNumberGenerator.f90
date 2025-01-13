!//////////////////////////////////////////////////////
#include "Includes.h"
Module RandomNumberGenerator_OpenACC
  use SMConstants
  use HexMeshClass
  use Utilities
  use FTValueDictionaryClass
  use PhysicsStorage
#ifdef _OPENACC
    use cudafor
    use openacc
    use curand
    use curand_device
    use openacc_curand
#endif
  Implicit None

  private
  public rndGeneratorOpenacc, rndnum_openacc_Addnoise

  type(curandGenerator) :: generator

  !definition of the complete tripClass 
  type RandomNumberGeneratorClass

    contains
    procedure :: construct      => rndnum_openacc_Construct
    procedure :: destruct       => rndnum_openacc_Destruct

  end type RandomNumberGeneratorClass

  type(RandomNumberGeneratorClass)  :: rndGeneratorOpenacc

  contains
!/////////////////////////////////////////////////////////////////////////

  Subroutine rndnum_openacc_Construct(self, mesh, controlVariables)
    use FileReadingUtilities, only: getRealArrayFromString, getCharArrayFromString
    use Utilities,            only: toLower
    use MPI_Process_Info
    use Headers
implicit none

    class(RandomNumberGeneratorClass)       :: self
    type(HexMesh), intent(in)               :: mesh
    type(FTValueDictionary), intent(in)     :: controlVariables

  end Subroutine rndnum_openacc_Construct
!
  subroutine rndnum_openacc_Destruct(self)

    class(RandomNumberGeneratorClass), intent(inout)    :: self


  end Subroutine rndnum_openacc_Destruct

  subroutine rndnum_openacc_Addnoise(self, mesh)
    implicit none
    class(RandomNumberGeneratorClass)       :: self
    type(HexMesh), intent(inout)            :: mesh

    integer :: numStates, seed
    type(curandStateXORWOW) :: h
    real(kind=RP) :: r, state
    integer :: i,j,k,eID, rand_int 

    call random_number(r)         
    rand_int = 1 + int(r * 4500)  

    !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(rand_int) async(1)
    do eID = 1 , size(mesh % elements)
        !$acc loop vector collapse(3) private(h,state,seed)
        do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
            seed = rand_int + eID + (k*(mesh % elements(eID) % Nxyz(1)+1)*(mesh % elements(eID) % Nxyz(2)+1) + j*(mesh % elements(eID) % Nxyz(2)+1) + i)  ! Unique seed
            ! Set the seed for the curandGenerator
            call curand_init(seed,0,0,h)
            state = curand_uniform(h)
            call addnoise(state, mesh % elements(eID) % geom % x(1,i,j,k), mesh % elements(eID) % storage % S_NS(:,i,j,k))

        end do ;  end do ; end do
    end do

  end subroutine rndnum_openacc_Addnoise

  subroutine addnoise(state, x, source_term)
    !$acc routine seq
    use ElementClass
    implicit none
    real(kind=RP) , intent(in)    :: state
    real(kind=RP) , intent(in)    :: x
    real(kind=RP) , intent(inout) :: source_term(NCONS)

    real(kind=RP) ::  d, La, dummy1

    d = 1.0_RP
    La = 6.0_RP

    if (x .le. La) then
      dummy1 = (state - 0.5_RP)*pi*2.0_RP
      source_term(2:4) = d*sin(dummy1)
    end if

  end subroutine addnoise


End Module RandomNumberGenerator_OpenACC