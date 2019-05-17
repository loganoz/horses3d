!
!////////////////////////////////////////////////////////////////////////
!
!      Jacobian.f90
!      Created: 2011-09-27 16:21:15 +0200 
!      By: Gonzalo Rubio Calzado
!          Carlos Redondo
!          Andrés Rueda  
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module Jacobian 

   use SMConstants
   use HexMeshClass
   use PhysicsStorage 
   use MPI_Process_Info    , only: MPI_Process
   use Utilities           , only: QsortWithFriend
   use DGSEMClass          , only: DGSem, ComputeTimeDerivative_f
   use GenericMatrixClass  , only: Matrix_t
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none 

   private
   public Look_for_neighbour, ijk2local, local2ijk
   public JACEPS, Jacobian_t
   
   real(kind=RP), parameter :: JACEPS = 1.e-9_RP ! Minimum value of a Jacobian entry (smaller values are considered as 0._RP)
   
!  ------------------
!  Main Jacobian type
!  ------------------
   type :: Jacobian_t
      integer, allocatable :: ndofelm_l(:)   ! Array containing the number of degrees of freedom for all the LOCAL elements
      integer, allocatable :: ndofelm (:)    ! Array containing the number of degrees of freedom for all the elements
      integer, allocatable :: firstIdx(:)    ! Array containing the first index of the matrix for every element (supposing ordered numbering)
      integer              :: nEqn
      contains
         procedure :: Construct  => Jacobian_Construct
         procedure :: Destruct   => Jacobian_Destruct
         procedure :: Compute    => Jacobian_Compute
   end type Jacobian_t
   
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------
!     Construct the JacobianInfo variables
!     ------------------------------------
      subroutine Jacobian_Construct(this, mesh, nEqn)
         implicit none
         !-arguments-----------------------------------------
         class(Jacobian_t), intent(inout) :: this
         type(HexMesh)    , intent(in)    :: mesh
         integer          , intent(in)    :: nEqn
         !-local-variables-----------------------------------
         integer :: eID, ierr
         integer :: globIDs_l  (mesh % no_of_elements)
         integer :: all_globIDs(mesh % no_of_allElements)
         integer :: all_ndofelm(mesh % no_of_allElements)
         integer :: no_of_elems(MPI_Process % nProcs)
         integer :: displs     (MPI_Process % nProcs)
         !---------------------------------------------------
         
         this % nEqn = nEqn
         
!
!        ***************
!        Allocate memory
!        ***************
!
         allocate ( this % ndofelm_l(mesh % no_of_elements       ) )
         allocate ( this % ndofelm  (mesh % no_of_allElements    ) )
         allocate ( this % firstIdx (mesh % no_of_allElements + 1) )
!
!        *********************
!        Get local block sizes
!        *********************
!
         do eID=1, mesh % no_of_elements
            this % ndofelm_l(eID)  = nEqn * (mesh % Nx(eID)+1) * (mesh % Ny(eID)+1) * (mesh % Nz(eID)+1)
         end do
!
!        ***********
!        MPI sharing
!        ***********
!
         if (MPI_Process % doMPIAction) then
#ifdef _HAS_MPI_
!
!           Get global element IDs in all partitions
!           ----------------------------------------
            do eID=1, mesh % no_of_elements
               globIDs_l(eID) = mesh % elements(eID) % globID
            end do
!
!           Share info with other processes
!           -------------------------------

            call mpi_allgather(mesh % no_of_elements, 1, MPI_INT, no_of_elems, 1, MPI_INT, MPI_COMM_WORLD, ierr)
            
            displs(1) = 0
            do eID = 1, MPI_Process % nProcs-1
               displs(eID+1) = displs(eID) + no_of_elems(eID)
            end do
            
            call mpi_allgatherv(globIDs_l, mesh % no_of_elements, MPI_INT, all_globIDs, no_of_elems , displs, MPI_INT, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(this % ndofelm_l, mesh % no_of_elements, MPI_INT, all_ndofelm, no_of_elems , displs, MPI_INT, MPI_COMM_WORLD, ierr)
!
!           Reorganize the ndofelm
!           ----------------------
            call QsortWithFriend(all_globIDs, all_ndofelm)
            
            this % ndofelm = all_ndofelm
#endif
         else
            this % ndofelm = this % ndofelm_l
         end if
         
!
!        ************
!        Get firstIdx
!        ************
!
         this % firstIdx(1) = 1
         do eID=1, mesh % no_of_allElements
            this % firstIdx(eID+1) = this % firstIdx(eID) + this % ndofelm(eID)
         end do
         
      end subroutine Jacobian_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------
!     Destruct the JacobianInfo variables
!     ------------------------------------
      subroutine Jacobian_Destruct(this)
         implicit none
         !-arguments-----------------------------------------
         class(Jacobian_t), intent(inout) :: this
         !---------------------------------------------------
         
         deallocate (this % ndofelm)
         deallocate (this % firstIdx)
         
      end subroutine Jacobian_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------
!     Generic subroutine for computing the Jacobian matrix
!     ----------------------------------------------------
      subroutine Jacobian_Compute(this, sem, nEqn, time, Matrix, TimeDerivative, eps_in, BlockDiagonalized)
         implicit none
         !-arguments----------------------------------
         class(Jacobian_t)        , intent(inout)     :: this
         type(DGSem)              , intent(inout)     :: sem
         integer,                   intent(in)        :: nEqn
         real(kind=RP)            , intent(in)        :: time
         class(Matrix_t)          , intent(inout)     :: Matrix
         procedure(ComputeTimeDerivative_f), optional :: TimeDerivative      !   
         real(kind=RP)  , optional, intent(in)        :: eps_in
         logical        , optional, intent(in)        :: BlockDiagonalized  !<? Construct only the block diagonal? (Only for AnJacobian_t)
         !--------------------------------------------
         
         ERROR stop 'Jacobian_t must be cast to NumJacobian_t or AnJacobian_t'
         
      end subroutine Jacobian_Compute
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Look_for_neighbour(this, mesh)
         IMPLICIT NONE 
!
!     -----------------------------------------------------------------------------
!     This subroutine finds the neighbours of all elements (conforming meshes only)
!        current implementation only finds neighbour 3D elements, but not faces
!     -----------------------------------------------------------------------------
!
!
!     ----------------
!     Input parameters
!     ----------------
!
         TYPE(Neighbor_t)                  :: this(:)
         TYPE(HexMesh)                     :: mesh
!
!     ---------------
!     Local Variables
!     ---------------
!
         INTEGER                         :: i,j,iEl
         
         DO iEl = 1, SIZE(mesh%elements)
            this(iEl)%elmnt(7) = iEl  ! The last one is itself
            DO j = 1, 6
               IF (mesh%elements(iEl)%NumberOfConnections(j) == 0) THEN
                  this(iEl)%elmnt(j) = 0
               ELSE
                  this(iEl)%elmnt(j) = mesh%elements(iEl)%Connection(j)%ElementIDs(1)
               ENDIF
            ENDDO
         ENDDO
    
      END SUBROUTINE 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
!  Returns the local index relative to an element from the local coordinates: i(lagrange node x), j(lagrange node y), 
!  k(lagrange node z), l(equation number)
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
   
      FUNCTION ijk2local(i,j,k,l,N_EQN,Nx,Ny,Nz) RESULT(idx)
         IMPLICIT NONE
         
         INTEGER, INTENT(IN)   :: i, j, k, l, Nx, Ny, Nz, N_EQN
         INTEGER               :: idx
         
         IF (l < 1 .OR. l > N_EQN)  STOP 'error in ijk2local, l has a wrong value'
         IF (i < 0 .OR. i > Nx)     STOP 'error in ijk2local, i has a wrong value'
         IF (j < 0 .OR. j > Ny)     STOP 'error in ijk2local, j has a wrong value'
         IF (k < 0 .OR. k > Nz)     STOP 'error in ijk2local, k has a wrong value'
         
         idx = k*(Nx+1)*(Ny+1)*N_EQN + j*(Nx+1)*N_EQN + i*N_EQN + l
      END FUNCTION
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Returns the coordinates relative to an element: l(equation number), i(lagrange node x), j(lagrange node y), k(lagrange node z)
!  from the local index  
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
   
      FUNCTION local2ijk(idx,N_EQN,Nx,Ny,Nz) RESULT (indices)
      
         INTEGER, INTENT(IN)   :: idx, Nx, Ny, Nz, N_EQN
         INTEGER               :: indices(4)
         INTEGER               :: tmp1, tmp2
         
         IF (idx < 1 .OR. idx > (Nx+1)*(Ny+1)*(Nz+1)*N_EQN) STOP 'error in local2ijk, idx has wrong value'
         
         indices(4) = (idx-1) / ((Nx+1)*(Ny+1) * N_EQN)
         tmp1       = MOD((idx-1),((Nx+1)*(Ny+1) * N_EQN) )
         indices(3) = tmp1 / ((Nx+1)*N_EQN)
         tmp2       = MOD(tmp1,((Nx+1)*N_EQN))
         indices(2) = tmp2 / (N_EQN)
         indices(1) = MOD(tmp2, N_EQN) + 1
      END FUNCTION

END MODULE  


