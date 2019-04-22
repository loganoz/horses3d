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
Module Jacobian 

USE SMConstants
USE HexMeshClass
USE PhysicsStorage 
IMPLICIT NONE 

   private
   public Look_for_neighbour, ijk2local, local2ijk
   public JACEPS
   
   real(kind=RP), parameter :: JACEPS = 1.e-8_RP ! Minimum value of a Jacobian entry (smaller values are considered as 0._RP)
   
!========
 CONTAINS
!========
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


