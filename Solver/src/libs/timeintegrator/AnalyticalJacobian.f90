!
!//////////////////////////////////////////////////////
!
!   @File:    AnalyticalJacobian.f90
!   @Author:  Andrés Rueda (a.rueda@upm.es)
!   @Created: Tue Oct 31 14:00:00 2017
!   @Last revision date: Tue Feb 13 19:37:39 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: c01958bbb74b2de9252027cd1c501fe081a58ef2
!
!//////////////////////////////////////////////////////
!
!  This module provides the routines for computing the analytical Jacobian matrix
!  -> Implemented for inviscid terms
!  -> TODO: MPI implementation is missing
!  -> The Jacobian of the BCs is temporarily computed numerically
!
!//////////////////////////////////////////////////////
#include "Includes.h"
module AnalyticalJacobian
   use SMConstants
   use ElementClass
   use HexMeshClass
   use NodalStorageClass
   use PhysicsStorage
   use Physics
   use Jacobian
   use MatrixClass
   use DGSEMClass
   use StopWatchClass
   use MeshTypes
   implicit none
   
   private
   public AnalyticalJacobian_Compute
   
!
!  ----------
!  Parameters
!  ----------
!
   real(kind=RP), parameter :: jaceps = 1.e-8_RP ! Minimum value of a Jacobian entry (smaller values are considered as 0._RP)
   
   ! Variables to be moved to Jacobian Storage
   integer, allocatable :: ndofelm(:)
   integer, allocatable :: firstIdx(:)
   
contains

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  Subroutine for computing the analytical Jacobian matrix
!  -------------------------------------------------------
   subroutine AnalyticalJacobian_Compute(sem,time,Matrix,BlockDiagonalized)
      implicit none
      !--------------------------------------------
      type(DGSem)              , intent(inout) :: sem
      real(kind=RP)            , intent(in)    :: time
      class(Matrix_t)          , intent(inout) :: Matrix
      logical        , optional, intent(in)    :: BlockDiagonalized  !<? Construct only the block diagonal?
      !--------------------------------------------
      integer :: eID, fID  ! Element and face counters
      integer :: nnz
      integer :: nelem
      logical :: BlockDiagonal
      logical, save :: IsFirst = .TRUE.
      !--------------------------------------------
      
      if (IsFirst) then
         call Stopwatch % CreateNewEvent("Analytical Jacobian construction")
         IsFirst = .FALSE.
      end if
      
      call Stopwatch % Start("Analytical Jacobian construction")
!
!     Initializations
!     ---------------
      
      if ( present(BlockDiagonalized) ) then
         BlockDiagonal = BlockDiagonalized
      else
         BlockDiagonal = .FALSE.
      end if
      
      nelem = size(sem % mesh % elements)
      
!
!     Get block sizes and position in matrix. This can be moved elsewhere since it's needed by both the numerical and analytical Jacobians
!     ---------------------------------------
      
      safedeallocate(ndofelm)  ; allocate (ndofelm(nelem))
      safedeallocate(firstIdx) ; allocate (firstIdx(nelem+1))
      
      firstIdx = 1
      DO eID=1, nelem
         ndofelm(eID)  = NCONS * (sem % Nx(eID)+1) * (sem % Ny(eID)+1) * (sem % Nz(eID)+1)
         if (eID>1) firstIdx(eID) = firstIdx(eID-1) + ndofelm(eID-1)
      END DO
      firstIdx(nelem+1) = firstIdx(nelem) + ndofelm(nelem)
         
      nnz = MAXVAL(ndofelm) ! 
      
!
!     Preallocate Jacobian matrix
!     ---------------------------
      
      call Matrix % Preallocate(nnz)
      call Matrix % Reset
!$omp parallel
!
!     ***************
!     Diagonal blocks
!     ***************
!
      call AnalyticalJacobian_DiagonalBlocks(sem % mesh, time, sem % externalState, Matrix)
!
!     *******************
!     Off-Diagonal blocks
!     *******************
!
      if (.not. BlockDiagonal) call AnalyticalJacobian_OffDiagonalBlocks(sem % mesh,Matrix)
!$omp end parallel
      
!
!     Finish assembling matrix
!     ------------------------
      
      call Matrix % Assembly(firstIdx,ndofelm)
      
      call Stopwatch % Pause("Analytical Jacobian construction")
      
      write(STD_OUT,'(A,ES10.3,A)') "Analytical Jacobian construction: ", Stopwatch % Elapsedtime("Analytical Jacobian construction"), ' seconds'
      
      call Stopwatch % Reset("Analytical Jacobian construction")
   end subroutine AnalyticalJacobian_Compute
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the off-diagonal blocks of the Jacobian matrix
!  -----------------------------------------------------------------------------------------------
   subroutine AnalyticalJacobian_DiagonalBlocks(mesh,time,externalStateProcedure,Matrix)
      use FaceClass
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout) :: mesh
      real(kind=RP)            , intent(in)    :: time
      procedure(BCState_FCN)                   :: externalStateProcedure
      class(Matrix_t)          , intent(inout) :: Matrix
      !--------------------------------------------
      integer :: eID, fID
      !--------------------------------------------
#if defined(NAVIERSTOKES)
      call ComputeNumericalFluxJacobian(mesh,time,externalStateProcedure)
#endif

!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID)) 
         call f % ProjectFluxJacobianToElements(LEFT ,LEFT )   ! dF/dQL to the left element 
         if (.not. (f % faceType == HMESH_BOUNDARY)) call f % ProjectFluxJacobianToElements(RIGHT,RIGHT)   ! dF/dQR to the right element
                                                               ! check if right element's rotation is ok for the Jacobian
         end associate
      end do
!$omp end do

!$omp do schedule(runtime)
      do eID = 1, size(mesh % elements)
         associate (e=> mesh % elements(eID))
         call computeBlock(e,mesh,Matrix)
         end associate
      end do
!$omp end do
      
   end subroutine AnalyticalJacobian_DiagonalBlocks   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the off-diagonal blocks of the Jacobian matrix
!  -----------------------------------------------------------------------------------------------
   subroutine AnalyticalJacobian_OffDiagonalBlocks(mesh,Matrix)
      use FaceClass
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout) :: mesh
      class(Matrix_t)          , intent(inout) :: Matrix
      !--------------------------------------------
      
      ERROR stop ':: Analytical Jacobian not implemented for off-diagonal blocks'
      
   end subroutine AnalyticalJacobian_OffDiagonalBlocks
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  
!  -----------------------------------------------------------------------------------------------
#if defined(NAVIERSTOKES)
   subroutine ComputeNumericalFluxJacobian(mesh,time,externalStateProcedure)
      use RiemannSolvers
      use FaceClass
      implicit none
      !--------------------------------------------
      type(HexMesh), intent(inout)    :: mesh
      real(kind=RP), intent(in)       :: time
      procedure(BCState_FCN)          :: externalStateProcedure
      !--------------------------------------------
      integer :: fID
      !--------------------------------------------
      
      call mesh % ProlongSolutionToFaces
      
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID) )
         select case (f % faceType)
            case (HMESH_INTERIOR)
               call ComputeInterfaceFluxJacobian(f)
            case (HMESH_BOUNDARY)
               call ComputeBoundaryFluxJacobian(f,time,externalStateProcedure)
         end select
         end associate
      end do
!$omp end do
   
   end subroutine ComputeNumericalFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeBoundaryFluxJacobian(f,time,externalStateProcedure)
      use RiemannSolvers
      use FaceClass
      implicit none
      !--------------------------------------------
      type(Face), intent(inout) :: f
      real(kind=RP), intent(in) :: time
      procedure(BCState_FCN)    :: externalStateProcedure
      !--------------------------------------------
      integer :: i,j
      real(kind=RP) :: BCjac(NCONS,NCONS)
      real(kind=RP), dimension(NCONS,NCONS) :: dFStar_dqL, dFStar_dqR
      !--------------------------------------------
      
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
!
!        Get external state
!        ------------------
         
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL externalStateProcedure( f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      f % boundaryType )
!
!        Get numerical flux jacobian on the face point (i,j)
!        ---------------------------------------------------

         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = dFStar_dqL, & ! this is dFStar/dqL
                                 side = LEFT)
         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = dFStar_dqR, & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         dFStar_dqL = dFStar_dqL * f % geom % jacobian(i,j)
         dFStar_dqR = dFStar_dqR * f % geom % jacobian(i,j)
!
!        Correct dFstar/dQL with the Jacobian of the boundary condition
!        --------------------------------------------------------------
         call ExternalStateJacobian( f % geom % x(:,i,j), &
                                     time, &
                                     f % geom % normal(:,i,j), &
                                     f % storage(1) % Q(:,i,j),&
                                     f % storage(2) % Q(:,i,j),&
                                     f % boundaryType, &
                                     externalStateProcedure, &
                                     BCjac )
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = dFStar_dqL !&
!~                                            + matmul(dFStar_dqR,BCjac)  ! TODO: Why is this not working anymore???
      end do             ; end do
      
   end subroutine ComputeBoundaryFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeInterfaceFluxJacobian(f)
      use RiemannSolvers
      use FaceClass
      implicit none
      !--------------------------------------------
      type(Face), intent(inout) :: f
      !--------------------------------------------
      integer :: i,j
      !--------------------------------------------
      
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
!
!        Get numerical flux jacobian on the face point (i,j)
!        ---------------------------------------------------

         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = f % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = f % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!           Scale with the mapping Jacobian
!           -------------------------------
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         f % storage(RIGHT) % dFStar_dqF (:,:,i,j) = f % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         
      end do             ; end do
      
   end subroutine ComputeInterfaceFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
   subroutine ExternalStateJacobian(x,time,nHat,Qin,Qex,boundaryType,externalStateProcedure,BCjac)
      implicit none
      !--------------------------------------------
      real(kind=RP), intent(in)       :: x(3)
      real(kind=RP), intent(in)       :: time
      real(kind=RP), intent(in)       :: nHat(3)
      real(kind=RP), intent(in)       :: Qin(NCONS)
      real(kind=RP), intent(in)       :: Qex(NCONS)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
      procedure(BCState_FCN)          :: externalStateProcedure
      real(kind=RP), intent(out)      :: BCjac(NCONS,NCONS)
      !--------------------------------------------
      real(kind=RP) :: newQext (NCONS)
      real(kind=RP) :: q(NCONS), buffer
      real(kind=RP),parameter :: eps = 1.e-8_RP
      integer :: i
      !--------------------------------------------
      
      q = Qin
      
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         newQext = q
         CALL externalStateProcedure( x, time, nHat, newQext, boundaryType )
         
         BCjac(:,i) = (newQext-Qex)/eps
         
         q(i) = buffer
      end do
      
   end subroutine ExternalStateJacobian
#endif
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeBlock(e,mesh,Matrix)
      implicit none
      type(Element)            , intent(in) :: e
      type(HexMesh)            , intent(in)    :: mesh
      class(Matrix_t), intent(inout) :: Matrix
      !-------------------------------------
      real(kind=RP) :: dFdQ(NCONS,NCONS,0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3),NDIM) 
      real(kind=RP) :: LocalMat(ndofelm(e % eID),ndofelm(e % eID)) 
      integer :: irow_0(ndofelm(e % eID))         ! Row indexes for the matrix
      integer :: irow(ndofelm(e % eID))           ! Formatted row indexes for the column assignment (small values with index -1)
      integer              :: j, NDOFEL
      
      NDOFEL = ndofelm(e % eID)
      irow_0 (1:NDOFEL) = (/ (firstIdx(e % eID) + j, j=0,NDOFEL-1) /)   
         
         
         call ComputeContravariantFluxJacobian( e, dFdQ) 
         call Local_GetDiagonalBlock(e, dFdQ, &
               mesh % faces(e % faceIDs(EFRONT )) % storage(e %faceSide(EFRONT )) % dFStar_dqEl(:,:,:,:,e %faceSide(EFRONT )), &
               mesh % faces(e % faceIDs(EBACK  )) % storage(e %faceSide(EBACK  )) % dFStar_dqEl(:,:,:,:,e %faceSide(EBACK  )), &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e %faceSide(EBOTTOM)) % dFStar_dqEl(:,:,:,:,e %faceSide(EBOTTOM)), &
               mesh % faces(e % faceIDs(ERIGHT )) % storage(e %faceSide(ERIGHT )) % dFStar_dqEl(:,:,:,:,e %faceSide(ERIGHT )), &
               mesh % faces(e % faceIDs(ETOP   )) % storage(e %faceSide(ETOP   )) % dFStar_dqEl(:,:,:,:,e %faceSide(ETOP   )), &
               mesh % faces(e % faceIDs(ELEFT  )) % storage(e %faceSide(ELEFT  )) % dFStar_dqEl(:,:,:,:,e %faceSide(ELEFT  )), &
               LocalMat  )
      
      ! Dump to mat
      
      do j=1, NDOFEL
         irow = irow_0
         where (LocalMat(:,j) < jaceps) irow = -1
         
         call Matrix % SetColumn (NDOFEL, irow_0, firstIdx(e % eID) + j-1, LocalMat(:,j) )
      end do
      
   end subroutine ComputeBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  
!  -----------------------------------------------------------------------------------------------
   subroutine Local_GetDiagonalBlock(e,dFdQ,dfdq_fr,dfdq_ba,dfdq_bo,dfdq_ri,dfdq_to,dfdq_le,LocalMat)
      implicit none
      !-------------------------------------------
      type(Element)                    , intent(in) :: e
      real(kind=RP), dimension(NCONS,NCONS,0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3),3), intent(in) :: dFdQ
      real(kind=RP), dimension(NCONS,NCONS,0:e%Nxyz(1),0:e%Nxyz(3)), intent(in) :: dfdq_fr, dfdq_ba
      real(kind=RP), dimension(NCONS,NCONS,0:e%Nxyz(1),0:e%Nxyz(2)), intent(in) :: dfdq_bo, dfdq_to
      real(kind=RP), dimension(NCONS,NCONS,0:e%Nxyz(2),0:e%Nxyz(3)), intent(in) :: dfdq_ri, dfdq_le
      real(kind=RP) :: LocalMat(ndofelm(e % eID),ndofelm(e % eID))
      !-------------------------------------------
      integer :: i, j             ! Matrix indices
      integer :: i1, j1, k1, eq1  ! variable counters
      integer :: i2, j2, k2, eq2  ! variable counters
      integer :: nXi, nEta, nZeta ! Number of nodes in every direction
      real(kind=RP) :: di, dj, dk       ! Kronecker deltas
      !-------------------------------------------
      
      nXi   = e % Nxyz(1) + 1
      nEta  = e % Nxyz(2) + 1
      nZeta = e % Nxyz(3) + 1
      
      LocalMat = 0._RP
      do k2 = 0, e % Nxyz(3) ; do j2 = 0, e % Nxyz(2) ; do i2 = 0, e % Nxyz(1) ; do eq2 = 1, NCONS 
         do k1 = 0, e % Nxyz(3) ; do j1 = 0, e % Nxyz(2) ; do i1 = 0, e % Nxyz(1) ; do eq1 = 1, NCONS 
            
!           Kronecker deltas
!           -----------------

            if (i1 == i2) then
               di = 1._RP
            else
               di = 0._RP
            end if
            if (j1 == j2) then
               dj = 1._RP
            else
               dj = 0._RP
            end if
            if (k1 == k2) then
               dk = 1._RP
            else
               dk = 0._RP
            end if
            


            i = eq1 + i1*NCONS + j1*NCONS*nXi + k1*NCONS*nEta*nXi ! row index (1-based)
            j = eq2 + i2*NCONS + j2*NCONS*nXi + k2*NCONS*nEta*nXi ! column index (1-based)
            
            LocalMat(i,j) = &
            
!           Volumetric contribution
!           -----------------------
                           ( dFdQ(eq1,eq2,i2,j2,k2,1) * e % spAXi   % hatD(i1,i2) * dj * dk &
                           + dFdQ(eq1,eq2,i2,j2,k2,2) * e % spAEta  % hatD(j1,j2) * di * dk &
                           + dFdQ(eq1,eq2,i2,j2,k2,3) * e % spAZeta % hatD(k1,k2) * di * dj &
!           Faces contribution
!           ------------------
                           -   dfdq_fr(eq1,eq2,i1,k1) * e % spAeta % b(j1,FRONT ) * e % spAeta % v(j2,FRONT ) * di * dk   & ! 1 Front
                           -   dfdq_ba(eq1,eq2,i1,k1) * e % spAeta % b(j1,BACK  ) * e % spAeta % v(j2,BACK  ) * di * dk   & ! 2 Back
                           -   dfdq_bo(eq1,eq2,i1,j1) * e % spAZeta% b(k1,BOTTOM) * e % spAZeta% v(k2,BOTTOM) * di * dj   & ! 3 Bottom
                           -   dfdq_to(eq1,eq2,i1,j1) * e % spAZeta% b(k1,TOP   ) * e % spAZeta% v(k2,TOP   ) * di * dj   & ! 5 Top
                           -   dfdq_ri(eq1,eq2,j1,k1) * e % spAXi  % b(i1,RIGHT ) * e % spAXi  % v(i2,RIGHT ) * dj * dk   & ! 4 Right
                           -   dfdq_le(eq1,eq2,j1,k1) * e % spAXi  % b(i1,LEFT  ) * e % spAXi  % v(i2,LEFT  ) * dj * dk ) & ! 6 Left
                                                                                       / e % geom % jacobian(i1,j1,k1) ! Scale with Jacobian from mass matrix
            
         end do                ; end do                ; end do                ; end do
      end do                ; end do                ; end do                ; end do
      
   end subroutine Local_GetDiagonalBlock
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------
!  Subroutine to get the Jacobian of the contravariant fluxes
!  -> dFdQ (:,:,i,j,k,dim)
!              |     |
!           jac|coord|flux in cartesian direction dim 
!  ----------------------------------------------------------
   subroutine ComputeContravariantFluxJacobian( e, dFdQ) 
      implicit none
      !--------------------------------------------
      type(Element)              :: e
      real(kind=RP), intent(out) :: dFdQ( NCONS, NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3), NDIM )
      !--------------------------------------------
      real(kind=RP), DIMENSION(NCONS,NCONS)  :: dfdq_,dgdq_,dhdq_
      integer                                :: i,j,k
      !--------------------------------------------
#if defined(NAVIERSTOKES)      
      do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

         call InviscidJacobian(e % storage % Q(:,i,j,k),dfdq_,dgdq_,dhdq_)
         
         
         dFdQ(:,:,i,j,k,IX) = e % geom % jGradXi  (1,i,j,k) * dfdq_ + &
                              e % geom % jGradXi  (2,i,j,k) * dgdq_ + &
                              e % geom % jGradXi  (3,i,j,k) * dhdq_ 

         dFdQ(:,:,i,j,k,IY) = e % geom % jGradEta (1,i,j,k) * dfdq_ + &
                              e % geom % jGradEta (2,i,j,k) * dgdq_ + &
                              e % geom % jGradEta (3,i,j,k) * dhdq_ 

         dFdQ(:,:,i,j,k,IZ) = e % geom % jGradZeta(1,i,j,k) * dfdq_ + &
                              e % geom % jGradZeta(2,i,j,k) * dgdq_ + &
                              e % geom % jGradZeta(3,i,j,k) * dhdq_ 
      end do                ; end do                ; end do
#endif
   end subroutine ComputeContravariantFluxJacobian
   

end module AnalyticalJacobian
