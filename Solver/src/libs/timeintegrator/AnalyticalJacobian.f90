!
!//////////////////////////////////////////////////////
!
!   @File:    AnalyticalJacobian.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Oct 31 14:00:00 2017
!   @Last revision date: Thu Feb 14 16:16:40 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 17d60e4e57235a57aa406023ebe4c26157bc211a
!
!//////////////////////////////////////////////////////
!
!  This module provides the routines for computing the analytical Jacobian matrix
!  -> Only for p-conforming representations (TODO: make general)
!  -> TODO: Implement as class to prevent memory leaking and additional computations
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
   use Jacobian                        , only: JACEPS
   use MatrixClass
   use DGSEMClass
   use StopWatchClass
   use MeshTypes
   use EllipticDiscretizations
   use ElementConnectivityDefinitions  , only: axisMap, normalAxis
   use BoundaryConditions              , only: BCs
   use FaceClass                       , only: Face
   implicit none
   
   private
   public AnalyticalJacobian_Compute
   
   integer, parameter :: other(2) = [2, 1]
   
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
   subroutine AnalyticalJacobian_Compute(sem,nEqn, time,Matrix,BlockDiagonalized)
      implicit none
      !--------------------------------------------
      type(DGSem)              , intent(inout) :: sem
      integer,                   intent(in)    :: nEqn
      real(kind=RP)            , intent(in)    :: time
      class(Matrix_t)          , intent(inout) :: Matrix
      logical        , optional, intent(in)    :: BlockDiagonalized  !<? Construct only the block diagonal?
#if defined(NAVIERSTOKES)
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
         
         ! TODO: Add conformity check!
         
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
         ndofelm(eID)  = nEqn * (sem % mesh % Nx(eID)+1) * (sem % mesh % Ny(eID)+1) * (sem % mesh % Nz(eID)+1)
         if (eID>1) firstIdx(eID) = firstIdx(eID-1) + ndofelm(eID-1)
      END DO
      firstIdx(nelem+1) = firstIdx(nelem) + ndofelm(nelem)
      
!
!     ***************************
!     Preallocate Jacobian matrix
!     ***************************
!
      select type(Matrix_p => Matrix)
!
!        If block-diagonal matrix, construct with size of blocks
!        -------------------------------------------------------
         type is(DenseBlockDiagMatrix_t)
            call Matrix_p % Preallocate(nnzs=ndofelm)
         type is(SparseBlockDiagMatrix_t)
            call Matrix_p % Preallocate(nnzs=ndofelm)
            
!
!        If matrix is CSR, standard preallocate with LinkedListMatrix
!        ------------------------------------------------------------
         type is(csrMat_t)
            call Matrix_p % Preallocate()
!
!        Otherwise, construct with nonzeros in each row
!        ----------------------------------------------
         class default 
            if (BlockDiagonal) then
               nnz = maxval(ndofelm)
            else
               nnz = 7*maxval(ndofelm) ! 20180201: hard-coded to 7, since the analytical Jacobian can only compute off-diagonal blocks for compact schemes (neighbors' effect)
            end if
            call Matrix_p % Preallocate(nnz)
      end select
      
      call Matrix % Reset (ForceDiagonal = .TRUE.)
      call Matrix % SpecifyBlockInfo(firstIdx,ndofelm)
      
!$omp parallel
!
!     **************************************************
!     Compute the Jacobian of the Numerical Flux (FStar)
!     **************************************************
!
      call ComputeNumericalFluxJacobian(sem % mesh,nEqn,time)
!
!     ***************
!     Diagonal blocks
!     ***************
!
      call AnalyticalJacobian_DiagonalBlocks(sem % mesh, nEqn,time, Matrix)
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
      
      call Matrix % Assembly()
      
      call Stopwatch % Pause("Analytical Jacobian construction")
      
      write(STD_OUT,'(A,ES10.3,A)') "Analytical Jacobian construction: ", Stopwatch % lastElapsedtime("Analytical Jacobian construction"), ' seconds'
      
#else
      ERROR stop ':: Analytical Jacobian only for NS'
#endif
   end subroutine AnalyticalJacobian_Compute
#if defined(NAVIERSTOKES)
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the diagonal blocks of the Jacobian matrix
!  -------------------------------------------------------------------------------------------
   subroutine AnalyticalJacobian_DiagonalBlocks(mesh,nEqn,time,Matrix)
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout) :: mesh
      integer,                   intent(in)    :: nEqn
      real(kind=RP)            , intent(in)    :: time
      class(Matrix_t)          , intent(inout) :: Matrix
      !--------------------------------------------
      integer :: eID, fID
      !--------------------------------------------
!
!     Project flux Jacobian to corresponding element
!     ----------------------------------------------
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID)) 
         call f % ProjectFluxJacobianToElements(nEqn,LEFT ,LEFT )   ! dF/dQL to the left element 
         if (.not. (f % faceType == HMESH_BOUNDARY)) call f % ProjectFluxJacobianToElements(nEqn,RIGHT,RIGHT)   ! dF/dQR to the right element
         end associate
      end do
!$omp end do
!
!     Project flux Jacobian with respect to gradients to corresponding elements
!     -------------------------------------------------------------------------
      if (flowIsNavierStokes) then
!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            associate (f => mesh % faces(fID)) 
            call f % ProjectGradJacobianToElements(LEFT, LEFT)   ! dF/dQL to the left element 
            if (.not. (f % faceType == HMESH_BOUNDARY)) call f % ProjectGradJacobianToElements(RIGHT,RIGHT)   ! dF/dQR to the right element
            end associate
         end do
!$omp end do
      end if
!
!     Compute each element's diagonal block
!     -------------------------------------
!$omp do schedule(runtime)
      do eID = 1, size(mesh % elements)
         associate (e=> mesh % elements(eID))
         call Local_SetDiagonalBlock( e, &
                                      mesh % faces( e % faceIDs(EFRONT ) ), &
                                      mesh % faces( e % faceIDs(EBACK  ) ), &
                                      mesh % faces( e % faceIDs(EBOTTOM) ), &
                                      mesh % faces( e % faceIDs(ERIGHT ) ), &
                                      mesh % faces( e % faceIDs(ETOP   ) ), &
                                      mesh % faces( e % faceIDs(ELEFT  ) ), &
                                      Matrix )
         end associate
      end do
!$omp end do
      
   end subroutine AnalyticalJacobian_DiagonalBlocks   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the off-diagonal blocks of the Jacobian matrix
!  -> Note: Only the interior faces contribute to the off-diagonal blocks
!  -> only for p-conforming meshes
!  -----------------------------------------------------------------------------------------------
   subroutine AnalyticalJacobian_OffDiagonalBlocks(mesh,Matrix)
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout) :: mesh
      class(Matrix_t)          , intent(inout) :: Matrix
      !--------------------------------------------
      integer :: eID, fID, elSide, side
      type(Element), pointer :: e_plus
      type(Element), pointer :: e_minus
      type(Face)   , pointer :: f
      !--------------------------------------------
      
!
!     Project flux Jacobian to opposed elements (RIGHT to LEFT and LEFT to RIGHT)
!     ---------------------------------------------------------------------------
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID)) 
         if (f % faceType == HMESH_INTERIOR) then
            call f % ProjectFluxJacobianToElements(NCONS, LEFT ,RIGHT)   ! dF/dQR to the left element
            call f % ProjectFluxJacobianToElements(NCONS, RIGHT,LEFT )   ! dF/dQL to the right element 
         end if
         end associate
      end do
!$omp end do

!
!     Project flux Jacobian with respect to gradients to opposed elements (RIGHT to LEFT and LEFT to RIGHT)
!     -----------------------------------------------------------------------------------------------------
      if (flowIsNavierStokes) then
!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            associate (f => mesh % faces(fID)) 
            if (f % faceType == HMESH_INTERIOR) then
               call f % ProjectGradJacobianToElements(LEFT ,RIGHT)   ! dF/dGradQR to the left element
               call f % ProjectGradJacobianToElements(RIGHT,LEFT )   ! dF/dGradQL to the right element 
            end if
            end associate
         end do
!$omp end do
      end if
      
!
!     Compute the off-diagonal blocks for each element's equations
!     ------------------------------------------------------------
!$omp do schedule(runtime) private(e_plus,elSide,fID,side,f,e_minus)
      do eID = 1, mesh % no_of_elements
         e_plus => mesh % elements(eID)
!
!        One block for every neighbor element
!        ------------------------------------
         do elSide = 1, 6
            if (e_plus % NumberOfConnections(elSide) == 0) cycle
            
            fID  = e_plus % faceIDs(elSide)
            side = e_plus % faceSide(elSide)
            
            f => mesh % faces(fID)
            e_minus => mesh % elements(f % elementIDs( other(side) ))
            
            call Local_GetOffDiagonalBlock(f,e_plus,e_minus,side,Matrix)
         end do
      end do
!$omp end do
      nullify (f, e_plus, e_minus)
   end subroutine AnalyticalJacobian_OffDiagonalBlocks
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeNumericalFluxJacobian(mesh,nEqn,time)
      use RiemannSolvers_NS
      implicit none
      !--------------------------------------------
      type(HexMesh), intent(inout)    :: mesh
      integer,       intent(in)       :: nEqn
      real(kind=RP), intent(in)       :: time
      !--------------------------------------------
      integer :: fID
      !--------------------------------------------
      
      call mesh % ProlongSolutionToFaces(NCONS)
      if (flowIsNavierStokes) call mesh % ProlongGradientsToFaces(NGRAD)
      
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID) )
         select case (f % faceType)
            case (HMESH_INTERIOR)
               call ComputeInterfaceFluxJacobian(f)
            case (HMESH_BOUNDARY)
               call ComputeBoundaryFluxJacobian(f,time)
         end select
         end associate
      end do
!$omp end do
   
   end subroutine ComputeNumericalFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Computes the Jacobian of the numerical flux on a boundary as
!     dF*    dF*(Q⁺,Q⁻,n^)     dF*(Q⁺,Q⁻,n^)     dQ⁻
!    ---- = --------------- + --------------- * ----,
!     dQ⁺      dQ⁺               dQ⁻             dQ⁺
!  where the last term is the Jacobian of the boundary condition.... +: internal state
!                                                                    -: external state
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeBoundaryFluxJacobian(f,time)
      use RiemannSolvers_NS
      implicit none
      !--------------------------------------------
      type(Face), intent(inout) :: f
      real(kind=RP), intent(in) :: time
      !--------------------------------------------
      integer :: i,j
      real(kind=RP) :: BCjac(NCONS,NCONS)
      real(kind=RP), dimension(NCONS,NCONS) :: dFStar_dqL, dFStar_dqR
      !--------------------------------------------
      
!
!     *********************
!     Inviscid contribution
!     *********************
!
      
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
!
!        Get external state
!        ------------------
         
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % StateForEqn( NCONS, &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))
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
!        Scale with the mapping Jacobian
!        -------------------------------
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         f % storage(RIGHT) % dFStar_dqF (:,:,i,j) = f % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)

      end do             ; end do
!
!     ********************
!     Viscous contribution
!     ********************
!
      if (flowIsNavierStokes) call ViscousDiscretization % RiemannSolver_Jacobians(f)  ! TODO: Check if external gradient has to be taken into account
!
!     **************************************************************
!     Correct dFstar/dQL with the Jacobian of the boundary condition
!     **************************************************************
!
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
         
         call ExternalStateJacobian( f % geom % x(:,i,j), &
                                     time, &
                                     f % geom % normal(:,i,j), &
                                     f % storage(1) % Q(:,i,j),&
                                     f % storage(2) % Q(:,i,j),&
                                     f % zone, &
                                     BCjac )
         
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT  ) % dFStar_dqF (:,:,i,j) &
                                            + matmul(f % storage(RIGHT ) % dFStar_dqF (:,:,i,j),BCjac)
      end do             ; end do
      
   end subroutine ComputeBoundaryFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Computes the Jacobians of the numerical flux for an inner face as
!                         dF*(Q⁺,Q⁻,n^)                        dF*(Q⁺,Q⁻,n^) 
!           dFStar/dqL = ---------------         dFStar/dqR = ---------------
!                          dQ⁺                                  dQ⁻         
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeInterfaceFluxJacobian(f)
      use RiemannSolvers_NS
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
!        Scale with the mapping Jacobian
!        -------------------------------
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         f % storage(RIGHT) % dFStar_dqF (:,:,i,j) = f % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         
      end do             ; end do
      
      if (flowIsNavierStokes) call ViscousDiscretization % RiemannSolver_Jacobians(f)
      
   end subroutine ComputeInterfaceFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  This routine obtains the Jacobian of the boundary condition numerically.
!  This can be optimized introducing the analytical Jacobian of every single implemented BC...
!  -----------------------------------------------------------------------------------------------
   subroutine ExternalStateJacobian(x,time,nHat,Qin,Qex,zone,BCjac)
      implicit none
      !--------------------------------------------
      real(kind=RP), intent(in)       :: x(3)
      real(kind=RP), intent(in)       :: time
      real(kind=RP), intent(in)       :: nHat(3)
      real(kind=RP), intent(in)       :: Qin(NCONS)
      real(kind=RP), intent(in)       :: Qex(NCONS)
      integer,       intent(in)       :: zone
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
         CALL BCs(zone) % bc % StateForEqn( NCONS, x, time, nHat, newQext)
         
         BCjac(:,i) = (newQext-Qex)/eps
         
         q(i) = buffer
      end do
      
   end subroutine ExternalStateJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------
!  Local_SetDiagonalBlock:
!     Subroutine to set a diagonal block in the Jacobian
!  -----------------------------------------------------
   subroutine Local_SetDiagonalBlock(e, fF, fB, fO, fR, fT, fL, Matrix)
      use HyperbolicDiscretizations
      implicit none
      !-------------------------------------------
      type(Element)  , intent(inout) :: e
      type(Face)     , intent(in)    :: fF, fB, fO, fR, fT, fL !< The six faces of the element
      class(Matrix_t), intent(inout) :: Matrix
      !-------------------------------------------
      real(kind=RP) :: MatEntry
      real(kind=RP) :: dFdQ      (NCONS,NCONS,     0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3),NDIM)
      real(kind=RP) :: dF_dgradQ (NCONS,NCONS,NDIM,NDIM,0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3))
      integer :: i, j             ! Matrix indices
      integer :: r                ! Additional index for counting
      integer :: i1, j1, k1, eq1  ! variable counters
      integer :: i2, j2, k2, eq2  ! variable counters
      integer :: nXi, nEta        ! Number of nodes in every direction
      integer :: EtaSpa, ZetaSpa  ! Spacing for these two coordinate directions
      real(kind=RP) :: di, dj, dk ! Kronecker deltas
      integer :: Deltas           ! A variable to know if enough deltas are zero, in which case this is a zero entry of the matrix..
      real(kind=RP), dimension(:,:,:,:)    , pointer :: dfdq_fr, dfdq_ba, dfdq_bo, dfdq_to, dfdq_le, dfdq_ri
      real(kind=RP), dimension(:,:,:,:,:,:), pointer :: dfdGradQ_fr, dfdGradQ_ba, dfdGradQ_bo, dfdGradQ_to, dfdGradQ_ri, dfdGradQ_le
      real(kind=RP) :: nF( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(3) )   ! Normal vecs for FRONT  face
      real(kind=RP) :: nB( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(3) )   ! Normal vecs for BACK   face
      real(kind=RP) :: nO( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2) )   ! Normal vecs for BOTTOM face
      real(kind=RP) :: nR( NDIM, 0:e % Nxyz(2), 0:e % Nxyz(3) )   ! Normal vecs for RIGHT  face
      real(kind=RP) :: nT( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2) )   ! Normal vecs for TOP    face
      real(kind=RP) :: nL( NDIM, 0:e % Nxyz(2), 0:e % Nxyz(3) )   ! Normal vecs for LEFT   face
      real(kind=RP) :: xiAux  (4), Gvec_xi  (3)
      real(kind=RP) :: etaAux (4), Gvec_eta (3)
      real(kind=RP) :: zetaAux(4), Gvec_zeta(3)
      real(kind=RP) :: a_plus, a_minus, temp
      !-------------------------------------------
!
!     *******************
!     Initial definitions
!     *******************
!
      nXi   = e % Nxyz(1) + 1
      nEta  = e % Nxyz(2) + 1
      EtaSpa  = NCONS*nXi
      ZetaSpa = NCONS*nXi*nEta
      
      a_plus  = 0.5_RP  ! Temp... TODO: read from ViscousDiscretization
      a_minus = 0.5_RP  ! Temp... TODO: read from ViscousDiscretization
      
!     Get normal vectors in element indexing
!     --------------------------------------      
      call AnJac_GetNormalVec(fF, e % faceSide(EFRONT ), nF)
      call AnJac_GetNormalVec(fB, e % faceSide(EBACK  ), nB)
      call AnJac_GetNormalVec(fO, e % faceSide(EBOTTOM), nO)
      call AnJac_GetNormalVec(fR, e % faceSide(ERIGHT ), nR)
      call AnJac_GetNormalVec(fT, e % faceSide(ETOP   ), nT)
      call AnJac_GetNormalVec(fL, e % faceSide(ELEFT  ), nL)
      
!
!     *********************
!     Inviscid contribution
!     *********************
!
      call HyperbolicDiscretization % ComputeInnerFluxJacobian( e, dFdQ) 
      if (flowIsNavierStokes) call ViscousDiscretization % ComputeInnerFluxJacobian( e, dF_dgradQ, dFdQ)
      
!
!     Pointers to the flux Jacobians with respect to q on the faces
!     -------------------------------------------------------------
      dfdq_fr(1:,1:,0:,0:) => e % Storage % dfdq_fr(1:,1:,0:,0:,e %faceSide(EFRONT ))
      dfdq_ba(1:,1:,0:,0:) => e % Storage % dfdq_ba(1:,1:,0:,0:,e %faceSide(EBACK  ))
      dfdq_bo(1:,1:,0:,0:) => e % Storage % dfdq_bo(1:,1:,0:,0:,e %faceSide(EBOTTOM))
      dfdq_to(1:,1:,0:,0:) => e % Storage % dfdq_to(1:,1:,0:,0:,e %faceSide(ETOP   ))
      dfdq_ri(1:,1:,0:,0:) => e % Storage % dfdq_ri(1:,1:,0:,0:,e %faceSide(ERIGHT ))
      dfdq_le(1:,1:,0:,0:) => e % Storage % dfdq_le(1:,1:,0:,0:,e %faceSide(ELEFT  ))
      
      do k2 = 0, e % Nxyz(3) ; do j2 = 0, e % Nxyz(2) ; do i2 = 0, e % Nxyz(1) ; do eq2 = 1, NCONS 
         do k1 = 0, e % Nxyz(3) ; do j1 = 0, e % Nxyz(2) ; do i1 = 0, e % Nxyz(1) ; do eq1 = 1, NCONS 
            Deltas = 0
!           Kronecker deltas
!           -----------------

            if (i1 == i2) then
               di = 1._RP
               Deltas = Deltas + 1
            else
               di = 0._RP
            end if
            if (j1 == j2) then
               dj = 1._RP
               Deltas = Deltas + 1
            else
               dj = 0._RP
            end if
            if (k1 == k2) then
               dk = 1._RP
               Deltas = Deltas + 1
            else
               dk = 0._RP
            end if
            
            if (Deltas < 2) cycle

            i = eq1 + i1*NCONS + j1*EtaSpa + k1*ZetaSpa ! row index (1-based)
            j = eq2 + i2*NCONS + j2*EtaSpa + k2*ZetaSpa ! column index (1-based)
            
            
            MatEntry = &
            
!           Volumetric contribution (inner fluxes)
!           **************************************
                           (  dFdQ(eq1,eq2,i2,j2,k2,1) * e % spAXi   % hatD(i1,i2) * dj * dk &
                            + dFdQ(eq1,eq2,i2,j2,k2,2) * e % spAEta  % hatD(j1,j2) * di * dk &
                            + dFdQ(eq1,eq2,i2,j2,k2,3) * e % spAZeta % hatD(k1,k2) * di * dj &
!           Faces contribution (numerical fluxes)
!           *************************************
                            -   dfdq_fr(eq1,eq2,i1,k1) * e % spAeta % b(j1,FRONT ) * e % spAeta % v(j2,FRONT ) * di * dk   & ! 1 Front
                            -   dfdq_ba(eq1,eq2,i1,k1) * e % spAeta % b(j1,BACK  ) * e % spAeta % v(j2,BACK  ) * di * dk   & ! 2 Back
                            -   dfdq_bo(eq1,eq2,i1,j1) * e % spAZeta% b(k1,BOTTOM) * e % spAZeta% v(k2,BOTTOM) * di * dj   & ! 3 Bottom
                            -   dfdq_to(eq1,eq2,i1,j1) * e % spAZeta% b(k1,TOP   ) * e % spAZeta% v(k2,TOP   ) * di * dj   & ! 5 Top
                            -   dfdq_ri(eq1,eq2,j1,k1) * e % spAXi  % b(i1,RIGHT ) * e % spAXi  % v(i2,RIGHT ) * dj * dk   & ! 4 Right
                            -   dfdq_le(eq1,eq2,j1,k1) * e % spAXi  % b(i1,LEFT  ) * e % spAXi  % v(i2,LEFT  ) * dj * dk ) & ! 6 Left
                                                                                       * e % geom % invJacobian(i1,j1,k1) ! Scale with Jacobian from mass matrix
            
            call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntry)
            
         end do                ; end do                ; end do                ; end do
      end do                ; end do                ; end do                ; end do
      nullify (dfdq_fr, dfdq_ba, dfdq_bo, dfdq_to, dfdq_le, dfdq_ri)
      
!
!     ********************
!     Viscous contribution
!     ********************
!      
      if (flowIsNavierStokes) then
!
!        Pointers to the flux Jacobians with respect to grad(q) on the faces
!        -------------------------------------------
         dfdGradQ_fr(1:,1:,1:,1:,0:,0:) => e % Storage % dfdGradQ_fr(1:,1:,1:,1:,0:,0:,e %faceSide(EFRONT  ))
         dfdGradQ_ba(1:,1:,1:,1:,0:,0:) => e % Storage % dfdGradQ_ba(1:,1:,1:,1:,0:,0:,e %faceSide(EBACK   ))
         dfdGradQ_bo(1:,1:,1:,1:,0:,0:) => e % Storage % dfdGradQ_bo(1:,1:,1:,1:,0:,0:,e %faceSide(EBOTTOM ))
         dfdGradQ_to(1:,1:,1:,1:,0:,0:) => e % Storage % dfdGradQ_to(1:,1:,1:,1:,0:,0:,e %faceSide(ETOP    ))
         dfdGradQ_ri(1:,1:,1:,1:,0:,0:) => e % Storage % dfdGradQ_ri(1:,1:,1:,1:,0:,0:,e %faceSide(ERIGHT  ))
         dfdGradQ_le(1:,1:,1:,1:,0:,0:) => e % Storage % dfdGradQ_le(1:,1:,1:,1:,0:,0:,e %faceSide(ELEFT   ))
         
         do k2 = 0, e % Nxyz(3) ; do j2 = 0, e % Nxyz(2) ; do i2 = 0, e % Nxyz(1) ; do eq2 = 1, NCONS 
            do k1 = 0, e % Nxyz(3) ; do j1 = 0, e % Nxyz(2) ; do i1 = 0, e % Nxyz(1) ; do eq1 = 1, NCONS 
               Deltas = 0
!              Kronecker deltas
!              -----------------

               if (i1 == i2) then
                  di = 1._RP
                  Deltas = Deltas + 1
               else
                  di = 0._RP
               end if
               if (j1 == j2) then
                  dj = 1._RP
                  Deltas = Deltas + 1
               else
                  dj = 0._RP
               end if
               if (k1 == k2) then
                  dk = 1._RP
                  Deltas = Deltas + 1
               else
                  dk = 0._RP
               end if
               
               if (Deltas < 1) cycle
               
               i = eq1 + i1*NCONS + j1*EtaSpa + k1*ZetaSpa ! row index (1-based)
               j = eq2 + i2*NCONS + j2*EtaSpa + k2*ZetaSpa ! column index (1-based)
               
!              Get Xi auxiliar terms
!              ---------------------
               xiAux   = 0._RP
               if (dj*dk > 0.5_RP) then
                  do r=0, e % Nxyz(1)
                     
                     temp = e % spAxi   % hatD(i1,r) * e % geom % invJacobian(r,j1,k1)
                     
                     Gvec_xi   = [dF_dgradQ(eq1,eq2,1,1,r,j1,k1), dF_dgradQ(eq1,eq2,2,1,r,j1,k1), dF_dgradQ(eq1,eq2,3,1,r,j1,k1)]
                     
                     xiAux(4)   = xiAux(4)   + temp * ( e % spAXi   % b(r,LEFT  ) * e % spAXi   % v(i2,LEFT  ) * dot_product( Gvec_xi  , nL(:,j2,k2) ) &
                                                       +e % spAXi   % b(r,RIGHT ) * e % spAXi   % v(i2,RIGHT ) * dot_product( Gvec_xi  , nR(:,j2,k2) ) )
                     temp = temp * e % spAxi   % hatD(r,i2)
                     
                     xiAux  (1:3) = xiAux  (1:3) + temp * Gvec_xi
                  end do
                  xiAux(4) = xiAux(4) * a_plus
               end if
               
!              Get Eta auxiliar terms
!              ----------------------
               etaAux   = 0._RP
               if (di*dk > 0.5_RP) then
                  do r=0, e % Nxyz(2)
                     
                     temp = e % spAEta  % hatD(j1,r) * e % geom % invJacobian(i1,r,k1)
                     
                     Gvec_eta  = [dF_dgradQ(eq1,eq2,1,2,i1,r,k1), dF_dgradQ(eq1,eq2,2,2,i1,r,k1), dF_dgradQ(eq1,eq2,3,2,i1,r,k1)]
                     
                     etaAux(4)  = etaAux(4)  + temp * ( e % spAEta  % b(r,FRONT ) * e % spAEta  % v(j2,FRONT ) * dot_product( Gvec_eta , nF(:,i2,k2) ) &
                                                       +e % spAEta  % b(r,BACK  ) * e % spAEta  % v(j2,BACK  ) * dot_product( Gvec_eta , nB(:,i2,k2) ) )
                     temp = temp * e % spAEta  % hatD(r,j2)
                     
                     etaAux (1:3) = etaAux (1:3) + temp * Gvec_eta
                  end do
                  etaAux(4) = etaAux(4) * a_plus
               end if
               
!              Get Zeta auxiliar terms
!              -----------------------
               zetaAux   = 0._RP
               if (di*dj > 0.5_RP) then
                  do r=0, e % Nxyz(3)
                     
                     temp = e % spAZeta % hatD(k1,r) * e % geom % invJacobian(i1,j1,r)
                     
                     Gvec_zeta = [dF_dgradQ(eq1,eq2,1,3,i1,j1,r), dF_dgradQ(eq1,eq2,2,3,i1,j1,r), dF_dgradQ(eq1,eq2,3,3,i1,j1,r)]
                     
                     zetaAux(4) = zetaAux(4) + temp * ( e % spAZeta % b(r,BOTTOM) * e % spAZeta % v(k2,BOTTOM) * dot_product( Gvec_zeta, nO(:,i2,j2) ) &
                                                       +e % spAZeta % b(r,TOP   ) * e % spAZeta % v(k2,TOP   ) * dot_product( Gvec_zeta, nT(:,i2,j2) ) )
                     temp = temp * e % spAZeta % hatD(r,k2)
                     
                     zetaAux(1:3) = zetaAux(1:3) + temp * Gvec_zeta
                  end do
                  zetaAux(4) = zetaAux(4) * a_plus
               end if
               
               Gvec_xi   = [dF_dgradQ(eq1,eq2,1,1,i2,j1,k1), dF_dgradQ(eq1,eq2,2,1,i2,j1,k1), dF_dgradQ(eq1,eq2,3,1,i2,j1,k1)]
               Gvec_eta  = [dF_dgradQ(eq1,eq2,1,2,i1,j2,k1), dF_dgradQ(eq1,eq2,2,2,i1,j2,k1), dF_dgradQ(eq1,eq2,3,2,i1,j2,k1)]
               Gvec_zeta = [dF_dgradQ(eq1,eq2,1,3,i1,j1,k2), dF_dgradQ(eq1,eq2,2,3,i1,j1,k2), dF_dgradQ(eq1,eq2,3,3,i1,j1,k2)]
               
               MatEntry = 0._RP
!
!              Volumetric contribution (inner fluxes)
!              **************************************
!                 Xi-component of the flux
!                 ------------------------
!                   Xi:
               MatEntry = MatEntry + dj * dk * ( dot_product(xiAux  (1:3), e % geom % jGradXi  (:,i2,j2,k2) ) - xiAux(4) )
               
               MatEntry = MatEntry + e % spAXi   % hatD(i1,i2) * e % geom % invJacobian(i2,j1,k1) *                                       &
!                   Eta:
                  (+ dk *  (  e % spAEta  % hatD(j1,j2) * dot_product( Gvec_xi  , e % geom % jGradEta (:,i2,j2,k2) )                      &
                            - a_plus * (  e % spAEta  % b(j1,FRONT ) * e % spAEta  % v(j2,FRONT ) * dot_product( Gvec_xi  , nF(:,i2,k2) ) &
                                        + e % spAEta  % b(j1,BACK  ) * e % spAEta  % v(j2,BACK  ) * dot_product( Gvec_xi  , nB(:,i2,k2) ) &
                                       )                                                                                                  &
                           )                                                                                                              &
!                   Zeta:
                   + dj *  (  e % spAZeta % hatD(k1,k2) * dot_product( Gvec_xi  , e % geom % jGradZeta(:,i2,j2,k2) )                      &
                            - a_plus * (  e % spAZeta % b(k1,BOTTOM) * e % spAZeta % v(k2,BOTTOM) * dot_product( Gvec_xi  , nO(:,i2,j2) ) &
                                        + e % spAZeta % b(k1,TOP   ) * e % spAZeta % v(k2,TOP   ) * dot_product( Gvec_xi  , nT(:,i2,j2) ) &
                                       )                                                                                                  &
                           )                                                                                                              &
                  )
!                 Eta-component of the flux
!                 -------------------------
!                   Eta:
               MatEntry = MatEntry + di * dk * ( dot_product(etaAux (1:3), e % geom % jGradEta (:,i2,j2,k2) ) - etaAux(4) )
               
               MatEntry = MatEntry + e % spAEta  % hatD(j1,j2) * e % geom % invJacobian(i1,j2,k1) *                                       &
!                   Xi:
                  (+ dk *  (  e % spAXi   % hatD(i1,i2) * dot_product( Gvec_eta , e % geom % jGradXi  (:,i2,j2,k2) )                      &
                            - a_plus * (  e % spAXi   % b(i1,LEFT  ) * e % spAXi   % v(i2,LEFT  ) * dot_product( Gvec_eta , nL(:,j2,k2) ) &
                                        + e % spAXi   % b(i1,RIGHT ) * e % spAXi   % v(i2,RIGHT ) * dot_product( Gvec_eta , nR(:,j2,k2) ) &
                                       )                                                                                                  &
                           )                                                                                                              &
!                   Zeta:
                   + di *  (  e % spAZeta % hatD(k1,k2) * dot_product( Gvec_eta , e % geom % jGradZeta(:,i2,j2,k2) )                      &
                            - a_plus * (  e % spAZeta % b(k1,BOTTOM) * e % spAZeta % v(k2,BOTTOM) * dot_product( Gvec_eta , nO(:,i2,j2) ) &
                                        + e % spAZeta % b(k1,TOP   ) * e % spAZeta % v(k2,TOP   ) * dot_product( Gvec_eta , nT(:,i2,j2) ) &
                                       )                                                                                                  &
                           )                                                                                                              &
                  )
!                 Zeta-component of the flux
!                 -------------------------
!                   Zeta:
               MatEntry = MatEntry + di * dj * ( dot_product(zetaAux(1:3), e % geom % jGradZeta(:,i2,j2,k2) ) - zetaAux(4))
               
               MatEntry = MatEntry + e % spAZeta % hatD(k1,k2) * e % geom % invJacobian(i1,j1,k2) *                                       &
!                   Xi:
                  (+ dj *  (  e % spAXi   % hatD(i1,i2) * dot_product( Gvec_zeta, e % geom % jGradXi  (:,i2,j2,k2) )                      &
                            - a_plus * (  e % spAXi   % b(i1,LEFT  ) * e % spAXi   % v(i2,LEFT  ) * dot_product( Gvec_zeta, nL(:,j2,k2) ) &
                                        + e % spAXi   % b(i1,RIGHT ) * e % spAXi   % v(i2,RIGHT ) * dot_product( Gvec_zeta, nR(:,j2,k2) ) &
                                       )                                                                                                  &
                           )                                                                                                              &
!                   Eta:
                   + di *  (  e % spAEta  % hatD(j1,j2) * dot_product( Gvec_zeta, e % geom % jGradEta (:,i2,j2,k2) )                      &
                            - a_plus * (  e % spAEta  % b(j1,FRONT ) * e % spAEta  % v(j2,FRONT ) * dot_product( Gvec_zeta, nF(:,i2,k2) ) &
                                        + e % spAEta  % b(j1,BACK  ) * e % spAEta  % v(j2,BACK  ) * dot_product( Gvec_zeta, nB(:,i2,k2) ) &
                                       )                                                                                                  &
                           )                                                                                                              &
                  )
               
               
               MatEntry = MatEntry + ( &
!
!              Faces contribution (surface integrals from the outer equation) - PENALTY TERM IS BEING CONSIDERED IN THE INVISCID PART - TODO: Reorganize storage to put it explicitely in another place (needed for purely viscous equations)
!              ***********************************************|**************
!                 Front face                 _________________|
!                 ----------                 |
                   +   dfdGradQ_fr(eq1,eq2,1,2,i1,k1) * e % spAEta  % b(j1,FRONT ) * e % spAXi   % D(i1,i2) * e % spAEta  % v(j2,FRONT ) * dk   & ! Xi
                   +   dfdGradQ_fr(eq1,eq2,2,2,i1,k1) * e % spAEta  % b(j1,FRONT ) * e % spAEta  % vd(j2,FRONT ) * di * dk                      & ! Eta
                   +   dfdGradQ_fr(eq1,eq2,3,2,i1,k1) * e % spAEta  % b(j1,FRONT ) * e % spAZeta % D(k1,k2) * e % spAEta  % v(j2,FRONT ) * di   & ! Zeta
!                 Back face
!                 ---------
                   +   dfdGradQ_ba(eq1,eq2,1,2,i1,k1) * e % spAEta  % b(j1,BACK  ) * e % spAXi   % D(i1,i2) * e % spAEta  % v(j2,BACK  ) * dk   & ! Xi
                   +   dfdGradQ_ba(eq1,eq2,2,2,i1,k1) * e % spAEta  % b(j1,BACK  ) * e % spAEta  % vd(j2,BACK  ) * di * dk                      & ! Eta
                   +   dfdGradQ_ba(eq1,eq2,3,2,i1,k1) * e % spAEta  % b(j1,BACK  ) * e % spAZeta % D(k1,k2) * e % spAEta  % v(j2,BACK  ) * di   & ! Zeta
!                 Bottom face
!                 -----------
                   +   dfdGradQ_bo(eq1,eq2,1,2,i1,j1) * e % spAZeta % b(k1,BOTTOM) * e % spAXi   % D(i1,i2) * e % spAZeta % v(k2,BOTTOM) * dj   & ! Xi
                   +   dfdGradQ_bo(eq1,eq2,2,2,i1,j1) * e % spAZeta % b(k1,BOTTOM) * e % spAEta  % D(j1,j2) * e % spAZeta % v(k2,BOTTOM) * di   & ! Eta
                   +   dfdGradQ_bo(eq1,eq2,3,2,i1,j1) * e % spAZeta % b(k1,BOTTOM) * e % spAZeta % vd(k2,BOTTOM) * di * dj                      & ! Zeta 
!                 Top face
!                 --------
                   +   dfdGradQ_to(eq1,eq2,1,2,i1,j1) * e % spAZeta % b(k1,TOP   ) * e % spAXi   % D(i1,i2) * e % spAZeta % v(k2,TOP   ) * dj   & ! Xi
                   +   dfdGradQ_to(eq1,eq2,2,2,i1,j1) * e % spAZeta % b(k1,TOP   ) * e % spAeta  % D(j1,j2) * e % spAZeta % v(k2,TOP   ) * di   & ! Eta
                   +   dfdGradQ_to(eq1,eq2,3,2,i1,j1) * e % spAZeta % b(k1,TOP   ) * e % spAZeta % vd(k2,TOP   ) * di * dj                      & ! Zeta
!                 Right face
!                 ----------
                   +   dfdGradQ_ri(eq1,eq2,1,2,j1,k1) * e % spAXi   % b(i1,RIGHT ) * e % spAXi   % vd(i2,RIGHT ) * dj * dk                      & ! Xi
                   +   dfdGradQ_ri(eq1,eq2,2,2,j1,k1) * e % spAXi   % b(i1,RIGHT ) * e % spAeta  % D(j1,j2) * e % spAXi   % v(i2,RIGHT ) * dk   & ! Eta
                   +   dfdGradQ_ri(eq1,eq2,3,2,j1,k1) * e % spAXi   % b(i1,RIGHT ) * e % spAZeta % D(k1,k2) * e % spAXi   % v(i2,RIGHT ) * dj   & ! Zeta
!                 Left face
!                 ---------
                   +   dfdGradQ_le(eq1,eq2,1,2,j1,k1) * e % spAXi   % b(i1,LEFT  ) * e % spAXi   % vd(i2,LEFT  ) * dj * dk                      & ! Xi
                   +   dfdGradQ_le(eq1,eq2,2,2,j1,k1) * e % spAXi   % b(i1,LEFT  ) * e % spAeta  % D(j1,j2) * e % spAXi   % v(i2,LEFT  ) * dk   & ! Eta
                   +   dfdGradQ_le(eq1,eq2,3,2,j1,k1) * e % spAXi   % b(i1,LEFT  ) * e % spAZeta % D(k1,k2) * e % spAXi   % v(i2,LEFT  ) * dj   & ! Zeta
                  )
               
               MatEntry = MatEntry * e % geom % invJacobian(i1,j1,k1) ! Scale with Jacobian from mass matrix
               
               call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntry)
               
            end do                ; end do                ; end do                ; end do
         end do                ; end do                ; end do                ; end do
         nullify (dfdGradQ_fr, dfdGradQ_ba, dfdGradQ_bo, dfdGradQ_to, dfdGradQ_ri, dfdGradQ_le )
      end if
      
   end subroutine Local_SetDiagonalBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Local_GetOffDiagonalBlock:
!     Routine to compute the off-diagonal block that results of connecting e_plus with e_minus for 
!     the equations that correspond to e_plus and inserting it in the Matrix
!
!     -> Currently, this routine only works for p-conforming representations
!           TODO: Add routine for p-nonconforming representations (block computation is more expensive and delta cycling is ruled out)
!  -----------------------------------------------------------------------------------------------
   subroutine Local_GetOffDiagonalBlock(f,e_plus,e_minus,side,Matrix)
      implicit none
      !-arguments----------------------------------------------------------------------
      type(Face), target, intent(in)    :: f       !<  Face connecting elements
      type(Element)     , intent(in)    :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
      type(Element)     , intent(in)    :: e_minus !<  Element that connects with e_plus through "f"
      integer           , intent(in)    :: side    !<  side of face where element (1) is
      class(Matrix_t)   , intent(inout) :: Matrix  !<> Jacobian matrix
      !-local-variables----------------------------------------------------------------
      integer :: i, j                     ! Matrix indexes
      integer :: i1, j1, k1, eq1          ! variable counters
      integer :: i2, j2, k2, eq2          ! variable counters
      integer :: nXi1, nEta1              ! Number of nodes in every direction
      integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
      integer :: nXi2, nEta2              ! Number of nodes in every direction
      integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
      integer :: elSide_plus              ! Element side where f is on e⁻
      integer :: elSide_minus             ! Element side where f is on e⁺
      integer :: elInd_plus(3)            ! Element indexes on e⁺
      integer :: elInd_minus(3)           ! Element indexes on e⁻
      integer :: normAx_plus              ! Normal axis to f on e⁺
      integer :: normAx_minus             ! Normal axis to f on e⁻
      integer :: tanAx_plus (2)           ! Tangent axes to f on e⁺
      integer :: tanAx_minus(2)           ! Tangent axes to f on e⁻
      integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
      integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
      integer :: faceInd_plus (2)         ! Face indexes on e⁺
      integer :: faceInd_minus(2)         ! Face indexes on e⁻
      integer :: NxyFace_plus (2)         ! Polynomial orders of face on element e⁺
      integer :: NxyFace_minus(2)         ! Polynomial orders of face on element e⁻ (only needed for viscous fluxes)
      integer :: faceInd_plus2minus(2)    ! Face indexes on e⁺ passed to the reference frame of e⁻
      integer :: faceInd_minus2plus(2)    ! Face indexes on e⁻ passed to the reference frame of e⁺ (only needed for viscous fluxes)
      integer :: dtan1, dtan2             ! Kronecker deltas in the tangent directions
      integer :: dtan1_minus, dtan2_minus ! Kronecker deltas in the tangent directions in the reference frame of e⁻ (only needed for viscous fluxes)
      integer :: Deltas                   ! Number of Kronecker deltas /= 0
      real(kind=RP) :: MatrixEntry        ! Value of the matrix entry
      type(NodalStorage_t), target  :: spA_plus (3)      ! Nodal storage in the different directions for e_plus  - local copy
      type(NodalStorage_t), target  :: spA_minus(3)      ! Nodal storage in the different directions for e_minus - local copy
      type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
      type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAnorm_minus     ! Nodal storage in the direction that is normal to the face for e⁻
      type(NodalStorage_t), pointer :: spAtan1_minus     ! Nodal storage in the tangent direction "1" to the face for e⁻ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAtan2_minus     ! Nodal storage in the tangent direction "2" to the face for e⁻ (only needed for viscous fluxes)
      real(kind=RP)       , pointer :: dfdq(:,:,:,:)     ! 
      real(kind=RP)       , pointer :: dfdGradQ(:,:,:,:,:,:)
      !--------------------------------------------------------------------------------
      
!
!     ***********
!     Definitions
!     ***********
!
      ! Entry spacing for element e⁺
      nXi1     = e_plus % Nxyz(1) + 1
      nEta1    = e_plus % Nxyz(2) + 1
      EtaSpa1  = NCONS*nXi1
      ZetaSpa1 = NCONS*nXi1*nEta1
      
      ! Entry spacing for element e⁻
      nXi2     = e_minus % Nxyz(1) + 1
      nEta2    = e_minus % Nxyz(2) + 1
      EtaSpa2  = NCONS*nXi2
      ZetaSpa2 = NCONS*nXi2*nEta2
      
      ! Element sides
      elSide_plus  = f % elementSide(side)
      elSide_minus = f % elementSide(other(side))
      
      ! Normal and tangent axes
      normAx_plus  = normalAxis (elSide_plus)
      normAx_minus = normalAxis (elSide_minus)
      tanAx_plus   = axisMap(:,f % elementSide(    side     ) )
      tanAx_minus  = axisMap(:,f % elementSide( other(side) ) )
      
      ! Side of axis where f is
      if (normAx_plus < 0) then
         normAxSide_plus = LEFT
      else
         normAxSide_plus = RIGHT
      end if
      normAx_plus = abs(normAx_plus)
      if (normAx_minus < 0) then
         normAxSide_minus = LEFT
      else
         normAxSide_minus = RIGHT
      end if
      normAx_minus = abs(normAx_minus)
      
      ! Nodal storage
      spA_plus  = NodalStorage(e_plus  % Nxyz)
      spA_minus = NodalStorage(e_minus % Nxyz)
      spAnorm_plus  => spA_plus( normAx_plus   )
      spAtan1_plus  => spA_plus( tanAx_plus(1) )
      spAtan2_plus  => spA_plus( tanAx_plus(2) )
      
      spAnorm_minus => spA_minus( normAx_minus   )
      spAtan1_minus => spA_minus( tanAx_minus(1) )
      spAtan2_minus => spA_minus( tanAx_minus(2) )
      
      ! Polynomial orders
      NxyFace_plus  = e_plus  % Nxyz ( tanAx_plus  )
      NxyFace_minus = e_minus % Nxyz ( tanAx_minus )
      
!
!     *********************
!     Inviscid contribution
!     *********************
!
      
!
!     Pointers to the flux Jacobians with respect to q on the faces
!     -------------------------------------------------------------
      dfdq(1:,1:,0:,0:) => f % storage(side) % dFStar_dqEl(1:,1:,0:,0:,other(side))
      
      do k2 = 0, e_minus % Nxyz(3) ; do j2 = 0, e_minus % Nxyz(2) ; do i2 = 0, e_minus % Nxyz(1) ; do eq2 = 1, NCONS 
         do k1 = 0, e_plus % Nxyz(3) ; do j1 = 0, e_plus % Nxyz(2) ; do i1 = 0, e_plus % Nxyz(1) ; do eq1 = 1, NCONS 
            
            elInd_plus  = [i1, j1, k1]
            elInd_minus = [i2, j2, k2]
            
            faceInd_plus  = elInd_plus ( tanAx_plus  )
            faceInd_minus = elInd_minus( tanAx_minus )
            
            call indexesOnOtherFace(faceInd_plus(1),faceInd_plus(2), NxyFace_plus(1),NxyFace_plus(2), f % rotation, side, faceInd_plus2minus(1),faceInd_plus2minus(2))
            
            ! "delta" cycling
            if ( any(faceInd_plus2minus /= faceInd_minus) ) cycle

            i = eq1 + i1*NCONS + j1*EtaSpa1 + k1*ZetaSpa1 ! row index (1-based)
            j = eq2 + i2*NCONS + j2*EtaSpa2 + k2*ZetaSpa2 ! column index (1-based)
            
            MatrixEntry = -   dfdq(eq1,eq2,faceInd_plus(1),faceInd_plus(2)) &
                            * spAnorm_plus  % b(elInd_plus ( normAx_plus  ), normAxSide_plus ) &
                            * spAnorm_minus % v(elInd_minus( normAx_minus ), normAxSide_minus )  & 
                            * e_plus % geom % invJacobian(i1,j1,k1)
            
            call Matrix % AddToBlockEntry (e_plus % GlobID, e_minus % GlobID, i, j, MatrixEntry)
         end do                ; end do                ; end do                ; end do
      end do                ; end do                ; end do                ; end do
      nullify(dfdq)
!
!     ********************
!     Viscous contribution
!     ********************
!
      if (flowIsNavierStokes) then
!
!        Pointers to the flux Jacobians with respect to grad(q) on the faces
!        -------------------------------------------
         dfdGradQ(1:,1:,1:,1:,0:,0:) => f % Storage(side) % dFv_dGradQEl (1:,1:,1:,1:,0:,0:, other(side) )
         
         do k2 = 0, e_minus % Nxyz(3) ; do j2 = 0, e_minus % Nxyz(2) ; do i2 = 0, e_minus % Nxyz(1) ; do eq2 = 1, NCONS 
            do k1 = 0, e_plus % Nxyz(3) ; do j1 = 0, e_plus % Nxyz(2) ; do i1 = 0, e_plus % Nxyz(1) ; do eq1 = 1, NCONS 
               
               elInd_plus  = [i1, j1, k1]
               elInd_minus = [i2, j2, k2]
               
               faceInd_plus  = elInd_plus ( tanAx_plus  )
               faceInd_minus = elInd_minus( tanAx_minus )
               
               call indexesOnOtherFace(faceInd_plus (1),faceInd_plus (2), NxyFace_plus (1),NxyFace_plus (2), f % rotation, side       , faceInd_plus2minus(1),faceInd_plus2minus(2))
               call indexesOnOtherFace(faceInd_minus(1),faceInd_minus(2), NxyFace_minus(1),NxyFace_minus(2), f % rotation, other(side), faceInd_minus2plus(1),faceInd_minus2plus(2))
               
!              Delta computation on "plus" side (used for delta cycling)
!              ---------------------------------------------------------
               Deltas = 0
               if ( faceInd_plus(1) == faceInd_minus2plus(1) ) then
                  dtan1 = 1._RP
                  Deltas = Deltas + 1
               else
                  dtan1 = 0._RP
               end if
               if ( faceInd_plus(2) == faceInd_minus2plus(2) ) then
                  dtan2 = 1._RP
                  Deltas = Deltas + 1
               else
                  dtan2 = 0._RP
               end if
               
               ! delta cycling
               if (Deltas < 1) cycle
               
!              Delta computation on "minus" side
!              ---------------------------------
               if ( faceInd_minus(1) == faceInd_plus2minus(1) ) then
                  dtan1_minus = 1._RP
               else
                  dtan1_minus = 0._RP
               end if
               if ( faceInd_minus(2) == faceInd_plus2minus(2) ) then
                  dtan2_minus = 1._RP
               else
                  dtan2_minus = 0._RP
               end if
               
!              Computation of the matrix entry
!              -------------------------------
               
               i = eq1 + i1*NCONS + j1*EtaSpa1 + k1*ZetaSpa1 ! row index (1-based)
               j = eq2 + i2*NCONS + j2*EtaSpa2 + k2*ZetaSpa2 ! column index (1-based)
               
               MatrixEntry = &
!~!
!~!              Faces contribution (surface integrals from the inner equation)
!~!              ***********************************************|**************
!~!                                                      _______|
!~!                                                      |
!~                 +(-   dfdGradQ(eq1,eq2,tanAx_plus(1) ,1,faceInd_minus2plus(1),faceInd_minus2plus(2))        &
!~                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  )  &
!~                     * spAtan1_plus  % hatD(faceInd_plus(1),faceInd_minus2plus(1))              &
!~                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan2         & ! Tangent direction 1
                     
!~                   -   dfdGradQ(eq1,eq2, normAx_plus  ,1,faceInd_minus2plus(1),faceInd_minus2plus(2))        &
!~                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus )  &
!~                     * spAnorm_plus  % bd  (elInd_plus ( normAx_plus  ), normAxSide_plus  ) * dtan1 * dtan2 & ! Normal direction
                   
!~                   -   dfdGradQ(eq1,eq2,tanAx_plus(2) ,1,faceInd_minus2plus(1),faceInd_minus2plus(2))        &
!~                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  )  &
!~                     * spAtan2_plus  % hatD(faceInd_plus(2),faceInd_minus2plus(2))              &
!~                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan1         & ! Tangent direction 2
!
!              Faces contribution (surface integrals from the outer equation) - PENALTY TERM IS BEING CONSIDERED IN THE INVISCID PART - TODO: Reorganize storage to put it explicitely in another place (needed for purely viscous equations)
!                 The tangent directions here are taken in the|reference frame of e⁻
!              ***********************************************|*********************
!                                                      _______|
!                                                      |
                 +(+   dfdGradQ(eq1,eq2,tanAx_minus(1),2,faceInd_plus(1),faceInd_plus(2)) &
                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  ) &
                     * spAtan1_minus % D   (faceInd_plus2minus(1),faceInd_minus(1)) &
                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan2_minus                  &  ! Tangent direction 1
                     
                   +   dfdGradQ(eq1,eq2, normAx_minus ,2,faceInd_plus(1),faceInd_plus(2)) &
                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  ) &
                     * spAnorm_minus % vd  (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan1_minus * dtan2_minus    & ! Normal direction
                     
                   +   dfdGradQ(eq1,eq2,tanAx_minus(2),2,faceInd_plus(1),faceInd_plus(2)) &
                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  ) &
                     * spAtan2_minus % D   (faceInd_plus2minus(2),faceInd_minus(2)) &
                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan1_minus                  & ! Tangent direction 2
                                                                                    ) * e_plus % geom % invJacobian(i1,j1,k1) ! Scale with Jacobian from mass matrix
                   
               call Matrix % AddToBlockEntry (e_plus % GlobID, e_minus % GlobID, i, j, MatrixEntry)
            end do                ; end do                ; end do                ; end do
         end do                ; end do                ; end do                ; end do
         nullify(dfdGradQ)
      end if
!
!     *********
!     Finish up
!     *********
!
      nullify(spAnorm_plus, spAtan1_plus, spAtan2_plus, spAnorm_minus, spAtan1_minus, spAtan2_minus)
      
   end subroutine Local_GetOffDiagonalBlock
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  AnJac_GetNormalVec: 
!  Get normal vector for element indexing and scale with surface Jacobian
!     (Only for p-conforming meshes)
!  ----------------------------------------------------------------------
   subroutine AnJac_GetNormalVec(f, faceSide, nEl) 
      implicit none
      !-arguments-----------------------------------------------------
      type(Face)     , intent(in)    :: f
      integer        , intent(in)    :: faceSide
      real(kind=RP)  , intent(inout) :: nEl (:,0:,0:)
      !-local-variables-----------------------------------------------
      integer :: i,j, ii, jj
      !---------------------------------------------------------------
      
      select case (faceSide)
         case (LEFT)
            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
               nEl(:,i,j) = f % geom % normal(:,i,j) * f % geom % jacobian(i,j)
            end do                        ; end do
         case (RIGHT)
            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
               call leftIndexes2Right(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
               
               nEl(1:NDIM,ii,jj) = - f % geom % normal (1:NDIM,i,j) * f % geom % jacobian(i,j)
            end do                        ; end do
      end select
      
   end subroutine AnJac_GetNormalVec   
   
#endif
end module AnalyticalJacobian
