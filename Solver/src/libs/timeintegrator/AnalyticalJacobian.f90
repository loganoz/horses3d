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
#include "Includes.h"
module AnalyticalJacobian
   use SMConstants
   use ElementClass
   use HexMeshClass
   use NodalStorageClass
   use PhysicsStorage
   use Physics
   use Jacobian
   use LinearSolverClass
   use DGSEMClass
   implicit none
   
   private
   public AnalyticalJacobian_Compute
   
   type :: ElementDivMatrix_t
      real(kind=RP), allocatable :: Dx(:,:)
      real(kind=RP), allocatable :: Dy(:,:)
      real(kind=RP), allocatable :: Dz(:,:)
      logical       :: created =.false.
   end type ElementDivMatrix_t
   
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
   subroutine AnalyticalJacobian_Compute(sem,time,linsolver,BlockDiagonalized) !,Jacobian
      implicit none
      !--------------------------------------------
      type(DGSem)              , intent(inout) :: sem
      real(kind=RP)            , intent(in)    :: time
      class(GenericLinSolver_t), intent(inout) :: linsolver          !
      logical        , optional, intent(in)    :: BlockDiagonalized  !<? Construct only the block diagonal?
      !--------------------------------------------
      integer :: eID, fID  ! Element and face counters
      integer :: nnz
      integer :: nelem
      logical :: BlockDiagonal
      !--------------------------------------------
      
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
!     Get block sizes and position in matrix (using 0-based index because of PETSc) TODO:: change!! ... This can be moved elsewhere since it's needed by both the numerical and analytical Jacobians
!     ---------------------------------------------------------------------------------------------
      
      safedeallocate(ndofelm)  ; allocate (ndofelm(nelem))
      safedeallocate(firstIdx) ; allocate (firstIdx(nelem+1))
      
      firstIdx = 0
      DO eID=1, nelem
         ndofelm(eID)  = NCONS * (sem % Nx(eID)+1) * (sem % Ny(eID)+1) * (sem % Nz(eID)+1)
         if (eID>1) firstIdx(eID) = firstIdx(eID-1) + ndofelm(eID-1)
      END DO
      firstIdx(nelem+1) = firstIdx(nelem) + ndofelm(nelem)
         
      nnz = MAXVAL(ndofelm) ! 
      
!
!     Preallocate Jacobian matrix
!     ---------------------------

      call linsolver%PreallocateA(nnz)
      call linsolver%ResetA
      
!$omp parallel
!
!     ***********************
!     Volumetric contribution
!     ***********************
!

!$omp do schedule(runtime)
      do eID = 1, nelem
         call Local_VolumetricMatrix(sem % mesh % elements(eID),linsolver)
      end do
!$omp end do

!
!     ******************
!     Faces contribution
!     ******************
!
#if defined(NAVIERSTOKES)
      call ComputeNumericalFluxJacobian(sem % mesh,time,sem % externalState)
#endif
      
!
!     Add faces contribution to off-diagonal blocks
!     ---------------------------------------------
      if (.not. BlockDiagonal) call OffDiagonalBlocks(sem % mesh,linsolver)
      
!     Pre-assembly the matrix with the volumetric contribution
!     --------------------------------------------------------
      call linsolver % PreAssemblyA 
!
!     Add faces contribution to diagonal blocks
!     -----------------------------------------

      call BlockDiagonal_FacesContribution(sem % mesh,linsolver)

!$omp end parallel
      
!
!     Finish assembling matrix
!     ------------------------
      
      call linsolver%AssemblyA
      
      
   end subroutine AnalyticalJacobian_Compute

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
!  Subroutine for computing the elemental volumetric contribution to the Jacobian matrix
!
!  -> Takes the linear contribution in ElementDivMatrix and multiplies it by the Jacobian
!     of the flux terms
!  -------------------------------------------------------------------------------------
   subroutine Local_VolumetricMatrix(e,linsolver)
      implicit none
      !--------------------------------------------
      type(Element)            , intent(in)    :: e
      class(GenericLinSolver_t), intent(inout) :: linsolver
      !--------------------------------------------
      real(kind=RP), allocatable, target :: LocalMatrix(:,:)
      real(kind=RP), allocatable         :: dFdQ(:,:)
      real(kind=RP), allocatable         :: dGdQ(:,:)
      real(kind=RP), allocatable         :: dHdQ(:,:)
      type(ElementDivMatrix_t)      :: ElDivMatrix       ! Matrix containing the terms for computing the volumetric terms of the DG div
      integer              :: N(3)              ! Polynomial orders of element
      integer              :: NDOFEL            ! Number of degrees of freedom in element
      integer, allocatable :: irow_0(:)         ! Row indexes for the matrix
      integer, allocatable :: irow(:)           ! Formatted row indexes for the column assignment (small values with index -1)
      integer              :: j
      
      !--------------------------------------------
      
!
!     Definitions and allocations
!     ---------------------------
      
      N = e % Nxyz
      NDOFEL = ndofelm(e % eID)
      
      allocate(irow_0(NDOFEL),irow(NDOFEL))
      irow_0 (1:NDOFEL) = (/ (firstIdx(e % eID) + j, j=0,NDOFEL-1) /)
      
      allocate (LocalMatrix(NDOFEL,NDOFEL), &
                       dFdQ(NDOFEL,NDOFEL), &
                       dGdQ(NDOFEL,NDOFEL), &
                       dHdQ(NDOFEL,NDOFEL))
      
      ElDivMatrix = CreateDivergenceMatrices(e)
      
      call CreateElementFluxJacobians(e,dFdQ,dGdQ,dHdQ)
      
      LocalMatrix = MATMUL(ElDivMatrix % Dx,dFdQ) + &
                    MATMUL(ElDivMatrix % Dy,dGdQ) + &
                    MATMUL(ElDivMatrix % Dz,dHdQ)
      
      ! Setting matrix by columns
      do j=1, NDOFEL
         irow = irow_0
         where (LocalMatrix(:,j) < jaceps) irow = -1
         
         call linsolver % setAColumn(NDOFEL, irow_0, firstIdx(e % eID) + j-1, LocalMatrix(:,j) )
      end do
      
      deallocate (LocalMatrix,irow_0,irow)
      
   end subroutine Local_VolumetricMatrix
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
      use MeshTypes
      implicit none
      !--------------------------------------------
      type(HexMesh), intent(inout)    :: mesh
      real(kind=RP), intent(in)       :: time
      procedure(BCState_FCN)          :: externalStateProcedure
      !--------------------------------------------
      integer :: i,j, fID
      !--------------------------------------------
      
      call mesh % ProlongSolutionToFaces
      
!$omp do private(i,j) schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID) )
         
!        Get external states
!        -------------------
         
         if (f % faceType == HMESH_BOUNDARY) then
            do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
               f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
               CALL externalStateProcedure( f % geom % x(:,i,j), &
                                            time, &
                                            f % geom % normal(:,i,j), &
                                            f % storage(2) % Q(:,i,j),&
                                            f % boundaryType )
            end do               ; end do
         end if
         
         do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
!
!           Get numerical flux jacobian on the face point (i,j)
!           ---------------------------------------------------

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
            f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
            f % storage(RIGHT) % dFStar_dqF (:,:,i,j) = f % storage(LEFT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
            
         end do             ; end do
         
         end associate
      end do
!$omp end do
   
   end subroutine ComputeNumericalFluxJacobian
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the off-diagonal blocks of the Jacobian matrix
!  -----------------------------------------------------------------------------------------------
   subroutine OffDiagonalBlocks(mesh,linsolver)
      use FaceClass
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout)    :: mesh
      class(GenericLinSolver_t), intent(inout) :: linsolver
      !--------------------------------------------
      
      
      
   end subroutine OffDiagonalBlocks

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the diagonal blocks of the Jacobian matrix
!  -------------------------------------------------------------------------------------------
   subroutine BlockDiagonal_FacesContribution(mesh,linsolver)
      use FaceClass
      use MeshTypes
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout)    :: mesh
      class(GenericLinSolver_t), intent(inout) :: linsolver
      !--------------------------------------------
      integer :: fID, eID
      
      !--------------------------------------------
      
      !--------------------------------------------
      
!     Project flux Jacobian to element face
!        -> For the block diagonal, it is only needed to transfer dF/dQL to the left element and
!                                                                 dF/dQR to the right element
!     -------------------------------------
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         associate (f => mesh % faces(fID)) 
         call f % ProjectFluxJacobianToElements(LEFT ,LEFT )   ! dF/dQL to the left element 
         if (.not. f % faceType == HMESH_BOUNDARY) call f % ProjectFluxJacobianToElements(RIGHT,RIGHT)   ! dF/dQR to the right element
                                                               ! check if right element's rotation is ok for the Jacobian
         end associate
      end do
!$omp end do
      !private(irow_0,irow,NDOFEL,LocalMat,j)
!$omp do schedule(runtime)
      do eID = 1, size(mesh % elements)
         associate (e => mesh % elements(eID)) 
         
         call AddBlock(e,mesh,linsolver)
         end associate
      end do
!$omp end do
      
      
   end subroutine BlockDiagonal_FacesContribution
   
   subroutine AddBlock(e,mesh,linsolver)
      use MeshTypes
      implicit none
      type(Element)                    , intent(in) :: e
      type(HexMesh), target    , intent(inout)    :: mesh
      class(GenericLinSolver_t), intent(inout) :: linsolver
      !--------------------------------------------
      
      integer :: irow_0(ndofelm(e % eID))         ! Row indexes for the matrix
      integer :: irow(ndofelm(e % eID))           ! Formatted row indexes for the column assignment (small values with index -1)
      integer              :: j, NDOFEL
      real(kind=RP) :: LocalMat(ndofelm(e % eID),ndofelm(e % eID))
      
      NDOFEL = ndofelm(e % eID)
      irow_0 (1:NDOFEL) = (/ (firstIdx(e % eID) + j, j=0,NDOFEL-1) /)   
      
      call Local_GetDiagonalBlock(e, &
               mesh % faces(e % faceIDs(EFRONT )) % storage(e %faceSide(EFRONT )) % dFStar_dqEl(:,:,:,:,e %faceSide(EFRONT )), &
               mesh % faces(e % faceIDs(EBACK  )) % storage(e %faceSide(EBACK  )) % dFStar_dqEl(:,:,:,:,e %faceSide(EBACK  )), &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e %faceSide(EBOTTOM)) % dFStar_dqEl(:,:,:,:,e %faceSide(EBOTTOM)), &
               mesh % faces(e % faceIDs(ERIGHT )) % storage(e %faceSide(ERIGHT )) % dFStar_dqEl(:,:,:,:,e %faceSide(ERIGHT )), &
               mesh % faces(e % faceIDs(ETOP   )) % storage(e %faceSide(ETOP   )) % dFStar_dqEl(:,:,:,:,e %faceSide(ETOP   )), &
               mesh % faces(e % faceIDs(ELEFT  )) % storage(e %faceSide(ELEFT  )) % dFStar_dqEl(:,:,:,:,e %faceSide(ELEFT  )), &
               LocalMat  )
      
      do j=1, NDOFEL
         irow = irow_0
         where (LocalMat(:,j) < jaceps) irow = -1
         
         call linsolver % AddToAColumn(NDOFEL, irow_0, firstIdx(e % eID) + j-1, LocalMat(:,j) )
      end do
      
      
   end subroutine AddBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------
!  
!  ----------------------------------------------------------
   subroutine Local_GetDiagonalBlock(e,dfdq_fr,dfdq_ba,dfdq_bo,dfdq_ri,dfdq_to,dfdq_le,LocalMat)
      implicit none
      !-------------------------------------------
      type(Element)                    , intent(in) :: e
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
            
!           Add faces contribution
!           ----------------------

            i = eq1 + i1*NCONS + j1*NCONS*nXi + k1*NCONS*nEta*nXi ! row index (1-based)
            j = eq2 + i2*NCONS + j2*NCONS*nXi + k2*NCONS*nEta*nXi ! column index (1-based)
            
            LocalMat(i,j) = LocalMat(i,j) &
                                    + ( dfdq_fr(eq1,eq2,i1,k1) * e % spAeta % b(j1,FRONT ) * e % spAeta % v(j2,FRONT ) * di * dk   & ! 1 Front
                                    +   dfdq_ba(eq1,eq2,i1,k1) * e % spAeta % b(j1,BACK  ) * e % spAeta % v(j2,BACK  ) * di * dk   & ! 2 Back
                                    +   dfdq_bo(eq1,eq2,i1,j1) * e % spAZeta% b(k1,BOTTOM) * e % spAZeta% v(k2,BOTTOM) * di * dj   & ! 3 Bottom
                                    +   dfdq_to(eq1,eq2,i1,j1) * e % spAZeta% b(k1,TOP   ) * e % spAZeta% v(k2,TOP   ) * di * dj   & ! 5 Top
                                    +   dfdq_ri(eq1,eq2,j1,k1) * e % spAXi  % b(i1,RIGHT ) * e % spAXi  % v(i2,RIGHT ) * dj * dk   & ! 4 Right
                                    +   dfdq_le(eq1,eq2,j1,k1) * e % spAXi  % b(i1,LEFT  ) * e % spAXi  % v(i2,LEFT  ) * dj * dk ) & ! 6 Left
                                                                                                / e % geom % jacobian(i1,j1,k1) ! Scale with Jacobian from mass matrix
            
         end do                ; end do                ; end do                ; end do
      end do                ; end do                ; end do                ; end do
      
   end subroutine Local_GetDiagonalBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------
!  Subroutine for computing the global flux Jacobian matrices
!  ----------------------------------------------------------
   subroutine CreateElementFluxJacobians(e,dFdQ,dGdQ,dHdQ)
      implicit none
      !-------------------------------------------
      type(Element), intent(in)    :: e
      real(kind=RP), intent(inout) :: dFdQ(:,:) ! Global Jacobian of the flux in the x-direction
      real(kind=RP), intent(inout) :: dGdQ(:,:) ! Global Jacobian of the flux in the y-direction
      real(kind=RP), intent(inout) :: dHdQ(:,:)
      !-------------------------------------------
      real(kind=RP) :: dfdq_(NCONS,NCONS) ! Local Jacobian of the flux in the x-direction
      real(kind=RP) :: dgdq_(NCONS,NCONS) ! Local Jacobian of the flux in the y-direction
      real(kind=RP) :: dhdq_(NCONS,NCONS) ! Local Jacobian of the flux in the z-direction
      integer       :: Idx0               ! Start index to add a block matrix
      integer       :: i,j,k              ! Coordinate counters
      integer       :: N(3)               ! Polynomial orders of element
      !-------------------------------------------
      
      N = e % Nxyz
      dFdQ = 0._RP
      dGdQ = 0._RP
      dHdQ = 0._RP
#if defined(NAVIERSTOKES)
      do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
         Idx0 = ijk2local(i,j,k,1,NCONS,N(1),N(2),N(3)) ! Beginning index of local block matrix
         
         call InviscidJacobian(e % storage % Q(:,i,j,k),dfdq_,dgdq_,dhdq_)
         
         dFdQ(Idx0:Idx0+NCONS-1,Idx0:Idx0+NCONS-1) = dfdq_
         dGdQ(Idx0:Idx0+NCONS-1,Idx0:Idx0+NCONS-1) = dgdq_
         dHdQ(Idx0:Idx0+NCONS-1,Idx0:Idx0+NCONS-1) = dhdq_
         
      end do ; end do ; end do
#endif
   end subroutine CreateElementFluxJacobians
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------------
!  Subroutine for computing and storing a divergence matrix for the element level
!  ------------------------------------------------------------------------------
   function CreateDivergenceMatrices(e) result(ElementDivMatrix)
      implicit none
      !-------------------------------------------
      type(Element)           , intent(in)  :: e
      type(ElementDivMatrix_t)              :: ElementDivMatrix
      !-------------------------------------------
      integer :: DIMi
      !-------------------------------------------
      
!
!     Create a derivative matix for each direction
!     --------------------------------------------

      call CreateMatrix_StdVolumeGreen(ElementDivMatrix % Dx, e, 1)
      call CreateMatrix_StdVolumeGreen(ElementDivMatrix % Dy, e, 2)
      call CreateMatrix_StdVolumeGreen(ElementDivMatrix % Dz, e, 3)
      
!
!     All done
!     --------

      ElementDivMatrix % created = .true.
      
   end function CreateDivergenceMatrices
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
!  Subroutine for computing the elemental volumetric contribution to the Jacobian matrix
!
!  -> This matrix can be created at the beginning of the simulation since it gathers all
!     the linear contributions to the weak Green integral. The nonlinearities are introduced
!     in runtime taking into account the Jacobians specified in libs/physics/Physics.f90
!  -> Currently, the subroutine takes advantage of the tensor product expansion and loads
!     directional matrices by blocks.
!  -------------------------------------------------------------------------------------
   subroutine CreateMatrix_StdVolumeGreen(LocalMatrix,e,DIMi)
      implicit none
      !-------------------------------------------
      real(kind=RP), allocatable, target, intent(inout) :: LocalMatrix(:,:)
      type(Element)                     , intent(in)    :: e
      integer                           , intent(in)    :: DIMi
      !-------------------------------------------
      real(kind=RP), pointer  :: Mat_p(:,:) ! Pointer to blocks of the elemental matrix
      
      integer :: N(3)              ! Polynomial orders of element
      integer :: NDOFEL            ! Number of degrees of freedom in element
      integer :: iEQ               ! Equation counter
      integer :: i,j,k             ! Coordinate counters
      integer :: Idx0              ! Start index to add a block matrix
      
      integer :: dIdx_Xi   ! Spacing between matrix entries for Xi bases
      integer :: dIdx_Eta  ! Spacing between matrix entries for Eta bases
      integer :: dIdx_Zeta ! Spacing between matrix entries for Zeta bases
      integer :: Size_Xi   ! Size of block matrix for Xi bases
      integer :: Size_Eta  ! Size of block matrix for Eta bases
      integer :: Size_Zeta ! Size of block matrix for Zeta bases
      !--------------------------------------------
      
      if (DIMi < 1 .or. DIMi>3) ERROR STOP 'Wrong dimension in CreateMatrix_StdVolumeGreen'
      
      N = e % Nxyz
      NDOFEL = NCONS * (N(1) + 1) * (N(2) + 1) * (N(3) + 1)
      
      allocate (LocalMatrix(NDOFEL,NDOFEL))
      LocalMatrix = 0._RP
      
!
!     General definitions for the local Jacobian matrix
!     -------------------------------------------------
      
      dIdx_Xi   = NCONS
      dIdx_Eta  = NCONS * (N(1) + 1)
      dIdx_Zeta = NCONS * (N(1) + 1) * (N(2) + 1)
      
      Size_Xi   = NCONS * (N(1) + 1)
      Size_Eta  = NCONS * (N(1) + 1) * (N(2) + 1)
      Size_Zeta = NCONS * (N(1) + 1) * (N(2) + 1) * (N(3) + 1)
      
!
!      Xi-contribution
!     --------------------------------------------
      
      do k = 0, N(3) ; do j = 0, N(2) ; do iEQ = 1, NCONS
         
         Idx0 = ijk2local(0,j,k,iEQ,NCONS,N(1),N(2),N(3)) ! Beginning index of local block matrix
         
         Mat_p =>  LocalMatrix(Idx0:(Idx0-1) + Size_Xi:dIdx_Xi, &
                               Idx0:(Idx0-1) + Size_Xi:dIdx_Xi)
         
         Mat_p = Mat_p + MultiplyRowsByVec (e % spAxi % hatD, e % geom % jGradXi(DIMi,:,j,k)) 
         
      end do ; end do ; end do
      
!
!     Eta-contribution
!     --------------------------------------------
      
      do k = 0, N(3) ; do i = 0, N(1) ; do iEQ = 1, NCONS
         
         Idx0 = ijk2local(i,0,k,iEQ,NCONS,N(1),N(2),N(3)) ! Beginning index of local block matrix
         
         Mat_p =>  LocalMatrix(Idx0:(Idx0-1) + Size_Eta:dIdx_Eta, &
                               Idx0:(Idx0-1) + Size_Eta:dIdx_Eta)
         
         Mat_p = Mat_p + MultiplyRowsByVec (e % spAeta % hatD, e % geom % jGradEta(DIMi,i,:,k)) 
         
      end do ; end do ; end do
      
!
!     Zeta-contribution
!     --------------------------------------------

      do j = 0, N(2) ; do i = 0, N(1) ; do iEQ = 1, NCONS
         
         Idx0 = ijk2local(i,j,0,iEQ,NCONS,N(1),N(2),N(3)) ! Beginning index of local block matrix
         
         Mat_p =>  LocalMatrix(Idx0:(Idx0-1) + Size_Zeta:dIdx_Zeta, &
                               Idx0:(Idx0-1) + Size_Zeta:dIdx_Zeta)
         
         Mat_p = Mat_p + MultiplyRowsByVec (e % spAzeta % hatD, e % geom % jGradZeta(DIMi,i,j,:)) 
         
      end do ; end do ; end do
      
!
!     Scale with Jacobian (from mass matrix terms)
!     --------------------------------------------
      
      do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1) ; do iEQ = 1, NCONS
         Idx0 = ijk2local(i,j,k,iEQ,NCONS,N(1),N(2),N(3))
         
         LocalMatrix (Idx0,:) = LocalMatrix (Idx0,:) / e % geom % jacobian(i,j,k)
      end do ; end do ; end do ; end do
      
      nullify (Mat_p)
      
   end subroutine CreateMatrix_StdVolumeGreen
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
!  Subroutine for multipĺyng the rows of a matrix by an array... TODO: Move to DenseMat_Utils
!  -------------------------------------------------------------------------------------
   function MultiplyRowsByVec(Mat,Vec) result(NewMat)
      implicit none
      !-------------------------------------------
      real(kind=RP), intent(in)  :: Mat(:,:)
      real(kind=RP), intent(in)  :: Vec(size(Mat,2))
      real(kind=RP)              :: NewMat(size(Mat,1),size(Mat,2))
      !-------------------------------------------
      integer :: j
      !-------------------------------------------
      
      do j = 1, size(Mat,2)
         NewMat(:,j) = Mat(:,j) * Vec(j)
      end do
   end function MultiplyRowsByVec
   

end module AnalyticalJacobian
