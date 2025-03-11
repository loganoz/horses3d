!//////////////////////////////////////////////////////
!
!  This module provides the routines for computing the analytical Jacobian matrix
!  -> Only for p-conforming representations (TODO: make general)
!  -> The Jacobian of the BCs is temporarily computed numerically
!
!//////////////////////////////////////////////////////


! /////////////////////////////////////////////////////////////////////////////
! ----------------------ZhangYu--------------------
! For irragular mesh with the physical boundarys of element are not ortharganal, 
! the transform Jacobian should provide scaling and rotating effects, 
! then the Jacobian is not symmertic.  
! If the Jacobian is not symmetic, how to make the coefficient matrix symmetric?
! If the Jacobian matrix is unsymmetric, 
! the coefficient matrix should not be symmetric.
! /////////////////////////////////////////////////////////////////////////////

#include "Includes.h"
module AnalyticalJacobianPoisson
    use SMConstants
    use ElementClass
    use HexMeshClass
    use NodalStorageClass
    use PhysicsStorage
    use Physics
    use JacobianDefinitions, only: JACEPS
    use JacobianComputerClass, only: JacobianComputer_t
    use MatrixClass
    use DGSEMClass 
    ! use DGSEMClass, only: DGSem, ComputeTimeDerivative_f
    use StopWatchClass
    use MeshTypes
    use EllipticDiscretizations
    use MPI_Process_Info, only: MPI_Process
    use ElementConnectivityDefinitions, only: axisMap, normalAxis
    use BoundaryConditions, only: BCs
    use FaceClass, only: Face
    use Utilities, only: dot_product
    use ConnectivityClass, only: Connectivity
    use FTValueDictionaryClass
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
    use RiemannSolvers_NS, only: RiemannSolver_dFdQ, RiemannSolver
    use HyperbolicDiscretizations, only: HyperbolicDiscretization
    use VariableConversion, only: NSGradientVariables_STATE
    use FluidData_NS, only: dimensionless
#elif defined(SPALARTALMARAS)
    use RiemannSolvers_NSSA, only: RiemannSolver_dFdQ, RiemannSolver
    use HyperbolicDiscretizations, only: HyperbolicDiscretization
    use VariableConversion, only: NSGradientVariables_STATE
    use FluidData_NSSA, only: dimensionless
#elif defined(SCALAR)
    use RiemannSolvers_SLR 
    ! use HyperbolicDiscretizations, only: HyperbolicDiscretization
    use VariableConversion
    use FluidData_SLR
#elif defined(SCALAR_INS_V04)
    use RiemannSolvers_SLR_INS_V04
    ! use HyperbolicDiscretizations, only: HyperbolicDiscretization
    use VariableConversion
    use FluidData_SLR_INS_V04
#endif
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

    private
    public AnPoissonJacobian_t

    integer, parameter :: other(2) = [2, 1]

!
!  **************************************************
!  Main type for the analytical Jacobian computations
!  **************************************************
    type, extends(JacobianComputer_t) :: AnPoissonJacobian_t

    contains
        procedure :: Construct => AnPoissonJacobian_Construct
        procedure :: Compute => AnPoissonJacobian_Compute
        
#if defined(SCALAR_INS_V04)
        procedure :: TimeDerivative_ComputeRHS
        procedure :: TimeDerivative_ComputeRHS_pressure_specific
        procedure :: TimeDerivative_ComputeRHS_velocity_specific
#endif
        ! ===================================================================
        ! procedure :: Local_Get_BC_Dir_RHS
    end type AnPoissonJacobian_t

#if defined(NAVIERSTOKES)
    real(kind=RP) :: Identity(NCONS, NCONS) ! identity matrix. TODO: Define only once in the constructor (When this is a class!!)
#elif  defined(SCALAR) || defined(SCALAR_INS_V04)
    real(kind=RP) :: Identity(NCONS, NCONS) ! identity matrix. TODO: Define only once in the constructor (When this is a class!!)
#endif
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  AnPoissonJacobian_Construct:
!  Main class constructor
!  -------------------------------------------------------
    subroutine AnPoissonJacobian_Construct(this, mesh, nEqn, controlVariables)
        implicit none
        !-arguments-----------------------------------------
        class(AnPoissonJacobian_t), intent(inout) :: this
        type(HexMesh), intent(inout) :: mesh
        integer, intent(in)    :: nEqn
        type(FTValueDictionary), intent(in)    :: controlVariables
        !-local-variables-----------------------------------
        integer :: fID, eID
        !---------------------------------------------------

!
!     Construct parent
!     ----------------
        call this%JacobianComputer_t%construct(mesh, nEqn, controlVariables)
!
!     Create stopwatch event
!     ----------------------
        call Stopwatch%CreateNewEvent("Analytical Jacobian Poisson construction")
!
!     Construct specific storage
!     --------------------------
        ! Global storage
        mesh%storage%AnPoissonJacobian = .TRUE.
        ! Elements:
!$omp parallel do schedule(runtime)
        do eID = 1, mesh%no_of_elements
            ! write (*,*) "bug in analytic ", 1./(eID-1)
            call mesh%storage%elements(eID)%ConstructAnJac()
            
        end do

!$omp end parallel do


#if defined(SCALAR_INS_V04)
        write (*,*) "NCONS, N_INS, NDIM, NDIM, self % Nxyz(1), self % Nxyz(2), self % Nxyz(3)",&
                     NCONS, N_INS, NDIM, NDIM
#endif
        
        ! Faces:
!$omp parallel do schedule(runtime)
        do fID = 1, size(mesh%faces)
            call mesh%faces(fID)%storage%ConstructAnJac(NDIM)
        end do
!$omp end parallel do

        !TODO: Add conformity check

    end subroutine AnPoissonJacobian_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  Subroutine for computing the analytical Jacobian matrix
!  -------------------------------------------------------
   subroutine AnPoissonJacobian_Compute(this, sem, nEqn, time, Matrix,               &
    TimeDerivative, TimeDerivativeIsolated, &
    eps_in, BlockDiagonalized, mode  &
#if defined(SCALAR_INS_V04)
    ,startNum                     &        
#endif
    )
        implicit none
        !-arguments----------------------------------
        class(AnPoissonJacobian_t), intent(inout)     :: this
        type(DGSem), intent(inout)     :: sem
        integer, intent(in)        :: nEqn
        real(kind=RP), intent(in)        :: time
        class(Matrix_t), intent(inout)     :: Matrix

        real(kind=RP), optional, intent(in)        :: eps_in            ! Not needed here...
        logical, optional, intent(in)        :: BlockDiagonalized !<? Construct only the block diagonal?
        integer, optional, intent(in)        :: mode
#if defined(SCALAR_INS_V04)
        integer, optional, intent(in)        :: startNum   
        procedure(ComputeNonlinearStep1_f), optional :: TimeDerivative    ! Not needed here...
        procedure(ComputeNonlinearStep1_f), optional :: TimeDerivativeIsolated
#else            
    procedure(ComputeTimeDerivative_f), optional :: TimeDerivative    ! Not needed here...
    procedure(ComputeTimeDerivative_f), optional :: TimeDerivativeIsolated
#endif

        ! =======================================================
#if defined(SCALAR)
        !--------------------------------------------
        integer :: nnz
        integer :: nelem, ierr
        logical :: BlockDiagonal
        integer :: eID, i
        integer :: j, k!--------------------------------------------

        call Stopwatch%Start("Analytical Jacobian construction")
        !
        !     Initializations
        !     ---------------
        if (present(BlockDiagonalized)) then
            BlockDiagonal = BlockDiagonalized
        else
            BlockDiagonal = .FALSE.
        end if

        nelem = size(sem%mesh%elements)

        ! write (*,*) "AnPoissonJacobian_Compute, nelem = size(sem%mesh%elements), this%ndofelm_l ====== ",&
                    !  nelem,this%ndofelm_l

        if (.not. this%preAllocate) then
            select type (Matrix_p => Matrix)
                !
                !           If block-diagonal matrix, construct with size of blocks
                !           -------------------------------------------------------
            type is (DenseBlockDiagMatrix_t)
                call Matrix_p%Preallocate(nnzs=this%ndofelm_l)
            type is (SparseBlockDiagMatrix_t)
                call Matrix_p%Preallocate(nnzs=this%ndofelm_l)

            type is (csrMat_t)
                call Matrix_p%Preallocate()
                !           ----------------------------------------------
            class default
                if (BlockDiagonal) then
                    nnz = maxval(this%ndofelm_l)
                else
                    nnz = 7*maxval(this%ndofelm_l) ! 20180201: hard-coded to 7, since the analytical Jacobian can only compute off-diagonal blocks for compact schemes (neighbors' effect)
                end if
                call Matrix_p%Preallocate(nnz)
            end select

            call Matrix%SpecifyBlockInfo(this%firstIdx, this%ndofelm)
        end if

        call Matrix%Reset(ForceDiagonal=.TRUE.)

!$omp parallel
        !     ------------------
        !     Prolong Q to faces
        !     ------------------
    !   call sem%mesh%ProlongSolutionToFaces(NCONS)
    !   call IP_Poisson_ComputeGradientCoe(nEqn, nEqn, sem % mesh)
    !   write (*,*) "=self % dF_dgradQ=  ---=-=-", size(self % dF_dgradQ)


!     ***************
!     Diagonal blocks
!     ***************
!
      call AnalyticalJacobianPoisson_DiagonalBlocks(sem % mesh, nEqn,time, Matrix)

      if (.not. BlockDiagonal) call AnalyticalJacobianPoisson_OffDiagonalBlocks(sem % mesh,nEqn ,Matrix)
      !$omp end parallel
!
!     ************************
!     Finish assembling matrix
!     ************************
!
        call Matrix%Assembly()

#elif defined(SCALAR_INS_V04)
! ============================================================================================
        !--------------------------------------------
        integer :: nnz
        integer :: nelem, ierr
        logical :: BlockDiagonal
        integer :: eID, i
        integer :: j, k!--------------------------------------------

        call Stopwatch%Start("Analytical Jacobian construction")
        !
        !     Initializations
        !     ---------------
        if (present(BlockDiagonalized)) then
            BlockDiagonal = BlockDiagonalized
        else
            BlockDiagonal = .FALSE.
        end if

        nelem = size(sem%mesh%elements)

        ! write (*,*) "AnPoissonJacobian_Compute, nelem = size(sem%mesh%elements), this%ndofelm_l ====== ",&
                    !  nelem,this%ndofelm_l

        if (.not. this%preAllocate) then
            select type (Matrix_p => Matrix)
                !
                !           If block-diagonal matrix, construct with size of blocks
                !           -------------------------------------------------------
            type is (DenseBlockDiagMatrix_t)
                call Matrix_p%Preallocate(nnzs=this%ndofelm_l)
            type is (SparseBlockDiagMatrix_t)
                call Matrix_p%Preallocate(nnzs=this%ndofelm_l)

            type is (csrMat_t)
                call Matrix_p%Preallocate()
                !           ----------------------------------------------
            class default
                if (BlockDiagonal) then
                    nnz = maxval(this%ndofelm_l)
                else
                    nnz = 7*maxval(this%ndofelm_l) ! 20180201: hard-coded to 7, since the analytical Jacobian can only compute off-diagonal blocks for compact schemes (neighbors' effect)
                end if
                call Matrix_p%Preallocate(nnz)
            end select

            call Matrix%SpecifyBlockInfo(this%firstIdx, this%ndofelm)
        end if

        call Matrix%Reset(ForceDiagonal=.TRUE.)

!$omp parallel
 
!     ***************
!     Diagonal blocks
!     ***************
!
      call AnalyticalJacobianPoisson_DiagonalBlocks(sem % mesh, nEqn,time, Matrix)

      if (.not. BlockDiagonal) call AnalyticalJacobianPoisson_OffDiagonalBlocks(sem % mesh,nEqn ,Matrix, startNum)
!$omp end parallel
!
    if (startNum == 4) then
        ! call AnalyticalJacobianPoisson_Neu_Per_setOne(sem % mesh,nEqn ,Matrix, startNum)
      call AnalyticalJacobianPoisson_Neu_Per_setOne(sem % mesh, nEqn,time, Matrix)
    endif


!     ************************
!     Finish assembling matrix
!     ************************
!
        call Matrix%Assembly()
!
#elif defined(NAVIERSTOKES)
        
#else
        error stop ':: Analytical Jacobian only for NS, SCALAR_INS_V04'
#endif

        call Matrix%Visualize('Jacobian.txt')
    end subroutine AnPoissonJacobian_Compute


#if defined(SCALAR_INS_V04)

    subroutine TimeDerivative_ComputeRHS(this,  mesh, nEqn,t, RHS, AllNDOF )
        implicit none
        class(AnPoissonJacobian_t), intent(inout)     :: this
        type(HexMesh), target    , intent(inout) :: mesh
     !   type(HexMesh)              :: mesh
        integer, intent(in)    :: AllNDOF    !
        real(kind=RP)              :: t
        integer, intent(in)        :: nEqn
        ! real(kind=RP)                       :: RHS( size(e_plus % storage % Q) )
    real(kind=RP)                       :: RHS(nEqn * AllNDOF)
        !
!        ---------------
!        Local variables
!        ---------------
!
        integer     :: eID , i, j, k, ierr, fID

        type(Element), pointer :: e
        type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
        integer ::  elSide, side
        ! type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Face)   , pointer :: f

!$omp do schedule(runtime) private(i,j,k)
        do eID = 1, size(mesh % elements) 
           associate(e => mesh % elements(eID)) 
           if ( e % hasSharedFaces ) cycle

           spAxi => NodalStorage(e%Nxyz(1))
           spAeta => NodalStorage(e%Nxyz(2))
           spAzeta => NodalStorage(e%Nxyz(3))

  !!        One block for every neighbor element
   !        ------------------------------------
           do elSide = 1, 6
              ! write (*,*) "eID, fe_plus % NumberOfConnections(elSide) ======= ", eID, e_plus % NumberOfConnections(elSide)
              if (e % NumberOfConnections(elSide) == 0) then 

                  fID  = e % faceIDs(elSide)
                  side = e % faceSide(elSide)

                  write (*,*) "fID, side ====BCBCBCBC RHSRHSRHS====**********==**********========= ", fID, side

                  f => mesh % faces(fID)


                  if(f%FaceType == HMESH_BOUNDARY) then
                      call Local_Get_BC_Dir_RHS_Poisson(f,nEqn,e,t,side, RHS, AllNDOF, 0)

                      
                  endif
          
              endif
         
          end do
           
           end associate 

        end do
!$omp end do
    end subroutine TimeDerivative_ComputeRHS


    subroutine TimeDerivative_ComputeRHS_pressure_specific(this,  mesh, nEqn,t, RHS, AllNDOF,  startVarNum )
        implicit none
        class(AnPoissonJacobian_t), intent(inout)     :: this
        type(HexMesh), target    , intent(inout) :: mesh
     !   type(HexMesh)              :: mesh
        integer, intent(in)    :: AllNDOF    !
        real(kind=RP)              :: t
        integer, intent(in)        :: nEqn
        ! real(kind=RP)                       :: RHS( size(e_plus % storage % Q) )
        real(kind=RP)               :: RHS(nEqn * AllNDOF)
        ! REAL(kind=RP), intent(inout)       :: variable(nEqn, dim1, dim2,dim3)
        integer, intent(in)             :: startVarNum
        !
!        ---------------
!        Local variables
!        ---------------
!
        integer     :: eID , i, j, k, ierr, fID
        integer     :: dim1, dim2,dim3

        type(Element), pointer :: e
        integer ::  elSide, side
        ! type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Face)   , pointer :: f

!$omp do schedule(runtime) private(i,j,k)
        do eID = 1, size(mesh % elements) 
           associate(e => mesh % elements(eID)) 
           if ( e % hasSharedFaces ) cycle

           dim1 = e%Nxyz(1)
           dim2 = e%Nxyz(2)
           dim3 = e%Nxyz(3)

  !!        One block for every neighbor element
   !        ------------------------------------
           do elSide = 1, 6
                ! write (*,*) "eID, fe_plus % NumberOfConnections(elSide) ======= ", eID, e_plus % NumberOfConnections(elSide)
                if (e % NumberOfConnections(elSide) == 0) then 

                    fID  = e % faceIDs(elSide)
                    side = e % faceSide(elSide)

                    ! write (*,*) "fID, side ====BCBCBCBC RHSRHSRHS====**********==**********========= ", fID, side
                    ! write (*,*) "fID, side ====BCBCBCBC dim 1, 2, 3*******==**** ==== ",dim1, dim2,dim3

                    f => mesh % faces(fID)

                    if(f%FaceType == HMESH_BOUNDARY) then

                        if (TRIM(BCs(f % zone) % bc % BCType) == "zydiyDirichletBC") then

                            call Local_Get_BC_Neumann_RHS_Poisson_Specific_variable(f,nEqn,e,t,side, RHS, AllNDOF, &
                                                                                e % storage % pre_source(:,:,:,:),  & 
                                                                                startVarNum, dim1, dim2,dim3)
                            
                            ! print *, "Boundary condition matches: zydiyDirichletBC"

                        else if (TRIM(BCs(f % zone) % bc % BCType) == "zydiyfreeslipwall") then

                            call Local_Get_BC_Dir_RHS_Poisson_Specific_variable(f,nEqn,e,t,side, RHS, AllNDOF, &
                                                                            e % storage % pre_source(:,:,:,:),  & 
                                                                            startVarNum, dim1, dim2,dim3)

                        else if (TRIM(BCs(f % zone) % bc % BCType) == "user-defined") then

                            call Local_Get_BC_Dir_RHS_Poisson_Specific_variable(f,nEqn,e,t,side, RHS, AllNDOF, &
                                                                            e % storage % pre_source(:,:,:,:),  & 
                                                                            startVarNum, dim1, dim2,dim3)
                        
                        end if


                        
                    endif
                endif
            end do
           
            end associate 
        end do
!$omp end do
     end subroutine TimeDerivative_ComputeRHS_pressure_specific

    subroutine TimeDerivative_ComputeRHS_velocity_specific(this,  mesh, nEqn,t, RHS, AllNDOF,  startVarNum )
        implicit none
        class(AnPoissonJacobian_t), intent(inout)     :: this
        type(HexMesh), target    , intent(inout) :: mesh
     !   type(HexMesh)              :: mesh
        integer, intent(in)    :: AllNDOF    !
        real(kind=RP)              :: t
        integer, intent(in)        :: nEqn
        real(kind=RP)               :: RHS(nEqn * AllNDOF)
        integer, intent(in)             :: startVarNum
        !
!        ---------------
!        Local variables
!        ---------------
!
        integer     :: eID , i, j, k, ierr, fID
        integer     :: dim1, dim2,dim3

        type(Element), pointer :: e
        integer ::  elSide, side
        type(Face)   , pointer :: f

!$omp do schedule(runtime) private(i,j,k)
        do eID = 1, size(mesh % elements) 
           associate(e => mesh % elements(eID)) 
           if ( e % hasSharedFaces ) cycle

           dim1 = e%Nxyz(1)
           dim2 = e%Nxyz(2)
           dim3 = e%Nxyz(3)

  !!        One block for every neighbor element
   !        ------------------------------------
           do elSide = 1, 6
                if (e % NumberOfConnections(elSide) == 0) then 

                    fID  = e % faceIDs(elSide)
                    side = e % faceSide(elSide)


                    f => mesh % faces(fID)

                    if(f%FaceType == HMESH_BOUNDARY) then

                        if (TRIM(BCs(f % zone) % bc % BCType) == "zydiyDirichletBC") then

                            call Local_Get_BC_Dir_RHS_Poisson_3velocity(f,nEqn,e,t,side, RHS, AllNDOF, &
                                                                        e % storage % vel_source(:,:,:,:),  & 
                                                                        startVarNum, dim1, dim2,dim3)
                            
                            ! print *, "Boundary condition matches: zydiyDirichletBC"

                        else if (TRIM(BCs(f % zone) % bc % BCType) == "zydiyfreeslipwall") then

                            call Local_Get_BC_Nuemann_RHS_Poisson_3velocity(f,nEqn,e,t,side, RHS, AllNDOF, &
                                                                        e % storage % vel_source(:,:,:,:),  & 
                                                                        startVarNum, dim1, dim2,dim3)
                            ! print *, "Boundary condition does match --- zydiyfreeslipwall."

                        else if (TRIM(BCs(f % zone) % bc % BCType) == "user-defined") then

                            call Local_Get_BC_Dir_RHS_Poisson_3velocity(f,nEqn,e,t,side, RHS, AllNDOF, &
                                                                        e % storage % vel_source(:,:,:,:),  & 
                                                                        startVarNum, dim1, dim2,dim3)
                            ! print *, "Boundary condition does match --- user-defined."
                        
                        end if

                    endif
                endif
            end do
           
            end associate 
        end do
!$omp end do
     end subroutine TimeDerivative_ComputeRHS_velocity_specific



#endif

#if defined(SCALAR)

 


#endif

! ===========================================================================================================
! ===========================================================================================================
! ===========================================================================================================
! ===========================================================================================================

#if defined(SCALAR_INS_V04)

! -------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the diagonal blocks of the Jacobian matrix
!  -------------------------------------------------------------------------------------------
    subroutine AnalyticalJacobianPoisson_DiagonalBlocks(mesh, nEqn, time, Matrix)
        implicit none
        !--------------------------------------------
        type(HexMesh), target, intent(inout) :: mesh
        integer, intent(in)    :: nEqn
        real(kind=RP), intent(in)    :: time
        class(Matrix_t), intent(inout) :: Matrix
        !--------------------------------------------
        integer :: eID, fID
        type(Element), pointer :: e
        !--------------------------------------------
!
!     Compute each element's diagonal block
!     -------------------------------------
        write (*,*) "N_INS ====================== ", N_INS
!$omp do schedule(runtime) private(e)
        do eID = 1, size(mesh%elements)
            ! write (*,*) "eID ============= ", eID
            e => mesh%elements(eID)
            ! write (*, *) "F B B R T L ",e%faceIDs(EFRONT), e%faceIDs(EBACK), e%faceIDs(EBOTTOM), & 
            !                             e%faceIDs(ERIGHT), e%faceIDs(ETOP), e%faceIDs(ELEFT)
            call Local_SetDiagonalBlock(e, neqn,&
                                        mesh%faces(e%faceIDs(EFRONT)), &
                                        mesh%faces(e%faceIDs(EBACK)), &
                                        mesh%faces(e%faceIDs(EBOTTOM)), &
                                        mesh%faces(e%faceIDs(ERIGHT)), &
                                        mesh%faces(e%faceIDs(ETOP)), &
                                        mesh%faces(e%faceIDs(ELEFT)), &
                                        Matrix, &
                                        mesh%IBM)
        end do
!$omp end do
        nullify (e)
    end subroutine AnalyticalJacobianPoisson_DiagonalBlocks


    !  -----------------------------------------------------
!  Local_SetDiagonalBlock:
!     Subroutine to set a diagonal block in the Jacobian
!  -----------------------------------------------------
    subroutine Local_SetDiagonalBlock(e, neqn, fF, fB, fO, fR, fT, fL, Matrix, IBM)
        use IBMClass
        implicit none
        !-------------------------------------------
        type(Element), intent(inout) :: e
        type(Face), target, intent(in)    :: fF, fB, fO, fR, fT, fL !< The six faces of the element
        class(Matrix_t), intent(inout) :: Matrix
        type(IBM_type), intent(inout) :: IBM
        integer, intent(in)    :: nEqn
        !-------------------------------------------
        real(kind=RP) :: MatEntries(neqn, neqn)
        integer :: i, j             ! Matrix indices
        integer :: i1, j1, k1, eq1  ! variable counters for row
        integer :: i2, j2, k2, eq2  ! variable counters for column
        integer :: i12, j12, k12    ! variable counters when (row index) = (column index)
        integer :: baseRow, baseCol ! Position of neqn by neqn miniblock of Jacobian
        integer :: nXi, nEta        ! Number of nodes in every direction
        integer :: EtaSpa, ZetaSpa  ! Spacing for these two coordinate directions
        real(kind=RP) :: temp


        integer :: index  ! variable counters for column

        type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

        optional :: IBM

        !-------------------------------------------
        spAxi => NodalStorage(e%Nxyz(1))
        spAeta => NodalStorage(e%Nxyz(2))
        spAzeta => NodalStorage(e%Nxyz(3))
!
!     *******************
!     Initial definitions
!     *******************
!
        nXi = e%Nxyz(1) + 1
        nEta = e%Nxyz(2) + 1
        EtaSpa = neqn*nXi
        ZetaSpa = neqn*nXi*nEta


        ! write(*,*) "e%geom%Jacobian ============== ",  e%geom%Jacobian
        ! write(*,*) "e%geom%jGradXi ============== ",  e%geom%jGradXi
        ! write(*,*) "e%geom%jGradEta ============== ",  e%geom%jGradEta
        ! write(*,*) "e%geom%jGradzeta ============== ",  e%geom%jGradzeta



        !     Xi contribution (dj*dk)
        do k12 = 0, e%Nxyz(3)
            do j12 = 0, e%Nxyz(2)
                do i2 = 0, e%Nxyz(1)
                    do i1 = 0, e%Nxyz(1)
                        MatEntries = 0._RP

                        temp = 0.0_RP
                        ! temp = 10.0_RP/(i1-4)   ! a bug made by ZY on purpose

                        do index = 0, e%Nxyz(1)
                            temp = temp &
                                        + spAxi%DT(i1   ,index)   &
                                        * spAxi%DT(i2   ,index)   &
                                        * spAxi%w(index)       &
                                        * DOT_PRODUCT(                     &
                                            e%geom%jGradXi(:,index,j12,k12)/ e%geom%Jacobian(index, j12, k12),&
                                            e%geom%jGradXi(:,index,j12,k12) &
                                            ! e%geom%jGradXi(:,index,j12,k12)/ e%geom%Jacobian(index, j12, k12) &
                                        )
                                                    
                        end do

                        do index = 1, neqn
                            MatEntries(index,index)=temp &
                                                / spAXi%w(i1) / e%geom%Jacobian(i1, j12, k12)
                        end do


                        MatEntries = MatEntries 

            
                        baseRow = i1*neqn + j12*EtaSpa + k12*ZetaSpa
                        baseCol = i2*neqn + j12*EtaSpa + k12*ZetaSpa

                        do eq2 = 1, neqn
                            do eq1 = 1, neqn
                                i = eq1 + baseRow  ! row index (1-based)
                                j = eq2 + baseCol  ! column index (1-based)

                                call Matrix%AddToBlockEntry(e%GlobID, e%GlobID, i, j, MatEntries(eq1, eq2))  
                                ! TODO: This can be improved by setting the whole matrix at once

                            end do
                        end do
                        ! =========--------===*****************************======================
                    end do
                end do
            end do
        end do


        !     Eta contribution (di*dk)
!     ------------------------
        do k12 = 0, e % Nxyz(3) 
            do j2 = 0, e % Nxyz(2) 
                do j1 = 0, e % Nxyz(2) 
                    do i12 = 0, e % Nxyz(1)
         
                        temp = 0._RP
                        

                        do index = 0, e%Nxyz(2)
                            temp = temp &
                                    + spAeta%DT(j1   ,index)  &
                                    * spAeta%DT(j2   ,index)  &
                                    * spAeta%w(index)                &
                                    * DOT_PRODUCT(                     &
                                            e%geom%jGradEta(:,i12,index,k12)/e%geom%Jacobian(i12, index, k12),&
                                            ! e%geom%jGradEta(:,i12,index,k12)/e%geom%Jacobian(i12, index, k12) &
                                            e%geom%jGradEta(:,i12,index,k12) &
                                    )                    
                        end do

                        do index = 1, neqn
                            MatEntries(index,index)=temp
                        end do

                        MatEntries = MatEntries  &
                                                / spAEta%w(j1) / e%geom%Jacobian(i12, j1, k12)  


                        baseRow = i12*neqn + j1*EtaSpa + k12*ZetaSpa
                        baseCol = i12*neqn + j2*EtaSpa + k12*ZetaSpa


                        do eq2 = 1, neqn
                            do eq1 = 1, neqn 
                                ! write (*,*) "bug devide 0--------------------", 1./(eq2-1.0)
                                i = eq1 + baseRow  ! row index (1-based)
                                j = eq2 + baseCol  ! column index (1-based)
                            
                                call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                            
                            end do
                        end do
                    end do
                end do
            end do
        end do

           
!     Zeta contribution (di*dj)
!     -------------------------
        do k2 = 0, e % Nxyz(3)  
            do k1 = 0, e % Nxyz(3)   
                do j12 = 0, e % Nxyz(2)  
                    do i12 = 0, e % Nxyz(1)
         
                        MatEntries = 0._RP

                        temp = 0._RP
                        

                        do index = 0, e%Nxyz(3)
                            temp = temp &
                                    + spAzeta%DT(k1   ,index)  &
                                    * spAzeta%DT(k2   ,index)  &
                                    * spAzeta%w(index)                &
                                    * DOT_PRODUCT(                     &
                                            e%geom%jGradZeta(:,i12,j12,index)/e%geom%Jacobian(i12,j12,index),&
                                            e%geom%jGradZeta(:,i12,j12,index)/e%geom%Jacobian(i12,j12,index) &
                                    )                    
                        end do

                        do index = 1, neqn
                            MatEntries(index,index)=temp
                        end do

                        MatEntries = MatEntries  &
                                                / spAzeta%w(k1)


                        baseRow = i12*neqn + j12*EtaSpa + k1*ZetaSpa
                        baseCol = i12*neqn + j12*EtaSpa + k2*ZetaSpa

                        do eq2 = 1, neqn 
                            do eq1 = 1, neqn 
                                i = eq1 + baseRow  ! row index (1-based)
                                j = eq2 + baseCol  ! column index (1-based)
                            
                                call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                            
                            end do
                        end do
                    end do
                end do
            end do
        end do


!     Xi, Eta contribution (dk)
!     -------------------------
        do k12 = 0, e % Nxyz(3)
            do j2 = 0, e % Nxyz(2)  
                do j1 = 0, e % Nxyz(2)  
                    do i2 = 0, e % Nxyz(1)  
                        do i1 = 0, e % Nxyz(1)  


                            MatEntries = 0._RP
                            temp = 0.0_RP
                                                    
                            ! ===================================
                            ! D(s1,i1)
                            ! D(s2,j2)
                            ! ===================================

                            ! ===================================
                            ! D(s1,i2)
                            ! D(s2,j1)
                            ! ===================================
                            temp = temp + DOT_PRODUCT(           &
                                            e%geom%jGradXi (:,i1,j2,k12)/ e %geom%Jacobian(i1, j2, k12),  &
                                            e%geom%jGradEta(:,i1,j2,k12)/ e %geom%Jacobian(i1, j2, k12)   &
                                        )                                                            & 
                                        * spAxi   %DT(i2,i1) * spAxi%w(i1)                           &
                                        * spAeta  %DT(j1,j2) * spAeta%w(j2)                          &
                                        * spAzeta%w(k12)                                              &
                                        * e%geom%jacobian(i1,j2,k12)                &
                                        / spAxi   %w(i1) / spAeta%w(j1)/spAzeta%w(k12)/e%geom%jacobian(i1,j1,k12)
                            ! write(*,*) "!  Xi, Eta contribution (dk) ================== ", temp
                                         
                                                     
                            temp = temp + DOT_PRODUCT(             &
                                            e%geom%jGradXi (:,i2,j1,k12)/ e %geom%Jacobian(i2, j1, k12),  &
                                            e%geom%jGradEta(:,i2,j1,k12)/ e %geom%Jacobian(i2, j1, k12)   &
                                        )                                                            &
                                        * spAxi   %DT(i1,i2) * spAxi%w(i2)                                &
                                        * spAeta  %DT(j2,j1) * spAeta%w(j1)                           &
                                        * spAzeta%w(k12)                                              &
                                        * e%geom%jacobian(i2,j1,k12)                &
                                        / spAxi   %w(i1) / spAeta%w(j1)/spAzeta%w(k12)/e%geom%jacobian(i1,j1,k12)
                            ! write(*,*) "!  Xi, Eta contribution (dk) ================== ", temp

                            do index = 1, neqn
                                MatEntries(index,index)=temp
                            end do

                          
                            ! write(*,*) "!  Xi, Eta contribution (dk) ================== ", MatEntries
                            baseRow = i1*neqn + j1*EtaSpa + k12*ZetaSpa
                            baseCol = i2*neqn + j2*EtaSpa + k12*ZetaSpa

                            do eq2 = 1, neqn 
                                do eq1 = 1, neqn 
                                    i = eq1 + baseRow  ! row index (1-based)
                                    j = eq2 + baseCol  ! column index (1-based)
                            
                                    call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                            
                                end do
                            end do
                            ! ---------===========*****************************======================
                            ! ---------===========*****************************======================
                        end do
                    end do
                end do
            end do
        end do

! !     Xi, Zeta contribution (dj)
! !     -------------------------
!         do k2 = 0, e % Nxyz(3)
!             do k1 = 0, e % Nxyz(3)
!                 do j12 = 0, e % Nxyz(2)  
!                     do i2 = 0, e % Nxyz(1)  
!                         do i1 = 0, e % Nxyz(1)  
         
!                             ! MatEntries = -16._RP

!                             MatEntries = 0._RP
!                             temp = 0.0_RP
!                             ! ===================================
!                             ! D(s1,i1)
!                             ! D(s3,k2)
!                             ! ===================================

!                             ! ===================================
!                             ! D(s1,i2)
!                             ! D(s3,k1)
!                             ! ===================================

        
!                             temp = temp + DOT_PRODUCT(                                                   &
!                                                         e%geom%jGradXi  (:,i2,j12,k1)/ e %geom%Jacobian(i2, j12, k1), &
!                                                         e%geom%jGradZeta(:,i1,j12,k1)/ e %geom%Jacobian(i1, j12, k1)  &
!                                                     )                                                                   &
!                                                     * spAxi   %DT(i1,i2) * spAxi  %w(i2)                             &
!                                                     * spAzeta %DT(k2,k1)                                     &
!                                                     / spAxi   %w(i1)  
                                                    
!                             temp = temp + DOT_PRODUCT(                                                   &
!                                                         e%geom%jGradXi  (:,i1,j12,k1)/ e %geom%Jacobian(i1, j12, k1), &
!                                                         e%geom%jGradZeta(:,i1,j12,k2)/ e %geom%Jacobian(i1, j12, k2)  &
!                                                     )                                                                   &
!                                                     * spAZeta %DT(k1,k2) * spAZeta%w(k2)                           &
!                                                     * spAxi   %DT(i2,i1)                            &
!                                                     / spAZeta %w(k1) 

                            
!                             do index = 1, neqn
!                                 MatEntries(index,index)=temp
!                             end do
!                             write(*,*) "!   Xi, Zeta contribution (dj) ================== ", MatEntries

!                             baseRow = i1*neqn + j12*EtaSpa + k1*ZetaSpa
!                             baseCol = i2*neqn + j12*EtaSpa + k2*ZetaSpa
                            

!                             do eq2 = 1, neqn 
!                                 do eq1 = 1, neqn 
!                                     i = eq1 + baseRow  ! row index (1-based)
!                                     j = eq2 + baseCol  ! column index (1-based)
                            
!                                     call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                            
!                                 end do
!                             end do
!                             ! ---------===========*****************************======================
!                         end do
!                     end do
!                 end do
!             end do
!         end do


! !     Eta, Zeta contribution (di)
! !     -------------------------
!         do k2 = 0, e % Nxyz(3)
!             do k1 = 0, e % Nxyz(3)
!                 do j2 = 0, e % Nxyz(2)  
!                     do j1 = 0, e % Nxyz(2)  
!                         do i12 = 0, e % Nxyz(1)  
         
!                             MatEntries = 0._RP
!                             temp = 0.0_RP

!                             temp = temp + DOT_PRODUCT(                                         &
!                                                         e%geom%jGradEta(:,i12,j2,k1) / e %geom%Jacobian(i12, j2, k1) ,  &
!                                                         e%geom%jGradZeta(:,i12,j1,k1)/ e %geom%Jacobian(i12, j1, k1)   &
!                                                     )                                                                  &
!                                                     * spAeta  %DT(j1,j2) * spAeta%w(j2)                             &
!                                                     * spAzeta %DT(k2,k1)                           &
!                                                     / spAeta  %w(j1)     
                                                    
!                             temp = temp + DOT_PRODUCT(                                         &
!                                                         e%geom%jGradEta (:,i12,j1,k1) / e %geom%Jacobian(i12, j1, k1) ,  &
!                                                         e%geom%jGradZeta(:,i12,j1,k2) / e %geom%Jacobian(i12, j1, k2)   &
!                                                     )                                                                  &
!                                                     * spAZeta %DT(k1,k2) * spAZeta%w(k2)                           &
!                                                     * spAEta  %DT(j2,j1)                             &
!                                                     / spAZeta %w(k1)   

!                             do index = 1, neqn
!                                 MatEntries(index,index)=temp
!                             end do


!                             write(*,*) "!  Eta, Zeta  contribution (dj) ================== ", MatEntries

!                             baseRow = i12*neqn + j1*EtaSpa + k1*ZetaSpa
!                             baseCol = i12*neqn + j2*EtaSpa + k2*ZetaSpa

!                             do eq2 = 1, neqn 
!                                 do eq1 = 1, neqn 
!                                     i = eq1 + baseRow  ! row index (1-based)
!                                     j = eq2 + baseCol  ! column index (1-based)
                            
!                                     call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                            
!                                 end do
!                             end do
!                             ! ---------===========*****************************======================
!                         end do
!                     end do
!                 end do
!             end do
!         end do


        nullify (spAxi, spAeta, spAzeta)
    end subroutine Local_SetDiagonalBlock


    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  -----------------------------------------------------------------------------------------------
    !  Subroutine for adding the faces' contribution to the off-diagonal blocks of the Jacobian matrix
    !  -> Note: Only the interior faces contribute to the off-diagonal blocks
    !  -> only for p-conforming meshes
    !  -----------------------------------------------------------------------------------------------
    subroutine AnalyticalJacobianPoisson_OffDiagonalBlocks(mesh,nEqn, Matrix, startNum)
        implicit none
        !--------------------------------------------
        type(HexMesh), target    , intent(inout) :: mesh
        class(Matrix_t)          , intent(inout) :: Matrix
        integer, intent(in)    :: nEqn

        integer, optional, intent(in)    :: startNum
        
        !--------------------------------------------
        integer :: eID, fID, elSide, side
        type(Element), pointer :: e_plus
        type(Element), pointer :: e_minus
        type(Face)   , pointer :: f
        !--------------------------------------------

    ! write (*,*) "mesh % no_of_elements,",  mesh % no_of_elements

    !
    !     Compute the off-diagonal blocks for each element's equations
    !     ------------------------------------------------------------
    !$omp do schedule(runtime) private(e_plus,e_minus,elSide,fID,side,f)
        do eID = 1, mesh % no_of_elements
            e_plus => mesh % elements(eID)

    !
    !        One block for every neighbor element
    !        ------------------------------------
            do elSide = 1, 6
                ! write (*,*) "eID, fe_plus % NumberOfConnections(elSide) ======= ", eID, e_plus % NumberOfConnections(elSide)
                if (e_plus % NumberOfConnections(elSide) == 0) then 

                    ! if (e_plus % NumberOfConnections(elSide) == 0) cycle

                    fID  = e_plus % faceIDs(elSide)
                    side = e_plus % faceSide(elSide)

                    ! write (*,*) "fID, side ====BCBCBCBC====**********==**********========= ", fID, side

                    f => mesh % faces(fID)

                    ! write (*,*) "e_plus % Connection(elSide) % globID, e_plus % Connection(elSide) % Nxyz***========= ", e_plus % Connection(elSide) % globID, e_plus % Connection(elSide) % Nxyz

                    if(f%FaceType == HMESH_BOUNDARY) then
                        write(*,*) "startNum in off blocks ============== ", startNum

                        if (startNum == 4) then
                            ! *********************************************************
                            ! ### for the pressure, 
                            ! if velocity dirichlet, 
                            !                   pressure zero gradient
                            !                   Neumann
                            ! if velocity zero gradient
                            !                   pressure has dirichlet BC
                            !                   dirichlet
                            ! *********************************************************
                            if      (TRIM(BCs(f % zone) % bc % BCType) == "zydiyDirichletBC" )  then
                                ! call Local_Get_BC_DiagonalBlock(f, neqn, e_plus,side,Matrix)
                            else if (TRIM(BCs(f % zone) % bc % BCType) == "zydiyfreeslipwall") then
                                call Local_Get_BC_DiagonalBlock(f, neqn, e_plus,side,Matrix)
                            else if (TRIM(BCs(f % zone) % bc % BCType) == "user-defined") then
                                call Local_Get_BC_DiagonalBlock(f, neqn, e_plus,side,Matrix)
                            endif

                        else if (startNum == 5) then
                            ! *********************************************************
                            ! ### for the velocity, 
                            ! if velocity dirichlet, 
                            !                   dirichlet
                            ! if velocity zero gradient
                            !                   Neumann
                            ! *********************************************************

                            if      (TRIM(BCs(f % zone) % bc % BCType) == "zydiyDirichletBC" )  then
                                call Local_Get_BC_DiagonalBlock(f, neqn, e_plus,side,Matrix)

                            else if (TRIM(BCs(f % zone) % bc % BCType) == "zydiyfreeslipwall") then

                            else if (TRIM(BCs(f % zone) % bc % BCType) == "user-defined") then

                                call Local_Get_BC_DiagonalBlock(f, neqn, e_plus,side,Matrix)

                            endif
                                ! call Local_Get_BC_DiagonalBlock(f, neqn, e_plus,side,Matrix)
                        endif

                    endif
            
                else 

                    ! if (e_plus % NumberOfConnections(elSide) == 0) cycle

                    fID  = e_plus % faceIDs(elSide)
                    side = e_plus % faceSide(elSide)

                    ! write (*,*) "fID, side ========**********==**********========= ", fID, side

                    f => mesh % faces(fID)


                    e_minus => mesh % elements(e_plus % Connection(elSide)%globID)
                    call Local_GetOffDiagonalBlock(f, nEqn,e_plus,e_plus % Connection(elSide),e_minus,side,Matrix)
                endif
           
            end do
        end do
    !$omp end do
        nullify (f, e_plus)
     end subroutine AnalyticalJacobianPoisson_OffDiagonalBlocks



     subroutine Local_Get_BC_DiagonalBlock(f, nEqn, e_plus, side, Matrix)
        implicit none
        !-arguments----------------------------------------------------------------------
        type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Element), intent(in)    :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
        integer, intent(in)    :: side    !<  side of face where e_plus is
        class(Matrix_t), intent(inout) :: Matrix  !<> Jacobian matrix
        integer, intent(in)    :: nEqn
        !-local-variables----------------------------------------------------------------
        integer :: i, j                     ! Matrix indexes
        integer :: i1, j1, k1, eq1          ! variable counters
        integer :: i2, j2, k2, eq2          ! variable counters
        integer :: baseRow, baseCol         ! Position of neqn by neqn miniblock of Jacobian
        integer :: r                        ! Additional index for counting in the normal direction (only needed for viscous fluxes)
        integer :: nXi1, nEta1              ! Number of nodes in every direction
        integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
        integer :: nXi2, nEta2              ! Number of nodes in every direction
        integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
        integer :: elSide_plus              ! Element side where f is on e⁻
        integer :: elInd_plus(3)            ! Element indexes on e⁺
        integer :: normAx_plus              ! Normal axis to f on e⁺
        integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
        integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
        integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
        integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
        real(kind=RP) :: MatEntries(neqn, neqn)           ! Values of the matrix entries
        type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
        type(NodalStorage_t), target  :: spA_minus(3)      ! Nodal storage in the different directions for e_minus - local copy
        type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
        type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
        !--------------------------------------------------------------------------------

        ! *********************************************************************************************
        integer :: i12, j12, k12 , index   ,index1          ! i1=i2, j1=j2, k1=k2

        real(kind=RP) :: epsilon, sigma      

        integer :: elInd_plus2(3)            ! Element indexes on e⁺
        integer :: elInd_minus2(3)           ! Element indexes on e⁻
        real(kind=RP) :: JGradXi0(3), JGradXi02(3)  ,jJacobian     ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: MatEntries1(neqn, neqn)           ! Values of the matrix entries
        real(kind=RP) :: MatEntries2(neqn, neqn)           ! Values of the matrix entries
        real(kind=RP) :: normCartesian(3)     ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: temp00, temp01,temp02    ! the variables in the face after interpolation
        integer :: faceInd_plus2(2)         ! Face indexes on e⁺
        integer :: eNorm(3)            ! Element indexes on e⁺
        integer :: Nd1, Nd2, Nd3       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
        real(kind=RP) :: Grad1(3), Grad2(3), Grad3(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face

        ! *********************************************************************************************

        !
!     ***********
!     Definitions
!     ***********
!

        epsilon =-0.0_RP

        sigma =  ( maxval(f % Nf) )**2_RP *10000_RP / (f % geom % h)**1


        ! Entry spacing for element e⁺
        nXi1 = e_plus%Nxyz(1) + 1
        nEta1 = e_plus%Nxyz(2) + 1
        EtaSpa1 = neqn*nXi1
        ZetaSpa1 = neqn*nXi1*nEta1

        ! Entry spacing for element e⁻
        nXi2 = e_plus%Nxyz(1) + 1
        nEta2 = e_plus%Nxyz(2) + 1
        EtaSpa2 = neqn*nXi2
        ZetaSpa2 = neqn*nXi2*nEta2

        ! Element sides
        elSide_plus = f%elementSide(side)


        ! Normal and tangent axes
        normAx_plus = normalAxis(elSide_plus)
        tanAx_plus = axisMap(:, f%elementSide(side))

        ! Side of axis where f is
        if (normAx_plus < 0) then
            normAxSide_plus = LEFT
            normAxSide_minus = RIGHT
        else
            normAxSide_plus = RIGHT
            normAxSide_minus = LEFT
        end if
        normAx_plus = abs(normAx_plus)

     

        ! Nodal storage
        ! --------------------------------------
        ! TODO: Why this doesn't work since ifort ver. 19.1?
        ! --------------------------------------
        ! spA_plus  = NodalStorage(e_plus  % Nxyz)
        ! spA_minus = NodalStorage(e_minus % Nxyz)

        spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
        spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
        spA_plus(3) = NodalStorage(e_plus%Nxyz(3))
 

        spAnorm_plus => spA_plus(normAx_plus)
        spAtan1_plus => spA_plus(tanAx_plus(1))
        spAtan2_plus => spA_plus(tanAx_plus(2))

        
    if (e_plus%Nxyz(normAx_plus) /= 0 ) then

        ! Polynomial orders
        NxyFace_plus = e_plus%Nxyz(tanAx_plus)

        
        if (side == LEFT) then
            normLR =  1.0_RP
        else
            normLR =      1.0_RP
        end if

            ! =====================================M11===========================================================================================
            ! =====================================M11===========================================================================================
            if(normAx_plus == 1) then
                ! ---------------- normal direction of computational space is in the Xi dir 
                Nd3=e_plus%Nxyz(1)
                Nd1=e_plus%Nxyz(2)
                Nd2=e_plus%Nxyz(3)   
            elseIF(normAx_plus == 2) then
                ! ---------------- normal direction of computational space is in the Eta dir 
                Nd3=e_plus%Nxyz(2)
                Nd1=e_plus%Nxyz(1)
                Nd2=e_plus%Nxyz(3)
            else 
                ! ---------------- normal direction of computational space is in the Zeta dir 
                Nd3=e_plus%Nxyz(3)
                Nd1=e_plus%Nxyz(1)
                Nd2=e_plus%Nxyz(2)
            endif

            do j2=0,Nd2
                do j1=0,Nd2
                    do i2=0,Nd1
                        do i1=0,Nd1
                            do k2=0,Nd3
                                do k1=0,Nd3

                            MatEntries = 0
                            MatEntries1 = 0
                            MatEntries2 = 0

                            temp00 = 0
                            temp01 = 0
                            temp02 = 0
                
                        if(i1==i2 .and. j1 == j2) then
                            do index = 0, Nd3
                                if(normAx_plus == 1) then
                                    ! ---------------- normal direction of computational space is in the Xi dir 
                                    Grad3=e_plus%geom%jGradXi(:, index, i2, j2)/e_plus%geom%jacobian(index, i2, j2)
                                    elseIF(normAx_plus == 2) then
                                    ! ---------------- normal direction of computational space is in the Eta dir 
                                    Grad3=e_plus%geom%jGradEta(:, i2, index, j2)/e_plus%geom%jacobian(i2, index, j2)
                                else 
                                    ! ---------------- normal direction of computational space is in the Zeta dir 
                                    Grad3=e_plus%geom%jGradZeta(:, i2, j2, index)/e_plus%geom%jacobian(i2, j2, index)
                                endif

                                temp00  = temp00  + spAtan1_plus%w ( i1)   &
                                                *   spAtan2_plus%w ( j1)   &
                                                *   spAnorm_plus%v ( k1, normAxSide_plus )  &
                                                *   spAnorm_plus%v ( index, normAxSide_plus ) &
                                                *   (   0                                                      &
                                                   +   Grad3(1) *   f%geom%normal(1,i1,j1)   &
                                                   +   Grad3(2) *   f%geom%normal(2,i1,j1)   &              
                                                   +   Grad3(3) *   f%geom%normal(3,i1,j1)   &
                                                    )                             &
                                                *   spAnorm_plus%DT( k2,   index )        &
                                                *   f%geom%jacobian(i1,j1)                                 

                            end do
                        endif

                            ! do index = 0, Nd1
                        if(j1 == j2) then

                                if(normAx_plus == 1) then
                                    ! ---------------- normal direction of computational space is in the Xi dir 
                                    Grad1=e_plus%geom%jGradEta(:, k2, i1, j2)/e_plus%geom%jacobian(k2, i1, j2)
                                    elseIF(normAx_plus == 2) then
                                    ! ---------------- normal direction of computational space is in the Eta dir 
                                    Grad1=e_plus%geom%jGradXi(:, i1, k2, j2)/e_plus%geom%jacobian(i1, k2, j2)
                                else 
                                    ! ---------------- normal direction of computational space is in the Zeta dir 
                                    Grad1=e_plus%geom%jGradXi(:, i1, j2, k2)/e_plus%geom%jacobian(i1, j2, k2)
                                endif

                                temp00  = temp00  + spAtan1_plus%w ( i1)   &
                                                *   spAtan2_plus%w ( j1)   &
                                                *   spAnorm_plus%v ( k1                      , normAxSide_plus )  &
                                                *   spAnorm_plus%v ( k2                      , normAxSide_plus ) &
                                                *   (   0                                                      &
                                                   +    Grad1(1) *   f%geom%normal(1,i1,j1)   &
                                                   +    Grad1(2) *   f%geom%normal(2,i1,j1)   &
                                                   +    Grad1(3) *   f%geom%normal(3,i1,j1)   &
                                                    )                                   &
                                                *   spAtan1_plus%DT( i2,   i1 )  &
                                                *   f%geom%jacobian(i1,j1)                                 
                        endif
                         
                            ! -------------------------------------------------------------
                        if(i1==i2 .and. j1 == j2) then
                            temp02 = temp02 +   spAtan1_plus  %w( i1)  &
                                            * spAtan2_plus  %w( j1)  &
                                            * spAnorm_plus%v ( k1, normAxSide_plus ) &
                                            * spAnorm_plus%v ( k2, normAxSide_plus ) &
                                            * f%geom%jacobian(i1,j1)                                 

                        endif
                            
        
                            if(normAx_plus == 1) then
                                ! ---------------- normal direction of computational space is in the Xi dir 
                                baseRow = k1*neqn + i1*EtaSpa1 + j1*ZetaSpa1
                                baseCol = k2*neqn + i2*EtaSpa1 + j2*ZetaSpa1
                                jJacobian = e_plus%geom%jacobian(k1,i1,j1)
                                elseIF(normAx_plus == 2) then
                                ! ---------------- normal direction of computational space is in the Eta dir 
                                baseRow = i1*neqn + k1*EtaSpa1 + j1*ZetaSpa1
                                baseCol = i2*neqn + k2*EtaSpa1 + j2*ZetaSpa1
                                jJacobian = e_plus%geom%jacobian(i1,k1,j1)
                            else 
                                ! ---------------- normal direction of computational space is in the Zeta dir 
                                baseRow = i1*neqn + j1*EtaSpa1 + k1*ZetaSpa1
                                baseCol = i2*neqn + j2*EtaSpa1 + k2*ZetaSpa1
                                jJacobian = e_plus%geom%jacobian(i1,j1,k1)
                            endif

                            do index = 1,neqn
                                MatEntries(index,index) = 0&
                                    - 2.0_RP * 0.500_RP* temp00 &
                                    + 2.0_RP * 0.500_RP* epsilon * temp01  &
                                    + 2.0_RP * sigma * temp02  
                            end do
                            
                          

                            MatEntries = MatEntries &
                                        /  spAtan1_plus %w(i1) / spAtan2_plus %w(j1) / spAnorm_plus%w(k1) / jJacobian   

                            
                            do eq2 = 1, neqn; 
                                do eq1 = 1, neqn
                                    i = eq1 + baseRow ! row index (1-based)
                                        j = eq2 + baseCol ! column index (1-based)

                                        call Matrix%AddToBlockEntry(e_plus%globID, e_plus%globID, i, j, MatEntries(eq1, eq2))
                                end do; 
                            end do

                        end do
                    end do
                end do
            end do
            end do
            end do

    endif

!     *********
!     Finish up
!     *********
!
    nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus)

end subroutine Local_Get_BC_DiagonalBlock !





!  -----------------------------------------------------------------------------------------------
!  Local_GetOffDiagonalBlock:
!     Routine to compute the off-diagonal block that results of connecting e_plus with e_minus for
!     the equations that correspond to e_plus and inserting it in the Matrix
!
!     -> Currently, this routine only works for p-conforming representations
!           TODO: Add routine for p-nonconforming representations (block computation is more expensive and delta cycling is ruled out)
!  -----------------------------------------------------------------------------------------------
    subroutine Local_GetOffDiagonalBlock(f, nEqn, e_plus, e_minus, e_minus0,side, Matrix)
        implicit none
        !-arguments----------------------------------------------------------------------
        type(Face), target, intent(in)    :: f       !<  Face connecting elements
        integer, intent(in)    :: nEqn
        
        type(Element), intent(in)    :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
        type(Connectivity), intent(in)    :: e_minus !<  Element that connects with e_plus through "f"
        type(Element), intent(in)    :: e_minus0  !<  The off-diagonal block is the contribution to this element's equations
        integer, intent(in)    :: side    !<  side of face where e_plus is
        class(Matrix_t), intent(inout) :: Matrix  !<> Jacobian matrix
        !-local-variables----------------------------------------------------------------
        integer :: i, j                     ! Matrix indexes
        integer :: i1, j1, k1, eq1          ! variable counters
        integer :: i2, j2, k2, eq2          ! variable counters
        integer :: baseRow, baseCol         ! Position of neqn by neqn miniblock of Jacobian
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
        integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
        integer :: tanAx_minus(2)           ! Tangent axes to f on e⁻
        integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
        integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
        integer :: faceInd_plus(2)         ! Face indexes on e⁺
        integer :: faceInd_minus(2)         ! Face indexes on e⁻
        integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
        integer :: NxyFace_minus(2)         ! Polynomial orders of face on element e⁻ (only needed for viscous fluxes)
        integer :: faceInd_plus2minus(2)    ! Face indexes on e⁺ passed to the reference frame of e⁻
        integer :: faceInd_minus2plus(2)    ! Face indexes on e⁻ passed to the reference frame of e⁺ (only needed for viscous fluxes)
        integer :: elInd_1(3)       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
        integer :: elInd_2(3)       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
        integer :: elInd_plus_tan1(3)       ! Element indexes on e⁺, but the index in the first tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
        real(kind=RP) :: MatEntries(neqn, neqn)           ! Values of the matrix entries
        type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
        type(NodalStorage_t), target  :: spA_minus(3)      ! Nodal storage in the different directions for e_minus - local copy
        type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
        type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAnorm_minus     ! Nodal storage in the direction that is normal to the face for e⁻
        type(NodalStorage_t), pointer :: spAtan1_minus     ! Nodal storage in the tangent direction "1" to the face for e⁻ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_minus     ! Nodal storage in the tangent direction "2" to the face for e⁻ (only needed for viscous fluxes)
        !--------------------------------------------------------------------------------

        ! *********************************************************************************************
        ! integer :: elInd_minus(3)            ! Element indexes on e⁺
        ! integer :: elInd_minus2(3)            ! Element indexes on e⁻
        integer :: i12, j12, k12 ,index, index1             ! i1=i2, j1=j2, k1=k2
        integer :: i12other, j12other, k12other, i1other, i2other, j1other, j2other

        real(kind=RP) :: epsilon, sigma      

        integer :: elInd_plus2(3)            ! Element indexes on e⁺
        integer :: elInd_minus2(3)           ! Element indexes on e⁻
        real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: MatEntries1(neqn, neqn)           ! Values of the matrix entries
        real(kind=RP) :: MatEntries2(neqn, neqn)           ! Values of the matrix entries
        real(kind=RP) :: normCartesian(3)     ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: normLR , tmp01   ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        integer :: faceInd_plus2(2)         ! Face indexes on e⁺
        integer :: faceInd_minus2(2)         ! Face indexes on e⁻
        integer :: faceInd_plus2minus2(2)    ! Face indexes on e⁺ passed to the reference frame of e⁻
        integer :: faceInd_minus2plus2(2)         ! Face indexes on e⁻
        integer :: eNorm(3)       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
        integer :: Nd1, Nd2, Nd3       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
        
        real(kind=RP) :: Grad1(3), Grad2(3), Grad3(3), jJacobian      ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face

        real(kind=RP) :: temp00, temp01,temp02    ! the variables in the face after interpolation
        ! e_minus0 = 
        ! *********************************************************************************************

        !
!     ***********
!     Definitions
!     ***********
!

        epsilon =-0.0_RP

        sigma =  ( maxval(f % Nf) )**2_RP *10000_RP / (f % geom % h)**1


      
        ! Entry spacing for element e⁺
        nXi1 = e_plus%Nxyz(1) + 1
        nEta1 = e_plus%Nxyz(2) + 1
        EtaSpa1 = neqn*nXi1
        ZetaSpa1 = neqn*nXi1*nEta1
        ! write (*,*) "nXi1, nEta1, EtaSpa1, ZetaSpa1 === ", nXi1, nEta1, EtaSpa1, ZetaSpa1

        ! Entry spacing for element e⁻
        nXi2 = e_minus%Nxyz(1) + 1
        nEta2 = e_minus%Nxyz(2) + 1
        EtaSpa2 = neqn*nXi2
        ZetaSpa2 = neqn*nXi2*nEta2
        ! write (*,*) "nXi2, nEta2, EtaSpa2, ZetaSpa2 === ", nXi2, nEta2, EtaSpa2, ZetaSpa2

        ! Element sides
        elSide_plus = f%elementSide(side)
        elSide_minus = f%elementSide(other(side))

        ! write (*,*) "elSide_plus, elSide_minus === ", elSide_plus, elSide_minus


        ! Normal and tangent axes
        normAx_plus = normalAxis(elSide_plus)
        normAx_minus = normalAxis(elSide_minus)
        tanAx_plus = axisMap(:, f%elementSide(side))
        tanAx_minus = axisMap(:, f%elementSide(other(side)))

        ! write (*,*) "normAx_plus, normAx_minus, tanAx_plus, tanAx_minus=== ", normAx_plus, normAx_minus, tanAx_plus, tanAx_minus

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
        ! --------------------------------------
        ! TODO: Why this doesn't work since ifort ver. 19.1?
        ! --------------------------------------
        ! spA_plus  = NodalStorage(e_plus  % Nxyz)
        ! spA_minus = NodalStorage(e_minus % Nxyz)

        spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
        spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
        spA_plus(3) = NodalStorage(e_plus%Nxyz(3))
        spA_minus(1) = NodalStorage(e_minus%Nxyz(1))
        spA_minus(2) = NodalStorage(e_minus%Nxyz(2))
        spA_minus(3) = NodalStorage(e_minus%Nxyz(3))

        spAnorm_plus => spA_plus(normAx_plus)
        spAtan1_plus => spA_plus(tanAx_plus(1))
        spAtan2_plus => spA_plus(tanAx_plus(2))

        spAnorm_minus => spA_minus(normAx_minus)
        spAtan1_minus => spA_minus(tanAx_minus(1))
        spAtan2_minus => spA_minus(tanAx_minus(2))


        ! Polynomial orders
        NxyFace_plus = e_plus%Nxyz(tanAx_plus)
        NxyFace_minus = e_minus%Nxyz(tanAx_minus)

        
!     ********************
!     Viscous contribution
!     ********************
!
        if (flowIsScalar_ins_v04) then
            if (e_plus%Nxyz(normAx_plus) /= 0 ) then


            if (side == LEFT) then
                normLR =  1.0_RP
            else
                normLR =      1.0_RP
            end if

            ! =====================================M11===========================================================================================
            ! =====================================M11===========================================================================================
            if(normAx_plus == 1) then
                ! ---------------- normal direction of computational space is in the Xi dir 
                Nd3=e_plus%Nxyz(1)
                Nd1=e_plus%Nxyz(2)
                Nd2=e_plus%Nxyz(3) 
                ! write (*,*) "if(normAx_plus == 1) then", normAx_plus

            elseIF(normAx_plus == 2) then
                ! ---------------- normal direction of computational space is in the Eta dir 
                Nd3=e_plus%Nxyz(2)
                Nd1=e_plus%Nxyz(1)
                Nd2=e_plus%Nxyz(3)
                ! write (*,*) "if(normAx_plus == 2) then", normAx_plus
            else 
                ! ---------------- normal direction of computational space is in the Zeta dir 
                Nd3=e_plus%Nxyz(3)
                Nd1=e_plus%Nxyz(1)
                Nd2=e_plus%Nxyz(2)
                ! write (*,*) "if(normAx_plus == 3) then", normAx_plus
            endif

          
                    do k2 = 0, e_plus%Nxyz(3); 
                        do k1 = 0, e_plus%Nxyz(3); 
                            do j2 = 0, e_plus%Nxyz(2); 
                                do j1 = 0, e_plus%Nxyz(2); 
                                do i2 = 0, e_plus%Nxyz(1)
                                    do i1 = 0, e_plus%Nxyz(1)

                                    elInd_plus    = [i1, j1, k1]
                                    elInd_plus2   = [i2, j2, k2]
                                    elInd_minus   = [i1, j1, k1]
                                    elInd_minus2  = [i2, j2, k2]

                                    faceInd_plus     = elInd_plus  ( tanAx_plus )
                                    faceInd_plus2    = elInd_plus2 ( tanAx_plus )
                                    faceInd_minus    = elInd_minus ( tanAx_minus )
                                    faceInd_minus2   = elInd_minus2( tanAx_minus )

                                    elInd_plus_tan1 = elInd_plus
                                    elInd_plus_tan1(normAx_plus) = elInd_minus2(normAx_minus)

                                    call indexesOnOtherFace( faceInd_plus  (1),faceInd_plus  (2), NxyFace_plus (1), NxyFace_plus (2), f % rotation, side       , faceInd_plus2minus (1),faceInd_plus2minus (2) )
                                    call indexesOnOtherFace( faceInd_plus2 (1),faceInd_plus2 (2), NxyFace_plus (1), NxyFace_plus (2), f % rotation, side       , faceInd_plus2minus2(1),faceInd_plus2minus2(2) )
                                    call indexesOnOtherFace( faceInd_minus (1),faceInd_minus (2), NxyFace_minus(1), NxyFace_minus(2), f % rotation, other(side), faceInd_minus2plus (1),faceInd_minus2plus (2) )
                                    call indexesOnOtherFace( faceInd_minus2(1),faceInd_minus2(2), NxyFace_minus(1), NxyFace_minus(2), f % rotation, other(side), faceInd_minus2plus2(1),faceInd_minus2plus2(2) )
                                    

                            MatEntries = 0
                            MatEntries1 = 0
                            MatEntries2 = 0

                            temp00 = 0
                            temp01 = 0
                            temp02 = 0

                        if (( elInd_plus(tanAx_plus(1)) == elInd_plus2(tanAx_plus(1)) ) .and. &
                            ( elInd_plus(tanAx_plus(2)) == elInd_plus2(tanAx_plus(2)) ) ) then
                           
                           
                            do index = 0, e_plus%Nxyz(normAx_plus)

                                if(normAx_plus == 1) then 
                                    Grad3 = e_plus%geom%jGradXi  (:,index, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))/e_plus%geom%jacobian(index, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))
                                elseIF(normAx_plus == 2) then
                                    Grad3 = e_plus%geom%jGradEta  (:,elInd_plus(tanAx_plus(1)), index, elInd_plus(tanAx_plus(2)))/e_plus%geom%jacobian(elInd_plus(tanAx_plus(1)), index, elInd_plus(tanAx_plus(2)))
                                else 
                                    Grad3 = e_plus%geom%jGradZeta  (:,elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)), index)/e_plus%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)), index)
                                endif

                                temp00  = temp00  + spAtan1_plus%w ( elInd_plus(tanAx_plus(1)) )   &
                                                *   spAtan2_plus%w ( elInd_plus(tanAx_plus(2)) )   &
                                                *   spAnorm_plus%v ( elInd_plus (normAx_plus) , normAxSide_plus )  &
                                                *   spAnorm_plus%v ( index                    , normAxSide_plus ) &
                                                *   (   0                                                      &
                                                   +  Grad3(1) *   f%geom%normal(1, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))  ) &
                                                   +  Grad3(2) *   f%geom%normal(2, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))  ) &
                                                   +  Grad3(3) *   f%geom%normal(3, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))  ) &
                                                    )                                         &  
                                                *   spAnorm_plus%DT( elInd_plus2(normAx_plus),   index )                      &
                                                *   f%geom%jacobian( elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )                                 

                            enddo
                                    
                            temp02 = temp02 +   spAtan1_plus  %w( elInd_plus(tanAx_plus(1)) )  &
                                            *   spAtan2_plus  %w( elInd_plus(tanAx_plus(2)) )  &
                                            *   spAnorm_plus%v ( elInd_plus (normAx_plus) , normAxSide_plus ) &
                                            *   spAnorm_plus%v ( elInd_plus2(normAx_plus) , normAxSide_plus ) &
                                            *   f%geom%jacobian( elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )                                 
                                    
                        endif

                        if( elInd_plus(tanAx_plus(2)) == elInd_plus2(tanAx_plus(2)) ) then

                            if(normAx_plus == 1) then 
                                Grad1 = e_plus%geom%jGradEta  (:,elInd_plus2(normAx_plus), elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))/e_plus%geom%jacobian(elInd_plus2(normAx_plus), elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))
                            elseIF(normAx_plus == 2) then
                                Grad1 = e_plus%geom%jGradXi  (:,elInd_plus(tanAx_plus(1)), elInd_plus2(normAx_plus), elInd_plus(tanAx_plus(2)))/e_plus%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus2(normAx_plus), elInd_plus(tanAx_plus(2)))
                            else 
                                Grad1 = e_plus%geom%jGradXi  (:,elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)), elInd_plus2(normAx_plus))/e_plus%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)), elInd_plus2(normAx_plus))
                            endif

                            temp00  = temp00  + spAtan1_plus%w ( elInd_plus(tanAx_plus(1)) )   &
                                            *   spAtan2_plus%w ( elInd_plus(tanAx_plus(2)) )   &
                                            *   spAnorm_plus%v ( elInd_plus (normAx_plus) , normAxSide_plus )  &
                                            *   spAnorm_plus%v ( elInd_plus2(normAx_plus) , normAxSide_plus ) &
                                            *   (   0                                                      &
                                               +    Grad1(1) *  f%geom%normal(1, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))  )  &
                                               +    Grad1(2) *  f%geom%normal(2, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))  )  &
                                               +    Grad1(3) *  f%geom%normal(3, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))  )  &
                                                )                                   &
                                            *   spAtan1_plus%DT(  elInd_plus2(tanAx_plus(1)) ,    elInd_plus(tanAx_plus(1) ) )  &
                                            *   f%geom%jacobian( elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )                                 

                        endif

                    
                        if (side == LEFT) then
                            do index = 1, neqn
                                MatEntries(index,index) = 0&
                                                        - 0.2500_RP* temp00           &
                                                        + 0.2500_RP* epsilon * temp01  &
                                                        + 0.5_RP*sigma * temp02
                            enddo
                        
                        else
                            do index = 1, neqn
                                MatEntries(index,index) = 0&
                                                        + 0.2500_RP* temp00           &
                                                        - 0.2500_RP* epsilon * temp01  &
                                                        + 0.5_RP*sigma * temp02
                            enddo
                        
                        endif
                        
                        MatEntries = MatEntries &
                                        /   spAtan1_plus %w( elInd_plus(tanAx_plus(1)) ) &
                                        /   spAtan2_plus %w( elInd_plus(tanAx_plus(2)) ) &
                                        /   spAnorm_plus %w( elInd_plus(normAx_plus)  ) &
                                        /   e_plus%geom%jacobian(i1, j1, k1)   
                            
                    
                        baseRow = i1*neqn + j1*EtaSpa1 + k1*ZetaSpa1
                        baseCol = i2*neqn + j2*EtaSpa2 + k2*ZetaSpa2
                    
                        do eq2 = 1, neqn; 
                            do eq1 = 1, neqn
                                i = eq1 + baseRow ! row index (1-based)
                                j = eq2 + baseCol ! column index (1-based)
                                call Matrix%AddToBlockEntry(e_plus%globID, e_plus%globID, i, j, MatEntries(eq1, eq2))
                            end do; 
                        end do
                    

            ! =====================================M22===========================================================================================
            ! =====================================M22===========================================================================================

                            MatEntries = 0
                            MatEntries1 = 0
                            MatEntries2 = 0

                            temp00 = 0
                            temp01 = 0
                            temp02 = 0

                            elInd_1 = elInd_minus
                            elInd_1( tanAx_plus(1) ) =  faceInd_minus2plus(1)
                            elInd_1( tanAx_plus(2) ) =  faceInd_minus2plus(2)
                                
                            elInd_2 = elInd_minus2
                            elInd_2( tanAx_plus(1) ) =  faceInd_minus2plus2(1)
                            elInd_2( tanAx_plus(2) ) =  faceInd_minus2plus2(2)
                            
                        if ( ( elInd_minus(tanAx_minus(1)) == elInd_minus2(tanAx_minus(1)) )  .and. &
                             ( elInd_minus(tanAx_minus(2)) == elInd_minus2(tanAx_minus(2)) )  ) then



                            do index = 0, e_minus0%Nxyz(normAx_minus)
                                if(normAx_minus == 1) then
                                    ! ---------------- normal direction of computational space is in the Xi dir 
                                    Grad3=e_minus0%geom%jGradXi  (:,index     , elInd_minus(2), elInd_minus(3) ) / e_minus0%geom%jacobian(index     , elInd_minus(2), elInd_minus(3) )
                                elseIF(normAx_minus == 2) then
                                    ! ---------------- normal direction of computational space is in the Eta dir 
                                    Grad3=e_minus0%geom%jGradEta (:,elInd_minus(1), index     , elInd_minus(3) ) / e_minus0%geom%jacobian(elInd_minus(1), index     , elInd_minus(3) )
                                else 
                                    ! ---------------- normal direction of computational space is in the Zeta dir 
                                    Grad3=e_minus0%geom%jGradZeta(:,elInd_minus(1), elInd_minus(2), index      ) / e_minus0%geom%jacobian(elInd_minus(1), elInd_minus(2), index      )
                                endif


                                temp00  = temp00  + spAtan1_minus %w ( elInd_minus( tanAx_minus(1) ) )   &
                                                *   spAtan2_minus %w ( elInd_minus( tanAx_minus(2) ) )   &
                                                *   spAnorm_minus %v ( elInd_minus( normAx_minus  )  , normAxSide_minus )  &
                                                *   spAnorm_minus %v ( index                         , normAxSide_minus ) &
                                                *   (   0                                                      &
                                                   +    Grad3(1) *   f%geom%normal(1, elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )   &
                                                   +    Grad3(2) *   f%geom%normal(2, elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )   &
                                                   +    Grad3(3) *   f%geom%normal(3, elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )   &
                                                    )                                           &
                                                *   spAnorm_minus%DT( elInd_minus2( normAx_minus  ),   index )             &
                                                *   f%geom%jacobian( elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )

                            enddo

                            temp02 = temp02 +   spAtan1_minus %w ( elInd_minus( tanAx_minus(1) ) )  &
                                            *   spAtan2_minus %w ( elInd_minus( tanAx_minus(2) ) )  &
                                            *   spAnorm_minus%v ( elInd_minus( normAx_minus  ) , normAxSide_minus ) &
                                            *   spAnorm_minus%v ( elInd_minus2( normAx_minus  ) , normAxSide_minus ) &
                                            *   f%geom%jacobian( elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )
                                            
                                            
                        endif

                        if(( elInd_minus(tanAx_minus(2)) == elInd_minus2(tanAx_minus(2)) ) ) then
                            if(normAx_minus == 1) then
                                ! ---------------- normal direction of computational space is in the Xi dir 
                                Grad1 = e_minus0%geom%jGradEta(:,elInd_minus2( normAx_minus  ) , elInd_minus(2), elInd_minus(3) ) / e_minus0%geom%jacobian(elInd_minus2( normAx_minus  ) , elInd_minus(2), elInd_minus(3) )
                            elseIF(normAx_minus == 2) then
                                ! ---------------- normal direction of computational space is in the Eta dir 
                                Grad1 = e_minus0%geom%jGradXi (:,elInd_minus(1), elInd_minus2( normAx_minus  ) , elInd_minus(3) ) / e_minus0%geom%jacobian(elInd_minus(1), elInd_minus2( normAx_minus  ) , elInd_minus(3) )
                            else 
                                ! ---------------- normal direction of computational space is in the Zeta dir 
                                Grad1 = e_minus0%geom%jGradXi(:,elInd_minus(1), elInd_minus(2), elInd_minus2( normAx_minus  )  ) / e_minus0%geom%jacobian(elInd_minus(1), elInd_minus(2), elInd_minus2( normAx_minus  )  )
                            endif



                            temp00  = temp00  + spAtan1_minus%w ( elInd_minus ( tanAx_minus(1) ) )   &
                                            *   spAtan2_minus%w ( elInd_minus ( tanAx_minus(2) ) )   &
                                            *   spAnorm_minus%v ( elInd_minus ( normAx_minus )   , normAxSide_minus )  &
                                            *   spAnorm_minus%v ( elInd_minus2( normAx_minus )   , normAxSide_minus ) &
                                            *   (   0                                                      &
                                               +  Grad1(1) *   f%geom%normal(1, elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )  &
                                               +  Grad1(2) *   f%geom%normal(2, elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )  &
                                               +  Grad1(3) *   f%geom%normal(3, elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )  &
                                                )                                    &
                                            *   spAtan1_minus%DT( elInd_minus2(tanAx_minus(1)), elInd_minus( tanAx_minus(1) ) )   &
                                            *   f%geom%jacobian( elInd_1( tanAx_plus(1) ), elInd_1( tanAx_plus(2) ) )
                        endif

                      
                        if (side == LEFT) then
                            do index = 1,neqn
                                MatEntries(index,index) = 0&
                                                        + 0.2500_RP* temp00           &
                                                        - 0.2500_RP* epsilon * temp01  &
                                                        + 0.5_RP*sigma * temp02
                            enddo
                        else
                            do index = 1,neqn
                                MatEntries(index,index) = 0&
                                                        - 0.2500_RP* temp00           &
                                                        + 0.2500_RP* epsilon * temp01  &
                                                        + 0.5_RP* sigma * temp02
                            enddo
                        endif

                        MatEntries = MatEntries &
                                    /   spAtan1_minus %w( elInd_minus( tanAx_minus(1) ) ) &
                                    /   spAtan2_minus %w( elInd_minus( tanAx_minus(2) ) ) &
                                    /   spAnorm_minus %w( elInd_minus( normAx_minus   ) ) &
                                    ! /   e_minus0 % geom%jacobian(elInd_minus(1), elInd_minus(2), elInd_minus(3))   
                                    /   e_minus0 % geom%jacobian(elInd_minus(1), elInd_minus(2), elInd_minus(3))   

                        baseRow = elInd_minus (1)*neqn + elInd_minus (2)*EtaSpa1 + elInd_minus (3)*ZetaSpa1
                        baseCol = elInd_minus2 (1)*neqn + elInd_minus2 (2)*EtaSpa2 + elInd_minus2 (3)*ZetaSpa2
            
                        do eq2 = 1, neqn; 
                            do eq1 = 1, neqn
                                i = eq1 + baseRow ! row index (1-based)
                                j = eq2 + baseCol ! column index (1-based)
                            
                                call Matrix%AddToBlockEntry(e_minus0%globID, e_minus0%globID, i, j, MatEntries(eq1, eq2))
                            end do; 
                        end do
                  

            ! =====================================M12===========================================================================================
            ! =====================================M12===========================================================================================

                            MatEntries = 0
                            MatEntries1 = 0
                            MatEntries2 = 0

                            temp00 = 0
                            temp01 = 0
                            temp02 = 0

                            elInd_1 = elInd_plus
                            elInd_1( tanAx_plus(1) ) =  faceInd_minus2plus2(1)
                            elInd_1( tanAx_plus(2) ) =  faceInd_minus2plus2(2)
                                    
                            elInd_2 =elInd_minus2
                            elInd_2( tanAx_minus(1) ) =  faceInd_plus2minus(1)
                            elInd_2( tanAx_minus(2) ) =  faceInd_plus2minus(2)


                        if (( elInd_plus(tanAx_plus(1)) == elInd_1( tanAx_plus(1) ) ) .and. &
                            ( elInd_plus(tanAx_plus(2)) == elInd_1( tanAx_plus(2) ) ) ) then

                            do index = 0, e_minus0%Nxyz(normAx_minus)
                                if(normAx_minus == 1) then
                                    Grad3=e_minus0%geom%jGradXi  (:,index     , elInd_minus2(2), elInd_minus2(3) ) / e_minus0%geom%jacobian(index     , elInd_minus2(2), elInd_minus2(3) )
                                elseIF(normAx_minus == 2) then
                                    Grad3=e_minus0%geom%jGradEta (:,elInd_minus2(1), index     , elInd_minus2(3) ) / e_minus0%geom%jacobian(elInd_minus2(1), index     , elInd_minus2(3) )
                                else 
                                    Grad3=e_minus0%geom%jGradZeta(:,elInd_minus2(1), elInd_minus2(2), index      ) / e_minus0%geom%jacobian(elInd_minus2(1), elInd_minus2(2), index      )
                                endif

                                temp00  = temp00  + spAtan1_plus %w ( elInd_plus(tanAx_plus(1)) )   &
                                                *   spAtan2_plus %w ( elInd_plus(tanAx_plus(2)) )   &
                                                *   spAnorm_plus  %v ( elInd_plus( normAx_plus  ) , normAxSide_plus )  &
                                                *   spAnorm_minus %v ( index                      , normAxSide_minus ) &
                                                *   (   0                                                      &
                                                   +    Grad3(1) *   f%geom%normal(1, elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )   &
                                                   +    Grad3(2) *   f%geom%normal(2, elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )   &
                                                   +    Grad3(3) *   f%geom%normal(3, elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )   &
                                                    )                                           &
                                                *   spAnorm_minus%DT( elInd_minus2( normAx_minus  ),   index )             &
                                                *   f%geom%jacobian( elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )    

                            enddo
                            temp02 = temp02 +   spAtan1_plus %w ( elInd_plus(tanAx_plus(1)) )    &
                                            *   spAtan2_plus %w ( elInd_plus(tanAx_plus(2)) )    &
                                            *   spAnorm_plus %v ( elInd_plus  ( normAx_plus  )      , normAxSide_plus ) &
                                            *   spAnorm_minus%v ( elInd_minus2( normAx_minus )      , normAxSide_minus ) &
                                            *   f%geom%jacobian( elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )    
                            

                        endif

                        if( ( elInd_plus(tanAx_plus(2)) == elInd_1( tanAx_plus(2) ) ) ) then
                            if(normAx_minus == 1) then
                                Grad1=e_minus0%geom%jGradEta (:,elInd_minus2(normAx_minus), elInd_2(tanAx_minus(1)),      elInd_2(tanAx_minus(2)))/e_minus0%geom%jacobian(elInd_minus2( normAx_minus ), elInd_2(tanAx_minus(1)), elInd_2(tanAx_minus(2)))
                            elseIF(normAx_minus == 2) then
                                Grad1=e_minus0%geom%jGradXi  (:, elInd_2(tanAx_minus(1)),   elInd_minus2( normAx_minus ), elInd_2(tanAx_minus(2)))/e_minus0%geom%jacobian(elInd_2(tanAx_minus(1)), elInd_minus2( normAx_minus ), elInd_2(tanAx_minus(2)))
                            else 
                                Grad1=e_minus0%geom%jGradXi  (:, elInd_2(tanAx_minus(1)),   elInd_2(tanAx_minus(2)),      elInd_minus2( normAx_minus ))/e_minus0%geom%jacobian(elInd_2(tanAx_minus(1)), elInd_2(tanAx_minus(2)), elInd_minus2( normAx_minus ))
                            endif


                            temp00  = temp00  + spAtan1_plus %w ( elInd_plus(tanAx_plus(1)) )   &
                                            *   spAtan2_plus %w ( elInd_plus(tanAx_plus(2)) )   &
                                            *   spAnorm_plus %v ( elInd_plus  ( normAx_plus  ), normAxSide_plus  )  &
                                            *   spAnorm_minus%v ( elInd_minus2( normAx_minus ), normAxSide_minus ) &
                                            *   (   0                                                      &
                                               +  Grad1(1) *   f%geom%normal(1, elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )   &
                                               +  Grad1(2) *   f%geom%normal(2, elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )   &
                                               +  Grad1(3) *   f%geom%normal(3, elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )   &
                                                )                                          &
                                            *   spAtan1_minus%DT( elInd_minus2(tanAx_minus(1)),   elInd_2( tanAx_minus(1) ) )      &
                                            *   f%geom%jacobian( elInd_plus(tanAx_plus(1)) , elInd_plus(tanAx_plus(2)) )    

                    endif
                    

                        if (side == LEFT) then
                            do index = 1,neqn
                                MatEntries(index,index) = 0&
                                            - 0.2500_RP* temp00           &
                                            - 0.2500_RP* epsilon * temp01  &
                                            - 0.5_RP*sigma * temp02
                            enddo
                        else
                            do index = 1,neqn
                                MatEntries(index,index) = 0&
                                        + 0.2500_RP* temp00           &
                                        + 0.2500_RP* epsilon * temp01  &
                                        - 0.5_RP*sigma * temp02
                            enddo

                        endif

                        MatEntries = MatEntries &
                                        / spAtan1_plus %w( elInd_plus( tanAx_plus(1) ) ) &
                                        / spAtan2_plus %w( elInd_plus( tanAx_plus(2) ) ) &
                                        / spAnorm_plus %w( elInd_plus( normAx_plus   ) ) &
                                        / e_plus%geom%jacobian(elInd_plus(1) , elInd_plus(2) , elInd_plus(3) )   
                                        
                        baseRow = elInd_plus(1) *neqn + elInd_plus(2) *EtaSpa1 + elInd_plus(3)*ZetaSpa1
                        baseCol = elInd_minus2(1)*neqn + elInd_minus2(2)*EtaSpa2 + elInd_minus2(3)*ZetaSpa2

                            
                        do eq2 = 1, neqn; 
                            do eq1 = 1, neqn
                                i = eq1 + baseRow ! row index (1-based)
                                j = eq2 + baseCol ! column index (1-based)
                                call Matrix%AddToBlockEntry(e_plus%globID, e_minus0%globID, i, j, MatEntries(eq1, eq2))
                            end do; 
                        end do


            ! =====================================M21===========================================================================================
            ! =====================================M21===========================================================================================

                            MatEntries = 0
                            MatEntries1 = 0
                            MatEntries2 = 0

                            temp00 = 0
                            temp01 = 0
                            temp02 = 0

                            elInd_2 = elInd_plus2
                            elInd_2( tanAx_plus(1) ) =  faceInd_minus2plus(1)
                            elInd_2( tanAx_plus(2) ) =  faceInd_minus2plus(2)
                           
                        if (( elInd_plus2( tanAx_plus(1) )  == elInd_2( tanAx_plus(1) ) ) .and. &
                            ( elInd_plus2( tanAx_plus(2) )  == elInd_2( tanAx_plus(2) ) ) ) then

                                    
                            do index = 0, e_plus%Nxyz(normAx_plus)

                                if(normAx_plus == 1) then 
                                    Grad3 = e_plus%geom%jGradXi  (:,index, elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2)) )/e_plus%geom%jacobian(index, elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2)))
                                elseIF(normAx_plus == 2) then
                                    Grad3 = e_plus%geom%jGradEta (:,elInd_plus2(tanAx_plus(1)), index, elInd_plus2(tanAx_plus(2)))/e_plus%geom%jacobian(elInd_plus2(tanAx_plus(1)), index, elInd_plus2(tanAx_plus(2)))
                                else 
                                    Grad3 = e_plus%geom%jGradZeta(:,elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2)), index)/e_plus%geom%jacobian(elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2)), index)
                                endif

                                temp00  = temp00  + spAtan1_plus %w ( elInd_plus2   (tanAx_plus(1)) )   &
                                                *   spAtan2_plus %w ( elInd_plus2   (tanAx_plus(2)) )   &
                                                *   spAnorm_minus %v ( elInd_minus   (normAx_minus) , normAxSide_minus )  &
                                                *   spAnorm_plus  %v ( index                    , normAxSide_plus ) &
                                                *   (   0                                                      &
                                                   +  Grad3(1) *   f%geom%normal(1, elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2))  ) &
                                                   +  Grad3(2) *   f%geom%normal(2, elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2))  ) &
                                                   +  Grad3(3) *   f%geom%normal(3, elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2))  ) &
                                                    )                                         &  
                                                *   spAnorm_plus%DT( elInd_plus2(normAx_plus),   index )                      &
                                                *   f%geom%jacobian( elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2)) )                                 


                            enddo

                            temp02 = temp02 +   spAtan1_plus %w ( elInd_plus2   (tanAx_plus(1)) )   &
                                            *   spAtan2_plus %w ( elInd_plus2   (tanAx_plus(2)) )   &
                                            *   spAnorm_minus %v ( elInd_minus   (normAx_minus) , normAxSide_minus) &
                                            *   spAnorm_plus  %v ( elInd_plus2(normAx_plus) , normAxSide_plus ) &
                                            *   f%geom%jacobian( elInd_plus2(tanAx_plus(1)), elInd_plus2(tanAx_plus(2)) )                                 

                        endif


                        if( ( elInd_plus2( tanAx_plus(2) )  == elInd_2( tanAx_plus(2) ) )) then

                            if(normAx_plus == 1) then 
                                Grad1 = e_plus%geom%jGradEta  (:,elInd_plus2(normAx_plus), elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2)) )/e_plus%geom%jacobian(elInd_plus2(normAx_plus), elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2)))
                            elseIF(normAx_plus == 2) then
                                Grad1 = e_plus%geom%jGradXi (:,elInd_2(tanAx_plus(1)), elInd_plus2(normAx_plus), elInd_2(tanAx_plus(2)))/e_plus%geom%jacobian(elInd_2(tanAx_plus(1)), elInd_plus2(normAx_plus), elInd_2(tanAx_plus(2)))
                            else 
                                Grad1 = e_plus%geom%jGradXi(:,elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2)), elInd_plus2(normAx_plus))/e_plus%geom%jacobian(elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2)), elInd_plus2(normAx_plus))
                            endif


                            temp00  = temp00  + spAtan1_plus %w ( elInd_2   (tanAx_plus(1)) )   &
                                            *   spAtan2_plus %w ( elInd_2   (tanAx_plus(2)) )   &
                                            *   spAnorm_minus %v ( elInd_minus   (normAx_minus) , normAxSide_minus )  &
                                            *   spAnorm_plus  %v ( elInd_plus2   (normAx_plus) , normAxSide_plus ) &
                                            *   (   0                                                      &
                                               + Grad1(1) *    f%geom%normal(1, elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2))  ) &
                                               + Grad1(2) *    f%geom%normal(2, elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2))  ) &
                                               + Grad1(3) *    f%geom%normal(3, elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2))  ) &
                                                )                                            &
                                            *   spAnorm_plus%DT( elInd_plus2   (tanAx_plus(1)) ,   elInd_2(tanAx_plus(1)) )                     &
                                            *   f%geom%jacobian( elInd_2(tanAx_plus(1)), elInd_2(tanAx_plus(2)))                                 

                        endif


                        if (side == LEFT) then
                            do index=1,neqn
                                MatEntries(index,index) = 0&
                                        + 0.2500_RP* temp00           &
                                        + 0.2500_RP* epsilon * temp01 &
                                        - 0.5_RP*sigma * temp02
                            enddo
                        else
                            do index=1,neqn
                                MatEntries(index,index) = 0&
                                        - 0.2500_RP* temp00           &
                                        - 0.2500_RP* epsilon * temp01  &
                                        - 0.5_RP*sigma * temp02
                            enddo
                        endif


                        MatEntries = MatEntries &
                                    /   spAtan1_minus %w( elInd_minus( tanAx_minus(1) ) ) &
                                    /   spAtan2_minus %w( elInd_minus( tanAx_minus(2) ) ) &
                                    /   spAnorm_minus  %w( elInd_minus( normAx_minus   ) ) &
                                    /   e_minus0 % geom%jacobian(elInd_minus(1), elInd_minus(2), elInd_minus(3))   

                        baseRow = elInd_minus(1)*neqn + elInd_minus (2)*EtaSpa1 + elInd_minus (3)*ZetaSpa1
                        baseCol = elInd_plus2(1) *neqn + elInd_plus2(2) *EtaSpa2 + elInd_plus2(3) *ZetaSpa2
                        

                        do eq2 = 1, neqn; 
                            do eq1 = 1, neqn
                                i = eq1 + baseRow ! row index (1-based)
                                j = eq2 + baseCol ! column index (1-based)
                                call Matrix%AddToBlockEntry(e_minus0%globID, e_plus%globID, i, j, MatEntries(eq1, eq2))
                            end do; 
                        end do




                        end do
                        end do
                    end do
                    end do
                end do
            end do

         
            
        end if


            

   
        end if

!     *********
!     Finish up
!     *********
!
        nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus, spAnorm_minus, spAtan1_minus, spAtan2_minus)

    end subroutine Local_GetOffDiagonalBlock






   subroutine Local_Get_BC_Dir_RHS_Poisson(f, nEqn, e_plus,time, side, RHS, AllNDOF, startVarNum)
    implicit none
    !-arguments----------------------------------------------------------------------
    type(Face), target, intent(in)    :: f       !<  Face connecting elements
    type(Element), intent(in)    :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
    integer, intent(in)    :: side    !<  side of face where e_plus is
    integer, intent(in)    :: AllNDOF    !
    REAL(KIND=RP)                :: time
    integer, intent(in)    :: nEqn
    integer, intent(in)    :: startVarNum
    ! real(kind=RP)                       :: RHS( size(e_plus % storage % Q) )

    real(kind=RP)                       :: RHS(nEqn * AllNDOF)


    !-local-variables----------------------------------------------------------------
    integer :: i1, j1, k1, eq1          ! variable counters
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
    integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
    integer :: tanAx_minus(2)           ! Tangent axes to f on e⁻
    integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
    integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
    integer :: faceInd_plus(2)         ! Face indexes on e⁺
    integer :: faceInd_minus(2)         ! Face indexes on e⁻
    integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
    integer :: NxyFace_minus(2)         ! Polynomial orders of face on element e⁻ (only needed for viscous fluxes)
    integer :: faceInd_minus2plus(2)    ! Face indexes on e⁻ passed to the reference frame of e⁺ (only needed for viscous fluxes)
    integer :: Deltas                   ! Number of Kronecker deltas /= 0
    integer :: elInd_plus_norm(3)       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
    integer :: elInd_plus_tan1(3)       ! Element indexes on e⁺, but the index in the first tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
    integer :: elInd_plus_tan2(3)       ! Element indexes on e⁺, but the index in the second tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
    type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
    type(NodalStorage_t), target  :: spA_minus(3)      ! Nodal storage in the different directions for e_minus - local copy
    type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
    type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
    type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
    type(NodalStorage_t), pointer :: spAnorm_minus     ! Nodal storage in the direction that is normal to the face for e⁻
    type(NodalStorage_t), pointer :: spAtan1_minus     ! Nodal storage in the tangent direction "1" to the face for e⁻ (only needed for viscous fluxes)
    type(NodalStorage_t), pointer :: spAtan2_minus     ! Nodal storage in the tangent direction "2" to the face for e⁻ (only needed for viscous fluxes)
    real(kind=RP), pointer :: dfdq(:, :, :, :)     !
    real(kind=RP), pointer :: df_dGradQ_f(:, :, :, :, :)       ! Pointer to the Jacobian with respect to gradQ on face
    real(kind=RP), pointer :: df_dGradQ_e(:, :, :, :, :, :, :) ! Pointer to the Jacobian with respect to gradQ on face
    real(kind=RP) :: a_minus
    real(kind=RP) :: normAux(NCONS, NCONS)
    real(kind=RP), pointer :: Gvec_norm(:, :, :)       ! Auxiliary vector containing values of dFv_dgradQ in the direction normal to the face
    real(kind=RP), pointer :: Gvec_tan1(:, :, :)       ! Auxiliary vector containing values of dFv_dgradQ in the first tangent direction to the face
    real(kind=RP), pointer :: Gvec_tan2(:, :, :)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
    real(kind=RP), allocatable :: nHat(:, :, :)

    ! QDot(:,i1,j1,k1)
    !--------------------------------------------------------------------------------

    ! *********************************************************************************************
    ! integer :: elInd_minus(3)            ! Element indexes on e⁺
    ! integer :: elInd_minus2(3)            ! Element indexes on e⁻
    integer :: i12, j12, k12              ! i1=i2, j1=j2, k1=k2

    real(kind=RP) :: epsilon, sigma      

    integer :: elInd_plus2(3)            ! Element indexes on e⁺
    integer :: elInd_minus2(3)           ! Element indexes on e⁻
    real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
    real(kind=RP) :: normCartesian(3)     ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
    real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
    real(kind=RP) :: tmp01,tmp02    ! the variables in the face after interpolation
    integer :: faceInd_plus2(2)         ! Face indexes on e⁺
    real(kind=RP) :: BcDirichlet(nEqn)     ! Boundary condition values
    real(kind=RP) :: MatEntries_RHS(nEqn)           ! Values of the matrix entries

    ! *********************************************************************************************
    
    
    !
!     ***********
!     Definitions
!     ***********
!

    epsilon = -0.0_RP

    sigma =  ( maxval(f % Nf) )**2_RP *10000_RP / (f % geom % h)**1
    ! write (*,*) "sigma ===BCBC==============", sigma
    BcDirichlet = 3.0

    ! Entry spacing for element e⁺
    nXi1 = e_plus%Nxyz(1) + 1
    nEta1 = e_plus%Nxyz(2) + 1
    EtaSpa1 = nEqn*nXi1
    ZetaSpa1 = nEqn*nXi1*nEta1
    ! write (*,*) "nXi1, nEta1, EtaSpa1, ZetaSpa1 === ", nXi1, nEta1, EtaSpa1, ZetaSpa1

    ! Entry spacing for element e⁻
    nXi2 = e_plus%Nxyz(1) + 1
    nEta2 = e_plus%Nxyz(2) + 1
    EtaSpa2 = nEqn*nXi2
    ZetaSpa2 = nEqn*nXi2*nEta2
    ! write (*,*) "nXi2, nEta2, EtaSpa2, ZetaSpa2 === ", nXi2, nEta2, EtaSpa2, ZetaSpa2

    ! Element sides
    elSide_plus = f%elementSide(side)


    ! Normal and tangent axes
    normAx_plus = normalAxis(elSide_plus)
    tanAx_plus = axisMap(:, f%elementSide(side))

    ! Side of axis where f is
    if (normAx_plus < 0) then
        normAxSide_plus = LEFT
        normAxSide_minus = RIGHT
    else
        normAxSide_plus = RIGHT
        normAxSide_minus = LEFT
       end if
    normAx_plus = abs(normAx_plus)

    ! write(*,*) "((fsideBC---normAx_plus))normAxSide_plus))", normAx_plus,normAxSide_plus


    ! Nodal storage
    ! --------------------------------------
    ! TODO: Why this doesn't work since ifort ver. 19.1?
    ! --------------------------------------
    ! spA_plus  = NodalStorage(e_plus  % Nxyz)
    ! spA_minus = NodalStorage(e_minus % Nxyz)

    spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
    spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
    spA_plus(3) = NodalStorage(e_plus%Nxyz(3))


    spAnorm_plus => spA_plus(normAx_plus)
    spAtan1_plus => spA_plus(tanAx_plus(1))
    spAtan2_plus => spA_plus(tanAx_plus(2))



    ! if (elInd_plus(tanAx_plus(1))== elInd_plus(tanAx_plus(2) ) ) then
    if (e_plus%Nxyz(normAx_plus) /= 0 ) then
    
        ! Polynomial orders
        NxyFace_plus = e_plus%Nxyz(tanAx_plus)
    
        ! write (*,*) "NxyFace_plus, NxyFace_minus=== ", NxyFace_plus, NxyFace_minus
        ! write (*,*) "elInd_plus(normAx_plus) --==--------------=== ",elInd_plus(normAx_plus)
        ! write (*,*) "Hello BC ----------------=== "
    

        if (side == LEFT) then
            normLR =  1.0_RP
        else
            normLR =    1.0_RP
        end if
        ! write (*,*) "normLR ================, ", normLR
        ! write (*,*) "f%geom%jacobian BC ================, ", f%geom%jacobian
        ! write (*,*) "f%geom%normal BC ================, ", f%geom%normal
        ! write (*,*) "normLR BC ================, ", normLR
    
        do k1 = 0, e_plus % Nxyz(3)
            do j1 = 0, e_plus % Nxyz(2)
               do i1 = 0, e_plus % Nxyz(1) 
                  ! write (*,*) ", i,j,k, e % storage % QDot(:,i,j,k)", i,j,k, e % storage % QDot(:,i,j,k)
                  ! write(*,'(A1 1F12.8 A2)', advance="no") " ", e % storage % QDot(:,i,j,k), ", "
            
                   elInd_plus = [i1, j1, k1]
                   elInd_plus2 = [i1, j1, k1]
                   faceInd_plus  = elInd_plus ( tanAx_plus  )
                   faceInd_plus2 = elInd_plus ( tanAx_plus )

      
                   MatEntries_RHS = 0

                   CALL BCs(f % zone) % bc % ScalarState(                                              &
                               f % geom % x(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                               time,                                                                     &
                               f % geom % normal(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                               f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ) &
                            )

                   BcDirichlet = f % storage(2) % Q(startVarNum:startVarNum + nEqn-1,  &
                                                    elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) )

                   write(*,*) "BcDirichlet ==========", BcDirichlet


            
                    if(normAx_plus == 1) then 
                        JGradXi0 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                        JGradXi02 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                    elseIF(normAx_plus == 2) then
                        JGradXi0 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                        JGradXi02 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                    else 
                        JGradXi0 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                        JGradXi02 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                    endif

                
                   MatEntries_RHS = MatEntries_RHS                                          &
                                 + epsilon                                                   &
                                 * spAtan1_plus %w (elInd_plus(tanAx_plus(1)))                            &
                                 * spAtan2_plus %w (elInd_plus(tanAx_plus(2)))   &
                                 * spAnorm_plus %vd(elInd_plus2(normAx_plus), normAxSide_plus) &
                                 * f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))) &
                                 * normLR                                                      &  
                                 * BcDirichlet                                                &
                                 * DOT_PRODUCT(JGradXi02                                         &
                                             / e_plus%geom%jacobian(i1, j1, k1),                &
                                             f%geom%normal(:, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )  &
                                 )

                   MatEntries_RHS = MatEntries_RHS                   &
                                  + 2_RP * sigma                           &
                                  * (f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))) &
                                  * spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                  * spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                  * spAnorm_plus  %v (elInd_plus  (normAx_plus), normAxSide_plus) &
                                  * BcDirichlet

                   MatEntries_RHS = MatEntries_RHS &
                                  / spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                  / spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                  / spAnorm_plus  %w(elInd_plus (normAx_plus)  )  &
                                  / e_plus%geom%jacobian(i1, j1, k1)   

                
                !    RHS(:,i1,j1,k1) = RHS(:,i1,j1,k1) + MatEntries_RHS
                e_plus % storage % slrDot1(:,i1,j1,k1) = e_plus % storage % slrDot1(:,i1,j1,k1) + MatEntries_RHS(1:N_INS)

                ! write (*,*) " ===e_plus % storage % slrDot1(:,i1,j1,k1)==== ",e_plus % storage % slrDot1(:,i1,j1,k1)
                ! write (*,*) " ===e_plus % storage % QDot(:,i1,j1,k1)==== ",e_plus % storage % QDot(:,i1,j1,k1)
                e_plus % storage % QDot(:,i1,j1,k1) = e_plus % storage % QDot(:,i1,j1,k1) + MatEntries_RHS
                ! write (*,*) " ===e_plus % storage % QDot(:,i1,j1,k1)222==== ",e_plus % storage % QDot(:,i1,j1,k1)

               end do
            !    write(*,*)
            end do       
            ! write(*,*)
        end do 
    end if

    !     *********
    !     Finish up
    !     *********
    !
    nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus, spAnorm_minus, spAtan1_minus, spAtan2_minus)

end subroutine Local_Get_BC_Dir_RHS_Poisson !



   subroutine Local_Get_BC_Dir_RHS_Poisson_Specific_variable(f, nEqn, e_plus,time, side, RHS, AllNDOF, variable, startVarNum, dim1, dim2,dim3)
        implicit none
        !-arguments----------------------------------------------------------------------
        type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Element), intent(in)       :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
        integer, intent(in)             :: side    !<  side of face where e_plus is
        integer, intent(in)             :: AllNDOF    !
        integer, intent(in)             :: dim1,dim2,dim3    !
        REAL(KIND=RP)                   :: time
        integer, intent(in)             :: nEqn
        REAL(kind=RP), intent(inout)    :: variable(nEqn, dim1, dim2,dim3)
        integer, intent(in)             :: startVarNum
        ! real(kind=RP)                       :: RHS( size(e_plus % storage % Q) )

        real(kind=RP)                       :: RHS(nEqn * AllNDOF)

        !-local-variables----------------------------------------------------------------
        integer :: i1, j1, k1, eq1          ! variable counters
        integer :: nXi1, nEta1              ! Number of nodes in every direction
        integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
        integer :: nXi2, nEta2              ! Number of nodes in every direction
        integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
        integer :: elSide_plus              ! Element side where f is on e⁻
        integer :: elInd_plus(3)            ! Element indexes on e⁺
        integer :: normAx_plus              ! Normal axis to f on e⁺
        integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
        integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
        integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
        integer :: faceInd_plus(2)         ! Face indexes on e⁺
        integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
        type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
        type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
        type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      

        ! QDot(:,i1,j1,k1)
        !--------------------------------------------------------------------------------
        real(kind=RP) :: epsilon, sigma      

        integer :: elInd_plus2(3)            ! Element indexes on e⁺
        real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        integer :: faceInd_plus2(2)         ! Face indexes on e⁺
        real(kind=RP) :: BcDirichlet(nEqn)     ! Boundary condition values
        real(kind=RP) :: MatEntries_RHS(nEqn)           ! Values of the matrix entries

        ! *********************************************************************************************
    !     ***********
    !     Definitions
    !     ***********
    !

        epsilon = -0.0_RP

        sigma =  ( maxval(f % Nf) )**2_RP *10000_RP / (f % geom % h)**1
        ! write (*,*) "sigma ===BCBC==============", sigma
        BcDirichlet = 3.0

        ! Entry spacing for element e⁺
        nXi1 = e_plus%Nxyz(1) + 1
        nEta1 = e_plus%Nxyz(2) + 1
        EtaSpa1 = nEqn*nXi1
        ZetaSpa1 = nEqn*nXi1*nEta1

        ! Entry spacing for element e⁻
        nXi2 = e_plus%Nxyz(1) + 1
        nEta2 = e_plus%Nxyz(2) + 1
        EtaSpa2 = nEqn*nXi2
        ZetaSpa2 = nEqn*nXi2*nEta2

        ! Element sides
        elSide_plus = f%elementSide(side)

        ! Normal and tangent axes
        normAx_plus = normalAxis(elSide_plus)
        tanAx_plus = axisMap(:, f%elementSide(side))

        ! Side of axis where f is
        if (normAx_plus < 0) then
            normAxSide_plus = LEFT
            normAxSide_minus = RIGHT
        else
            normAxSide_plus = RIGHT
            normAxSide_minus = LEFT
           end if
        normAx_plus = abs(normAx_plus)

        ! Nodal storage
        ! --------------------------------------
        ! TODO: Why this doesn't work since ifort ver. 19.1?
        ! --------------------------------------

        spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
        spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
        spA_plus(3) = NodalStorage(e_plus%Nxyz(3))


        spAnorm_plus => spA_plus(normAx_plus)
        spAtan1_plus => spA_plus(tanAx_plus(1))
        spAtan2_plus => spA_plus(tanAx_plus(2))



        ! if (elInd_plus(tanAx_plus(1))== elInd_plus(tanAx_plus(2) ) ) then
        if (e_plus%Nxyz(normAx_plus) /= 0 ) then
        
            ! Polynomial orders
            NxyFace_plus = e_plus%Nxyz(tanAx_plus)
        

            if (side == LEFT) then
                normLR =  1.0_RP
            else
                normLR =    1.0_RP
            end if
        
            do k1 = 0, e_plus % Nxyz(3)
                do j1 = 0, e_plus % Nxyz(2)
                   do i1 = 0, e_plus % Nxyz(1) 
                
                       elInd_plus = [i1, j1, k1]
                       elInd_plus2 = [i1, j1, k1]
                       faceInd_plus  = elInd_plus ( tanAx_plus  )
                       faceInd_plus2 = elInd_plus ( tanAx_plus )

                
                       MatEntries_RHS = 0

                       CALL BCs(f % zone) % bc % ScalarState(                                              &
                                   f % geom % x(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   time,                                                                     &
                                   f % geom % normal(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ) &
                                )

                       BcDirichlet = f % storage(2) % Q(startVarNum:startVarNum + nEqn-1,  &
                                                        elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) )

                    !    write(*,*) "startVarNum Dir Var4==========", startVarNum
                    !    write(*,*) "BcDirichlet ==========", BcDirichlet

                
                        if(normAx_plus == 1) then 
                            JGradXi0 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                        elseIF(normAx_plus == 2) then
                            JGradXi0 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                        else 
                            JGradXi0 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                        endif

                    
                       MatEntries_RHS = MatEntries_RHS                                          &
                                     + epsilon                                                   &
                                     * spAtan1_plus %w (elInd_plus(tanAx_plus(1)))                            &
                                     * spAtan2_plus %w (elInd_plus(tanAx_plus(2)))   &
                                     * spAnorm_plus %vd(elInd_plus2(normAx_plus), normAxSide_plus) &
                                     * f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))) &
                                     * normLR                                                      &  
                                     * BcDirichlet                                                &
                                     * DOT_PRODUCT(JGradXi02                                         &
                                                 / e_plus%geom%jacobian(i1, j1, k1),                &
                                                 f%geom%normal(:, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )  &
                                     )

                       MatEntries_RHS = MatEntries_RHS                   &
                                      + 2_RP * sigma                           &
                                      * (f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))) &
                                      * spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      * spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      * spAnorm_plus  %v (elInd_plus  (normAx_plus), normAxSide_plus) &
                                      * BcDirichlet

                       MatEntries_RHS = MatEntries_RHS &
                                      / spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      / spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      / spAnorm_plus  %w(elInd_plus (normAx_plus)  )  &
                                      / e_plus%geom%jacobian(i1, j1, k1)   

                    ! if (.not. allocated(variable)) then
                    !     write(*,*) "Error: variable is not allocated!"
                    !     stop
                    ! end if
                                    

                    ! write(*,*) "******** variable ******", variable(:,i1, j1, k1)
                    ! write(*,*) "******** MatEntries_RHS ******", MatEntries_RHS

                    
                    ! e_plus % storage % slrDot1(:,i1,j1,k1) = e_plus % storage % slrDot1(:,i1,j1,k1) + MatEntries_RHS(1:N_INS)

                    ! e_plus % storage % QDot(:,i1,j1,k1)    = e_plus % storage % QDot   (:,i1,j1,k1) + MatEntries_RHS
                    ! e_plus % storage % pre_source(:,i1,j1,k1)    = MatEntries_RHS
                    e_plus % storage % pre_source(:,i1,j1,k1)    = e_plus % storage % pre_source   (:,i1,j1,k1) + MatEntries_RHS
                    ! variable     (:,i1,j1,k1)    = variable   (:,i1,j1,k1) + MatEntries_RHS

                   end do
                !    write(*,*)
                end do       
                ! write(*,*)
            end do 
        end if

        !     *********
        !     Finish up
        !     *********
        !
        nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus)

    end subroutine Local_Get_BC_Dir_RHS_Poisson_Specific_variable !


   subroutine Local_Get_BC_Neumann_RHS_Poisson_Specific_variable(f, nEqn, e_plus,time, side, RHS, AllNDOF, variable, startVarNum, dim1, dim2,dim3)
        implicit none
        !-arguments----------------------------------------------------------------------
        type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Element), intent(in)       :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
        integer, intent(in)             :: side    !<  side of face where e_plus is
        integer, intent(in)             :: AllNDOF    !
        integer, intent(in)             :: dim1,dim2,dim3    !
        REAL(KIND=RP)                   :: time
        integer, intent(in)             :: nEqn
        REAL(kind=RP), intent(inout)    :: variable(nEqn, dim1, dim2,dim3)
        integer, intent(in)             :: startVarNum
        ! real(kind=RP)                       :: RHS( size(e_plus % storage % Q) )

        real(kind=RP)                       :: RHS(nEqn * AllNDOF)

        !-local-variables----------------------------------------------------------------
        integer :: i1, j1, k1, eq1          ! variable counters
        integer :: nXi1, nEta1              ! Number of nodes in every direction
        integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
        integer :: nXi2, nEta2              ! Number of nodes in every direction
        integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
        integer :: elSide_plus              ! Element side where f is on e⁻
        integer :: elInd_plus(3)            ! Element indexes on e⁺
        integer :: normAx_plus              ! Normal axis to f on e⁺
        integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
        integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
        integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
        integer :: faceInd_plus(2)         ! Face indexes on e⁺
        integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
        type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
        type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
        type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      

        ! QDot(:,i1,j1,k1)
        !--------------------------------------------------------------------------------
        real(kind=RP) :: epsilon, sigma      

        integer :: elInd_plus2(3)            ! Element indexes on e⁺
        real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        integer :: faceInd_plus2(2)         ! Face indexes on e⁺
        real(kind=RP) :: BcNeumann(nEqn)     ! Boundary condition values
        real(kind=RP) :: MatEntries_RHS(nEqn)           ! Values of the matrix entries

        ! *********************************************************************************************
    !     ***********
    !     Definitions
    !     ***********
    !

        epsilon = -0.0_RP

        BcNeumann = 0.0

        ! Entry spacing for element e⁺
        nXi1 = e_plus%Nxyz(1) + 1
        nEta1 = e_plus%Nxyz(2) + 1
        EtaSpa1 = nEqn*nXi1
        ZetaSpa1 = nEqn*nXi1*nEta1

        ! Entry spacing for element e⁻
        nXi2 = e_plus%Nxyz(1) + 1
        nEta2 = e_plus%Nxyz(2) + 1
        EtaSpa2 = nEqn*nXi2
        ZetaSpa2 = nEqn*nXi2*nEta2

        ! Element sides
        elSide_plus = f%elementSide(side)

        ! Normal and tangent axes
        normAx_plus = normalAxis(elSide_plus)
        tanAx_plus = axisMap(:, f%elementSide(side))

        ! Side of axis where f is
        if (normAx_plus < 0) then
            normAxSide_plus = LEFT
            normAxSide_minus = RIGHT
        else
            normAxSide_plus = RIGHT
            normAxSide_minus = LEFT
           end if
        normAx_plus = abs(normAx_plus)

        ! Nodal storage
        ! --------------------------------------
        ! TODO: Why this doesn't work since ifort ver. 19.1?
        ! --------------------------------------

        spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
        spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
        spA_plus(3) = NodalStorage(e_plus%Nxyz(3))


        spAnorm_plus => spA_plus(normAx_plus)
        spAtan1_plus => spA_plus(tanAx_plus(1))
        spAtan2_plus => spA_plus(tanAx_plus(2))



        ! if (elInd_plus(tanAx_plus(1))== elInd_plus(tanAx_plus(2) ) ) then
        if (e_plus%Nxyz(normAx_plus) /= 0 ) then
        
            ! Polynomial orders
            NxyFace_plus = e_plus%Nxyz(tanAx_plus)
        

            if (side == LEFT) then
                normLR =  1.0_RP
            else
                normLR =    1.0_RP
            end if
        
            do k1 = 0, e_plus % Nxyz(3)
                do j1 = 0, e_plus % Nxyz(2)
                   do i1 = 0, e_plus % Nxyz(1) 
                
                       elInd_plus = [i1, j1, k1]
                       elInd_plus2 = [i1, j1, k1]
                       faceInd_plus  = elInd_plus ( tanAx_plus  )
                       faceInd_plus2 = elInd_plus ( tanAx_plus )

                
                       MatEntries_RHS = 0

                       CALL BCs(f % zone) % bc % ScalarState(                                              &
                                   f % geom % x(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   time,                                                                     &
                                   f % geom % normal(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ) &
                                )

                       BcNeumann = f % storage(2) % Q(startVarNum:startVarNum + nEqn-1,  &
                                                        elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) )

                    !    write(*,*) "startVarNum NUE 4==========", startVarNum
                    !    write(*,*) "BcNeumann ==========", BcNeumann

                
                        if(normAx_plus == 1) then 
                            JGradXi0 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                        elseIF(normAx_plus == 2) then
                            JGradXi0 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                        else 
                            JGradXi0 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                        endif

                       MatEntries_RHS = MatEntries_RHS                   &
                                      + 1.0_RP                      &
                                      * (f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))) &
                                      * spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      * spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      * spAnorm_plus  %v (elInd_plus  (normAx_plus), normAxSide_plus) &
                                      * BcNeumann

                       MatEntries_RHS = MatEntries_RHS &
                                      / spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      / spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      / spAnorm_plus  %w(elInd_plus (normAx_plus)  )  &
                                      / e_plus%geom%jacobian(i1, j1, k1)   

                        
                        ! write(*,*) "startVarNum ==========", startVarNum
                        ! write(*,*) "bname, BCType ==========", TRIM(BCs(f % zone) % bc % bname), " ",TRIM(BCs(f % zone) % bc % BCType)
                        ! write(*,*) "It is BCNuemann ==========", BcNeumann
               
                               
                    ! e_plus % storage % slrDot1(:,i1,j1,k1) = e_plus % storage % slrDot1(:,i1,j1,k1) + MatEntries_RHS(1:N_INS)

                    ! e_plus % storage % QDot(:,i1,j1,k1)    = e_plus % storage % QDot   (:,i1,j1,k1) + MatEntries_RHS
                    ! e_plus % storage % pre_source(:,i1,j1,k1)    = MatEntries_RHS
                    e_plus % storage % pre_source(:,i1,j1,k1)    = e_plus % storage % pre_source   (:,i1,j1,k1) + MatEntries_RHS
                    ! variable     (:,i1,j1,k1)    = variable   (:,i1,j1,k1) + MatEntries_RHS

                   end do
                !    write(*,*)
                end do       
                ! write(*,*)
            end do 
        end if

        !     *********
        !     Finish up
        !     *********
        !
        nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus)

    end subroutine Local_Get_BC_Neumann_RHS_Poisson_Specific_variable !




    subroutine Local_Get_BC_Dir_RHS_Poisson_3velocity(f, nEqn, e_plus,time, side, RHS, AllNDOF, variable, startVarNum, dim1, dim2,dim3)
        implicit none
        !-arguments----------------------------------------------------------------------
        type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Element), intent(in)       :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
        integer, intent(in)             :: side    !<  side of face where e_plus is
        integer, intent(in)             :: AllNDOF    !
        integer, intent(in)             :: dim1,dim2,dim3    !
        REAL(KIND=RP)                   :: time
        integer, intent(in)             :: nEqn
        REAL(kind=RP), intent(inout)    :: variable(nEqn, dim1, dim2,dim3)
        integer, intent(in)             :: startVarNum
        real(kind=RP)                       :: RHS(nEqn * AllNDOF)

        !-local-variables----------------------------------------------------------------
        integer :: i1, j1, k1, eq1          ! variable counters
        integer :: nXi1, nEta1              ! Number of nodes in every direction
        integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
        integer :: nXi2, nEta2              ! Number of nodes in every direction
        integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
        integer :: elSide_plus              ! Element side where f is on e⁻
        integer :: elInd_plus(3)            ! Element indexes on e⁺
        integer :: normAx_plus              ! Normal axis to f on e⁺
        integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
        integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
        integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
        integer :: faceInd_plus(2)         ! Face indexes on e⁺
        integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
        type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
        type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
        type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      

        ! QDot(:,i1,j1,k1)
        !--------------------------------------------------------------------------------
        real(kind=RP) :: epsilon, sigma      

        integer :: elInd_plus2(3)            ! Element indexes on e⁺
        real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        integer :: faceInd_plus2(2)         ! Face indexes on e⁺
        real(kind=RP) :: BcDirichlet(nEqn)     ! Boundary condition values
        real(kind=RP) :: MatEntries_RHS(nEqn)           ! Values of the matrix entries

        ! *********************************************************************************************
    !     ***********
    !     Definitions
    !     ***********
    !

        epsilon = -0.0_RP

        sigma =  ( maxval(f % Nf) )**2_RP *10000_RP / (f % geom % h)**1
        ! write (*,*) "sigma ===BCBC==============", sigma
        BcDirichlet = 3.0

        ! Entry spacing for element e⁺
        nXi1 = e_plus%Nxyz(1) + 1
        nEta1 = e_plus%Nxyz(2) + 1
        EtaSpa1 = nEqn*nXi1
        ZetaSpa1 = nEqn*nXi1*nEta1

        ! Entry spacing for element e⁻
        nXi2 = e_plus%Nxyz(1) + 1
        nEta2 = e_plus%Nxyz(2) + 1
        EtaSpa2 = nEqn*nXi2
        ZetaSpa2 = nEqn*nXi2*nEta2

        ! Element sides
        elSide_plus = f%elementSide(side)

        ! Normal and tangent axes
        normAx_plus = normalAxis(elSide_plus)
        tanAx_plus = axisMap(:, f%elementSide(side))

        ! Side of axis where f is
        if (normAx_plus < 0) then
            normAxSide_plus = LEFT
            normAxSide_minus = RIGHT
        else
            normAxSide_plus = RIGHT
            normAxSide_minus = LEFT
           end if
        normAx_plus = abs(normAx_plus)

        ! Nodal storage
        ! --------------------------------------
        ! TODO: Why this doesn't work since ifort ver. 19.1?
        ! --------------------------------------

        spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
        spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
        spA_plus(3) = NodalStorage(e_plus%Nxyz(3))


        spAnorm_plus => spA_plus(normAx_plus)
        spAtan1_plus => spA_plus(tanAx_plus(1))
        spAtan2_plus => spA_plus(tanAx_plus(2))



        ! if (elInd_plus(tanAx_plus(1))== elInd_plus(tanAx_plus(2) ) ) then
        if (e_plus%Nxyz(normAx_plus) /= 0 ) then
        
            ! Polynomial orders
            NxyFace_plus = e_plus%Nxyz(tanAx_plus)
        

            if (side == LEFT) then
                normLR =  1.0_RP
            else
                normLR =    1.0_RP
            end if
        
            do k1 = 0, e_plus % Nxyz(3)
                do j1 = 0, e_plus % Nxyz(2)
                   do i1 = 0, e_plus % Nxyz(1) 
                
                       elInd_plus = [i1, j1, k1]
                       elInd_plus2 = [i1, j1, k1]
                       faceInd_plus  = elInd_plus ( tanAx_plus  )
                       faceInd_plus2 = elInd_plus ( tanAx_plus )

                
                       MatEntries_RHS = 0

                       CALL BCs(f % zone) % bc % ScalarState(                                              &
                                   f % geom % x(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   time,                                                                     &
                                   f % geom % normal(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ) &
                                )

                       BcDirichlet = f % storage(2) % Q(startVarNum:startVarNum + nEqn-1,  &
                                                        elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) )


                        
                    !    write(*,*) "startVarNum ==========", startVarNum
                    !    write(*,*) "bname, BCType ==========", TRIM(BCs(f % zone) % bc % bname), " ",TRIM(BCs(f % zone) % bc % BCType)

                    !    character(len=LINE_LENGTH) :: bname
                    !    character(len=LINE_LENGTH) :: BCType
                    

                    !    write(*,*) "It is BcDirichlet ==========", BcDirichlet

                
                        if(normAx_plus == 1) then 
                            JGradXi0 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                        elseIF(normAx_plus == 2) then
                            JGradXi0 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                        else 
                            JGradXi0 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                        endif

                    
                       MatEntries_RHS = MatEntries_RHS                                          &
                                     + epsilon                                                   &
                                     * spAtan1_plus %w (elInd_plus(tanAx_plus(1)))                            &
                                     * spAtan2_plus %w (elInd_plus(tanAx_plus(2)))   &
                                     * spAnorm_plus %vd(elInd_plus2(normAx_plus), normAxSide_plus) &
                                     * f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))) &
                                     * normLR                                                      &  
                                     * BcDirichlet                                                &
                                     * DOT_PRODUCT(JGradXi02                                         &
                                                 / e_plus%geom%jacobian(i1, j1, k1),                &
                                                 f%geom%normal(:, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )  &
                                     )

                       MatEntries_RHS = MatEntries_RHS                   &
                                      + 2_RP * sigma                           &
                                      * (f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))) &
                                      * spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      * spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      * spAnorm_plus  %v (elInd_plus  (normAx_plus), normAxSide_plus) &
                                      * BcDirichlet

                       MatEntries_RHS = MatEntries_RHS &
                                      / spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      / spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      / spAnorm_plus  %w(elInd_plus (normAx_plus)  )  &
                                      / e_plus%geom%jacobian(i1, j1, k1)   
                    
                    ! write(*,*) "e_plus % storage % vel_source(:,i1,j1,k1) =2=========", e_plus % storage % vel_source(:,i1,j1,k1)
                    e_plus % storage % vel_source(:,i1,j1,k1)    = e_plus % storage % vel_source   (:,i1,j1,k1) + MatEntries_RHS
                    ! write(*,*) "e_plus % storage % vel_source(:,i1,j1,k1) =3=========", e_plus % storage % vel_source(:,i1,j1,k1)

                   end do
                !    write(*,*)
                end do       
                ! write(*,*)
            end do 
        end if

        !     *********
        !     Finish up
        !     *********
        !
        nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus)

    end subroutine Local_Get_BC_Dir_RHS_Poisson_3velocity !

    subroutine Local_Get_BC_Nuemann_RHS_Poisson_3velocity(f, nEqn, e_plus,time, side, RHS, AllNDOF, variable, startVarNum, dim1, dim2,dim3)
        implicit none
        !-arguments----------------------------------------------------------------------
        type(Face), target, intent(in)    :: f       !<  Face connecting elements
        type(Element), intent(in)       :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
        integer, intent(in)             :: side    !<  side of face where e_plus is
        integer, intent(in)             :: AllNDOF    !
        integer, intent(in)             :: dim1,dim2,dim3    !
        REAL(KIND=RP)                   :: time
        integer, intent(in)             :: nEqn
        REAL(kind=RP), intent(inout)    :: variable(nEqn, dim1, dim2,dim3)
        integer, intent(in)             :: startVarNum
        real(kind=RP)                       :: RHS(nEqn * AllNDOF)

        !-local-variables----------------------------------------------------------------
        integer :: i1, j1, k1, eq1          ! variable counters
        integer :: nXi1, nEta1              ! Number of nodes in every direction
        integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
        integer :: nXi2, nEta2              ! Number of nodes in every direction
        integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
        integer :: elSide_plus              ! Element side where f is on e⁻
        integer :: elInd_plus(3)            ! Element indexes on e⁺
        integer :: normAx_plus              ! Normal axis to f on e⁺
        integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
        integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
        integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
        integer :: faceInd_plus(2)         ! Face indexes on e⁺
        integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
        type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
        type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
        type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
        type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      

        ! QDot(:,i1,j1,k1)
        !--------------------------------------------------------------------------------
        real(kind=RP) :: epsilon, sigma      

        integer :: elInd_plus2(3)            ! Element indexes on e⁺
        real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
        integer :: faceInd_plus2(2)         ! Face indexes on e⁺
        real(kind=RP) :: BCNuemann(nEqn)     ! Boundary condition values
        real(kind=RP) :: MatEntries_RHS(nEqn)           ! Values of the matrix entries

        ! *********************************************************************************************
    !     ***********
    !     Definitions
    !     ***********
    !


        BCNuemann = 0.0

        ! Entry spacing for element e⁺
        nXi1 = e_plus%Nxyz(1) + 1
        nEta1 = e_plus%Nxyz(2) + 1
        EtaSpa1 = nEqn*nXi1
        ZetaSpa1 = nEqn*nXi1*nEta1

        ! Entry spacing for element e⁻
        nXi2 = e_plus%Nxyz(1) + 1
        nEta2 = e_plus%Nxyz(2) + 1
        EtaSpa2 = nEqn*nXi2
        ZetaSpa2 = nEqn*nXi2*nEta2

        ! Element sides
        elSide_plus = f%elementSide(side)

        ! Normal and tangent axes
        normAx_plus = normalAxis(elSide_plus)
        tanAx_plus = axisMap(:, f%elementSide(side))

        ! Side of axis where f is
        if (normAx_plus < 0) then
            normAxSide_plus = LEFT
            normAxSide_minus = RIGHT
        else
            normAxSide_plus = RIGHT
            normAxSide_minus = LEFT
           end if
        normAx_plus = abs(normAx_plus)

        ! Nodal storage
        ! --------------------------------------
        ! TODO: Why this doesn't work since ifort ver. 19.1?
        ! --------------------------------------

        spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
        spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
        spA_plus(3) = NodalStorage(e_plus%Nxyz(3))


        spAnorm_plus => spA_plus(normAx_plus)
        spAtan1_plus => spA_plus(tanAx_plus(1))
        spAtan2_plus => spA_plus(tanAx_plus(2))



        ! if (elInd_plus(tanAx_plus(1))== elInd_plus(tanAx_plus(2) ) ) then
        if (e_plus%Nxyz(normAx_plus) /= 0 ) then
        
            ! Polynomial orders
            NxyFace_plus = e_plus%Nxyz(tanAx_plus)
        

            if (side == LEFT) then
                normLR =  1.0_RP
            else
                normLR =    1.0_RP
            end if
        
            do k1 = 0, e_plus % Nxyz(3)
                do j1 = 0, e_plus % Nxyz(2)
                   do i1 = 0, e_plus % Nxyz(1) 
                
                       elInd_plus = [i1, j1, k1]
                       elInd_plus2 = [i1, j1, k1]
                       faceInd_plus  = elInd_plus ( tanAx_plus  )
                       faceInd_plus2 = elInd_plus ( tanAx_plus )

                
                       MatEntries_RHS = 0

                       CALL BCs(f % zone) % bc % ScalarState(                                              &
                                   f % geom % x(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   time,                                                                     &
                                   f % geom % normal(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                   f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ) &
                                )

                       BCNuemann = f % storage(2) % Q(startVarNum:startVarNum + nEqn-1,  &
                                                        elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) )


                        
                    !    write(*,*) "startVarNum ==========", startVarNum
                    !    write(*,*) "bname, BCType ==========", TRIM(BCs(f % zone) % bc % bname), " ",TRIM(BCs(f % zone) % bc % BCType)

                    !    character(len=LINE_LENGTH) :: bname
                    !    character(len=LINE_LENGTH) :: BCType
                    

                    !    write(*,*) "It is BCNuemann ==========", BCNuemann

                
                        if(normAx_plus == 1) then 
                            JGradXi0 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                        elseIF(normAx_plus == 2) then
                            JGradXi0 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                        else 
                            JGradXi0 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                            JGradXi02 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                        endif


                       MatEntries_RHS = MatEntries_RHS                   &
                                      + 1.0_RP                           &
                                      * (f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))) &
                                      * spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      * spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      * spAnorm_plus  %v (elInd_plus  (normAx_plus), normAxSide_plus) &
                                      * BCNuemann

                       MatEntries_RHS = MatEntries_RHS &
                                      / spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                      / spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                      / spAnorm_plus  %w(elInd_plus (normAx_plus)  )  &
                                      / e_plus%geom%jacobian(i1, j1, k1)   
                    
                    e_plus % storage % vel_source(:,i1,j1,k1)    = e_plus % storage % vel_source   (:,i1,j1,k1) + MatEntries_RHS

                   end do
                !    write(*,*)
                end do       
                ! write(*,*)
            end do 
        end if

        !     *********
        !     Finish up
        !     *********
        !
        nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus)

    end subroutine Local_Get_BC_Nuemann_RHS_Poisson_3velocity !
 


! -------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the diagonal blocks of the Jacobian matrix
!  -------------------------------------------------------------------------------------------
    subroutine AnalyticalJacobianPoisson_Neu_Per_setOne(mesh, nEqn, time, Matrix)
        implicit none
        !--------------------------------------------
        type(HexMesh), target, intent(inout) :: mesh
        integer, intent(in)    :: nEqn
        real(kind=RP), intent(in)    :: time
        class(Matrix_t), intent(inout) :: Matrix
        !--------------------------------------------
        integer :: eID, fID
        type(Element), pointer :: e
        !--------------------------------------------
!
!     Compute each element's diagonal block
!     -------------------------------------
!$omp do schedule(runtime) private(e)
        do eID = 1, size(mesh%elements)
            ! write (*,*) "eID ============= ", eID
            e => mesh%elements(eID)
            ! write (*, *) "F B B R T L ",e%faceIDs(EFRONT), e%faceIDs(EBACK), e%faceIDs(EBOTTOM), & 
            !                             e%faceIDs(ERIGHT), e%faceIDs(ETOP), e%faceIDs(ELEFT)
            call Local_Get_1Line_SetOne(e, neqn,&
                                        mesh%faces(e%faceIDs(EFRONT)), &
                                        mesh%faces(e%faceIDs(EBACK)), &
                                        mesh%faces(e%faceIDs(EBOTTOM)), &
                                        mesh%faces(e%faceIDs(ERIGHT)), &
                                        mesh%faces(e%faceIDs(ETOP)), &
                                        mesh%faces(e%faceIDs(ELEFT)), &
                                        Matrix, &
                                        mesh%IBM)
        end do
!$omp end do
        nullify (e)
    end subroutine AnalyticalJacobianPoisson_Neu_Per_setOne



    !  -----------------------------------------------------
!  Local_SetDiagonalBlock:
!     Subroutine to set a diagonal block in the Jacobian
!  -----------------------------------------------------
subroutine Local_Get_1Line_SetOne(e, neqn, fF, fB, fO, fR, fT, fL, Matrix, IBM)
    use IBMClass
    implicit none
    !-------------------------------------------
    type(Element), intent(inout) :: e
    type(Face), target, intent(in)    :: fF, fB, fO, fR, fT, fL !< The six faces of the element
    class(Matrix_t), intent(inout) :: Matrix
    type(IBM_type), intent(inout) :: IBM
    integer, intent(in)    :: nEqn
    !-------------------------------------------
    real(kind=RP) :: MatEntries(neqn, neqn)
    integer :: i, j             ! Matrix indices
    integer :: i1, j1, k1, eq1  ! variable counters for row
    integer :: i2, j2, k2, eq2  ! variable counters for column
    integer :: i12, j12, k12    ! variable counters when (row index) = (column index)
    integer :: baseRow, baseCol ! Position of neqn by neqn miniblock of Jacobian
    integer :: nXi, nEta        ! Number of nodes in every direction
    integer :: EtaSpa, ZetaSpa  ! Spacing for these two coordinate directions
    real(kind=RP) :: temp


    integer :: index  ! variable counters for column

    type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

    optional :: IBM

    !-------------------------------------------
    spAxi => NodalStorage(e%Nxyz(1))
    spAeta => NodalStorage(e%Nxyz(2))
    spAzeta => NodalStorage(e%Nxyz(3))
!
!     *******************
!     Initial definitions
!     *******************
!
    nXi = e%Nxyz(1) + 1
    nEta = e%Nxyz(2) + 1
    EtaSpa = neqn*nXi
    ZetaSpa = neqn*nXi*nEta

 
    !     Xi contribution (dj*dk)
    do k12 = 0, e%Nxyz(3)
        do j12 = 0, e%Nxyz(2)
            do i12 = 0, e%Nxyz(1)
                    ! MatEntries = 10
                    MatEntries = 1e-8

                    baseRow = i12*neqn + j12*EtaSpa + k12*ZetaSpa
                    baseCol = i12*neqn + j12*EtaSpa + k12*ZetaSpa

                    do eq2 = 1, neqn
                        do eq1 = 1, neqn
                            i = eq1 + baseRow  ! row index (1-based)
                            ! i = 1  ! row index (1-based)
                            j = eq2 + baseCol  ! column index (1-based)

                            call Matrix%SetBlockEntry(1, e%GlobID, 13, j, 0.0_RP)  
                            ! call Matrix%SetBlockEntry(23, e%GlobID, 13, j, 0.0_RP)  
                            ! call Matrix%SetBlockEntry(1, e%GlobID, 1, j, 2e-9_RP)  
                            ! call Matrix%AddToBlockEntry(1, e%GlobID, 1, j, 10.0_RP)  
                            ! call Matrix%SetBlockEntry(1, e%GlobID, 1, j, MatEntries(eq1, eq2))  
                            ! TODO: This can be improved by setting the whole matrix at once

                        end do
                    end do
                    ! =========--------===*****************************======================
            end do
        end do
    end do

    call Matrix%SetBlockEntry(1, 1, 13, 13, 1e9_RP)


    nullify (spAxi, spAeta, spAzeta)
end subroutine Local_Get_1Line_SetOne




#endif

#if defined(NAVIERSTOKES)  


!-------------------------------------===000level NS --------------------------------------------------------
#endif
end module AnalyticalJacobianPoisson
