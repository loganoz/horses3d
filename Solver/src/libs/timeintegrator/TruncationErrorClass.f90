module TruncationErrorClass
   use SMConstants
   use MultigridTypes
   use DGSEMClass
   use FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use PhysicsStorage
   use FluidData
   use HexMeshClass
   use NodalStorageClass
#if defined(CAHNHILLIARD)
   use BoundaryConditionFunctions, only: C_BC, MU_BC
#endif
   implicit none
   
   private
   public TruncationError_t, EstimateTruncationError, InitializeForTauEstimation, GenerateExactTEmap
   public NON_ISOLATED_TE, ISOLATED_TE
   
   !---------------------------------------------------------------------------------------------------------
   ! Type for storing the truncation error for one element in one direction according to the polynomial order
   !  -> Currently, we are only storing the infinity norm. This needs less storage and is indeed  a conservative 
   !     creterion |\tau| <= |\tau_1| + |\tau_2| (see paper "Truncation Error Estimation in the p-Anisotropic Discontinuous Galerkin Spectral Element Method")
   !---------------------------------------------------------------------------------------------------------
   type :: TruncErrorPol_t
      real(kind=RP), allocatable  :: maxTE(:)   ! |\tau|∞
      integer                     :: P          ! Polynomial order in this direction
   end type TruncErrorPol_t
   
   !-----------------------------------------------------------------------------
   ! Type for storing the truncation error of one element in the three directions 
   !-----------------------------------------------------------------------------
   type :: TruncationError_t
      type(TruncErrorPol_t)       :: Dir(3)
      integer                     :: TruncErrorType
      contains
         procedure :: construct => ConstructTruncationError
         procedure :: destruct  => DestructTruncationError
   end type TruncationError_t

!
!  -------------------
!  External procedures
!  -------------------
!
#if defined(NAVIERSTOKES)
   procedure(BCState_FCN)   :: externalStateForBoundaryName_NS
   procedure(BCGradients_FCN)   :: ExternalGradientForBoundaryName_NS
#elif defined(CAHNHILLIARD)
   procedure(BCState_FCN)   :: externalCHStateForBoundaryName
   procedure(BCGradients_FCN)   :: ExternalChemicalPotentialGradientForBoundaryName
   procedure(BCGradients_FCN)   :: ExternalConcentrationGradientForBoundaryName
#endif
   
#if defined(NAVIERSTOKES)
   interface
      subroutine UserDefinedSourceTermNS(x, time, S, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         USE HexMeshClass
         use PhysicsStorage
         use FluidData
         IMPLICIT NONE
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(out)  :: S(NCONS)
         type(Thermodynamics_t),    intent(in)  :: thermodynamics_
         type(Dimensionless_t),     intent(in)  :: dimensionless_
         type(RefValues_t),         intent(in)  :: refValues_
      end subroutine UserDefinedSourceTermNS
   end interface
#endif
!
!  ----------------
!  Module variables
!  ----------------
!
   procedure(ComputeQDot_FCN), pointer :: TimeDerivative
   
   !! Parameters
   integer, parameter :: ISOLATED_TE = 0
   integer, parameter :: NON_ISOLATED_TE = 1
   
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ROUTINES FOR TRUNCATION ERROR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that constructs the variable that stores tau in one element
!  ----------------------------------------------------------------------
   subroutine ConstructTruncationError(TE,Nmin,Nmax)
      implicit none
      !------------------------------------------
      class(TruncationError_t), intent(inout) :: TE         !<> Variable that stores the truncation error
      integer                 , intent(in)    :: Nmin       !<  Minimum polynomial order for estimation
      integer                 , intent(in)    :: Nmax(3)    !<  Maximum polynomial order for estimation
      !------------------------------------------
      
      !------------------------------------------
      
      ! Allocate maxTE
      allocate (TE % Dir(1) % maxTE(NMin:Nmax(1)))
      allocate (TE % Dir(2) % maxTE(NMin:Nmax(2)))
      allocate (TE % Dir(3) % maxTE(NMin:Nmax(3)))
!
!     --------------------------
!     Initialize maxTE to zero 
!        Remarks:
!           a. This value will be overwritten when the TE is estimated.
!           b. If the polynomial order specified by the user is NMin, no TE estimation can be done and therefore the value TE = 0 is assumed.
!     --------------------------
!
      TE % Dir(1) % maxTE = 0._RP
      TE % Dir(2) % maxTE = 0._RP
      TE % Dir(3) % maxTE = 0._RP
   end subroutine ConstructTruncationError
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that destructs the variable that stores tau in one element
!  ----------------------------------------------------------------------
   subroutine DestructTruncationError(TE)
      implicit none
      !------------------------------------------
      class(TruncationError_t) :: TE     !<> Variable that stores the truncation error
      !------------------------------------------
      
      deallocate (TE % Dir(1) % maxTE)
      deallocate (TE % Dir(2) % maxTE)
      deallocate (TE % Dir(3) % maxTE)
      
   end subroutine DestructTruncationError
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that sets P for all elements in TE
!  ----------------------------------------------------------------------
   subroutine InitializeForTauEstimation(TE,sem,TruncErrorType, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      implicit none
      !------------------------------------------
      type(TruncationError_t) :: TE(:)
      type(DGSem), intent(in) :: sem
      integer    , intent(in) :: TruncErrorType !<  Either NON_ISOLATED_TE or ISOLATED_TE
      procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      procedure(ComputeQDot_FCN) :: ComputeTimeDerivativeIsolated
      !------------------------------------------
      integer                 :: eID
      !------------------------------------------
      
      do eID=1, size(sem % mesh % elements)
         TE(eID) % TruncErrorType = TruncErrorType
         TE(eID) % Dir(1) % P = sem % Nx(eID)
         TE(eID) % Dir(2) % P = sem % Ny(eID)
         TE(eID) % Dir(3) % P = sem % Nz(eID)
      end do
      
      if (TruncErrorType == ISOLATED_TE) then
         TimeDerivative => ComputeTimeDerivativeIsolated
      else
         TimeDerivative => ComputeTimeDerivative
      end if
      
   end subroutine InitializeForTauEstimation
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Routine for computing the truncation error associated with every element
!  in a specified direction for a specific evaluation of the fine solution
!     (a-posteriori or quasi a-priori without correction)
!  ------------------------------------------------------------------------
   subroutine EstimateTruncationError(TE,sem,t,Var,Dir)
      implicit none
      !--------------------------------------------------------
      type(TruncationError_t)  :: TE(:)       !<> Type containing the truncation error estimation
      type(DGSem)              :: sem         !<> sem (to evaluate Qdot in a given coarser mesh)
      real(kind=RP)            :: t           !<  time 
      type(MGSolStorage_t)     :: Var(:)      !<  Type containing the source term in the mesh where the truncation error is being estimated (to be deprecated.. maintained just to keep consistency with manufactured solutions module)
      integer                  :: Dir         !<  Direction in which the truncation error is being estimated
      !--------------------------------------------------------
      integer                  :: iEl           !   Element counter
      integer                  :: iEQ           !   Equation counter
      integer                  :: i,j,k         !   Coordinate index counters
      real(kind=RP)            :: wx, wy, wz    ! 
      real(kind=RP)            :: Jac
      real(kind=RP)            :: maxTE
      real(kind=RP)            :: S(NTOTALVARS)      !   Source term
      !--------------------------------------------------------
      
      call TimeDerivative(sem % mesh, sem % particles, t, sem % BCFunctions)
      
      S = 0._RP ! Initialize source term
      
      do iEl = 1, size(sem % mesh % elements)
         associate (e => sem % mesh % elements(iEl))
         
         if (e % Nxyz(Dir) >= TE(iEl) % Dir(Dir) % P) cycle ! it is actually never going to be greater than... but for security
         
         maxTE = 0._RP
         
         ! loop over all the degrees of freedom of the element
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
#if defined(NAVIERSTOKES)            
            call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), t, S, thermodynamics, dimensionless, refValues)
#endif
            
            do iEQ = 1, NTOTALVARS
               wx  = NodalStorage(e % Nxyz(1)) % w (i)
               wy  = NodalStorage(e % Nxyz(2)) % w (j)
               wz  = NodalStorage(e % Nxyz(3)) % w (k)
               Jac = e % geom % jacobian(i,j,k)
               
               maxTE =  MAX(maxTE , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,i,j,k) + S(iEQ) + Var(iEl) % Scase(iEQ,i,j,k) )  )
            end do
         end do         ; end do         ; end do
         
         TE(iEl) % Dir(Dir) % maxTE(e % Nxyz(Dir)) = maxTE
         
         end associate
      end do
   end subroutine EstimateTruncationError
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for printing the TE map(s) of one element
!  -----------------------------------------------------------------------
   subroutine PrintTEmap(TEmap,iEl)
      implicit none
      !-------------------------------------------
      real(kind=RP)  :: TEmap(:,:,:)
      integer        :: iEl
      !-------------------------------------------
      integer                :: k, i, l
      integer                :: fd
      character(LINE_LENGTH) :: TEmapfile
      !-------------------------------------------
      
      do k = 1, size(TEmap,3)
         write(TEmapfile,'(A,I7.7,A,I2.2,A)') 'RegressionFiles/TEmapXY-Elem_',iEl,'-Nz_',k,'.dat'
   
         open(newunit = fd, file=TRIM(TEmapfile), action='write')
            do i = 1, size(TEmap, 1)
               write(fd,*) (TEmap(i,l,k),l=1,size(TEmap,2))
            end do
         close(fd)
      end do
   end subroutine PrintTEmap
!
!  ------------------------------------------------------------------------
!  Subroutine for generating the exact  truncation error map in one element
!  -> The exact solution must be coded down in UserDefinedState1 (ProblemFile.f90)
!  ------------------------------------------------------------------------
!
   subroutine GenerateExactTEmap(sem, NMIN, NMAX, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, iEl, TruncErrorType)
      implicit none
      !-------------------------------------------------------------------------
      type(DGSem)                :: sem
      integer, intent(in)        :: NMIN
      integer, intent(in)        :: NMAX(NDIM)
      real(kind=RP)              :: t
      procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      procedure(ComputeQDot_FCN) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables
      integer, intent(in)        :: iEl
      integer, intent(in)        :: TruncErrorType
      !-------------------------------------------------------------------------
      real(kind=RP)              :: TEmap (NMIN:NMAX(1),NMIN:NMAX(2),NMIN:NMAX(3))
      integer                    :: i,j,k
      integer                    :: nelem      ! Number of elements
      type(BCFunctions_t)        :: BCFunctions(no_of_BCsets)
!~       type(SolStorage_t), allocatable :: Var(:)  
      logical                              :: success   
      integer, dimension(size(sem % mesh % elements) ) :: Nx, Ny, Nz
!~       integer                              :: eID
!~       integer                              :: ii,jj,kk
      
!~       real(kind=RP)  :: maxTE, wx, wy, wz, Jac, maxTEtemp
!~       integer        :: iEQ, ElMax
      !-------------------------------------------------------------------------
      
      ! Initializations
      nelem   = SIZE(sem % mesh % elements)
      
      if (TruncErrorType == ISOLATED_TE) then
         TimeDerivative => ComputeTimeDerivativeIsolated
      else
         TimeDerivative => ComputeTimeDerivative
      end if
      
#if defined(NAVIERSTOKES)
      BCFunctions(1) % externalState => externalStateForBoundaryName_NS
      BCFunctions(1) % externalGradients => externalGradientForBoundaryName_NS
#elif defined(CAHNHILLIARD)
      BCFunctions(C_BC) % externalState      => externalCHStateForBoundaryName
      BCFunctions(C_BC) % externalGradients  => externalConcentrationGradientForBoundaryName

      BCFunctions(MU_BC) % externalState     => externalCHStateForBoundaryName
      BCFunctions(MU_BC) % externalGradients => externalChemicalPotentialGradientForBoundaryName
#endif
      
      do k = NMIN, NMAX(3)
         do j = NMIN, NMAX(2)
            do i = NMIN, NMAX(1)
               
               Nx = i
               Ny = j
               Nz = k
               
               !-------------------------------------------
               !Destruct previous sem and construct new one
               !-------------------------------------------
               CALL sem % destruct()
               
               sem % ManufacturedSol = controlVariables % containsKey("manufactured solution")
               
               call sem % construct (  controlVariables  = controlVariables                              ,  &
                           BCFunctions       = BCFunctions, &
                           Nx_ = Nx ,     Ny_ = Ny,     Nz_ = Nz,                   &
                           success           = success)
               
               print*, 'Computing TE for N=',i,j,k,'. success=', success
               
               if(.NOT. success)   ERROR STOP ":: problem creating sem"
               
#if defined(NAVIERSTOKES)
               CALL UserDefinedFinalSetup(sem % mesh , thermodynamics, dimensionless, refValues)
#endif
               
               
               TEmap(i,j,k) = EstimateTauOfElem(sem,t,controlVariables,iEl)
               
               print*, 'Done for N=',i,j,k
            end do
         end do
      end do
      
      CALL PrintTEmap(TEmap,iEl)
      stop
      
   end subroutine GenerateExactTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----
!     Estimates the infinity norm of the truncation error for a given mesh given a restart solution or a manufactured solution
!     ---
      function EstimateTauOfElem(sem,t,controlVariables,iEl) result(maxTE)
         implicit none
         !-------------------------------------------------------
         type(DGSem), target    , intent(inout) :: sem              !<> sem class (inout cause' we compute Qdot)
         real(kind=RP)          , intent(in)    :: t                !>  Time
         type(FTValueDictionary), intent(in)    :: controlVariables !<
         integer                , intent(in)    :: iEl              !<  Present if the result is wanted for a certain element
         real(kind=RP)                          :: maxTE            !>  |\tau|_{\infty}
         !-------------------------------------------------------
         integer          :: nelem
         integer          :: eID             ! element counter
         integer          :: iEQ             ! Equation counter
         integer          :: ii,jj,kk        ! doF counters
         real(kind=RP)    :: wx, wy, wz      ! Quadrature weights
         real(kind=RP)    :: Jac             ! Jacobian (mapping)
         !-------------------------------------------------------
         
         nelem = SIZE(sem % mesh % elements)
         
         !-------------------------------------------
         !Get exact solution (from ProblemFile.f90)
         !-------------------------------------------
         
         do eID = 1, nelem
            associate (e => sem % mesh % elements(eID))
#if defined(NAVIERSTOKES)            
            do kk = 0, e % Nxyz(3) ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
               call UserDefinedState1(e % geom % x(:,ii,jj,kk), t, [0._RP, 0._RP, 0._RP], e % storage % Q(:,ii,jj,kk), thermodynamics, dimensionless, refValues)
            end do                 ; end do                 ; end do
#endif
            
            end associate
         end do
         
         call TimeDerivative(sem % mesh, sem % particles, t, sem % BCFunctions)
         
         maxTE = 0._RP ! Initialization
         
         associate ( e => sem % mesh % elements(iEl) )
         
         ! loop over all the degrees of freedom of the element
         do kk = 0, e % Nxyz(3) ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            
            do iEQ = 1, NTOTALVARS
               wx  = NodalStorage(e % Nxyz(1)) % w (ii)
               wy  = NodalStorage(e % Nxyz(2)) % w (jj)
               wz  = NodalStorage(e % Nxyz(3)) % w (kk)
               Jac = e % geom % jacobian(ii,jj,kk)
      
               maxTE =  MAX(maxTE , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,ii,jj,kk) ) )
            end do
         end do          ; end do          ; end do
          
         end associate
      end function EstimateTauOfElem
end module TruncationErrorClass
