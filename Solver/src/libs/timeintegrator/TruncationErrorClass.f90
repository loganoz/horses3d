module TruncationErrorClass
   use SMConstants
   use MultigridTypes            , only: MGSolStorage_t
   use DGSEMClass                , only: DGSem, ComputeTimeDerivative_f
   use FTValueDictionaryClass    , only: FTValueDictionary
   use PhysicsStorage            , only: NCONS, CTD_IGNORE_MODE
   use HexMeshClass              , only: HexMesh
#if defined(NAVIERSTOKES) || defined(INCNS)
   use FluidData                 , only: Thermodynamics, RefValues, Dimensionless
   use ProblemFileFunctions      , only: UserDefinedState_f, UserDefinedFinalSetup_f, UserDefinedSourceTermNS_f
#endif
   use NodalStorageClass         , only: NodalStorage
   use FileReadingUtilities      , only: RemovePath
   use Utilities                 , only: AlmostEqual, LeastSquaresLinRegression
   use ElementClass              , only: Element
   implicit none
   
   private
   public TruncationError_t, EstimateTruncationError, InitializeForTauEstimation, PrintTEmap, AssignTimeDerivative, GenerateExactTEmap, EstimateAndPlotTruncationError
   public NON_ISOLATED_TE, ISOLATED_TE, NMINest, OLD_TE, NEW_TE
   
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
      integer                     :: TruncErrorForm
      contains
         procedure :: construct           => TruncationError_Construct
         procedure :: destruct            => TruncationError_Destruct
         procedure :: GenerateTEmap       => TruncationError_GenerateTEmap
         procedure :: ExportToFile        => TruncationError_ExportToFile
         procedure :: ReadFromFile        => TruncationError_ReadFromFile
         procedure :: ExtrapolateInOneDir => TruncationError_ExtrapolateInOneDir
         procedure :: reset               => TruncationError_reset
   end type TruncationError_t
!
!  ----------------
!  Module variables
!  ----------------
!
   procedure(ComputeTimeDerivative_f), pointer :: TimeDerivative
   
   !! Parameters
   integer, parameter :: ISOLATED_TE = 0
   integer, parameter :: NON_ISOLATED_TE = 1
   integer, parameter :: NMINest    = 1      ! Minimum polynomial order used for estimation 
   integer, parameter :: OLD_TE = 2
   integer, parameter :: NEW_TE = 3
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
   subroutine TruncationError_Construct(TE,Nmax)
      implicit none
      !------------------------------------------
      class(TruncationError_t), intent(inout) :: TE         !<> Variable that stores the truncation error
      integer                 , intent(in)    :: Nmax(3)    !<  Maximum polynomial order for estimation
      !------------------------------------------
      
      !------------------------------------------
      
      ! Allocate maxTE
      allocate (TE % Dir(1) % maxTE(NMINest:Nmax(1)))
      allocate (TE % Dir(2) % maxTE(NMINest:Nmax(2)))
      allocate (TE % Dir(3) % maxTE(NMINest:Nmax(3)))
!
!     --------------------------
!     Initialize maxTE to zero 
!        Remarks:
!           a. This value will be overwritten when the TE is estimated.
!           b. If the polynomial order specified by the user is NMINest, no TE estimation can be done and therefore the value TE = 0 is assumed.
!     --------------------------
!
      TE % Dir(1) % maxTE = 0._RP
      TE % Dir(2) % maxTE = 0._RP
      TE % Dir(3) % maxTE = 0._RP
   end subroutine TruncationError_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that destructs the variable that stores tau in one element
!  ----------------------------------------------------------------------
   elemental subroutine TruncationError_Destruct(TE)
      implicit none
      !------------------------------------------
      class(TruncationError_t), intent(inout) :: TE     !<> Variable that stores the truncation error
      !------------------------------------------
      
      deallocate (TE % Dir(1) % maxTE)
      deallocate (TE % Dir(2) % maxTE)
      deallocate (TE % Dir(3) % maxTE)
      
   end subroutine TruncationError_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------
!  Export truncation error to file
!  -------------------------------
   subroutine TruncationError_ExportToFile(this, fName)
      implicit none
      !-arguments-------------------------------------------------
      class(TruncationError_t)   , intent(in) :: this
      character(len=LINE_LENGTH) , intent(in) :: fName
      !-local-variables-------------------------------------------
      integer :: fd
      !-----------------------------------------------------------
      
      open (newunit=fd, file = trim(fName), position='append') 
            
      write(fd,*) this % Dir(1) % P, this % Dir(2) % P, this % Dir(3) % P
      write(fd,*) this % Dir(1) % maxTE
      write(fd,*) this % Dir(2) % maxTE
      write(fd,*) this % Dir(3) % maxTE
      
      close(fd)
      
   end subroutine TruncationError_ExportToFile
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------
!  Read truncation error from file
!  -------------------------------
   subroutine TruncationError_ReadFromFile(this, fName, fd_in)
      implicit none
      !-arguments-------------------------------------------------
      class(TruncationError_t)            , intent(inout)   :: this
      character(len=LINE_LENGTH), optional, intent(in)      :: fName
      integer                   , optional, intent(in)      :: fd_in 
      !-local-variables-------------------------------------------
      integer :: fd, Pxyz(3)
      !-----------------------------------------------------------
      
      if (present(fName) ) then
         open (newunit=fd, file = trim(fName)) 
      elseif (present(fd_in) ) then
         fd = fd_in
      else
         error stop 'TruncationError_ReadFromFile: fName of fd_in must be provided'
      end if
      
      read(fd,*) Pxyz
      read(fd,*) this % Dir(1) % maxTE(NMINest:Pxyz(1)-1)
      read(fd,*) this % Dir(2) % maxTE(NMINest:Pxyz(2)-1)
      read(fd,*) this % Dir(3) % maxTE(NMINest:Pxyz(3)-1)
      
      if (present(fName) ) then
         close(fd)
      end if
      
      this % Dir(1) % P = Pxyz(1)
      this % Dir(2) % P = Pxyz(2)
      this % Dir(3) % P = Pxyz(3)
      
   end subroutine TruncationError_ReadFromFile
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine that extrapolates the behavior of the directional components
!  of the truncation error.
!  -----------------------------------------------------------------------
   subroutine TruncationError_ExtrapolateInOneDir(this,P_1,NMax,Dir,notenough,error)
      implicit none
      !---------------------------------------
      class(TruncationError_t), intent(inout)   :: this
      integer                    :: P_1               !<  P-1 (max polynomial order with tau estimation for regression)
      integer                    :: NMax
      integer                    :: Dir
      logical                    :: notenough         !>  .TRUE. if there are not enough points in every direction for regression 
      integer                    :: error             !>  error=1 if line behavior is not as expected
      !---------------------------------------
      
      real(kind=RP)              :: x   (P_1-NMINest+1)
      real(kind=RP)              :: y   (P_1-NMINest+1)
      integer                    :: N
      integer                    :: i
      real(kind=RP)              :: C,eta,r             ! Regression variables
      integer                    :: fd
      !---------------------------------------
      
      ! Initializations
      error     = 0
      notenough = .FALSE.
      
      ! Check if there are enough points for regression
      if (P_1 < NMINest + 1) then
         notenough = .TRUE.
         return
      end if
      
      ! Check if last point is NMIN=2
      if ( AlmostEqual(this % Dir(Dir) % maxTE (P_1),0._RP) ) then
         return ! nothing to do here
      end if
      
      ! Perform regression analysis   
      N = P_1 - NMINest + 1
      y = LOG10(this % Dir(Dir) % maxTE (NMINest:P_1))
      x = (/ (real(i,RP), i=NMINest,P_1) /)
      call LeastSquaresLinRegression(N,x,y,C,eta,r)
      
      ! Extrapolate the TE
      do i = P_1+1, NMax
         this % Dir(Dir) % maxTE(i) = 10**(C + eta*i)
      end do
      
      ! Check if there is an unexpected behavior
      if (eta >= 0) then
         Error = 1 
         return
      end if
      
   end subroutine TruncationError_ExtrapolateInOneDir
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  TruncationError_GenerateTEmap:
!  Routine to generate TEmap out of directional components
!  -----------------------------------------------------------------------
   subroutine TruncationError_GenerateTEmap(this,Nmax,TEmap)
      implicit none
      !-arguments-------------------------------------------------
      class(TruncationError_t), intent(in)    :: this
      integer                 , intent(in)    :: Nmax(NDIM)
      real(kind=RP)           , intent(inout) :: TEmap(NMINest:Nmax(1),NMINest:Nmax(2),NMINest:Nmax(3))
      !-local-variables-------------------------------------------
      integer :: i,j,k
      !-----------------------------------------------------------
      
      do k = NMINest, Nmax(3) ; do j = NMINest, Nmax(2)  ; do i = NMINest, Nmax(1)
         TEmap(i,j,k) = this % Dir(1) % maxTE(i) + &  !xi   contribution
                        this % Dir(2) % maxTE(j) + &  !eta  contribution
                        this % Dir(3) % maxTE(k)      !zeta contribution
      end do                  ; end do                   ; end do
      
   end subroutine TruncationError_GenerateTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------
!  Reset values of the truncation error
!  ------------------------------------
   subroutine TruncationError_reset(this)
      implicit none
      !-arguments-------------------------------------------------
      class(TruncationError_t), intent(inout)  :: this
      !-----------------------------------------------------------
      
      this % dir(1) % maxTE = 0._RP
      this % dir(2) % maxTE = 0._RP
      this % dir(3) % maxTE = 0._RP
   end subroutine TruncationError_reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Subroutine that sets P for all elements in TE
!  ---------------------------------------------
   subroutine InitializeForTauEstimation(TE,sem,TruncErrorType,TruncErrorForm,ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      implicit none
      !------------------------------------------
      type(TruncationError_t) :: TE(:)
      type(DGSem), intent(in) :: sem
      integer    , intent(in) :: TruncErrorType !<  Either NON_ISOLATED_TE or ISOLATED_TE
      integer    , intent(in) :: TruncErrorForm !<  Either NON_ISOLATED_TE or ISOLATED_TE
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      !------------------------------------------
      integer                 :: eID
      !------------------------------------------
      
!$omp parallel do schedule(runtime)
      do eID=1, size(sem % mesh % elements)
         TE(eID) % TruncErrorType = TruncErrorType
         TE(eID) % TruncErrorForm = TruncErrorForm
         TE(eID) % Dir(1) % P = sem % mesh % Nx(eID)
         TE(eID) % Dir(2) % P = sem % mesh % Ny(eID)
         TE(eID) % Dir(3) % P = sem % mesh % Nz(eID)
      end do
!$omp end parallel do
      
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
      type(TruncationError_t), intent(inout) :: TE(:)       !<> Type containing the truncation error estimation
      type(DGSem), target    , intent(inout) :: sem         !<> sem (to evaluate Qdot in a given coarser mesh)
      real(kind=RP)          , intent(in)    :: t           !<  time 
      type(MGSolStorage_t)   , intent(in)    :: Var(:)      !<  Type containing the source term in the mesh where the truncation error is being estimated (to be deprecated.. maintained just to keep consistency with manufactured solutions module)
      integer                , intent(in)    :: Dir         !<  Direction in which the truncation error is being estimated
      !--------------------------------------------------------
      integer                  :: iEl           !   Element counter
      integer                  :: iEQ           !   Equation counter
      integer                  :: i,j,k         !   Coordinate index counters
      real(kind=RP)            :: wx, wy, wz    ! 
      real(kind=RP)            :: Jac
      real(kind=RP)            :: maxTE
      real(kind=RP)            :: S(NCONS)      !   Source term
      type(Element), pointer   :: e
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))         
      procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
#endif
      !--------------------------------------------------------
      
      call TimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      
      S = 0._RP ! Initialize source term
      
!$omp parallel do private(iEl,maxTE,i,j,k,iEQ,wx,wy,wz,Jac,S,e) schedule(runtime)
      do iEl = 1, size(sem % mesh % elements)
         
         e => sem % mesh % elements(iEl)  ! An associate does not work well here
         
         if (e % Nxyz(Dir) >= TE(iEl) % Dir(Dir) % P) cycle ! it is actually never going to be greater than... but for security
         
         maxTE = 0._RP
         
         ! loop over all the degrees of freedom of the element
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))         
            call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, S, thermodynamics, dimensionless, refValues)
#endif
            wx  = NodalStorage(e % Nxyz(1)) % w (i)
            wy  = NodalStorage(e % Nxyz(2)) % w (j)
            wz  = NodalStorage(e % Nxyz(3)) % w (k)
            Jac = e % geom % jacobian(i,j,k)
            do iEQ = 1, NCONS
               if (TE(1) % TruncErrorForm .EQ. OLD_TE) then
                  maxTE =  MAX(maxTE , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,i,j,k) + S(iEQ) + Var(iEl) % Scase(iEQ,i,j,k) - Var(iEl) % S(iEQ,i,j,k) )  )  ! The last term is included to do time-accurate p-adaptation...For steady-state it can be neglected (original formulation)
               elseif (TE(1) % TruncErrorForm .EQ. NEW_TE) then
                  maxTE =  MAX(maxTE , (e % storage % Qdot (iEQ,i,j,k) + S(iEQ) + Var(iEl) % Scase(iEQ,i,j,k) - Var(iEl) % S(iEQ,i,j,k) )  ) 
               endif
            end do
         end do         ; end do         ; end do
         
         TE(iEl) % Dir(Dir) % maxTE(e % Nxyz(Dir)) = maxTE
         
      end do
!$omp end parallel do
   end subroutine EstimateTruncationError
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for printing the TE map(s) of one element
!  -----------------------------------------------------------------------
   subroutine PrintTEmap(NMIN,TEmap,iEl,FileName)
      implicit none
      !-------------------------------------------
      integer         , intent(in) :: NMIN(3)
      real(kind=RP)   , intent(in) :: TEmap(NMIN(1):,NMIN(2):,NMIN(3):)
      integer         , intent(in) :: iEl
      character(len=*), intent(in) :: FileName
      !-------------------------------------------
      integer                :: k, i, l
      integer                :: fd
      character(LINE_LENGTH) :: TEmapfile
      !-------------------------------------------
      
      do k = lbound(TEmap,3), ubound(TEmap,3)
         write(TEmapfile,'(3A,I7.7,A,I2.2,A)') 'RegressionFiles/TEmap_', trim(RemovePath(FileName)), '-Elem_',iEl,'-Nz_',k,'.dat'
   
         open(newunit = fd, file=TRIM(TEmapfile), action='write')
            do i = lbound(TEmap, 1), ubound(TEmap, 1)
               write(fd,*) (TEmap(i,l,k),l=lbound(TEmap,2),ubound(TEmap,2))
            end do
         close(fd)
      end do
   end subroutine PrintTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for printing the TE map(s) of one element
!  -----------------------------------------------------------------------
   subroutine AssignTimeDerivative(ComputeTimeDerivative)
      implicit none
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      
      TimeDerivative => ComputeTimeDerivative
   end subroutine AssignTimeDerivative
!
!  ------------------------------------------------------------------------
!  Subroutine for generating the exact  truncation error map in one element
!  -> The exact solution must be coded down in UserDefinedState1 (ProblemFile.f90)
!  ------------------------------------------------------------------------
!
   subroutine GenerateExactTEmap(sem, NMIN, NMAX, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, iEl, TruncErrorType, TruncErrorForm)
      implicit none
      !-arguments---------------------------------------------------------------
      type(DGSem)                :: sem
      integer, intent(in)        :: NMIN(NDIM)
      integer, intent(in)        :: NMAX(NDIM)
      real(kind=RP)              :: t
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables
      integer, intent(in)        :: iEl
      integer, intent(in)        :: TruncErrorType
      integer, intent(in)        :: TruncErrorForm

      !-external-routines------------------------------------------------------
#if defined(NAVIERSTOKES) && !defined(CAHNHILLIARD)
      procedure(UserDefinedFinalSetup_f) :: UserDefinedFinalSetup
#endif
      !-local variables---------------------------------------------------------
      real(kind=RP)              :: TEmap (NMIN(1):NMAX(1),NMIN(2):NMAX(2),NMIN(3):NMAX(3))
      integer                    :: i,j,k
      integer                    :: nelem      ! Number of elements
      logical                              :: success   
      integer, dimension(size(sem % mesh % elements) ) :: Nx, Ny, Nz
      !-------------------------------------------------------------------------
      
      ! Initializations
      nelem   = SIZE(sem % mesh % elements)
      
      if (TruncErrorType == ISOLATED_TE) then
         TimeDerivative => ComputeTimeDerivativeIsolated
      else
         TimeDerivative => ComputeTimeDerivative
      end if
      
      do k = NMIN(3), NMAX(3)
         do j = NMIN(2), NMAX(2)
            do i = NMIN(1), NMAX(1)
               
               Nx = i
               Ny = j
               Nz = k
               
               !-------------------------------------------
               !Destruct previous sem and construct new one
               !-------------------------------------------
               CALL sem % destruct()
               
               sem % ManufacturedSol = controlVariables % containsKey("manufactured solution")
               
               call sem % construct (  controlVariables  = controlVariables                              ,  &
                           Nx_ = Nx ,     Ny_ = Ny,     Nz_ = Nz,                   &
                           success           = success)
               
               print*, 'Computing TE for N=',i,j,k,'. success=', success
               
               if(.NOT. success)   error stop ":: problem creating sem"
               
#if defined(NAVIERSTOKES) && !defined(CAHNHILLIARD) && !defined(SPALARTALMARAS)
               CALL UserDefinedFinalSetup(sem % mesh , thermodynamics, dimensionless, refValues)
#endif
               
               
               TEmap(i,j,k) = EstimateTauOfElem(sem,t,TruncErrorForm,iEl)
               
               print*, 'Done for N=',i,j,k
            end do
         end do
      end do
      
      CALL PrintTEmap(NMIN,TEmap,iEl,"Exact")
      error stop
      
   end subroutine GenerateExactTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----
!     Estimates the infinity norm of the truncation error for a given mesh given a restart solution or a manufactured solution
!     ---
      function EstimateTauOfElem(sem,t,iEl,TruncErrorForm, UseLoadedSol) result(maxTE)
         implicit none
         !-arguments---------------------------------------------
         type(DGSem), target    , intent(inout) :: sem              !<> sem class (inout cause' we compute Qdot)
         real(kind=RP)          , intent(in)    :: t                !>  Time
         integer                , intent(in)    :: iEl              !<  Present if the result is wanted for a certain element
         real(kind=RP)                          :: maxTE            !>  |\tau|_{\infty}
         integer                , intent(in)    :: TruncErrorForm
         logical      , optional, intent(in)    :: UseLoadedSol

         !-external-routines------------------------------------
#if defined(NAVIERSTOKES)  
         procedure(UserDefinedState_f) :: UserDefinedState1
#endif
         !-local-variables--------------------------------------
         integer          :: nelem
         integer          :: eID             ! element counter
         integer          :: iEQ             ! Equation counter
         integer          :: ii,jj,kk        ! doF counters
         real(kind=RP)    :: wx, wy, wz      ! Quadrature weights
         real(kind=RP)    :: Jac             ! Jacobian (mapping)
         logical          :: UserDefinedSol
         !-------------------------------------------------------
         
         nelem = SIZE(sem % mesh % elements)
         
         UserDefinedSol = .TRUE.
         if ( present(UseLoadedSol) ) then
            if (UseLoadedSol) UserDefinedSol = .FALSE.
         end if
         
         !-------------------------------------------
         !Get exact solution (from ProblemFile.f90)
         !-------------------------------------------
         
         do eID = 1, nelem
            associate (e => sem % mesh % elements(eID))
            
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
            if (UserDefinedSol) then
               do kk = 0, e % Nxyz(3) ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
                  call UserDefinedState1(e % geom % x(:,ii,jj,kk), t, [0._RP, 0._RP, 0._RP], e % storage % Q(:,ii,jj,kk), thermodynamics, dimensionless, refValues)
               end do                 ; end do                 ; end do
            end if
#endif
            
            end associate
         end do
         
         call TimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
         
         maxTE = 0._RP ! Initialization
         
         associate ( e => sem % mesh % elements(iEl) )
         
         ! loop over all the degrees of freedom of the element
         do kk = 0, e % Nxyz(3) ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            
            do iEQ = 1, NCONS
               wx  = NodalStorage(e % Nxyz(1)) % w (ii)
               wy  = NodalStorage(e % Nxyz(2)) % w (jj)
               wz  = NodalStorage(e % Nxyz(3)) % w (kk)
               Jac = e % geom % jacobian(ii,jj,kk)
            
              if (TruncErrorForm .EQ. OLD_TE) then 
                  maxTE =  MAX(maxTE , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,ii,jj,kk) ) )
              elseif (TruncErrorForm .EQ. NEW_TE) then
                  maxTE =  MAX(maxTE , ABS  (e % storage % Qdot (iEQ,ii,jj,kk) ) )
              endif
            
            end do
         end do          ; end do          ; end do
          
         end associate
      end function EstimateTauOfElem
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  This subroutine estimates and plot the truncation error of a representation
!  -----------------------------------------------------------------------    
   subroutine EstimateAndPlotTruncationError(sem, t, controlVariables, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      implicit none
      !-arguments------------------------------------------------
      type(DGSem), target    , intent(inout) :: sem              !<> 
      real(kind=RP)          , intent(in)    :: t
      type(FTValueDictionary), intent(in)    :: controlVariables
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      !-local-variables------------------------------------------
      integer       :: eID, ii, jj, kk, iEQ
      real(kind=RP) :: Tau    (sem % mesh % no_of_elements)
      real(kind=RP) :: IsolTau(sem % mesh % no_of_elements), wx, wy, wz, Jac
      type(HexMesh) :: auxMesh
      integer       :: NDOF
      !----------------------------------------------------------
      print*, 'Estimating truncerror'
      
      Tau = 0._RP
      IsolTau = 0._RP
!
!     Non-isolated truncation error
!     -----------------------------
      
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      do eID = 1, sem % mesh % no_of_elements
         associate (e => sem % mesh % elements(eID))
         
         do kk = 0, e % Nxyz(3) ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            
            wx  = NodalStorage(e % Nxyz(1)) % w (ii)
            wy  = NodalStorage(e % Nxyz(2)) % w (jj)
            wz  = NodalStorage(e % Nxyz(3)) % w (kk)
            Jac = e % geom % jacobian(ii,jj,kk)
            
            do iEQ = 1, NCONS
               Tau(eID) =  MAX(Tau(eID) , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,ii,jj,kk) ) )
            end do
         end do          ; end do          ; end do
          
         end associate
      end do
      
!
!     Isolated truncation error
!     -------------------------
      
      call ComputeTimeDerivativeIsolated(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      do eID = 1, sem % mesh % no_of_elements
         associate (e => sem % mesh % elements(eID) )
         
         do kk = 0, e % Nxyz(3) ; do jj = 0, e % Nxyz(2) ; do ii = 0, e % Nxyz(1)
            
            wx  = NodalStorage(e % Nxyz(1)) % w (ii)
            wy  = NodalStorage(e % Nxyz(2)) % w (jj)
            wz  = NodalStorage(e % Nxyz(3)) % w (kk)
            Jac = e % geom % jacobian(ii,jj,kk)
            
            do iEQ = 1, NCONS
               IsolTau(eID) =  MAX(IsolTau(eID) , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,ii,jj,kk) ) )
            end do
         end do          ; end do          ; end do
          
         end associate
      end do
      
!     Create auxMesh, save the truncation errors there as Q(1) and Q(2), and export it
!     ---------------------------------------------------------------------------------
      
      ! 1. Prepare mesh
      auxMesh % nodeType = sem % mesh % nodeType
      auxMesh % no_of_elements = sem % mesh % no_of_elements
      allocate ( auxMesh % elements (sem % mesh % no_of_elements) )
      
      NDOF = 0
      do eID = 1, sem % mesh % no_of_elements
         associate ( e_aux => auxMesh % elements(eID), &
                     e     =>    sem % mesh % elements(eID) )
         e_aux % globID = e % globID
         e_aux % Nxyz = e % Nxyz
         NDOF = NDOF + product(e % Nxyz + 1)
         end associate
      end do
      
      call auxMesh % PrepareForIO
      call auxMesh % AllocateStorage (NDOF, controlVariables,.FALSE.)
      
      ! 2. Save the error there
      do eID = 1, sem % mesh % no_of_elements
         associate (e => auxMesh % elements(eID) )
         e % Storage % Q(1,:,:,:) = Tau(eID)
         e % Storage % Q(2,:,:,:) = IsolTau(eID)
         end associate
      end do
      
!~      print*, 'Plotting truncerror'
!~      call PlotTruncationError(sem % mesh, Tau, IsolTau)
      
      call auxMesh % SaveSolution(sem % numberOfTimeSteps,t,'RESULTS/TruncError.hsol',.FALSE.)
      
   end subroutine EstimateAndPlotTruncationError
      
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Subroutine for plotting the adaptation information:
!  Plots the truncation error contained in array Tau
!  -----------------------------------------------------------------------
   subroutine PlotTruncationError(mesh,Tau,IsolTau)
      use HexMeshClass, only: HexMesh
      implicit none
      !--------------------------------------
      TYPE(HexMesh)          , intent(in) :: mesh       !<  mesh...
      real(kind=RP)          , intent(in) :: Tau(:)     !<  Type containing the element-wise truncation error estimation
      real(kind=RP)          , intent(in) :: IsolTau(:) !<  Type containing the element-wise truncation error estimation
      !--------------------------------------
      integer                            :: fd
      integer                            :: eID
      integer                            :: i,j,k
      CHARACTER(len=LINE_LENGTH)         :: plotFileName
      !--------------------------------------
      
      write(plotFileName,'(A)') 'RESULTS/TruncationError.tec'
      
      open(newunit = fd, file=plotFileName, action='WRITE')
         
         write(fd,*) 'TITLE = "Truncation error information (HORSES3D)"'
         write(fd,*) 'VARIABLES = "X","Y","Z","Tau","IsolTau"'
         
         do eID = 1, size(mesh % elements)
            associate (N => mesh % elements(eID) % Nxyz)
            write(fd,*) 'ZONE I=', N(1)+1, ",J=",N(2)+1, ",K=",N(3)+1,", F=POINT"
            
            !-------
            ! Plot!
            !-------
            do k = 0, N(3) ; do j= 0, N(2) ; do i = 0, N(1)
                     write(fd,'(5E13.5)') &
                                mesh % elements(eID) % geom % x(1,i,j,k), &
                                mesh % elements(eID) % geom % x(2,i,j,k), &
                                mesh % elements(eID) % geom % x(3,i,j,k), &
                                Tau(eID), IsolTau(eID)
            end do ; end do ; end do
            
            end associate
            
         end do
      
      close(fd)
      
   end subroutine PlotTruncationError
end module TruncationErrorClass