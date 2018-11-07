!
!//////////////////////////////////////////////////////
!
!   @File:    TruncationErrorClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Feb 28 14:00:00 2018
!   @Last revision date: Tue Oct 30 18:24:03 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: e995d28c5dd9fe9c7286ec6c7e53405ab11a7c14
!
!//////////////////////////////////////////////////////
!
module TruncationErrorClass
   use SMConstants
   use MultigridTypes            , only: MGSolStorage_t
   use DGSEMClass                , only: DGSem, ComputeTimeDerivative_f
   use FTValueDictionaryClass    , only: FTValueDictionary
   use PhysicsStorage            , only: NTOTALVARS, CTD_IGNORE_MODE
   use HexMeshClass              , only: HexMesh
#if defined(NAVIERSTOKES)   || defined(INCNS)
   use FluidData_NS              , only: Thermodynamics, RefValues, Dimensionless
   use ProblemFileFunctions      , only: UserDefinedState_f, UserDefinedFinalSetup_f, UserDefinedSourceTermNS_f
#endif
   use NodalStorageClass         , only: NodalStorage
   use FileReadingUtilities      , only: RemovePath
   implicit none
   
   private
   public TruncationError_t, EstimateTruncationError, InitializeForTauEstimation, PrintTEmap, AssignTimeDerivative, GenerateExactTEmap, EstimateAndPlotTruncationError
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
!  ----------------
!  Module variables
!  ----------------
!
   procedure(ComputeTimeDerivative_f), pointer :: TimeDerivative
   
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
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      !------------------------------------------
      integer                 :: eID
      !------------------------------------------
      
      do eID=1, size(sem % mesh % elements)
         TE(eID) % TruncErrorType = TruncErrorType
         TE(eID) % Dir(1) % P = sem % mesh % Nx(eID)
         TE(eID) % Dir(2) % P = sem % mesh % Ny(eID)
         TE(eID) % Dir(3) % P = sem % mesh % Nz(eID)
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
#if defined(NAVIERSTOKES)            
      procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
#endif
      !--------------------------------------------------------
      
      call TimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      
      S = 0._RP ! Initialize source term
      
      do iEl = 1, size(sem % mesh % elements)
         
         associate (e => sem % mesh % elements(iEl))
         
         if (e % Nxyz(Dir) >= TE(iEl) % Dir(Dir) % P) cycle ! it is actually never going to be greater than... but for security
         
         maxTE = 0._RP
         
         ! loop over all the degrees of freedom of the element
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
#if defined(NAVIERSTOKES)            
            call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, S, thermodynamics, dimensionless, refValues)
#endif
            
            do iEQ = 1, NTOTALVARS
               wx  = NodalStorage(e % Nxyz(1)) % w (i)
               wy  = NodalStorage(e % Nxyz(2)) % w (j)
               wz  = NodalStorage(e % Nxyz(3)) % w (k)
               Jac = e % geom % jacobian(i,j,k)
               
               maxTE =  MAX(maxTE , wx * wy * wz * Jac * ABS  (e % storage % Qdot (iEQ,i,j,k) + S(iEQ) + Var(iEl) % Scase(iEQ,i,j,k) - Var(iEl) % S(iEQ,i,j,k) )  )  ! The last term is included to do time-accurate p-adaptation...For steady-state it can be neglected (original formulation)
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
   subroutine GenerateExactTEmap(sem, NMIN, NMAX, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, iEl, TruncErrorType)
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
               
               if(.NOT. success)   ERROR STOP ":: problem creating sem"
               
#if defined(NAVIERSTOKES) && !defined(CAHNHILLIARD)
               CALL UserDefinedFinalSetup(sem % mesh , thermodynamics, dimensionless, refValues)
#endif
               
               
               TEmap(i,j,k) = EstimateTauOfElem(sem,t,iEl)
               
               print*, 'Done for N=',i,j,k
            end do
         end do
      end do
      
      CALL PrintTEmap(NMIN,TEmap,iEl,"Exact")
      stop
      
   end subroutine GenerateExactTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----
!     Estimates the infinity norm of the truncation error for a given mesh given a restart solution or a manufactured solution
!     ---
      function EstimateTauOfElem(sem,t,iEl, UseLoadedSol) result(maxTE)
         implicit none
         !-arguments---------------------------------------------
         type(DGSem), target    , intent(inout) :: sem              !<> sem class (inout cause' we compute Qdot)
         real(kind=RP)          , intent(in)    :: t                !>  Time
         integer                , intent(in)    :: iEl              !<  Present if the result is wanted for a certain element
         real(kind=RP)                          :: maxTE            !>  |\tau|_{\infty}
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
            
#if defined(NAVIERSTOKES)      
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
            
            do iEQ = 1, NTOTALVARS
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
            
            do iEQ = 1, NTOTALVARS
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
      call auxMesh % AllocateStorage (NDOF, controlVariables,.FALSE.,.FALSE.)
      
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
