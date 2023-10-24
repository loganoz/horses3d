!
!//////////////////////////////////////////////////////
!
!  Module that provides routines to estimate the truncation error every n number of 
!  interactions and store it in files.
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module MultiTauEstimationClass
   use SMConstants
   use TruncationErrorClass   , only: TruncationError_t, NON_ISOLATED_TE, ISOLATED_TE, NMINest
   use AnisFASMultigridClass  , only: AnisFASMultigrid_t
   use FTValueDictionaryClass , only: FTValueDictionary
   use DGSEMClass             , only: DGSem, ComputeTimeDerivative_f
   use ParamfileRegions       , only: readValueInRegion, getSquashedLine
   use Utilities              , only: toLower
   use MPI_Process_Info       , only: MPI_Process
   use Headers                , only: Subsection_Header
   implicit none
   
   private
   public MultiTauEstim_t
   
!
!  Main type for multi tau-estimation
!  ----------------------------------
   type MultiTauEstim_t
      logical, private                    :: Active = .FALSE.
      integer                             :: interval
      integer                             :: TruncErrorType
      integer                             :: stage
      integer                             :: num_of_elements
      character(len=LINE_LENGTH)          :: FolderName
      type(TruncationError_t),allocatable :: TE(:)
      type(AnisFASMultigrid_t)            :: AnisFAS
      
      contains
         procedure :: construct              => MultiTau_Construct
         procedure :: constructForPAdaptator => MultiTau_ConstructForPAdaptator
         procedure :: destruct               => MultiTau_Destruct
         procedure :: describe               => MultiTau_Describe
         procedure :: estimate               => MultiTau_Estimate
         procedure :: GetTEmap               => MultiTau_GetTEmap
   end type MultiTauEstim_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MultiTau_Construct(this, controlVariables, sem)
         implicit none
         !-arguments----------------------------------------------------
         class(MultiTauEstim_t)  , intent(inout)   :: this
         type(FTValueDictionary) , intent(in)      :: controlVariables
         type(DGSem)             , intent(in)      :: sem
         !-local-variables----------------------------------------------
         integer :: eID
         character(LINE_LENGTH) :: paramFile
         character(LINE_LENGTH) :: in_label
         character(LINE_LENGTH) :: R_TruncErrorType, R_Folder
         integer, allocatable   :: R_interval
         !--------------------------------------------------------------
!
!        Initialize
!        **********
         this % Active = MultiTauEstimationIsDefined()
         if (.not. this % Active) return
!
!        Read the input
!        **************
         write(in_label , '(A)') "#define multi tau-estimation"
      
         call get_command_argument(1, paramFile) !
      
         call readValueInRegion ( trim ( paramFile )  , "truncation error type"  , R_TruncErrorType     , in_label , "# end" )
         call readValueInRegion ( trim ( paramFile )  , "interval"               , R_interval           , in_label , "# end" )
         call readValueInRegion ( trim ( paramFile )  , "folder"                 , R_Folder             , in_label , "# end" )
         
!        Truncation error type
!        ---------------------
         call toLower (R_TruncErrorType)
         select case ( trim(R_TruncErrorType) )
            case ('isolated')
               this % TruncErrorType = ISOLATED_TE
            case ('non-isolated')
               this % TruncErrorType = NON_ISOLATED_TE
            case default
               this % TruncErrorType = NON_ISOLATED_TE
         end select
         
!        Estimation interval
!        -------------------
         if ( allocated(R_interval) ) then
            this % interval = R_interval
         else
            error stop 'Multi tau-estimation :: An interval must be specified'
         end if
         
!        Output folder
!        -------------
         if ( R_Folder == "") then
            this % FolderName = "MultiTau"
         else
            this % FolderName = trim(R_Folder)
         end if
!
!        Allocate memory to store the truncation error
!        *********************************************
         this % num_of_elements = sem % mesh % no_of_elements
         allocate (this % TE(this % num_of_elements))
         do eID = 1, this % num_of_elements
            call this % TE(eID) % construct(sem % mesh % elements(eID) % Nxyz-1)
            this % TE(eID) % Dir(1) % P = sem % mesh % elements(eID) % Nxyz(1) - 1
            this % TE(eID) % Dir(2) % P = sem % mesh % elements(eID) % Nxyz(2) - 1
            this % TE(eID) % Dir(3) % P = sem % mesh % elements(eID) % Nxyz(3) - 1
         end do
         this % stage = 0
!
!        Construct anisotropic FAS for error estimation
!        **********************************************
         call this % AnisFAS % construct(controlVariables,sem,estimator=.TRUE.,NMINestim = NMINest)
!
!        Create folder for storing the truncation error
!        **********************************************
         call execute_command_line ('mkdir -p ' // this % FolderName)
         
         call this % describe
         
      end subroutine MultiTau_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MultiTau_ConstructForPAdaptator(this, TruncErrorType, NxyzMax, num_of_elements, Folder)
         implicit none
         !-arguments----------------------------------------------------
         class(MultiTauEstim_t)  , intent(inout)   :: this
         integer                 , intent(in)      :: TruncErrorType
         integer                 , intent(in)      :: NxyzMax(NDIM)
         integer                 , intent(in)      :: num_of_elements
         character(LINE_LENGTH)  , intent(in)      :: Folder
         !-local-variables----------------------------------------------
         integer :: eID
         !--------------------------------------------------------------
!
!        Initialize
!        **********
         this % Active = .FALSE.
         
!        Truncation error type
!        ---------------------
         this % TruncErrorType = TruncErrorType
         
!        Output folder
!        -------------
         this % FolderName = Folder
!
!        Allocate memory to store the truncation error
!        *********************************************
         this % num_of_elements = num_of_elements
         this % stage = 0
         
!~          call this % describe
         
      end subroutine MultiTau_ConstructForPAdaptator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MultiTau_Destruct(this)
         implicit none
         !-arguments----------------------------------------------------
         class(MultiTauEstim_t), intent(inout) :: this
         !--------------------------------------------------------------
         
         if (.not. this % Active) return
         
!        Destruct truncation error storage
!        ---------------------------------
         call this % TE % destruct
         deallocate (this % TE)
         
!        Destruct multigrid
!        ------------------
         call this % AnisFAS % destruct
         
      end subroutine MultiTau_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MultiTau_Describe(this)
         implicit none
         !-arguments----------------------------------------------------
         class(MultiTauEstim_t), intent(inout) :: this
         !--------------------------------------------------------------
         
         if (.not. MPI_Process % isRoot) return
         
         write(STD_OUT,'(/)')
         call Subsection_Header("Multi tau-estimation")
         
         write(STD_OUT,'(30X,A,A30,I6)') "->","Iteration interval: ", this % interval
         select case (this % TruncErrorType)
            case (NON_ISOLATED_TE)  ; write(STD_OUT,'(30X,A,A30,A)') "->", "Truncation error type: ", "non-isolated"
            case (ISOLATED_TE)      ; write(STD_OUT,'(30X,A,A30,A)') "->", "Truncation error type: ", "isolated"
         end select
         write(STD_OUT,'(30X,A,A30,A)') "->", "Storing in folder: ", this % FolderName
         
      end subroutine MultiTau_Describe
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MultiTau_Estimate(this, sem, iter, t, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
         implicit none
         !-arguments----------------------------------------------------
         class(MultiTauEstim_t)  , intent(inout)   :: this
         type(DGSem)             , intent(inout)   :: sem
         integer                 , intent(in)      :: iter
         real(kind=RP)           , intent(in)      :: t
         procedure(ComputeTimeDerivative_f)        :: ComputeTimeDerivative
         procedure(ComputeTimeDerivative_f)        :: ComputeTimeDerivativeIsolated
         !-local-variables----------------------------------------------
         integer :: eID
         character(len=LINE_LENGTH) :: fName
         !--------------------------------------------------------------
         
!        Initial checks
!        **************
         if (.not. this % Active) return
         if (iter == 1 .or. mod(iter,this % interval) /= 0) return
         
         this % stage = this % stage + 1
!
!        Estimate the truncation error
!        *****************************
         call this % AnisFAS % solve(iter,t,computeTimeDerivative,ComputeTimeDerivativeIsolated,this % TE, this % TruncErrorType)
!
!        Write the estimation to files
!        *****************************
         do eID = 1, this % num_of_elements
            write(fName,'(A,A,I8.8,A)') trim(this % FolderName), "/Elem", sem % mesh % elements(eID) % globID, ".dat"
            call this % TE(eID) % ExportToFile(fName)
         end do
!
!        Reset the truncation error
!        **************************
         do eID = 1, this % num_of_elements
            
            this % TE(eID) % Dir(1) % maxTE = 0._RP
            this % TE(eID) % Dir(2) % maxTE = 0._RP
            this % TE(eID) % Dir(3) % maxTE = 0._RP
            
         end do
         
      end subroutine MultiTau_Estimate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------
!     Get maximum values for the TEmap of an element
!     ----------------------------------------------
      subroutine MultiTau_GetTEmap(this,stages,glob_eID,Nmax,Nmin,P_1,TEmap)
         implicit none
         !-arguments----------------------------------------------------
         class(MultiTauEstim_t)  , intent(in)      :: this
         integer                 , intent(in)      :: stages(2) ! Initial and final estimation stage
         integer                 , intent(in)      :: glob_eID
         integer                 , intent(in)      :: Nmax(NDIM)  ! Maximum polynomial order for adaptation
         integer                 , intent(in)      :: Nmin(NDIM)  ! Minimum polynomial order for adaptation
         integer                 , intent(inout)   :: P_1 (NDIM)
         real(kind=RP)           , intent(out)     :: TEmap(NMINest:Nmax(1),NMINest:Nmax(2),NMINest:Nmax(3))
         !-local-variables----------------------------------------------
         integer        :: stage, fd
         integer        :: i, j, k, dir, error
         logical        :: notenough
         real(kind=RP)  :: TEmap_temp(NMINest:Nmax(1),NMINest:Nmax(2),NMINest:Nmax(3))
         character(len=LINE_LENGTH) :: fName
         type(TruncationError_t)    :: TE
         !--------------------------------------------------------------
         
         write(fName,'(A,A,I8.8,A)') trim(this % FolderName), "/Elem", glob_eID, ".dat"
         open (newunit=fd, file = trim(fName)) 
         
         call TE % construct( Nmax )
         
         TEmap = 0._RP
         TEmap_temp = 0._RP
         do stage=stages(1), stages(2)
            call TE % reset
            call TE % ReadFromFile(fd_in = fd)
            
            do dir=1, 3
               ! Adjust P-1
               P_1(dir)  = TE % Dir(dir) % P - 1
               if ( (P_1(dir) < NMIN(dir)) .and. (P_1(dir) < NMINest+1) ) P_1(dir) = NMIN(dir)
               
               ! Extrapolate
               call TE % ExtrapolateInOneDir (P_1(dir),NMax(dir),Dir,notenough,error)
               if (notenough .or. error==1) TE % Dir(dir) % maxTE(P_1(dir)+1:Nmax(dir)) = huge(1._RP) / 4 ! Divide by 4 to prevent overflow
            end do
            
            call TE % GenerateTEmap(Nmax,TEmap_temp)
         
            do k = NMINest, Nmax(3) ; do j = NMINest, Nmax(2)  ; do i = NMINest, Nmax(1)
               if (TEmap_temp(i,j,k) > TEmap(i,j,k)) TEmap(i,j,k) = TEmap_temp(i,j,k)
            end do                  ; end do                   ; end do
         end do
         
         call TE % destruct
         close (fd)
         
      end subroutine MultiTau_GetTEmap
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Auxiliary routines
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   logical function MultiTauEstimationIsDefined()
      implicit none
      !-local-variables----------------------------------------------
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
      !--------------------------------------------------------------
!
!     Initialize
!     ----------
      MultiTauEstimationIsDefined = .FALSE.
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)
!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )
!
!     Read the whole file to find overenriching regions
!     -------------------------------------------------
readloop:do 
         read ( fID , '(A)' , iostat = io ) line

         if ( io .lt. 0 ) then
!
!           End of file
!           -----------
            line = ""
            exit readloop

         elseif ( io .gt. 0 ) then
!
!           Error
!           -----
            errorMessage(STD_OUT)
            error stop "Stopped."
         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#definemultitau-estimation' ) .gt. 0 ) then
               MultiTauEstimationIsDefined = .TRUE.
               close(fID)  
               return
            end if
         end if
      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

   end function MultiTauEstimationIsDefined

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MultiTauEstimationClass