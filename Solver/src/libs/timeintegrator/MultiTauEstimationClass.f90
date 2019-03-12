!
!//////////////////////////////////////////////////////
!
!   @File:    MultiTauEstimationClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Mar 12 15:43:41 2019
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
!  Module that provides routines to estimate the truncation error every n number of 
!  iteractions and store it in files.
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module MultiTauEstimationClass
   use SMConstants
   use TruncationErrorClass   , only: TruncationError_t, NON_ISOLATED_TE, ISOLATED_TE
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
      logical, private                    :: Active
      integer                             :: interval
      integer                             :: TruncErrorType
      integer                             :: stage
      integer                             :: num_of_elements
      character(len=LINE_LENGTH)          :: FolderName
      type(TruncationError_t),allocatable :: TE(:)
      type(AnisFASMultigrid_t)            :: AnisFAS
      
      contains
         procedure :: construct  => MultiTau_Construct
         procedure :: destruct   => MultiTau_Destruct
         procedure :: describe   => MultiTau_Describe
         procedure :: estimate   => MultiTau_Estimate
   end type MultiTauEstim_t
!
!  Parameters
!  ----------
   integer, parameter :: NMINest = 1
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
            call this % TE(eID) % construct(NMINest,sem % mesh % elements(eID) % Nxyz-1)
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
         integer :: eID, fd
         character(len=LINE_LENGTH) :: CurrentFolder
         character(len=LINE_LENGTH) :: fName
         !--------------------------------------------------------------
         
!        Initial checks
!        **************
         if (.not. this % Active) return
         if (mod(iter,this % interval) /= 0) return
         
         this % stage = this % stage + 1
!
!        Estimate the truncation error
!        *****************************
         call this % AnisFAS % solve(iter,t,computeTimeDerivative,ComputeTimeDerivativeIsolated,this % TE, this % TruncErrorType)
!
!        Write the estimation to files
!        *****************************
!        Create new folder
!        -----------------
         write(CurrentFolder,'(A,A,I5.5)') trim(this % FolderName), "/", this % stage
         call execute_command_line('mkdir -p '// CurrentFolder)

!        Write general info file
!        -----------------------
         if (MPI_Process % isRoot) then
            open (newunit=fd, file = trim(CurrentFolder) // "/info")
            
            write(fd,'(A24,ES24.16)') 'Time ='                 , t
            write(fd,'(A24,I24)')     'Iteration ='            , iter
            select case (this % TruncErrorType)
               case (ISOLATED_TE)     ; write(fd,'(A24,A24)')     'Truncation error type =', 'isolated'
               case (NON_ISOLATED_TE) ; write(fd,'(A24,A24)')     'Truncation error type =', 'non-isolated'
            end select
            write(fd,'(A24,I24)')     'NMIN      ='            , NMINest
            
            close (fd)
         end if
         
!        Write truncation error
!        ----------------------
         
         do eID = 1, this % num_of_elements
            write(fName,'(A,A,I8.8,A)') trim(CurrentFolder), "/Elem", sem % mesh % elements(eID) % globID, ".dat"
            open (newunit=fd, file = trim(fName)) 
            
            write(fd,*) sem % mesh % elements(eID) % Nxyz
            write(fd,*) this % TE(eID) % Dir(1) % maxTE
            write(fd,*) this % TE(eID) % Dir(2) % maxTE
            write(fd,*) this % TE(eID) % Dir(3) % maxTE
            
            close (fd)
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
!  Auxiliar routines
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
            stop "Stopped."
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
