module StatsAveragingModule
   use SMConstants
   use SolutionFile
   use StatisticsMonitor, only: NO_OF_VARIABLES
   implicit none
   
   private
   public Initialize_StatsAveragingModule, Finalize_StatsAveragingModule, PerformAveraging
   
   type Element_t
!                                /* Solution quantities */
      integer                    :: Nsol(NDIM)
      real(kind=RP), pointer     :: stats(:,:,:,:)
   end type Element_t
   
   type Mesh_t
      integer                          :: num_of_elements
      integer                          :: nodeType
      integer                          :: iter
      real(kind=RP)                    :: time
      real(kind=RP)                    :: refs(NO_OF_SAVED_REFS)
      type(Element_t),   allocatable   :: elements(:)
      character(len=LINE_LENGTH)       :: meshName
      character(len=LINE_LENGTH)       :: solutionName
      contains
         procedure   :: construct               => Mesh_Construct
         procedure   :: destruct                => Mesh_Destruct
         procedure   :: ReadStatistics          => Mesh_ReadStatistics
         procedure   :: SaveStatistics          => Mesh_SaveStatistics
         procedure   :: AddWeightedContribution => Mesh_AddWeightedContribution
   end type Mesh_t
   
   type(Mesh_t)   :: AveragedSol, CurrentSol
   
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_StatsAveragingModule(num_of_statFiles,fileNames)
         implicit none
         !-arguments--------------------------------------------------
         integer                    , intent(in) :: num_of_statFiles
         character(len=LINE_LENGTH) , intent(in) :: filenames(num_of_statFiles)
         !-local-variables--------------------------------------------
         integer                    :: eID
         integer                    :: fid
         integer                    :: num_of_elements
         integer                    :: arrayDimensions(4)
         integer, allocatable       :: Nsol(:,:)
         character(len=LINE_LENGTH) :: solutionName
         real(kind=RP), allocatable :: tempStats(:,:,:,:)
         !------------------------------------------------------------
         
!        Get information from last file 
!           (the rest must be of the same size)
!        --------------------------------------
         solutionName = fileNames(num_of_statFiles)
         num_of_elements = getSolutionFileNoOfElements( solutionName )
         
         allocate ( Nsol(3,num_of_elements) )
         
         
         fid = putSolutionFileInReadDataMode(solutionName)
         do eID = 1, num_of_elements
            call getSolutionFileArrayDimensions(fid,arrayDimensions)
            Nsol(1:3,eID) = arrayDimensions(2:4) - 1
            allocate (tempStats (1:NO_OF_VARIABLES,0:Nsol(1,eID),0:Nsol(3,eID),0:Nsol(2,eID)) )
            read(fid) tempStats
            deallocate (tempStats)
         end do
         close(fid)
         
         call AveragedSol % construct (num_of_elements, Nsol)
         call CurrentSol  % construct (num_of_elements, Nsol)
         
!
!        Get time and iteration
!        ----------------------
         call getSolutionFileTimeAndIteration( trim( solutionName ),AveragedSol % iter, AveragedSol % time)
!
!        Read reference values
!        ---------------------
         AveragedSol % refs = getSolutionFileReferenceValues( trim( solutionName ) )
!
!        Get node type
!        -------------
         AveragedSol % nodeType = getSolutionFileNodeType(solutionName)
         
         
      end subroutine Initialize_StatsAveragingModule
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine PerformAveraging(num_of_statFiles,fileNames,weights,fileName)
         use Headers
         implicit none
         !-arguments--------------------------------------------------
         integer                    , intent(in) :: num_of_statFiles
         character(len=LINE_LENGTH) , intent(in) :: filenames  (num_of_statFiles)
         real(kind=RP)              , intent(in) :: weights    (num_of_statFiles)
         character(len=LINE_LENGTH) , intent(in) :: fileName
         !-local-variables--------------------------------------------
         integer :: fileNum
         !------------------------------------------------------------
         
         write(STD_OUT,'(/,/)')
         call Section_Header("Averaging files")
         write(STD_OUT,'(/,/)')
         
         do fileNum = 1, num_of_statFiles
            call CurrentSol  % ReadStatistics(fileNames(fileNum))
            call AveragedSol % AddWeightedContribution (CurrentSol, weights(fileNum))
         end do
         
         write(STD_OUT,'(A20,A,A)') 'Done.', ' Writing merged statistics to: ', trim(fileName)
         
         call AveragedSol % SaveStatistics( trim(fileName) )
         write(STD_OUT,'(/,/)')
         write(STD_OUT,'(/,/)')
      end subroutine PerformAveraging
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_StatsAveragingModule
         implicit none
         
         call AveragedSol % destruct
         call CurrentSol % destruct
         
      end subroutine Finalize_StatsAveragingModule
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Mesh_Construct(this, num_of_elements, Nsol)
         implicit none
         !-arguments---------------------------------------------------
         class(Mesh_t)  , intent(inout)   :: this
         integer        , intent(in)      :: num_of_elements
         integer        , intent(in)      :: Nsol(3,num_of_elements)
         !-local-variables---------------------------------------------
         integer :: eID
         !-------------------------------------------------------------
         
         this % num_of_elements = num_of_elements
         allocate ( this % elements(num_of_elements) )
         
         do eID = 1, this % num_of_elements
            associate (e => this % elements(eID) )
            allocate( e % stats(1:NO_OF_VARIABLES,0:Nsol(1,eID), 0:Nsol(2,eID), 0:Nsol(3,eID)) )
            e % stats = 0._RP
            end associate
         end do
      end subroutine Mesh_Construct
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Mesh_Destruct(this)
         implicit none
         !-arguments---------------------------------------------------
         class(Mesh_t)  , intent(inout)   :: this
         !-------------------------------------------------------------
         
         deallocate (this % elements)
         
      end subroutine Mesh_Destruct
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Mesh_AddWeightedContribution(this,contribution,weight)
         implicit none
         !-arguments---------------------------------------------------
         class(Mesh_t)         :: this
         type(Mesh_t)          :: contribution
         real(kind=RP)         :: weight
         !-local-variables---------------------------------------------
         integer       :: eID         
         !-------------------------------------------------------------
         
!$omp parallel do
         do eID = 1, this % num_of_elements
            this % elements(eID) % stats = this % elements(eID) % stats + weight * contribution % elements(eID) % stats            
         end do
!$omp end parallel do
         
      end subroutine Mesh_AddWeightedContribution
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Mesh_ReadStatistics(self,solutionName)
         implicit none
         !-arguments---------------------------------------------------
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: solutionName
         !-local-variables---------------------------------------------
         integer       :: arrayDimensions(4)
         integer       :: fid, eID         
         !-------------------------------------------------------------
         
!~         self % solutionName = trim(solutionName)
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(solutionName)
         
         do eID = 1, self % num_of_elements
            associate ( e => self % elements(eID) )
            call getSolutionFileArrayDimensions(fid,arrayDimensions) ! needed to get to the position of stats
            read(fid) e % stats
            end associate
         end do
!
!        Close file
!        ----------
         close(fid)

      end subroutine Mesh_ReadStatistics
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Mesh_SaveStatistics(self, name)
         implicit none
         !-arguments---------------------------------------------------
         class(Mesh_t),       intent(in)        :: self
         character(len=*),    intent(in)        :: name
         !-local-variables---------------------------------------------
         integer  :: fid, eID
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
         !-------------------------------------------------------------
!
!        Create new file
!        ---------------
         call CreateNewSolutionFile(trim(name),STATS_FILE, self % nodeType, self % num_of_elements, self % iter, self % time, self % refs)
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % num_of_elements
            associate( e => self % elements(eID) )
            call writeArray(fid, e % stats)
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine Mesh_SaveStatistics
end module StatsAveragingModule