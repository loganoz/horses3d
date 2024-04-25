module StatsAveragingModule
   use SMConstants
   use SolutionFile
   use StatisticsMonitor, only: NO_OF_VARIABLES_Sij
   use PhysicsStorage_NS, only: NCONS, NGRAD
   implicit none
   
   private
   public Initialize_StatsAveragingModule, Finalize_StatsAveragingModule, PerformAveraging, hasGradients
   
   type Element_t
!                                /* Solution quantities */
      integer                    :: Nsol(NDIM)
      real(kind=RP), pointer     :: stats(:,:,:,:)
      real(kind=RP), pointer     :: Q(:,:,:,:)
      real(kind=RP), pointer     :: Q_x(:,:,:,:)
      real(kind=RP), pointer     :: Q_y(:,:,:,:)
      real(kind=RP), pointer     :: Q_z(:,:,:,:)
      real(kind=RP)              :: offsetIO
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
         procedure   :: PrepareForIO            => Mesh_PrepareForIO
   end type Mesh_t
   
   type(Mesh_t)   :: AveragedSol, CurrentSol

   logical        :: hasGradients
   
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
         real(kind=RP), allocatable :: tempQStats(:,:,:,:)
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
            allocate (tempStats (1:NO_OF_VARIABLES_Sij,0:Nsol(1,eID),0:Nsol(2,eID),0:Nsol(3,eID)), tempQStats (1:NCONS,0:Nsol(1,eID),0:Nsol(2,eID),0:Nsol(3,eID)) )
            read(fid) tempStats
            read(fid) tempQStats
            if (hasGradients) then
                read(fid) tempQStats
                read(fid) tempQStats
                read(fid) tempQStats
            end if 
            deallocate (tempQStats,tempStats)
         end do
         close(fid)
         
         call AveragedSol % construct (num_of_elements, Nsol)
         call CurrentSol  % construct (num_of_elements, Nsol)
         call AveragedSol % PrepareForIO()
         
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
         integer        , intent(in)      :: Nsol(NDIM,num_of_elements)
         !-local-variables---------------------------------------------
         integer :: eID
         !-------------------------------------------------------------
         
         this % num_of_elements = num_of_elements
         allocate ( this % elements(num_of_elements) )
         
         do eID = 1, this % num_of_elements
            associate (e => this % elements(eID) )
            allocate( e % Q(1:NCONS,0:Nsol(1,eID), 0:Nsol(2,eID), 0:Nsol(3,eID)) , e % stats(1:NO_OF_VARIABLES_Sij,0:Nsol(1,eID), 0:Nsol(2,eID), 0:Nsol(3,eID)) )
            e % stats = 0._RP
            e % Q = 0._RP
            if (hasGradients) then
                allocate( e % Q_x(1:NCONS,0:Nsol(1,eID), 0:Nsol(2,eID), 0:Nsol(3,eID)), &
                          e % Q_y(1:NCONS,0:Nsol(1,eID), 0:Nsol(2,eID), 0:Nsol(3,eID)), &
                          e % Q_z(1:NCONS,0:Nsol(1,eID), 0:Nsol(2,eID), 0:Nsol(3,eID)) )
                e % Q_x = 0._RP
                e % Q_y = 0._RP
                e % Q_z = 0._RP
            endif
            e % Nsol(:) = Nsol(:,eID)
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
            this % elements(eID) % Q = this % elements(eID) % Q + weight * contribution % elements(eID) % Q            
            if (hasGradients) then
                this % elements(eID) % Q_x = this % elements(eID) % Q_x + weight * contribution % elements(eID) % Q_x            
                this % elements(eID) % Q_y = this % elements(eID) % Q_y + weight * contribution % elements(eID) % Q_y            
                this % elements(eID) % Q_z = this % elements(eID) % Q_z + weight * contribution % elements(eID) % Q_z            
            end if 
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
            read(fid) e % Q
            if (hasGradients) then
                read(fid) e % Q_x
                read(fid) e % Q_y
                read(fid) e % Q_z
            end if 
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
         integer  :: fid, eID, no_of_stats_variables
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
         integer(kind=AddrInt)            :: pos
         !-------------------------------------------------------------

         if (hasGradients) then
             no_of_stats_variables = NO_OF_VARIABLES_Sij + NCONS + NCONS*NDIM
         else
             no_of_stats_variables = NO_OF_VARIABLES_Sij
         end if 
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
                pos = POS_INIT_DATA + (eID-1)*5_AddrInt*SIZEOF_INT + no_of_stats_variables*e % offsetIO*SIZEOF_RP
                call writeArray(fid, e % stats, position=pos)
                write(fid) e % Q
                if (hasGradients) then
                    write(fid) e % Q_x
                    write(fid) e % Q_y
                    write(fid) e % Q_z
                end if 
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine Mesh_SaveStatistics
!
      Subroutine Mesh_PrepareForIO(self)
          Implicit None
         !-arguments---------------------------------------------------
         class(Mesh_t),       intent(inout)        :: self
         !-local-variables---------------------------------------------
         integer  :: eID
         integer, allocatable :: elementSizes(:)
!
!        Get each element size
!        ---------------------
         allocate(elementSizes(self % num_of_elements))
         do eID = 1, self % num_of_elements
            associate(e => self % elements(eID))
                elementSizes(eID) = product( e % Nsol + 1)
            end associate
         end do
!
         self % elements(1) % offsetIO = 0
         do eID = 2, self % num_of_elements
                self % elements(eID) % offsetIO = self % elements(eID-1) % offsetIO + elementSizes(eID-1)
         end do
!
         deallocate(elementSizes)
!
      End Subroutine Mesh_PrepareForIO
!
end module StatsAveragingModule
