#include "Includes.h"
module LoadBalancingMonitorClass
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use MPI_Process_Info

#ifdef _HAS_MPI_
   use mpi
#endif

   implicit none

   private

   public LoadBalancingMonitor_t

!
!  ***************************************
!  Load balancing monitor class definition
!  ***************************************
!
   type LoadBalancingMonitor_t
      logical                         :: active
      integer                         :: ID
      integer                         :: bufferLine
      integer                         :: num_of_vars
      real(kind=RP), allocatable      :: values(:,:)
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: variable
      contains
         procedure   :: Initialization => LoadBalancingMonitor_Initialization
         procedure   :: Update         => LoadBalancingMonitor_Update
         procedure   :: WriteLabel     => LoadBalancingMonitor_WriteLabel
         procedure   :: WriteValues    => LoadBalancingMonitor_WriteValue
         procedure   :: WriteToFile    => LoadBalancingMonitor_WriteToFile
         procedure   :: getLast        => LoadBalancingMonitor_GetLast
         procedure   :: destruct       => LoadBalancingMonitor_Destruct
         procedure   :: copy           => LoadBalancingMonitor_Assign
         generic     :: assignment(=)  => copy
   end type LoadBalancingMonitor_t
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////
!
!           LOAD BALANCING MONITOR PROCEDURES
!           ---------------------------------
!///////////////////////////////////////////////////////////////////////////
!
      subroutine LoadBalancingMonitor_Initialization( self , mesh , ID, solution_file , FirstCall)
!
!        *****************************************************************************
!              This subroutine initializes the load balancing monitor. The following
!           data is obtained from the case file:
!              -> Name: The monitor name (10 characters maximum)
!              -> Variable: The variable to be monitorized.
!        *****************************************************************************
!
         use ParamfileRegions
         use Utilities, only: toLower
         implicit none
         class(LoadBalancingMonitor_t) :: self
         class(HexMesh)         :: mesh
         integer                :: ID
         character(len=*)       :: solution_file
         logical, intent(in)    :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         integer                          :: fID
!
!        Get monitor ID
!        --------------
         self % ID = ID
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define load balancing monitor " , self % ID

         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "name"              , self % monitorName      , in_label , "# end" )
         call readValueInRegion ( trim ( paramFile )  , "variable"          , self % variable         , in_label , "# end" )
!
!        Enable the monitor
!        ------------------
         self % active = .true.
         self % bufferLine = 1
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
         self % num_of_vars = 1
         call toLower(self % variable)

         select case ( trim ( self % variable ) )
         case ("max dof per partition")
         case ("max dof")
         case ("min dof per partition")
         case ("min dof")
         case ("avg dof per partition")
         case ("avg dof")
         case ("absolute dof unbalancing")
         case ("abs unbalancing")
         case ("relative dof unbalancing")
         case ("rel unbalancing")

         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for load balancing monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" load balancing monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Max DOF per partition or Max DOF"
               print*, "   * Min DOF per partition or Min DOF"
               print*, "   * Avg DOF per partition or Avg DOF"
               print*, "   * Absolute DOF unbalancing or Abs unbalancing"
               print*, "   * Relative DOF unbalancing or Rel unbalancing"
               error stop "error stopped."
            end if
         end select

         allocate ( self % values(self % num_of_vars,BUFFER_SIZE) )

!
!        Prepare the file in which the monitor is exported
!        -------------------------------------------------
         write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % monitorName) , ".loadbalancing"
!
!        Create file
!        -----------
         if (FirstCall) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" )
!
!        Write the file headers
!        ----------------------
            write( fID , '(A20,A  )') "#Monitor name:      ", trim(self % monitorName)
            write( fID , '(A20,A  )') "#Selected variable: " , trim(self % variable)

            write( fID , * )

            write( fID , '(A10,2X,A24,2X,A24)'   ) "#Iteration" , "Time" , trim(self % variable)

            close ( fID )
         end if

      end subroutine LoadBalancingMonitor_Initialization

      subroutine LoadBalancingMonitor_Update ( self , mesh , bufferPosition )
!
!        *******************************************************************
!           This subroutine updates the monitor value computing it from
!           the mesh. It is stored in the "bufferPosition" position of the
!           buffer.
!        *******************************************************************
!
         implicit none
         class   (  LoadBalancingMonitor_t ) :: self
         class   (  HexMesh       )          :: mesh
         integer                             :: bufferPosition
!        --- local-variables ---------------------------------
         integer                             :: maxDOF, minDOF, avgDOF
         integer                             :: partitionDOF(MPI_Process % nProcs)

         self % bufferLine = bufferPosition
!
!        Compute the volume integral
!        ---------------------------
         select case ( trim(self % variable) )
#ifdef _HAS_MPI_
         case ("max dof per partition")
            self % values(1,bufferPosition) = maxval(CumputePartitionDOF(mesh))
         case ("max dof")
            self % values(1,bufferPosition) = maxval(CumputePartitionDOF(mesh))
         case ("min dof per partition")
            self % values(1,bufferPosition) = minval(CumputePartitionDOF(mesh))
         case ("min dof")
            self % values(1,bufferPosition) = minval(CumputePartitionDOF(mesh))
         case ("avg dof per partition")
            self % values(1,bufferPosition) = sum(CumputePartitionDOF(mesh)) / MPI_Process % nProcs
         case ("avg dof")
            self % values(1,bufferPosition) = sum(CumputePartitionDOF(mesh)) / MPI_Process % nProcs
         case ("absolute dof unbalancing")
            partitionDOF = CumputePartitionDOF(mesh)
            maxDOF = maxval(partitionDOF)
            minDOF = minval(partitionDOF)
            self % values(1,bufferPosition) = maxDOF - minDOF
         case ("abs unbalancing")
            partitionDOF = CumputePartitionDOF(mesh)
            maxDOF = maxval(partitionDOF)
            minDOF = minval(partitionDOF)
            self % values(1,bufferPosition) = maxDOF - minDOF
         case ("relative dof unbalancing")
            partitionDOF = CumputePartitionDOF(mesh)
            maxDOF = maxval(partitionDOF)
            minDOF = minval(partitionDOF)
            avgDOF = sum(partitionDOF) / MPI_Process % nProcs
            self % values(1,bufferPosition) = 100.0_RP * (maxDOF - minDOF) / avgDOF
         case ("rel unbalancing")
            partitionDOF = CumputePartitionDOF(mesh)
            maxDOF = maxval(partitionDOF)
            minDOF = minval(partitionDOF)
            avgDOF = sum(partitionDOF) / MPI_Process % nProcs
            self % values(1,bufferPosition) = 100.0_RP * (maxDOF - minDOF) / avgDOF
#endif
         end select


      end subroutine LoadBalancingMonitor_Update

      subroutine LoadBalancingMonitor_WriteLabel ( self )
!
!        *************************************************************
!              This subroutine writes the label for the load balancing
!           monitor, when invoked from the time integrator Display
!           procedure.
!        *************************************************************
!
         implicit none
         class(LoadBalancingMonitor_t)             :: self

         write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))

      end subroutine LoadBalancingMonitor_WriteLabel
!
!        *************************************************************
!        WriteValue: This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
      subroutine LoadBalancingMonitor_WriteValue ( self , bufferLine )
         implicit none
         !-arguments-----------------------------------------------
         class(LoadBalancingMonitor_t)     :: self
         integer                           :: bufferLine
         !-local-variables-----------------------------------------
         integer :: i
         !---------------------------------------------------------

         do i=1, self % num_of_vars
            write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values (i, bufferLine )
         end do

      end subroutine LoadBalancingMonitor_WriteValue

      subroutine LoadBalancingMonitor_WriteToFile ( self , iter , t , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none
         class(LoadBalancingMonitor_t) :: self
         integer                       :: iter(:)
         real(kind=RP)                 :: t(:)
         integer                       :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i
         integer                       :: fID
         character(len=LINE_LENGTH)    :: fmt
!        ---------------------------------------------------------

         if ( MPI_Process % isRoot ) then
            open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )

            write(fmt,'(A,I0,A)') '(I10,2X,ES24.16,', size(self % values,1), '(2X,ES24.16))'
            do i = 1 , no_of_lines
               write( fID , fmt ) iter(i) , t(i) , self % values(:,i)

            end do

            close ( fID )
         end if

         if ( no_of_lines .ne. 0 ) self % values(:,1) = self % values(:,no_of_lines)

      end subroutine LoadBalancingMonitor_WriteToFile
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      function LoadBalancingMonitor_GetLast(self) result(lastValues)
         implicit none
         class(LoadBalancingMonitor_t), intent(in) :: self
         real(kind=RP)                             :: lastValues ( size(self % values,1) )

         lastValues(:) = self % values (:,self % bufferLine)

      end function LoadBalancingMonitor_GetLast
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine LoadBalancingMonitor_Destruct (self)
         implicit none
         class(LoadBalancingMonitor_t), intent(inout) :: self

         deallocate (self % values)
      end subroutine LoadBalancingMonitor_Destruct

      elemental subroutine LoadBalancingMonitor_Assign (to, from)
         implicit none
         class(LoadBalancingMonitor_t), intent(inout)  :: to
         type(LoadBalancingMonitor_t) , intent(in)     :: from

         to % active       = from % active
         to % ID           = from % ID

         safedeallocate (to % values)
         allocate ( to % values( size(from % values,1) , size(from % values,2) ) )
         to % values       = from % values

         to % monitorName  = from % monitorName
         to % fileName     = from % fileName
         to % variable     = from % variable
      end subroutine LoadBalancingMonitor_Assign

!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!     Additional functions
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      function CumputePartitionDOF(mesh) result(partitionDOF)
         implicit none
         class(HexMesh)                :: mesh
         integer                       :: partitionDOF(MPI_Process % nProcs)
!        ---- local-variables ----------------------------------------------------
         integer                       :: eID, ierr
         integer                       :: Nx(mesh % no_of_elements), Ny(mesh % no_of_elements), Nz(mesh % no_of_elements)
         integer                       :: local_DOFs
!        ---------------------------------------------------------------------------
!
!        Get polynomial order in all partitions
!        ----------------------------------------
!$omp parallel do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            Nx(eID) = mesh % elements(eID) % Nxyz(1)
            Ny(eID) = mesh % elements(eID) % Nxyz(2)
            Nz(eID) = mesh % elements(eID) % Nxyz(3)
         enddo
!$omp end parallel do 

!
!     --------------------------------------------------------------------------
!     Perform a reduction to know how many DOFs are in each process
!     --------------------------------------------------------------------------
      local_DOFs = SUM((Nx(:)+1)*(Ny(:)+1)*(Nz(:)+1))

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
      call mpi_allgather(local_DOFs, 1, MPI_INT, partitionDOF, 1, MPI_INT, MPI_COMM_WORLD, ierr)
#endif
      else
         partitionDOF(1) = local_DOFs
      end if

      end function CumputePartitionDOF

end module LoadBalancingMonitorClass