#include "Includes.h"
module Samplings
   use SMConstants
   use NodalStorageClass
   use HexMeshClass
   use MonitorDefinitions
   use FileReadingUtilities      , only: getFileName
   use SurfaceSampling
   use PlaneSampling
   use SpatialMeanNode
   implicit none
!

   private
   public      Sampling_t
!
!  *****************************
!  Main sampling class definition
!  *****************************
!  
   type Sampling_t
      character(len=LINE_LENGTH)           :: solution_file
      integer                              :: no_of_surfaceSamplings
	  integer  							   :: no_of_planeSamplings
	  integer  							   :: no_of_spatialMeanNodes
      integer                              :: dt_restriction
      logical                              :: write_dt_restriction
      class(SurfaceSampling_t), allocatable :: surfaceSamplings(:)
	  class(PlaneSampling_t),   allocatable :: planeSamplings(:)
	  class(SpatialMeanNode_t),   allocatable :: spatialMeanNodes(:)
      contains
         procedure   :: Construct       => Samplings_Construct
		 procedure   :: UpdateInterp    => Samplings_UpdateLagrangeInterp
         procedure   :: UpdateValues    => Sampling_UpdateValues
         procedure   :: WriteToFile     => Sampling_WriteToFile
         procedure   :: destruct        => Sampling_Destruct
         procedure   :: copy            => Sampling_Assign
         generic     :: assignment(=)   => copy
   end type Sampling_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Samplings_Construct( Samplings, mesh, controlVariables )
		 use Headers
         use FTValueDictionaryClass
         use mainKeywordsModule
         implicit none
         class(Sampling_t)                    :: Samplings
         class(HexMesh), intent(in)           :: mesh
         class(FTValueDictionary), intent(in) :: controlVariables
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: fID , io
         integer                         :: i
         character(len=STR_LEN_MONITORS) :: line
         character(len=STR_LEN_MONITORS) :: solution_file                                            
         logical, save                   :: FirstCall = .TRUE.
!
!        Setup the buffer
!        ----------------
         if (controlVariables % containsKey("monitors flush interval") ) then
            BUFFER_SIZE = controlVariables % integerValueForKey("monitors flush interval")
         end if
!
!        Get the solution file name
!        --------------------------
         solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_MONITORS )
!
!        Remove the *.hsol termination
!        -----------------------------
         solution_file = trim(getFileName(solution_file))
         Samplings % solution_file = trim(solution_file)
!
!        Search in case file for samplings
!        ---------------------------------------------------------------------
         if (mesh % child) then ! Return doing nothing if this is a child mesh
            Samplings % no_of_surfaceSamplings = 0
			Samplings % no_of_planeSamplings   = 0
			Samplings % no_of_spatialMeanNodes   = 0
         else
            call getNoOfSamplings( Samplings % no_of_surfaceSamplings, Samplings % no_of_planeSamplings, Samplings % no_of_spatialMeanNodes)
         end if
!
!        Initialize
!        ----------


         allocate ( Samplings % surfaceSamplings ( Samplings % no_of_surfaceSamplings )  )
		 allocate ( Samplings % planeSamplings   ( Samplings % no_of_planeSamplings   )  )
		 allocate ( Samplings % spatialMeanNodes ( Samplings % no_of_spatialMeanNodes   )  )
		 
		 if (Samplings % no_of_surfaceSamplings .GT. 0) then
			call Section_Header("Initialize Surface Samplings")
		 end if
         do i = 1 , Samplings % no_of_surfaceSamplings
            call Samplings % surfaceSamplings(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do

		 if (Samplings % no_of_planeSamplings .GT. 0) then
			call Section_Header("Initialize Plane Samplings")
		 end if
		 do i = 1 , Samplings % no_of_planeSamplings
            call Samplings % planeSamplings(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
		 
		 if (Samplings % no_of_spatialMeanNodes .GT. 0) then
			call Section_Header("Initialize Spatial Mean Node")
		 end if
		 do i = 1 , Samplings % no_of_spatialMeanNodes
            call Samplings % spatialMeanNodes(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
		
         FirstCall = .FALSE.
		 
      end subroutine Samplings_Construct
	  
	  subroutine Samplings_UpdateLagrangeInterp (self, mesh)
!
!        *******************************************************************
!        This subroutine updates the Lagrange interpolants after pAdaptation
!        *******************************************************************
!        
         implicit none
         class(Sampling_t)   :: self
         class(HexMesh)      :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i 
		 
!
!        Update interpolants plane Samplings
!        -----------------------------------
         do i = 1 , self % no_of_planeSamplings
            call self % planeSamplings(i) % UpdateInterp( mesh )
         end do

      end subroutine Samplings_UpdateLagrangeInterp

	  subroutine Sampling_UpdateValues (self, mesh, t)
!
!        ***************************************************************
!              This subroutine updates the values for the Samplings.
!        ***************************************************************
!        
         use PhysicsStorage
         use StopwatchClass
         implicit none
         class(Sampling_t)   :: self
         class(HexMesh)      :: mesh
		 real(kind=RP)       :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i 
		 
!
!        Update surface Samplings
!        ------------------------
         do i = 1 , self % no_of_surfaceSamplings
            call self % surfaceSamplings(i) % Update( mesh , self % surfaceSamplings(i) % bufferLine, t )
         end do
!
!        Update plane Samplings
!        ----------------------
         do i = 1 , self % no_of_planeSamplings
            call self % planeSamplings(i) % Update( mesh , self % planeSamplings(i) % bufferLine, t )
         end do
!
!        Update spatial mean node
!        ------------------------
         do i = 1 , self % no_of_spatialMeanNodes
            call self % spatialMeanNodes(i) % Update( mesh , self % spatialMeanNodes(i) % bufferLine, t )
         end do
		 
      end subroutine Sampling_UpdateValues

      subroutine Sampling_WriteToFile ( self , mesh, force) 
!
!        ******************************************************************
!              This routine has a double behaviour:
!           force = .true.  -> Writes to file and resets buffers
!           force = .false. -> Just writes to file if the buffer is full
!        ******************************************************************
!
         use MPI_Process_Info
         implicit none
         class(Sampling_t)       :: self
         class(HexMesh)          :: mesh
         logical, optional       :: force
!        ------------------------------------------------
         integer                 :: i 
         logical                 :: forceVal
         if ( present ( force ) ) then
            forceVal = force

         else
            forceVal = .false.

         end if

         if ( forceVal ) then 
!
!           In this case the Samplings are exported to their files and the buffer is reseted
!           -------------------------------------------------------------------------------
  
            do i = 1 , self % no_of_surfaceSamplings
               call self % surfaceSamplings(i) % WriteToFile ( self % surfaceSamplings(i) % bufferLine )
            end do
			
			do i = 1 , self % no_of_planeSamplings
               call self % planeSamplings(i) % WriteToFile ( self % planeSamplings(i) % bufferLine )
            end do
			
			do i = 1 , self % no_of_spatialMeanNodes
               call self % spatialMeanNodes(i) % WriteToFile ( self % spatialMeanNodes(i) % bufferLine )
            end do

         else
!
!           The Samplings are exported just if the buffer is full
!           ----------------------------------------------------

            do i = 1 , self % no_of_surfaceSamplings
				if ( self % surfaceSamplings(i) % bufferLine .eq. self % surfaceSamplings(i) % bufferSize) then
					call self % surfaceSamplings(i) % WriteToFile ( self % surfaceSamplings(i) % bufferLine )
				end if
            end do
            do i = 1 , self % no_of_planeSamplings
				if ( self % planeSamplings(i) % bufferLine .eq. self % planeSamplings(i) % bufferSize) then
					call self % planeSamplings(i) % WriteToFile ( self % planeSamplings(i) % bufferLine )
				end if
            end do
			
			do i = 1 , self % no_of_spatialMeanNodes
				if ( self % spatialMeanNodes(i) % bufferLine .eq. self % spatialMeanNodes(i) % bufferSize) then
					call self % spatialMeanNodes(i) % WriteToFile ( self % spatialMeanNodes(i) % bufferLine )
				end if
            end do

         end if
      end subroutine Sampling_WriteToFile
      
      subroutine Sampling_Destruct (self)
         implicit none
         class(Sampling_t)        :: self
         
		 if ( self % no_of_surfaceSamplings .gt. 0)   call self % surfaceSamplings % destruct
		 if ( self % no_of_planeSamplings   .gt. 0)   call self % planeSamplings % destruct
		 if ( self % no_of_spatialMeanNodes .gt. 0)   call self % spatialMeanNodes % destruct

         safedeallocate (self % surfaceSamplings)
		 safedeallocate (self % planeSamplings)
		 safedeallocate (self % spatialMeanNodes)
		 
      end subroutine
      
      elemental subroutine Sampling_Assign ( to, from )
         implicit none
         !-arguments--------------------------------------
         class(Sampling_t), intent(inout)  :: to
         type(Sampling_t) , intent(in)     :: from
         !-local-variables--------------------------------
         !------------------------------------------------
         
         to % solution_file          = from % solution_file
         to % no_of_surfaceSamplings = from % no_of_surfaceSamplings
		 to % no_of_planeSamplings   = from % no_of_planeSamplings
		 to % no_of_spatialMeanNodes = from % no_of_spatialMeanNodes
         
         to % dt_restriction        = from % dt_restriction
         to % write_dt_restriction  = from % write_dt_restriction
         

         safedeallocate ( to % surfaceSamplings )
		 safedeallocate ( to % planeSamplings )
		 safedeallocate ( to % spatialMeanNodes )
         allocate ( to % surfaceSamplings ( size(from % surfaceSamplings) ) )
		 allocate ( to % planeSamplings   ( size(from % planeSamplings  ) ) )
		 allocate ( to % spatialMeanNodes ( size(from % spatialMeanNodes  ) ) )
         to % surfaceSamplings = from % surfaceSamplings
		 to % planeSamplings   = from % planeSamplings
		 to % spatialMeanNodes = from % spatialMeanNodes
         
      end subroutine Sampling_Assign
      
!
!//////////////////////////////////////////////////////////////////////////////
!
!        Auxiliars
!
!//////////////////////////////////////////////////////////////////////////////
!
   subroutine getNoOfSamplings(no_of_surfaceSamplings, no_of_planeSamplings, no_of_spatialMeanNodes)
      use ParamfileRegions
      implicit none
      integer, intent(out)    :: no_of_surfaceSamplings
	  integer, intent(out)    :: no_of_planeSamplings
	  integer, intent(out)    :: no_of_spatialMeanNodes
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
!
!     Initialize
!     ----------
      no_of_surfaceSamplings = 0
	  no_of_planeSamplings   = 0
	  no_of_spatialMeanNodes = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find Samplings
!     ------------------------------------
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

            if ( index ( line , '#definesurfacesampling' ) .gt. 0 ) then
               no_of_surfaceSamplings = no_of_surfaceSamplings + 1 
			
			else if ( index ( line , '#defineplanesampling' ) .gt. 0 ) then
               no_of_planeSamplings = no_of_planeSamplings + 1 
			   
			else if ( index ( line , '#definespatialmeannode' ) .gt. 0 ) then
               no_of_spatialMeanNodes = no_of_spatialMeanNodes + 1 
			   
            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

end subroutine getNoOfSamplings

end module Samplings
!
!///////////////////////////////////////////////////////////////////////////////////
!