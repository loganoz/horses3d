module ResidualsMonitorClass
   use SMConstants
   use PhysicsStorage
   use HexMeshClass
   use MonitorDefinitions

   private
   public   Residuals_t
!
!  **************************
!  Residuals class definition
!  **************************
!
   type Residuals_t
      logical                         :: active
      real(kind=RP)                   :: values(NCONS,BUFFER_SIZE)
      character(len=STR_LEN_MONITORS) :: fileName
      contains
         procedure   :: Initialization => Residuals_Initialization
         procedure   :: Update         => Residuals_Update
         procedure   :: WriteLabel     => Residuals_WriteLabel
         procedure   :: WriteValues    => Residuals_WriteValue
         procedure   :: WriteToFile    => Residuals_WriteToFile
   end type Residuals_t
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           RESIDUALS ROUTINES
!           ------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Residuals_Initialization( self, solution_file ) 
!
!        *******************************************************************
!              This subroutine initializes the residuals structure
!        *******************************************************************
!
         implicit none
         class(Residuals_t) :: self
         character(len=*)   :: solution_file
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: fileName
         integer                          :: fID
         integer                          :: pos
!
!        Enable the monitor
!        ------------------
         self % active = .true.
!
!        Get monitor file name
!        ---------------------
         write( self % fileName , '(A,A)') trim(solution_file) , ".residuals"  
!
!        Create file to write the residuals
!        ----------------------------------
         open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
         write ( fID , ' ( A                                      ) ' ) "Residuals file"
#if defined(NAVIERSTOKES)
         write ( fID , ' ( A10,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24 ) ' ) "Iteration" , "Time" , "continuity" , &
                                                              "x-momentum" , "y-momentum" , "z-momentum", "energy"
#elif defined(CAHNHILLIARD)
         write ( fID , ' ( A10,2X,A24,2X,A24) ' ) "Iteration" , "Time" , "concentration"

#endif
!
!        Close file
!        ----------
         close ( fID ) 
              
      end subroutine Residuals_Initialization

      subroutine Residuals_Update ( self, mesh, maxResiduals, bufferPosition)
!
!        *********************************************************
!              This subroutine updates the residuals values from
!           those computed in the Monitor procedure
!        *********************************************************
!
         implicit none
         class(Residuals_t)         :: self
         class(HexMesh), intent(in) :: mesh
         real(kind=RP)              :: maxResiduals(NCONS)
         integer                    :: bufferPosition
!
!        Update buffer values
!        --------------------      
         self % values( 1:NCONS , bufferPosition ) = maxResiduals

      end subroutine Residuals_Update

      subroutine Residuals_WriteLabel ( self )
!
!        ************************************************************
!              This subroutine displays the residuals labels for the
!           time integrator Display procedure.
!        ************************************************************
!
         implicit none
         class(Residuals_t)             :: self
#if defined(NAVIERSTOKES)
         write(STD_OUT , '(3X,A10)' , advance = "no") "continuity"
         write(STD_OUT , '(3X,A10)' , advance = "no") "x-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "y-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "z-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "energy"

#elif defined(CAHNHILLIARD)
         write(STD_OUT , '(3X,A10)' , advance = "no") "concentration"

#endif

      end subroutine Residuals_WriteLabel
   
      subroutine Residuals_WriteValue ( self , bufferLine ) 
!
!        ***************************************************************
!              This subroutine displays the residuals values for the 
!           time integrator Display procedure
!        ***************************************************************
!
         implicit none
         class(Residuals_t) :: self
         integer            :: bufferLine
!        ---------------------------------------------------------
         integer            :: eq
      
         do eq = 1 , NCONS
            write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values(eq , bufferLine)
         end do

      end subroutine Residuals_WriteValue 

      subroutine Residuals_WriteToFile ( self , iter , t , no_of_lines)
!
!        *********************************************************************
!              This subroutine exports the results to the monitor file.
!           Just "no_of_lines" buffer lines are written.
!        *********************************************************************
!
         implicit none  
         class(Residuals_t)             :: self
         integer                    :: iter(:)
         real(kind=RP)              :: t(:)
         integer                    :: no_of_lines
!        -------------------------------------------
         integer                    :: i
         integer                    :: fID
!
!        Open file
!        ---------
         open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
!
!        Write values
!        ------------         
         do i = 1 , no_of_lines
            write( fID , '(I10,2X,ES24.16,5(2X,ES24.16))' ) iter(i) , t(i) , self % values(1:NCONS,i)
         end do
!
!        Close file
!        ----------        
         close ( fID )

         if ( no_of_lines .ne. 0 ) then
            self % values(1:NCONS,1) = self % values(1:NCONS,no_of_lines)
         end if
      
      end subroutine Residuals_WriteToFile

end module ResidualsMonitorClass
