#include "Includes.h"
module SurfaceSampling
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use MPI_Process_Info
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString, GetRealValue, getCharArrayFromString
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
 
   private
   public   SurfaceSampling_t


!
!  ********************************
!  Surface Sampling class definition
!  ********************************
!
   type SurfaceSampling_t
      logical                         :: active=.false.
      logical                         :: isDimensionless
      integer                         :: ID
	  integer						  :: nVariables
      integer                         :: marker
	  integer                         :: rank
	  integer 						  :: interval
	  integer                         :: bufferSize
	  integer  						  :: bufferLine
	  integer						  :: intervalCount
	  integer, allocatable      	  :: nData(:)
      real(kind=RP), allocatable      :: values(:,:,:)
      character(len=STR_LEN_MONITORS) :: SamplingName
      character(len=STR_LEN_MONITORS), allocatable :: fileName(:)
      character(len=STR_LEN_MONITORS), allocatable :: variable(:)
      contains
         procedure   :: Initialization => SurfaceSampling_Initialization
         procedure   :: Update         => SurfaceSampling_Update
         procedure   :: WriteToFile    => SurfaceSampling_WriteToFile
         procedure   :: destruct       => SurfaceSampling_Destruct
         procedure   :: copy           => SurfaceSampling_Assign
         generic     :: assignment(=)  => copy
   end type SurfaceSampling_t

   contains
!
!/////////////////////////////////////////////////////////////////////////
!
!           SURFACE Sampling PROCEDURES
!           --------------------------
!/////////////////////////////////////////////////////////////////////////
!
      subroutine SurfaceSampling_Initialization( self , mesh , ID, solution_file , FirstCall)
!
!        *****************************************************************************
!              This subroutine initializes the surface Sampling. The following
!           data is obtained from the case file:
!              -> Name: The Sampling name (10 characters maximum)
!              -> Marker: The surface marker in which the Sampling will be computed.
!              -> Variable: The variable to be Samplingized.
!        *****************************************************************************
!  
         use Headers
		 use ParamfileRegions
		 use MPI_Process_Info
         implicit none
         class(SurfaceSampling_t):: self
         class(HexMesh)          :: mesh
         integer                 :: ID
         character(len=*)        :: solution_file
         logical, intent(in)     :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         integer, allocatable             :: marker
         character(len=STR_LEN_MONITORS)  :: markerName
		 character(len=STR_LEN_MONITORS)  :: bufferInt
		 character(len=STR_LEN_MONITORS)  :: formatFile
		 character(len=STR_LEN_MONITORS)  :: variables
		 character(len=STR_LEN_MONITORS)  :: writeInterval
         integer                          :: pos
         integer                          :: fID
         integer                          :: zoneID
		 integer						  :: Nf(2)
		 integer						  :: nData
		 integer						  :: allnData(MPI_Process % nProcs), ierr
		 integer						  :: zonefID
		 integer						  :: i,j,k
		 real(kind=RP), allocatable		  :: X(:,:)
		 real(kind=RP), allocatable		  :: X1(:,:)
		 real(kind=RP), allocatable		  :: X2(:,:)
		 real(kind=RP), allocatable		  :: X3(:,:)
		
!
!        Get Sampling ID, assign zero to bufferLine and intervalCount
!        --------------
         self % ID = ID
		 self % bufferLine = 0
		 self % intervalCount = 0
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define surface sampling " , self % ID
         
         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "name"              , self % SamplingName      , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "surface"           , markerName               , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "variables"         , variables                , in_label , "# end" ) 
		 call readValueInRegion ( trim ( paramFile )  , "sampling interval" , bufferInt                , in_label , "# end" ) 
		 call readValueInRegion(trim(paramFile), "write interval", writeInterval      , in_label, "#end" )
!
!        Enable the Sampling
!        ------------------
         self % active = .true.
!
!        Get the variables and its size
!        ------------------------------		 
		 call getCharArrayFromString(variables, STR_LEN_MONITORS, self % variable)
		 self % nVariables = size(self % variable) 
!
!        Get the surface marker
!        ----------------------
         self % marker = -1
         do zoneID = 1, size(mesh % zones)
            if ( trim(mesh % zones(zoneID) % name) .eq. trim(markerName) ) then
               self % marker = zoneID
               exit
            end if
         end do

         if ( self % marker .eq. -1 ) then
            self % active = .false.
			if (MPI_Process % isRoot ) then 
            write(*,'(A,I0)') "Warning: Marker not specified for surface sampling ", self % ID
            write(*,'(A,I0,A)') "     Surface Sampling ", self % ID, " disabled."
			end if 
         end if
		 
		 if (mesh % zones(self % marker) % no_of_faces .eq. 0) then
			self % active = .false. 
		 end if 
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
!		
		DO i=1, self % nVariables
!        ****************************************
         select case ( trim ( self % variable(i) ) )
!        ****************************************
!
			case ("shearstress-tangent")
               self % isDimensionless = .false.
			   
			case ("shearstress-x")
               self % isDimensionless = .false.
			
			case ("shearstress-y")
               self % isDimensionless = .false.
			   
			case ("shearstress-z")
               self % isDimensionless = .false.
			   
			case ("pressure")
               self % isDimensionless = .false.
			
			case ("q1")
               self % isDimensionless = .true.
			   
			case ("q2")
               self % isDimensionless = .true.
			   
			case ("q3")
               self % isDimensionless = .true.
			
			case ("q4")
               self % isDimensionless = .true.
			   
			case ("q5")
               self % isDimensionless = .true.   
			   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            case default
				if (MPI_Process % isRoot ) then 
				   if ( len_trim (self % variable(i)) .eq. 0 ) then
					  print*, "Variable was not specified for surface Sampling " , self % ID , "."
				   else
					  print*, 'Variable "',trim(self % variable(i)),'" surface Sampling ', self % ID, ' not implemented yet.'
					  print*, "Options available are:"
					  print*, "   * shearstress-tangent"
					  print*, "   * shearstress-x"
					  print*, "   * shearstress-y"
					  print*, "   * shearstress-z"
					  print*, "   * pressure"
					  stop "Stopped."

				   end if
				end if 
!
!        **********
         end select
!        **********
		END DO
!
!        Get the number of Order and the sampling interval
!        -------------------------------------------------
		 if (mesh % zones(self % marker) % no_of_faces .gt.0) then
			Nf(1) = mesh % faces (mesh % zones(self % marker) % faces(1))%Nf(1)+1
			Nf(2) = mesh % faces (mesh % zones(self % marker) % faces(1))%Nf(2)+1
			nData = mesh % zones(self % marker) % no_of_faces * Nf(1) * Nf(2)
		 else
			self % active = .false. 
			Nf(1)=0
			Nf(2)=0
			nData=0
		 end if 
		 self % interval = GetRealValue(bufferInt)
		 
		 ALLOCATE (self % nData(MPI_Process % nProcs))
		 
         if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
			call mpi_allgather(nData, 1, MPI_INT, allnData, 1, MPI_INT, MPI_COMM_WORLD, ierr)
			nData=sum(allnData)
			self % nData = allnData
			self % rank = maxloc(allnData,1)-1
			call mpi_allgather(Nf(1), 1, MPI_INT, allnData, 1, MPI_INT, MPI_COMM_WORLD, ierr)
			Nf(1)=maxval(allnData)
			call mpi_allgather(Nf(2), 1, MPI_INT, allnData, 1, MPI_INT, MPI_COMM_WORLD, ierr)
			Nf(2)=maxval(allnData)
#endif  
		end if 
			
!
!       Get the max. number of timestep in the buffer file before being written
!       -----------------------------------------------------------------------			
		IF (LEN(TRIM(writeInterval)) .EQ. 0) THEN
			self % bufferSize = 1;
		ELSE
			self % bufferSize = GetRealValue(writeInterval)
!
!           Failsafe to prevent too many data being written at one time
!           -----------------------------------------------------------
			IF (nData .GT. 200000) THEN   
				self % bufferSize = 1;
			END IF
		END IF

		 
		 ALLOCATE( X(nDATA,3), self % fileName(self % nVariables), X1(nData, MPI_Process % nProcs), &
					X2(nData, MPI_Process % nProcs), X3(nData, MPI_Process % nProcs) )
		
		if ( MPI_Process % isRoot ) then
			ALLOCATE(self % values(nData+1,self % bufferSize , self % nVariables))
		end if 
						
		 k=0
		 do zonefID = 1, mesh % zones(self % marker) % no_of_faces
!
!           Face global ID
!           --------------
            fID = mesh % zones(self % marker) % faces(zonefID)
			do j = 0, Nf(2)-1;    do i = 0, Nf(1)-1
				k=k+1
				X(k,1)=mesh % faces(fID) % geom % x(1,i,j)
				X(k,2)=mesh % faces(fID) % geom % x(2,i,j)
				X(k,3)=mesh % faces(fID) % geom % x(3,i,j)
			end do;			end do
			
		end do
		
		 if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
			call mpi_allgather(X(:,1), nData, MPI_DOUBLE, X1, nData, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
			call mpi_allgather(X(:,2), nData, MPI_DOUBLE, X2, nData, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
			call mpi_allgather(X(:,3), nData, MPI_DOUBLE, X3, nData, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
			if ( MPI_Process % isRoot ) then
				X(:,1)=X1(:,1)
				X(:,2)=X2(:,1)
				X(:,3)=X3(:,1)
				do i=1, MPI_Process % nProcs-1
					X(sum(self % nData(1:i))+1:sum(self % nData(1:i+1)),1)=X1(1:self % nData(i+1),i+1)
					X(sum(self % nData(1:i))+1:sum(self % nData(1:i+1)),2)=X2(1:self % nData(i+1),i+1)
					X(sum(self % nData(1:i))+1:sum(self % nData(1:i+1)),3)=X3(1:self % nData(i+1),i+1)
				end do
			end if
#endif  
		end if 
		
		
		DO i=1, self % nVariables
!
!        Prepare the file in which the Sampling is exported
!        -------------------------------------------------
         write( self % fileName(i) , '(A,A,A,A,A,A,I0,A)') trim(solution_file) , "_" , trim(markerName) , "_" , trim(self % variable(i)) &
									, "_surface_" , self % ID,".sampling"  
									
#if defined(NAVIERSTOKES) && (!(INCNS))
!
!        Create file
!        -----------
         if (FirstCall) then
		    if (MPI_Process % isRoot ) then 
				open ( newunit = fID , file = trim(self % fileName(i)) , action = "write" , access = "stream" , status = "replace", position='append' )
				
!
!        Write the file headers
!        ----------------------
				write( fID) self % ID
				write( fID) Nf(1)
				write( fID) Nf(2)
				write( fID) nData
				write( fID ) self % interval
				write( fID ) refValues % rho
				write( fID ) refValues % V
				write( fID ) refValues % p
				write( fID ) refValues % T
				write( fID ) refValues % mu
				write( fID ) refValues % AoATheta
				write( fID ) refValues % AoAPhi
				write( fID ) X(:,1)
				write( fID ) X(:,2)
				write( fID ) X(:,3)
				close ( fID )
			end if 
         end if
#endif
		END DO
		
!
!        Write Information into log
!        --------------------------
         write( formatFile , '(A,A,A,A,I0,A)') trim(solution_file) , "_" , trim(markerName) , "_'variable'_surface_", self % ID, ".sampling"  
		 
		 if ( .not. MPI_Process % isRoot ) return
!
!        Write Information 
!        -----------------------------------------------
		 
		 write(STD_OUT,'(/)')
		 call SubSection_Header("Surface Samplings")
			write(STD_OUT,'(30X,A,A27,I4)') "->" , "Surface Sampling ID: " , self % ID
			write(STD_OUT,'(30X,A,A27,A27)') "->" , "Surface Name: " , markerName
			write(STD_OUT,'(30X,A,A27,A128)') "->" , "Variables: " , variables
			write(STD_OUT,'(30X,A,A27,I4)') "->" , "Samplings Interval: ", self % interval
			write(STD_OUT,'(30X,A,A27,A128)') "->" , "Filename: ", formatFile
			write(STD_OUT,'(30X,A,A27,A25)') "->" ,"Note: ","Extracted with dimension"
		 
	     deallocate(X,X1,X2,X3)
		 
      end subroutine SurfaceSampling_Initialization

      subroutine SurfaceSampling_Update ( self, mesh, bufferPosition, t )
!
!        *******************************************************************
!           This subroutine updates the Sampling value computing it from
!           the mesh. It is stored in the "bufferPosition" position of the 
!           buffer.
!        *******************************************************************
!
         use SamplingOperator
         implicit none
         class   (  SurfaceSampling_t )  :: self
         class   (  HexMesh       )      :: mesh
         integer                         :: bufferPosition
		 real(kind=RP)                   :: t
		 integer                         :: i, j, k, recv_req(MPI_Process % nProcs-1), send_req, ierr
         real(kind=RP)                   :: F(NDIM)
		 real(kind=RP), allocatable      :: data_out(:)

		 
		 if (self % intervalCount .EQ. 0 ) then
		    self % bufferLine = self % bufferLine + 1
			
			if ((self % nData(MPI_Process % rank +1) .gt.0 ).or. MPI_Process % isRoot) then
			DO i=1, self % nVariables
			select case ( trim ( self % variable(i) ) )

				case ("shearstress-tangent")
					Call VectorSurfaceSampling(mesh, self % marker, SHEAR_STRESS_TANGENT, self % SamplingName, data_out)

				case ("shearstress-x")
					Call VectorSurfaceSampling(mesh, self % marker, SHEAR_STRESS_X, self % SamplingName, data_out)
					
				case ("shearstress-y")
					Call VectorSurfaceSampling(mesh, self % marker, SHEAR_STRESS_Y, self % SamplingName, data_out)
					
				case ("shearstress-z")
					Call VectorSurfaceSampling(mesh, self % marker, SHEAR_STRESS_Z, self % SamplingName, data_out)
					
				case ("pressure")
					Call VectorSurfaceSampling(mesh, self % marker, PRESSURE_SURF, self % SamplingName, data_out)
					
				case ("q1")
					Call VectorSurfaceSampling(mesh, self % marker, Q1, self % SamplingName, data_out)
					
				case ("q2")
					Call VectorSurfaceSampling(mesh, self % marker, Q2, self % SamplingName, data_out)
					
				case ("q3")
					Call VectorSurfaceSampling(mesh, self % marker, Q3, self % SamplingName, data_out)
					
				case ("q4")
					Call VectorSurfaceSampling(mesh, self % marker, Q4, self % SamplingName, data_out)
					
				case ("q5")
					Call VectorSurfaceSampling(mesh, self % marker, Q5, self % SamplingName, data_out)

			end select
			
				if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
					if ( MPI_Process % isRoot ) then
					
						self % values(:,bufferPosition,i)=0_RP
						self % values(1,bufferPosition,i)=t
						if (self % nData(1) .gt.0 ) then
							self % values(2:self % nData(1) +1,bufferPosition,i)=data_out
						end if
						k=0
						DO j=1, MPI_Process % nProcs -1
							if (self % nData(j+1) .gt.0 ) then
								k=k+1
								call mpi_irecv(self % values(sum(self % nData(1:j))+1:sum(self % nData(1:j))+1+self % nData(j+1),bufferPosition,i), &
									self % nData(j+1), MPI_DOUBLE, j, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(k), ierr)
							    call mpi_wait(recv_req(k), MPI_STATUS_IGNORE, ierr) 
							end if 
						END DO 
				    else 
					
						if (self % nData(MPI_Process % rank +1) .gt.0 ) then
						   call mpi_isend(data_out, self % nData(MPI_Process % rank +1), MPI_DOUBLE, 0, &
								DEFAULT_TAG, MPI_COMM_WORLD, send_req, ierr)
						   call mpi_wait(send_req, MPI_STATUS_IGNORE, ierr) 
						   send_req=0
						end if 
					end if 
#endif  
				else 
					self % values(1,bufferPosition,i)=t
					self % values(2:size(data_out,1)+1,bufferPosition,i)=data_out
				end if 
				
			END DO
			
			end if 
					
					
		 end if
		 
		 self % intervalCount = self % intervalCount + 1
		 
		 if (self % intervalCount .EQ. self % interval) then
			self % intervalCount = 0 
		 end if 

      end subroutine SurfaceSampling_Update

      subroutine SurfaceSampling_WriteToFile ( self , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(SurfaceSampling_t) :: self
         integer                  :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i,k
         integer                    :: fID
		 
		 if ( MPI_Process % isRoot ) then
		 
		 DO k=1, self % nVariables

            open( newunit = fID , file = trim ( self % fileName(k) ) , action = "write" , access = "stream" , status = "old", position='append' )
         
            do i = 1 , no_of_lines
               write( fID ) self % values(:,i,k)
   
            end do
        
            close ( fID )
         
			if ( no_of_lines .ne. 0 ) self % values(:,1,k) = self % values(:,no_of_lines,k)
		 
		 END DO
		 
		 end if
		 
		 self % bufferLine = 0
      
      end subroutine SurfaceSampling_WriteToFile
!
      elemental subroutine SurfaceSampling_Destruct (self)
         implicit none
         class(SurfaceSampling_t), intent(inout) :: self
         
         safedeallocate (self % values)
		 safedeallocate (self % nData)
		 safedeallocate (self % fileName)
		 safedeallocate (self % variable)
		 
      end subroutine SurfaceSampling_Destruct
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine SurfaceSampling_Assign(to, from)
         implicit none
         class(SurfaceSampling_t), intent(inout) :: to
         type(SurfaceSampling_t),  intent(in)    :: from
         
         if ( from % active) then
            to % active          = from % active
            to % isDimensionless = from % isDimensionless
            to % ID              = from % ID
			to % nVariables      = from % nVariables
            to % marker          = from % marker
			to % interval        = from % interval
			to % bufferSize      = from % bufferSize
			to % bufferLine      = from % bufferLine
			to % intervalCount   = from % intervalCount
			
            safedeallocate(to % values)
            allocate ( to % values( size(from % values,1),size(from % values,2),size(from % values,3) ) )
            to % values          = from % values
			
            to % SamplingName    = from % SamplingName
			
			safedeallocate(to % fileName)
            allocate ( to % fileName( size(from % fileName) ) )
            to % fileName         = from % fileName
			
			safedeallocate(to % variable)
            allocate ( to % variable( size(from % variable) ) )
            to % variable         = from % variable

         else
            to % active = .FALSE.
         end if
         
      end subroutine SurfaceSampling_Assign
end module SurfaceSampling

