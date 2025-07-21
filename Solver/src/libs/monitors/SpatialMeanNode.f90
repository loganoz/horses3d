#include "Includes.h"
module SpatialMeanNode
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use VariableConversion
   use MPI_Process_Info
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString, GetRealValue, getCharArrayFromString, GetIntValue, GetLogicalValue
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
   private
   public   SpatialMeanNode_t
!
!  *********************************
!  SpatialMeanNodes class definition
!  *********************************
!

   type SpatialMeanNode_t
      integer                         :: ID
      integer                         :: nVariables
      integer                         :: interval
      integer                         :: bufferSize
      integer                         :: bufferLine
      integer                         :: intervalCount
      integer                         :: nActive
      integer                         :: dirAxis
      integer                         :: nUniqueAll
      integer                         :: iVarU, iVarV, iVarW
      integer      , allocatable      :: activeLoc(:,:)
      integer      , allocatable      :: nMultiply(:)
      integer      , allocatable      :: nMultiplyAll(:)
      logical                         :: meanData = .false.
      real(kind=RP)                   :: pmin(3), pmax(3)
      real(kind=RP)                   :: error=0.000001 ! tolerance of coordinate
      real(kind=RP), allocatable      :: geom(:)        ! size nUnique
	  real(kind=RP), allocatable      :: meanU(:), meanV(:), meanW(:)
      real(kind=RP), allocatable      :: values(:,:,:)  ! (nUnique, bufferSize, nVariables)
      character(len=STR_LEN_MONITORS),  allocatable :: fileName (:)
      character(len=STR_LEN_MONITORS)               :: spatialMeanName
      character(len=STR_LEN_MONITORS),  allocatable :: variable (:)
      contains
         procedure   :: Initialization          => SpatialMeanNode_Initialization
         procedure   :: Update                  => SpatialMeanNode_Update
         procedure   :: WriteToFile             => SpatialMeanNode_WriteToFile
         procedure   :: LookForUniqueCoordinate => SpatialMeanNode_LookForUniqueCoordinate
         procedure   :: destruct              => SpatialMeanNode_Destruct
         procedure   :: copy                  => SpatialMeanNode_Assign
         generic     :: assignment(=)  => copy
   end type SpatialMeanNode_t

   contains

      subroutine SpatialMeanNode_Initialization(self, mesh, ID, solution_file, FirstCall)
         use Headers
         use ParamfileRegions
         use MPI_Process_Info
         use Utilities, only: toLower
         implicit none
         class(SpatialMeanNode_t)  :: self
         class(HexMesh)            :: mesh
         integer                   :: ID
         character(len=*)          :: solution_file
         logical, intent(in)       :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: i, j, k, fID
         real(kind=RP)                    :: point(2,3)
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         character(len=STR_LEN_MONITORS)  :: interval, direction
         character(len=STR_LEN_MONITORS)  :: point1_char,point2_char,point3_char
         character(len=STR_LEN_MONITORS)  :: variables
         character(len=STR_LEN_MONITORS)  :: fileFormat
         character(len=STR_LEN_MONITORS)  :: writeInterval
         character(len=STR_LEN_MONITORS)  :: meanData
#if defined(NAVIERSTOKES) && (!(INCNS))
         if (FirstCall) then
!
!           Get monitor ID, assign zero to bufferLine and intervalCount
!           --------------
            self % ID = ID
            self % bufferLine = 0
            self % intervalCount = 0
!
!           Search for the parameters in the case file
!           ------------------------------------------
            write(in_label , '(A,I0)') "#define spatial mean node " , self % ID
         
            call get_command_argument(1, paramFile)
            call readValueInRegion(trim(paramFile), "name"              , self % spatialMeanName  , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "variables"         , variables         , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "direction axis"    , direction         , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "xrange [xmin,xmax]", point1_char       , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "yrange [ymin,ymax]", point2_char       , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "zrange [zmin,zmax]", point3_char       , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "sampling interval" , interval          , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "write interval"    , writeInterval     , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "mean data"         , meanData          , in_label, "#end" )
!
!           Get the variables, points, N discretization, and interval
!           ----------------------------------------------
            call getCharArrayFromString(variables, STR_LEN_MONITORS, self % variable)
            self % nVariables = size(self % variable)
            point(:,1) = getRealArrayFromString(point1_char)
            point(:,2) = getRealArrayFromString(point2_char)
            point(:,3) = getRealArrayFromString(point3_char)
            
            self % pmin(1) = point(1,1)
            self % pmin(2) = point(1,2) 
            self % pmin(3) = point(1,3) 
            
            self % pmax(1) = point(2,1)
            self % pmax(2) = point(2,2) 
            self % pmax(3) = point(2,3)
            
            self % dirAxis  = GetIntValue(direction)
            self % interval = GetIntValue(interval)
            if (len(trim(meanData)).lt.3) then 
                self % meanData = .false.
            else    
                self % meanData = GetLogicalValue(meanData)
            end if 
            
            if ( .not. allocated(mesh % elements(1) % storage % stats % data ) ) then
                self % meanData = .false.
            end if 
            
            interval=TRIM(ADJUSTL(interval))
!
!           Get the max. number of timestep in the buffer file before being written
!           -----------------------------------------------------------------------         
            IF (LEN(TRIM(writeInterval)) .EQ. 0) THEN
                self % bufferSize = 1;
            ELSE
                self % bufferSize = GetIntValue(writeInterval)
!
!               Failsafe to prevent too many data being written at one time
!               -----------------------------------------------------------             
                IF (self % bufferSize .GT. 10000) THEN   
                    self % bufferSize = 10000;
                END IF
            END IF
!
!           Look for unique data point in the range
!           ---------------------------------------
            call self % LookForUniqueCoordinate (mesh)
!
!           Allocate Variables
!           ------------------          
            ALLOCATE(self % fileName (self % nVariables))
            self % values = 0_RP
         end if
         if (self % meanData) then
			allocate(self % meanU(self % nUniqueAll), self % meanV(self % nUniqueAll), self % meanW(self % nUniqueAll))
			self % meanU=0_RP
			self % meanV=0_RP
			self % meanW=0_RP
         end if 
!
!        Check Variables, Create Files, and Write Header Files
!        -----------------------------------------------------
         do i=1,self % nVariables
            if (self % meanData) then
!
!           Prepare the file in which the SpatialMeanNode is exported
!           ---------------------------------------------------------
            write( self % fileName (i) , '(A,A,A,A,A,A,I0,A)') trim(solution_file) ,"_", trim(self % spatialMeanName) ,"_", trim(self % variable(i)) &
                                    , "_spatialTemporalMean_" , self % ID,".node"  
!
!           Check the variable
!           ------------------
            call tolower(self % variable (i))
#ifdef NAVIERSTOKES
			select case ( trim(self % variable (i)) )
            case ("density")
            case ("pressure")
            case ("ptotal")
            case ("velocity")
			case ("viscosity")
            case ("u")
				self % iVarU=i
            case ("v")
				self % iVarV=i
            case ("w")
				self % iVarW=i
            case ("uu")
            case ("vv")
            case ("ww")
            case ("uv")
            case ("uw")
            case ("vw")
			case ("uprime2")
            case ("vprime2")
            case ("wprime2")
            case ("uvprime")
            case ("uwprime")
            case ("vwprime")
            case ("mach")
            case ("k")
            case ("q1")
            case ("q2")
            case ("q3")
            case ("q4")
            case ("q5")
            case default
                if ( MPI_Process % isRoot ) then 
                   print*, 'SpatialMeanNode from temporal mean data, variable "',trim(self % variable(i)),'" not implemented.'
                   print*, "Options available are:"
                   print*, "   * density"
                   print*, "   * pressure"
                   print*, "   * velocity"
				   print*, "   * viscosity"
                   print*, "   * u"
                   print*, "   * v"
                   print*, "   * w"
                   print*, "   * uu"
                   print*, "   * vv"
                   print*, "   * ww"
                   print*, "   * uv"
                   print*, "   * uw"
                   print*, "   * vw"
				   print*, "   * uprime2"
				   print*, "   * vprime2"
				   print*, "   * wprime2"
				   print*, "   * uvprime"
				   print*, "   * uwprime"
				   print*, "   * vwprime"
                   print*, "   * Mach"
                   print*, "   * K"
                   print*, "   * q1"
                   print*, "   * q2"
                   print*, "   * q3"
                   print*, "   * q4"
                   print*, "   * q5"
                end if 
            end select
#endif          
            else
!
!           Prepare the file in which the SpatialMeanNode is exported
!           ---------------------------------------------------------
            write( self % fileName (i) , '(A,A,A,A,A,A,I0,A)') trim(solution_file) ,"_", trim(self % spatialMeanName) ,"_", trim(self % variable(i)) &
                                    , "_spatialmean_" , self % ID,".node"  
!
!           Check the variable
!           ------------------
            call tolower(self % variable (i))

            select case ( trim(self % variable (i)) )
#ifdef NAVIERSTOKES
            case ("density")
            case ("pressure")
            case ("ptotal")
            case ("velocity")
			case ("viscosity")
            case ("u")
            case ("v")
            case ("w")
            case ("mach")
            case ("k")
            case ("omegax")
            case ("omegay")
            case ("omegaz")
            case ("q1")
            case ("q2")
            case ("q3")
            case ("q4")
            case ("q5")
            case default
                if ( MPI_Process % isRoot ) then 
                   print*, 'SpatialMeanNode variable "',trim(self % variable(i)),'" not implemented.'
                   print*, "Options available are:"
                   print*, "   * density"
                   print*, "   * pressure"
                   print*, "   * velocity"
				   print*, "   * viscosity"
                   print*, "   * u"
                   print*, "   * v"
                   print*, "   * w"
                   print*, "   * Mach"
                   print*, "   * K"
                   print*, "   * q1"
                   print*, "   * q2"
                   print*, "   * q3"
                   print*, "   * q4"
                   print*, "   * q5"
                end if 
#endif
#ifdef INCNS
            case default
               print*, "SpatialMeanNodes are not implemented for the incompressible NSE"
#endif
#ifdef MULTIPHASE
            case default
               print*, 'SpatialMeanNodes are not implemented.'
#endif
		 end select
         end if 
		 end do 
		 

         if ( .not. MPI_Process % isRoot ) return

         do i=1,self % nVariables
#if defined(NAVIERSTOKES) && (!(INCNS))
!
!        Create file
!        -----------
         if (FirstCall) then
            if (MPI_Process % isRoot ) then 
                open ( newunit = fID , file = trim(self % fileName (i)) , action = "write" , access = "stream" , status = "replace", position='append' )
!
!           Write the file headers
!           ----------------------
                write( fID) self % ID
                write( fID) sum(self % nMultiplyAll)
                write( fID) self % nUniqueAll
                write( fID) self % dirAxis
                write( fID ) self % interval
                write( fID ) refValues % rho
                write( fID ) refValues % V
                write( fID ) refValues % p
                write( fID ) refValues % T
                write( fID ) refValues % mu
                write( fID ) refValues % AoATheta
                write( fID ) refValues % AoAPhi
                write( fID ) self % geom
                close ( fID )
            end if
         end if
#endif      
        end do 
#if defined(NAVIERSTOKES) && (!(INCNS))
!
!           File Format 
!           -----------------------------------------------
            write( fileFormat , '(A,A,A,A,I0,A)') trim(solution_file) ,"_", trim(self % spatialMeanName) ,"_'variable'_spatialmean_", self % ID, ".node"
!
!        Write Information 
!        -----------------------------------------------         
         write(STD_OUT,'(/)')
         call SubSection_Header("SpatialMean Nodes")
            write(STD_OUT,'(30X,A,A27,I4)') "->" , "SpatialMeanNode ID: " , self % ID
            write(STD_OUT,'(30X,A,A27,A128)') "->" , "Variables: " , variables
            write(STD_OUT,'(30X,A,A27,A,F6.4,A,F6.4,A)') "->" , "xrange [xmin,xmax] (m): ","[", &
                                                   self % pmin(1), ", ", self % pmax(1), "]"
            write(STD_OUT,'(30X,A,A27,A,F6.4,A,F6.4,A)') "->" , "yrange [ymin,ymax] (m): ","[", &
                                                   self % pmin(2), ", ", self % pmax(2), "]"
            write(STD_OUT,'(30X,A,A27,A,F6.4,A,F6.4,A)') "->" , "zrange [zmin,zmax] (m): ","[", &
                                                   self % pmin(3), ", ", self % pmax(3), "]"
            write(STD_OUT,'(30X,A,A27,I5)') "->" , "Number of unique nodes: ", self % nUniqueAll
            write(STD_OUT,'(30X,A,A27,I4)') "->" , "Samplings Interval: ", self % interval
            write(STD_OUT,'(30X,A,A27,I4)') "->" , "Write Interval: ", self % bufferSize
            write(STD_OUT,'(30X,A,A27,L1)') "->" , "Meanflow data: ", self % meanData
            write(STD_OUT,'(30X,A,A27,A128)') "->" , "Filename: ", fileFormat
#endif  

#endif
      end subroutine SpatialMeanNode_Initialization

      subroutine SpatialMeanNode_Update(self, mesh, bufferPosition, t)
         use Physics
         use MPI_Process_Info
         use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
         implicit none
         class(SpatialMeanNode_t)         :: self
         class(HexMesh)                   :: mesh
         integer                          :: bufferPosition
         real(kind=RP)                    :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k, l, m, ierr, eID, pID
         integer, parameter :: NO_OF_VARIABLES_Sij = 9
         integer, parameter ::  U  = 1
         integer, parameter ::  V  = 2
         integer, parameter ::  W  = 3
         integer, parameter ::  UU = 4
         integer, parameter ::  VV = 5
         integer, parameter ::  WW = 6
         integer, parameter ::  UV = 7
         integer, parameter ::  UW = 8
         integer, parameter ::  VW = 9
         real(kind=RP)  :: value, rhoInv, kappa
         real(kind=RP)  , allocatable   :: buff(:,:)
#ifdef NAVIERSTOKES
!
!        Update data based on interval
!        -----------------------------       
         if (self % intervalCount .EQ. 0 ) then
            self % bufferLine = self % bufferLine + 1
            DO m=1,self % nVariables
              self % values(:, bufferPosition, m) = 0_RP
              self % values(1, bufferPosition, m) = t
              DO l=1, self % nActive 
!
!                   Update the Node
!                   ----------------
                    associate( e => mesh % elements(self % activeLoc(1,l)) )
                    associate( Q => e % storage % Q,  S => e % storage, meanQ => e % storage % stats % data)
                    eID=self % activeLoc(1,l)
                    i = self % activeLoc(2,l)
                    j = self % activeLoc(3,l)
                    k = self % activeLoc(4,l)
                    pID=self % activeLoc(5,l)

                select case (trim(self % variable(m)))
                case("density")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHO,i,j,k)
                case("pressure")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Pressure(Q(:,i,j,k)) 
                case("ptotal")
                    value = POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))/POW2(Q(IRHO,i,j,k))     ! Vabs**2
                    value =  value / ( thermodynamics % gamma*(thermodynamics % gamma-1.0_RP)*(Q(IRHOE,i,j,k)/Q(IRHO,i,j,k)-0.5_RP * value) ) ! Mach ^2
                    value = Pressure(Q(:,i,j,k))*(1.0_RP+0.5_RP*(thermodynamics % gamma-1.0_RP)*value)**( thermodynamics % gamma/(thermodynamics % gamma-1.0_RP))
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value
                case("velocity")
                    value = sqrt(POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/Q(IRHO,i,j,k)
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value
			    case("viscosity")
					call get_laminar_mu_kappa(Q(:,i,j,k),value,kappa)
					self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value
                case("omegax")
                    value = (1/Q(IRHO,i,j,k) * S % U_y(IRHOW,i,j,k) - Q(IRHOW,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_y(IRHO,i,j,k)) &
                                    - (1/Q(IRHO,i,j,k) * S % U_z(IRHOV,i,j,k) - Q(IRHOV,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_z(IRHO,i,j,k))
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value
                case("omegay")
                    value = (1/Q(IRHO,i,j,k) * S % U_z(IRHOU,i,j,k) - Q(IRHOU,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_z(IRHO,i,j,k)) &
                                    - (1/Q(IRHO,i,j,k) * S % U_x(IRHOW,i,j,k) - Q(IRHOW,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_x(IRHO,i,j,k))
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value               
                case("omegaz")
                    value = (1/Q(IRHO,i,j,k) * S % U_x(IRHOV,i,j,k) - Q(IRHOV,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_x(IRHO,i,j,k)) &
                                    - (1/Q(IRHO,i,j,k) * S % U_y(IRHOU,i,j,k) - Q(IRHOU,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_y(IRHO,i,j,k))
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value                   
                case("u")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
                case("v")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
                case("w")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
                case("mach")
                    value = POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))/POW2(Q(IRHO,i,j,k))     ! Vabs**2
                    value = sqrt( value / ( thermodynamics % gamma*(thermodynamics % gamma-1.0_RP)*(Q(IRHOE,i,j,k)/Q(IRHO,i,j,k)-0.5_RP * value) ) )
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value       
                case("k")
                    value = 0.5_RP * (POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/Q(IRHO,i,j,k)
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ value       
                case("q1")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHO,i,j,k)   
                case("q2")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOU,i,j,k)  
                case("q3")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOV,i,j,k)  
                case("q4")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOW,i,j,k)   
                case("q5")
                    self % values ( pID+1,bufferPosition, m) = self % values (pID+1,bufferPosition, m)+ Q(IRHOE,i,j,k)  
                end select
                    end associate
                    end associate
              END DO
!           Combine all result from all proc MPI into root
!           ----------------------------------------------          
#ifdef _HAS_MPI_            
                call mpi_barrier(MPI_COMM_WORLD, ierr)  
                ALLOCATE(buff(self % nUniqueAll,MPI_Process % nProcs))
                call mpi_allgather(self % values(2:self % nUniqueAll+1, bufferPosition, m), self % nUniqueAll, MPI_DOUBLE, buff, self % nUniqueAll, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
                self % values(2:self % nUniqueAll+1, bufferPosition, m) = SUM(buff,DIM=2)
                DEALLOCATE(buff)
#endif
!           Average with the number of data
!           -------------------------------
                do i = 1, self % nUniqueAll
                    self % values(i+1, bufferPosition, m) = self % values(i+1, bufferPosition, m) / self % nMultiplyAll(i)
                end do 
            end do
			if (self % meanData) then
				self % meanU = self % values (2: self % nUniqueAll+1,bufferPosition,self % iVarU)
				self % meanV = self % values (2: self % nUniqueAll+1,bufferPosition,self % iVarV)
				self % meanW = self % values (2: self % nUniqueAll+1,bufferPosition,self % iVarW)
			end if 
        end if
        
         self % intervalCount = self % intervalCount + 1
         
         if (self % intervalCount .EQ. self % interval) then
            self % intervalCount = 0 
         end if 
#endif
      end subroutine SpatialMeanNode_Update

      subroutine SpatialMeanNode_WriteToFile ( self, no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(SpatialMeanNode_t)  :: self
         integer                   :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i, j, ierr
         integer                    :: fID
        
         if ( MPI_Process % isRoot ) then
         DO j=1, self % nVariables
            open( newunit = fID , file = trim ( self % fileName (j) ) , action = "write" , access = "stream" , status = "old", position='append' )
            do i = 1 , no_of_lines
                write( fID ) self % values(:,i,j)
            end do
            close ( fID )
         END DO
         end if
         if ( no_of_lines .ne. 0 ) self % values(:,1,:) = self % values(:,no_of_lines,:)
         
         self % bufferLine = 0
         
         self % values = 0_RP
      
      end subroutine SpatialMeanNode_WriteToFile

      subroutine SpatialMeanNode_LookForUniqueCoordinate(self, mesh)
         use MPI_Process_Info
         implicit none
         class(SpatialMeanNode_t)       :: self
         class(HexMesh)                 :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: eID, i, j, k ,l , nPotentialUnique, nUnique, nCount
         integer, allocatable       :: uniqueProcs(:), activeLoc(:,:)
         integer                    :: ierr
         real(kind=RP), allocatable :: yUnique(:), yFinal(:), yUniqueProc(:,:)
         real(kind=RP)              :: error=0.000001
         real(kind=RP)              :: buff, y
         logical                    :: unique, finish
!
!        Obtain the direction at which spatial mean is not applied and store the coordinate
!        ----------------------------------------------------------------------------------
         nPotentialUnique=0
         do eID=1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
            nPotentialUnique=nPotentialUnique+(e % Nxyz(self % dirAxis) +1) * (e % Nxyz(1) +1) !nPotentialUnique Points
            end associate
         end do 
         ALLOCATE(activeLoc(1:5, mesh % NDOF), yUnique(nPotentialUnique))
         nUnique=0
         nCount =0
         do eID=1, mesh % no_of_elements
                if (.not.((any(mesh % elements(eID) % geom % x(1,:,:,:).gt.(self % pmin(1)-error))) & 
                            .and.(any(mesh % elements(eID) % geom % x(1,:,:,:).lt.(self % pmax(1)+error))))) then
                    cycle
                end if 
                if (.not.((any(mesh % elements(eID) % geom % x(2,:,:,:).gt.(self % pmin(2)-error))) & 
                            .and.(any(mesh % elements(eID) % geom % x(2,:,:,:).lt.(self % pmax(2)+error))))) then
                    cycle
                end if 
                if (.not.((any(mesh % elements(eID) % geom % x(3,:,:,:).gt.(self % pmin(3)-error))) & 
                            .and.(any(mesh % elements(eID) % geom % x(3,:,:,:).lt.(self % pmax(3)+error))))) then
                    cycle
                end if 
                
                associate(e => mesh % elements(eID))
                do i=0,e % Nxyz(3); do j=0,e % Nxyz(2); do k=0,e % Nxyz(1)
                    if (.not.(((e % geom % x(1,k,j,i)).gt.(self % pmin(1)-error)) &
                         .and.((e % geom % x(1,k,j,i)).lt.(self % pmax(1)+error)))) cycle
                    if (.not.(((e % geom % x(2,k,j,i)).gt.(self % pmin(2)-error)) &
                         .and.((e % geom % x(2,k,j,i)).lt.(self % pmax(2)+error)))) cycle
                    if (.not.(((e % geom % x(3,k,j,i)).gt.(self % pmin(3)-error)) &
                         .and.((e % geom % x(3,k,j,i)).lt.(self % pmax(3)+error)))) cycle
!
!                   Store the location of the active node ( within the bounded x,y,z ) 
!                   ------------------------------------------------------------------
                    nCount = nCount+1
                    activeLoc(1, nCount) = eID
                    activeLoc(2, nCount) = k
                    activeLoc(3, nCount) = j
                    activeLoc(4, nCount) = i
!
!                   Find the unique points
!                   ----------------------                  
                    y = e % geom % x( self % dirAxis, k,j,i)
                    unique=.true.
                    do l=1,nUnique
                        if (abs(yUnique(l)-y).lt.error) then
                            unique=.false.
                            exit
                        end if
                    end do 
                    
                    if (unique) then
                        nUnique = nUnique + 1
                        yUnique(nUnique)=y
                    end if
                end do ; end do ; end do 
            end associate
         end do 
         
         ALLOCATE(yFinal(nUnique), self % activeLoc(5,nCount))
         self % nActive = nCount
         yFinal=yUnique(1:nUnique)
         self % activeLoc(1:4,:) = activeLoc(1:4,1:nCount)
         DEALLOCATE(yUnique, activeLoc)
!
!        MPI Operation to combine operation from different MPI
!        -----------------------------------------------------
#ifdef _HAS_MPI_
            ALLOCATE(uniqueProcs(MPI_Process % nProcs))
!
!           Gather all data from all processes
!           ----------------------------------
            call mpi_allgather(nUnique, 1, MPI_INT, uniqueProcs, 1, MPI_INT, MPI_COMM_WORLD, ierr)
            nPotentialUnique = sum(uniqueProcs)
            call mpi_barrier(MPI_COMM_WORLD, ierr) 
!
!           Send unique coordinate on each proc to every proc with allgather
!           ----------------------------------------------------------------
            ALLOCATE(yUnique(maxval(uniqueProcs)))
            ALLOCATE(yUniqueProc(maxval(uniqueProcs),MPI_Process % nProcs))
            yUnique = 0_RP
            yUnique(1:nUnique)=yFinal
            call mpi_allgather(yUnique, maxval(uniqueProcs), MPI_DOUBLE, yUniqueProc, maxval(uniqueProcs), MPI_DOUBLE, MPI_COMM_WORLD, ierr)
            call mpi_barrier(MPI_COMM_WORLD, ierr) 
            
            nCount = uniqueProcs(1)
            DEALLOCATE(yUnique)
            ALLOCATE(yUnique(nPotentialUnique))
            yUnique(1:nCount) = yUniqueProc(1:nCount,1)
            
            do i=2, MPI_Process % nProcs
                do j=1, uniqueProcs(i)
                    unique=.true. 
                    if (any((yUnique-yUniqueProc(j,i)).lt.error)) then
                        unique = .false.
                        cycle
                    end if
                    if (unique)then
                        nCount = nCount+1
                        yUnique(nCount)=yUniqueProc(j,i)
                    end if 
                end do 
            end do

            DEALLOCATE(yFinal, yUniqueProc)
            ALLOCATE(yFinal(1:nCount))
            yFinal=yUnique(1:nCount)
            DEALLOCATE (yUnique)
!
!           Sort the coordinate in ascending order
!           --------------------------------------          
            CALL sortDoubleArrayMinMax(nCount, yFinal, 0)
            self % nUniqueAll = nCount

            call mpi_barrier(MPI_COMM_WORLD, ierr) 
#else           
            self % nUniqueAll = nUnique
#endif          
            ALLOCATE(self % values (self % nUniqueAll+1, self % bufferSize, self % nVariables), self % geom (self % nUniqueAll), &
                    self % nMultiply ( self % nUniqueAll), self % nMultiplyAll ( self % nUniqueAll))
            self % geom= yFinal
!
!           Assign location of all active node w.r.t. unique coordinate location
!           --------------------------------------------------------------------
            self % nMultiply = 0
            do i=1, self % nActive
                eID = self % activeLoc(1,i)
                do j= 1, self % nUniqueAll
                    buff = mesh % elements(eID) % geom % x (self % dirAxis, self % activeLoc(2,i), self % activeLoc(3,i), self % activeLoc(4,i))
                    if (abs(buff - self % geom (j)).lt.error) then
                        self % nMultiply(j) = self % nMultiply(j) +1 ! Local for each proc
                        self % activeLoc(5,i) = j
                        exit
                    end if 
                end do 
            end do 
#ifdef _HAS_MPI_
!
!           Gather nMultiply to be combined into nMultiplyAll in root
!           ---------------------------------------------------------
            self % nMultiplyAll = 0
            ALLOCATE(activeLoc(self % nUniqueAll, MPI_Process % nProcs))
            call mpi_allgather(self % nMultiply, self % nUniqueAll, MPI_INT, activeLoc, self % nUniqueAll, MPI_INT, MPI_COMM_WORLD, ierr)
            
            self % nMultiplyAll = sum(activeLoc, DIM=2)
            
            DEALLOCATE(activeLoc)
            
#else   
            self % nMultiplyAll = self % nMultiply

#endif
      end subroutine SpatialMeanNode_LookForUniqueCoordinate
!
!////////////////////////////////////////////////////////////////////////
!
! This Recursive Subroutine sort Double ArraySet at nColumn from its Minimum to Maximum Value
! When call for the first time k must be set to 0 -- WARNING DO NOT USE FOR HUGE ARRAY
!
        RECURSIVE SUBROUTINE sortDoubleArrayMinMax(sizeArray, ArraySet, k)
            IMPLICIT NONE
            INTEGER                                             ,INTENT(IN)       :: sizeArray
            REAL(kind=RP)   ,DIMENSION(1:sizeArray)             ,INTENT(INOUT)    :: ArraySet
            INTEGER                                             ,INTENT(IN)       :: k
!
!        ---------------
!        Local variables
!        ---------------
!
            REAL(kind=RP)     :: BufferSet
            INTEGER           :: i
!
!        ---------------
!        Perform Sorting Value
!        ---------------
!
            SELECT CASE(k)
            CASE(0)
                DO i=1,sizeArray-1
                    IF (ArraySet(i).GT.ArraySet(i+1)) THEN
                        BufferSet     = ArraySet(i)
                        ArraySet(i)   = ArraySet(i+1)
                        ArraySet(i+1) = BufferSet
                        IF(i.GT.1) THEN
                            CALL sortDoubleArrayMinMax(sizeArray, ArraySet, i)
                        END IF
                    END IF
                END DO
            CASE DEFAULT ! For Recursive
                IF (ArraySet(k-1).GT.ArraySet(k)) THEN
                    BufferSet     = ArraySet(k-1)
                    ArraySet(k-1) = ArraySet(k)
                    ArraySet(k)   = BufferSet
                    IF(k.GT.2) THEN
                        CALL sortDoubleArrayMinMax(sizeArray, ArraySet, k-1)
                    END IF
                END IF
            END SELECT

        END SUBROUTINE sortDoubleArrayMinMax
      
      elemental subroutine SpatialMeanNode_Destruct (self)
         implicit none
         class(SpatialMeanNode_t), intent(inout) :: self
         
         safedeallocate (self % activeLoc)
         safedeallocate (self % nMultiply)
         safedeallocate (self % nMultiplyAll)
         safedeallocate (self % geom)
		 safedeallocate (self % values)
		 safedeallocate (self % meanU)
		 safedeallocate (self % meanV)
         safedeallocate (self % meanW)
         safedeallocate (self % fileName)
         safedeallocate (self % variable)

      end subroutine SpatialMeanNode_Destruct
      
      elemental subroutine SpatialMeanNode_Assign (to, from)
         implicit none
         class(SpatialMeanNode_t), intent(inout) :: to
         type(SpatialMeanNode_t) , intent(in) :: from
        
         
         to % ID = from %  ID
         to % nVariables      = from % nVariables
         to % interval        = from % interval
         to % bufferSize      = from % bufferSize
         to % bufferLine      = from % bufferLine
         to % intervalCount   = from % intervalCount    
         to % nActive         = from % nActive
         to % dirAxis         = from % dirAxis
         to % nUniqueAll      = from % nUniqueAll   
         to % pmin            = from % pmin
         to % pmax            = from % pmax
         to % error           = from % error    
         
         safedeallocate ( to % activeLoc )
         allocate ( to % activeLoc( size(from % activeLoc,1),size(from % activeLoc,2) ) )
         to % activeLoc = from % activeLoc
         
         safedeallocate ( to % nMultiply )
         allocate ( to % nMultiply( size(from % nMultiply) ) )
         to % nMultiply = from % nMultiply
         
         safedeallocate ( to % nMultiplyAll )
         allocate ( to % nMultiplyAll( size(from % nMultiplyAll) ) )
         to % nMultiplyAll = from % nMultiplyAll
         
         safedeallocate ( to % geom )
         allocate ( to % geom( size(from % geom) ) )
         to % geom = from % geom
         
         safedeallocate ( to % values )
         allocate ( to % values( size(from % values,1),size(from % values,2),size(from % values,3) ) )
         to % values = from % values

         safedeallocate ( to % fileName )
         allocate ( to % fileName ( size(from % fileName) ) )
         to % fileName = from % fileName
         
         to % spatialMeanName = from % spatialMeanName

         safedeallocate ( to % variable )
         allocate ( to % variable ( size(from % variable) ) )
         to % variable = from % variable         

         
         
      end subroutine SpatialMeanNode_Assign
      
end module SpatialMeanNode