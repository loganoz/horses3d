#include "Includes.h"

module PlaneSampling
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use VariableConversion
   use MPI_Process_Info
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString, GetRealValue, getCharArrayFromString
   use NodalStorageClass   , only: NodalStorage
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
   private
   public   PlaneSampling_t
!
!  *******************************
!  PlaneSamplings class definition
!  *******************************
!

   type PlaneSampling_t
      logical      , allocatable      :: active(:)
      integer      , allocatable      :: rank (:)
      integer                         :: ID
      integer                         :: nVariables
      integer                         :: interval
      integer                         :: bufferSize
      integer                         :: bufferLine
      integer                         :: intervalCount
      integer                         :: N(2)
      integer                         :: nNodes
      integer      , allocatable      :: eID (:)
      real(kind=RP), allocatable      :: lxi(:,:) , leta(:,:), lzeta(:,:)
      real(kind=RP), allocatable      :: values(:,:,:)
      real(kind=RP), allocatable      :: x(:,:)
      real(kind=RP), allocatable      :: xi(:,:)
      logical                         :: disturbanceData =.false.
      character(len=STR_LEN_MONITORS),  allocatable :: fileName (:)
      character(len=STR_LEN_MONITORS)               :: planeName
      character(len=STR_LEN_MONITORS)               :: fileInput
      character(len=STR_LEN_MONITORS),  allocatable :: variable (:)
      contains
         procedure   :: Initialization => Plane_Initialization
         procedure   :: Update         => Plane_Update
		 procedure   :: UpdateInterp   => Plane_UpdateLagrangeInterp
         procedure   :: WriteToFile    => Plane_WriteToFile
         procedure   :: LookInOtherPartitions => Plane_LookInOtherPartitions
         procedure   :: destruct       => Plane_Destruct
         procedure   :: copy           => Plane_Assign
         generic     :: assignment(=)  => copy
   end type PlaneSampling_t

   contains

      subroutine Plane_Initialization(self, mesh, ID, solution_file, FirstCall)
         use Headers
         use ParamfileRegions
         use MPI_Process_Info
         use Utilities, only: toLower
         implicit none
         class(PlaneSampling_t)  :: self
         class(HexMesh)          :: mesh
         integer                 :: ID
         character(len=*)        :: solution_file
         logical, intent(in)     :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: i, j, k, l, fid, nNodes,inputType, ierr
         logical                          :: firstOutDomain
         real(kind=RP)                    :: point_1(NDIM,1),point_2(NDIM,1),point_3(NDIM,1)
         real(kind=RP), allocatable       :: pStart(:,:), pEnd(:,:)
         real(kind=RP)                    :: x(NDIM), xi(NDIM)
         real(kind=RP)                    :: delta(NDIM,2), invNodes(2)
         real(kind=RP), allocatable       :: deltaArray1(:,:)
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         character(len=STR_LEN_MONITORS)  :: interval
         character(len=STR_LEN_MONITORS)  :: point1_char,point2_char,point3_char
         character(len=STR_LEN_MONITORS)  :: N12,N23, inputType_char
         character(len=STR_LEN_MONITORS)  :: variables
         character(len=STR_LEN_MONITORS)  :: fileFormat
         character(len=STR_LEN_MONITORS)  :: writeInterval
         character(len=STR_LEN_MONITORS)  :: distData
         
         firstOutDomain = .FALSE.
         
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
            write(in_label , '(A,I0)') "#define plane sampling " , self % ID
            
            call get_command_argument(1, paramFile)
            call readValueInRegion(trim(paramFile), "name"           , self % planeName  , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "variables"      , variables         , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "node input type", inputType_char    , in_label, "#end" )
			if (len(trim(inputType_char))>0) then
				inputType  = GetRealValue(inputType_char)
			else
				inputType  = 0
			end if 
            if (inputType.ne.0) then
                call readValueInRegion(trim(paramFile), "input file"           , self % fileInput  , in_label, "#end" )
            else
                call readValueInRegion(trim(paramFile), "point 1 [x,y,z]", point1_char       , in_label, "#end" )
                call readValueInRegion(trim(paramFile), "point 2 [x,y,z]", point2_char       , in_label, "#end" )
                call readValueInRegion(trim(paramFile), "point 3 [x,y,z]", point3_char       , in_label, "#end" )
                call readValueInRegion(trim(paramFile), "n12"            , N12               , in_label, "#end" )
                call readValueInRegion(trim(paramFile), "n23"            , N23               , in_label, "#end" )
            end if 
            call readValueInRegion(trim(paramFile), "sampling interval", interval        , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "write interval", writeInterval      , in_label, "#end" )
            call readValueInRegion(trim(paramFile), "data type", distData, in_label, "#end" )
            
!
!           Get the variables, points, N discretization, and interval
!           ----------------------------------------------
            call getCharArrayFromString(variables, STR_LEN_MONITORS, self % variable)
            self % nVariables = size(self % variable)
            if (inputType.eq.0) then
                point_1(:,1) = getRealArrayFromString(point1_char)
                point_2(:,1) = getRealArrayFromString(point2_char)
                point_3(:,1) = getRealArrayFromString(point3_char)
                self % N(1)  = GetRealValue(N12)
                self % N(2)  = GetRealValue(N23)
            else 
!               Open file input for nodes
                open(newunit=fid, file= trim(self % fileInput),status="old",action="read", iostat=ierr)
                if (ierr /= 0) then
                    write(*,*) "ERROR-opening input file for plane sampling: ", trim(self % fileInput)
                    stop
                end if
                read(fid, *, iostat=ierr) self % N(1), self % N(2)
            end if 
            self % interval = GetRealValue(interval)
            
            
            interval=TRIM(ADJUSTL(interval))
            writeInterval=TRIM(ADJUSTL(writeInterval))

!
!           Failsafe if the number of node is less than 2
!           ---------------------------------------------
            IF ((self % N(1).LT.2)) THEN
                self % N(1)=2
            END IF
            IF ((self % N(2).LT.2)) THEN
                self % N(2)=2
            END IF
            
            self % nNodes = self % N(1) * self % N(2)
            nNodes = self % nNodes
!
!           Get the max. number of timestep in the buffer file before being written
!           -----------------------------------------------------------------------         
            IF (LEN(TRIM(writeInterval)) .EQ. 0) THEN
                self % bufferSize = 1;
            ELSE
                self % bufferSize = GetRealValue(writeInterval)
!
!               Failsafe to prevent too many data being written at one time
!               -----------------------------------------------------------             
                IF (nNodes * self % bufferSize .GT. 500000) THEN   
                    self % bufferSize = 1;
                END IF
            END IF
!
!           Allocate Variables
!           ------------------          
            ALLOCATE(self % x (NDIM, nNodes), self % xi (NDIM, nNodes), self % active (nNodes) &
                    , self % eID (nNodes), self % rank (nNodes), self % fileName (self % nVariables) &
                    , self % values(nNodes+1,self % bufferSize,self % nVariables), deltaArray1(NDIM,self % N(1))  &
                    , pStart(NDIM, self % N(1)))
!
!           Allocate Lagrange interpolants - Assumed identical polynomial for all elements
!           -----------------------------------------------------------------------------
            ALLOCATE(self % lxi( 0 : mesh % elements (1) % Nxyz(1), 1:nNodes), &
                        self % leta( 0 : mesh % elements (1) % Nxyz(2), 1:nNodes), &
                        self % lzeta( 0 : mesh % elements (1) % Nxyz(3), 1:nNodes))
!
!           Nodes spacing
!           -------------           
            invNodes(1)=1.0_RP/(self % N(1)-1)
            invNodes(2)=1.0_RP/(self % N(2)-1)
            if (inputType.eq.0) then 
                delta(1,1)=(point_2(1,1)-point_1(1,1))*invNodes(1)
                delta(2,1)=(point_2(2,1)-point_1(2,1))*invNodes(1)
                delta(3,1)=(point_2(3,1)-point_1(3,1))*invNodes(1)
                delta(1,2)=(point_3(1,1)-point_2(1,1))*invNodes(2)
                delta(2,2)=(point_3(2,1)-point_2(2,1))*invNodes(2)
                delta(3,2)=(point_3(3,1)-point_2(3,1))*invNodes(2)
!
!           Create spacing Array for Nodes in direction 1
!           ---------------------------------------------               
                DO l=1,self % N(1)
                    deltaArray1(1,l)=(delta(1,1))*(l-1)
                    deltaArray1(2,l)=(delta(2,1))*(l-1)
                    deltaArray1(3,l)=(delta(3,1))*(l-1)
                END DO
!
!           Get Nodes location of the plane
!           -------------------------------             
                DO l=1,self % N(2) 
                    DO k=1, self % N(1)
                        DO j=1,3
                            pStart(j,k) = (point_1(j,1)+delta(j,2)*(l-1))
                            self % x(j,(l-1)*self % N(1)+k) = pStart(j,k) + deltaArray1(j,k)
                        END DO
                    END DO 
                END DO
            else  
                ALLOCATE(pEnd(NDIM, self % N(1)))
!
!           Read first points location from input file              
!           ----------------------------------      
                DO l=1,self % N(1)
                    read(fid, *, iostat=ierr) pStart(1,l), pStart(2,l), pStart(3,l)
                    if (ierr /= 0) exit
                END DO 
!
!           Read last points location from input file               
!           ----------------------------------      
                DO l=1,self % N(1)
                    read(fid, *, iostat=ierr) pEnd(1,l), pEnd(2,l), pEnd(3,l)
                    if (ierr /= 0) exit
                END DO 
                close(fid)
!
!           Create spacing Array for Nodes in sweep direction
!           -------------------------------------------------               
                DO l=1,self % N(1)
                    deltaArray1(1,l)=(pEnd(1,l)-pStart(1,l))*invNodes(2)
                    deltaArray1(2,l)=(pEnd(2,l)-pStart(2,l))*invNodes(2)
                    deltaArray1(3,l)=(pEnd(3,l)-pStart(3,l))*invNodes(2)
                END DO
!
!           Get Nodes location of the plane
!           -------------------------------             
                DO l=1,self % N(2) 
                    DO j=1,NDIM
                        DO k=1, self % N(1)
                            self % x(j,(l-1)*self % N(1)+k) = pStart(j,k) + deltaArray1(j,k)*(l-1)
                        END DO
                    END DO 
                END DO       
                DEALLOCATE(pEnd)				
            end if 

            
            DEALLOCATE (deltaArray1,pStart)

!
!           Find node location in the element ID
!           ------------------------------------                
            DO i=1,self % nNodes
        
                x(1)=self % x(1,i)
                x(2)=self % x(2,i)
                x(3)=self % x(3,i)
!
!               Find the requested point in the mesh
!               ------------------------------------
                self % active (i) = mesh % FindPointWithCoords(x, self % eID(i), xi)
                
                self % xi(1,i) = xi(1)
                self % xi(2,i) = xi(2)
                self % xi(3,i) = xi(3)
!
!               Check whether the Plane is located in other partition
!               -----------------------------------------------------
                call self % LookInOtherPartitions (i)
!
!               Disable the nodes if the point is not found
!               -------------------------------------------
                IF ( .not. self % active (i) ) then
                    IF ( .not. firstOutDomain) then
                        IF ( MPI_Process % isRoot ) then
                            firstOutDomain = .TRUE.
                        END IF
                    END IF
                END IF
!
!               Get the Lagrange interpolants
!               -----------------------------               
                if ( (MPI_Process % rank .ne. self % rank (i)).OR.( .not. self % active (i) ) ) then
                    self % lxi (:,i)= 0.0_RP
                    self % leta (:,i)= 0.0_RP
                    self % lzeta (:,i)= 0.0_RP
                ELSE
                    associate(e => mesh % elements(self % eID(i)))
                    associate( spAxi   => NodalStorage(e % Nxyz(1)), &
                                spAeta  => NodalStorage(e % Nxyz(2)), &
                                spAzeta => NodalStorage(e % Nxyz(3)) )
                    self % lxi  (:,i)  = spAxi % lj(self % xi(1,i))
                    self % leta  (:,i) = spAeta % lj(self % xi(2,i))
                    self % lzeta  (:,i)= spAzeta % lj(self % xi(3,i))
                    end associate
                    end associate
                END IF
                
            END DO
            
         end if
!
!        Check Variables, Create Files, and Write Header Files
!        -----------------------------------------------------
         do i=1,self % nVariables
!
!           Prepare the file in which the Plane is exported
!           -----------------------------------------------
            write( self % fileName (i) , '(A,A,A,A,A,A,I0,A)') trim(solution_file) ,"_", trim(self % planeName) ,"_", trim(self % variable(i)) &
                                    , "_plane_" , self % ID,".sampling"  
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
                   print*, 'Plane variable "',trim(self % variable(i)),'" not implemented.'
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
               print*, "Planes are not implemented for the incompressible NSE"
#endif
#ifdef MULTIPHASE
            case ("density")
            case ("pressure")
			case ("static-pressure")
            case ("concentration")
            case ("velocity")
            case ("u")
            case ("v")
            case ("w")
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
                   print*, 'Plane variable "',trim(self % variable(i)),'" not implemented.'
                   print*, "Options available are:"
                   print*, "   * density"
                   print*, "   * pressure"
				   print*, "   * static-pressure"
                   print*, "   * velocity"
                   print*, "   * concentration"
                   print*, "   * u"
                   print*, "   * v"
                   print*, "   * w"
                   print*, "   * omegax"
				   print*, "   * omegay"
				   print*, "   * omegaz"
                   print*, "   * q1"
                   print*, "   * q2"
                   print*, "   * q3"
                   print*, "   * q4"
                   print*, "   * q5"
                end if 
#endif
		end select

        if ( .not. MPI_Process % isRoot ) return


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
                write( fID) self % N(1)
                write( fID) self % N(2)
                write( fID) self % nNodes
                write( fID ) self % interval
#if defined(NAVIERSTOKES) && (!(INCNS))
                write( fID ) refValues % rho
                write( fID ) refValues % V
                write( fID ) refValues % p
                write( fID ) refValues % T
                write( fID ) refValues % mu
                write( fID ) refValues % AoATheta
                write( fID ) refValues % AoAPhi
#elif defined(MULTIPHASE)
                write( fID ) refValues % rho
                write( fID ) refValues % V
                write( fID ) refValues % p
                write( fID ) thermodynamics % rho(2)
                write( fID ) thermodynamics % mu(1)
                write( fID ) thermodynamics % mu(2)
                write( fID ) refValues % g0
#endif   
                write( fID ) transpose(self % x)
                close ( fID )
            end if
         end if
   
        end do 
!
!           File Format 
!           -----------------------------------------------
            write( fileFormat , '(A,A,A,A,I0,A)') trim(solution_file) ,"_", trim(self % planeName) ,"_'variable'_plane_", self % ID, ".sampling"
!
!        Write Information 
!        -----------------------------------------------         
         write(STD_OUT,'(/)')
         call SubSection_Header("Plane Samplings")
            write(STD_OUT,'(30X,A,A27,I4)') "->" , "Plane ID: " , self % ID
            write(STD_OUT,'(30X,A,A27,A128)') "->" , "Variables: " , variables
			if (inputType.eq.0) then
			    write(STD_OUT,'(30X,A,A27,A36)') "->" , "Input Type: " , "Specified 3 points location"
				write(STD_OUT,'(30X,A,A27,A,F7.2,A,F7.2,A,F7.2,A)') "->" , "Point x1 [m]: ","[", &
													   point_1(1,1), ", ", point_1(2,1), ", ", point_1(3,1), "]"
				write(STD_OUT,'(30X,A,A27,A,F7.2,A,F7.2,A,F7.2,A)') "->" , "Point x2 [m]: ","[", &
													   point_2(1,1), ", ", point_2(2,1), ", ", point_2(3,1), "]"
				write(STD_OUT,'(30X,A,A27,A,F7.2,A,F7.2,A,F7.2,A)') "->" , "Point x3 [m]: ","[", &
													   point_3(1,1), ", ", point_3(2,1), ", ", point_3(3,1), "]"
	        else
			    write(STD_OUT,'(30X,A,A27,A27)') "->" , "Input Type: " , "Specified Input File"
				write(STD_OUT,'(30X,A,A27,A27)') "->" , "Input File: " , self % fileInput
				write(STD_OUT,'(30X,A,A27,A,F7.2,A,F7.2,A,F7.2,A)') "->" , "First Point [m]: ","[", &
													   self % x(1,1), ", ", self % x(2,1), ", ", self % x(3,1), "]"
				write(STD_OUT,'(30X,A,A27,A,F7.2,A,F7.2,A,F7.2,A)') "->" , "Last Point [m]: ","[", &
													   self % x(1,self % nNodes), ", ", self % x(2,self % nNodes), ", ", self % x(3,self % nNodes), "]"
			end if 
            write(STD_OUT,'(30X,A,A27,A,I4,A,I4,A)') "->" , "Nodes [N12,N23]: ","[", &
                                                   self % N(1), ", ", self % N(2), "]"
            write(STD_OUT,'(30X,A,A27,I4)') "->" , "Samplings Interval: ", self % interval
            write(STD_OUT,'(30X,A,A27,A128)') "->" , "Filename: ", fileFormat
            if (firstOutDomain) then    
                write(STD_OUT,'(30X,A,A27,A61)') "->" ,"Note: ","Node/Nodes are located out of domain. A NaN will be assigned"
            end if
    
      end subroutine Plane_Initialization
!
!     Update Lagrange Interpolants After pAdaptation
!     -----------------------------------------------------
      subroutine Plane_UpdateLagrangeInterp(self, mesh)
         use Physics
         use MPI_Process_Info
         use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
         implicit none
         class(PlaneSampling_t)           :: self
         class(HexMesh)                   :: mesh
	  
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: i
!
!           Find node location in the element ID
!           ------------------------------------                
            DO i=1,self % nNodes
!
!               Get the Lagrange interpolants
!               -----------------------------               
                if ( (MPI_Process % rank .ne. self % rank (i)).OR.( .not. self % active (i) ) ) then
                    self % lxi (:,i)= 0.0_RP
                    self % leta (:,i)= 0.0_RP
                    self % lzeta (:,i)= 0.0_RP
                ELSE
                    associate(e => mesh % elements(self % eID(i)))
                    associate( spAxi   => NodalStorage(e % Nxyz(1)), &
                                spAeta  => NodalStorage(e % Nxyz(2)), &
                                spAzeta => NodalStorage(e % Nxyz(3)) )
                    self % lxi  (:,i)  = spAxi % lj(self % xi(1,i))
                    self % leta  (:,i) = spAeta % lj(self % xi(2,i))
                    self % lzeta  (:,i)= spAzeta % lj(self % xi(3,i))
                    end associate
                    end associate
                END IF
                
            END DO
	  end subroutine Plane_UpdateLagrangeInterp

      subroutine Plane_Update(self, mesh, bufferPosition, t)
         use Physics
         use MPI_Process_Info
         use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
         implicit none
         class(PlaneSampling_t)           :: self
         class(HexMesh)                   :: mesh
         integer                          :: bufferPosition
         real(kind=RP)                    :: t
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k, l, m, ierr
         real(kind=RP)  :: value, kappa
		 real(kind=RP)  :: c, sqrtRho
         real(kind=RP)  :: NaN
         real(kind=RP)  :: Sym, Asym
         real(kind=RP)  :: invRhoBase, invRhoF
         real(kind=RP), allocatable   :: var(:,:,:)   
                              
         if (self % intervalCount .EQ. 0 ) then
            self % bufferLine = self % bufferLine + 1

            DO m=1,self % nVariables
            
              self % values(1, bufferPosition, m) = t
              
              DO l=1, self % nNodes 
!
!               Assign NaN if the node located outside of domain
!               ------------------------------------------------              
                IF ( .not. self % active (l) ) THEN  
                    self % values(l, bufferPosition,m) = IEEE_VALUE(NaN, IEEE_QUIET_NAN)
                ELSE
                    if ( MPI_Process % rank .eq. self % rank (l) ) then
!
!                   Update the Node
!                   ----------------
                    associate( e => mesh % elements(self % eID(l)) )
					allocate (var(0:e % Nxyz(1),0:e % Nxyz(2),0:e % Nxyz(3)  ))	 
                    associate( Q => e % storage % Q,  S => e % storage )
   
#ifdef NAVIERSTOKES
            select case (trim(self % variable(m)))
            case("density")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHO,i,j,k)
               end do            ; end do             ; end do
            case("pressure")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Pressure(Q(:,i,j,k))
               end do            ; end do             ; end do
               
            case("viscosity")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  call get_laminar_mu_kappa(Q(:,i,j,k),var(i,j,k),kappa)
               end do            ; end do             ; end do
               
            case("ptotal")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))/POW2(Q(IRHO,i,j,k))     ! Vabs**2
                  var(i,j,k) =  var(i,j,k) / ( thermodynamics % gamma*(thermodynamics % gamma-1.0_RP)*(Q(IRHOE,i,j,k)/Q(IRHO,i,j,k)-0.5_RP * var(i,j,k)) ) ! Mach ^2
                  var(i,j,k) = Pressure(Q(:,i,j,k))*(1.0_RP+0.5_RP*(thermodynamics % gamma-1.0_RP)*var(i,j,k))**( thermodynamics % gamma/(thermodynamics % gamma-1.0_RP))
               end do            ; end do             ; end do
   
            case("velocity")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  var(i,j,k) = sqrt(POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/Q(IRHO,i,j,k)
               end do         ; end do         ; end do
            
            case("omegax")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = (1/Q(IRHO,i,j,k) * S % U_y(IRHOW,i,j,k) - Q(IRHOW,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_y(IRHO,i,j,k)) &
                                - (1/Q(IRHO,i,j,k) * S % U_z(IRHOV,i,j,k) - Q(IRHOV,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_z(IRHO,i,j,k))
               end do            ; end do             ; end do
               
            case("omegay")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = (1/Q(IRHO,i,j,k) * S % U_z(IRHOU,i,j,k) - Q(IRHOU,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_z(IRHO,i,j,k)) &
                                - (1/Q(IRHO,i,j,k) * S % U_x(IRHOW,i,j,k) - Q(IRHOW,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_x(IRHO,i,j,k))
               end do            ; end do             ; end do
               
            case("omegaz")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = (1/Q(IRHO,i,j,k) * S % U_x(IRHOV,i,j,k) - Q(IRHOV,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_x(IRHO,i,j,k)) &
                                - (1/Q(IRHO,i,j,k) * S % U_y(IRHOU,i,j,k) - Q(IRHOU,i,j,k)/(Q(IRHO,i,j,k)**2) * S % U_y(IRHO,i,j,k))
               end do            ; end do             ; end do
   
            case("u")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
               end do            ; end do             ; end do
   
            case("v")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
               end do            ; end do             ; end do
   
            case("w")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
               end do            ; end do             ; end do
   
            case("mach")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  var(i,j,k) = POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))/POW2(Q(IRHO,i,j,k))     ! Vabs**2
                  var(i,j,k) = sqrt( var(i,j,k) / ( thermodynamics % gamma*(thermodynamics % gamma-1.0_RP)*(Q(IRHOE,i,j,k)/Q(IRHO,i,j,k)-0.5_RP * var(i,j,k)) ) )
               end do         ; end do         ; end do
      
            case("k")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  var(i,j,k) = 0.5_RP * (POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/Q(IRHO,i,j,k)
               end do         ; end do         ; end do
               
            case("q1")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHO,i,j,k)
               end do            ; end do             ; end do
               
            case("q2")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOU,i,j,k)
               end do            ; end do             ; end do
               
            case("q3")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOV,i,j,k)
               end do            ; end do             ; end do
               
            case("q4")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOW,i,j,k)
               end do            ; end do             ; end do
               
            case("q5")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IRHOE,i,j,k)
               end do            ; end do             ; end do
               
            end select
#endif
#ifdef MULTIPHASE
            select case (trim(self % variable(m)))
            case("static-pressure")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  var(i,j,k) = Q(IMP,i,j,k) + Q(IMC,i,j,k)*e % storage % mu(1,i,j,k) - 12.0_RP*multiphase%sigma*multiphase%invEps*(POW2(Q(IMC,i,j,k)*(1.0_RP-Q(IMC,i,j,k)))) &
                               - 0.25_RP*3.0_RP*multiphase % sigma * multiphase % eps * (POW2(e % storage % c_x(1,i,j,k))+POW2(e % storage % c_y(1,i,j,k))+POW2(e % storage % c_z(1,i,j,k)))
               end do         ; end do         ; end do
			   
            case("concentration")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
               end do            ; end do             ; end do
			   
            case("density")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  var(i,j,k) = c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2)
               end do            ; end do             ; end do
			   
            case("pressure")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(5,i,j,k) 
               end do            ; end do             ; end do
   
            case("velocity")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = sqrt(POW2(Q(2,i,j,k)) + POW2(Q(3,i,j,k)) + POW2(Q(4,i,j,k)))/sqrtRho
               end do         ; end do         ; end do
            
            case("omegax")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = (1/sqrtRho * S % U_y(IMSQRHOW,i,j,k) - Q(IMSQRHOW,i,j,k)/(sqrtRho**2) * S % U_y(IMC,i,j,k)) &
                                - (1/sqrtRho * S % U_z(IMSQRHOV,i,j,k) - Q(IMSQRHOV,i,j,k)/(sqrtRho**2) * S % U_z(IMC,i,j,k))
               end do            ; end do             ; end do
               
            case("omegay")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = (1/sqrtRho * S % U_z(IMSQRHOU,i,j,k) - Q(IMSQRHOU,i,j,k)/(sqrtRho**2) * S % U_z(IMC,i,j,k)) &
                                - (1/sqrtRho * S % U_x(IMSQRHOW,i,j,k) - Q(IMSQRHOW,i,j,k)/(sqrtRho**2) * S % U_x(IMC,i,j,k))
               end do            ; end do             ; end do
               
            case("omegaz")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = (1/sqrtRho * S % U_x(IMSQRHOV,i,j,k) - Q(IMSQRHOV,i,j,k)/(sqrtRho**2) * S % U_x(IMC,i,j,k)) &
                                - (1/sqrtRho * S % U_y(IMSQRHOU,i,j,k) - Q(IMSQRHOU,i,j,k)/(sqrtRho**2) * S % U_y(IMC,i,j,k))
               end do            ; end do             ; end do
   
            case("u")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = Q(IMSQRHOU,i,j,k) / sqrtRho
               end do            ; end do             ; end do
   
            case("v")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = Q(IMSQRHOV,i,j,k) / sqrtRho
               end do            ; end do             ; end do
   
            case("w")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
			      c = min(max(Q(1,i,j,k),0.0_RP),1.0_RP)
                  sqrtRho = sqrt(c * dimensionless % rho (1) + (1.0_RP-c) * dimensionless % rho(2))
                  var(i,j,k) = Q(IMSQRHOW,i,j,k) / sqrtRho
               end do            ; end do             ; end do
   
            case("q1")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IMC,i,j,k)
               end do            ; end do             ; end do
               
            case("q2")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IMSQRHOU,i,j,k)
               end do            ; end do             ; end do
               
            case("q3")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IMSQRHOV,i,j,k)
               end do            ; end do             ; end do
               
            case("q4")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IMSQRHOW,i,j,k)
               end do            ; end do             ; end do
               
            case("q5")
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
                  var(i,j,k) = Q(IMP,i,j,k)
               end do            ; end do             ; end do
			   
            end select
#endif 
   
                    value = 0.0_RP
                    do k = 0, e % Nxyz(3)    ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
                        value = value + var(i,j,k) * self % lxi(i,l) * self % leta(j,l) * self % lzeta(k,l)
                    end do               ; end do             ; end do
   
                    self % values(l+1, bufferPosition, m) = value
   
                    end associate
					deallocate(var)
                    end associate
#ifdef _HAS_MPI_            
                    if ( MPI_Process % doMPIAction ) then

!                   Share the result with the rest of the processes
!                   -----------------------------------------------         
                        call mpi_bcast(value, 1, MPI_DOUBLE, self % rank (l), MPI_COMM_WORLD, ierr)
                    end if
#endif
                    else

!                   Receive the result from the rank that contains the Plane
!                   --------------------------------------------------------
#ifdef _HAS_MPI_
                        if ( MPI_Process % doMPIAction ) then
                            call mpi_bcast(self % values(l+1, bufferPosition, m), 1, MPI_DOUBLE, self % rank (l), MPI_COMM_WORLD, ierr)
                        end if
#endif
                    end if
                END IF  
              END DO
            END DO
         end if
         self % intervalCount = self % intervalCount + 1
        
         
         if (self % intervalCount .EQ. self % interval) then
            self % intervalCount = 0 
         end if 

      end subroutine Plane_Update

      subroutine Plane_WriteToFile ( self, no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(PlaneSampling_t)  :: self
         integer                 :: no_of_lines
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
      
      end subroutine Plane_WriteToFile 

      subroutine Plane_LookInOtherPartitions(self, i)
         use MPI_Process_Info
         implicit none
         class(PlaneSampling_t)     :: self
         integer ,intent(in)        :: i
         integer                    :: allActives(MPI_Process % nProcs)
         integer                    :: active, ierr

         if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
!
!           Cast the logicals onto integers
!           -------------------------------
            if ( self % active (i) ) then
               active = 1
            else
               active = 0
            end if
            allActives=0
!
!           Gather all data from all processes
!           ----------------------------------
            call mpi_allgather(active, 1, MPI_INT, allActives, 1, MPI_INT, MPI_COMM_WORLD, ierr)
!
!           Check if any of them found the Plane
!           ------------------------------------
            if ( any(allActives .eq. 1) ) then
!
!              Assign the domain of the partition that contains the Plane
!              ----------------------------------------------------------
               self % active (i) = .true.
               self % rank (i) = maxloc(allActives, dim = 1) - 1

            else
!
!              Disable the Plane
!              -----------------
               self % active (i) = .false.
               self % rank (i)  = -1
               self % eID (i) = 1

            end if
#endif
         else
!
!           Without MPI select the rank 0 as default
!           ----------------------------------------
            self % rank (i)= 0

         end if

      end subroutine Plane_LookInOtherPartitions
      
      elemental subroutine Plane_Destruct (self)
         implicit none
         class(PlaneSampling_t), intent(inout) :: self
         
         safedeallocate (self % values)
         safedeallocate (self % lxi)
         safedeallocate (self % leta)
         safedeallocate (self % lzeta)
         safedeallocate (self % active)
         safedeallocate (self % rank)
         safedeallocate (self % eID)
         safedeallocate (self % x)
         safedeallocate (self % fileName)
         safedeallocate (self % variable)
      end subroutine Plane_Destruct
      
      elemental subroutine Plane_Assign (to, from)
         implicit none
         class(PlaneSampling_t), intent(inout) :: to
         type(PlaneSampling_t) , intent(in) :: from
         
         safedeallocate ( to % active )
         allocate ( to % active ( size(from % active) ) )
         to % active = from % active
         
         safedeallocate ( to % rank )
         allocate ( to % rank ( size(from % rank) ) )
         to % rank = from % rank
         
         to % ID = from %  ID
         to % nVariables      = from % nVariables
         to % interval        = from % interval
         to % bufferSize      = from % bufferSize
         to % bufferLine      = from % bufferLine
         to % intervalCount   = from % intervalCount
         to % N               = from % N
         to % nNodes          = from % nNodes
         to % disturbanceData = from % disturbanceData
         
         safedeallocate ( to % eID )
         allocate ( to % eID ( size(from % eID) ) )
         to % eID = from % eID      

         safedeallocate ( to % lxi )
         allocate ( to % lxi ( size(from % lxi,1), size(from % lxi,2) ) )
         to % lxi = from % lxi
         
         safedeallocate ( to % leta )
         allocate ( to % leta ( size(from % leta,1), size(from % leta,2) ) )
         to % leta = from % leta
         
         safedeallocate ( to % lzeta )
         allocate ( to % lzeta ( size(from % lzeta,1), size(from % lzeta,2) ) )
         to % lzeta = from % lzeta
         
         safedeallocate ( to % values )
         allocate ( to % values( size(from % values,1),size(from % values,2),size(from % values,3) ) )
         to % values = from % values
         
         safedeallocate ( to % x )
         allocate ( to % x( size(from % x,1),size(from % x,2) ) )
         to % x = from % x       
         
         safedeallocate ( to % xi )
         allocate ( to % xi( size(from % xi,1),size(from % xi,2) ) )
         to % xi = from % xi    

         safedeallocate ( to % fileName )
         allocate ( to % fileName ( size(from % fileName) ) )
         to % fileName = from % fileName
         
         to % planeName = from % planeName

         safedeallocate ( to % variable )
         allocate ( to % variable ( size(from % variable) ) )
         to % variable = from % variable         

         
         
      end subroutine Plane_Assign
      
end module PlaneSampling

