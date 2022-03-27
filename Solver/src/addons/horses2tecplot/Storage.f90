!
!//////////////////////////////////////////////////////
!
!   @File:    Storage.f90
!   @Author:  UNK
!   @Created: UNK
!   @Last revision date: Mon Sep  6 22:44:49 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: 3334a040b8cdf3201850a2deec9950c84f2dc21f
!
! TODO: Create destructors
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module Storage
   use SMConstants
   use PhysicsStorage
   use SolutionFile
   use VariableConversion
   implicit none
   
   private
   public Mesh_t, Element_t, Boundary_t
   public NVARS, NGRADVARS, hasMPIranks, hasBoundaries

   integer     :: NVARS, NGRADVARS
   logical     :: hasMPIranks, hasBoundaries

   type Element_t
!                                /* Mesh quantities */
      integer                    :: eID
      integer                    :: Nmesh(NDIM)
      integer                    :: mpi_rank = 0
      real(kind=RP), pointer     :: x(:,:,:,:)
!                                /* Solution quantities */
      integer                    :: Nsol(NDIM)
      real(kind=RP), pointer     :: Q(:,:,:,:)
      real(kind=RP), pointer     :: QDot(:,:,:,:)
      real(kind=RP), pointer     :: U_x(:,:,:,:)
      real(kind=RP), pointer     :: U_y(:,:,:,:)
      real(kind=RP), pointer     :: U_z(:,:,:,:)
      real(kind=RP), pointer     :: Q_x(:,:,:,:)
      real(kind=RP), pointer     :: Q_y(:,:,:,:)
      real(kind=RP), pointer     :: Q_z(:,:,:,:)
      real(kind=RP), pointer     :: mu_NS(:,:,:,:)
      real(kind=RP), pointer     :: stats(:,:,:,:)
!                                /* Output quantities */
      integer                    :: Nout(NDIM)
      real(kind=RP), pointer     :: xOut(:,:,:,:)
      real(kind=RP), pointer     :: Qout(:,:,:,:)
      real(kind=RP), pointer     :: QDot_out(:,:,:,:)
      real(kind=RP), pointer     :: U_xout(:,:,:,:)
      real(kind=RP), pointer     :: U_yout(:,:,:,:)
      real(kind=RP), pointer     :: U_zout(:,:,:,:)
      real(kind=RP), pointer     :: mu_NSout(:,:,:,:)
      real(kind=RP), pointer     :: statsout(:,:,:,:)
      
      real(kind=RP), allocatable :: outputVars(:,:,:,:)
   end type Element_t
   
   type Boundary_t
      character(len=LINE_LENGTH) :: Name
      integer                    :: no_of_faces
      integer, allocatable       :: elements(:)
      integer, allocatable       :: elementSides(:)
   end type Boundary_t
   
   type Mesh_t
      integer  :: no_of_elements
      integer  :: nodeType
      type(Element_t),   allocatable    :: elements(:)
      type(Boundary_t),  allocatable    :: boundaries(:)
      character(len=LINE_LENGTH) :: meshName
      character(len=LINE_LENGTH) :: solutionName
      real(kind=RP)              :: refs(NO_OF_SAVED_REFS)
      logical                    :: hasGradients
      logical                    :: isStatistics
      contains
         procedure   :: ReadMesh     => Mesh_ReadMesh
         procedure   :: ReadSolution => Mesh_ReadSolution
   end type Mesh_t

   contains
      subroutine Mesh_ReadMesh(self,meshName)
         use Headers
         implicit none
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: meshName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, dimension(:), allocatable     :: arrayDimensions
         integer                                :: fileType, dimensionsSize
         integer                                :: fid, eID, i, pos
         character(len=1024)                    :: msg
         character(len=LINE_LENGTH)             :: flag

         self % meshName = trim(meshName)
!
!        Get mesh node type and type of file
!        -----------------------------------
         self % nodeType = getSolutionFileNodeType(meshName)
         fileType = getSolutionFileType(meshName)
!
!        Get number of elements
!        ----------------------
         self % no_of_elements = getSolutionFileNoOfElements(meshName)
!
!        Allocate elements
!        -----------------
         allocate(self % elements(self % no_of_elements))
!
!        Allocate dimension based on type of file
!        ----------------------------------------
         select case (fileType)
         case (ZONE_MESH_FILE)
             dimensionsSize = 3
         case default
             dimensionsSize = 4
         end select

         allocate(arrayDimensions(dimensionsSize))
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(meshName)
      
         do eID = 1, self % no_of_elements
            associate ( e => self % elements(eID) ) 
            e % eID = eID
            
            call getSolutionFileArrayDimensions(fid,arrayDimensions)
!
!           Allocate memory for the coordinates
!           -----------------------------------            
            ! e % Nmesh(1:3) = arrayDimensions(2:4) - 1 
            e % Nmesh(1:dimensionsSize-1) = arrayDimensions(2:dimensionsSize) - 1 
!           Use a 0 index for the case of a surface mesh
            if (dimensionsSize .eq. 3) e % Nmesh(3) = 0
            allocate( e % x(NDIM,0:e % Nmesh(1),0:e % Nmesh(2),0:e % Nmesh(3)) )
!
!           Read data
!           ---------
            read(fid) e % x

            end associate
         end do
!
!        Close file
!        ----------
         close(fid)
         
!        Read boundary file (if present)
!        -------------------------------
         do i = 1, command_argument_count()
            call get_command_argument(i, flag)
            
            pos = index(trim(flag),"--boundary-file")
            if ( pos .ne. 0 ) then
               hasBoundaries = .TRUE.
               
               call readBoundaryFile(self % boundaries, flag)
               
            end if
         end do
         
!
!        Describe the mesh
!        -----------------
         write(msg,'(A,A,A)') 'Mesh file "',trim(meshName),'":'
         write(STD_OUT,'(/)')
         call SubSection_Header(trim(msg))

         write(STD_OUT,'(30X,A,A30,I0)') "->", "Number of elements: ", self % no_of_elements
         select case ( self % nodeType )
         case(1)
            write(STD_OUT,'(30X,A,A30,A)') "->","Discretization nodes: ","Gauss"
         case(2)
            write(STD_OUT,'(30X,A,A30,A)') "->","Discretization nodes: ","Gauss-Lobatto"
         end select

!
!        Describe the boundaries
!        -----------------
         if (hasBoundaries) then
             pos = index(trim(flag),"=")
             flag = flag(pos+1:len_trim(flag))
             write(msg,'(A,A,A)') 'Boundary file "',trim(flag),'":'
             write(STD_OUT,'(/)')
             call SubSection_Header(trim(msg))
             write(STD_OUT,'(30X,A,A30,I0)') "->","Number of Boundaries: ", size(self % boundaries)
             
         end if 

      end subroutine Mesh_ReadMesh

      subroutine Mesh_ReadSolution(self,solutionName)
         use Headers
         implicit none
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: solutionName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: no_of_elements
         integer       :: arrayDimensions(4)
         integer       :: fid, eID, pos
         integer       :: i,j,k
         integer       :: iter
         real(kind=RP) :: time
         character(len=1024)  :: msg
         character(len=LINE_LENGTH)    :: flag

         self % solutionName = trim(solutionName)
!
!        Get the solution file type
!        --------------------------
         select case ( getSolutionFileType(trim(solutionName)) )

         case (SOLUTION_FILE)
            self % hasGradients = .false.
            self % isStatistics = .false.

         case (SOLUTION_AND_GRADIENTS_FILE)
            self % hasGradients = .true.
            self % isStatistics = .false.

         case (STATS_FILE)
            self % hasGradients = .false.
            self % isStatistics = .true.

         case default
            print*, "File expected to be a solution file"
            errorMessage(STD_OUT)
            stop
         end select
!
!        Get node type
!        -------------
         if ( getSolutionFileNodeType(solutionName) .ne. self % nodeType ) then
            print*, "Solution and Mesh node type differs"
            errorMessage(STD_OUT)
            stop
         end if
!
!        Get number of elements
!        ----------------------
         no_of_elements = getSolutionFileNoOfElements(solutionName)
         if ( self % no_of_elements .ne. no_of_elements ) then
            write(STD_OUT,'(30X,A,I0,A,I0,A)') "The number of elements in the mesh (",self % no_of_elements,&
                                           ") differs to that of the solution (",no_of_elements,")."
            errorMessage(STD_OUT)
            stop 
         end if
!
!        Get time and iteration
!        ----------------------
         call getSolutionFileTimeAndIteration(trim(solutionName),iter,time)
!
!        Read reference values
!        ---------------------
         self % refs = getSolutionFileReferenceValues(trim(solutionName))
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(solutionName)
         
         ! call set_getVelocityGradients(GRADVARS_STATE) ! FIXME: MIGHT BE NEEDED FOR HORSES2PLT
         ! write(STD_OUT,'(15X,A)') " WARNING horses2tecplot.90 :: Velocity Gradients set to default (GRADVARS_STATE)"
      
         if ( .not. self % isStatistics ) then
            do eID = 1, self % no_of_elements
               associate ( e => self % elements(eID) ) 
               call getSolutionFileArrayDimensions(fid,arrayDimensions)

               NVARS = arrayDimensions(1)
               NGRADVARS = NVARS     ! TODO: Read NCONS and NGRAD from physics!
               select case(NVARS)
               case(1)      ; NGRADVARS = 1
               case(6)      ; NGRADVARS = 3
               case default ; NGRADVARS = 3
               end select
!   
!              Allocate memory for the coordinates
!              -----------------------------------            
               e % Nsol(1:3) = arrayDimensions(2:4) - 1
               allocate( e % Q(1:arrayDimensions(1),0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
!   
!              Read data
!              ---------
               read(fid) e % Q
   
               if ( self % hasGradients ) then
!   
!                 Allocate memory for the gradients
!                 ---------------------------------
                  allocate( e % Q_x(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) ) ! TODO: Allocate depending on physics
                  allocate( e % Q_y(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  allocate( e % Q_z(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  
                  allocate( e % U_x(1:NGRADVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  allocate( e % U_y(1:NGRADVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  allocate( e % U_z(1:NGRADVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
!   
!                 Read data
!                 ---------
                  read(fid) e % Q_x
                  read(fid) e % Q_y
                  read(fid) e % Q_z

!                 Call set_getVelocityGradients to make the pointer to the actual subroutine, is needed only for the NS
!                 Set state as is the default option TODO point to the correct one if its posible (oscar note)
!                 ---------------------------
                  call set_getVelocityGradients(GRADVARS_STATE)
                  
                  ! Following block works for NS, CH, NSCH and iNS .... but not iNSCH: change 5 by 6 to use iNSCH (NS won't work) 
                  if (NVARS .ge. 5) then 
                     do k = 0,e % Nsol(3) ; do j = 0, e % Nsol(2) ; do i = 0, e % Nsol(1)
                        call getVelocityGradients(e % Q(:,i,j,k), e % Q_x(1:5,i,j,k), e % Q_y(1:5,i,j,k), e % Q_z(1:5,i,j,k), &
                                                  e % U_x(1:3,i,j,k), e % U_y(1:3,i,j,k), e % U_z(1:3,i,j,k) )
                     end do               ; end do                ; end do
                     !if (NVARS == 6) then
                     !   e % U_x(6,:,:,:) = e % Q_x(6,:,:,:)
                     !   e % U_y(6,:,:,:) = e % Q_y(6,:,:,:)
                     !   e % U_z(6,:,:,:) = e % Q_z(6,:,:,:)
                     !end if
                  else
                     e % U_x = e % Q_x(1:NGRADVARS,:,:,:)
                     e % U_y = e % Q_y(1:NGRADVARS,:,:,:)
                     e % U_z = e % Q_z(1:NGRADVARS,:,:,:)
                  end if
                  
               end if
               end associate
            end do

         else

            do eID = 1, self % no_of_elements
               associate ( e => self % elements(eID) )
               call getSolutionFileArrayDimensions(fid,arrayDimensions)
!
!              Allocate memory for the statistics
!              ----------------------------------
               e % Nsol(1:3) = arrayDimensions(2:4) - 1
               allocate( e % stats(1:9,0:e % Nsol(1), 0:e % Nsol(2), 0:e % Nsol(3)) )
!
!              Read data
!              ---------
               read(fid) e % stats
               end associate
            end do
         end if
!
!        Close file
!        ----------
         close(fid)
         
!        Read mpi ranks (if present)
!        ---------------------------
         do i = 1, command_argument_count()
            call get_command_argument(i, flag)
            
            pos = index(trim(flag),"--partition-file")
            if ( pos .ne. 0 ) then
               hasMPIranks = .TRUE.
               
               call readPartitionFile(self,flag)
               
            end if
         end do
         
!
!        Describe the solution
!        ---------------------
         write(msg,'(A,A,A)') 'Solution file "',trim(solutionName),'":'
         write(STD_OUT,'(/)')
         call SubSection_Header(trim(msg))

         if ( self % isStatistics ) then
            write(STD_OUT,'(30X,A,A30)') "->","File is statistics file."
         
         else
            if ( self % hasGradients ) then
               write(STD_OUT,'(30X,A,A40,A)') "->","Solution file contains gradients: ", "yes"
            else
               write(STD_OUT,'(30X,A,A40,A)') "->","Solution file contains gradients: ", "no"
            end if
         end if

         write(STD_OUT,'(30X,A,A30,I0)') "->","Iteration: ", iter
         write(STD_OUT,'(30X,A,A30,ES10.3)') "->","Time: ", time
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference velocity: ", self % refs(V_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference density: ", self % refs(RHO_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference Temperature: ", self % refs(T_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference Mach number: ", self % refs(MACH_REF)
      end subroutine Mesh_ReadSolution
      
      
      subroutine readPartitionFile(self,flag)
         implicit none
         !-arguments-----------------------------------------------
         class(Mesh_t)   , intent(inout)  :: self
         character(len=*), intent(in)     :: flag
         !-local-variables-----------------------------------------
         integer                    :: pos, nelem, eID, fID
         character(len=LINE_LENGTH) :: partitionFileName
         !---------------------------------------------------------
         
         pos = index(trim(flag),"=")
               
         if ( pos .eq. 0 ) then
            print*, 'Missing "=" operator in --partition-file flag'
            errorMessage(STD_OUT)
            stop
         end if
         
         partitionFileName = flag(pos+1:len_trim(flag))
         
         open(newunit = fID, file=trim(partitionFileName),action='read')
         
            read(fID,*) nelem
            
            do eID = 1, nelem
               read(fID,*) self % elements(eID) % mpi_rank
            end do
            
         close(fID)
      end subroutine readPartitionFile
      
      subroutine readBoundaryFile(boundaries, flag)
         implicit none
         !-arguments-----------------------------------------------
         type(Boundary_t), allocatable, intent(inout)  :: boundaries(:)
         character(len=*)             , intent(in)     :: flag
         !-local-variables-----------------------------------------
         integer                    :: pos, fd, no_of_boundaries,bID
         character(len=LINE_LENGTH) :: boundaryFileName
         !---------------------------------------------------------
         
         pos = index(trim(flag),"=")
               
         if ( pos .eq. 0 ) then
            print*, 'Missing "=" operator in --boundary-file flag'
            errorMessage(STD_OUT)
            stop
         end if
         
         boundaryFileName = flag(pos+1:len_trim(flag))
         
         open(newunit = fd, file=trim(boundaryFileName),action='read')
         
            read(fd,*) no_of_boundaries
            allocate ( boundaries(no_of_boundaries) )
            
            do bID = 1, no_of_boundaries
               
               read(fd,*) boundaries(bID) % Name
               read(fd,*) boundaries(bID) % no_of_faces
               
               allocate ( boundaries(bID) % elements    (boundaries(bID) % no_of_faces) )
               allocate ( boundaries(bID) % elementSides(boundaries(bID) % no_of_faces) )
               
               read(fd,*) boundaries(bID) % elements
               read(fd,*) boundaries(bID) % elementSides
            end do
            
         close(fd)
         
      end subroutine readBoundaryFile
end module Storage
