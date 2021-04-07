!
!   @File:    ObserverClass.f90
!   @Author:  Oscar Marino (oscar.marino@upm.es)
!   @Created: Mar 22 2020
!   @Last revision date: 
!   @Last revision author: 
!   @Last revision commit: 
!
!//////////////////////////////////////////////////////
!
!This class represents the observer for the fW-H accoustic analogy, including the relations with several observers

#include "Includes.h"
Module  FWHObseverClass  !

   use SMConstants
   use FaceClass
   use Physics
   use PhysicsStorage
   use NodalStorageClass
   use MonitorDefinitions, only: BUFFER_SIZE_DEFAULT, BUFFER_SIZE, STR_LEN_MONITORS
   use ZoneClass
   use HexMeshClass
   use MPI_Process_Info
   Implicit None

!
!  *****************************
!  Observer source pair class definition
!  class for the coupling of each pair of observer and source(face)
!  mainly accoustic geometrical relations and face link
!  *****************************
   type ObserverSourcePairClass
       real(kind=RP), dimension(:,:,:), allocatable        :: rVect
       real(kind=RP), dimension(:,:),   allocatable        :: r
       real(kind=RP), dimension(:,:),   allocatable        :: re
       real(kind=RP), dimension(:,:,:), allocatable        :: reUnitVect 
       real(kind=RP), dimension(:,:),   allocatable        :: reStar
       real(kind=RP), dimension(:,:,:), allocatable        :: reStarUnitVect 
       real(kind=RP)                                       :: tDelay
       integer                                             :: faceIDinMesh    ! ID of the source (face) at the Mesh array (linked list)

       contains

           procedure :: construct       => ObserverSourcePairConstruct
           procedure :: destruct        => ObserverSourcePairDestruct
           procedure :: FWHSurfaceIntegral

   end type ObserverSourcePairClass
!
!  *****************************
!  General observer class definition
!   (similar to a monitor, mostly surface monitor in many behaviours)
!  *****************************
!  
   type ObserverClass

       integer                                                         :: ID
       real(kind=RP), dimension(NDIM)                                  :: x        ! position of the observer at global coordinates
       integer                                                         :: numberOfFaces
       class(ObserverSourcePairClass), dimension(:), allocatable       :: sourcePair
       real(kind=RP), dimension(:,:), allocatable                      :: Pac      ! accoustic pressure, two componenets and the total (sum)
       real(kind=RP)                                                   :: tDelay
       logical                                                         :: active
       character(len=STR_LEN_MONITORS)                                 :: observerName
       character(len=STR_LEN_MONITORS)                                 :: fileName

       contains

           procedure :: construct      => ObserverConstruct
           procedure :: destruct       => ObserverDestruct
           procedure :: update         => ObserverUpdate
           procedure :: writeToFile    => ObserverWriteToFile
           procedure :: updateTdelay   => ObserverUpdateTdelay

   end type ObserverClass

   contains

!/////////////////////////////////////////////////////////////////////////
!           OBSERVER CLASS PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

   Subroutine ObserverConstruct(self, sourceZone, mesh, ID, solution_file, FirstCall)

!        *****************************************************************************
!              This subroutine initializes the observer similar to a monitor. The following
!           data is obtained from the case file:
!              -> Name: The observer name (10 characters maximum)
!              -> x: The observer position
!        *****************************************************************************

       use ParamfileRegions
       use FileReadingUtilities, only: getRealArrayFromString
       implicit none

       class(ObserverClass)                                 :: self
       class(Zone_t), intent(in)                            :: sourceZone
       class(HexMesh), intent(in)                           :: mesh
       integer, intent(in)                                  :: ID
       character(len=*), intent(in)                         :: solution_file
       logical, intent(in)                                  :: FirstCall

       ! local variables
       character(len=STR_LEN_MONITORS)  :: in_label
       character(len=STR_LEN_MONITORS)  :: fileName
       character(len=STR_LEN_MONITORS)  :: paramFile
       character(len=STR_LEN_MONITORS)  :: coordinates
       integer                          :: fID
       integer                          :: MeshFaceID, zoneFaceID
!
!      Get observer ID
!      --------------
       self % ID = ID
!
!      Search for the parameters in the case file
!      ------------------------------------------
       write(in_label , '(A,I0)') "#define accoustic observer " , self % ID

       call get_command_argument(1, paramFile)
       call readValueInRegion(trim ( paramFile), "name",   self % observerName, in_label, "# end" ) 
       call readValueInRegion(trim(paramFile), "position", coordinates        , in_label, "# end" )

!      Get the coordinates
!      -------------------
      ! print *, "Observer: ", trim(self%observerName)
       self % x = getRealArrayFromString(coordinates)

!     Enable the observer
!     ------------------
      self % active = .true.
      allocate ( self % Pac(BUFFER_SIZE,3) )

      !     ------------------
!     Get source information
      self % numberOfFaces = sourceZone % no_of_faces

!     Construct each pair observer-source
!     ------------------
      allocate( self % sourcePair(self % numberOfFaces) )
!     Loop the zone to get faces
      do zoneFaceID = 1, self % numberOfFaces
!         Face global ID
          MeshFaceID = sourceZone % faces(zoneFaceID)
          call self % sourcePair(zoneFaceID) % construct(self % x, mesh % faces(MeshFaceID), MeshFaceID, FirstCall)
      end do  

!     Set the average time delay of the observer
!     -------------------------------------------------
      call self % updateTdelay()

!     Prepare the file in which the observer is exported
!     -------------------------------------------------
      write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % observerName) , ".observer"
!
!     Create file
!     -----------
      if (FirstCall) then
         open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 

!        Write the file headers
!        ----------------------
         write( fID , '(A20,A  )') "Observer name:      ", trim(self % observerName)
         write( fID , '(A20,ES24.10)') "x coordinate: ", self % x(1)
         write( fID , '(A20,ES24.10)') "y coordinate: ", self % x(2)
         write( fID , '(A20,ES24.10)') "z coordinate: ", self % x(3)

         write( fID , * )
         write( fID , '(A10,5(2X,A24))' ) "Iteration" , "Time" , "Observer Time", "P'T", "P'L", "P'"

         close ( fID )
      end if

   End Subroutine ObserverConstruct

   Subroutine ObserverUpdate(self, mesh, BufferPosition, isSolid)

!     *******************************************************************
!        This subroutine updates the observer accoustic pressure computing it from
!        the mesh storage. It is stored in the "bufferPosition" position of the 
!        buffer.
!         TODO: use openmp (commented)
!         TODO: use mpi (see surface integral)
!     *******************************************************************
!
      implicit none
      class (ObserverClass)                                :: self
      class (HexMesh), intent(in)                          :: mesh
      integer,intent(in)                                   :: bufferPosition
      logical, intent(in)                                  :: isSolid

      ! local variables
      real(kind=RP)                                        :: Pt, Pl  ! pressure of each pair
      real(kind=RP), dimension(3)                          :: localPacc, Pacc   ! temporal variable to store the sum of the pressure
      real(kind=RP)                                        :: valx, valy, valz
      integer                                              :: zoneFaceID, meshFaceID !,ierr

!     Initialization
!     --------------            
      self % Pac(bufferPosition,:) = 0.0_RP
      Pacc = 0.0_RP
      valx = 0.0_RP
      valy = 0.0_RP
      valz = 0.0_RP


!        Loop the pairs (equivalent to loop the zone) and get the values
!        ---------------------------------------
!!$omp do private(fID,localVal) reduction(+:valx,valy,valz) schedule(runtime)
      ! print *, "obs: ", trim(self%observerName)
      do zoneFaceID = 1, self % numberOfFaces
!        Compute the integral
!        --------------------
         meshFaceID = self % sourcePair(zoneFaceID) % faceIDinMesh
         localPacc = self % sourcePair(zoneFaceID) % FWHSurfaceIntegral( mesh % faces(meshFaceID), isSolid )

         ! sum without interpolate: supose little change of each tDelay
         valx = valx + localPacc(1)
         valy = valy + localPacc(2)
         valz = valz + localPacc(3)
      end do  
!!$omp end do
!!$omp end parallel

      Pacc = (/valx, valy, valz/)

! #ifdef _HAS_MPI_
!       localPacc = Pacc
!       call mpi_allreduce(localPacc, Pacc, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
! #endif

      self % Pac(bufferPosition,:) = Pacc
      ! print *, "obs: ", trim(self%observerName), "id: ", self%id
      ! print *, "Pacc obs: ", self%Pac(bufferPosition,:)

   End Subroutine ObserverUpdate

   Subroutine ObserverUpdateTdelay(self)

!     *******************************************************************
!        This subroutine updates the observer time delay. For static surfaces it 
!        doesn't need to be updated at every iteration
!     *******************************************************************
!
      implicit none
      class   (ObserverClass)                              :: self
      
      ! local variables
      integer                                              :: i
      real(kind=RP)                                        :: t

      t = 0.0_RP
      do i = 1, self % numberOfFaces
          ! for moving surfaces, each pair of observer-source should updated the tDelay first
          ! print *, "time of source ", i, "of observer ", self % ID, "is: ", self % sourcePair(i)%tDelay
          t = t + self % sourcePair(i) % tDelay
      end do  

      self % tDelay = t / real(self % numberOfFaces,RP)
          ! print *, "time of average of observer ", self % ID, "is: ", self % tDelay

   End Subroutine ObserverUpdateTdelay

   Subroutine ObserverWriteToFile(self, iter, tsource, no_of_lines)
!
!     *************************************************************
!           This subroutine writes the buffer to the file.
!     *************************************************************
!
      implicit none  
      class(ObserverClass)                                 :: self
      integer, dimension(:)                                :: iter
      real(kind=RP), dimension(:)                          :: tsource
      integer                                              :: no_of_lines
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: i
      integer                    :: fID

      if ( MPI_Process % isRoot ) then

         open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
      
         do i = 1 , no_of_lines
            write( fID , '(I10,5(2X,ES24.16))' ) iter(i) , tsource(i), tsource(i) + self % tDelay, self % Pac(i,:)
         end do
      
         close ( fID )
      end if
      
      if ( no_of_lines .ne. 0 ) self % Pac(1,:) = self % Pac(no_of_lines,:)
      
   End Subroutine ObserverWriteToFile

   Subroutine ObserverDestruct(self)

        implicit none
        class(ObserverClass), intent(inout)               :: self

        safedeallocate (self % Pac)
        call self % sourcePair % destruct
        safedeallocate (self % sourcePair)

   End Subroutine ObserverDestruct

!/////////////////////////////////////////////////////////////////////////
!           OBSERVER SOURCE PAIR CLASS PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

   Subroutine  ObserverSourcePairConstruct(self, x, f, fID, FirstCall)

       ! use fluiddata
       use FWHDefinitions, only: rho0, P0, c0, U0, M0, fwGamma2
       implicit none

       class(ObserverSourcePairClass)                      :: self
       real(kind=RP), dimension(NDIM), intent(in)          :: x       ! observer position
       type(face), intent(in)                              :: f    ! source
       integer                                             :: fID
       logical, intent(in)                                 :: FirstCall

       ! local variables
       integer                                             :: Nx,Ny
       integer                                             :: i, j
       real(kind=RP)                                       :: fwGammaInv

       self % faceIDinMesh = fID

       Nx = f % Nf(1)
       Ny = f % Nf(2)

       allocate( self % r(0:Nx,0:Ny), self % re(0:Nx,0:Ny), self % reStar(0:Nx,0:Ny) )
       allocate( self % rVect(NDIM,0:Nx,0:Ny), self % reUnitVect(NDIM,0:Nx,0:Ny) ,self % reStarUnitVect(NDIM,0:Nx,0:Ny) )

       fwGammaInv = 1.0_RP / sqrt(fwGamma2)
       ! source position, for each node of the face
       associate (y => f % geom % x)
           do j= 0, Ny; do i = 0,Nx
               ! store geometrical accoustic relations for each node
               self % rVect(:,i,j) = x(:) - y(:,i,j)
               self % r(i,j) = norm2(self % rVect(:,i,j))
               self % reStar(i,j) = fwGammaInv*sqrt( self%r(i,j)**2 + fwGamma2*( dot_product(M0, self%rVect(:,i,j)) )**2 )
               self % reStarUnitVect(:,i,j) = ( self%rVect(:,i,j) + fwGamma2*dot_product(M0, self%rVect(:,i,j))*M0(:) ) / &
                                            (fwGamma2*self%reStar(i,j))
               self % re(i,j) = fwGamma2*( self%reStar(i,j) - dot_product(M0, self%rVect(:,i,j)) )
               self % reUnitVect(:,i,j) = fwGamma2*( self%reStarUnitVect(:,i,j) - M0(:) )
               self % tDelay = (sum(self%re))/real(size(self%re),RP) / c0
           end do; end do
       end associate

       ! print *, "r: ", self%r(1,1)
       ! print *, "rv: ", self%rVect(:,1,1)
       ! print *, "res: ", self%reStar(1,1)
       ! print *, "reStarUnitVect: ", self%reStarUnitVect(:,1,1), "norm: ",norm2(self%reStarUnitVect(:,1,1))
       ! print *, "re: ", self%re(1,1)
       ! print *, "reUnitVect: ", self%reUnitVect(:,1,1), "norm: ",norm2(self%reUnitVect(:,1,1))

   End Subroutine ObserverSourcePairConstruct 

  elemental Subroutine ObserverSourcePairDestruct(self)

      Class(ObserverSourcePairClass), intent(inout)       :: self
 
      safedeallocate(self % rVect)
      safedeallocate(self % r)
      safedeallocate(self % re)
      safedeallocate(self % reUnitVect)
      safedeallocate(self % reStar)
      safedeallocate(self % reStarUnitVect)

  End Subroutine ObserverSourcePairDestruct 

   ! calculate the surface integrals of the FW-H analogy for stacionary surfaces (permable or impermeable) with a general flow
   ! direction of the medium
   ! the integrals are for a single face (pane in FWH terminology) for a single observer
!         TODO: check if is more efficient to store FWHvariables for each face instead of calculating it always
!               for many observers, its being recomputed as many as observers

   Function FWHSurfaceIntegral(self, f, isSolid) result(Pacc)

       use FWHDefinitions, only: rho0, P0, c0, U0, M0
       implicit none

       class(ObserverSourcePairClass)                      :: self
       class(Face), intent(in)                             :: f
       real(kind=RP),dimension(3)                          :: Pacc  ! accoustic pressure values
       logical, intent(in)                                 :: isSolid

       ! local variables
       integer                                             :: i, j  ! face indexes
       real(kind=RP), dimension(NDIM)                      :: Qi,QiDot, n
       real(kind=RP), dimension(NDIM,NDIM)                 :: Lij, LijDot
       type(NodalStorage_t), pointer                       :: spAxi, spAeta
       real(kind=RP)                                       :: Pt, Pl

       ! Initialization
       Pt = 0.0_RP
       Pl = 0.0_RP
       spAxi  => NodalStorage(f % Nf(1))
       spAeta => NodalStorage(f % Nf(2))


       associate( Q => f % storage(1) % Q )
           associate( Qdot => f % storage(1) % Qdot )

    !           **********************************
    !           Computes the surface integral
    !              I = \int vec{f}·vec{n} * vec{g}·vec{r} dS
    !           **********************************
       ! print *, "reStarUnitVect: ", self%reStarUnitVect
    !
                do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)

                   n = f % geom % normal(:,i,j)
                   call calculateFWHVariables(Q(:,i,j), Qdot(:,i,j), isSolid, Qi, QiDot, Lij, LijDot)
                   ! print *, "Q: ", Q(:,i,j)
                   ! print *, "Qdot: ", Qdot(:,i,j)
                   ! print *, "Qi: ", Qi
                   ! print *, "QiDot: ", QiDot
                   ! print *, "Lij: ", Lij
                   ! print *, "LijDot: ", LijDot

                   ! loading term integrals
                   Pl = Pl +  dot_product(matmul(LijDot,n(:)),self%reUnitVect(:,i,j)) / (self%reStar(i,j) * c0) * &
                             spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
                   Pl = Pl +  dot_product(matmul(Lij,n(:)),self%reStarUnitVect(:,i,j)) / (self%reStar(i,j)**2) * &
                             spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)

                   ! thickness term integrals, only for permable surfaces
                   if (.not. isSolid) then
                       Pt = Pt + (1 - dot_product(M0(:),self%reUnitVect(:,i,j))) * dot_product(QiDot(:),n(:)) / (self%reStar(i,j)) * &
                                 spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
                       Pt = Pt -  dot_product(U0(:),self%reStarUnitVect(:,i,j)) * dot_product(Qi(:),n(:)) / (self%reStar(i,j)**2) * &
                                 spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
                   end if  
                end do          ;    end do
           end associate
       end associate

       Pt = Pt / (4.0_RP * PI)
       Pl = Pl / (4.0_RP * PI)
                   ! print *, "Pl: ", pl

      ! get total accoustic pressure as the sum of the two components (the quadrapol terms are being ignored)
       Pacc = (/Pt, Pl, Pt+Pl/)
       ! print *, "Pacc: ", Pacc

   End Function FWHSurfaceIntegral 

   Subroutine SourceProlongSolution(observer, mesh)

!     *******************************************************************
!        This subroutine prolong the solution from the mesh storage to the faces (source).
!         TODO: use openmp (commented)
!         TODO: use mpi (see surface integral)
!     *******************************************************************
!
      use ElementClass
      implicit none
      class (ObserverClass), intent(in)                    :: observer
      class (HexMesh), intent(inout), target               :: mesh

      ! local variables
      integer                                              :: zoneFaceID, meshFaceID, eID !,ierr
      integer, dimension(6)                                :: meshFaceIDs
      class(Element), pointer                              :: elements(:)

!     *************************
!     Perform the interpolation
!     *************************
!
      elements => mesh % elements
!!$omp parallel private(fID, eID, fIDs, localVal) shared(elements,mesh,NodalStorage,zoneID,integralType,val,&
!!$omp&                                        valx,valy,valz,computeGradients)
!!$omp single

!        Loop the zone to get faces and elements
!        ---------------------------------------
      do zoneFaceID = 1, observer % numberOfFaces
         meshFaceID = observer % sourcePair(zoneFaceID) % faceIDinMesh

         eID = mesh % faces(meshFaceID) % elementIDs(1)
         meshFaceIDs = mesh % elements(eID) % faceIDs

!!$omp task depend(inout:elements(eID))
         call elements(eID) % ProlongSolutionToFaces(NCONS,&
                                                            mesh % faces(meshFaceIDs(1)),&
                                                            mesh % faces(meshFaceIDs(2)),&
                                                            mesh % faces(meshFaceIDs(3)),&
                                                            mesh % faces(meshFaceIDs(4)),&
                                                            mesh % faces(meshFaceIDs(5)),&
                                                            mesh % faces(meshFaceIDs(6)),&
                                                            computeQdot = .TRUE.)

         ! if ( computeGradients ) then
         !    call elements(eID) % ProlongGradientsToFaces(NGRAD, mesh % faces(meshFaceIDs(1)),&
         !                                     mesh % faces(meshFaceIDs(2)),&
         !                                     mesh % faces(meshFaceIDs(3)),&
         !                                     mesh % faces(meshFaceIDs(4)),&
         !                                     mesh % faces(meshFaceIDs(5)),&
         !                                     mesh % faces(meshFaceIDs(6)) )
         ! end if
!!$omp end task
      end do
!!$omp end single

   End Subroutine SourceProlongSolution

   Subroutine calculateFWHVariables(Q, Qdot, isSolid, Qi, QiDot, Lij, LijDot)

       use VariableConversion, only: Pressure, PressureDot
       use FWHDefinitions,     only: rho0, P0, c0, U0, M0
       implicit none

       real(kind=RP), dimension(NCONS), intent(in)         :: Q        ! horses variables array
       real(kind=RP), dimension(NCONS), intent(in)         :: Qdot     ! horses time derivatives array
       logical, intent(in)                                 :: isSolid
       real(kind=RP), dimension(NDIM), intent(out)         :: Qi       ! fwh Qi array, related with the accoustic pressure thickness
       real(kind=RP), dimension(NDIM), intent(out)         :: Qidot
       real(kind=RP), dimension(NDIM,NDIM), intent(out)    :: Lij      ! fwh Lij tensor: related with the accoustic pressure loading
       real(kind=RP), dimension(NDIM,NDIM), intent(out)    :: LijDot

       !local variables
       real(kind=RP)                                       :: P, pDot
       real(kind=RP), dimension(NDIM,NDIM)                 :: Pij      ! fwh perturbation stress tensor
       ! real(kind=RP), dimension(NDIM:NDIM)                 :: tau
       integer                                             :: i, j, ii, jj

       P = Pressure(Q)
       pDot = PressureDot(Q,Qdot)
       ! print *, "P0: ", P0
       ! print *, "P: ", P
       ! print *, "QDot: ", QDot
       ! print *, "pDot: ", pDot

       LijDot = 0.0_RP
       do i=1,NDIM
           Pij(i,i) = P - P0
           !pressure derivative of LijDot
           LijDot(i,i) = pDot
       end do

       !TODO use the stress tensor and the time derivative for Lij and LijDot respectively
       ! call getStressTensor(Q, U_x, U_y, U_z, tau)
       ! Pij = Pij - tau
       ! LijDot = LijDot - tauDot

       ! set values for solid (impermeable) surface
       Qi(:) = -rho0*U0(:)
       Qidot = 0.0_RP
       Lij = Pij

       !calculate terms for permable surface
       if (.not. isSolid) then
           Qi(:) = Qi(:) + Q(2:4)
           ! convert to complete velocity instead of perturbation velocity
           QiDot(:) = QiDot(:) + Qdot(2:4) - Qdot(1)*U0(:)

           do j = 1, NDIM
               jj = j + 1
               do i = 1, NDIM
                   ! one index is added since rhoV1 = Q(2), rhoV2 = Q(3) ...
                   ii = i + 1
                   Lij(i,j) = Lij(i,j) + (Q(ii) - Q(1)*U0(i))*(Q(jj)/Q(1))
                   LijDot(i,j) = LijDot(i,j) + (1/Q(1)) * &
                                 ( (Qdot(ii) - Qdot(1)) * (Q(jj) - Q(1)*U0(j)) + &
                                 ( Qdot(jj) - Q(jj)/Q(1)*Qdot(1) ) * (Q(ii) - Q(1)*U0(i)) )
               end do  
           end do  
       end if

   End Subroutine calculateFWHVariables

End Module  FWHObseverClass 
