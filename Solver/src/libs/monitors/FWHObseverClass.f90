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
       ! class(Zone_t), pointer                                          :: sourceZone
       logical                                                         :: active
       character(len=STR_LEN_MONITORS)                                 :: observerName
       character(len=STR_LEN_MONITORS)                                 :: fileName

       contains

           procedure :: construct      => ObserverConstruct
           procedure :: destruct       => ObserverDestruct
           procedure :: update         => ObserverUpdate
           ! procedure :: writeToFile    => ObserverWriteToFile
           procedure :: updateTdelay   => ObserverUpdateTdelay

   end type ObserverClass

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
       ! type(Face)                                          :: faceSource
       integer                                             :: faceIDinMesh
       logical                                             :: isSolid

       !TODO: once it works, remove comments of faceSource
       !TODO: remove isSolid as an attribute and passid it as a parameter for FWHSurfaceIntegral (only need there), could be an
       !attribute of the Observer or better of the FW class, and pass as a parameter when the update is called

       contains

           procedure :: construct       => ObserverSourcePairConstruct
           procedure :: destruct        => ObserverSourcePairDestruct
           procedure :: FWHSurfaceIntegral

   end type ObserverSourcePairClass

   contains

!/////////////////////////////////////////////////////////////////////////
!           OBSERVER CLASS PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

   Subroutine ObserverConstruct(self, sourceZone, mesh, ID, solution_file, FirstCall, isSolid)

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
       ! class(Zone_t), intent(in), target                            :: sourceZone
       class(Zone_t), intent(in)                            :: sourceZone
       class(HexMesh), intent(in)                           :: mesh
       integer, intent(in)                                  :: ID
       character(len=*), intent(in)                         :: solution_file
       logical, intent(in)                                  :: FirstCall, isSolid

       ! local variables
       character(len=STR_LEN_MONITORS)  :: in_label
       character(len=STR_LEN_MONITORS)  :: fileName
       character(len=STR_LEN_MONITORS)  :: paramFile
       character(len=STR_LEN_MONITORS)  :: coordinates
       integer                          :: fID
       integer                          :: MeshFaceID, zoneFaceID
       ! integer                          :: number
!
!      Get observer ID
!      --------------
       self % ID = ID
!
!      Get observer zone (surface for integration)
!      --------------
       ! self % sourceZone => sourceZone
!
!      Search for the parameters in the case file
!      ------------------------------------------
       write(in_label , '(A,I0)') "#define accoustic observer " , self % ID

       call get_command_argument(1, paramFile)
       call readValueInRegion(trim ( paramFile), "name",   self % observerName, in_label, "# end" ) 
       call readValueInRegion(trim(paramFile), "position", coordinates        , in_label, "# end" )

!      Get the coordinates
!      -------------------
      print *, "Observer: ", self%observerName
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
          call self % sourcePair(zoneFaceID) % construct(self % x, mesh % faces(MeshFaceID), MeshFaceID, isSolid, FirstCall)
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

   Subroutine ObserverUpdate(self, mesh, BufferPosition)

!     *******************************************************************
!        This subroutine updates the observer accoustic pressure computing it from
!        the mesh storage. It is stored in the "bufferPosition" position of the 
!        buffer.
!         TODO: use openmp (commented)
!         TODO: use mpi (see surface integral)
!     *******************************************************************
!
      use ElementClass
      implicit none
      class (ObserverClass)                                :: self
      class (HexMesh), intent(inout), target               :: mesh
      integer,intent(in)                                   :: bufferPosition

      ! local variables
      real(kind=RP)                                        :: Pt, Pl  ! pressure of each pair
      real(kind=RP), dimension(3)                          :: localPacc, Pacc   ! temporal variable to store the sum of the pressure
      real(kind=RP)                                        :: valx, valy, valz
      integer                                              :: zoneFaceID, meshFaceID, eID !,ierr
      integer, dimension(6)                                :: meshFaceIDs
      class(Element), pointer                              :: elements(:)

!     Initialization
!     --------------            
      self % Pac = 0.0_RP
      Pacc = 0.0_RP
      valx = 0.0_RP
      valy = 0.0_RP
      valz = 0.0_RP

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
      do zoneFaceID = 1, self % numberOfFaces
         meshFaceID = self % sourcePair(zoneFaceID) % faceIDinMesh

         eID = mesh % faces(meshFaceID) % elementIDs(1)
         meshFaceIDs = mesh % elements(eID) % faceIDs

!!$omp task depend(inout:elements(eID))
         call elements(eID) % ProlongSolutionToFaces(NCONS, mesh % faces(meshFaceIDs(1)),&
                                         mesh % faces(meshFaceIDs(2)),&
                                         mesh % faces(meshFaceIDs(3)),&
                                         mesh % faces(meshFaceIDs(4)),&
                                         mesh % faces(meshFaceIDs(5)),&
                                         mesh % faces(meshFaceIDs(6)) )
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

!        Loop the pairs (equivalent to lood the zone) and get the values
!        ---------------------------------------
!!$omp do private(fID,localVal) reduction(+:valx,valy,valz) schedule(runtime)
      do zoneFaceID = 1, self % numberOfFaces
!        Compute the integral
!        --------------------
         ! localPacc = self % sourcePair(zoneFaceID) % FWHSurfaceIntegral()
         meshFaceID = self % sourcePair(zoneFaceID) % faceIDinMesh
         localPacc = self % sourcePair(zoneFaceID) % FWHSurfaceIntegral( mesh % faces(meshFaceID) )
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

      self%Pac(bufferPosition,:) = Pacc

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
          ! for moving surfaces each pair of observer-source should updated the tDelay first
          ! print *, "time of source ", i, "of observer ", self % ID, "is: ", self % sourcePair(i)%tDelay
          t = t + self % sourcePair(i) % tDelay
      end do  

      self % tDelay = t / real(self % numberOfFaces,RP)
          ! print *, "time of average of observer ", self % ID, "is: ", self % tDelay

   End Subroutine ObserverUpdateTdelay

   Subroutine ObserberWriteToFile(self, iter, tsource, no_of_lines)
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
      
   End Subroutine ObserberWriteToFile

   Subroutine ObserverDestruct(self)

        implicit none
        class(ObserverClass)                               :: self

        safedeallocate (self % Pac)
        call self % sourcePair % destruct
        safedeallocate (self % sourcePair)

   End Subroutine ObserverDestruct

!/////////////////////////////////////////////////////////////////////////
!           OBSERVER SOURCE PAIR CLASS PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

   Subroutine  ObserverSourcePairConstruct(self, x, f, fID, isSolid, FirstCall)

       use fluiddata
       implicit none

       class(ObserverSourcePairClass)                      :: self
       real(kind=RP), dimension(NDIM), intent(in)          :: x       ! observer position
       type(face), intent(in)                              :: f    ! source
       integer                                             :: fID
       logical, intent(in)                                 :: FirstCall, isSolid

       ! local variables
       integer                                             :: Nx,Ny
       integer                                             :: i, j
       real(kind=RP)                                       :: fwGamma2, fwGammaInv
       real(kind=RP)                                       :: rho0, c0, P0
       real(kind=RP), dimension(NDIM)                      :: U0, M0

       ! self % faceSource = f
       self % faceIDinMesh = fID
       self % isSolid = isSolid

       Nx = f % Nf(1)
       Ny = f % Nf(2)

       allocate( self % r(Nx,Ny), self % re(Nx,Ny), self % reStar(Nx,Ny) )
       allocate( self % rVect(NDIM,Nx,Ny), self % reUnitVect(NDIM,Nx,Ny) ,self % reStarUnitVect(NDIM,Nx,Ny) )

       call getMeanStreamValues(rho0, U0, M0, c0, P0)

       fwGamma2 = 1.0_RP / (1.0_RP - dimensionless % Mach**2)
       fwGammaInv = sqrt(1.0_RP - dimensionless % Mach**2)
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
       ! print *, "resUnitVect: ", self%reStarUnitVect(:,1,1), "norm: ",norm2(self%reStarUnitVect(:,1,1))
       ! print *, "res: ", self%re(1,1)
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

   ! Function FWHSurfaceIntegral(self) result(Pacc)
   Function FWHSurfaceIntegral(self, f) result(Pacc)

       class(ObserverSourcePairClass)                      :: self
       class(Face), intent(in)                           :: f
       real(kind=RP),dimension(3)                          :: Pacc  ! accoustic pressure values

       ! local variables
       integer                                             :: i, j  ! face indexes
       real(kind=RP), dimension(NDIM)                      :: Qi,QiDot, n
       real(kind=RP), dimension(NDIM,NDIM)                 :: Lij, LijDot
       type(NodalStorage_t), pointer                       :: spAxi, spAeta
       real(kind=RP)                                       :: rho0, c0, P0
       real(kind=RP), dimension(NDIM)                      :: U0, M0
       real(kind=RP)                                       :: Pt, Pl

       ! Initialization
       Pt = 0.0_RP
       Pl = 0.0_RP
       ! spAxi  => NodalStorage(self % faceSource % Nf(1))
       ! spAeta => NodalStorage(self % faceSource % Nf(2))
       spAxi  => NodalStorage(f % Nf(1))
       spAeta => NodalStorage(f % Nf(2))

       call getMeanStreamValues(rho0, U0, M0, c0, P0)

       associate( Q => f % storage(1) % Q )
           associate( Qdot => f % storage(1) % Qdot )

    !           **********************************
    !           Computes the surface integral
    !              I = \int vec{f}·vec{n} * vec{g}·vec{r} dS
    !           **********************************
    !
                do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)

                   n = f % geom % normal(:,i,j)
                   call calculateFWHVariables(Q, Qdot, self%isSolid, Qi, QiDot, Lij, LijDot)

                   ! loading term integrals
                   Pl = Pl +  dot_product(matmul(LijDot,n(:)),self%reUnitVect(:,i,j)) / (self%reStar(i,j) * c0) * &
                             spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
                   Pl = Pl +  dot_product(matmul(Lij,n(:)),self%reStarUnitVect(:,i,j)) / (self%reStar(i,j)**2) * &
                             spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)

                   ! thickness term integrals, only for permable surfaces
                   if (.not. self % isSolid) then
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

      ! get total accoustic pressure as the sum of the two components (the quadraplo terms are being ignored)
       Pacc = (/Pt, Pl, Pt+Pl/)

   End Function FWHSurfaceIntegral 

   Subroutine calculateFWHVariables(Q, Qdot, isSolid, Qi, QiDot, Lij, LijDot)

       use VariableConversion, only: Pressure, PressureDot
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
       real(kind=RP)                                       :: rho0, c0, P0
       real(kind=RP), dimension(NDIM)                      :: U0, M0
       integer                                             :: i, j, ii, jj

       call getMeanStreamValues(rho0, U0, M0, c0, P0)
       P = Pressure(Q)
       pDot = PressureDot(Q,Qdot)

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

   Subroutine  getMeanStreamValues(rho0, U0, M0, c0, P0)

       use fluiddata

       real(kind=RP), intent(out)                          :: rho0, P0, c0
       real(kind=RP), dimension(NDIM), intent(out)         :: U0, M0

       ! local variables
       real(kind=RP)                                       :: theta, phi, U0Magnitud


       theta = refvalues % AOAtheta*(pi/180.0_RP)
       phi   = refvalues % AOAphi*(pi/180.0_RP)

       ! set 1 by default
       ! TODO use values of boundary conditions (inflow if exists or outflow, or set this defaults if not exists)
       U0Magnitud = 1.0_RP
       rho0 = 1.0_RP

       U0(1)  = U0Magnitud*cos(theta)*cos(phi)
       U0(2)  = U0Magnitud*sin(theta)*cos(phi)
       U0(3)  = U0Magnitud*sin(phi)

       M0 = U0 * dimensionless % Mach
       c0 = U0Magnitud / dimensionless % Mach

       ! default initial condition and outflow BC for energy without external pressure
       ! TODO include external pressure
       P0 = 1.0_RP / (dimensionless % gammaM2)
       ! rhoe0 = P0 / thermodynamics%gammaMinus1 + 0.5_RP*rho0*(U0Magnitud**2)


   End Subroutine getMeanStreamValues

End Module  FWHObseverClass 
