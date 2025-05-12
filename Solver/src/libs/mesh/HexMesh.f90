#include "Includes.h"
MODULE HexMeshClass
      use Utilities                       , only: toLower, almostEqual, AlmostEqualRelax
      use SMConstants
      USE MeshTypes
      USE NodeClass
      USE ElementClass
      USE FaceClass
      use FacePatchClass
      USE TransfiniteMapClass
      use SharedBCModule
      use ElementConnectivityDefinitions
      use ZoneClass                       , only: Zone_t, ConstructZones, ReassignZones
      use PhysicsStorage
      use NodalStorageClass
      use MPI_Process_Info
      use MPI_Face_Class
      use FluidData
      use StorageClass
      use FileReadingUtilities            , only: RemovePath, getFileName
      use FTValueDictionaryClass          , only: FTValueDictionary
      use SolutionFile
      use BoundaryConditions,               only: BCs
      use IntegerDataLinkedList           , only: IntegerDataLinkedList_t
      use PartitionedMeshClass            , only: mpi_partition
      use IBMClass
#if defined(NAVIERSTOKES)
      use WallDistance
#endif
#ifdef _HAS_MPI_
      use mpi
#endif
      IMPLICIT NONE

      private
      public      HexMesh
      public      Neighbor_t, NUM_OF_NEIGHBORS

      public      GetOriginalNumberOfFaces
      public      ConstructFaces, ConstructPeriodicFaces
      public      DeletePeriodicMinusFaces, GetElementsFaceIDs
      public      no_of_stats_variables
!
!     ---------------
!     Mesh definition
!     ---------------
!
      type HexMesh
         integer                                   :: numberOfFaces
         integer                                   :: nodeType
         integer                                   :: no_of_elements
         integer                                   :: no_of_allElements
         integer                                   :: no_of_faces
         integer                                   :: dt_restriction = DT_FIXED     ! Time step restriction of last step (DT_FIXED -initial value-, DT_DIFF or DT_CONV)
         integer      , dimension(:), allocatable  :: Nx, Ny, Nz
         integer                                   :: NDOF
         integer,                     allocatable  :: faces_interior(:)
         integer,                     allocatable  :: faces_mpi(:)
         integer,                     allocatable  :: faces_boundary(:)
         integer,                     allocatable  :: elements_sequential(:)
         integer,                     allocatable  :: elements_mpi(:)
         integer, allocatable                      :: HOPRnodeIDs(:)
         character(len=LINE_LENGTH)                :: meshFileName
         type(SolutionStorage_t)                   :: storage              ! Here the solution and its derivative are stored
         type(Node)   , dimension(:), allocatable  :: nodes
         type(Face)   , dimension(:), allocatable  :: faces
         type(Face)   , dimension(:), allocatable  :: mortar_faces         !Mortars for sliding meshes, 4:1 mortars are embeded in faces
         type(Element), dimension(:), allocatable  :: elements
         type(MPI_FacesSet_t)                      :: MPIfaces
         type(MPI_FacesSet_t)                      :: MPImortar_faces
         type(IBM_type)                            :: IBM
         class(Zone_t), dimension(:), allocatable  :: zones
         logical                                   :: child       = .FALSE.         ! Is this a (multigrid) child mesh? default .FALSE.
         logical                                   :: meshIs2D    = .FALSE.         ! Is this a 2D mesh? default .FALSE.
         integer                                   :: dir2D       = 0               ! If it is in fact a 2D mesh, dir 2D stores the global direction IX, IY or IZ
         integer                                   :: dir2D_ctrl  = 0               ! dir2D as in the control file
         logical                                   :: anisotropic = .FALSE.         ! Is the mesh composed by elements with anisotropic polynomial orders? default false
         logical                                   :: ignoreBCnonConformities = .FALSE.
         integer,                     allocatable  :: HO_Elements(:)           !List of elements with polynomial order greater than 1
         integer,                     allocatable  :: LO_Elements(:)           !List of elements with polynomial order less or equal than 1
         integer,                     allocatable  :: HO_FacesInterior(:)      !List of interior faces with polynomial order greater than 1
         integer,                     allocatable  :: HO_FacesBoundary(:)      !List of boundary faces with polynomial order greater than 1
         integer,                     allocatable  :: HO_ElementsMPI(:)        !List of MPI elements with polynomial order greater than 1
         integer,                     allocatable  :: HO_ElementsSequential(:) !List of sequential elements with polynomial order greater than 1
         logical                                   :: nonconforming= .FALSE. 
         logical                                   :: sliding= .FALSE.
         real(kind=RP)                             :: omega 
         !!!!!
         integer :: numBFacePoints
         integer :: n_sliding
         integer :: n_slidingnewnodes 
         integer, allocatable       :: arr1(:)
         integer, allocatable       :: arr2(:)
         integer, allocatable       :: arr3(:)
         integer, allocatable       :: mortararr1(:,:)
         integer, allocatable       :: mortararr2(:,:)
         integer, allocatable       ::face_nodes(:,:)
         integer, allocatable       ::face_othernodes(:,:)
         integer, allocatable       :: Mat(:,:)
         integer, allocatable       :: Connect(:,:,:)
         integer, allocatable       :: rotmortars(:)
         contains
            procedure :: destruct                      => HexMesh_Destruct
            procedure :: Describe                      => HexMesh_Describe
            procedure :: DescribePartition             => DescribeMeshPartition
            procedure :: AllocateStorage               => HexMesh_AllocateStorage
            procedure :: ConstructZones                => HexMesh_ConstructZones
            procedure :: DefineAsBoundaryFaces         => HexMesh_DefineAsBoundaryFaces
            procedure :: CheckIfMeshIs2D               => HexMesh_CheckIfMeshIs2D
            procedure :: CorrectOrderFor2DMesh         => HexMesh_CorrectOrderFor2DMesh
            procedure :: SetConnectivitiesAndLinkFaces => HexMesh_SetConnectivitiesAndLinkFaces
            procedure :: UpdateFacesWithPartition      => HexMesh_UpdateFacesWithPartition
            procedure :: ConstructGeometry             => HexMesh_ConstructGeometry
            procedure :: ProlongSolutionToFaces        => HexMesh_ProlongSolutionToFaces
            procedure :: ProlongGradientsToFaces       => HexMesh_ProlongGradientsToFaces
            procedure :: PrepareForIO                  => HexMesh_PrepareForIO
            procedure :: Export                        => HexMesh_Export
            procedure :: ExportOrders                  => HexMesh_ExportOrders
            procedure :: ExportBoundaryMesh            => HexMesh_ExportBoundaryMesh
            procedure :: SaveSolution                  => HexMesh_SaveSolution
            procedure :: pAdapt                        => HexMesh_pAdapt
            procedure :: pAdapt_MPI                    => HexMesh_pAdapt_MPI
            procedure :: UpdateHOArrays                => HexMesh_UpdateHOArrays
#if defined(NAVIERSTOKES) || defined(INCNS)
            procedure :: SaveStatistics                => HexMesh_SaveStatistics
            procedure :: ResetStatistics               => HexMesh_ResetStatistics
#endif
            procedure :: LoadSolution                  => HexMesh_LoadSolution
            procedure :: LoadSolutionForRestart        => HexMesh_LoadSolutionForRestart
            procedure :: WriteCoordFile
#if defined(ACOUSTIC)
            procedure :: SetUniformBaseFlow            => HexMesh_SetUniformBaseFlow
            ! procedure :: LoadBaseFlowSolution          => HexMesh_LoadBaseFlowSolution
            procedure :: ProlongBaseSolutionToFaces    => HexMesh_ProlongBaseSolutionToFaces
            procedure :: UpdateMPIFacesBaseSolution    => HexMesh_UpdateMPIFacesBaseSolution
            procedure :: GatherMPIFacesBaseSolution    => HexMesh_GatherMPIFacesBaseSolution
#endif
            procedure :: UpdateMPIFacesPolynomial      => HexMesh_UpdateMPIFacesPolynomial
            procedure :: UpdateMPIFacesSolution        => HexMesh_UpdateMPIFacesSolution
            procedure :: UpdateMPIFacesGradients       => HexMesh_UpdateMPIFacesGradients
            procedure :: UpdateMPIFacesAviscflux       => HexMesh_UpdateMPIFacesAviscflux
            procedure :: UpdateMPIFacesMortarflux      => HexMesh_UpdateMPIFacesMortarflux
            procedure :: UpdateMPIFacesGradMortarflux  => HexMesh_UpdateMPIFacesGradMortarflux
            procedure :: GatherMPIFacesSolution        => HexMesh_GatherMPIFacesSolution
            procedure :: GatherMPIFacesGradients       => HexMesh_GatherMPIFacesGradients
            procedure :: GatherMPIFacesAviscFlux       => HexMesh_GatherMPIFacesAviscFlux
            procedure :: GatherMPIFacesMortarFlux      => HexMesh_GatherMPIFacesMortarFlux
            procedure :: GatherMPIFacesGradMortarFlux  => HexMesh_GatherMPIFacesGradMortarFlux
            procedure :: FindPointWithCoords           => HexMesh_FindPointWithCoords
            procedure :: FindPointWithCoordsInNeighbors=> HexMesh_FindPointWithCoordsInNeighbors
            procedure :: ComputeWallDistances          => HexMesh_ComputeWallDistances
            procedure :: ConformingOnZone              => HexMesh_ConformingOnZone
            procedure :: SetStorageToEqn          => HexMesh_SetStorageToEqn
#if defined(INCNS) && defined(CAHNHILLIARD)
            procedure :: ConvertDensityToPhaseFIeld    => HexMesh_ConvertDensityToPhaseField
            procedure :: ConvertPhaseFieldToDensity    => HexMesh_ConvertPhaseFieldToDensity
#endif
            procedure :: RotateMesh                    => HexMesh_RotateMesh      
            procedure :: MarkSlidingElementsRadius     => HexMesh_MarkSlidingElementsRadius
            procedure :: MarkRadius                    => HexMesh_MarkRadius
            procedure :: RotateNodes                   => HexMesh_RotateNodes
            procedure :: Modifymesh                    => HexMesh_Modifymesh
            procedure :: ConstructSlidingMortarsConforming =>  HexMesh_ConstructSlidingMortarsConforming
            procedure :: ConstructMortars              =>HexMesh_ConstructMortars

            procedure :: copy                          => HexMesh_Assign
            generic   :: assignment(=)                 => copy
      end type HexMesh

      integer, parameter :: NUM_OF_NEIGHBORS = 6 ! Hardcoded: Hexahedral conforming meshes
      integer            :: no_of_stats_variables

      TYPE Neighbor_t         ! added to introduce colored computation of numerical Jacobian (is this the best place to define this type??) - only usable for conforming meshes
         INTEGER :: elmnt(NUM_OF_NEIGHBORS+1) ! "7" hardcoded for 3D hexahedrals in conforming meshes (the last one is itself)... This definition must change if the code is expected to be more general
      END TYPE Neighbor_t
!
!     ========
      CONTAINS
!     ========
!
!     -----------
!     Destructors
!     -----------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE HexMesh_Destruct( self )
         IMPLICIT NONE
         CLASS(HexMesh) :: self

         safedeallocate (self % Nx)
         safedeallocate (self % Ny)
         safedeallocate (self % Nz)
!
!        -----
!        Nodes
!        -----
!
         call self % nodes % destruct
         DEALLOCATE( self % nodes )
         safedeallocate (self % HOPRnodeIDs)
!
!        --------
!        Elements
!        --------
!
         call self % elements % destruct
         DEALLOCATE( self % elements )
!
!        -----
!        Faces
!        -----
!
         call self % faces % Destruct
         DEALLOCATE( self % faces )
!
!        -----
!        Zones
!        -----
!
         if (allocated(self % zones)) DEALLOCATE( self % zones )
!
!        ----------------
!        Solution storage
!        ----------------
!
         call self % storage % destruct

         safedeallocate(self % elements_sequential)
         safedeallocate(self % elements_mpi)
         safedeallocate(self % faces_interior)
         safedeallocate(self % faces_mpi)
         safedeallocate(self % faces_boundary)

         safedeallocate(self % HO_Elements)
         safedeallocate(self % LO_Elements)
         safedeallocate(self % HO_FacesInterior)
         safedeallocate(self % HO_FacesBoundary)
         safedeallocate(self % HO_ElementsMPI)
         safedeallocate(self % HO_ElementsSequential)

!
!        ----------------
!        IBM storage
!        ----------------
!
         if( self% IBM% active ) then
            if( self% child ) then
               call self% IBM% destruct( .true. )
            else
               call self% IBM% destruct( .false. )
            end if
         end if
         
      END SUBROUTINE HexMesh_Destruct
!
!     -------------
!     Print methods
!     -------------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintMesh( self )
      IMPLICIT NONE
      TYPE(HexMesh) :: self
      INTEGER ::  k

      PRINT *, "Nodes..."
      DO k = 1, SIZE(self % nodes)
         CALL PrintNode( self % nodes(k), k )
      END DO
      PRINT *, "Elements..."
      DO k = 1, SIZE(self % elements)
         CALL PrintElement( self % elements(k), k )
      END DO
      PRINT *, "Faces..."
      DO k = 1, SIZE(self % faces)
         CALL  self % faces(k) % Print
      END DO

      END SUBROUTINE PrintMesh
!
!////////////////////////////////////////////////////////////////////////
!
      integer function GetOriginalNumberOfFaces(self)
         USE FTMultiIndexTableClass
         USE FTValueClass

         IMPLICIT NONE
         TYPE(HexMesh) :: self

         INTEGER                 :: eID, faceNumber
         INTEGER                 :: faceID
         INTEGER                 :: nodeIDs(8), faceNodeIDs(4), j

         CLASS(FTMultiIndexTable), POINTER  :: table
         CLASS(FTObject), POINTER :: obj
         CLASS(FTValue) , POINTER :: v

         ALLOCATE(table)
         CALL table % initWithSize( SIZE( self % nodes) )

         GetOriginalNumberOfFaces = 0
         DO eID = 1, SIZE( self % elements )

            nodeIDs = self % elements(eID) % nodeIDs
            DO faceNumber = 1, 6
               DO j = 1, 4
                  faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
               END DO

               IF (.not. table % containsKeys(faceNodeIDs) )     THEN
                  GetOriginalNumberOfFaces = GetOriginalNumberOfFaces + 1

                  ALLOCATE(v)
                  CALL v % initWithValue(GetOriginalNumberOfFaces)
                  obj => v
                  CALL table % addObjectForKeys(obj,faceNodeIDs)
                  CALL release(v)
               END IF
            END DO
         END DO

         CALL release(table)

      end function GetOriginalNumberOfFaces

      SUBROUTINE ConstructFaces( self, success, numberOfElements, HorsesMortars, globalToLocalElementID)  !mod
        !
        !     -------------------------------------------------------------
        !     Go through the elements and find the unique faces in the mesh
        !     -------------------------------------------------------------
        !
                 use IntegerArrayLinkedListTable       
                 IMPLICIT NONE
                 TYPE(HexMesh), target   :: self
                 LOGICAL                 :: success
                 INTEGER,optional        :: numberOfElements
                 INTEGER, optional       :: HorsesMortars(6, 6*SIZE( self % elements ))
                 INTEGER, optional       :: globalToLocalElementID(self % no_of_allElements)
        
                 INTEGER                 :: eID, eIDM, faceNumber, faceNumberM
                 INTEGER                 :: faceID
                 INTEGER                 :: nodeIDs(8), faceNodeIDs(4), j, i, k, l
                 INTEGER                 :: MnodeIDs(8), MfaceNodeIDs(4)
                 type(Table_t)           :: table
                 logical                 :: ConformingMesh 
                 INTEGER                 :: nbface, nintface, nmaster, nslave, nintfacec, mID
        
                 ConformingMesh=.TRUE.
                 nbface=0
                 nintface=0
                 nmaster=0
                 nslave=0
                 nintfacec=0
                 if (present(numberOfElements)) ConformingMesh=.FALSE.
                 table = Table_t(size(self % nodes))
        
                 self % numberOfFaces = 0
                 DO eID = 1, SIZE( self % elements )
        
                    nodeIDs = self % elements(eID) % nodeIDs

                    DO faceNumber = 1, 6
        
                     IF (self%elements(eID)%MortarFaces(faceNumber)==0 .OR. &
                     self%elements(eID)%MortarFaces(faceNumber)>=2)   then ! keep the same

                     IF(self%elements(eID)%MortarFaces(faceNumber)>=2 .AND. &
                     .not.present(globalToLocalElementID))  CYCLE  

                     IF(self%elements(eID)%MortarFaces(faceNumber)>=2 .AND. &
                     present(globalToLocalElementID))  THEN  
                        IF(globalToLocalElementID(HorsesMortars(3,(eID*6)-5+faceNumber-1)).NE.-1)  CYCLE
                     END IF
        
                          DO j = 1, 4
                             faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
                          END DO
        
                          faceID = table % ContainsEntry(faceNodeIDs)
                          IF ( faceID .ne. 0 )     THEN
                           nintfacec=nintfacec+1
           !
           !                 --------------------------------------------------------------
           !                 Add this element to the slave side of the face associated with
           !                 these nodes.
           !                 --------------------------------------------------------------
           
                             self % faces(faceID) % elementIDs(2)  = eID
                             self % faces(faceID) % elementSide(2) = faceNumber
                             self % faces(faceID) % FaceType       = HMESH_INTERIOR
                             self % faces(faceID) % rotation       = faceRotation(masterNodeIDs = self % faces(faceID) % nodeIDs, &
                                                                                slaveNodeIDs  = faceNodeIDs                      )
                             !self % faces(faceID) % rotation       = 0
                                                               
                           !write(*,*) 'rotation', self % faces(faceID) % rotation
                          ELSE!
           !                 ------------------
           !                 Construct new face
           !                 ------------------
           !
                             self % numberOfFaces = self % numberOfFaces + 1
                             nintface=nintface+1
                             IF(self % numberOfFaces > SIZE(self % faces))     THEN
        
                                call table % Destruct
                                PRINT *, "Too many faces for # of elements (interior):", self % numberOfFaces, " vs ", SIZE(self % faces)
                                write(*,*) 'nintface=', nintface
                                write(*,*) 'nmaster=' , nmaster 
                                write(*,*) 'nslave=' , nslave
                                write(*,*) 'nint_constructed=' ,nintfacec
                                success = .FALSE.
                                RETURN
                             END IF
        
                             CALL self % faces(self % numberOfFaces) % Construct(ID  = self % numberOfFaces, &
                                                                                nodeIDs = faceNodeIDs, &
                                                                                elementID = eID,       &
                                                                                side = faceNumber)
        
                             self % faces(self % numberOfFaces) % boundaryName = &
                                      self % elements(eID) % boundaryName(faceNumber)
        
                             self % faces(self % numberOfFaces) % IsMortar=0
                             IF(self%elements(eID)%MortarFaces(faceNumber)>=2 .AND. present(globalToLocalElementID))then
                              if (globalToLocalElementID(HorsesMortars(3,(eID*6)-5+faceNumber-1))==-1) THEN 
                                 self % faces(self % numberOfFaces) % FaceType       = HMESH_INTERIOR
                                 self % faces(self % numberOfFaces) % IsMortar=2
                                 self % faces(self % numberOfFaces) % Mortarpos=MODULO(self%elements(eID)%MortarFaces(faceNumber),20)
                                 self % faces(self % numberOfFaces) %elementIDs(2)=self % faces(self % numberOfFaces) %elementIDs(1)
                                 self % faces(self % numberOfFaces) %elementIDs(1)=0
                                 self % faces(self % numberOfFaces) % elementSide(2) = faceNumber
                                 !write(*,*)'masterface in otherprocess'
                              end if
                              END IF 
           
           !                 ----------------------------------------------
           !                 Mark which face is associated with these nodes
           !                 ----------------------------------------------
           
                             call table % AddEntry(faceNodeIDs)
                          END IF
                       END IF 
                       IF (self%elements(eID)%MortarFaces(faceNumber)==1) THEN    !we construct the new big face and the 4 small faces
                        write(*,*) 'big mortar face construction line 431'
                        DO j = 1, 4
                           faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
                        END DO

                          self % numberOfFaces = self % numberOfFaces + 1

                           nmaster=nmaster+1


                          IF(self % numberOfFaces > SIZE(self % faces))     THEN
                             call table % Destruct
                             PRINT *, "Too many faces for # of elements (master):", self % numberOfFaces, " vs ", SIZE(self % faces)
                             write(*,*) 'nintface=', nintface
                             write(*,*) 'nmaster=' , nmaster 
                             write(*,*) 'nslave=' , nslave
                             write(*,*) 'nint_constructed=' ,nintfacec
                             success = .FALSE.
                             RETURN
                          END IF
                          
        
                          CALL self % faces(self % numberOfFaces) % Construct(ID  = self % numberOfFaces, &
                          nodeIDs = faceNodeIDs, &
                          elementID = eID,       &
                          side = faceNumber)

                          self % faces(self % numberOfFaces) % FaceType       = HMESH_INTERIOR
        
                          self % faces(self % numberOfFaces) % boundaryName = &
                                   self % elements(eID) % boundaryName(faceNumber)

                          call table % AddEntry(faceNodeIDs)
   
                          self % faces(self % numberOfFaces) % IsMortar=1
                          mID=self % numberOfFaces
                          allocate(self % faces(mID) % Mortar(4))
                          write(*,*)' allocated mortar 4'
                          self % faces(mID) % Mortar = 0
                          DO l=1, 4  
                             eIDM=HorsesMortars(l + 2, (eID*6)-5 + faceNumber-1)
                             IF ((present(globalToLocalElementID))) THEN 
                              IF (globalToLocalElementID(eIDM)==-1) then
                              !write(*,*)'slaveface in other process'
                               self % faces(mID) % n_mpi_mortar = self % faces(mID) % n_mpi_mortar + 1 
                              cycle
                              end if
                             END IF
                             self % numberOfFaces = self % numberOfFaces + 1
                              nslave=nslave+1
                              self % faces(mID) % Mortar(l) = self % numberOfFaces
                             IF(self % numberOfFaces > SIZE(self % faces))     THEN
                                call table % Destruct
                                PRINT *, "Too many faces for # of elements (slaves):", self % numberOfFaces, " vs ", SIZE(self % faces)
                                write(*,*) 'nintface=', nintface
                                write(*,*) 'nmaster=' , nmaster 
                                write(*,*) 'nslave=' , nslave
                                write(*,*) 'nint_constructed=' ,nintfacec
                                success = .FALSE.
                                RETURN
                             END IF
        
                             IF(eIDM==0)     THEN
                                call table % Destruct
                                PRINT *, "Mortar error:"
                                success = .FALSE.
                                RETURN
                             END IF
        
                             if (faceNumber==1) faceNumberM=2
                             if (faceNumber==2) faceNumberM=1
                             if (faceNumber==3) faceNumberM=5
                             if (faceNumber==4) faceNumberM=6
                             if (faceNumber==5) faceNumberM=3
                             if (faceNumber==6) faceNumberM=4


                             IF (present(globalToLocalElementID)) then
                              MnodeIDs = self % elements(globalToLocalElementID(eIDM)) % nodeIDs
                           ELSE         
                             MnodeIDs = self % elements(eIDM) % nodeIDs
                           END IF
                             DO j = 1, 4
                                MfaceNodeIDs(j) = MnodeIDs(localFaceNode(j,faceNumberM))
                             END DO
                             !construct small slave mortar 
                             !the left side (1) of a small slave mortar is always the big mortar
                             CALL self % faces(self % numberOfFaces) % Construct(ID  = self % numberOfFaces, &
                             nodeIDs = MfaceNodeIDs, &
                             elementID = eID,       &
                             side = faceNumber)
                             IF (present(globalToLocalElementID)) then
                              self % faces(self % numberOfFaces) % boundaryName = &
                                       self % elements(globalToLocalElementID(eIDM)) % boundaryName(faceNumberM)
                              self % faces(self % numberOfFaces) % elementIDs(2)  = globalToLocalElementID(eIDM)
                           ELSE  
                             self % faces(self % numberOfFaces) % boundaryName = &
                                      self % elements(eIDM) % boundaryName(faceNumberM)
                             self % faces(self % numberOfFaces) % elementIDs(2)  = eIDM
                           END IF
                             self % faces(self % numberOfFaces) % elementSide(2) = faceNumberM
                             self % faces(self % numberOfFaces) % FaceType       = HMESH_INTERIOR
                            ! self % faces(self % numberOfFaces) % rotation       = faceRotation(masterNodeIDs = self % faces(faceID) % nodeIDs, &
                             !                                                            slaveNodeIDs  = faceNodeIDs                      )         
                             self % faces(self % numberOfFaces) % rotation       =0
                             self % faces(self % numberOfFaces) % IsMortar=2
                             self % faces(self % numberOfFaces) % Mortarpos=l
                             call table % AddEntry(MfaceNodeIDs)
                             self % faces(mID) % Mortar(l) = self % numberOfFaces
                          END DO !l
                       END IF  
                       
                    END DO !faceNumber
        
                 END DO !eID 
        
                 call table % Destruct

                 write(*,*)"construction done"
        
              END SUBROUTINE ConstructFaces
        
              subroutine GetElementsFaceIDs(self)       !mod
                 implicit none
                 type(HexMesh), intent(inout)  :: self
        !
        !        ---------------
        !        Local variables
        !        ---------------
        !
                 integer  :: fID, eL, eR, e, side, k
                  k=0
                 do fID = 1, size(self % faces)
                    select case (self % faces(fID) % faceType)
                    case (HMESH_INTERIOR)
                       !select case (self % faces(fID) % IsMortar)  
                       if (self % faces(fID) % faceType==HMESH_INTERIOR .AND. self%faces(fID)%IsMortar==0) then !conforming
                          eL = self % faces(fID) % elementIDs(1)
                          eR = self % faces(fID) % elementIDs(2)
                        !  if (self%elements(el)%sliding_newnodes) then 
                        !   write(*,*) 'element', el,'sliding new nodes line 565 hexmesh'
                        !   write(*,*) 'right element eR', eR 
                        !  end if 

                         ! if (.not.(self%elements(el)%sliding_newnodes) .and. (self%elements(el)%sliding)) then 
                         !  write(*,*) 'element', el,'just sliding  line 570 hexmesh'
                         ! write(*,*) 'right eR', eR 
                         ! end if 
                          !if (eR==0) then 
                           !write(*,*)'er=0; fID:', fID 
                           !k=k+1
                          !end if 
                        !write(*,*) eL, eR 
                        !write(*,*) fID 
                          self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
                          self % elements(eR) % faceIDs(self % faces(fID) % elementSide(2)) = fID
                          self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1
                          self % elements(eR) % faceSide(self % faces(fID) % elementSide(2)) = 2
                       end if 
                       if (self % faces(fID) % faceType==HMESH_INTERIOR .AND. self%faces(fID)%IsMortar==1) then !Big master
        
                          eL = self % faces(fID) % elementIDs(1)
                          self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
                          self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1
                       end if 
                       if (self % faces(fID) % faceType==HMESH_INTERIOR .AND. self%faces(fID)%IsMortar==2) then
        
                          eR = self % faces(fID) % elementIDs(2)
                          self % elements(eR) % faceIDs(self % faces(fID) % elementSide(2)) = fID
                          self % elements(eR) % faceSide(self % faces(fID) % elementSide(2)) = 2
        
                       end if  
        
                    case (HMESH_BOUNDARY)
                       eL = self % faces(fID) % elementIDs(1)
                       self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
                       self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1
        
                    case (HMESH_MPI)
                       side = maxloc(self % faces(fID) % elementIDs, 1)
        
                       e = self % faces(fID) % elementIDs(side)
                       self % elements(e) % faceIDs(self % faces(fID) % elementSide(side)) = fID
                       self % elements(e) % faceSide(self % faces(fID) % elementSide(side)) = side
                       self % elements(e) % hasSharedFaces = .true.
        
                    case (HMESH_UNDEFINED)
                       eL = self % faces(fID) % elementIDs(1)
                       self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
                       self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1
        
                    end select
                 end do
                ! write(*,*)'nfaces with er=0:', k
              end subroutine GetElementsFaceIDs
!
!////////////////////////////////////////////////////////////////////////
!
!
!---------------------------------------------------------------------
!! Element faces can be rotated with respect to each other. Orientation
!! gives the relative orientation of the master (1) face to the
!! slave (2) face . In this routine,
!! orientation is measured in 90 degree increments:
!!                   rotation angle = orientation*pi/2
!!
!! As an example, faceRotation = 1 <=> rotating master by 90 deg.
!
      INTEGER pure FUNCTION faceRotation(masterNodeIDs, slaveNodeIDs)
         IMPLICIT NONE
         INTEGER, DIMENSION(4), intent(in) :: masterNodeIDs, slaveNodeIDs !< Node IDs
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, dimension(4), parameter :: NEXTNODE = (/2,3,4,1/)
         INTEGER :: j,l

!
!        Rotate until both first nodes match (each j corresponds to a 90deg rotation)
!        -----------------------------------
         DO j = 1, 4
            IF(masterNodeIDs(1) == slaveNodeIDs(j)) EXIT
         END DO 
!
!        Check whether the orientation is same or opposite
!        -------------------------------------------------
         if ( masterNodeIDS(2) == slaveNodeIDs(NEXTNODE(j)) ) then
            faceRotation = j - 1
         else
            faceRotation = j + 3
         end if

      END FUNCTION faceRotation
!
!////////////////////////////////////////////////////////////////////////
!
   INTEGER pure FUNCTION facerotation4sliding(masterNodeIDs, slaveNodeIDs)
      IMPLICIT NONE 
      INTEGER, DIMENSION(4), intent(in) :: masterNodeIDs, slaveNodeIDs

      !if ((masterNodeIDs(1)==slaveNodeIDs(2)) .AND. (masterNodeIDs(4)==slaveNodeIDs(3))) then 
      !   facerotation4sliding=0
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(3)) .AND. (masterNodeIDs(4)==slaveNodeIDs(4))) then 
      !   facerotation4sliding=1
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(4)) .AND. (masterNodeIDs(4)==slaveNodeIDs(1))) then 
      !   facerotation4sliding=2
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(1)) .AND. (masterNodeIDs(4)==slaveNodeIDs(2))) then 
      !   facerotation4sliding=3
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(4)) .AND. (masterNodeIDs(4)==slaveNodeIDs(3))) then 
      !   facerotation4sliding=4
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(1)) .AND. (masterNodeIDs(4)==slaveNodeIDs(4))) then 
      !   facerotation4sliding=5
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(2)) .AND. (masterNodeIDs(4)==slaveNodeIDs(1))) then 
      !   facerotation4sliding=6
      !elseif ((masterNodeIDs(1)==slaveNodeIDs(3)) .AND. (masterNodeIDs(4)==slaveNodeIDs(2))) then 
      !   facerotation4sliding=7
      !end if 

      !if (((masterNodeIDs(1)==slaveNodeIDs(1)) .AND. (masterNodeIDs(2)==slaveNodeIDs(2))) .OR. &
      !((slaveNodeIDs(1)==masterNodeIDs(1)) .AND. (slaveNodeIDs(2)==masterNodeIDs(2)))) then 
      !   facerotation4sliding=0
     ! elseif (((masterNodeIDs(1)==slaveNodeIDs(2)) .AND. (masterNodeIDs(2)==slaveNodeIDs(3))) .OR. &
      !   ((slaveNodeIDs(1)==masterNodeIDs(2)) .AND. (slaveNodeIDs(2)==masterNodeIDs(3)))) then 
      !   facerotation4sliding=1
      !elseif (((masterNodeIDs(1)==slaveNodeIDs(3)) .AND. (masterNodeIDs(2)==slaveNodeIDs(4))) .OR. &
      !   ((slaveNodeIDs(1)==masterNodeIDs(3)) .AND. (slaveNodeIDs(2)==masterNodeIDs(4)))) then 
      !   facerotation4sliding=2
     ! elseif (((masterNodeIDs(1)==slaveNodeIDs(4)) .AND. (masterNodeIDs(2)==slaveNodeIDs(1))) .OR. &
      !!   ((slaveNodeIDs(1)==masterNodeIDs(4)) .AND. (slaveNodeIDs(2)==masterNodeIDs(1)))) then 
      !   facerotation4sliding=3
      !elseif (((masterNodeIDs(1)==slaveNodeIDs(1)) .AND. (masterNodeIDs(2)==slaveNodeIDs(4))) .OR. &
      !   ((slaveNodeIDs(1)==masterNodeIDs(1)) .AND. (slaveNodeIDs(2)==masterNodeIDs(4)))) then 
      !   facerotation4sliding=4
      !elseif (((masterNodeIDs(1)==slaveNodeIDs(2)) .AND. (masterNodeIDs(2)==slaveNodeIDs(1))) .OR. &
      !   ((slaveNodeIDs(1)==masterNodeIDs(2)) .AND. (slaveNodeIDs(2)==masterNodeIDs(1)))) then 
      !   facerotation4sliding=5
      !elseif (((masterNodeIDs(1)==slaveNodeIDs(3)) .AND. (masterNodeIDs(2)==slaveNodeIDs(2))) .OR. &
      !!   ((slaveNodeIDs(1)==masterNodeIDs(3)) .AND. (slaveNodeIDs(2)==masterNodeIDs(2)))) then 
      !   facerotation4sliding=6
      !elseif (((masterNodeIDs(1)==slaveNodeIDs(4)) .AND. (masterNodeIDs(2)==slaveNodeIDs(3))) .OR. &
      !   ((slaveNodeIDs(1)==masterNodeIDs(4)) .AND. (slaveNodeIDs(2)==masterNodeIDs(3)))) then 
      !   facerotation4sliding=7
      !end if 

      if (((masterNodeIDs(3)==slaveNodeIDs(2)) .AND. (masterNodeIDs(4)==slaveNodeIDs(1))) .OR. &
      ((slaveNodeIDs(2)==masterNodeIDs(3)) .AND. (slaveNodeIDs(1)==masterNodeIDs(4)))) then 
         facerotation4sliding=0
      elseif (((masterNodeIDs(1)==slaveNodeIDs(1)) .AND. (masterNodeIDs(2)==slaveNodeIDs(2))) .OR. &
         ((slaveNodeIDs(3)==masterNodeIDs(3)) .AND. (slaveNodeIDs(4)==masterNodeIDs(4)))) then 
            facerotation4sliding=7
      end if 
   END FUNCTION facerotation4sliding
!
!////////////////////////////////////////////////////////////////////////
!
   INTEGER FUNCTION faceRotationnodes(masterNodeIDs, slaveNodeIDs)
   IMPLICIT NONE
   REAL(kind=RP), DIMENSION(4,3), intent(in) :: masterNodeIDs, slaveNodeIDs !< Node IDs
!
!        ---------------
!        Local variables
!        ---------------
!
   integer, dimension(4), parameter :: NEXTNODE = (/2,3,4,1/)
   INTEGER :: j,l
   integer :: coord 
   logical :: success 
   real(kind=RP) :: tol 
!
!        Rotate until both first nodes match (each j corresponds to a 90deg rotation)
!        -----------------------------------

   coord=0
   success=.FALSE.
   tol=0.00000001_RP
   DO j = 1, 4
      coord=0
     ! write(*,*) ' in facerotationnodes j=', j
      call CompareTwoNodesALL(masterNodeIDs(1,:), slaveNodeIDs(j,:), success, coord, tol)
     ! write(*,*) 'success', success
     ! write(*,*) 'coord', coord
      IF (success ) then 
        ! write(*,*) 'its a success'
         EXIT
      END IF  
     ! IF(masterNodeIDs(1) == slaveNodeIDs(j)) EXIT
   END DO 
!
!        Check whether the orientation is same or opposite
!        -------------------------------------------------
  ! write(*,*) 'coord', coord 
   !if (success) write(*,*) 'success'
   coord=0
   success=.FALSE.
  ! write(*,*) 'line 752 j=',j
   call CompareTwoNodesALL(masterNodeIDs(2,:), slaveNodeIDs(NEXTNODE(j),:), success, coord, tol)
   IF (success ) THEN 
    !if ( masterNodeIDS(2) == slaveNodeIDs(NEXTNODE(j)) ) then
      faceRotationnodes = j - 1
   else
      faceRotationnodes = j + 3
   end if
! write(*,*) 'faceRotationNodes', faceRotationnodes
END FUNCTION faceRotationnodes   
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructPeriodicFaces(self, useRelaxTol)
      USE Physics
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! This subroutine looks for periodic boundary conditions. If they
! are found, periodic+ face is set as an interior face. The slave
! face is the periodic- face and will be deleted in the following
! step.
!-------------------------------------------------------------------
!
!
!--------------------
! External variables
!--------------------
!
      TYPE(HexMesh)              :: self
      LOGICAL, intent(in)        :: useRelaxTol

!
!--------------------
! Local variables
!--------------------
!
!
      REAL(KIND=RP)              :: x1(NDIM), x2(NDIM), edge_length(4), min_edge_length
      LOGICAL                    :: master_matched(4), slave_matched(4), success, found
      INTEGER                    :: coord, slaveNodeIDs(4), localCoord

      INTEGER                    :: i,j,k,l
      integer                    :: zIDplus, zIDMinus, iFace, jFace
      character(len=LINE_LENGTH) :: associatedBname
!
!     --------------------------------------------
!     Loop to find faces with the label "periodic"
!     --------------------------------------------
!
!     -----------------------------
!     Loop zones with BC "periodic"
!     -----------------------------
!

      do zIDPlus = 1, size(self % zones)
!
!        Cycle if the zone is not periodic
!        ---------------------------------
         if ( trim(BCs(zIDPlus) % bc % bcType) .ne. "periodic") cycle
!
!        Cycle if the zone has already marked to be deleted
!        --------------------------------------------------
         if ( self % zones(zIDPlus) % toBeDeleted ) cycle
!
!        Reset the coordinate (changes when changing zones)
!        --------------------------------------------------
         coord = 0
!
!        Get the marker of the associated zone
!        -------------------------------------
         found = .false.
         do zIDMinus = 1, size(self % zones)
            call BCs(zIDPlus) % bc % GetPeriodicPair(associatedBname)
            if ( trim(associatedBname) .eq. trim(self % zones(zIDMinus) % Name) ) then
               found = .true.
               self % zones(zIDMinus) % toBeDeleted = .true.
               exit
            end if
         end do

         if ( .not. found ) then
            print*, 'coupled boundary "',trim(associatedBname),' for boundary "',trim(self % zones(zIDPlus) % Name),'" not found.'
            errorMessage(STD_OUT)
            error stop
         end if
!
!        Loop faces in the periodic+ zone
!        --------------------------------
ploop:   do iFace = 1, self % zones(zIDPlus) % no_of_faces
!
!           Get the face ID
!           ---------------
            i = self % zones(zIDPlus) % faces(iFace)
!
!           Consider only HMESH_UNDEFINED faces
!           -----------------------------------
            if ( (self % faces(i) % faceType .ne. HMESH_UNDEFINED)) cycle ploop
!
!           Loop faces in the periodic- zone
!           --------------------------------
mloop:      do jFace = 1, self % zones(zIDMinus) % no_of_faces
!
!              Get the face ID
!              ---------------
               j = self % zones(zIDMinus) % faces(jFace)
!
!              Consider only HMESH_UNDEFINED faces
!              -----------------------------------
               if ( (self % faces(j) % faceType .ne. HMESH_UNDEFINED)) cycle mloop
!
!              ----------------------------------------------------------------------------------------
!              The index i is a periodic+ face
!              The index j is a periodic- face
!              We are looking for couples of periodic+ and periodic- faces where 2 of the 3 coordinates
!              in all the corners are shared. The non-shared coordinate has to be always the same one.
!              ---------------------------------------------------------------------------------------
!
               master_matched(:)   = .FALSE.     ! True if the master corner finds a partner
               slave_matched(:)    = .FALSE.     ! True if the slave corner finds a partner

               ! compute minimum edge length to make matching tolerance relative to element size. Assumes that the periodic
               ! coordinate of all the nodes not vary significativelly compared to the other two coordinates.
               if (useRelaxTol) then
                   edge_length(1)=NORM2(self % nodes(self % faces(i) % nodeIDs(1)) % x - self % nodes(self % faces(i) % nodeIDs(2)) % x)
                   edge_length(2)=NORM2(self % nodes(self % faces(i) % nodeIDs(2)) % x - self % nodes(self % faces(i) % nodeIDs(3)) % x)
                   edge_length(3)=NORM2(self % nodes(self % faces(i) % nodeIDs(3)) % x - self % nodes(self % faces(i) % nodeIDs(4)) % x)
                   edge_length(4)=NORM2(self % nodes(self % faces(i) % nodeIDs(4)) % x - self % nodes(self % faces(i) % nodeIDs(1)) % x)
                   min_edge_length=minval(edge_length)
               end if

               if ( coord .eq. 0 ) then
!
!                 Check all coordinates
!                 ---------------------
                  do localCoord = 1, 3
                     master_matched = .false.
                     slave_matched = .false.
mastercoord:         DO k = 1, 4
                        x1 = self%nodes(self%faces(i)%nodeIDs(k))%x
slavecoord:             DO l = 1, 4
                           IF (.NOT.slave_matched(l)) THEN
                              x2 = self%nodes(self%faces(j)%nodeIDs(l))%x
                              IF (useRelaxTol) THEN
                                  CALL CompareTwoNodesRelax(x1, x2, master_matched(k), localCoord, min_edge_length)
                              ELSE
                                  CALL CompareTwoNodes(x1, x2, master_matched(k), localCoord)
                              END IF
                              IF (master_matched(k)) THEN
                                 slave_matched(l) = .TRUE.
                                 EXIT  slavecoord
                              ENDIF
                           ENDIF
                        ENDDO    slavecoord
                        IF (.NOT.master_matched(k)) EXIT mastercoord
                     ENDDO mastercoord

                     if ( all(master_matched) ) exit
                  end do

               else
!
!                 Check only the shared coordinates
!                 ---------------------------------
                  DO k = 1, 4
                     x1 = self%nodes(self%faces(i)%nodeIDs(k))%x
                     DO l = 1, 4
                        IF (.NOT.slave_matched(l)) THEN
                           x2 = self%nodes(self%faces(j)%nodeIDs(l))%x
                           IF (useRelaxTol) THEN
                               CALL CompareTwoNodesRelax(x1, x2, master_matched(k), coord, min_edge_length)
                           ELSE
                               CALL CompareTwoNodes(x1, x2, master_matched(k), coord)
                           END IF
                           IF (master_matched(k)) THEN
                              slave_matched(l) = .TRUE.
                              EXIT
                           ENDIF
                        ENDIF
                     ENDDO
                     IF (.NOT.master_matched(k)) EXIT
                  ENDDO

               end if

               IF ( all(master_matched) ) THEN
                  if ( coord .eq. 0 ) coord = localCoord
                  self % faces(i) % boundaryName   = emptyBCName
                  self % faces(i) % elementIDs(2)  = self % faces(j) % elementIDs(1)
                  self % faces(i) % elementSide(2) = self % faces(j) % elementSide(1)
                  self % faces(i) % FaceType       = HMESH_INTERIOR
                  self % elements(self % faces(i) % elementIDs(1)) % boundaryName(self % faces(i) % elementSide(1)) = emptyBCName
                  self % elements(self % faces(i) % elementIDs(2)) % boundaryName(self % faces(i) % elementSide(2)) = emptyBCName
!
!                 To obtain the face rotation, we traduce the right element node IDs to the left
!                 ------------------------------------------------------------------------------
                  do k = 1, 4
                     x1 = self % nodes ( self % faces(i) % nodeIDs(k)) % x
                     do l = 1, 4
                        x2 = self % nodes ( self % faces(j) % nodeIDs(l) ) % x
                        IF (useRelaxTol) THEN
                            CALL CompareTwoNodesRelax(x1, x2, success, coord, min_edge_length)
                        ELSE
                            CALL CompareTwoNodes(x1, x2, success, coord)
                        END IF
                        if ( success ) then
                           slaveNodeIDs(l) = self % faces(i) % nodeIDs(k)
                        end if
                     end do
                  end do
                  self % faces(i) % rotation = faceRotation(self % faces(i) % nodeIDs, &
                                                            slaveNodeIDs)
                  cycle ploop

               ENDIF
            end do   mloop ! periodic- faces
!
!           If the code arrives here, the periodic+ face was not able to find a partner
!           ---------------------------------------------------------------------------
            print*, "When constructing periodic boundary conditions,"
            write(STD_OUT,'(A,I0,A,I0,A,I0)') "Face ",i," in zone ",zIDPlus, &
                  " was not able to find a partner. Element: ", self % faces(i) % elementIDs(1)
            errorMessage(STD_OUT)
            error stop

            end do   ploop    ! periodic+ faces
         end do               ! periodic+ zones

         if ( MPI_Process % isRoot .and. useRelaxTol) print *, "Success: when matching all periodic boundary conditions with relaxed comparison"

      END SUBROUTINE ConstructPeriodicFaces
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CompareTwoNodes(x1, x2, success, coord, tol)
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! Comparison of two nodes. If two of the three coordinates are the
! same, there is success. If there is success, the coordinate which
! is not the same is saved. If the initial value of coord is not 0,
! only that coordinate is checked.
!-------------------------------------------------------------------
!
!
!     --------------------
!     External variables
!     --------------------
!
      REAL(KIND=RP) :: x1(3)
      REAL(KIND=RP) :: x2(3)
      LOGICAL       :: success
      INTEGER       :: coord
      REAL(KIND=RP), OPTIONAL :: tol
!
!     --------------------
!     Local variables
!     --------------------
!
      INTEGER :: i
      INTEGER :: counter

      counter = 0

      IF (coord == 0) THEN

         DO i = 1,3
            write(*,*) 'comparing', x1(i), 'with', x2(i)
            IF ( AlmostEqual( x1(i), x2(i), tol ) ) THEN
               counter = counter + 1
            ELSE
               coord = i
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ELSE

         DO i = 1,3
            IF (i /= coord) THEN
               IF ( AlmostEqual( x1(i), x2(i) ) ) THEN
                  counter = counter + 1
               ENDIF
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ENDIF


      END SUBROUTINE CompareTwoNodes
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CompareTwoNodesALL(x1, x2, success, coord, tol)
         IMPLICIT NONE
   !
   !-------------------------------------------------------------------
   ! Comparison of two nodes. If two of the three coordinates are the
   ! same, there is success. If there is success, the coordinate which
   ! is not the same is saved. If the initial value of coord is not 0,
   ! only that coordinate is checked.
   !-------------------------------------------------------------------
   !
   !
   !     --------------------
   !     External variables
   !     --------------------
   !
         REAL(KIND=RP) :: x1(3)
         REAL(KIND=RP) :: x2(3)
         LOGICAL       :: success
         INTEGER       :: coord
         REAL(KIND=RP), OPTIONAL :: tol
   !
   !     --------------------
   !     Local variables
   !     --------------------
   !
         INTEGER :: i
         INTEGER :: counter
   
         counter = 0
   
         IF (coord == 0) THEN
   
            DO i = 1,3
               !write(*,*) 'comparing', x1(i), 'with', x2(i)
               IF ( AlmostEqual( x1(i), x2(i), tol ) ) THEN
                  counter = counter + 1
               ELSE
                  coord = i
               ENDIF
            ENDDO
   
            IF (counter==3) THEN
               success = .TRUE.
            ELSE
               success = .FALSE.
            ENDIF
   
         ELSE
   
            DO i = 1,3
               IF (i /= coord) THEN
                  IF ( AlmostEqual( x1(i), x2(i) ) ) THEN
                     counter = counter + 1
                  ENDIF
               ENDIF
            ENDDO
   
            IF (counter.ge.2) THEN
               success = .TRUE.
            ELSE
               success = .FALSE.
            ENDIF
   
         ENDIF
   
   
         END SUBROUTINE CompareTwoNodesALL


      SUBROUTINE CompareTwoNodesRelax(x1, x2, success, coord, min_edge_length)
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! Similar to CompareTwoNodes, but the comparison of the two nodes
! is done relaxed by the minimum edge length
!-------------------------------------------------------------------
!     --------------------
!     External variables
!     --------------------
!
      REAL(KIND=RP) :: x1(3)
      REAL(KIND=RP) :: x2(3)
      REAL(KIND=RP) :: min_edge_length
      LOGICAL       :: success
      INTEGER       :: coord
!
!     --------------------
!     Local variables
!     --------------------
!
      INTEGER :: i
      INTEGER :: counter

      counter = 0

      IF (coord == 0) THEN

         DO i = 1,3
            IF ( AlmostEqualRelax( x1(i), x2(i) , min_edge_length ) ) THEN
               counter = counter + 1
            ELSE
               coord = i
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ELSE

         DO i = 1,3
            IF (i /= coord) THEN
               IF ( AlmostEqualRelax( x1(i), x2(i) , min_edge_length ) ) THEN
                  counter = counter + 1
               ENDIF
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ENDIF

      END SUBROUTINE CompareTwoNodesRelax
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DeletePeriodicminusfaces(self)
      use MPI_Face_Class
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! This subroutine looks for periodic boundary conditions. If they
! are found, periodic+ face is set as an interior face. The slave
! face is the periodic- face and will be deleted in the following
! step.
!-------------------------------------------------------------------
!
!
!     --------------------
!     External variables
!     --------------------
!
      TYPE(HexMesh) :: self
!
!     --------------------
!     Local variables
!     --------------------
!
      TYPE(Face),ALLOCATABLE  :: dummy_faces(:)
      INTEGER                 :: i, domain
      INTEGER                 :: iFace, numberOfFaces
      character(len=LINE_LENGTH)    :: bName
      integer                 :: newFaceID(self % numberOfFaces)
!
!     This first loop marks which faces will not be deleted
!     -----------------------------------------------------
      newFaceID = -1
      iFace = 0
      ALLOCATE( dummy_faces(self % numberOfFaces) )
      DO i = 1, self%numberOfFaces
         if ( self % faces(i) % faceType .ne. HMESH_UNDEFINED ) then
            iFace = iFace + 1
            dummy_faces(iFace) = self%faces(i)
            dummy_faces(iFace) % ID = iFace
            newFaceID(i) = iFace
         elseif (.not. self % zones(self % faces(i) % zone) % toBeDeleted) then
            iFace = iFace + 1
            dummy_faces(iFace) = self%faces(i)
            dummy_faces(iFace) % ID = iFace
            newFaceID(i) = iFace
         ENDIF
      ENDDO

      numberOfFaces = iFace

      DEALLOCATE(self%faces)
      ALLOCATE(self%faces(numberOfFaces))

      self%numberOfFaces = numberOfFaces

      DO i = 1, self%numberOfFaces
         self%faces(i) = dummy_faces(i)
      ENDDO
!
!     Update MPI face IDs
!     -------------------
      if ( (MPI_Process % doMPIAction) .and. self % MPIfaces % Constructed ) then
         do domain = 1, MPI_Process % nProcs
            do iFace = 1, self % MPIfaces % faces(domain) % no_of_faces
               self % MPIfaces % faces(domain) % faceIDs(iFace) = newFaceID(self % MPIfaces % faces(domain) % faceIDs(iFace))
            end do
         end do
      end if
!
!     Reassign zones
!     -----------------
      CALL ReassignZones(self % faces, self % zones)

      END SUBROUTINE DeletePeriodicminusfaces
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ProlongSolutionToFaces(self, nEqn, HO_Elements)
         implicit none
         class(HexMesh),    intent(inout) :: self
         integer,           intent(in)    :: nEqn
         logical, optional, intent(in)    :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!     
         integer  :: fIDs(6)
         integer  :: eID, i
         logical :: HOElements

         if (present(HO_Elements)) then
            HOElements = HO_Elements
         else
            HOElements = .false.
         end if

         if (HOElements) then
!$omp do schedule(runtime) private(eID, fIDs)
            do i = 1, size(self % HO_Elements)
               eID = self % HO_Elements(i)
               fIDs = self % elements(eID) % faceIDs
               if (.not.self%sliding) then 
                  if (.not.self%nonconforming) then
                  call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                                     self % faces(fIDs(1)),&
                                                                     self % faces(fIDs(2)),&
                                                                     self % faces(fIDs(3)),&
                                                                     self % faces(fIDs(4)),&
                                                                     self % faces(fIDs(5)),&
                                                                     self % faces(fIDs(6)) )
                  else 
                  call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                                     fFR=self % faces(fIDs(1)),&
                                                                     fBK=self % faces(fIDs(2)),&
                                                                     fBOT=self % faces(fIDs(3)),&
                                                                     fR=self % faces(fIDs(4)),&
                                                                     fT=self % faces(fIDs(5)),&
                                                                     fL=self % faces(fIDs(6)),&
                                                                     faces=self % faces )   
                  end if 
               else 
                  call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                                     fFR=self % faces(fIDs(1)),&
                                                                     fBK=self % faces(fIDs(2)),&
                                                                     fBOT=self % faces(fIDs(3)),&
                                                                     fR=self % faces(fIDs(4)),&
                                                                     fT=self % faces(fIDs(5)),&
                                                                     fL=self % faces(fIDs(6)),&
                                                                     faces=self % mortar_faces )  
               end if 
            end do
!$omp end do
         else
!$omp do schedule(runtime) private(fIDs)
            do eID = 1, size(self % elements)
               fIDs = self % elements(eID) % faceIDs
               if (.not.self%sliding) then 
                  if (.not.self%nonconforming) then
                  call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                                     self % faces(fIDs(1)),&
                                                                     self % faces(fIDs(2)),&
                                                                     self % faces(fIDs(3)),&
                                                                     self % faces(fIDs(4)),&
                                                                     self % faces(fIDs(5)),&
                                                                     self % faces(fIDs(6)) )
                  else 
                  call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                                     fFR=self % faces(fIDs(1)),&
                                                                     fBK=self % faces(fIDs(2)),&
                                                                     fBOT=self % faces(fIDs(3)),&
                                                                     fR=self % faces(fIDs(4)),&
                                                                     fT=self % faces(fIDs(5)),&
                                                                     fL=self % faces(fIDs(6)),&
                                                                     faces=self % faces )
                  end if
               else 
                  call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                                     fFR=self % faces(fIDs(1)),&
                                                                     fBK=self % faces(fIDs(2)),&
                                                                     fBOT=self % faces(fIDs(3)),&
                                                                     fR=self % faces(fIDs(4)),&
                                                                     fT=self % faces(fIDs(5)),&
                                                                     fL=self % faces(fIDs(6)),&
                                                                     faces=self % mortar_faces )
               end if 
            end do
!$omp end do
         end if

      end subroutine HexMesh_ProlongSolutionToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ProlongGradientsToFaces(self, nGradEqn)
         implicit none
         class(HexMesh),   intent(inout)  :: self
         integer,          intent(in)     :: nGradEqn
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fIDs(6)
         integer  :: eID

!$omp do schedule(runtime)
         do eID = 1, size(self % elements)
            fIDs = self % elements(eID) % faceIDs
            if (.not.self%sliding) then 
               if (.not.self%nonconforming) then 
                  call self % elements(eID) % ProlongGradientsToFaces(nGradEqn, &
                                                                  fFR=self % faces(fIDs(1)),&
                                                                  fBK=self % faces(fIDs(2)),&
                                                                  fBOT=self % faces(fIDs(3)),&
                                                                  fR=self % faces(fIDs(4)),&
                                                                  fT=self % faces(fIDs(5)),&
                                                                  fL=self % faces(fIDs(6)))
               else 
                  call self % elements(eID) % ProlongGradientsToFaces(nGradEqn, &
                                                                  fFR=self % faces(fIDs(1)),&
                                                                  fBK=self % faces(fIDs(2)),&
                                                                  fBOT=self % faces(fIDs(3)),&
                                                                  fR=self % faces(fIDs(4)),&
                                                                  fT=self % faces(fIDs(5)),&
                                                                  fL=self % faces(fIDs(6)),&
                                                                  faces=self % faces )
               end if 
            else 
               call self % elements(eID) % ProlongGradientsToFaces(nGradEqn, &
                                                                  fFR=self % faces(fIDs(1)),&
                                                                  fBK=self % faces(fIDs(2)),&
                                                                  fBOT=self % faces(fIDs(3)),&
                                                                  fR=self % faces(fIDs(4)),&
                                                                  fT=self % faces(fIDs(5)),&
                                                                  fL=self % faces(fIDs(6)),&
                                                                  faces=self % mortar_faces )
            end if 
         end do
!$omp end do

      end subroutine HexMesh_ProlongGradientsToFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------------------------------
!     HexMesh_UpdateMPIFacesPolynomial:
!        Send the face polynomial orders to the duplicate faces in the neighbor partitions
!        The information sent per face is:
!              [fNxi, fNeta, eNxi, eNeta, eNzeta, eGlobID],
!        where:
!           fNxi:    Polynomial order 1 for neighbor face
!           fNeta:   Polynomial order 2 for neighbor face
!           eNxi:    Polynomial order 1 for neighbor element
!           eNeta:   Polynomial order 2 for neighbor element
!           eNzeta:  Polynomial order 3 for neighbor element
!           eGlobID: Global ID of neighbor element
!     ------------------------------------------------------------------------------------
      subroutine HexMesh_UpdateMPIFacesPolynomial(self)
         use MPI_Face_Class
         implicit none
         !-arguments----------------------------------------------------------
         class(HexMesh)         :: self
#ifdef _HAS_MPI_
         !-local-variables----------------------------------------------------
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         !--------------------------------------------------------------------

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvN(domain)
         end do
!
!        *************
!        Send solution
!        *************
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate( f => self % faces(fID))
               associate( e => self % elements(maxval(f % elementIDs)) )


               self % MPIfaces % faces(domain) % Nsend(counter:counter+1  ) = e % Nxyz(axisMap(:,f % elementSide(thisSide)))
               self % MPIfaces % faces(domain) % Nsend(counter+2:counter+4) = e % Nxyz
               self % MPIfaces % faces(domain) % Nsend(counter+5)           = e % globID

               counter = counter + 6

               end associate
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendN(domain)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesPolynomial
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateMPIFacesSolution(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvQ(domain, nEqn)
         end do
!
!        *************
!        Send solution
!        *************
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % Qsend(counter:counter+nEqn-1) = f % storage(thisSide) % Q(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendQ(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesSolution

      subroutine HexMesh_UpdateMPIFacesGradients(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvU_xyz(domain, nEqn)
         end do
!
!        ***************
!        Gather gradients
!        ***************
!
         do domain = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            counter = 1

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)

               associate(f => self % faces(fID))

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % U_xyzsend(counter:counter+nEqn-1) = f % storage(thisSide) % U_x(:,i,j)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % U_xyzsend(counter:counter+nEqn-1) = f % storage(thisSide) % U_y(:,i,j)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % U_xyzsend(counter:counter+nEqn-1) = f % storage(thisSide) % U_z(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do

            call self % MPIfaces % faces(domain) % SendU_xyz(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesGradients
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateMPIFacesAviscFlux(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvAviscFlux(domain, nEqn)
         end do
!
!        ***********
!        Send H flux
!        ***********
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % AviscFluxSend(counter:counter+nEqn-1) = &
                     f % storage(thisSide) % AviscFlux(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendAviscFlux(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesAviscFlux

!
!////////////////////////////////////////////////////////////////////////
!
#if defined(ACOUSTIC)
      subroutine HexMesh_UpdateMPIFacesBaseSolution(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvQ(domain, nEqn)
         end do
!
!        *************
!        Send solution
!        *************
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % Qsend(counter:counter+nEqn-1) = f % storage(thisSide) % Qbase(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendQ(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesBaseSolution
#endif
!
      !////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateMPIFacesMortarFlux(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvMortarFlux(domain, nEqn)
         end do
!
!        ***********
!        Send H flux
!        ***********
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               if (f % Ismortar==2) then 
                  do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                     self % MPIfaces % faces(domain) % Flux_M_Send(counter:counter+nEqn-1) = &
                        f % storage(1) % MortarFlux(:,i,j)
                     counter = counter + nEqn
                  end do               ; end do
                  !write(*,*) 'update', f % storage(1) % MortarFlux
               end if 
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendMortarFlux(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesMortarFlux
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateMPIFacesGradMortarFlux(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter, d 
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvMortarGradFlux(domain, nEqn)
         end do
!
!        ***********
!        Send H flux
!        ***********
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather 
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               if (f % Ismortar==2) then 
                  do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1) 
                     self % MPIfaces % faces(domain) % GradFlux_M_Send(counter:counter+nEqn-1) = &
                        f % storage(1) % GradMortarFlux(:,1,i,j)
                     counter = counter + nEqn
                  end do               ; end do            
                  do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1) 
                     self % MPIfaces % faces(domain) % GradFlux_M_Send(counter:counter+nEqn-1) = &
                        f % storage(1) % GradMortarFlux(:,2,i,j)
                     counter = counter + nEqn
                  end do               ; end do  
                  do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1) 
                     self % MPIfaces % faces(domain) % GradFlux_M_Send(counter:counter+nEqn-1) = &
                        f % storage(1) % GradMortarFlux(:,3,i,j)
                     counter = counter + nEqn
                  end do               ; end do  
               end if 
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendMortarGradFlux(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesGradMortarFlux
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_GatherMPIFacesSolution(self, nEqn)
         use FaceClass
         implicit none
         class(HexMesh)    :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForSolution

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % Q(:,i,j) = self % MPIfaces % faces(domain) % Qrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do
               if (f % IsMortar==2) then 
                  call f % Interpolatebig2small(nEqn, f)
               end if 
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesSolution

      subroutine HexMesh_GatherMPIFacesGradients(self, nEqn)
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter, l
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs

!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForGradients

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_x(:,i,j) = self % MPIfaces % faces(domain) % U_xyzrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_y(:,i,j) = self % MPIfaces % faces(domain) % U_xyzrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_z(:,i,j) = self % MPIfaces % faces(domain) % U_xyzrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do

               if (f % IsMortar==2) then 
                  !write(*,*) 'gradient mpi ux before interpolation', f % storage(otherSide(thisSide)) % U_x
                   !write(*,*) 'gradient mpi uy before interpolation', f % storage(otherSide(thisSide)) % U_y
                   !write(*,*) 'gradient mpi uz before interpolation', f % storage(otherSide(thisSide)) % U_z
                   do l=1,3
                      call f % Interpolatebig2small(nEqn, f,grad=l)
                   end do 
 
 
                end if 
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesGradients
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_GatherMPIFacesAviscflux(self, nEqn)
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = [2, 1]

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForAviscflux

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % Aviscflux(:,i,j) = &
                     self % MPIfaces % faces(domain) % AviscfluxRecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesAviscflux
!
      !////////////////////////////////////////////////////////////////////////
!
#if defined(ACOUSTIC)
      subroutine HexMesh_GatherMPIFacesBaseSolution(self, nEqn)
         implicit none
         class(HexMesh)    :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForSolution

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % Qbase(:,i,j) = self % MPIfaces % faces(domain) % Qrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesBaseSolution
#endif
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_GatherMPIFacesMortarFlux(self, nEqn)
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter,p
         integer, parameter :: otherSide(2) = [2, 1]

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForMortarFlux

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
                  if (f%Ismortar==1) then 
                     associate(fStar => f % storage(1) % Fstar)
                        do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                           fStar(:,i,j) = fStar(:,i,j) + &
                              self % MPIfaces % faces(domain) % Flux_M_Recv(counter:counter+nEqn-1)
                           counter = counter + nEqn
                        end do               ; end do
                        !write(*,*) 'fstar mpi', fstar 
                     end associate
            end if 
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesMortarFlux
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_GatherMPIFacesGradMortarFlux(self, nEqn)
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter, d
         integer, parameter :: otherSide(2) = [2, 1]

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForMortarGradFlux

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
                  if (f%Ismortar==1) then 
                     associate(unStar => f % storage(1) % unStar)
                        do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1) 
                           unStar(:,1,i,j) = unStar(:,1,i,j) + &
                              self % MPIfaces % faces(domain) % GradFlux_M_Recv(counter:counter+nEqn-1)
                           counter = counter + nEqn
                        end do               ; end do       
                        do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1) 
                           unStar(:,2,i,j) = unStar(:,2,i,j) + &
                              self % MPIfaces % faces(domain) % GradFlux_M_Recv(counter:counter+nEqn-1)
                           counter = counter + nEqn
                        end do               ; end do  
                        do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1) 
                           unStar(:,3,i,j) = unStar(:,3,i,j) + &
                              self % MPIfaces % faces(domain) % GradFlux_M_Recv(counter:counter+nEqn-1)
                           counter = counter + nEqn
                        end do               ; end do  
                        !write(*,*) unStar
                     end associate
            end if 
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesGradMortarFlux
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_PrepareForIO(self)
         implicit none
         class(HexMesh)    :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr, eID
         integer, allocatable :: elementSizes(:)
         integer, allocatable :: allElementSizes(:)
         integer, allocatable :: allElementsOffset(:)
!
!        Get each element size
!        ---------------------
         allocate(elementSizes(self % no_of_allElements))

         elementSizes = 0     ! Default value 0 to use allreduce with SUM
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            elementSizes(e % globID) = product( e % Nxyz + 1)
            end associate
         end do
!
!        Gather the rest of the mesh values if it is a partition
!        -------------------------------------------------------
         allocate(allElementSizes(self % no_of_allElements))

         allElementSizes = 0
         if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
            call mpi_allreduce(elementSizes, allElementSizes, &
                               self % no_of_allElements, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
         else
            allElementSizes = elementSizes

         end if
!
!        Get all elements offset: the accumulation of allElementSizes
!        -----------------------
         allocate(allElementsOffset(self % no_of_allElements))

         allElementsOffset(1) = 0
         do eID = 2, self % no_of_allElements
            allElementsOffset(eID) = allElementsOffset(eID-1) + allElementSizes(eID-1)
         end do
!
!        Assign the results to the elements
!        ----------------------------------
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            e % offsetIO = allElementsOffset(e % globID)
            end associate
         end do
!
!        Free memory
!        -----------
         deallocate(elementSizes, allElementSizes, allElementsOffset)

      end subroutine HexMesh_PrepareForIO
!
!////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------
!     This subroutine describes the loaded mesh
!     --------------------------------------------------------------------
      SUBROUTINE HexMesh_Describe( self , fileName, bFaceOrder )
      USE Headers
      IMPLICIT NONE
      !-arguments------------------------------------------
      CLASS(HexMesh)      :: self
      CHARACTER(LEN=*)    :: fileName
      integer, intent(in) :: bFaceOrder
      !-local-variables------------------------------------
      integer           :: ierr
      integer           :: zoneID
      integer           :: no_of_bdry_faces
      integer           :: no_of_faces
      integer, allocatable :: facesPerZone(:)
      character(len=LINE_LENGTH) :: str
      !----------------------------------------------------

      allocate ( facesPerZone(size(self % zones)) )

!     Gather information
!     ------------------

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         do zoneID = 1, size(self % zones)
            call mpi_reduce ( self % zones(zoneID) % no_of_faces, facesPerZone(zoneID) , 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
         end do

         no_of_bdry_faces = sum(facesPerZone)
         no_of_faces      = (6*self % no_of_allElements + no_of_bdry_faces)/2
#endif
      else
         do zoneID = 1, size(self % zones)
            facesPerZone(zoneID) = self % zones(zoneID) % no_of_faces
         end do

         no_of_bdry_faces = sum(facesPerZone)
         no_of_faces = size ( self % faces )
      end if


!     Describe the mesh
!     -----------------

      if ( .not. MPI_Process % isRoot ) return

      write(STD_OUT,'(/)')
      call Section_Header("Mesh information")
      write(STD_OUT,'(/)')

      call SubSection_Header('Mesh file "' // trim(fileName) // '".')

      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of elements: " , self % no_of_allElements
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of faces: " , no_of_faces

      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of boundary faces: " , no_of_bdry_faces
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Order of curved faces: " , bFaceOrder
      write(STD_OUT,'(30X,A,A28,L10)') "->" , "2D extruded mesh: " , self % meshIs2D

!     Describe the zones
!     ------------------
      write(STD_OUT,'(/)')
      call Section_Header("Creating zones")
      write(STD_OUT,'(/)')

      do zoneID = 1, size(self % zones)
         write(str,'(A,I0,A,A)') "Zone ", zoneID, " for boundary: ",trim(self % zones(zoneID) % Name)
         call SubSection_Header(trim(str))
         write(STD_OUT,'(30X,A,A28,I0)') "->", ' Number of faces: ', facesPerZone(zoneID)
         call BCs(zoneID) % bc % Describe
         write(STD_OUT,'(/)')
      end do

!     Finish up
!     ---------
      deallocate ( facesPerZone )

      END SUBROUTINE HexMesh_Describe
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DescribeMeshPartition( self )
      USE Headers
      use PartitionedMeshClass
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!  This subroutine describes the loaded mesh partition
!--------------------------------------------------------------------
!
!
!     ------------------
!     External variables
!     ------------------
!
      CLASS(HexMesh)    :: self
#ifdef _HAS_MPI_
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER           :: fID, zoneID, rank, ierr
      INTEGER           :: no_of_bdryfaces, no_of_mpifaces
      integer           :: no_of_elementsP(MPI_Process % nProcs)
      integer           :: no_of_facesP(MPI_Process % nProcs)
      integer           :: no_of_bfacesP(MPI_Process % nProcs)
      integer           :: no_of_mpifacesP(MPI_Process % nProcs)
	  integer           :: no_of_mpiSuggest(2)
	  real(KIND=RP)     :: maxRatio_mpifacesShared
      character(len=64) :: partitionID

      if ( .not. MPI_Process % doMPIAction ) return

      no_of_bdryfaces = 0
      no_of_mpifaces  = 0

      do fID = 1 , size ( self % faces )
         if ( self % faces(fID) % faceType .eq. HMESH_BOUNDARY) then
            no_of_bdryfaces = no_of_bdryfaces + 1
         elseif ( self % faces(fID) % faceType .eq. HMESH_MPI ) then
            no_of_mpifaces = no_of_mpifaces + 1
         end if
      end do
!
!     Share all quantities to the root process
!     ----------------------------------------
      call mpi_gather(size(self % elements) , 1 , MPI_INT , no_of_elementsP , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(size(self % faces)    , 1 , MPI_INT , no_of_facesP    , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(no_of_bdryfaces       , 1 , MPI_INT , no_of_bfacesP   , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(no_of_mpifaces        , 1 , MPI_INT , no_of_mpifacesP , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)

      if ( .not. MPI_Process % isRoot ) return

      write(STD_OUT,'(/)')
      call Section_Header("Mesh partitions")
      write(STD_OUT,'(/)')

      write(STD_OUT,'(10X,A21)', advance='no') "Partitioning method: "

      select case (MPI_Partitioning)
         case (METIS_PARTITIONING)
            write(STD_OUT,'(A)') 'METIS'
         case (SFC_PARTITIONING)
            write(STD_OUT,'(A)') 'Space-filling curve'
      end select
      write(STD_OUT,*)

      do rank = 1, MPI_Process % nProcs

         write(partitionID,'(A,I0)') "Partition ", rank
         call SubSection_Header(trim(partitionID))

         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of elements: " , no_of_elementsP(rank)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of faces: " , no_of_facesP(rank)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of boundary faces: " , no_of_bfacesP(rank)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of mpi faces: " , no_of_mpifacesP(rank)

      end do
	  
	  maxRatio_mpifacesShared = maxval(REAL(no_of_mpifacesP)/REAL(no_of_facesP))
	  no_of_mpiSuggest(1) = FLOOR(0.4*MPI_Process % nProcs/(maxRatio_mpifacesShared))
	  no_of_mpiSuggest(2) = CEILING(0.6*MPI_Process % nProcs/(maxRatio_mpifacesShared))
	  
	  write(STD_OUT,'(/)')
	  call SubSection_Header("MPI Suggestion: ")
	  
	  if (self % NDOF.gt.5000000) then
		  write(STD_OUT,'(30X,A,A28,F4.2)') "->" , "Max mpi faces ratio: " ,maxRatio_mpifacesShared
		  write(STD_OUT,'(30X,A,A28,A10)')  "->" , "Optimum mpi faces ratio: " ,"0.4-0.6"
		  write(STD_OUT,'(30X,A,A28,I4,A3,I4)')"->" , "Suggested number of mpi: ",no_of_mpiSuggest(1)," - ", no_of_mpiSuggest(2)
      else 
	      write(STD_OUT,'(30X,A,A28,A28)') "->" , "NDOF < 5000000 : " ,"Maximize number of OpenMP"
	  end if 
	  
#endif

      END SUBROUTINE DescribeMeshPartition

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE WriteCoordFile(self,nEqn, FileName)
         USE PhysicsStorage
         IMPLICIT NONE
!
!        -----------------------------------------------------------------
!        This subroutine writes a *.coo file containing all the mesh nodes
!        that can be used for eigenvalue analysis using the TAUev code
!        -----------------------------------------------------------------
!
         !--------------------------------------------------------
         CLASS(HexMesh)       :: self        !<  this mesh
         integer              :: nEqn
         CHARACTER(len=*)     :: FileName    !<  ...
         !--------------------------------------------------------
         INTEGER              :: NumOfElem
         INTEGER              :: i, j, k, el, Nx, Ny, Nz, ndof, cooh
         !--------------------------------------------------------

         NumOfElem = SIZE(self % elements)
!
!        ------------------------------------------------------------------------
!        Determine the number of degrees of freedom
!           TODO: Move this to another place if needed in other parts of the code
!        ------------------------------------------------------------------------
!
         ndof = 0
         DO el = 1, NumOfElem
            Nx = self % elements(el) % Nxyz(1)
            Ny = self % elements(el) % Nxyz(2)
            Nz = self % elements(el) % Nxyz(3)
            ndof = ndof + (Nx + 1)*(Ny + 1)*(Nz + 1)*nEqn
         END DO

         OPEN(newunit=cooh, file=FileName, action='WRITE')

         WRITE(cooh,*) ndof, ndim   ! defined in PhysicsStorage
         DO el = 1, NumOfElem
            Nx = self % elements(el) % Nxyz(1)
            Ny = self % elements(el) % Nxyz(2)
            Nz = self % elements(el) % Nxyz(3)
            DO k = 0, Nz
               DO j = 0, Ny
                  DO i = 0, Nx
                     WRITE(cooh,*) self % elements(el) % geom % x(1,i,j,k), &
                                   self % elements(el) % geom % x(2,i,j,k), &
                                   self % elements(el) % geom % x(3,i,j,k)
                  END DO
               END DO
            END DO
         END DO

         CLOSE(cooh)

      END SUBROUTINE WriteCoordFile
!
!//////////////////////////////////////////////////////////////////////////////
!
!     This subroutine checks if the mesh is a 2D extruded mesh and gets the 2D direction
!     in the local frame for each element
!     --------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_CheckIfMeshIs2D(self)
         implicit none
         !-arguments---------------------------------------
         class(HexMesh),   intent(inout) :: self
         !-local-variables---------------------------------
         integer  :: eID, nID, no_of_orientedNodes
         integer  :: dir
         integer  :: ierr                             ! Error for MPI calls
         integer  :: face1Nodes(NODES_PER_FACE)
         integer  :: face2Nodes(NODES_PER_FACE)
         integer  :: no_of_orientedElems(NDIM)
         logical  :: meshExtrudedIn(NDIM)
         logical  :: meshExtrudedInLocal(NDIM)
         real(kind=RP)  :: xNodesF1(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: xNodesF2(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: dx(NDIM,NODES_PER_FACE)
         real(kind=RP), parameter   :: d2D(NDIM,NDIM) = reshape(  (/1._RP, 0._RP, 0._RP, &
                                                                    0._RP, 1._RP, 0._RP, &
                                                                    0._RP, 0._RP, 1._RP /) , (/3,3/) )
         !-------------------------------------------------

         no_of_orientedElems = 0

         elem_loop: do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
!
!           *****************************************
!           Check if the direction is xi (Left,Right)
!           *****************************************
!
!           Get both face node IDs
!           ----------------------
            face1Nodes = e % nodeIDs(localFaceNode(:, ELEFT))
            face2Nodes = e % nodeIDs(localFaceNode(:, ERIGHT))
!
!           Get the nodes coordinates
!           -------------------------
            do nID = 1, NODES_PER_FACE
               xNodesF1(:,nID) = self % nodes(face1Nodes(nID)) % x
               xNodesF2(:,nID) = self % nodes(face2Nodes(nID)) % x

            end do
!
!           Compute the delta x vectors
!           ---------------------------
            dx = xNodesF2 - xNodesF1
!
!           Normalize!
!           ----------
            do nID = 1, NODES_PER_FACE
               dx(:,nID) = dx(:,nID) / norm2( dx(:,nID) )
            end do
!
!           Check how many delta x vectors are parallel to the the x, y or z axis
!           ---------------------------------------------------------------

            do dir = 1, NDIM
               no_of_orientedNodes = 0
               do nID = 1, NODES_PER_FACE
                  if ( almostEqual(abs(dot_product(dx(:,nID),d2D(:,dir))),1.0_RP) ) then
                     no_of_orientedNodes = no_of_orientedNodes + 1
                  end if
               end do

               if ( no_of_orientedNodes .eq. 4 ) then
!
!                 This is (at least one of) the 2D direction(s)
!                 ---------------------------------------------
                  e % dir2D = IX
                  no_of_orientedElems(dir) = no_of_orientedElems(dir) + 1
                  e % globDir(dir) = IX

               end if
            end do

!
!           *****************************************
!           Check if the direction is eta (Front,Back)
!           *****************************************
!
!           Get both face node IDs
!           ----------------------
            face1Nodes = e % nodeIDs(localFaceNode(:, EFRONT))
            face2Nodes = e % nodeIDs(localFaceNode(:, EBACK))
!
!           Get the nodes coordinates
!           -------------------------
            do nID = 1, NODES_PER_FACE
               xNodesF1(:,nID) = self % nodes(face1Nodes(nID)) % x
               xNodesF2(:,nID) = self % nodes(face2Nodes(nID)) % x
            end do
!
!           Compute the delta x vectors
!           ---------------------------
            dx = xNodesF2 - xNodesF1
!
!           Normalize!
!           ----------
            do nID = 1, NODES_PER_FACE
               dx(:,nID) = dx(:,nID) / norm2( dx(:,nID) )
            end do
!
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            do dir = 1, NDIM
               no_of_orientedNodes = 0
               do nID = 1, NODES_PER_FACE
                  if ( almostEqual(abs(dot_product(dx(:,nID),d2D(:,dir))),1.0_RP) ) then
                     no_of_orientedNodes = no_of_orientedNodes + 1
                  end if
               end do

               if ( no_of_orientedNodes .eq. 4 ) then
!
!                 This is (at least one of) the 2D direction(s)
!                 ---------------------------------------------
                  if (e % dir2D == IX) then
                     e % dir2D = IXY
                  else
                     e % dir2D = IY
                  end if
                  no_of_orientedElems(dir) = no_of_orientedElems(dir) + 1
                  e % globDir(dir) = IY

               end if

            end do
!
!           *****************************************
!           Check if the direction is zeta (Bottom,Top)
!           *****************************************
!
!           Get both face node IDs
!           ----------------------
            face1Nodes = e % nodeIDs(localFaceNode(:, EBOTTOM))
            face2Nodes = e % nodeIDs(localFaceNode(:, ETOP))
!
!           Get the nodes coordinates
!           -------------------------
            do nID = 1, NODES_PER_FACE
               xNodesF1(:,nID) = self % nodes(face1Nodes(nID)) % x
               xNodesF2(:,nID) = self % nodes(face2Nodes(nID)) % x
            end do
!
!           Compute the delta x vectors
!           ---------------------------
            dx = xNodesF2 - xNodesF1
!
!           Normalize!
!           ----------
            do nID = 1, NODES_PER_FACE
               dx(:,nID) = dx(:,nID) / norm2( dx(:,nID) )
            end do
!
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            do dir = 1, NDIM
               no_of_orientedNodes = 0
               do nID = 1, NODES_PER_FACE
                  if ( almostEqual(abs(dot_product(dx(:,nID),d2D(:,dir))),1.0_RP) ) then
                     no_of_orientedNodes = no_of_orientedNodes + 1
                  end if
               end do

               if ( no_of_orientedNodes .eq. 4 ) then
!
!                 This is (at least one of) the 2D direction(s)
!                 ---------------------------------------------
                  select case (e % dir2D)
                     case (IX)
                        e % dir2D = IXZ
                     case (IY)
                        e % dir2D = IYZ
                     case (IXY)
                        e % dir2D = IXYZ
                     case default
                        e % dir2D = IZ
                  end select

                  no_of_orientedElems(dir) = no_of_orientedElems(dir) + 1
                  e % globDir(dir) = IZ

               end if
            end do

            end associate
         end do elem_loop

         meshExtrudedIn = ( no_of_orientedElems == self % no_of_elements )

!        MPI communication
!        -----------------
#if _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            meshExtrudedInLocal = meshExtrudedIn
            call mpi_allreduce ( meshExtrudedInLocal, meshExtrudedIn, NDIM, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr )
         end if
#endif

         if ( any(meshExtrudedIn) ) then
            self % meshIs2D = .TRUE.

            if ( all(meshExtrudedIn) ) then
               self % dir2D = IXYZ
            elseif ( meshExtrudedIn(IX) .and. meshExtrudedIn(IY) ) then
               self % dir2D = IXY
            elseif ( meshExtrudedIn(IX) .and. meshExtrudedIn(IZ) ) then
               self % dir2D = IXZ
            elseif ( meshExtrudedIn(IY) .and. meshExtrudedIn(IZ) ) then
               self % dir2D = IYZ
            elseif ( meshExtrudedIn(IX) ) then
               self % dir2D = IX
            elseif ( meshExtrudedIn(IY) ) then
               self % dir2D = IY
            elseif ( meshExtrudedIn(IZ) ) then
               self % dir2D = IZ
            end if
         end if

      end subroutine HexMesh_CheckIfMeshIs2D
!
!//////////////////////////////////////////////////////////////////////////////
!
!     If the mesh is a 2D extruded mesh, this subroutine sets the polynomial order
!     to zero in the corresponding direction (HexMesh_CheckIfMeshIs2D must have been called beforehand)
!     --------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_CorrectOrderFor2DMesh(self, dir2D,order_)
         implicit none
         class(HexMesh),   intent(inout) :: self
         integer,          intent(in)    :: dir2D
         integer, intent(in), optional   :: order_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, nID, no_of_orientedNodes, order
         integer  :: face1Nodes(NODES_PER_FACE)
         integer  :: face2Nodes(NODES_PER_FACE)
         logical  :: rightDir
         real(kind=RP)  :: d2D(NDIM)
         real(kind=RP)  :: xNodesF1(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: xNodesF2(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: dx(NDIM,NODES_PER_FACE)

         if ( present(order_) ) then
            order = order_
         else
            order = 0
         end if

         if (self % meshIs2D) then
            select case (dir2D)
               case (IX)
                  rightDir = any (self % dir2D == [IX, IXY, IXZ, IXYZ] )
               case (IY)
                  rightDir = any (self % dir2D == [IY, IXY, IYZ, IXYZ] )
               case (IZ)
                  rightDir = any (self % dir2D == [IZ, IXZ, IYZ, IXYZ] )
            end select
            if (.not. rightDir) then
               print*, "The mesh does not seem to be 2D for the selected direction"
               errorMessage(STD_OUT)
               return
            end if
         else
            print*, "The mesh does not seem to be 2D"
            errorMessage(STD_OUT)
         end if

         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))

            select case (e % globDir(dir2D))
               case (IX)
                  e % Nxyz(1) = order
                  self % Nx(eID) = order
               case (IY)
                  e % Nxyz(2) = order
                  self % Ny(eID) = order
               case (IZ)
                  e % Nxyz(3) = order
                  self % Nz(eID) = order
            end select

            end associate
         end do

      end subroutine HexMesh_CorrectOrderFor2DMesh
!
!////////////////////////////////////////////////////////////////////////
!
!        Set element connectivities
!        --------------------------
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_SetConnectivitiesAndLinkFaces(self,nodes,facesList)   !******
        implicit none
        !-arguments----------------------------------------------------
        class(HexMesh), target, intent(inout) :: self
        integer               , intent(in)    :: nodes
        integer, optional     , intent(in)    :: facesList(:)
        !-local-variables----------------------------------------------
        integer  :: fID, SideL, SideR, side, counter, l
        integer  :: NelL(2), NelR(2)    ! Polynomial orders on left and right of a face
        integer  :: Nel(2,2)            ! Polynomial orders on left and right of a face - for MPI
        integer  :: globID
        integer  :: Nxyz(NDIM)
        integer  :: domain, MPI_NDOFS(MPI_Process % nProcs), mpifID
        integer  :: MPI_MNDOFS(MPI_Process % nProcs)
        integer  :: num_of_Faces, ii
        integer, parameter :: other(2) = [2, 1]
        real(kind=RP) :: offset(2)
        !--------------------------------------------------------------

        if ( present(facesList) ) then
           num_of_Faces = size(facesList)
        else
           num_of_Faces = size(self % faces)
        end if
!
!        ---------------------------------------------------------------------------
!        Send polynomial orders across MPI faces
!           -> Currently doing for all MPI faces (not taking into account facesList)
!        ---------------------------------------------------------------------------
!
        if (mpi_partition % Constructed) call self % UpdateMPIFacesPolynomial
!
!        ------------------------
!        Link faces with elements
!        ------------------------
!
        do ii = 1, num_of_Faces

           if ( present(facesList) ) then
              fID = facesList(ii)
           else
              fID = ii
           end if
           associate (  f => self % faces(fID)   )

           select case (f % faceType)
           case (HMESH_INTERIOR)
              select case (f % IsMortar)
              case (0, 2)   
               !write(*,*) 'line 2664 f % IsMortar', f % IsMortar,'fID', f%ID, 'f % elementIDs(1)',f % elementIDs(1),'f % elementIDs(2)',f % elementIDs(2)
               if (f % elementIDs(1)==0) then 
                  associate(eR => self % elements(f % elementIDs(2))   )
                     NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
                     NelL = NelR
                  end associate 
               end if 
               if (f % elementIDs(2)==0) then 
                  associate(eL => self % elements(f % elementIDs(1))   )
                     NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
                     NelR = NelL
                  end associate 
                  end if 
               if ((f % elementIDs(1) .ne. 0) .and.  (f % elementIDs(2) .ne. 0) ) then 
                 associate(eL => self % elements(f % elementIDs(1)), &
                             eR => self % elements(f % elementIDs(2))   )
!
!                 Get polynomial orders of elements
!                 ---------------------------------
                 NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
                 NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
!
!                 Fill connectivity of element type
!                 ---------------------------------
                 SideL = f % elementSide(1)
                 SideR = f % elementSide(2)
!
!                 Construct connectivity
!                 ----------------------

                 if (NelL(1)==-1 .or. NelL(2)==-1) then 
                  write(*,*) 'line 2681; NelL of face', f%id,'=', NelL,'el', el%eID,'ftype', f%facetype 
                 end if 
                 if (NelR(1)==-1 .or. NelR(2)==-1) then 
                  write(*,*) 'line 2683; NelR of face', f%id,'=', NelR,'er', er%eID,'ftype', f%facetype
                 end if 
                 eL % NumberOfConnections(SideL) = 1
                 call eL % Connection(SideL) % Construct(eR % GlobID, eR % Nxyz)

                 eR % NumberOfConnections(SideR) = 1
                 call eR % Connection(SideR) % Construct(eL % GlobID, eL % Nxyz)
                 end associate
               end if 
                 if (f % Ismortar==2) then 
                 offset(1)=0.5_RP!!!!!!!if
                 offset(2)=0.5_RP
                  call f % LinkWithElements(NelL, NelR, nodes, offset)!!!!!!!
                   !write(*,*) NelL, NelR 
                  !write(*,*) 'after construction :', f % NelLeft, f % NelLeft
                                    !write(*,*) 'li,ki,g mortar slave face'
                 elseif (f % Ismortar ==0) then 
                 call f % LinkWithElements(NelL, NelR, nodes)!!!!!!!
                 end if 
               
              case (1, 3) 
               !if (f % Ismortar==3) write(*,*) 'Ismorta=3 2699 hex'
                 associate(eL => self % elements(f % elementIDs(1)))
!
!                 Get polynomial orders of elements
!                 ---------------------------------
                 NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
                 NelR = NelL
                  !if (f %Ismortar==3) write(*,*) 'mortar 3, face', f%ID
                 ! NumberOfConnections = 4 ????
                 ! write(*,*)'NELL 2708', NelL 
                 end associate
                  if (NelL(1)==-1 .or. NelL(2)==-1) then 
                     if (f%IsMortar==3) then 
                        associate(eL => self % elements(f % elementIDs(2)))
                           NelL = eL % Nxyz(axisMap(:, f % elementSide(2)))
                           NelR = NelL
                        end associate 
                           write(*,*) 'now Nell and NelR',NelL, NelR 
                     end if 
                  end if 
                 call f % LinkWithElements(NelL, NelR, nodes)   
                 !if (f %Ismortar==3) write(*,*) ' line 2720 hexmesh mortar 3, face', f%ID, 'Nf',f%Nf


              end select  

           case (HMESH_BOUNDARY)
              associate(eL => self % elements(f % elementIDs(1)))
!
!              Get polynomial orders of elements
!              ---------------------------------
              NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
              NelR = NelL

              ! Default NumberOfConnections = 0
              end associate

              call f % LinkWithElements(NelL, NelR, nodes)

!           case (HMESH_MPI): Do nothing. MPI faces are constructed in the next step.

           end select


           end associate
        end do
      
        !!!!only for sliding mortars 
       ! if (self%sliding) then 
        ! do ii = 1, size(self%mortar_faces)
           ! associate (  f => self % mortar_faces(ii)   )
              ! associate(eL => self % elements(f % elementIDs(1)), &
                 ! eR => self % elements(f % elementIDs(2))   )
   !
   !                 Get polynomial orders of elements
   !                 ---------------------------------
              ! NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
              ! NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
              ! call f % LinkWithElements(NelL, NelR, nodes, f%offset, f%s)
               !end associate 
            !end associate 
         !end do 
        !end if 

!
!        -----------------------------------------------
!        Gather faces polynomial and link MPI faces
!        -> All, not only the faces included in faceList
!        -----------------------------------------------
!
#ifdef _HAS_MPI_
        if ( MPI_Process % doMPIAction .and. mpi_partition % Constructed  )  then

           do domain = 1, MPI_Process % nProcs
!
!              Wait until messages have been received
!              --------------------------------------
!
              call self % MPIfaces % faces(domain) % WaitForN

              counter = 1
              do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
                 fID  = self % MPIfaces % faces(domain) % faceIDs(mpifID)
                 side = self % MPIfaces % faces(domain) % elementSide(mpifID)   ! face side 1/2

                 associate( f => self % faces(fID) )
                 associate( e => self % elements(maxval(f % elementIDs)) )
                  select case (f % IsMortar)
                  case (0, 2)   
                  sideL = f % elementSide(side)                                  ! element side 1/2/3/4/5/6
                  Nel(:,      side ) = e % Nxyz(axisMap(:,sideL))
                  Nel(:,other(side)) = self % MPIfaces % faces(domain) % Nrecv(counter:counter+1)
                  if (f % IsMortar == 2) then 
                   offset=0.5_RP
                   call f % LinkWithElements(Nel(:,1), Nel(:,2), nodes, offset)
                  elseif (f % IsMortar == 0 ) then 
                   call f % LinkWithElements(Nel(:,1), Nel(:,2), nodes)
                  end if 
                   Nxyz   = self % MPIfaces % faces(domain) % Nrecv(counter+2:counter+4)
                  globID = self % MPIfaces % faces(domain) % Nrecv(counter+5)
 
                  e % NumberOfConnections (sideL) = 1
                  call e % Connection(sideL) % construct (globID,Nxyz)
                  counter = counter + 6
               case (1)
                 sideL = f % elementSide(side)                                  ! element side 1/2/3/4/5/6

                 Nel(:,      side ) = e % Nxyz(axisMap(:,sideL))
                 Nel(:,other(side))= Nel(:,      side ) 
                 call f % LinkWithElements(Nel(:,1), Nel(:,2), nodes)

                 Nxyz   = self % MPIfaces % faces(domain) % Nrecv(counter+2:counter+4)
                 globID = self % MPIfaces % faces(domain) % Nrecv(counter+5)

                 e % NumberOfConnections (sideL) = 1
                 call e % Connection(sideL) % construct (globID,Nxyz)
               end select 
                 end associate
                 end associate
              end do
           end do
        end if
#endif
!
!        --------------------------
!        Allocate MPI Faces storage
!           TODO: This can be optimized when facesList is present
!        --------------------------
!
        if ( MPI_Process % doMPIAction ) then
           if ( .not. allocated(self % MPIfaces % faces) ) return
#if _HAS_MPI_
           MPI_NDOFS = 0
           MPI_MNDOFS = 0
           do domain = 1, MPI_Process % nProcs
              do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
                 fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
                 associate( fc => self % faces(fID) )
                 MPI_NDOFS(domain) = MPI_NDOFS(domain) + product(fc % Nf + 1)
                 if (fc % IsMortar == 1 .OR. fc % IsMortar == 2) then 
                  MPI_MNDOFS(domain) = MPI_MNDOFS(domain) + product(fc % Nf + 1)
              end if 
                 end associate
              end do
           end do

#if defined(NAVIERSTOKES)
if (.not.self % nonconforming) then 
           call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NGRAD, MPI_NDOFS)
         else 
            call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NGRAD, MPI_NDOFS, MPI_MNDOFS)
         end if 
#elif defined(INCNS)
if (.not.self % nonconforming) then 
           call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NCONS, MPI_NDOFS)
         else 
            call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NCONS, MPI_NDOFS, MPI_MNDOFS)
         end if
#elif defined(CAHNHILLIARD)
if (.not.self % nonconforming) then 
           call ConstructMPIFacesStorage(self % MPIfaces, NCOMP, NCOMP, MPI_NDOFS)
         else 
            call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NCONS, MPI_NDOFS, MPI_MNDOFS)
         end if 
#endif

#endif
        end if

     end subroutine HexMesh_SetConnectivitiesAndLinkFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateFacesWithPartition(self, partition, nAllElems, global2LocalIDs)
!
!        **************************************************
!        This subroutine casts HMESH_UNDEFINED faces which
!        are a partition boundary face as HMESH_MPI
!        **************************************************
!
         use PartitionedMeshClass
         implicit none
         class(HexMesh)    :: self
         class(PartitionedMesh_t),  intent(in)  :: partition
         integer,                   intent(in)  :: nAllElems
         integer,                   intent(in)  :: global2LocalIDs(nAllElems)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: k, eID, bFace, side, eSide, fID, domain
         integer  :: no_of_mpifaces(MPI_Process % nProcs)
         integer, parameter  :: otherSide(2) = (/2,1/)
         integer, parameter  :: invRot(1:4,0:7) = reshape( (/ 1, 2, 3, 4, &
                                                              4, 1, 2, 3, &
                                                              3, 4, 1, 2, &
                                                              2, 3, 4, 1, &
                                                              1, 4, 3, 2, &
                                                              2, 1, 4, 3, &
                                                              3, 2, 1, 4, &
                                                              4, 3, 2, 1 /), (/4,8/) )
!
!        First get how many faces are shared with each other partition
!        -------------------------------------------------------------
         no_of_mpifaces = 0
         do bFace = 1, partition % no_of_mpifaces
            domain = partition % mpiface_sharedDomain(bFace)
            no_of_mpifaces(domain) = no_of_mpifaces(domain) + 1
         end do
!
!        ---------------
!        Allocate memory
!        ---------------
!
                  !write(*,*) partition % no_of_mpifaces
         do domain = 1, MPI_Process % nProcs
            if ( no_of_mpifaces(domain) .ne. 0 ) then

               call self % MPIfaces % faces(domain) % Construct(no_of_mpifaces(domain))
            end if
         end do
                  !write(*,*) no_of_mpifaces

!
!        -------------
!        Assign values
!        -------------
!
         no_of_mpifaces = 0
         do bFace = 1, partition % no_of_mpifaces
!
!           Gather the face, and the relative position w.r.t. its element
!           -------------------------------------------------------------
            eID = global2LocalIDs(partition % mpiface_elements(bFace))
            side = partition % element_mpifaceSide(bFace)
            eSide = partition % mpiface_elementSide(bFace)
            fID = self % elements(eID) % faceIDs(side)

!
!           Change the face to a HMESH_MPI
!           ------------------------------
            associate(f => self % faces(fID))
            f % faceType = HMESH_MPI
            f % rotation = partition % mpiface_rotation(bFace)

            f % elementIDs(eSide) = eID

            f % elementIDs(otherSide(eSide)) = HMESH_NONE   ! This makes sense since elementIDs are in local numbering...
            f % elementSide(eSide) = side
            f % elementSide(otherSide(eSide)) = partition % element_mpifaceSideOther(bFace)
            if (eSide == RIGHT) f % nodeIDs = f % nodeIDs (invRot (:,f % rotation) )
            end associate

            self % elements(eID) % faceSide(side) = eSide

!
!           Create MPI Face
!           ---------------
            domain = partition % mpiface_sharedDomain(bFace)
            no_of_mpifaces(domain) = no_of_mpifaces(domain) + 1
            self % MPIfaces % faces(domain) % faceIDs(no_of_mpifaces(domain)) = fID
            self % MPIfaces % faces(domain) % elementSide(no_of_mpifaces(domain)) = eSide
               write(*,*) 'curent partition:', partition%ID 
               write(*,*) 'shared partition:', domain
                        !if (f % IsMortar==1) then 
              ! write(*,*) 'MPI Big Mortar fID', fID 
            !elseif (f % IsMortar==2) then 
              ! write(*,*) 'MPI slave Mortar fID', fID 
            !end if
         end do
                 write(*,*) 'number mpifaces', no_of_mpifaces
      
      end subroutine HexMesh_UpdateFacesWithPartition
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------
!     Construct geometry of faces and elements
!     -> This routine guarantees that the mapping is subparametric or isoparametric
!     -> TODO: Additional considerations are needed for p-nonconforming representations with inner curved faces
!     -----------------------------------------------------------------------------
      subroutine HexMesh_ConstructGeometry(self, facesList, elementList)
         implicit none
         !--------------------------------
         class(HexMesh)    , intent(inout) :: self
         integer, optional , intent(in)    :: facesList(:)
         integer, optional , intent(in)    :: elementList(:)
         !--------------------------------
         integer                       :: num_of_elements, num_of_faces
         integer                       :: i, ii
         integer                       :: fID, eID                        ! Face and element counters
         integer                       :: eIDLeft, eIDRight, e            ! Element IDs on left and right of a face
         integer                       :: SideIDL, SideIDR, side          ! Side of elements on left and right of a face
         integer                       :: buffer                          ! A temporal variable
         integer                       :: NSurfL(2), NSurfR(2), NSurf(2)  ! Polynomial order the face was constructef with
         integer                       :: Nelf(2), Nel(2), rot
         integer                       :: CLN(2)                          ! Chebyshev-Lobatto face orders
         REAL(KIND=RP)                 :: corners(NDIM,NODES_PER_ELEMENT) ! Variable just for initializing purposes
         real(kind=RP)   , allocatable :: faceCL(:,:,:)                   ! Coordinates of the Chebyshev-Lobatto nodes on the face
         type(SurfInfo_t), allocatable :: SurfInfo(:)                     ! Local copy of surf info that can be modified
         type(TransfiniteHexMap), pointer :: hexMap, hex8Map, genHexMap
         !--------------------------------
         logical                    :: isConforming   ! Is the representation conforming on a boundary?
         integer                    :: zoneID, zonefID, nZones
         integer, allocatable       :: bfOrder_local(:)  ! Polynomial order on a zone (partition wise)
         integer, allocatable       :: bfOrder(:)        ! Polynomial order on a zone (global)
         integer                    :: ierr              ! Error for MPI calls
         integer                    :: thisSide          ! The side of the MPI face that corresponds to current partition
         !--------------------------------

         corners = 0._RP

         if ( present(elementList) ) then
            num_of_elements = size(elementList)
         else
            num_of_elements = size(self % elements)
         end if
         if ( present(facesList) ) then
            num_of_faces = size(facesList)
         else
            num_of_faces = size(self % faces)
         end if

!
!        Generate a local copy of SurfInfo
!        ---------------------------------
         allocate ( SurfInfo(self % no_of_elements) )
         do ii=1, num_of_elements
            if ( present(elementList) ) then
               eID = elementList(ii)
            else
               eID = ii
            end if

            SurfInfo(eID) = self % elements(eID) % SurfInfo
         end do
!
!        ***********************************************************************
!        Find the polynomial order of the boundaries for 3D nonconforming meshes
!        -> In 3D meshes we force all the faces on a boundary to be mapped with the
!           same polynomial order, in order to have "water-tight" meshes. This condition
!           can be sometimes too strict...
!        -> When the representation is p-nonconforming on a boundary, the mapping
!           must be of order P <= min(N)/2. This is the general sufficient condition.
!        *************************************************************1**********
!
         if (self % anisotropic .and. (.not. self % meshIs2D) ) then
            nZones = size(self % zones)
            allocate ( bfOrder (nZones) )
            bfOrder = huge(bfOrder) ! Initialize to a big number

            do zoneID=1, nZones
               if (self % zones(zoneID) % no_of_faces == 0 ) cycle

!              Get the minimum polynomial order in this zone
!              ---------------------------------------------
               do zonefID = 1, self % zones(zoneID) % no_of_faces
                  fID = self % zones(zoneID) % faces(zonefID)

                  associate( f => self % faces(fID) )
                  bfOrder(zoneID) = min(bfOrder(zoneID),f % NfLeft(1),f % NfLeft(2))
                  end associate
               end do

!           MPI communication
!           -----------------
#if _HAS_MPI_
            end do

            if ( MPI_Process % doMPIAction ) then
               allocate ( bfOrder_local(nZones) )
               bfOrder_local = bfOrder
               call mpi_allreduce ( bfOrder_local, bfOrder, nZones, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr )
               deallocate ( bfOrder_local )
            end if

            do zoneID=1, nZones
#endif

!              Select the BC poynomial order
!              -----------------------------

               if ( self % ConformingOnZone(zoneID) .or. self % ignoreBCnonConformities) then
                  !bfOrder(zoneID) = bfOrder(zoneID)
               else
                  bfOrder(zoneID) = bfOrder(zoneID)/2
                  if ( bfOrder(zoneID) < 1 ) then
                     write(STD_OUT,*) 'ERROR :: The chosen polynomial orders are too low to represent the boundaries accurately'
                     write(STD_OUT,*) '      :: Nonconforming representations on boundaries need N>=2'
                     error stop
                  end if
               end if

               call NodalStorage (bfOrder(zoneID)) % construct (self % nodeType, bfOrder(zoneID))

            end do
         end if
!
!        **************************************************************
!        Check surfaces' integrity and adapt them to the solution order
!        **************************************************************
!
         do ii=1, num_of_faces
            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if
            associate( f => self % faces(fID) )

!
!           Check if the surfaces description in mesh file is consistent
!           ------------------------------------------------------------
            select case (f % faceType)

            case (HMESH_INTERIOR)
               if (f % IsMortar==0 .OR. f % IsMortar==2) then
                  eIDLeft  = f % elementIDs(1)
                  SideIDL  = f % elementSide(1)
                  !write(*,*)'line 3153 of HexMesh, eid:',eIDLeft
                  !if (self%elements(eIDLeft)%sliding) write(*,*)'its sliding element'
                  !if (self%elements(eIDLeft)%sliding_newnodes) write(*,*)'its sliding_newnodes element'
               !if (eIDLeft==0 ) cycle 
                  NSurfL   = SurfInfo(eIDLeft) % facePatches(SideIDL) % noOfKnots - 1
                  

                  eIDRight = f % elementIDs(2)
                  SideIDR  = f % elementSide(2)
                  !write(*,*)'line 3159 of HexMesh, eid:',eIDRight
                  !if (self%elements(eIDRight)%sliding) write(*,*)'its sliding element'
                  !if (self%elements(eIDRight)%sliding_newnodes) write(*,*)'its sliding_newnodes element'
                  !if (eIDRight==0 ) cycle 
                  NSurfR   = SurfInfo(eIDRight) % facePatches(SideIDR) % noOfKnots - 1

!              If both surfaces are of order 1.. There's no need to continue analyzing face
!              ----------------------------------------------------------------------------
                  if     ((SurfInfo(eIDLeft)  % IsHex8) .and. (SurfInfo(eIDRight) % IsHex8)) then
                     cycle
                  elseif ((SurfInfo(eIDLeft)  % IsHex8) .and. all(NSurfR == 1) ) then
                     cycle
                  elseif ((SurfInfo(eIDRight) % IsHex8) .and. all(NSurfL == 1) ) then
                     cycle
                  elseif (all(NSurfL == 1) .and. all(NSurfR == 1) ) then
                     cycle
                  elseif (any(NSurfL /= NSurfR)) then ! Only works for mesh files with isotropic boundary orders
                     write(STD_OUT,*) 'WARNING: Curved face definitions in mesh are not consistent.'
                     write(STD_OUT,*) '   Face:    ', fID
                     write(STD_OUT,*) '   Elements:', f % elementIDs
                     write(STD_OUT,*) '   N Left:  ', SurfInfo(eIDLeft)  % facePatches(SideIDL) % noOfKnots - 1
                     write(STD_OUT,*) '   N Right: ', SurfInfo(eIDRight) % facePatches(SideIDR) % noOfKnots - 1
                  end if

                  CLN(1) = min(f % NfLeft(1),f % NfRight(1))
                  CLN(2) = min(f % NfLeft(2),f % NfRight(2))
!
!              Adapt the curved face order to the polynomial order
!              ---------------------------------------------------
                  if ( any(CLN < NSurfL) ) then
                     allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                    !! if (size(SurfInfo(eIDLeft) % facePatches(SideIDL)%points(1,:,1))==2) then 
                     !write(*,*) 'size points(1,:,1)',size(SurfInfo(eIDLeft) % facePatches(SideIDL)%points(1,:,1))
                     !write(*,*) 'size points(1,1,:)',size(SurfInfo(eIDLeft) % facePatches(SideIDL)%points(1,1,:))
                     !   write(*,*)'here we are, eIDLeft=', eIDLeft 
                     !   write(*,*) 'isSliding?', self%elements(eIDLeft)%sliding
                    !    write(*,*) 'isSliding new nodes?', self%elements(eIDLeft)%sliding_newnodes
                     !   write(*,*) 'ID of the face', f%ID 
                     !   write(*,*) 'Mortar face', f%IsMortar
                     !end if 
                     call ProjectFaceToNewPoints(SurfInfo(eIDLeft) % facePatches(SideIDL), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                          CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                     call SurfInfo(eIDLeft) % facePatches(SideIDL) % Destruct()
                     call SurfInfo(eIDLeft) % facePatches(SideIDL) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                              NodalStorage(CLN(2)) % xCGL,faceCL)
                    ! write(*,*) 'line 3391 of hexmesh, NodalStorage(CLN(1)) % xCGL',NodalStorage(CLN(1)) % xCGL
                    ! write(*,*) 'line 3392 of hexmesh, NodalStorage(CLN(2)) % xCGL',NodalStorage(CLN(2)) % xCGL
                    ! write(*,*) 'line 3393 of hexmesh, faceCLL',faceCL

                     deallocate(faceCL)
                  end if

                  select case ( f % rotation )
                  case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
                     if (CLN(1) /= CLN(2)) then
                        buffer = CLN(1)
                        CLN(1) = CLN(2)
                        CLN(2) = buffer
                     end if
                  end select

                  if ( any(CLN < NSurfR) ) then       ! TODO JMT: I have added this.. is correct?      
                     allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                     call ProjectFaceToNewPoints(SurfInfo(eIDRight) % facePatches(SideIDR), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                          CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                     call SurfInfo(eIDRight) % facePatches(SideIDR) % Destruct()
                     call SurfInfo(eIDRight) % facePatches(SideIDR) % Construct(NodalStorage(CLN(1)) % xCGL,&
                                                                              NodalStorage(CLN(2)) % xCGL,faceCL)
                     deallocate(faceCL)
                  end if
               end if 
               if (f % IsMortar==1 .OR. f % IsMortar==3) then 
                 ! if (f%IsMortar==3) write(*,*)'line 3414 HexMesh mortar is 3'
                  eIDLeft  = f % elementIDs(1)
                  SideIDL  = f % elementSide(1)
                  NSurfL   = SurfInfo(eIDLeft)  % facePatches(SideIDL) % noOfKnots - 1
   
                  if     (SurfInfo(eIDLeft) % IsHex8 .or. all(NSurfL == 1)) cycle
   
                  if (self % anisotropic  .and. (.not. self % meshIs2D) ) then
                     CLN = bfOrder(f % zone)
                  else
                     CLN(1) = f % NfLeft(1)
                     CLN(2) = f % NfLeft(2)
                  end if
   !
   !              Adapt the curved face order to the polynomial order
   !              ---------------------------------------------------
                  if ( any(CLN < NSurfL) ) then
                     allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                     call ProjectFaceToNewPoints(SurfInfo(eIDLeft) % facePatches(SideIDL), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                           CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                     call SurfInfo(eIDLeft) % facePatches(SideIDL) % Destruct()
                     call SurfInfo(eIDLeft) % facePatches(SideIDL) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                               NodalStorage(CLN(2)) % xCGL,faceCL)
                     deallocate(faceCL)
                  end if
               end if 
            case (HMESH_BOUNDARY)
               eIDLeft  = f % elementIDs(1)
               SideIDL  = f % elementSide(1)
               NSurfL   = SurfInfo(eIDLeft)  % facePatches(SideIDL) % noOfKnots - 1

               if     (SurfInfo(eIDLeft) % IsHex8 .or. all(NSurfL == 1)) cycle

               if (self % anisotropic  .and. (.not. self % meshIs2D) ) then
                  CLN = bfOrder(f % zone)
               else
                  CLN(1) = f % NfLeft(1)
                  CLN(2) = f % NfLeft(2)
               end if
!
!              Adapt the curved face order to the polynomial order
!              ---------------------------------------------------
               if ( any(CLN < NSurfL) ) then
                  allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                  call ProjectFaceToNewPoints(SurfInfo(eIDLeft) % facePatches(SideIDL), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                        CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                  call SurfInfo(eIDLeft) % facePatches(SideIDL) % Destruct()
                  call SurfInfo(eIDLeft) % facePatches(SideIDL) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                            NodalStorage(CLN(2)) % xCGL,faceCL)
                  deallocate(faceCL)
               end if

            case (HMESH_MPI)
               eID = maxval(f % elementIDs)
               thisSide = maxloc(f % elementIDs, dim=1)
               side = f % elementSide(thisSide)
               NSurf = SurfInfo(eID) % facePatches(side) % noOfKnots - 1

               if ( SurfInfo(eID) % IsHex8 .or. all(NSurf == 1) ) cycle

               if (self % elements(eID) % faceSide(side) == LEFT) then
                  CLN(1) = f % NfLeft(1)  ! TODO in MPI faces, p-adaption has  
                  CLN(2) = f % NfLeft(2)  ! not been accounted yet.  
               else
                  CLN(1) = f % NfRight(1)  ! TODO in MPI faces, p-adaption has  
                  CLN(2) = f % NfRight(2)  ! not been accounted yet.  
               end if

               if ( side .eq. 2 ) then    ! Right faces need to be rotated
                  select case ( f % rotation )
                  case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular  ! TODO this is correct? 
                     if (CLN(1) /= CLN(2)) then
                        buffer = CLN(1)
                        CLN(1) = CLN(2)
                        CLN(2) = buffer
                     end if
                  end select
               end if

               if ( any(CLN < NSurf) ) then
                  allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                  call ProjectFaceToNewPoints(SurfInfo(eID) % facePatches(side), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                 CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                  call SurfInfo(eID) % facePatches(side) % Destruct()
                  call SurfInfo(eID) % facePatches(side) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                     NodalStorage(CLN(2)) % xCGL,faceCL)
                  deallocate(faceCL)
               end if
            end select
            end associate
         end do
         safedeallocate (bfOrder)

!
!        ----------------------------
!        Construct elements' geometry
!        ----------------------------
!
         allocate(hex8Map)
         call hex8Map % constructWithCorners(corners)
         allocate(genHexMap)

         do ii=1, num_of_elements
            if ( present(elementList) ) then
               eID = elementList(ii)
            else
               eID = ii
            end if

            if (SurfInfo(eID) % IsHex8) then
               !if (eID==170 ) then 
               !   write(*,*)'SurfInfo(eID) % corners of', eID, SurfInfo(eID) % corners
               !end if 
               !if (eID==167) then
               !   write(*,*)'SurfInfo(eID) % corners of', eID, SurfInfo(eID) % corners
               !end if 

               !if (self%elements(eID)%sliding_newnodes) then
                  !write(*,*)'sliding newnodes'
                  !write(*,*)'SurfInfo(eID) % corners of', eID, SurfInfo(eID) % corners
               !end if 
               call hex8Map % setCorners(SurfInfo(eID) % corners)
               hexMap => hex8Map
            else
               CALL genHexMap % destruct()
               CALL genHexMap % constructWithFaces(SurfInfo(eID) % facePatches)

               hexMap => genHexMap
            end if

            call self % elements(eID) % ConstructGeometry (hexMap)

         end do
!
!        -------------------------
!        Construct faces' geometry
!        -------------------------
!
         do ii=1, num_of_faces
            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if

            associate(f => self % faces(fID))
               !write(*,*) 'faceID= ', fID 
               !write(*,*) 'mortar type', f % IsMortar
            select case(f % faceType)
            case(HMESH_INTERIOR, HMESH_BOUNDARY)
               if (self%nonconforming) then
               select case (f % IsMortar)
               case (0, 1)
               associate(eL => self % elements(f % elementIDs(1)))
                  !write(*,*) f%elementIDs(1)
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         NodalStorage(f % Nf), NodalStorage(eL % Nxyz), &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0 )
               end associate
               case (2) 
               associate(eR => self % elements(f % elementIDs(2)))
                 ! write(*,*) f%elementIDs(2)
                 !write(*,*) f % NelRight
                 ! write(*,*)  f % NfRight
               call f % geom % construct(f % Nf, f % NelRight, f % NfRight, eR % Nxyz, &
                                         NodalStorage(f % Nf), NodalStorage(eR % Nxyz), &
                                         eR % geom, eR % hexMap, f % elementSide(2), &
                                         f % projectionType(2), 2, 0 )
               end associate
               end select 
            else 
               associate(eL => self % elements(f % elementIDs(1)))
                  !write(*,*) f%elementIDs(1)
                 ! write(*,*) 'fl',f % NfLeft
                  !write(*,*) 'constructing geometry if face:',fID
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         NodalStorage(f % Nf), NodalStorage(eL % Nxyz), &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0 )
               end associate
            end if
            case(HMESH_MPI)
               side = maxloc(f % elementIDs, dim=1)

               select case (side)
               case(1)
                  Nel = f % NelLeft
                  Nelf = f % NfLeft
                  rot = 0

               case(2)
                  Nel = f % NelRight
                  Nelf = f % NfRight
                  rot = f % rotation

               end select

               associate(e => self % elements(f % elementIDs(side)))
               call f % geom % construct(f % Nf, Nelf, Nel, e % Nxyz, &
                                         NodalStorage(f % Nf), NodalStorage(e % Nxyz), &
                                         e % geom, e % hexMap, f % elementSide(side), &
                                         f % projectionType(side), side, rot)

               end associate

            end select
            end associate

         end do
!
!        ------------------------------------------------------
!        Compute the faces minimum orthogonal distance estimate
!        ------------------------------------------------------
!
         do ii = 1, num_of_faces
            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if
            associate(f => self % faces(fID))
            select case(f % faceType)
            case(HMESH_INTERIOR)
               if (f % IsMortar ==0 .OR. f % IsMortar == 2) then 
               f % geom % h = min(minval(self % elements(f % elementIDs(1)) % geom % jacobian), &
                                  minval(self % elements(f % elementIDs(2)) % geom % jacobian)) &
                        / maxval(f % geom % jacobian)
               !if (f % geom % h==0.0_RP) then 
               !   write(*,*)'self % elements(f % elementIDs(1)) % geom % jacobian',self % elements(f % elementIDs(1)) % geom % jacobian
               !   write(*,*)'self % elements(f % elementIDs(2)) % geom % jacobian',self % elements(f % elementIDs(2)) % geom % jacobian
               !end if 
               elseif(f % IsMortar==1 .or. f % IsMortar==3) then 
                  f % geom % h = minval(self % elements(f % elementIDs(1)) % geom % jacobian) &
                  / maxval(f % geom % jacobian)
               end if 
            case(HMESH_BOUNDARY)
                              !write(*,*) 'jacface', f % geom % jacobian
               f % geom % h = minval(self % elements(f % elementIDs(1)) % geom % jacobian) &
                        / maxval(f % geom % jacobian)
            case(HMESH_MPI)
               f % geom % h = minval(self % elements(maxval(f % elementIDs)) % geom % jacobian) &
                        / maxval(f % geom % jacobian)
            end select
            end associate
         end do

         if ( MPI_Process % doMPIAction ) then
            call CommunicateMPIFaceMinimumDistance(self)
         end if

!
!        ---------
!        Finish up
!        ---------
!
         deallocate (SurfInfo)
         CALL hex8Map % destruct()
         DEALLOCATE(hex8Map)
         CALL genHexMap % destruct()
         DEALLOCATE(genHexMap)

      end subroutine HexMesh_ConstructGeometry

      subroutine CommunicateMPIFaceMinimumDistance(self)
         implicit none
         class(HexMesh) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
#ifdef _HAS_MPI_
         integer  :: no_of_max_faces, i, fID, mpifID, ierr
         real(kind=RP), allocatable  :: hsend(:,:)
         real(kind=RP), allocatable  :: hrecv(:,:)
         integer  :: sendReq(MPI_Process % nProcs)
         integer  :: recvReq(MPI_Process % nProcs)
         integer  :: sendSt(MPI_STATUS_SIZE, MPI_Process % nProcs)
         integer  :: recvSt(MPI_STATUS_SIZE, MPI_Process % nProcs)

         if ( .not. self % MPIfaces % Constructed ) return
!
!        Get the maximum number of faces
!        -------------------------------
         no_of_max_faces = 0
         do i = 1, MPI_Process % nProcs
            no_of_max_faces = max(no_of_max_faces, self % MPIfaces % faces(i) % no_of_faces)
         end do
!
!        Allocate an array to store the distances
!        ----------------------------------------
         allocate(hsend(no_of_max_faces, MPI_Process % nProcs))
         allocate(hrecv(no_of_max_faces, MPI_Process % nProcs))
!
!        Perform the receive calls
!        -------------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) then
               recvReq(i) = MPI_REQUEST_NULL
               cycle

            end if

            call mpi_irecv(hrecv(:,i), self % MPIfaces % faces(i) % no_of_faces, MPI_DOUBLE, i-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, recvReq(i), ierr)

         end do

!
!        Gather the distances to send
!        ------------------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(i) % no_of_faces
               fID = self % MPIfaces % faces(i) % faceIDs(mpifID)
               hsend(mpifID,i) = self % faces(fID) % geom % h
            end do
         end do
!
!        Send the distances
!        ------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) then
               sendReq(i) = MPI_REQUEST_NULL
               cycle

            end if

            call mpi_isend(hsend(:,i), self % MPIfaces % faces(i) % no_of_faces, MPI_DOUBLE, i-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, sendReq(i), ierr)

         end do

         call mpi_waitall(MPI_Process % nProcs, recvReq, recvSt, ierr)
!
!        Collect the distances
!        ---------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(i) % no_of_faces
               fID = self % MPIfaces % faces(i) % faceIDs(mpifID)
               self % faces(fID) % geom % h = min(self % faces(fID) % geom % h, hrecv(mpifID,i))
            end do
         end do

         call mpi_waitall(MPI_Process % nProcs, sendReq, sendSt, ierr)

         deallocate(hsend, hrecv)

#endif
      end subroutine CommunicateMPIFaceMinimumDistance

      subroutine HexMesh_DefineAsBoundaryFaces(self)
         implicit none
         class(HexMesh) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fID

         do fID = 1, size(self % faces)
            if ( self % faces(fID) % faceType .eq. HMESH_UNDEFINED ) then
               if ( trim(self % faces(fID) % boundaryName) .ne. emptyBCName ) then
                  self % faces(fID) % faceType = HMESH_BOUNDARY
               end if
            end if
         end do

      end subroutine HexMesh_DefineAsBoundaryFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------
!     Export mesh to a hmesh file
!     ---------------------------
      subroutine HexMesh_Export(self, fileName)
         use SolutionFile
         use MPI_Process_Info
         implicit none
         class(HexMesh),   intent(in)     :: self
         character(len=*), intent(in)     :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: fid, eID
         integer(kind=AddrInt)         :: pos
         character(len=LINE_LENGTH)    :: meshName
         real(kind=RP), parameter      :: refs(NO_OF_SAVED_REFS) = 0.0_RP

!
!        Create file: it will be contained in ./MESH
!        -------------------------------------------
         meshName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".hmesh"
         call CreateNewSolutionFile( trim(meshName), MESH_FILE, self % nodeType, &
                                     self % no_of_allElements, 0, 0.0_RP, refs)
!
!        Introduce all element nodal coordinates
!        ---------------------------------------
         fID = putSolutionFileInWriteDataMode(trim(meshName))
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 3_AddrInt*e % offsetIO*SIZEOF_RP
            call writeArray(fID, e % geom % x(:,0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3))*Lref, position=pos)
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(meshName))

      end subroutine HexMesh_Export
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ExportBoundaryMesh(self,fileName)
         implicit none
         !-arguments--------------------------------
         class(HexMesh),   intent(in)     :: self
         character(len=*), intent(in)     :: fileName          !<  Name of file containing polynomial orders to initialize
         !-local-variables--------------------------
         integer                          :: fd       ! File unit
         integer                          :: zoneID, zfID, fID
         character(len=LINE_LENGTH)       :: bMeshName
         !------------------------------------------


!
!        Create file: it will be contained in ./MESH
!        -------------------------------------------
         bMeshName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".bmesh"

         open( newunit = fd , file = trim(bMeshName), action = 'write')

!
!        Write file
!        ----------
         write(fd,*) size (self % zones)

         do zoneID = 1, size (self % zones)

            write (fd,*) self % zones(zoneID) % Name
            write (fd,*) self % zones(zoneID) % no_of_faces

            do zfID=1, self % zones(zoneID) % no_of_faces
               fID = self % zones(zoneID) % faces(zFID)
               write(fd,'(I8)',advance='no') self % faces(fID) % elementIDs(1)
            end do
            write(fd,*)

            do zfID=1, self % zones(zoneID) % no_of_faces
               fID = self % zones(zoneID) % faces(zFID)
               write(fd,'(I8)',advance='no') self % faces(fID) % elementSide(1)
            end do
            write(fd,*)

         end do

         close (fd)
      end subroutine HexMesh_ExportBoundaryMesh

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------
!     Export mesh orders to a file
!     ----------------------------
      subroutine HexMesh_ExportOrders(self,fileName)
         implicit none
         !------------------------------------------
         class(HexMesh),   intent(in)     :: self
         character(len=*), intent(in)     :: fileName          !<  Name of file containing polynomial orders to initialize
         !--local-variables-------------------------
         integer                          :: fd       ! File unit
         integer                          :: k
         character(len=LINE_LENGTH)       :: OrderFileName
         !--local-MPI-variables---------------------
         integer                          :: Nx(self % no_of_elements), total_Nx(self % no_of_allElements)
         integer                          :: Ny(self % no_of_elements), total_Ny(self % no_of_allElements)
         integer                          :: Nz(self % no_of_elements), total_Nz(self % no_of_allElements)
         integer                          :: globIDs(self % no_of_elements), all_globIDs(self % no_of_allElements), all_globIDs_copy(self % no_of_allElements)
         integer                          :: displ(MPI_Process % nProcs), no_of_elements_array(MPI_Process % nProcs)
         integer                          :: zID, eID, ierr

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
!
!        Share info with other processes
!        -------------------------------

         call mpi_allgather(self % no_of_elements, 1, MPI_INT, no_of_elements_array, 1, MPI_INT, MPI_COMM_WORLD, ierr)
!
!        Compute the displacements
!        -------------------------
         displ(1) = 0
         do zID = 1, MPI_Process % nProcs-1
            displ(zID+1) = displ(zID) + no_of_elements_array(zID)
         end do
!
!        Get global element IDs in all partitions
!        -----------------------------------------
         do eID=1, self % no_of_elements
            globIDs(eID) = self % elements(eID) % globID
         end do

         call mpi_allgatherv(globIDs, self % no_of_elements, MPI_INT, all_globIDs, no_of_elements_array , displ, MPI_INT, MPI_COMM_WORLD, ierr)
         all_globIDs_copy(:) = all_globIDs(:)
!
!        Get polynomial order in all partitions
!        ----------------------------------------
         do eID = 1, self % no_of_elements
            Nx(eID) = self % elements(eID) % Nxyz(1)
            Ny(eID) = self % elements(eID) % Nxyz(2)
            Nz(eID) = self % elements(eID) % Nxyz(3)
         enddo

         call mpi_allgatherv(Nx, self % no_of_elements, MPI_INT, total_Nx, no_of_elements_array, displ, MPI_INT, MPI_COMM_WORLD, ierr)

         call mpi_allgatherv(Ny, self % no_of_elements, MPI_INT, total_Ny, no_of_elements_array, displ, MPI_INT, MPI_COMM_WORLD, ierr)

         call mpi_allgatherv(Nz, self % no_of_elements, MPI_INT, total_Nz, no_of_elements_array, displ, MPI_INT, MPI_COMM_WORLD, ierr)

!
!        Reorganize polynomial order
!        ----------------------
         call QsortWithFriend(all_globIDs, total_Nx)
         all_globIDs(:) = all_globIDs_copy(:)

         call QsortWithFriend(all_globIDs, total_Ny)
         all_globIDs(:) = all_globIDs_copy(:)

         call QsortWithFriend(all_globIDs, total_Nz)
!
!        Create file: it will be contained in ./MESH
!        -------------------------------------------
         if ( MPI_Process % isRoot ) then
            OrderFileName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".omesh"
            open( newunit = fd , FILE = TRIM(OrderFileName), ACTION = 'write')

            write(fd,*) self % no_of_allElements

            do k=1, self % no_of_allElements
               write(fd,*) total_Nx(k), total_Ny(k), total_Nz(k)
            end do

            close (fd)
         end if
#endif
      else
         OrderFileName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".omesh"
         open( newunit = fd , FILE = TRIM(OrderFileName), ACTION = 'write')

         write(fd,*) self % no_of_elements

         do k=1, self % no_of_elements
            write(fd,*) self % elements(k) % Nxyz
         end do

         close (fd)

      end if

      end subroutine HexMesh_ExportOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!     ************************************************************************
!           Save solution subroutine for the Navier-Stokes solver. It saves
!        the state vector (Q), and optionally the gradients.
!     ************************************************************************
!
      subroutine HexMesh_SaveSolution(self, iter, time, name, saveGradients, saveSensor_, saveLES_)
         use SolutionFile
         use MPI_Process_Info
         implicit none
         class(HexMesh)                         :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
         logical,             intent(in)        :: saveGradients
         logical, optional,   intent(in)        :: saveSensor_
         logical, optional,   intent(in)        :: saveLES_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fid, eID, padding
         integer(kind=AddrInt)            :: pos
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS)
         real(kind=RP), allocatable       :: Q(:,:,:,:)
         logical                          :: saveSensor, saveLES
#if (!defined(NAVIERSTOKES))
         logical                          :: computeGradients = .true.
#endif
!
!        Gather reference quantities
!        ---------------------------
#if defined(NAVIERSTOKES)
         refs(GAMMA_REF) = thermodynamics % gamma
         refs(RGAS_REF)  = thermodynamics % R
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = refValues      % T
         refs(MACH_REF)  = dimensionless  % Mach
#elif defined(INCNS)
         refs(GAMMA_REF) = 0.0_RP
         refs(RGAS_REF)  = 0.0_RP
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = 0.0_RP
         refs(MACH_REF)  = 0.0_RP
#else
         refs = 0.0_RP
#endif
!
!        Create new file
!        ---------------
         if (present(saveSensor_)) then
            saveSensor = saveSensor_
         else
            saveSensor = .false.
         end if
         if (present(saveLES_)) then
            saveLES = saveLES_
         else
            saveLES = .false.
         end if

         if (saveGradients .and. computeGradients) then
            if (saveSensor) then
               call CreateNewSolutionFile(trim(name), SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE, &
                                          self % nodeType, self % no_of_allElements, iter, time, refs)
            else
               call CreateNewSolutionFile(trim(name), SOLUTION_AND_GRADIENTS_FILE, &
                                          self % nodeType, self % no_of_allElements, iter, time, refs)
            end if
            padding = NCONS + 3*NGRAD
         else
            if (saveSensor) then
               call CreateNewSolutionFile(trim(name), SOLUTION_AND_SENSOR_FILE, self % nodeType, &
                                          self % no_of_allElements, iter, time, refs)
            else
               call CreateNewSolutionFile(trim(name), SOLUTION_FILE, self % nodeType, &
                                          self % no_of_allElements, iter, time, refs)
            end if
            padding = NCONS
         end if

         if (saveLES) padding = padding + 2
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )

            allocate(Q(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
#ifdef FLOW
            Q(1:NCONS,:,:,:)  = e % storage % Q
#ifdef MULTIPHASE
            Q(IMP,:,:,:) = e % storage % Q(IMP,:,:,:) + e % storage % Q(IMC,:,:,:)*e % storage % mu(1,:,:,:)
#endif
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
            Q(NCONS,:,:,:) = e % storage % c(1,:,:,:)
#endif

            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*padding*e % offsetIO * SIZEOF_RP
            if (saveSensor) pos = pos + (e % globID - 1) * SIZEOF_RP
            call writeArray(fid, Q, position=pos)

            deallocate(Q)
            if ( saveGradients .and. computeGradients ) then

               allocate(Q(NGRAD,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))

#ifdef FLOW
               Q(1:NCONS,:,:,:) = e % storage % U_x
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(NGRAD,:,:,:) = e % storage % c_x(1,:,:,:)
#endif
               write(fid) Q

#ifdef FLOW
               Q(1:NCONS,:,:,:) = e % storage % U_y
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(NGRAD,:,:,:) = e % storage % c_y(1,:,:,:)
#endif
               write(fid) Q

#ifdef FLOW
               Q(1:NCONS,:,:,:) = e % storage % U_z
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(NGRAD,:,:,:) = e % storage % c_z(1,:,:,:)
#endif
               write(fid) Q

               deallocate(Q)
            end if

            if (saveSensor) then
               write(fid) e % storage % sensor
            end if

          if (saveLES) then
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
               allocate(Q(1,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               Q(1,:,:,:) = e % storage % mu_NS(1,:,:,:) ! total viscosity = mu + mu_sgs
               write(fid) Q
               Q(1,:,:,:) = e % storage % mu_turb_NS(:,:,:) !mu_sgs
               write(fid) Q
               deallocate(Q)
#endif
          end if 

            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveSolution

#if defined(NAVIERSTOKES)
      subroutine HexMesh_SaveStatistics(self, iter, time, name, saveGradients)
         use SolutionFile
         implicit none
         class(HexMesh),      intent(in)        :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
         logical,             intent(in)        :: saveGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fid, eID
         integer                          :: no_stat_s
         integer(kind=AddrInt)            :: pos
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
         real(kind=RP), allocatable       :: Q(:,:,:,:)
!
!        Gather reference quantities
!        ---------------------------
         refs(GAMMA_REF) = thermodynamics % gamma
         refs(RGAS_REF)  = thermodynamics % R
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = refValues      % T
         refs(MACH_REF)  = dimensionless  % Mach
         refs(RE_REF)    = dimensionless  % Re

!        Create new file
!        ---------------
         call CreateNewSolutionFile(trim(name),STATS_FILE, self % nodeType, self % no_of_allElements, iter, time, refs)
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_of_stats_variables*e % offsetIO*SIZEOF_RP
            no_stat_s = 9
            call writeArray(fid, e % storage % stats % data(1:no_stat_s,:,:,:), position=pos)
            allocate(Q(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            ! write(fid) e%storage%stats%data(7:,:,:,:)
            Q(1:NCONS,:,:,:) = e % storage % stats % data(no_stat_s+1:no_stat_s+NCONS,:,:,:)
            write(fid) Q
            deallocate(Q)
            if ( saveGradients .and. computeGradients ) then
               allocate(Q(NGRAD,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               ! UX
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1:no_stat_s+NCONS+NGRAD,:,:,:)
               write(fid) Q
               ! UY
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1+NGRAD:no_stat_s+NCONS+2*NGRAD,:,:,:)
               write(fid) Q
               ! UZ
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1+2*NGRAD:,:,:,:)
               write(fid) Q
               deallocate(Q)
            end if
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveStatistics

#endif

#if defined(INCNS) 
      subroutine HexMesh_SaveStatistics(self, iter, time, name, saveGradients)
         use SolutionFile
         implicit none
         class(HexMesh),      intent(in)        :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
         logical,             intent(in)        :: saveGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fid, eID
         integer                          :: no_stat_s
         integer(kind=AddrInt)            :: pos
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
         real(kind=RP), allocatable       :: Q(:,:,:,:)
!
!        Gather reference quantities (//# I COPIED THESE VERBATIM. Do I reaaly need these?)
!        ---------------------------

         refs(GAMMA_REF) = 0.0_RP
         refs(RGAS_REF)  = 0.0_RP
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = 0.0_RP
         refs(MACH_REF)  = 0.0_RP

!        Create new file
!        ---------------
         call CreateNewSolutionFile(trim(name),STATS_FILE, self % nodeType, self % no_of_allElements, iter, time, refs)
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + no_of_stats_variables*e % offsetIO*SIZEOF_RP
            no_stat_s = 9
            call writeArray(fid, e % storage % stats % data(1:no_stat_s,:,:,:), position=pos)
            allocate(Q(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
         !    ! write(fid) e%storage%stats%data(7:,:,:,:)
             Q(1:NCONS,:,:,:) = e % storage % stats % data(no_stat_s+1:no_stat_s+NCONS,:,:,:)
             write(fid) Q
             deallocate(Q)
            if ( saveGradients .and. computeGradients ) then
               allocate(Q(NGRAD,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               ! UX
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1:no_stat_s+NCONS+NGRAD,:,:,:)
               write(fid) Q
               ! UY
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1+NGRAD:no_stat_s+NCONS+2*NGRAD,:,:,:)
               write(fid) Q
               ! UZ
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1+2*NGRAD:,:,:,:)
               write(fid) Q
               deallocate(Q)
            end if
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveStatistics

#endif

#if defined(NAVIERSTOKES) || defined(INCNS)

      subroutine HexMesh_ResetStatistics(self)
         implicit none
         class(HexMesh)       :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID

         do eID = 1, self % no_of_elements
            self % elements(eID) % storage % stats % data = 0.0_RP
         end do

      end subroutine HexMesh_ResetStatistics
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------------
!     Subroutine to load a solution for restart using the information in the control file
!     -----------------------------------------------------------------------------------
      subroutine HexMesh_LoadSolutionForRestart( self, controlVariables, initial_iteration, initial_time, loadFromNSSA )
         use mainKeywordsModule, only: restartFileNameKey
         use FileReaders       , only: ReadOrderFile
         implicit none
         !-arguments-----------------------------------------------
         class(HexMesh)                       :: self
         type(FTValueDictionary), intent(in)  :: controlVariables
         logical                , intent(in)  :: loadFromNSSA
         integer                , intent(out) :: initial_iteration
         real(kind=RP)          , intent(out) :: initial_time
         !-local-variables-----------------------------------------
         character(len=LINE_LENGTH)           :: fileName
         character(len=LINE_LENGTH)           :: orderFileName
         integer, allocatable                 :: Nx(:), Ny(:), Nz(:)
         type(HexMesh), target                :: auxMesh
         integer                              :: NDOF, eID
         logical                              :: with_gradients
#if (!defined(NAVIERSTOKES)) || (!defined(INCNS)) || (!defined(MULTIPHASE))
         logical                          :: computeGradients = .true.
#endif
         !---------------------------------------------------------

         fileName = controlVariables % stringValueForKey(restartFileNameKey,requestedLength = LINE_LENGTH)

!
!        *****************************************************
!        The restart polynomial orders are different to self's
!        *****************************************************
!
         if ( controlVariables % containsKey("restart polorder")        .or. &
              controlVariables % containsKey("restart polorder file") ) then


!           Read the polynomial order of the solution to be loaded
!           ------------------------------------------------------

            if ( controlVariables % containsKey("restart polorder") ) then
               allocate ( Nx(self % no_of_allElements) , Ny(self % no_of_allElements) , Nz(self % no_of_allElements) )
               Nx = controlVariables % integerValueForKey ("restart polorder")
               Ny = Nx
               Nz = Nx
            elseif ( controlVariables % containsKey("restart polorder file") ) then
               orderFileName = controlVariables % stringValueForKey("restart polorder file",requestedLength = LINE_LENGTH)
               call ReadOrderFile(trim(orderFileName), Nx, Ny, Nz)
            end if

            if ( controlVariables % containsKey("restart nodetype") ) then
               select case ( trim(controlVariables % stringValueForKey("restart nodetype",requestedLength = LINE_LENGTH)) )
               case("Gauss")
                  auxMesh % nodeType = 1
               case("Gauss-Lobatto")
                  auxMesh % nodeType = 2
               end select
            else
               auxMesh % nodeType = self % nodeType
            end if
!
!           Construct an auxiliary mesh to read the solution
!           -----------------------------------------------

            auxMesh % no_of_elements = self % no_of_elements
            auxMesh % no_of_allElements = self % no_of_allElements
            allocate ( auxMesh % elements (self % no_of_elements) )
            allocate ( auxMesh % Nx (self % no_of_elements) )
            allocate ( auxMesh % Ny (self % no_of_elements) )
            allocate ( auxMesh % Nz (self % no_of_elements) )

            NDOF = 0
            do eID = 1, self % no_of_elements
               associate ( e_aux => auxMesh % elements(eID), &
                           e     =>    self % elements(eID) )
               auxMesh % Nx(eID) = Nx (e % globID)
               auxMesh % Ny(eID) = Ny (e % globID)
               auxMesh % Nz(eID) = Nz (e % globID)
               e_aux % globID = e % globID
               e_aux % Nxyz = [Nx(e % globID) , Ny(e % globID) , Nz(e % globID)]
               NDOF = NDOF + (Nx(e % globID) + 1) * (Ny(e % globID) + 1) * (Nz(e % globID) + 1)               ! TODO: change for new NDOF             
               end associate
            end do

            call auxMesh % PrepareForIO
            call auxMesh % AllocateStorage (NDOF, controlVariables,computeGradients,.FALSE.)
            call auxMesh % storage % pointStorage
            do eID = 1, auxMesh % no_of_elements
               auxMesh % elements(eID) % storage => auxMesh % storage % elements(eID)
            end do

!           Read the solution in the auxiliary mesh and interpolate to current mesh
!           ----------------------------------------------------------------------

            call auxMesh % LoadSolution ( fileName, initial_iteration, initial_time , with_gradients, loadFromNSSA=loadFromNSSA )
            do eID=1, self % no_of_elements
               call auxMesh % storage % elements (eID) % InterpolateSolution (self % storage % elements(eID), auxMesh % nodeType , with_gradients)
            end do

!           Clean up
!           --------

            do eID = 1, auxMesh % no_of_elements
               call auxMesh % elements(eID) % storage % destruct
            end do
            call auxMesh % storage % destruct
            deallocate (auxMesh % elements)

            if ( controlVariables % containsKey("get discretization error of") ) call GetDiscretizationError(self,controlVariables)
!
!        *****************************************************
!        The restart polynomial orders are the same as in self
!        *****************************************************
!
         else
            call self % LoadSolution ( fileName, initial_iteration, initial_time, loadFromNSSA=loadFromNSSA )
         end if

      end subroutine HexMesh_LoadSolutionForRestart
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------
!     Subroutine to load a solution (*.hsol) in self
!     ----------------------------------------------
      subroutine HexMesh_LoadSolution( self, fileName, initial_iteration, initial_time , with_gradients, loadFromNSSA)
         IMPLICIT NONE
         CLASS(HexMesh)                  :: self
         character(len=*)                :: fileName
         integer           , intent(out) :: initial_iteration
         real(kind=RP)     , intent(out) :: initial_time
         logical, optional , intent(out) :: with_gradients
         logical, optional , intent(in)  :: loadFromNSSA

!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                        :: fID, eID, fileType, no_of_elements, flag, nodetype
         integer                        :: padding, pos
         integer                        :: Nxp1, Nyp1, Nzp1, no_of_eqs, array_rank, expectedNoEqs
         real(kind=RP), allocatable     :: Q(:,:,:,:)
         character(len=SOLFILE_STR_LEN) :: rstName
         logical                        :: gradients
         logical                        :: NS_from_NSSA
         logical                        :: has_sensor

         gradients = .FALSE.
         has_sensor = .FALSE.
         if (present(loadFromNSSA)) then
             NS_from_NSSA = loadFromNSSA
         else
             NS_from_NSSA = .FALSE.
         end if 
         expectedNoEqs = NCONS
!
!        Get the file title
!        ------------------
         rstName = getSolutionFileName(trim(fileName))
!
!        Get the file type
!        -----------------
         fileType = getSolutionFileType(trim(fileName))

         select case (fileType)
         case(MESH_FILE)
            print*, "The selected restart file is a mesh file"
            errorMessage(STD_OUT)
            error stop

         case(SOLUTION_FILE)
            padding = 1*NCONS

         case(SOLUTION_AND_SENSOR_FILE)
            padding = 1*NCONS
            has_sensor = .TRUE.

         case(SOLUTION_AND_GRADIENTS_FILE)
            padding = NCONS + 3 * NGRAD
            gradients = .TRUE.

         case(SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
            padding = NCONS + 3 * NGRAD
            gradients = .TRUE.
            has_sensor = .TRUE.

         case(STATS_FILE)
            print*, "The selected restart file is a statistics file"
            errorMessage(STD_OUT)
            error stop
         case default
            print*, "Unknown restart file format"
            errorMessage(STD_OUT)
            error stop
         end select
         if (NS_from_NSSA) then
             expectedNoEqs = NCONS + 1
             if (gradients) then
                 ! add 1 as NNSA has one more NCONS, and 3 as has one more NGRAD, the whole will be padding = (NCONS + 1) + 3*(NGRAD+1)
                 padding = padding + 1 + 3
             else
                 padding = padding + 1
             end if
         end if
!
!        Get the node type
!        -----------------
         nodeType = getSolutionFileNodeType(trim(fileName))

         if ( nodeType .ne. self % nodeType ) then
            print*, "WARNING: Solution file uses a different discretization nodes than the mesh."
            print*, "Add restart polorder = (Pol order in your restart file) in the control file if you want interpolation routines to be used."
            print*, "If restart polorder is not specified the values in the original set of nodes are loaded into the new nodes without interpolation."
            errorMessage(STD_OUT)
         end if
!
!        Read the number of elements
!        ---------------------------
         no_of_elements = getSolutionFileNoOfElements(trim(fileName))

         if ( no_of_elements .ne. self % no_of_allElements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the restart file ", &
                                   "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
         end if
!
!        Read the initial iteration and time
!        -----------------------------------
         call getSolutionFileTimeAndITeration(trim(fileName), initial_iteration, initial_time)
!
!        Read the terminator indicator
!        -----------------------------
         flag = getSolutionFileDataInitFlag(trim(fileName))

         if ( flag .ne. BEGINNING_DATA ) then
            print*, "Beginning data flag was not found in the file."
            errorMessage(STD_OUT)
            error stop
         end if
!
!        Read elements data
!        ------------------
         fID = putSolutionFileInReadDataMode(trim(fileName))
         do eID = 1, size(self % elements)
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + 1_AddrInt*padding*e % offsetIO*SIZEOF_RP
            if (has_sensor) pos = pos + (e % globID - 1) * SIZEOF_RP
            read(fID, pos=pos) array_rank
            read(fID) no_of_eqs, Nxp1, Nyp1, Nzp1
            if (      ((Nxp1-1) .ne. e % Nxyz(1)) &
                 .or. ((Nyp1-1) .ne. e % Nxyz(2)) &
                 .or. ((Nzp1-1) .ne. e % Nxyz(3)) &
                 .or. (no_of_eqs .ne. expectedNoEqs )       ) then
               write(STD_OUT,'(A,I0,A)') "Error reading restart file: wrong dimension for element "&
                                           ,eID,"."

               write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Element dimensions: ", e % Nxyz(1), &
                                                                     " ,", e % Nxyz(2), &
                                                                     " ,", e % Nxyz(3), &
                                                                     "."

               write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Restart dimensions: ", Nxp1-1, &
                                                                     " ,", Nyp1-1, &
                                                                     " ,", Nzp1-1, &
                                                                     "."

               errorMessage(STD_OUT)
               error stop
            end if

            allocate(Q(expectedNoEqs, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            read(fID) Q
#ifdef FLOW
            e % storage % QNS = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
            e % storage % c(1,:,:,:) = Q(NCONS,:,:,:)
#endif

            deallocate(Q)

            if (gradients) then
               allocate(Q(NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               read(fID) Q
#ifdef FLOW
               e % storage % U_x = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               e % storage % c_x(1,:,:,:) = Q(NGRAD,:,:,:)
#endif

               read(fID) Q
#ifdef FLOW
               e % storage % U_y = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               e % storage % c_y(1,:,:,:) = Q(NGRAD,:,:,:)
#endif
               read(fID) Q
#ifdef FLOW
               e % storage % U_z = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               e % storage % c_z(1,:,:,:) = Q(NGRAD,:,:,:)
#endif
               deallocate(Q)
            end if

            if (has_sensor) then
               read(fID) e % storage % sensor
            end if

           end associate
         end do
!
!        Close the file
!        --------------
         close(fID)

         if (present(with_gradients) ) with_gradients = gradients

      END SUBROUTINE HexMesh_LoadSolution


#if defined(ACOUSTIC)
      Subroutine HexMesh_SetUniformBaseFlow(self,Q_in)
          Implicit None
           CLASS(HexMesh)                  :: self
           real(kind=RP), dimension(1:NCONS), intent(in)  :: Q_in
  !
  !        ---------------
  !        Local variables
  !        ---------------
           INTEGER                        :: eID, eq
  
           do eID = 1, size(self % elements)
              do eq = 1,NCONS
                  self % elements(eID) % storage % Qbase(eq,:,:,:) = Q_in(eq)
              end do
           end do
  !
      End Subroutine HexMesh_SetUniformBaseFlow
  !
  !////////////////////////////////////////////////////////////////////////
  !
        subroutine HexMesh_ProlongBaseSolutionToFaces(self, nEqn)
           implicit none
           class(HexMesh),    intent(inout) :: self
           integer,           intent(in)    :: nEqn
  !
  !        ---------------
  !        Local variables
  !        ---------------
  !
           integer  :: fIDs(6)
           integer  :: eID, i
  
  !$omp do schedule(runtime) private(fIDs)
           do eID = 1, size(self % elements)
              fIDs = self % elements(eID) % faceIDs
              call self % elements(eID) % ProlongBaseSolutionToFaces(nEqn, &
                                                                     self % faces(fIDs(1)),&
                                                                     self % faces(fIDs(2)),&
                                                                     self % faces(fIDs(3)),&
                                                                     self % faces(fIDs(4)),&
                                                                     self % faces(fIDs(5)),&
                                                                     self % faces(fIDs(6)) )
           end do
  !$omp end do
  
        end subroutine HexMesh_ProlongBaseSolutionToFaces
#endif
  !
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetDiscretizationError(mesh,controlVariables)
      implicit none
      !----------------------------------------------
      type(HexMesh) :: mesh
      type(FTValueDictionary) :: controlVariables
      !----------------------------------------------
      character(len=LINE_LENGTH) :: fileName
      type(HexMesh)              :: refMesh
      integer       :: NDOF, eID
      integer       :: initial_iteration
      real(kind=RP) :: initial_time
      !----------------------------------------------

      fileName = controlVariables % stringValueForKey("get discretization error of",LINE_LENGTH)

!     Construct an auxiliary mesh to read the solution
!     -----------------------------------------------

      refMesh % nodeType = mesh % nodeType
      refMesh % no_of_elements = mesh % no_of_elements
      allocate ( refMesh % elements (mesh % no_of_elements) )

      NDOF = 0
      do eID = 1, mesh % no_of_elements
         associate ( e_aux => refMesh % elements(eID), &
                     e     =>    mesh % elements(eID) )
         e_aux % globID = e % globID
         e_aux % Nxyz = e % Nxyz
         NDOF = NDOF + product(e % Nxyz + 1)
         end associate
      end do

      call refMesh % PrepareForIO
      call refMesh % AllocateStorage (NDOF, controlVariables,.FALSE.)

!     Read the solution in the auxiliary mesh and interpolate to current mesh
!     ----------------------------------------------------------------------

      call refMesh % LoadSolution ( fileName, initial_iteration, initial_time )

      refMesh % storage % Q = mesh % storage % Q - refMesh % storage % Q

      print*, '|disc_error|∞ = ', maxval(abs(refMesh % storage % Q))
      print*, '|disc_error|² = ', norm2(refMesh % storage % Q)

      call refMesh % SaveSolution(0, 0._RP, 'RESULTS/DiscError.hsol', .FALSE.)

!           Clean up
!           --------

      do eID = 1, refMesh % no_of_elements
         call refMesh % elements(eID) % storage % destruct
      end do
      call refMesh % storage % destruct
      deallocate (refMesh % elements)

   end subroutine GetDiscretizationError
!
!////////////////////////////////////////////////////////////////////////
!
!        AUXILIARY SUBROUTINES
!        --------------------
!
!////////////////////////////////////////////////////////////////////////
!
      logical function HexMesh_FindPointWithCoords(self, x, eID, xi, optionalElements)
         implicit none
         class(HexMesh), intent(in)         :: self
         real(kind=RP),    intent(in)       :: x(NDIM)
         integer,          intent(out)      :: eID
         real(kind=RP),    intent(out)      :: xi(NDIM)
         integer, optional,intent(in)       :: optionalElements(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: op_eID
         integer     :: zoneID, fID
         logical     :: success

         HexMesh_FindPointWithCoords = .false.
!
!        Search in optionalElements (if present)
!        ---------------------------------------
         if ( present(optionalElements) ) then
            do op_eID = 1, size(optionalElements)
               if ( optionalElements(op_eID) .eq. -1 ) cycle
               associate(e => self % elements(optionalElements(op_eID)))
               success = e % FindPointWithCoords(x, self % dir2D_ctrl, xi)
               if ( success ) then
                  eID = optionalElements(op_eID)
                  HexMesh_FindPointWithCoords = .true.
                  return
               end if
               end associate
            end do
         end if
!
!        Search in linear (not curved) mesh (faster and safer)
!        -----------------------------------------------------
         do eID = 1, self % no_of_elements
            success = self % elements(eID) % FindPointInLinElement(x, self % nodes)
            if ( success ) exit
         end do
!
!        If found in linear mesh, use FindPointWithCoords in that element and, if necessary, in neighbors...
!        ---------------------------------------------------------------------------------------------------
         if (eID <= self % no_of_elements) then
            success = self % FindPointWithCoordsInNeighbors(x, xi, eID, 2)
            if ( success ) then
               HexMesh_FindPointWithCoords = .true.
               return
            end if
         end if
!
!        As a last resource, search using FindPointWithCoords only in boundary elements
!        ------------------------------------------------------------------------------
         do zoneID=1, size(self % zones)
            do fID=1, self % zones(zoneID) % no_of_faces

               op_eID = self % faces ( self % zones(zoneID) % faces(fID) ) % elementIDs(1)
               success = self % elements (op_eID) % FindPointWithCoords(x, self % dir2D_ctrl, xi)
               if ( success ) then
                  HexMesh_FindPointWithCoords = .true.
                  return
               end if
            end do
         end do

      end function HexMesh_FindPointWithCoords
!
!////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------------------------------------------------------------
!     HexMesh_FindPointWithCoordsInNeighbors:
!     This subroutine looks for a point (defined by coordinates) in the neighbor elements of a specific element
!     Note: For MPI, this routine ONLY checks in neighbors that are in the same partition...
!     ---------------------------------------------------------------------------------------------------------
      logical recursive function HexMesh_FindPointWithCoordsInNeighbors(self, x, xi, eID, depth)
         implicit none
         !-arguments--------------------------------------------------
         class(HexMesh), intent(in)     :: self
         real(kind=RP) , intent(in)     :: x(NDIM)
         real(kind=RP) , intent(out)    :: xi(NDIM)
         integer       , intent(inout)  :: eID
         integer       , intent(in)     :: depth
         !-local-variables--------------------------------------------
         logical :: success
         integer :: fID, nID, new_eID
         !------------------------------------------------------------

         success = self % elements(eID) % FindPointWithCoords(x, self % dir2D_ctrl, xi)
         if ( success ) then
            HexMesh_FindPointWithCoordsInNeighbors = .TRUE.
            return
         end if

         if (depth > 1) then
            do fID=1, FACES_PER_ELEMENT

               new_eID = mpi_partition % global2localeID (self % elements(eID) % Connection(fID) % globID)
               if (new_eID == 0) cycle
               success = self % FindPointWithCoordsInNeighbors(x, xi, new_eID, depth-1)
               if ( success ) then
                  HexMesh_FindPointWithCoordsInNeighbors = .TRUE.
                  return
               end if

            end do
         end if

      end function HexMesh_FindPointWithCoordsInNeighbors
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ComputeWallDistances(self,facesList,elementList)
         implicit none
         class(HexMesh)     :: self
         integer, optional , intent(in)    :: facesList(:)
         integer, optional , intent(in)    :: elementList(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: eID, ii, i, j, k, no_of_wallDOFS, num_of_elems, num_of_faces
         real(kind=RP) :: currentDistance, minimumDistance
         integer       :: fID
         real(kind=RP) :: xP(NDIM)
         real(kind=RP), allocatable    :: Xwall(:,:)
!
!        Gather all walls coordinates
!        ----------------------------
         call HexMesh_GatherAllWallCoordinates(self, no_of_wallDOFS, Xwall)
!
!        Get the minimum distance to each elements nodal degree of freedom
!        -----------------------------------------------------------------
         if ( present(elementList) ) then
            num_of_elems = size (elementList)
         else
            num_of_elems = self % no_of_elements
         end if
         do ii = 1, num_of_elems
            if ( present(elementList) ) then
               eID = elementList (ii)
            else
               eID = ii
            end if
            associate(e => self % elements(eID))
            allocate(e % geom % dWall(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            if( self% IBM% active ) then
               allocate(e % geom % normal(NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               e % geom % dWall = huge(1.0_RP)
            endif

            if( .not. self% IBM% active ) then
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  xP = e % geom % x(:,i,j,k)

                  minimumDistance = HUGE(1.0_RP)
                  do fID = 1, no_of_wallDOFS
                     currentDistance = sum(POW2(xP - Xwall(:,fID)))
                     minimumDistance = min(minimumDistance, currentDistance)
                  end do

                  e % geom % dWall(i,j,k) = sqrt(minimumDistance)

               end do                  ; end do                ; end do
            end if
            
            end associate
         end do
!
!        Get the minimum distance to each face nodal degree of freedom
!        -------------------------------------------------------------
         if ( present(facesList) ) then
            num_of_faces = size (facesList)
         else
            num_of_faces = size (self % faces)
         end if
         do ii = 1, num_of_faces
            if ( present(facesList) ) then
               eID = facesList (ii)
            else
               eID = ii
            end if

            associate(fe => self % faces(eID))
            allocate(fe % geom % dWall(0:fe % Nf(1), 0:fe % Nf(2)))
            if( self% IBM% active ) then
               fe % geom % dWall = huge(1.0_RP)
               if (fe%IsMortar==3) then 
                  if(.not. allocated(self%mortar_faces(fe%mortar(1))%geom%dwall)) then 
                     allocate(self%mortar_faces(fe%mortar(1))% geom % dWall(0:fe % Nf(1), 0:fe % Nf(2)))
                     self%mortar_faces(fe%mortar(1))%geom%dwall=fe % geom % dWall
                  end if 
                  if (size(fe%mortar)==2) then 
                  if(.not. allocated(self%mortar_faces(fe%mortar(2))%geom%dwall)) then 
                     allocate(self%mortar_faces(fe%mortar(2))% geom % dWall(0:fe % Nf(1), 0:fe % Nf(2)))
                     self%mortar_faces(fe%mortar(2))%geom%dwall=fe % geom % dWall
                  end if 
                  end if 
               end if 
            endif
            
            if( .not. self% IBM% active ) then
               do j = 0, fe % Nf(2) ; do i = 0, fe % Nf(1)
                  
                  xP = fe % geom % x(:,i,j)

                  minimumDistance = HUGE(1.0_RP)
                  do fID = 1, no_of_wallDOFS
                     currentDistance = sum(POW2(xP - Xwall(:,fID)))
                     minimumDistance = min(minimumDistance, currentDistance)
                  end do

                  fe % geom % dWall(i,j) = sqrt(minimumDistance)

                end do                ; end do
                if (fe%IsMortar==3) then 
                  if(.not. allocated(self%mortar_faces(fe%mortar(1))%geom%dwall)) then 
                     allocate(self%mortar_faces(fe%mortar(1))%geom%dWall(0:fe % Nf(1), 0:fe % Nf(2)))
                     self%mortar_faces(fe%mortar(1))%geom%dwall=fe % geom % dWall
                  end if 
                  if (size(fe%mortar)==2) then 
                  if(.not. allocated(self%mortar_faces(fe%mortar(2))%geom%dwall)) then 
                     allocate(self%mortar_faces(fe%mortar(2))%geom%dWall(0:fe % Nf(1), 0:fe % Nf(2)))
                     self%mortar_faces(fe%mortar(2))%geom%dwall=fe % geom % dWall
                  end if 
               end if 
               end if 
            end if
            
            end associate
         end do

         deallocate(Xwall)

      end subroutine HexMesh_ComputeWallDistances

      subroutine HexMesh_GatherAllWallCoordinates(self, no_of_wallDOFS, Xwall)
         implicit none
         class(HexMesh),              intent(in)  :: self
         integer,                     intent(out) :: no_of_wallDOFS
         real(kind=RP),  allocatable, intent(out) :: Xwall(:,:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: no_of_localWallDOFS
         integer                    :: no_of_wallDOFS_perProcess(MPI_Process % nProcs)
         integer                    :: displ(MPI_Process % nProcs)
         real(kind=RP), allocatable :: localXwall(:,:)
         integer                    :: zID, ierr, fID, zonefID, i, j, displacement

         no_of_localWallDOFS = 0
         do zID = 1, size(self % zones)
            if ( (trim(BCs(zID) % bc % bcType) .ne. "noslipwall") ) then
               cycle
            end if

            do zonefID = 1, self % zones(zID) % no_of_faces
               fID = self % zones(zID) % faces(zonefID)
               no_of_localWallDOFS = no_of_localWallDOFS + product(self % faces(fID) % Nf + 1)
            end do
         end do
         allocate( localXwall(NDIM, no_of_localWallDOFS) )
!
!        Loop all faces to gather the local wall coordinates
!        ---------------------------------------------------
         no_of_localWallDOFS = 0
         do zID = 1, size(self % zones)
            if ( (trim(BCs(zID) % bc % bcType) .ne. "noslipwall") ) then
               cycle
            end if

            do zonefID = 1, self % zones(zID) % no_of_faces
               fID = self % zones(zID) % faces(zonefID)
               associate( f => self % faces(fID) )
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  no_of_localWalLDOFS = no_of_localWallDOFS + 1
                  localXWall(1:NDIM,no_of_localWallDOFS) = f % geom % x(1:NDIM,i,j)
               end do               ; end do
               end associate
            end do
         end do
!
!        ******************************************************************
!        MPI: To communicate each process walls, we use the allgatherv
!             function
!        ******************************************************************
!
         if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
!
!           Perform a reduction to how many walls are in each process
!           ---------------------------------------------------------
            call mpi_allgather(no_of_localWallDOFS, 1, MPI_INT, no_of_wallDOFS_perProcess, 1, MPI_INT, MPI_COMM_WORLD, ierr)
!
!           Compute the displacements
!           -------------------------
            displ(1) = 0
            do zID = 1, MPI_Process % nProcs - 1
               displ(zID+1) = displ(zID) + no_of_wallDOFS_perProcess(zID) * NDIM
            end do
!
!           Allocate the data
!           -----------------
            no_of_wallDOFS = sum(no_of_wallDOFS_perProcess)
            allocate( Xwall(NDIM, no_of_wallDOFS) )
!
!           Perform an allgatherv
!           ---------------------
            call mpi_allgatherv(localXwall, NDIM*no_of_localWallDOFS, MPI_DOUBLE, Xwall, NDIM*no_of_wallDOFS_perProcess, displ, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
#endif
         else
            no_of_wallDOFS = no_of_localWallDOFS
            allocate( Xwall(NDIM, no_of_wallDOFS) )
            Xwall = localXwall
         end if

         deallocate(localXwall)

      end subroutine HexMesh_GatherAllWallCoordinates
!
!////////////////////////////////////////////////////////////////////////
!
!        CONSTRUCT ZONES
!        ---------------
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ConstructZones( self )
      implicit none
      class(HexMesh)          :: self

      call ConstructZones ( self % faces , self % zones )

      end subroutine HexMesh_ConstructZones

!
!///////////////////////////////////////////////////////////////////////
!
!  HexMesh_AllocateStorage:
!  Allocates the storage for the simulation
!  -> Storage specific to the analytical Jacobian is constructed by the corresponding class (AnalyticalJacobian.f90)
!
   subroutine HexMesh_AllocateStorage(self,NDOF,controlVariables,computeGradients,Face_Storage)
      implicit none
      !-----------------------------------------------------------
      class(HexMesh), target                 :: self
      integer                 , intent(in)   :: NDOF
      class(FTValueDictionary), intent(in)   :: controlVariables
      logical                 , intent(in)   :: computeGradients
      logical, optional       , intent(in)   :: Face_Storage
      !-----------------------------------------------------------
      integer :: bdf_order, eID, fID, RKSteps_num
      logical :: Face_St, FaceComputeQdot
      character(len=LINE_LENGTH) :: time_int
      character(len=LINE_LENGTH) :: mg_smoother
      real(kind=RP), dimension(:,:,:),     pointer :: Qff
      !-----------------------------------------------------------

      if ( present(Face_Storage) ) then
         Face_St = Face_Storage
      else
         Face_St = .TRUE.
      end if
      FaceComputeQdot = controlVariables % containsKey("acoustic analogy")

      time_int = controlVariables % stringValueForKey("time integration",LINE_LENGTH)
      call toLower (time_int)

      if     ( controlVariables % containsKey("bdf order")) then
         bdf_order = controlVariables % integerValueForKey("bdf order")
         RKSteps_num = 0
      elseif ( trim(time_int) == "fas" ) then
         bdf_order = -1
         RKSteps_num = 0
        if ( controlVariables % containsKey("mg smoother")) then
          mg_smoother = controlVariables % stringValueForKey("mg smoother",LINE_LENGTH)
          call toLower (mg_smoother)
          if ( (trim(mg_smoother) .eq. "irk") .or. (trim(mg_smoother) .eq. "birk5") &
            .or. (trim(mg_smoother) .eq. "ilu") .or. (trim(mg_smoother) .eq. "sgs") ) then
            bdf_order = 1
            RKSteps_num = 0
          end if
        end if
      elseif ( trim(time_int) == "imex" ) then
         bdf_order = 1
#ifdef MULTIPHASE
         RKSteps_num = 3
!         RKSteps_num = 0
#else
         RKSteps_num = 0
#endif
      elseif ( trim(time_int) == "rosenbrock" ) then
         bdf_order = 0
         RKSteps_num = 0
      else
         bdf_order = -1
         RKSteps_num = 0
      end if

#ifdef MULTIPHASE
      ! This is a fix to prevent a seg fault in debug mode
      ! implemented by g.rubio@upm.es 09/2023
      if ( trim(time_int) == "explicit" ) then
         bdf_order = 1  
         RKSteps_num = 0   
      endif  
#endif 
!     Construct global and elements' storage
!     --------------------------------------
      call self % storage % construct (NDOF, self % Nx, self % Ny, self % Nz, computeGradients, .FALSE., bdf_order, RKSteps_num )

!     Construct faces' storage
!     ------------------------
      if (Face_St) then
         do fID = 1, size(self % faces)
            associate ( f => self % faces(fID) )
               if (f % IsMortar==2) then 
                  call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, .FALSE., FaceComputeQdot, Mortar=.TRUE.)
                  call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, .FALSE., FaceComputeQdot, Mortar=.TRUE.)
               else 
               call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, .FALSE., FaceComputeQdot)
               call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, .FALSE., FaceComputeQdot)
            end if 

            end associate
         end do
      end if

      if (self%sliding) then 
         do fID = 1, size(self % mortar_faces)
            associate ( f => self % mortar_faces(fID) )
                  call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, .FALSE., FaceComputeQdot, Mortar=.TRUE.)
                  call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, .FALSE., FaceComputeQdot, Mortar=.TRUE.)
            end associate
         end do 
      end if 

!     Point element storage
!     ---------------------
      DO eID = 1, SIZE(self % elements)
         associate (e => self % elements(eID))
         e % hn = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)  ! Also compute h/p here
         e % storage => self % storage % elements(eID)
         end associate
      END DO

   end subroutine HexMesh_AllocateStorage

   subroutine HexMesh_SetStorageToEqn(self, which)
      implicit none
      class(HexMesh), target :: self
      integer, intent(in)    :: which
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: off, ns, c, mu, nssa, caa
      integer  :: eID, fID

      call GetStorageEquations(off, ns, c, mu, nssa, caa)


      if ( which .eq. ns .or. which .eq. nssa .or. which .eq. caa) then
#ifdef FLOW
         self % storage % Q => self % storage % QNS
         self % storage % QDot => self % storage % QDotNS
         self % storage % PrevQ(1:,1:) => self % storage % PrevQNS(1:,1:)

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToNS
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToNS
            call self % faces(fID) % storage(2) % SetStorageToNS
         end do

#endif

      elseif ( which .eq. c ) then
#if defined(CAHNHILLIARD)
         self % storage % Q => self % storage % c
         self % storage % QDot => self % storage % cDot
         self % storage % PrevQ => self % storage % PrevC

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToCH_c
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToCH_c
            call self % faces(fID) % storage(2) % SetStorageToCH_c
         end do
#endif
      elseif ( which .eq. mu ) then
#if defined(CAHNHILLIARD)
         self % storage % Q => self % storage % c
         self % storage % QDot => self % storage % cDot
         self % storage % PrevQ => self % storage % PrevC

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToCH_mu
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToCH_mu
            call self % faces(fID) % storage(2) % SetStorageToCH_mu
         end do
#endif
      end if

   end subroutine HexMesh_SetStorageToEqn
!
!///////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------------
!  Checks if the representation is conforming on a zone (Boundary)
!  ---------------------------------------------------------------
   function HexMesh_ConformingOnZone (self,zoneID) result(conforming)
      implicit none
      !-----------------------------------------------------------
      class(HexMesh), intent(in) :: self
      integer       , intent(in) :: zoneID
      logical                    :: conforming
      !-----------------------------------------------------------
      integer :: fIdx   ! Local face index in zone
      integer :: fID    ! Face index
      integer :: eID    ! Element index
      integer :: eSide  ! Side of the element in contact with boundary
      integer :: nFace  ! Counter for neighbor faces
      !-----------------------------------------------------------

      if (zoneID < lbound(self % zones,1) .or. zoneID > ubound(self % zones,1) ) error stop 'HexMesh_ConformingOnZone :: Out of bounds'

      conforming = .TRUE.
      do fIdx = 1, self % zones(zoneID) % no_of_faces

         fID   = self % zones(zoneID) % faces(fIdx)
         eID   = self % faces(fID) % elementIDs(1)
         eSide = self % faces(fID) % elementSide(1)

         ! loop over the faces that are shared between boundary elements
         do nFace = 1, 4
            associate (f => self % faces ( self % elements(eID) % faceIDs (neighborFaces(nFace,eSide) ) ) )

            if (f % FaceType == HMESH_BOUNDARY) cycle

            if (any(f % NfLeft /= f % NfRight)) then
               conforming = .FALSE.
               return
            end if

            end associate
         end do

      end do

   end function
#if defined(INCNS) && defined(CAHNHILLIARD)
   subroutine HexMesh_ConvertDensityToPhaseField(self)
!
!     *************************************************************
!     Convert density to phase field only in element interior nodes
!     *************************************************************
!
      implicit none
      class(HexMesh),   intent(inout)  :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: eID, fID

      associate(rho1 => dimensionless % rho(1), &
                rho2 => dimensionless % rho(2))
      do eID = 1, self % no_of_elements
         associate(c => self % elements(eID) % storage % c, &
                   Q => self % elements(eID) % storage % QNS)
         c(1,:,:,:) = (-rho1 - rho2 + 2.0_RP * Q(INSRHO,:,:,:))/(rho2-rho1)
         end associate
      end do

      end associate
   end subroutine HexMesh_ConvertDensityToPhaseField

   subroutine HexMesh_ConvertPhaseFieldToDensity(self)
!
!     *************************************************************
!     Convert density to phase field only in element interior nodes
!     *************************************************************
!
      implicit none
      class(HexMesh),   intent(inout)  :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: eID, fID

      associate(rho1 => dimensionless % rho(1), &
                rho2 => dimensionless % rho(2))
      do eID = 1, self % no_of_elements
         associate(c => self % elements(eID) % storage % c, &
                   Q => self % elements(eID) % storage % QNS)
         Q(INSRHO,:,:,:) = 0.5_RP*(rho1*(1.0_RP-c(1,:,:,:)) + rho2*(1.0_RP + c(1,:,:,:)))
         end associate
      end do

      end associate

   end subroutine HexMesh_ConvertPhaseFieldToDensity
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  Adapts a mesh to new polynomial orders NNew
!  --------------------------------------------------------
subroutine HexMesh_pAdapt_MPI (self, NNew, controlVariables)
   implicit none
   !-arguments-----------------------------------------
   class(HexMesh), target  , intent(inout)   :: self
   integer                 , intent(in)      :: NNew(NDIM,self % no_of_elements)
   type(FTValueDictionary) , intent(in)      :: controlVariables
   !-local-variables-----------------------------------
   integer :: eID, fID, STLNum
   logical :: saveGradients, FaceComputeQdot
   logical :: analyticalJac   ! Do we need analytical Jacobian storage?
   type(Element)   , pointer :: e
   type(Face)      , pointer :: f

#if (!defined(NAVIERSTOKES))
   logical, parameter            :: computeGradients = .true.
#endif
   !---------------------------------------------------

!     **************************************
!     Check if resulting mesh is anisotropic
!     **************************************
   if ( maxval(NNew) /= minval(NNew) ) self % anisotropic = .TRUE.

   self % NDOF = 0
   do eID=1, self % no_of_elements
      self % NDOF = self % NDOF + product( NNew(:,eID) + 1 )
   end do

!     ********************
!     Some initializations
!     ********************
   saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
   FaceComputeQdot = controlVariables % containsKey("acoustic analogy")
   analyticalJac  = self % storage % anJacobian

!     ************************
!     Clean IBM Mask
!     ************************
   if (self % IBM% active) then
      do STLNum = 1, self% IBM% NumOfSTL
         call self% IBM% CleanMask( self % elements, self % no_of_elements, STLNum )
      end do
   end if

!     *********************************************
!     Adapt individual elements (geometry excluded)
!     *********************************************
!$omp parallel do schedule(runtime) private(e)
   do eID=1, self % no_of_elements
      e => self % elements(eID)   ! Associate fails(!) here
      if ( all( e % Nxyz == NNew(:,eID)) ) then
         cycle
      else
         call e % pAdapt ( NNew(:,eID), self % nodeType, saveGradients, self % storage % prevSol_num )
!$omp critical
         self % Nx(eID) = NNew(1,eID)
         self % Ny(eID) = NNew(2,eID)
         self % Nz(eID) = NNew(3,eID)
!$omp end critical
      end if
   end do
!$omp end parallel do    

!     *************************
!     Adapt corresponding faces
!     *************************

!     Destruct faces storage
!     ----------------------
!$omp parallel do schedule(runtime) 
   do fID=1, self % no_of_faces  !Destruct All faces storage
      call self % faces(fID) % storage % destruct
   end do
!$omp end parallel do

!     Set connectivities
!     ------------------
   call self % SetConnectivitiesAndLinkFaces (self % nodeType)

!     Construct faces storage
!     -----------------------
!$omp parallel do schedule(runtime) private(f)
   do fID=1, self % no_of_faces  !Construct All faces storage
      f => self % faces( fID )
      call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, analyticalJac, FaceComputeQdot)
      call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, analyticalJac, FaceComputeQdot)
   end do
!$omp end parallel do

!     ********************
!     Reconstruct geometry
!     ********************

   !* 1. Adapted elements
   !* 2. Surrounding faces of adapted elements
   !* 3. Neighbor elements of adapted elements whose intermediate face's geometry was adapted
   !* 4. Faces and elements that share a boundary with a reconstructed face (3D non-conforming representations)

!     Destruct old
!     ------------
   do eID=1, self % no_of_elements
      call self % elements (eID) % geom % destruct
   end do
   
   do eID=1, self % no_of_elements 
      e => self % elements(eID)
      do fID=1, 6
         call self % faces( e % faceIDs(fID) ) % geom % destruct
      end do
   end do

!     Construct new
!     ------------

   call self % ConstructGeometry()

!     ************************
!     Construct IBM Mask
!     ************************
   if (self % IBM% active) then
!$omp parallel do schedule(runtime) private(e)
      do eID=1, self % no_of_elements
         e => self % elements(eID) 
         if (allocated(e% isInsideBody)) deallocate(e% isInsideBody)
         if (allocated(e% isForcingPoint)) deallocate(e% isForcingPoint)
         if (allocated(e% STL)) deallocate(e% STL)

         call e % ConstructIBM(e% Nxyz(1), e% Nxyz(2), e% Nxyz(3), self% IBM% NumOfSTL)
      end do
!$omp end parallel do 

      do STLNum = 1, self% IBM% NumOfSTL
         call self% IBM% build(self % elements, self % no_of_elements, self % NDOF, .false.)
      end do
   end if

#if defined(NAVIERSTOKES)
   call self % ComputeWallDistances()
#endif

!     *********
!     Finish up
!     *********
   call self % PrepareForIO
   nullify (e)
   nullify (f)

end subroutine HexMesh_pAdapt_MPI
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Adapts a mesh to new polynomial orders NNew. 
!  Optimized version of HexMesh_Adapt_MPI, but this one does not work with MPI
!  This subroutine is required for pAdaptationClassTE
!  ------------------------------------------------------------------------
   subroutine HexMesh_pAdapt (self, NNew, controlVariables)
      implicit none
      !-arguments-----------------------------------------
      class(HexMesh), target  , intent(inout)   :: self
      integer                 , intent(in)      :: NNew(NDIM,self % no_of_elements)
      type(FTValueDictionary) , intent(in)      :: controlVariables
      !-local-variables-----------------------------------
      integer :: eID, fID, zoneID
      logical :: saveGradients, FaceComputeQdot
      logical :: analyticalJac   ! Do we need analytical Jacobian storage?
      type(IntegerDataLinkedList_t) :: elementList
      type(IntegerDataLinkedList_t) :: facesList
      type(IntegerDataLinkedList_t) :: zoneList
      integer         , allocatable :: zoneArray(:)
      integer         , allocatable :: facesArray(:)
      integer         , allocatable :: elementArray(:)
      type(Zone_t)    , pointer :: zone
      type(Element)   , pointer :: e
      type(Face)      , pointer :: f
#if (!defined(NAVIERSTOKES))
      logical, parameter            :: computeGradients = .true.
#endif
      !---------------------------------------------------

!     **************************************
!     Check if resulting mesh is anisotropic
!     **************************************
      if ( maxval(NNew) /= minval(NNew) ) self % anisotropic = .TRUE.

      self % NDOF = 0
      do eID=1, self % no_of_elements
         self % NDOF = self % NDOF + product( NNew(:,eID) + 1 )
      end do

!     ********************
!     Some initializations
!     ********************
      saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
      FaceComputeQdot = controlVariables % containsKey("acoustic analogy")

      facesList      = IntegerDataLinkedList_t(.FALSE.)
      elementList    = IntegerDataLinkedList_t(.FALSE.)
      zoneList = IntegerDataLinkedList_t(.FALSE.)
      analyticalJac  = self % storage % anJacobian

!     *********************************************
!     Adapt individual elements (geometry excluded)
!     *********************************************
!$omp parallel do schedule(runtime) private(fID, e)
      do eID=1, self % no_of_elements
         e => self % elements(eID)   ! Associate fails(!) here
         if ( all( e % Nxyz == NNew(:,eID)) ) then
            cycle
         else
            call e % pAdapt ( NNew(:,eID), self % nodeType, saveGradients, self % storage % prevSol_num )
!$omp critical
            self % Nx(eID) = NNew(1,eID)
            self % Ny(eID) = NNew(2,eID)
            self % Nz(eID) = NNew(3,eID)
            call elementList % add (eID)
            do fID=1, 6
               call facesList   % add (e % faceIDs(fID))
               if (self % faces(e % faceIDs(fID)) % FaceType  /= HMESH_BOUNDARY) then
                  call elementList % add ( mpi_partition % global2localeID (e % Connection(fID) % globID) )
               end if
            end do
!$omp end critical

         end if
!~         end associate
      end do
!$omp end parallel do

      call facesList % ExportToArray(facesArray, .TRUE.)

!     *************************
!     Adapt corresponding faces
!     *************************

!     Destruct faces storage
!     ----------------------
!$omp parallel do schedule(runtime)
      do fID=1, size(facesArray)
         call self % faces( facesArray(fID) ) % storage % destruct
      end do
!$omp end parallel do

!     Set connectivities
!     ------------------
      call self % SetConnectivitiesAndLinkFaces (self % nodeType, facesArray)

!     Construct faces storage
!     -----------------------
!$omp parallel do private(f) schedule(runtime)
      do fID=1, size(facesArray)
         f => self % faces( facesArray(fID) )  ! associate fails here in intel compilers
         call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, analyticalJac, FaceComputeQdot)
         call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, analyticalJac, FaceComputeQdot)
      end do
!$omp end parallel do

!     ********************
!     Reconstruct geometry
!     ********************

      !* 1. Adapted elements
      !* 2. Surrounding faces of adapted elements
      !* 3. Neighbor elements of adapted elements whose intermediate face's geometry was adapted
      !* 4. Faces and elements that share a boundary with a reconstructed face (3D non-conforming representations)


      if (self % anisotropic .and. (.not. self % meshIs2D) ) then

!        Check if any of the faces belongs to a boundary
!        -----------------------------------------------
         do fID=1, size(facesArray)
            associate (f => self % faces( facesArray(fID) ) )
            if ( f % FaceType == HMESH_BOUNDARY ) then
               call zoneList % add (f % zone)
            end if
            end associate
         end do

!        Add the corresponding faces and elements
!        ----------------------------------------
         call zoneList % ExportToArray (zoneArray)

         do zoneID=1, size(zoneArray)
            zone => self % zones( zoneArray(zoneID) )    ! Compiler bug(?): If zone was implemented as associate, gfortran would not compile
            do fID=1, zone % no_of_faces
               call facesList   % add ( zone % faces(fID) )
               call elementList % add ( self % faces(zone % faces(fID)) % elementIDs(1) )
            end do
         end do
         deallocate (zoneArray   )
      end if

      deallocate ( facesArray )

      call facesList   % ExportToArray(facesArray  , .TRUE.)
      call elementList % ExportToArray(elementArray, .TRUE.)

!     Destruct old
!     ------------
      do eID=1, size (elementArray)
         call self % elements (elementArray(eID)) % geom % destruct
      end do
      do fID=1, size (facesArray)
         call self % faces (facesArray(fID)) % geom % destruct
      end do

      call self % ConstructGeometry(facesArray, elementArray)

#if defined(NAVIERSTOKES)
      call self % ComputeWallDistances(facesArray, elementArray)
#endif

!     *********
!     Finish up
!     *********
      call self % PrepareForIO

      call facesList    % destruct
      call elementList  % destruct
      call zoneList     % destruct
      nullify (zone)
      nullify (e)
      nullify (f)
      deallocate (facesArray  )
      deallocate (elementArray)

   end subroutine HexMesh_pAdapt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  HexMesh_Assign:
!  Subroutine to assign a HexMesh to another.
!  It turns out that an "impure" procedure with explicit OMP is more efficient than pure or elemental in this case.
!
   subroutine HexMesh_Assign (to, from)
      implicit none
      !-arguments----------------------------------------
      class(HexMesh), intent(inout), target :: to
      type (HexMesh), intent(in)    :: from
      !-local-variables----------------------------------
      integer :: eID
      !--------------------------------------------------
      to % numberOfFaces      = from % numberOfFaces
      to % nodeType           = from % nodeType
      to % no_of_elements     = from % no_of_elements
      to % no_of_allElements  = from % no_of_allElements
      to % dt_restriction     = from % dt_restriction
      to % NDOF               = from % NDOF

      to % meshFileName       = from % meshFileName
      to % meshIs2D           = from % meshIs2D
      to % dir2D              = from % dir2D

      to % anisotropic        = from % anisotropic
      to % child              = from % child
      to % ignoreBCnonConformities = from % child

      safedeallocate (to % Nx)
      allocate ( to % Nx ( size(from % Nx) ) )
      to % Nx = from % Nx

      safedeallocate (to % Ny)
      allocate ( to % Ny ( size(from % Ny) ) )
      to % Ny = from % Ny

      safedeallocate (to % Nz)
      allocate ( to % Nz ( size(from % Nz) ) )
      to % Nz = from % Nz

      to % storage = from % storage

      safedeallocate (to % nodes)
      allocate ( to % nodes ( size(from % nodes) ) )
!$omp parallel do schedule(runtime)
      do eID=1, size(from % nodes)
         to % nodes(eID) = from % nodes(eID)
      end do
!$omp end parallel do

      safedeallocate (to % faces)
      allocate ( to % faces ( size(from % faces) ) )
!$omp parallel do schedule(runtime)
      do eID=1, size(from % faces)
         to % faces(eID) = from % faces(eID)
      end do
!$omp end parallel do
      safedeallocate (to % elements)
      allocate ( to % elements ( size(from % elements) ) )
!$omp parallel do schedule(runtime)
      do eID=1, from % no_of_elements
         to % elements(eID) = from % elements(eID)
      end do
!$omp end parallel do

      to % MPIfaces = from % MPIfaces

      safedeallocate (to % zones)
      allocate ( to % zones ( size(from % zones) ) )
      to % zones = from % zones

!
!     Point elements' storage
!     -----------------------
!$omp parallel do schedule(runtime)
      do eID = 1, to % no_of_elements
         to % elements(eID) % storage => to % storage % elements(eID)
      end do
!$omp end parallel do

      safedeallocate(to % elements_sequential)
      safedeallocate(to % elements_mpi       )
      safedeallocate(to % faces_interior     )
      safedeallocate(to % faces_mpi          )
      safedeallocate(to % faces_boundary     )

      allocate(to % elements_sequential(size(from % elements_sequential)))
      allocate(to % elements_mpi       (size(from % elements_mpi       )))
      allocate(to % faces_interior     (size(from % faces_interior     )))
      allocate(to % faces_mpi          (size(from % faces_mpi          )))
      allocate(to % faces_boundary     (size(from % faces_boundary     )))

      to % elements_sequential = from % elements_sequential
      to % elements_mpi        = from % elements_mpi
      to % faces_interior      = from % faces_interior
      to % faces_mpi           = from % faces_mpi
      to % faces_boundary      = from % faces_boundary


   end subroutine HexMesh_Assign
   subroutine HexMesh_UpdateHOArrays(self)
      implicit none
      !-arguments-----------------------------------------
      class(HexMesh), target  , intent(inout)   :: self
      !-local-variables-----------------------------------
      type(IntegerDataLinkedList_t)         :: HO_elementList
      type(IntegerDataLinkedList_t)         :: LO_elementList
      type(IntegerDataLinkedList_t)         :: faceInteriorList
      type(IntegerDataLinkedList_t)         :: faceBoundaryList
      type(IntegerDataLinkedList_t)         :: elementMPIList
      type(IntegerDataLinkedList_t)         :: elementSequentialList
      integer                               :: eID, fID, face, i
      !--------------------------------------------------
      
      HO_elementList            = IntegerDataLinkedList_t(.FALSE.)
      LO_elementList            = IntegerDataLinkedList_t(.FALSE.)
      faceInteriorList          = IntegerDataLinkedList_t(.FALSE.)
      faceBoundaryList          = IntegerDataLinkedList_t(.FALSE.)
      elementMPIList            = IntegerDataLinkedList_t(.FALSE.)
      elementSequentialList     = IntegerDataLinkedList_t(.FALSE.)

      if( allocated(self % HO_Elements) ) then
         deallocate(self % HO_Elements)
      endif

      if ( allocated(self % LO_Elements)) then
         deallocate(self % LO_Elements)
      endif

      if( allocated(self % HO_FacesInterior) ) then
         deallocate(self % HO_FacesInterior)
      endif

      if( allocated(self % HO_FacesBoundary) ) then
         deallocate(self % HO_FacesBoundary)
      endif

      if( allocated(self % HO_ElementsSequential) ) then
         deallocate(self % HO_ElementsSequential)
      endif

      ! All elements and faces
      do eID=1, self % no_of_elements
         if ( maxval( self % elements(eID) % Nxyz) > 1 ) then
            call HO_elementList % add(eID)
            do fID=1, 6
               face = self % elements(eID) % faceIDs(fID)
               if (self % faces(face) % FaceType  /= HMESH_BOUNDARY) then
                  call faceInteriorList % add (face)
               else
                  call faceBoundaryList % add (face)
               end if
            end do
         else
            call LO_elementList % add(eID)
         endif
      end do

      ! Sequential elements
      do i=1, size(self % elements_sequential)
         eID = self % elements_sequential(i)
         if ( maxval( self % elements(eID) % Nxyz) > 1 ) then
            call elementSequentialList % add(eID)
         endif
      end do


      call HO_elementList        % ExportToArray(self % HO_Elements, .TRUE.)
      call LO_elementList        % ExportToArray(self % LO_Elements, .TRUE.)
      call faceInteriorList      % ExportToArray(self % HO_FacesInterior, .TRUE.)
      call faceBoundaryList      % ExportToArray(self % HO_FacesBoundary, .TRUE.)
      call elementSequentialList % ExportToArray(self % HO_ElementsSequential, .TRUE.)

      call HO_elementList        % destruct
      call LO_elementList        % destruct
      call faceInteriorList      % destruct
      call faceBoundaryList      % destruct
      call elementSequentialList % destruct

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         if( allocated(self % HO_ElementsMPI) ) then
            deallocate(self % HO_ElementsMPI)
         endif

         do i=1, size(self % elements_mpi)
            eID = self % elements_mpi(i)
            if ( maxval( self % elements(eID) % Nxyz) > 1 ) then
               call elementMPIList % add(eID)
            endif
         end do

         call elementMPIList % ExportToArray(self % HO_ElementsMPI, .TRUE.)
      end if
#endif
call elementMPIList % destruct 

   end subroutine HexMesh_UpdateHOArrays

   subroutine HexMesh_RotateMesh(self, rad, center, numBFacePoints, nodes, mpi, angle )
      USE Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      IMPLICIT NONE 
      class(HexMesh), intent(inout)  :: self
      real(kind=RP), intent(inout) :: rad 
      real(kind=RP), intent(inout) :: center(2) 
      !integer, intent(inout)       :: dir2D
      integer, intent(inout)       :: numBFacePoints
      integer                          :: nodes
      logical, intent(in)          :: mpi 
      real(kind=RP), intent(in), optional :: angle

      integer :: n_sliding
      integer :: n_slidingnewnodes 
      integer, allocatable       :: arr1(:)
      integer, allocatable       :: arr2(:)
      integer, allocatable       :: arr3(:)
      integer, allocatable       :: mortararr1(:,:)
      integer, allocatable       :: mortararr2(:,:)
      integer, allocatable       ::face_nodes(:,:)
      integer, allocatable       ::face_othernodes(:,:)
      integer, allocatable       :: Mat(:,:)
      integer, allocatable       :: Connect(:,:,:)
      integer, allocatable       :: rotmortars(:)
      type(node), allocatable   :: tmp_nodes(:)
      integer ::  inter 
      real(kind=RP), allocatable :: o(:)
      real(kind=RP), allocatable :: s(:)
      real(kind=RP):: th
      integer :: oldnode
      integer :: numberOfNodes
      integer :: l, i, j 
      logical :: success
      integer :: oldnfaces
      real(kind=rp) :: PI 
      real(kind=rp) :: theta
      logical :: confor=.FALSE.

      call self % MarkRadius(rad, center, n_sliding, n_slidingnewnodes)

      !write(*,*) 'n_sliding',n_sliding
      !write(*,*) 'n_sliding',n_slidingnewnodes

      !reconstruct the nodes 
      if (.NOT. self%sliding) then 
         oldnode=size(self%nodes)
         allocate(tmp_nodes(oldnode))
         do i=1, oldnode
            CALL ConstructNode( tmp_nodes(i), self%nodes(i)%x, self%nodes(i)%globID )
         end do 
         do i=1, oldnode
            call self % nodes(i)%destruct
         end do 
         safedeallocate (self%nodes)
         numberOfNodes=(size(tmp_nodes)+2*n_slidingnewnodes)
         allocate(self%nodes(size(tmp_nodes)+2*n_slidingnewnodes))
         do i=1, oldnode
            CALL ConstructNode( self%nodes(i), tmp_nodes(i)%x, tmp_nodes(i)%globID )
         end do 
         safedeallocate (tmp_nodes)
      end if 
      allocate(arr1(n_slidingnewnodes))
      allocate(arr2(n_slidingnewnodes))
      allocate(mortararr1(n_slidingnewnodes,2))
      allocate(mortararr2(n_slidingnewnodes,2))
      allocate(arr3(n_sliding))
      allocate(Mat(n_slidingnewnodes,12))
      allocate(face_nodes(n_slidingnewnodes,4))
      allocate(face_othernodes(n_slidingnewnodes,4))
      allocate(Connect(n_slidingnewnodes, 9, 6))
      allocate(rotmortars(2*n_slidingnewnodes))
      arr1=0
      arr2=0
      arr3=0
      Mat=0
      Connect=0
      face_nodes=0
      face_othernodes=0
      allocate(o(4))
      allocate(s(4))
      PI=4.0_RP*DATAN(1.0_RP)
      theta=PI/2000.0_RP
      if (.not.allocated(self%arr1) ) allocate (self%arr1(n_slidingnewnodes))
      if (.not.allocated(self%arr2))  allocate (self%arr2(n_slidingnewnodes))
      if (.not.allocated(self%mortararr1))  allocate (self%mortararr1(n_slidingnewnodes,2))
      if (.not.allocated(self%mortararr2) ) allocate (self%mortararr2(n_slidingnewnodes,2))
      if (.not.allocated(self%arr3))  allocate (self%arr3(n_sliding))
      if (.not.allocated(self%Mat))  allocate (self%Mat(n_slidingnewnodes,12))
      if (.not.allocated(self%face_nodes))  allocate (self%face_nodes(n_slidingnewnodes,4))
      if (.not.allocated(self%face_othernodes))  allocate (self%face_othernodes(n_slidingnewnodes,4))
      if (.not.allocated(self%Connect))  allocate (self%Connect(n_slidingnewnodes, 9, 6))
      if (.not.allocated(self%rotmortars))  allocate (self%rotmortars(2*n_slidingnewnodes))

if (present(angle))theta=angle
      call self % Modifymesh(nodes,n_slidingnewnodes, arr1, arr2, arr3,Mat,center, o, s, face_nodes,face_othernodes, rotmortars, th, numberOfNodes,numBFacePoints,oldnode, theta)
      if (.not.self%sliding) then 
         self%arr1=arr1
         self%arr2=arr2 
         self%arr3=arr3
         self%Mat=Mat 
         self%mortararr1=mortararr1 
         self%mortararr2=mortararr2
         self%face_nodes=face_nodes
         self%face_othernodes=face_othernodes
         self%connect=connect 
         self%rotmortars=rotmortars
         self%n_slidingnewnodes=n_slidingnewnodes
         self%n_sliding=n_sliding
      end if 
      do l=1, size(arr2)
         do j=1,6
            if (self%elements(self%arr2(l))%MortarFaces(j)==1) then 
               self%mortararr2(l,1)=self%arr2(l)
               self%mortararr2(l,1)=self%Mat(l,1)
               self%mortararr2(l,2)=self%Mat(l,4)
               self%elements(self%arr2(l))%MortarFaces=0
            end if 
            if (self%elements(self%arr1(l))%MortarFaces(j)==1) then 
               self%mortararr1(l,1)=self%arr1(l)
               self%mortararr1(l,1)=self%Mat(l,3)
               self%mortararr1(l,2)=j
               self%mortararr1(l,2)=self%Mat(l,6)
               self%elements(self%arr1(l))%MortarFaces=0
            end if 
         end do 
      end do 

      oldnfaces=size(self%faces)
      do i=1,size(self%faces)
         call self%faces(i)%Destruct
      end do 
      if (.not. self%sliding) then 
         safedeallocate (self%faces)
         allocate(self%faces(oldnfaces + n_slidingnewnodes))
      end if 
      CALL ConstructFaces( self, success )
      if (allocated(self % zones) ) then 
         deallocate(self % zones)
      end if 
      call self % ConstructZones()

      CALL getElementsFaceIDs(self)
      call self % DefineAsBoundaryFaces()
      i=0
      
      if ( .not. MPI_Process % doMPIRootAction ) then
         call self % CheckIfMeshIs2D()
      end if

      !if ( dir2D .ne. 0 ) then
      !   call SetMappingsToCrossProduct
       !  call self % CorrectOrderFor2DMesh(dir2D)
      !end if

  

      do l=1, size(arr2) 
            self%faces(self%elements(self%arr2(l))%faceIDs( self%mortararr2(l,2)))%IsMortar=3
            self%faces(self%elements(self%arr2(l))%faceIDs( self%mortararr2(l,2)))%facetype=1

            self%faces(self%elements(self%arr1(l))%faceIDs( self%mortararr1(l,2)))%IsMortar=3
            self%faces(self%elements(self%arr1(l))%faceIDs( self%mortararr1(l,2)))%facetype=1
      end do 
      do l=1, size(self%faces)
         if (self % faces(l) % faceType==-1) then 
            self % faces(l) % faceType=1
         end if 
      end do

     call self % SetConnectivitiesAndLinkFaces(nodes) 

     call self % ConstructGeometry()

     PI=4.0_RP*DATAN(1.0_RP)

     if ( th .EQ. 0.0_RP ) confor=.TRUE. 
     write(*,*) 'confor is', confor 
      call self % ConstructMortars(nodes, self%n_slidingnewnodes, self%arr1, self%arr2,self%Mat,o, s, self%mortararr2,self%rotmortars, th, confor)

    ! if (th .EQ. PI/20_RP) then 
    !  write(*,*) 'sliding conforming'
    !  call self % ConstructSlidingMortarsConforming(nodes, n_slidingnewnodes, arr1, arr2,Mat, o, s, mortararr2,rotmortars)
    ! end if 
     self%sliding=.true.
     self%omega=self%omega+th
     deallocate(arr1)
     deallocate(arr2)
     deallocate(arr3)
     deallocate(Mat)
     deallocate(o)
     deallocate(s)
     deallocate(mortararr1)
     deallocate(mortararr2)
     deallocate(face_nodes)
     deallocate(face_othernodes)
     deallocate(Connect)
     deallocate(rotmortars)
   end subroutine HexMesh_RotateMesh


  subroutine HexMesh_MarkRadius (self, rad, center, n_sliding, n_slidingnewnodes)
   IMPLICIT NONE 
   class(HexMesh), intent(inout)  :: self
   real(kind=RP), intent(in) :: rad 
   real(kind=RP), intent(in) :: center(2) 
   integer, intent(inout) :: n_sliding
   integer, intent(inout) :: n_slidingnewnodes

   integer :: eID, eID2
   integer ::  i, f, j, ll, l
   integer :: new_nFaces

   f=0
   n_slidingnewnodes=0
   n_sliding=0
   new_nFaces=SIZE(self % faces)
   do i=1, size(self%elements)
      f=0
      do j=1,8
         if (((self%Nodes(self%elements(i)%nodeIDs(j))%X(1)-center(1))**2 +&
         (self%Nodes(self%elements(i)%nodeIDs(j))%X(3)-center(2))**2) .le. rad**2)  then 
            f=f+1
         end if 
      end do  
      if (f==8) then
         self%elements(i)%sliding=.true.
         n_sliding=n_sliding+1
      end if 
   end do 

   l=0
   do i=1, size(self%elements)
      if (self%elements(i)%sliding) then 
         eID=self%elements(i)%eID
         do j=1,6
           if (self%faces(self%elements(i)%faceIDs(j))%elementIDs(1)==eID) then 
            eID2=self%faces(self%elements(i)%faceIDs(j))%elementIDs(2)
           else   
            eID2=self%faces(self%elements(i)%faceIDs(j))%elementIDs(1)
           end if 
           if (eID2 .ne. 0) then 
               if (self%elements(eID2)%sliding) then 
                  cycle
               else 
                  self%elements(eID)%sliding_newnodes=.true.
                  n_slidingnewnodes=n_slidingnewnodes+1
                  new_nFaces=new_nFaces + 1
               end if 
           end if 
         end do 
      end if 
   end do 

  end subroutine HexMesh_MarkRadius
   subroutine HexMesh_MarkSlidingElementsRadius(self, rad, nelm, center, arr3, arr2, arr1, Mat, Connect,new_nFaces, face_nodes,face_othernodes,rotmortars)
   IMPLICIT NONE 
   class(HexMesh), intent(inout)  :: self

   real(kind=RP), intent(in) :: rad 
   real(kind=RP), intent(in) :: center(2)
   integer, intent(in)    :: nelm
   integer, intent(inout) :: arr3(nelm) !32 !384
   integer, intent(inout) :: arr2(nelm) !32   !48
   integer, intent(inout) :: arr1(nelm)  !32  !48
   integer, intent(inout) :: Mat(nelm,12) !32   !48 9
   integer, intent(inout) :: Connect(nelm, 9, 6)  !32  !48 9 6 
   !integer, intent(inout) :: Connect(48, 3, 6)
   integer, intent(inout) :: new_nFaces 
   integer, intent(inout) :: face_nodes(nelm,4)!32 
   integer, intent(inout) :: face_othernodes(nelm,4)!32 
   integer, intent(inout) :: rotmortars(nelm*2) !64 !96

   integer :: eID, eID2, eID3, e, m, fID2, faceNumber, faceNumber1, faceNumber2
   integer ::  i,f, ff ,fff, z, j, jj, jjj, ll,l, lll, k, kk, no, v,zz, f1,f2,ind1,ind2,ind3,sp
   integer ::ar2(nelm)
   integer :: faceNodeIDs(4)
   integer :: faceNodeIDs1(4)
   integer :: faceNodeIDs2(4)
   integer :: nodeIDs(8)
   integer :: nodeIDs1(8)
   integer :: nodeIDs2(8)

   real(kind=RP) :: NODESS(4,3)
   real(kind=RP) :: NODES2(4,3)

   integer :: rota 

   real(kind=RP)              :: XR(4,3)
   real(kind=RP) :: MNODES(nelm,4,3) !32   !48 5 3
   real(kind=RP) :: NODES(4,3)!32 
   real(kind=RP) :: theta 
   real(KIND=RP) :: ROT(3,3)
   integer :: Matrot(nelm,9) !32   !48 9 
   integer :: Matrott(nelm,9) !32   !48 9 
   integer :: facer(nelm,2)!32  !48 2
   integer :: y 
   y=3
   facer=0
   f=0
   l=0
   ll=0
   new_nFaces=SIZE(self % faces)
   do i=1, size(self%elements)
      f=0
      do j=1,8
         if (((self%Nodes(self%elements(i)%nodeIDs(j))%X(1)-center(1))**2 +&
         (self%Nodes(self%elements(i)%nodeIDs(j))%X(y)-center(2))**2) .le. rad**2)  then 
            f=f+1
         end if 
      end do  
      if (f==8) then
         self%elements(i)%sliding=.true.
         ll=ll+1
      end if 
   end do 

   l=0
   do i=1, size(self%elements)
      if (self%elements(i)%sliding) then 
         eID=self%elements(i)%eID
         do j=1,6
           if (self%faces(self%elements(i)%faceIDs(j))%elementIDs(1)==eID) then 
            eID2=self%faces(self%elements(i)%faceIDs(j))%elementIDs(2)
           else   
            eID2=self%faces(self%elements(i)%faceIDs(j))%elementIDs(1)
           end if 
           if (eID2 .ne. 0) then 
               if (self%elements(eID2)%sliding) then 
                  cycle
               else 
                  self%elements(eID)%MortarFaces(j)=1
                  self%elements(eID)%sliding_newnodes=.true.
                  l=l+1
                  arr2(l)=eID
                  new_nFaces=new_nFaces + 1
               end if 
           end if 
         end do 
      end if 
   end do 

   ll=0
   l=0
  
   !do i=1, size(arr2)
   !   write(*,*) 'arr2 i=',i, arr2(i)
   !end do 
   Mat=0
   do i=1,size(arr2)
      Mat(i,1)=arr2(i)
      ll=0
      lll=0
      no=0
      !check the position of the element
      do l=1,8
          if (self%nodes(self%elements(arr2(i))%nodeIDs(l))%X(y) .GT. 0.0_RP)  then 
              do z=1,8 
                 if (self%nodes(self%elements(arr2(i))%nodeIDs(z))%X(y) .LT. 0.0_RP) then 
                    no=1
                 end if 
              end do 
              if (no==0) then 
                 ll=1
                 exit 
              end if 
           end if 
  
           if (self%nodes(self%elements(arr2(i))%nodeIDs(l))%X(y) .LT. 0.0_RP) then 
              lll=1
              exit 
           end if 
      end do 
      !it's on the top of the circle 
      if (ll==1) then 
          do j=1,6
              if (self%faces(self%elements(arr2(i))%faceIDs(j))%elementIDs(1)==arr2(i)) then 
                  eID3=self%faces(self%elements(arr2(i))%faceIDs(j))%elementIDs(2)
               else   
                  eID3=self%faces(self%elements(arr2(i))%faceIDs(j))%elementIDs(1)
               end if
               if (eID3==0) cycle 
               if (.not.self%elements(eID3)%sliding) then 
                  Mat(i,4)=j 
                  Mat(i,7)=self%faces(self%elements(arr2(i))%faceIDs(j))%rotation 
                  select case (j)
                  case(1)
                      face_nodes(i,1)=1
                      face_nodes(i,2)=2
                      face_nodes(i,3)=5   !!!
                      face_nodes(i,4)=6
                      face_othernodes(i,1)=3
                      face_othernodes(i,2)=4
                      face_othernodes(i,3)=7
                      face_othernodes(i,4)=8
                  case(2)
                      face_nodes(i,1)=4
                      face_nodes(i,2)=3
                      face_nodes(i,3)=7
                      face_nodes(i,4)=8
                      face_othernodes(i,1)=1
                      face_othernodes(i,2)=2
                      face_othernodes(i,3)=6
                      face_othernodes(i,4)=5
                  case(3)
                      face_nodes(i,1)=1
                      face_nodes(i,2)=2
                      face_nodes(i,3)=3
                      face_nodes(i,4)=4
                      face_othernodes(i,1)=5
                      face_othernodes(i,2)=6
                      face_othernodes(i,3)=7
                      face_othernodes(i,4)=8
                  case(4)
                      face_nodes(i,1)=6
                      face_nodes(i,2)=2
                      face_nodes(i,3)=3
                      face_nodes(i,4)=7
                      face_othernodes(i,1)=5
                      face_othernodes(i,2)=1
                      face_othernodes(i,3)=4
                      face_othernodes(i,4)=8
                  case(5)
                      face_nodes(i,1)=5
                      face_nodes(i,2)=6
                      face_nodes(i,3)=7
                      face_nodes(i,4)=8
                      face_othernodes(i,1)=1
                      face_othernodes(i,2)=2
                      face_othernodes(i,3)=3
                      face_othernodes(i,4)=4
                  case(6)
                      face_nodes(i,1)=5
                      face_nodes(i,2)=1
                      face_nodes(i,3)=4
                      face_nodes(i,4)=8
                      face_othernodes(i,1)=6
                      face_othernodes(i,2)=2
                      face_othernodes(i,3)=3
                      face_othernodes(i,4)=7
                  end select 
               end if 
               if ((eID3 .NE. 0) .AND. (self%elements(eID3)%sliding_newnodes)) then
                  if (Mat(i,2)==0) then
                  do k=1,8
                     if (Mat(i,2)==0) then
                      do kk=1,8!!!!!!!!
                        if (Mat(i,2)==0) then
                        sp=0
               
                        if (self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1) .GE. center(1)) then 
                          if ((self%nodes(self%elements(eID3)%nodeIDs(k))%X(1) .LT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1)) .AND. &
                          (self%nodes(self%elements(eID3)%nodeIDs(kk))%X(y) .GE. center(2)) .AND. (self%nodes(self%elements(eID3)%nodeIDs(k))%X(y) .GT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(3))  .AND. (Mat(i,3)==0)) then 
                           Mat(i,2)=eID3 
                           exit
                          end if 
                        end if 
                        if (self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1) .LE. center(1)) then 
                           if ((self%nodes(self%elements(eID3)%nodeIDs(k))%X(1) .LT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1)) .AND. &
                           (self%nodes(self%elements(eID3)%nodeIDs(kk))%X(y) .GE. center(2)) .AND. (self%nodes(self%elements(eID3)%nodeIDs(k))%X(y) .LT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(3))  .AND. (Mat(i,3)==0)) then 
                            Mat(i,2)=eID3 
                            exit
                           end if 
                         end if 

                          end if 
                      end do 
                     end if 
                      if (Mat(i,2).NE. 0) exit
                  end do
               end if 
               end if 
               !if (Mat(i,1)==125) then 
               !   Mat(i,2)=10
               !   Mat(i,3)=9
               !end if 
          end do 
          !now we have the sliding neighbour element, we search for the non sliding neighbour (mortar)
          do j=1,6
            if(Mat(i,2) .NE. 0) then 
              if (self%faces(self%elements(Mat(i,2))%faceIDs(j))%elementIDs(1)==Mat(i,2)) then 
                  eID3=self%faces(self%elements(Mat(i,2))%faceIDs(j))%elementIDs(2)
               else   
                  eID3=self%faces(self%elements(Mat(i,2))%faceIDs(j))%elementIDs(1)
               end if
               if ((eID3 .NE. 0) ) then 
                  if  (.NOT.self%elements(eID3)%sliding) then 
                     Mat(i,3)=eID3
                     arr1(i)=eID3
                     Mat(i,5)=j 
                     Mat(i,8)=self%faces(self%elements(Mat(i,2))%faceIDs(j))%rotation
                  end if  
               end if    
            end if 
          end do 
          do j=1,6
            if (Mat(i,3) .NE. 0) then 
              if (self%faces(self%elements(Mat(i,3))%faceIDs(j))%elementIDs(1)==Mat(i,3)) then 
                  eID3=self%faces(self%elements(Mat(i,3))%faceIDs(j))%elementIDs(2)
               else   
                  eID3=self%faces(self%elements(Mat(i,3))%faceIDs(j))%elementIDs(1)
               end if
               if (eID3==Mat(i,2)) then 
                  Mat(i,6)=j 
                  self%elements(Mat(i,3))%MortarFaces(j)=1 

                  Mat(i,9)=self%faces(self%elements(Mat(i,3))%faceIDs(j))%rotation
               end if 
            end if 
          end do 
      end if 
      !it's on the bottom of the circle 
      if (lll==1) then
      do j=1,6
          if (self%faces(self%elements(arr2(i))%faceIDs(j))%elementIDs(1)==arr2(i)) then 
              eID3=self%faces(self%elements(arr2(i))%faceIDs(j))%elementIDs(2)
           else   
              eID3=self%faces(self%elements(arr2(i))%faceIDs(j))%elementIDs(1)
           end if
           if (eID3==0) cycle 
           if (.not.self%elements(eID3)%sliding) then 
              Mat(i,4)=j 
              Mat(i,7)=self%faces(self%elements(arr2(i))%faceIDs(j))%rotation
              select case (j)
              case(1)
                  face_nodes(i,1)=1
                  face_nodes(i,2)=2
                  face_nodes(i,3)=5   !!!
                  face_nodes(i,4)=6
                  face_othernodes(i,1)=3
                  face_othernodes(i,2)=4
                  face_othernodes(i,3)=7
                  face_othernodes(i,4)=8
              case(2)
                  face_nodes(i,1)=4
                  face_nodes(i,2)=3
                  face_nodes(i,3)=7
                  face_nodes(i,4)=8
                  face_othernodes(i,1)=1
                  face_othernodes(i,2)=2
                  face_othernodes(i,3)=6
                  face_othernodes(i,4)=5
              case(3)
                  face_nodes(i,1)=1
                  face_nodes(i,2)=2
                  face_nodes(i,3)=3
                  face_nodes(i,4)=4
                  face_othernodes(i,1)=5
                  face_othernodes(i,2)=6
                  face_othernodes(i,3)=7
                  face_othernodes(i,4)=8
              case(4)
                  face_nodes(i,1)=6
                  face_nodes(i,2)=2
                  face_nodes(i,3)=3
                  face_nodes(i,4)=7
                  face_othernodes(i,1)=5
                  face_othernodes(i,2)=1
                  face_othernodes(i,3)=4
                  face_othernodes(i,4)=8
              case(5)
                  face_nodes(i,1)=5
                  face_nodes(i,2)=6
                  face_nodes(i,3)=7
                  face_nodes(i,4)=8
                  face_othernodes(i,1)=1
                  face_othernodes(i,2)=2
                  face_othernodes(i,3)=3
                  face_othernodes(i,4)=4
              case(6)
                  face_nodes(i,1)=5
                  face_nodes(i,2)=1
                  face_nodes(i,3)=4
                  face_nodes(i,4)=8
                  face_othernodes(i,1)=6
                  face_othernodes(i,2)=2
                  face_othernodes(i,3)=3
                  face_othernodes(i,4)=7
              end select 
           end if 
           if ((eID3 .NE. 0) .AND. (self%elements(eID3)%sliding_newnodes)) then
              if(Mat(i,2)==0) then
              do k=1,8
               if(Mat(i,2)==0) then
                  do kk=1,8
                     if(Mat(i,2)==0) then 
                     if (self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1) .LE. center(1)) then 
                      if ((self%nodes(self%elements(eID3)%nodeIDs(k))%X(1) .GT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1)) .AND. &
                      (self%nodes(self%elements(eID3)%nodeIDs(kk))%X(y) .LE. center(2)) .AND. (self%nodes(self%elements(eID3)%nodeIDs(k))%X(y) .LT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(3)) ) then
                        Mat(i,2)=eID3
                        exit
                      end if  
                     end if
                     if (self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1) .GE. center(1)) then 
                        if ((self%nodes(self%elements(eID3)%nodeIDs(k))%X(1) .GT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(1)) .AND. &
                        (self%nodes(self%elements(eID3)%nodeIDs(kk))%X(y) .LE. center(2)) .AND. (self%nodes(self%elements(eID3)%nodeIDs(k))%X(y) .GT. self%nodes(self%elements(arr2(i))%nodeIDs(kk))%X(3)) ) then
                          Mat(i,2)=eID3
                          exit
                        end if  
                       end if  
                     end if 
                  end do 
               end if 
               if (Mat(i,2) .NE. 0) exit
              end do
            end if 
           end if 
      end do 
      do j=1,6
         if (Mat(i,2) .NE. 0) then 
          if (self%faces(self%elements(Mat(i,2))%faceIDs(j))%elementIDs(1)==Mat(i,2)) then 
              eID3=self%faces(self%elements(Mat(i,2))%faceIDs(j))%elementIDs(2)
           else   
              eID3=self%faces(self%elements(Mat(i,2))%faceIDs(j))%elementIDs(1)
           end if
            if ((eID3 .NE. 0)) then 
               if  (.NOT.self%elements(eID3)%sliding) then 
               Mat(i,3)=eID3
               arr1(i)=eID3
               Mat(i,5)=j 
               Mat(i,8)=self%faces(self%elements(Mat(i,2))%faceIDs(j))%rotation
               end if 
            end if     
         end if 
      end do 
      do j=1,6
         if (Mat(i,3) .NE. 0) then 
          if (self%faces(self%elements(Mat(i,3))%faceIDs(j))%elementIDs(1)==Mat(i,3)) then 
              eID3=self%faces(self%elements(Mat(i,3))%faceIDs(j))%elementIDs(2)
           else   
              eID3=self%faces(self%elements(Mat(i,3))%faceIDs(j))%elementIDs(1)
           end if
           if (eID3==Mat(i,2)) then 
              Mat(i,6)=j 
              self%elements(Mat(i,3))%MortarFaces(j)=1 
              Mat(i,9)=self%faces(self%elements(Mat(i,3))%faceIDs(j))%rotation
           end if 
         end if 
      end do 
      end if 
      !if (Mat(i,1)==125) then 
      !   Mat(i,2)=10
      !   arr1(i)=9
      !   Mat(i,3)=9
      !   self%elements(Mat(i,3))%MortarFaces(3)=1 
      !end if 
   end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do i=1,size(arr2)
      !write(*,*) 'Mat(i=', Mat(i,:)
      if ((Mat(i,2) .NE. 0) .AND. (Mat(i,3) .NE. 0)) then 
         if (.not.self%elements(Mat(i,1))%sliding_newnodes) write(*,*) 'for i=',i,'(i,1) is not sliding_newnodes'
         if (.not.self%elements(Mat(i,2))%sliding_newnodes) write(*,*) 'for i=',i,'(i,2) is not sliding_newnodes'
         if (self%elements(Mat(i,3))%sliding_newnodes) write(*,*) 'for i',i,'(i,3) is sliding_newnodes'
         if (self%elements(Mat(i,3))%sliding) write(*,*) 'for i',i,'(i,3) is sliding'
         eID=self%faces(self%elements(Mat(i,1))%faceIDs(Mat(i,4)))%elementIDs(1)
         eID2=self%faces(self%elements(Mat(i,1))%faceIDs(Mat(i,4)))%elementIDs(2)
         if (eID==Mat(i,1)) then 
            if (self%elements(eID2)%sliding) write(*,*) 'problem with (i,4)'
         else 
            if (self%elements(eID)%sliding) write(*,*) 'problem with (i,4)'
         end if 
         if ((eID .NE. Mat(i,1)) .AND.(eID2 .NE. Mat(i,1))) write(*,*) 'problem with (i,4) line 6195'
         eID=self%faces(self%elements(Mat(i,2))%faceIDs(Mat(i,5)))%elementIDs(1)
         eID2=self%faces(self%elements(Mat(i,2))%faceIDs(Mat(i,5)))%elementIDs(2)
         if (eID==Mat(i,2)) then 
            if (self%elements(eID2)%sliding) write(*,*) 'problem with (i,5)'
         else 
            if (self%elements(eID)%sliding) write(*,*) 'problem with (i,5)'
         end if 
         if ((eID .NE. Mat(i,2)) .AND.(eID2 .NE. Mat(i,2))) write(*,*) 'problem with (i,5) line 6203'
         eID=self%faces(self%elements(Mat(i,3))%faceIDs(Mat(i,6)))%elementIDs(1)
         eID2=self%faces(self%elements(Mat(i,3))%faceIDs(Mat(i,6)))%elementIDs(2)
         if (eID==Mat(i,3)) then 
            if (.not.self%elements(eID2)%sliding) write(*,*) 'problem with (i,6)'
         else 
            if (.not.self%elements(eID)%sliding) write(*,*) 'problem with (i,6)'
         end if 
         if ((eID .NE. Mat(i,3)) .AND.(eID2 .NE. Mat(i,3))) write(*,*) 'problem with (i,5) line 6211'
      end if 
   end do
   connect=0
   do i=1, size(arr2)
      do j=1,6
         if ((self%faces(self%Elements(Mat(i,1))%faceIDs(j))%elementIDs(1) .NE. 0) .AND. (self%faces(self%Elements(Mat(i,1))%faceIDs(j))%elementIDs(2) .NE. 0)) then 
            eID=self%faces(self%Elements(Mat(i,1))%faceIDs(j))%elementIDs(1)
            if ( .NOT. self%elements(eID)%sliding) then 
                  Mat(i,10)= self%faces(self%elements(Mat(i,1))%faceIDs(j))%elementIDs(1)
                  Mat(i,12)=self%faces(self%elements(Mat(i,1))%faceIDs(j))%rotation 
                  do jj=1,6
                     if (self%faces(self%elements(eID)%faceIDs(jj))%elementIDs(1)==Mat(i,1))  Mat(i,11)=jj 
                     if (self%faces(self%elements(eID)%faceIDs(jj))%elementIDs(2)==Mat(i,1))  Mat(i,11)=jj 
                  end do 
            else  
            eID=self%faces(self%Elements(Mat(i,1))%faceIDs(j))%elementIDs(2)
            if ( .NOT. self%elements(eID)%sliding) then
               Mat(i,10)= self%faces(self%elements(Mat(i,1))%faceIDs(j))%elementIDs(2)
               Mat(i,12)=self%faces(self%elements(Mat(i,1))%faceIDs(j))%rotation 
               do jj=1,6
                  if (self%faces(self%elements(eID)%faceIDs(jj))%elementIDs(1)==Mat(i,1))  Mat(i,11)=jj 
                  if (self%faces(self%elements(eID)%faceIDs(jj))%elementIDs(2)==Mat(i,1))  Mat(i,11)=jj 
               end do 
            end if 
         end if 
      end if 
      end do 
   end do 


   if ( MPI_Process % doMPIAction ) then
      if ( .NOT. mpi_partition % Constructed ) then
         open (unit=10,file="Mat.txt",action="write")
         do i=1,size(arr2)
            write (10,*)  (Mat(i,j), j=1,12)
         end do
         close (10)
      end if 
   end if 

   do i=1, size(arr2) 
      k=2
      Connect(i,1,1)=arr2(i)
      Connect(i,1,2)=i

      do l=1, size(arr2)
         do ll=1, 8
            do lll=1,8
            if ((self % elements (arr2(i)) % nodeIDs(ll)==self % elements (arr2(l)) % nodeIDs(lll)) .AND. &
               arr2(i) .ne. arr2(l) )then
               no=0
               do kk=2,9
                  if (Connect(i,kk,1)==arr2(l)) then 
                     no=1
                  endif 
               end do 
                  if (no .ne. 1) then 
                     Connect(i,k,1)=arr2(l)
                     Connect(i,k,2)=l
                     k=k+1
                  endif 
            end if 
            end do 
         end do 
      end do 
   end do  
   do i=1, size(arr2)
      eID=Connect(i,1,1)
      do z=2,9
         eID2=Connect(i,z,1)
         v=3
         if ((eID .ne. 0) .and. (eID2 .ne. 0)) then 
            do ll=1,8
               do lll=1,8
                  if ((self % elements (eID) % nodeIDs(ll)==self % elements (eID2) % nodeIDs(lll))  )then
                     Connect(i,z,v)=ll
                       v=v+1
                  end if 
               end do
            end do 
         end if 
      end do 
   end do   

   if ( MPI_Process % doMPIAction ) then
      if ( .NOT. mpi_partition % Constructed ) then
         open (unit=10,file="connect.txt",action="write")
         do i=1,size(arr2)
            !do j=1,9
               !do k=1,6
                  write (10,*)  (connect(i,:,:))
             !  end do 
            !end do 
         end do
         close (10)
      end if 
   end if 

   l=0
   do i=1, size(self%elements)
      if ((self%elements(i)%sliding) .and. .not.(self%elements(i)%sliding_newnodes)) then 
         l=l+1
         arr3(l)=i

      end if 
   end do 


   !theta=(4.0_RP*DATAN(1.0_RP))/16.0_RP
   ! ROT=0.0_RP
   ! ROT(1,1)=COS(theta)
   ! ROT(2,2)=COS(theta)
   ! ROT(3,3)=1.0_RP 
   ! ROT(1,2)=-SIN(theta)
   ! ROT(2,1)=SIN(theta)
   ! do i=1,size(arr2)
   !      nodeIDs1=self%elements(Mat(i,3))% nodeIDs    !3
   !      nodeIDs2=self%elements(Mat(i,1))% nodeIDs    !1
   !      faceNumber2=Mat(i,6)          !facenumber of 3
   !      faceNumber1=Mat(i,4)           !facenumber of 1
   !      DO j = 1, 4
   !          faceNodeIDs2(j) = nodeIDs1(localFaceNode(j,faceNumber2))
   !          NODESS(j,:)=self%nodes(faceNodeIDs2(j))%X
   !          faceNodeIDs1(j) = nodeIDs2(localFaceNode(j,faceNumber1))
   !          NODES2(j,:)=self%nodes(faceNodeIDs1(j))%X
   !          XR(j,:)=NODES2 (j,:)
   !          NODES2(j,:)= MATMUL(ROT, XR(j,:))
   !      END DO     
   !      DO j = 1, 4
   !            write(*,*) 'rotated nodes',j, NODES2(j,:)
   !      end do 
   !      DO j = 1, 4
   !         write(*,*) 'master nodes',j, NODESS(j,:)
   !   end do 
        ! rota=faceRotationnodes(masterNodeIDs=NODESS, slaveNodeIDs=NODES2)

   !      if (rota .NE. Mat(i,7)) write(*,*) 'problem rot .NE. MAT:, rot=', rota, 'MAT=', Mat(i,7)
   !      if (rota .EQ. Mat(i,7)) write(*,*) 'no problem', rota, 'MAT=', Mat(i,7)
   !      Mat(i,7)=rota
   !   end do 
  ! do i=1,size(arr2)
 !     write(*,*) 'element',Mat(i,1),'old rotation', Mat(i,7)
 !     write(*,*) 'element',Matrot(i,1),'new rotation', Matrot(i,7)
  !    write(*,*) '***********************************************'
  !!    write(*,*) 'Matrott i;',i, Matrott(i,:)
  !!    if (Matrott(i,7) .NE. Matrot(i,7)) write(*,*) 'it seems whos master and slave importe...'
   !end do

   do i=1, size(arr2)
         if ((Mat(i,2) .NE. 0) .AND. (Mat(i,3) .NE. 0))  then 
         if (self%faces(self%elements(Mat(i,1))%faceIDs(Mat(i,4)))%rotation .NE. Mat(i,7)) write(*,*) 'problem with rotation of mat(i7)'
         if ((self%faces(self%elements(Mat(i,1))%faceIDs(Mat(i,4)))%elementIDs(1) .NE. Mat(i,1)) .AND. (self%faces(self%elements(Mat(i,1))%faceIDs(Mat(i,4)))%elementIDs(2) .NE. Mat(i,1))) write(*,*) 'problem with face i1'
         if (self%faces(self%elements(Mat(i,2))%faceIDs(Mat(i,5)))%rotation .NE. Mat(i,8)) write(*,*) 'problem with rotation of mat(i8)'
         if ((self%faces(self%elements(Mat(i,2))%faceIDs(Mat(i,5)))%elementIDs(1) .NE. Mat(i,2)) .AND. (self%faces(self%elements(Mat(i,2))%faceIDs(Mat(i,5)))%elementIDs(2) .NE. Mat(i,2))) write(*,*) 'problem with face i2'
         if (self%faces(self%elements(Mat(i,3))%faceIDs(Mat(i,6)))%rotation .NE. Mat(i,9)) write(*,*) 'problem with rotation of mat(i9)'
         if ((self%faces(self%elements(Mat(i,3))%faceIDs(Mat(i,6)))%elementIDs(1) .NE. Mat(i,3)) .AND. (self%faces(self%elements(Mat(i,3))%faceIDs(Mat(i,6)))%elementIDs(2) .NE. Mat(i,3))) write(*,*) 'problem with face i3'
         end if          
   end do 
   !Mat=Matrot
   do i=1, size(arr2)
      if ((Mat(i,2) .NE. 0) .AND. (Mat(i,3) .NE. 0)) then 
         nodeIDs1=self%elements(Mat(i,3))% nodeIDs
         nodeIDs2=self%elements(Mat(i,2))% nodeIDs
         faceNumber1=Mat(i,6)
         faceNumber2=Mat(i,5)
         do j=1,4
            faceNodeIDs1(j)=nodeIDs1(localFaceNode(j,faceNumber1))
            faceNodeIDs2(j)=nodeIDs2(localFaceNode(j,faceNumber2))
         end do 
         facer(i,1)=facerotation(masterNodeIDs=faceNodeIDs1, slaveNodeIDs=faceNodeIDs2)
         facer(i,2)=facerotation(masterNodeIDs=faceNodeIDs2, slaveNodeIDs=faceNodeIDs1)
      end if 
   end do

   ind1=0
   ind2=0
   ind3=0
   do i=1,size(arr2)
      ind1=0
      ind2=0
      ind3=0
      do l=1,size(arr2)
         if (Mat(i,1)==Mat(l,1)) then 
            ind1=ind1+1
         end if 
         if (Mat(i,2)==Mat(l,2)) then 
            ind2=ind2+1
         end if 
         if (Mat(i,3)==Mat(l,3)) then 
            ind3=ind3+1
         end if 
      end do 
      if (ind1 .NE. 1) write(*,*) 'ind1 .NE 1 problem', ind1
      if (ind2 .NE. 1) write(*,*) 'ind2 .NE 1 problem' ,ind2
      if (ind3 .NE. 1) write(*,*) 'ind3 .NE 1 problem' ,ind3
   end do 

  ! ind1=0
  ! do i=1,size(arr2)
  !    ind1=0
  !    do l=1,size(arr2)
  !       if (Mat(i,1)==Mat(l,2)) then 
  !          ind1=ind1+1
  !          write(*,*) 'Mat(i,1)=',Mat(i,1),'Mat(l,2)=',Mat(l,2)
  !       end if 
  !    end do 
  !    write(*,*) 'line 6548, i=',i,'ind1=', ind1
  ! end do 
   do i=1,size(arr2)
      if (Mat(i,2) .NE. 0) then 
         Mat(i,8)=self%faces(self%elements(Mat(i,2))%faceIDs(Mat(i,5)))%rotation
      end if 
   end do 
  end subroutine HexMesh_MarkSlidingElementsRadius

  subroutine HexMesh_RotateNodes(self, theta,nelm, n, m , new_nNodes, new_nodes, arr1, arr2, arr3, Connect, o, s , face_nodes,face_othernodes, numBFacePoints, oldnnode)
   IMPLICIT NONE
   class(HexMesh), intent(inout)  :: self
    real(KIND=RP), intent(in)     :: theta
    integer, intent(in)           :: nelm
    integer,     intent(in)       :: n 
    integer,     intent(in)       :: m
    integer, intent(inout)        :: new_nNodes 
    type(Node), intent(inout)     :: new_nodes(new_nNodes)
    integer, intent(in)           :: arr1(nelm)
    integer, intent(in)           :: arr2(nelm)
    integer, intent(in)           :: arr3(nelm)
    integer, intent(in)           :: Connect(nelm, 9, 6)
    real(kind=RP), intent(inout)     :: o(4)
    real(kind=RP), intent(inout)     :: s(4)
    integer, intent(inout) :: face_nodes(nelm,4)
    integer, intent(inout) :: face_othernodes(nelm,4)
    integer, intent(inout) :: numBFacePoints
    integer, intent(inout) :: oldnnode

    real(KIND=RP) :: ROT(3,3)
    real(KIND=RP) :: XYZ(8,3)
    real(KIND=RP) :: XR(3)
    real(KIND=RP) :: Xflat(3,2,2)
    real(KIND=RP), allocatable :: Xpatch(:,:,:)
    integer :: i, l, eID, eID2, nm, j, k, kk, z , ll , lll, kkk, zz
    real(kind=RP) :: ss(2), oo(2), x, lb, ls, lss
    logical :: offset 
    integer :: inter(8)
    REAL(KIND=RP)           :: corners(3,8)
    REAL(KIND=RP)           :: points1(3,2,2)
    REAL(KIND=RP),allocatable   :: points2(:,:,:)
    real(kind=RP) ::  xxx(3, 0:1, 0:1)
    real(kind=RP) :: NODES(128,3)
    REAL(KIND=RP),allocatable   :: uNodes(:)
    REAL(KIND=RP),allocatable   :: vNodes(:)
    s=0.0_RP
    o=0.0_RP
    oo=0.0_RP
    ss=0.0_RP 
    XYZ=0.0_RP 
    XR=0.0_RP
    inter=0
    offset=.false. 
    nm=MOD(m,n)
    kkk=1
    points1=0.0_RP

    x=1.0_RP-nm*(2.0_RP/n)
    x=(-40.0_RP/(4.0_RP*DATAN(1.0_RP)))*theta + 1.0_RP
    write(*,*)'x=', x
    do i=1, oldnnode
        new_nodes(i) % X =self % nodes(i) % X 
        new_nodes(i) % globID =self % nodes(i) % globID 
       ! if (i .ne. new_nodes(i) % globID) write(*,*) 'we are at line 6131 of hex mesh...'
    end do 

    do i=1,size(self%elements)
      if (self%elements(i)%sliding) then 
         do j=1,8
            new_nodes(self%elements(i)%nodeIDs(j))%tbrotated=.true.
         end do 
      end if 
    end do 

    l=SIZE(self % nodes)+1
    !X
    ROT=0.0_RP
    ROT(1,1)=1.0_RP 
    ROT(2,2)=COS(theta)
    ROT(2,3)=-SIN(theta)
    ROT(3,2)=SIN(theta)
    ROT(3,3)=COS(theta)
    
    !Z
    ROT=0.0_RP
    ROT(1,1)=COS(theta)
    ROT(2,2)=COS(theta)
    ROT(3,3)=1.0_RP 
    ROT(1,2)=-SIN(theta)
    ROT(2,1)=SIN(theta)
    !Y
    ROT=0.0_RP
    ROT(1,1)=COS(-theta)
    ROT(2,2)=1.0_RP
    ROT(3,3)=COS(-theta) 
    ROT(1,3)=SIN(-theta)
    ROT(3,1)=-SIN(-theta)
    o(1)=(1.0_RP+x)/2.0_RP
    o(2)=-o(1)
    !oo(1)=-o(1)
    o(3)=(x-1.0_RP)/2.0_RP
    o(4)=-o(3)
    !oo(2)=-o(2)
    s(1)=1.0_RP-o(1)
    s(2)=s(1)
    !ss(1)=s(1)
    s(3)=x-o(2)
    s(4)=s(3)
    !ss(2)=s(2)

    o(1)=(x-1.0_rp)/2.0_rp !-0.5
    o(2)=(1.0_rp-x)/2.0_rp !0.5 
    o(3)=(1.0_rp+x)/2.0_rp !0.5
    o(4)=(-x-1.0_rp)/2.0_rp   !-0.5
    s(1)=o(1)+1.0_rp 
    s(2)=1.0_rp-o(2)
    s(3)=1.0_rp-o(3)
    s(4)=o(4)+1.0_rp 
    write(*,*)'offset o:',o
    write(*,*)'scale s',s
    !rotate_f ace_patchs 
    allocate (Xpatch(3,self%numBFacePoints,self%numBFacePoints))
    z=1
    do i=1, self % no_of_elements 
      if (self%elements(i)%sliding) then 
         if (self % elements(i) % SurfInfo % IsHex8 ) then 
            do zz=1,8
               corners(:,zz)=MATMUL(ROT, self % elements(i) % SurfInfo % corners(:,zz))
            end do !z
            self % elements(i) % SurfInfo % corners=corners
         else 
        do j=1, 6 
               
            if (allocated(self % elements(i) % SurfInfo % facePatches(j) % uKnots)) then 
               if (self % elements(i) % SurfInfo % facePatches(j) % noOfKnots(1)==2) then

                  points1= self % elements(i) % SurfInfo % facePatches(j) % points
                  if (self % elements(i) % sliding) then
                     Xflat(:,1,1)=MATMUL(ROT, points1(:,1,1))
                     Xflat(:,2,1)=MATMUL(ROT, points1(:,2,1))
                     Xflat(:,2,2)=MATMUL(ROT, points1(:,2,2))
                     Xflat(:,1,2)=MATMUL(ROT, points1(:,1,2))
                     self % elements(i) % SurfInfo % facePatches(j) % points=Xflat
                  end if 
               else 
                  allocate(uNodes(size(self % elements(i) % SurfInfo % facePatches(j) % uKnots)))
                  allocate(vNodes(size(self % elements(i) % SurfInfo % facePatches(j) % vKnots)))

                  uNodes=self % elements(i) % SurfInfo % facePatches(j) %uKnots
                  vNodes=self % elements(i) % SurfInfo % facePatches(j) %vKnots
                  allocate(points2(3,size(uNodes), size(vNodes)))
                  points2= self % elements(i) % SurfInfo % facePatches(j) % points
                     if (self % elements(i) % sliding) then
                     !rotate face patch curved
                        do k=1, numBFacePoints
                            do kk=1, numBFacePoints
                               Xpatch(:,kk,k)= MATMUL(ROT, points2(:,kk,k) )
                           end do 
                        end do 
                        call self % elements(i) % SurfInfo % facePatches(j) % destruct
                        call self % elements(i) % SurfInfo % facePatches(j) % construct(uNodes, vNodes, Xpatch)
                        
                     end if 
                     deallocate(points2)
                     deallocate(uNodes)
                     deallocate(vNodes)
               end if 
            z=z+1
            end if 
         end do!j
         end if 
      end if !self%elements(i)%sliding
    end do !i

    deallocate(Xpatch)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    l=oldnnode+1
    do i=1, size(arr2)
      if (.not. (self % elements(arr2(i)) % sliding)) write(*,*) "problem elem arr2"
      do j=1, 8
         if (self % elements(arr2(i)) % nodeIDs(j) .LE. oldnnode) then 
            if (self % elements(arr2(i)) % nodeIDs(j)==0) write(*,*) 'node ', j, 'of element', i, '=0 wtf'
            XR= self % nodes (self % elements(arr2(i)) % nodeIDs(j)) % X 
            XYZ(j,:)= MATMUL(ROT, XR)
         else 
            XYZ(j,:)=0.0_RP
         end if 
      end do 
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,1))) % X=XYZ(face_othernodes(i,1),:)
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,2))) % X=XYZ(face_othernodes(i,2),:)
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,3))) % X=XYZ(face_othernodes(i,3),:)
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,4))) % X=XYZ(face_othernodes(i,4),:)
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,1))) % rotated=.true.
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,2))) % rotated=.true.
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,3))) % rotated=.true.
      new_nodes(self % elements (arr2(i))% nodeIDs(face_othernodes(i,4))) % rotated=.true.
      inter=0
     ! do j=5, 8
      if (.not.  self%sliding) then 
      do j=1,4
         if (self % elements(arr2(i)) % nodeIDs(face_nodes(i,j)) .LE. oldnnode) then 
               if (self % elements(arr2(i)) % nodeIDs(face_nodes(i,j)) .GT. oldnnode) then 
                  write(*,*) 'problem logistic 0'
              end if 
               new_nodes(l) % globID=l
               new_nodes(l) % X = XYZ(face_nodes(i,j),:)
               new_nodes(l) % rotated=.true.
               inter(face_nodes(i,j))=l
               l=l+1
            do k=2,3
               do kkk=3,6
               eID=Connect(i,k,1)
               if ((eID .ne. 0)) then
               do ll=1,8
                  if ( (self % elements (arr2(i)) % nodeIDs(face_nodes(i,j))==self % elements (eID) % nodeIDs(ll)))then

                        self % elements (eID) % nodeIDs(ll)=inter(face_nodes(i,j))
                        if (inter(face_nodes(i,j))==0) write(*,*) 'element', eID, 'is receiving 0 for node', ll
                     
                  end if  
               end do !ll
            end if 
            end do !kk
            end do  !k
         end if 
      end do !j
         !do j=5,8
      do j=1,4
            if ((self % elements (arr2(i)) % nodeIDs(face_nodes(i,j))) .LE. oldnnode) then 
               if (inter(face_nodes(i,j))==0) write(*,*) 'element', arr2(i), 'is receiving 0 for node 5' ,self % elements (arr2(i)) % nodeIDs(face_nodes(i,j)) 
               self % elements (arr2(i)) % nodeIDs(face_nodes(i,j))=inter(face_nodes(i,j))
            end if 
         end do  !j
      end if 
   end do 
  
   do i=1, size(arr3)
      if (.not. (self % elements(arr3(i)) % sliding)) write(*,*) "problem elem arr2"
      do j=1, 8
            XR= self % nodes (self % elements(arr3(i)) % nodeIDs(j)) % X 
            XYZ(j,:)= MATMUL(ROT, XR)
      end do 
      do j=1,8
         if (.not.new_nodes(self % elements (arr3(i))% nodeIDs(j)) % rotated ) then 
      new_nodes(self % elements (arr3(i))% nodeIDs(j)) % X=XYZ(j,:)
      new_nodes(self % elements (arr3(i))% nodeIDs(j)) % rotated =.true.
         end if 
      end do 
   end do 

end subroutine HexMesh_RotateNodes


subroutine HexMesh_Modifymesh(self, nodes, nelm, arr1, arr2,arr3,Mat, center, o, s, face_nodes, face_othernodes,rotmortars, th,newnNodes,numBFacePoints,oldnnode, theta)
   IMPLICIT NONE 
   Class(HexMesh), intent(inout)    :: self 
   integer         , intent(in)    :: nodes
   integer, intent(in)    :: nelm
   integer, intent(inout) :: arr1(nelm)
   integer, intent(inout) :: arr2(nelm)
   integer, intent(inout) :: arr3(nelm)
   integer, intent(inout) :: Mat(nelm,12)
   real(kind=RP), intent(inout) :: center(2)
   real(kind=RP), intent(inout) :: o(4)
   real(kind=RP), intent(inout) :: s(4)
   integer, intent(inout) :: face_nodes(nelm,4)
   integer, intent(inout) :: face_othernodes(nelm,4)
   integer, intent(inout) :: rotmortars(nelm*2)
   real(kind=RP), intent(inout)  :: th 
   integer,intent(in) :: newnNodes
   integer,intent(inout) :: numBFacePoints
   integer,intent(inout) :: oldnnode
   real(kind=RP), intent(inout)   :: theta 

   type(Node), allocatable      :: new_nodes(:)
   !type(Node), allocatable  :: tmpNodes(:)
   integer :: l , i, j , new_nNodes, new_nFaces, n, m
   integer :: Connect(nelm, 9,6)
   !real(kind=RP) :: theta 
   real(kind=RP) :: PI
   logical :: success
   integer :: dealloc_status, sn, eID, fID

   allocate (new_nodes(newnNodes))
   arr1=0
   arr2=0
   new_nNodes=0
   new_nFaces=0
   !PI=4.0_RP*DATAN(1.0_RP)
   !theta=PI/3.0_RP
   !theta=PI/48.0_RP   !1/2
   !theta=PI/32.0_RP   !1/2
   !theta=PI/40.0_RP
   !theta=PI/60.0_RP
   !theta=PI/20.0_RP
   !theta=0.0_RP
   !theta=PI/20000.0_RP



   th=theta
   n=3
   m=1
   sn=SIZE(self % nodes)
   new_nNodes=SIZE(self % nodes) 
   o=0.0_RP
   s=0.0_RP
  
   if (.not.self%sliding) then 
      call self % MarkSlidingElementsRadius(1.01_RP,nelm, center, arr3, arr2, arr1, Mat, Connect, new_nFaces, face_nodes,face_othernodes, rotmortars)
   end if 
   new_nNodes=SIZE(self % nodes) 

   if (.not.self%sliding) then 
    call self % RotateNodes(theta,nelm, n, m , new_nNodes, new_nodes, arr1, arr2, arr3, Connect, o, s, face_nodes,face_othernodes, numBFacePoints,oldnnode)
   else 
      call self % RotateNodes(theta,nelm, n, m , new_nNodes, new_nodes, self%arr1, self%arr2, self%arr3, self%Connect, o, s, self%face_nodes,self%face_othernodes, numBFacePoints,oldnnode)
   end if 
   if (.not. self%sliding) then 
      do i=1, size(self % Nodes)
      self % nodes(i) % X =new_nodes(i) % X
      self % nodes(i) % GlobID=new_nodes(i) % GlobID
      end do 
   end if 


   do l=1, SIZE(self%faces)
       !call self % faces(l) % Destruct
            self % faces(l) % ID = -1
            self % faces(l) % FaceType = HMESH_NONE
            self % faces(l) % rotation = 0
            self % faces(l) % NelLeft = -1
            self % faces(l) % NelRight = -1       
            self % faces(l) % NfLeft = -1       
            self % faces(l) % NfRight = -1       
            self % faces(l) % Nf = -1         
            self % faces(l) % nodeIDs = -1             
            self % faces(l) % elementIDs = -1
            self % faces(l) % elementSide = -1
            self % faces(l) % projectionType = -1
            self % faces(l) % boundaryName = ""
   end do 
   if (allocated(self%mortar_faces)) then 
      do l=1, SIZE(self%mortar_faces)
         !call self % faces(l) % Destruct
              self % mortar_faces(l) % ID = -1
              self % mortar_faces(l) % FaceType = HMESH_NONE
              self % mortar_faces(l) % rotation = 0
              self % mortar_faces(l) % NelLeft = -1
              self % mortar_faces(l) % NelRight = -1       
              self % mortar_faces(l) % NfLeft = -1       
              self % mortar_faces(l) % NfRight = -1       
              self % mortar_faces(l) % Nf = -1         
              self % mortar_faces(l) % nodeIDs = -1             
              self % mortar_faces(l) % elementIDs = -1
              self % mortar_faces(l) % elementSide = -1
              self % mortar_faces(l) % projectionType = -1
              self % mortar_faces(l) % boundaryName = ""
     end do 
   end if 
   deallocate(new_nodes)
end subroutine HexMesh_Modifymesh

!*********************************************************************************************************
subroutine HexMesh_ConstructSlidingMortarsConforming(self, nodes, nelm, array1, array2,Mat, o, s, mortararr2,rotmortars)
   Class(HexMesh), intent(inout) :: self
    integer         , intent(in)    :: nodes
    integer,      intent(in)    :: nelm
    integer, intent(in)          :: array1(nelm)
    integer, intent(in)          :: array2(nelm)
    integer, intent(in)          :: Mat(nelm,12)
    real(kind=RP), intent(in)              :: o(4)
    real(kind=RP), intent(in)               :: s(4)
    integer, intent(in)            :: mortararr2(nelm, 2)
    integer, intent(in)            :: rotmortars(2*nelm)
   integer :: i, j, k, l ,m , fID1, fID2, eID1, eID2, faceNumber, faceNumber2
   integer :: faceNodeIDs(4)
   integer :: faceNodeIDs2(4)
   integer :: faceNodeIDs1(4)
   integer :: nodeIDs(8)
   integer  :: NelL(2), NelR(2)
   real(kind=RP) :: NODESS(4,3)
   real(kind=RP) :: NODES2(4,3)
   integer :: nodeIDs1(8)
   integer :: nodeIDs2(8)
   integer :: rot 
   integer :: el(nelm,2)
   do i=1,nelm 
      do l=1,nelm 
         if (array1(i)==array2(l)) write(*,*) 'tnaket array1(i)==array2(l)', i,l, array1(i),array2(l)
      end do 
   end do

   ALLOCATE (self % mortar_faces(nelm))
   l=1
   i=1
   do while (i .LE. nelm)
       eID1=Mat(i,3)
       nodeIDs = self % elements(eID1) % nodeIDs
       faceNumber=Mat(i,6)
       !write(*,*) 'mine 7521 of hexmesh'n Mat(i,:)
       DO j = 1, 4
           faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
       END DO

      fID1=self %elements(eID1) % faceIDs(Mat(i,6))!!!!!!!!!!!

             allocate(self % faces(fID1) % Mortar(1)) 
      ! write(*,*) 'computing the to mortars of element:', eID1
       !do m=1, 2
           CALL self % mortar_faces(l) % Construct(ID  = l, &
           nodeIDs = faceNodeIDs, &
           elementID = eID1,       &
           side = faceNumber)
           self%mortar_faces(l)=self % faces(fID1)
           !self%mortar_faces(l)=self % faces(self %elements(Mat(i,1) ) % faceIDs(Mat(i,4)))
           self%mortar_faces(l)%Mortarpos=1
           self%mortar_faces(l)%ID=l
         !write(*,*) 'ID of the mortar', m,'ID=', l
           self % mortar_faces(l) % FaceType       = HMESH_INTERIOR
           !!!!!!!!!self % mortar_faces(l) % boundaryName = self % elements(eID1) % boundaryName(faceNumber)
           allocate(self % mortar_faces(l) % Mortar(2)) 
           self % mortar_faces(l) % Mortar(1)=fID1
           !self % faces(fID1) %IsMortar=3
           self % faces(fID1) % Mortar(1)=l
          ! if (m==2)  self % faces(fID1) % Mortar(2)=l

           !if (m==1) eID2= Mat(i,2)   !2
           eID2= Mat(i,1)   !3
         write(*,*) 'second element of the mortar',l,'eID', eID2
           self % mortar_faces(l) % elementIDs(2)  = eID2
           self % mortar_faces(l) % elementSide(2) = Mat(i,4)
           self % mortar_faces(l) % FaceType       = HMESH_INTERIOR
    
           self%elements(eID1)%faceSide(Mat(i,6))=1

            self%elements(eID2)%faceSide(Mat(i,4))=2
            fID2=self %elements(eID2) % faceIDs(Mat(i,4))
           self % mortar_faces(l) % Mortar(2)=fID2
            self % mortar_faces(l) % rotation=Mat(i,7)
   
           if (.not.allocated(self % faces(fID2) % Mortar)) allocate(self % faces(fID2) % Mortar(1)) 
  
                self % mortar_faces(l)%offset(1)=0.0_RP
                self % mortar_faces(l)%offset(2)=0.0_RP
               self % mortar_faces(l)%s(1)=1.0_RP
               self % mortar_faces(l)%s(2)= 1.0_RP
               self % faces(fID2) % Mortar(1)=l

           !end if  
           !self % faces(fID2) % Mortar(m)=l
           l=l+1
      !end do     
       self % elements(array2(i))% MortarFaces(1)=3 
       self % elements(array1(i))% MortarFaces(2)=3                                                        
        i=i+1
   end do

   do i=1, size(self%mortar_faces)
       associate(f => self % mortar_faces(i))
         associate(eL => self % elements(f % elementIDs(1)), &
            eR => self % elements(f % elementIDs(2))   )
            NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
            NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
            call f % LinkWithElements(NelL, NelR, nodes, f%offset, f%s)
           end associate
    end associate
   end do
   

   do i=1, size(self%mortar_faces)
      associate(f => self % mortar_faces(i))
        associate(eL => self % elements(f % elementIDs(1)), &
           eR => self % elements(f % elementIDs(2))   )
           NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
           NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
           call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                   NodalStorage(f % Nf), NodalStorage(eL % Nxyz), &
                                     eL % geom, eL % hexMap, f % elementSide(1), &
                                     f % projectionType(1), 1, 0,.true.,1,f%s(1),f%IsMortar)

          end associate
          f % geom % h = minval(self % elements(f % elementIDs(2)) % geom % jacobian) &
          / maxval(f % geom % jacobian)
   end associate
  end do 
  do i=1, size(self%mortar_faces)
   self%mortar_faces(i)%geom%x=self%faces(self%mortar_faces(i)%Mortar(1))%geom%x
   self%mortar_faces(i)%geom%GradXi=self%faces(self%mortar_faces(i)%Mortar(1))%geom%GradXi
   self%mortar_faces(i)%geom%GradEta=self%faces(self%mortar_faces(i)%Mortar(1))%geom%GradEta
   self%mortar_faces(i)%geom%GradZeta=self%faces(self%mortar_faces(i)%Mortar(1))%geom%GradZeta
   self%mortar_faces(i)%geom%normal=self%faces(self%mortar_faces(i)%Mortar(1))%geom%normal
   self%mortar_faces(i)%geom%t1=self%faces(self%mortar_faces(i)%Mortar(1))%geom%t1 
   self%mortar_faces(i)%geom%t2=self%faces(self%mortar_faces(i)%Mortar(1))%geom%t2
   self%mortar_faces(i)%geom%surface=self%faces(self%mortar_faces(i)%Mortar(1))%geom%surface
   self%mortar_faces(i)%geom%jacobian=self%faces(self%mortar_faces(i)%Mortar(1))%geom%jacobian
   self%mortar_faces(i)%geom%h=self%faces(self%mortar_faces(i)%Mortar(1))%geom%h
   self%faces(self%mortar_faces(i)%Mortar(2))%elementIDs(2)=   self%faces(self%mortar_faces(i)%Mortar(2))%elementIDs(1)
   self%faces(self%mortar_faces(i)%Mortar(2))%elementIDs(1)=0
   write(*,*) 'mortar number', i, 'elements:', self%mortar_faces(i)%elementIDs
   write(*,*) 'mortar number', i, 'elemenfaces:', self%mortar_faces(i)%Mortar 
   write(*,*) 'element of first face:', self%faces(self%mortar_faces(i)%Mortar(1))%elementIDs
   write(*,*) 'element of second face:', self%faces(self%mortar_faces(i)%Mortar(2))%elementIDs

  end do 

  do i=1, size(array2)
   nodeIDs1=self%elements(Mat(i,3))% nodeIDs    !3
   nodeIDs2=self%elements(Mat(i,1))% nodeIDs    !1
   faceNumber2=Mat(i,6)          !rot of 3
   faceNumber=Mat(i,4)           !rot of 1
   DO j = 1, 4
       faceNodeIDs2(j) = nodeIDs1(localFaceNode(j,faceNumber2))
       NODESS(j,:)=self%nodes(faceNodeIDs2(j))%X
       faceNodeIDs(j) = nodeIDs2(localFaceNode(j,faceNumber))
       NODES2(j,:)=self%nodes(faceNodeIDs(j))%X
   END DO        
   rot=faceRotationnodes(masterNodeIDs=NODESS, slaveNodeIDs=NODES2)
   write(*,*) 'i=',i
   if (rot .NE. Mat(i,7)) write(*,*) 'problem rot .NE. MAT:, rot=', rot, 'MAT=', Mat(i,7)
   if (rot .EQ. Mat(i,7)) write(*,*) 'no problem', rot, 'MAT=', Mat(i,7)
   self % mortar_faces(i) % rotation=rot

  end do 

  do i=1,size(self%mortar_faces)
   write(*,*) 'mortar',i,'between element:',self%mortar_faces(i)%elementIDs
   write(*,*)'rotation:', self%mortar_faces(i)%rotation
end do 

end subroutine HexMesh_ConstructSlidingMortarsConforming

 subroutine HexMesh_ConstructMortars(self, nodes, nelm, array1, array2,Mat, o, s, mortararr2,rotmortars, th,confor)
   Class(HexMesh), intent(inout) :: self
    integer         , intent(in)    :: nodes
    integer,      intent(in)    :: nelm
    integer, intent(in)          :: array1(nelm)
    integer, intent(in)          :: array2(nelm)
    integer, intent(in)          :: Mat(nelm,12)
    real(kind=RP), intent(in)              :: o(4)
    real(kind=RP), intent(in)               :: s(4)
    integer, intent(in)            :: mortararr2(nelm, 2)
    integer, intent(in)            :: rotmortars(2*nelm)
   real(kind=RP), intent(in)               :: th
   logical, intent(in)             :: confor

   integer :: i, j, k, l ,m , fID1, fID2, eID1, eID2, faceNumber, faceNumber1,faceNumber2
   integer :: faceNodeIDs(4)
   integer :: nodeIDs(8)
   integer  :: NelL(2), NelR(2)
   real(kind=RP)              :: theta
   real(kind=RP)              :: XR(4,3)
   integer :: faceNodeIDs2(4)
   integer :: faceNodeIDs1(4)
   real(kind=RP) :: NODESS(4,3)
   real(kind=RP) :: NODES2(4,3)
   integer :: nodeIDs1(8)
   integer :: nodeIDs2(8)
   real(KIND=RP) :: ROT(3,3)
   integer :: rota 
   integer :: el(nelm,2)
   integer :: ell(nelm,2)
   integer :: nf 
  do i=1,nelm
     do l=1,nelm 
        if (array1(i)==array2(l)) write(*,*) ' array1(i)==array2(l)', i,l, array1(i),array2(l)
     end do 
  end do

  if (.not. allocated(self%mortar_faces)) ALLOCATE (self % mortar_faces(2*nelm))
  l=1
  i=1
  do i=1,size(array2)
      if (.not.confor) then 
         eID1=Mat(i,3)
         faceNumber=Mat(i,6)
      else 
         eID1=Mat(i,10)
         faceNumber=Mat(i,11)
      end if 
      nodeIDs = self % elements(eID1) % nodeIDs
      DO j = 1, 4
          faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
      END DO

     fID1=self %elements(eID1) % faceIDs(faceNumber)!!!!!!!!!!!

     if (.not. allocated(self % faces(fID1) %  Mortar)) allocate(self % faces(fID1) % Mortar(2)) 

          CALL self % mortar_faces(l) % Construct(ID  = l, &
          nodeIDs = faceNodeIDs, &
          elementID = eID1,       &
          side = faceNumber)
          !self%mortar_faces(l)=self % faces(fID1)

         ! self%mortar_faces(l)%ID=l
          self%mortar_faces(l)%Mortarpos=0
          self % mortar_faces(l) % FaceType       = HMESH_INTERIOR
          if (.not. allocated(self % mortar_faces(l) % Mortar))  allocate(self % mortar_faces(l) % Mortar(2)) 
          self % mortar_faces(l) % Mortar(1)=fID1
          self % faces(fID1) % Mortar(1)=l

          eID2= Mat(i,2)   

          self % mortar_faces(l) % elementIDs(2)  = eID2
   
          self % mortar_faces(l) % elementSide(2) = Mat(i,5)
          self % mortar_faces(l) % FaceType       = HMESH_INTERIOR

          self%elements(eID1)%faceSide(faceNumber)=1

           self%elements(eID2)%faceSide(Mat(i,5))=2

           fID2=self %elements(eID2) % faceIDs(Mat(i,5))

          self % mortar_faces(l) % Mortar(2)=fID2

           self % mortar_faces(l) % rotation=Mat(i,8)

          if (.not.allocated(self % faces(fID2) % Mortar)) allocate(self % faces(fID2) % Mortar(2)) 
               self % mortar_faces(l)%offset(1)=o(1)
               self % mortar_faces(l)%offset(2)=o(2)

              self % mortar_faces(l)%s(1)=s(1)
              self % mortar_faces(l)%s(2)= s(2)
              if (confor) then 
               self % mortar_faces(l)%offset(1)=0.0_RP
               self % mortar_faces(l)%offset(2)=0.0_RP
   
              self % mortar_faces(l)%s(1)=0.0_RP
              self % mortar_faces(l)%s(2)=0.0_RP
              end if 
              self % faces(fID2) % Mortar(2)=l

          l=l+1    
      self % elements(array2(i))% MortarFaces(1)=3 
      self % elements(array1(i))% MortarFaces(2)=3                                                        
       !i=i+1
  end do


l=nelm+1
i=1
do i=1,size(array2)
   if (.not.confor) then 
      eID1=Mat(i,3)
      faceNumber=Mat(i,6)
   else 
      eID1=Mat(i,10)
      faceNumber=Mat(i,11)
   end if 
   nodeIDs = self % elements(eID1) % nodeIDs
   DO j = 1, 4
       faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
   END DO

  fID1=self %elements(eID1) % faceIDs(faceNumber)!!!!!!!!!!!



       CALL self % mortar_faces(l) % Construct(ID  = l, &
       nodeIDs = faceNodeIDs, &
       elementID = eID1,       &
       side = faceNumber)
      ! self%mortar_faces(l)=self % faces(fID1)

      ! self%mortar_faces(l)%ID=l
       self%mortar_faces(l)%Mortarpos=1
       self % mortar_faces(l) % FaceType       = HMESH_INTERIOR
       if (.not.allocated(self % mortar_faces(l) % Mortar)) allocate(self % mortar_faces(l) % Mortar(2)) 
       self % mortar_faces(l) % Mortar(1)=fID1
       self % faces(fID1) % Mortar(2)=l

       eID2= Mat(i,1)   !3

       self % mortar_faces(l) % elementIDs(2)  = eID2

       self % mortar_faces(l) % elementSide(2) = Mat(i,4)
       self % mortar_faces(l) % FaceType       = HMESH_INTERIOR

       self%elements(eID1)%faceSide(faceNumber)=1

        self%elements(eID2)%faceSide(Mat(i,4))=2

        fID2=self %elements(eID2) % faceIDs(Mat(i,4))

       self % mortar_faces(l) % Mortar(2)=fID2
       self % mortar_faces(l) % rotation=Mat(i,8)
       if (.not.allocated(self % faces(fID2) % Mortar)) then 
         !write(*,*) 'line 7622 not allocated'
         if (.not.allocated(self % faces(fID2) % Mortar)) allocate(self % faces(fID2) % Mortar(2)) 
       end if 

            self % mortar_faces(l)%offset(1)=o(3)
            self % mortar_faces(l)%offset(2)=o(4)

           self % mortar_faces(l)%s(1)=s(3)
           self % mortar_faces(l)%s(2)= s(4)
           if (confor) then 
            self % mortar_faces(l)%offset(1)=0.0_RP
            self % mortar_faces(l)%offset(2)=0.0_RP

           self % mortar_faces(l)%s(1)=1.0_RP
           self % mortar_faces(l)%s(2)=1.0_RP
           end if 
           self % faces(fID2) % Mortar(1)=l

       l=l+1    
   self % elements(array2(i))% MortarFaces(1)=3 
   self % elements(array1(i))% MortarFaces(2)=3                                                        
    !i=i+1
end do

do i=1, size(self%mortar_faces)
   associate(f => self % mortar_faces(i))
     associate(eL => self % elements(f % elementIDs(1)), &
        eR => self % elements(f % elementIDs(2))   )
        NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
        NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
        call f % LinkWithElements(NelL, NelR, nodes, f%offset, f%s)
       end associate
end associate
end do


do i=1, size(self%mortar_faces)
  associate(f => self % mortar_faces(i))
    associate(eL => self % elements(f % elementIDs(1)), &
       eR => self % elements(f % elementIDs(2))   )
       NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
       NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
       call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                               NodalStorage(f % Nf), NodalStorage(eL % Nxyz), &
                                 eL % geom, eL % hexMap, f % elementSide(1), &
                                 f % projectionType(1), 1, 0,.true.,f%Mortarpos, f%s(1))


      end associate
      f % geom % h = minval(self % elements(f % elementIDs(2)) % geom % jacobian) &
      / maxval(f % geom % jacobian)
      nf=f % Nf(1)
end associate
end do 
l=nelm+1
do i=1,size(array2)
theta=(4.0_RP*DATAN(1.0_RP))/16.0_RP -th 
ROT=0.0_RP
ROT(1,1)=COS(theta)
ROT(2,2)=COS(theta)
ROT(3,3)=1.0_RP 
ROT(1,2)=-SIN(theta)
ROT(2,1)=SIN(theta)
   nodeIDs1=self%elements(Mat(i,3))% nodeIDs    !3
   nodeIDs2=self%elements(Mat(i,1))% nodeIDs    !1
   faceNumber2=Mat(i,6)          !facenumber of 3
   faceNumber1=Mat(i,4)           !facenumber of 1
   DO j = 1, 4
       faceNodeIDs2(j) = nodeIDs1(localFaceNode(j,faceNumber2))
       NODESS(j,:)=self%nodes(faceNodeIDs2(j))%X
       faceNodeIDs1(j) = nodeIDs2(localFaceNode(j,faceNumber1))
       NODES2(j,:)=self%nodes(faceNodeIDs1(j))%X
       XR(j,:)=NODES2 (j,:)
       NODES2(j,:)= MATMUL(ROT, XR(j,:))
   END DO        
   !rota=faceRotationnodes(masterNodeIDs=NODESS, slaveNodeIDs=NODES2)
   !write(*,*) 'i=,',i
   !write(*,*) 'l=,',l
   !if (rota .NE. Mat(i,7)) write(*,*) 'problem rot .NE. MAT:, rot=', rota, 'MAT=', Mat(i,7)
   !if (rota .EQ. Mat(i,7)) write(*,*) 'no problem', rota, 'MAT=', Mat(i,7)
   !  self % mortar_faces(l) % rotation=rota

   l=l+1
end do 




end subroutine HexMesh_ConstructMortars
END MODULE HexMeshClass
