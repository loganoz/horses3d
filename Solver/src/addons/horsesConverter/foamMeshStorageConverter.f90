!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!   HORSES3D to Foam Result - foamMeshStorage Module
!
!      This module convert Horses3D mesh storage to foamMesh Storage
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE foamMeshStorageConverter
    USE SMConstants
    USE Storage
    IMPLICIT NONE
	
	
    !
    TYPE foamFaces_t
        INTEGER                     :: facePoints    (4)
        INTEGER                     :: faceOwner
        INTEGER                     :: faceNeighbour
        INTEGER                     :: faceBoundary
	    LOGICAL                     :: master
    END TYPE foamFaces_t
    !
    TYPE foamBoundaries_t
        INTEGER                     :: faceStart
        INTEGER                     :: nFace
        TYPE(foamFaces_t)           ,ALLOCATABLE    :: faceBoundaries(:)
    END TYPE foamBoundaries_t
	!
	TYPE horsesPoints_t
		INTEGER                     :: ID_OLD			  ! Old Point ID 		- start from 0
		INTEGER						:: ID_NEW			  ! New Point ID for OF - start from 0
		INTEGER						:: copyID      		  ! Old ID of points with identical location
		INTEGER						:: eID				  ! ID of element this point represented
		INTEGER						:: i				  ! i N(1) of ID of element this point represented
		INTEGER						:: j				  ! j N(2) of ID of element this point represented
		INTEGER						:: k				  ! k N(3) of ID of element this point represented
		LOGICAL                     :: master			  ! if .true. than this is the master
		LOGICAL                     :: edgeline
		LOGICAL                     :: boundary
		LOGICAL                     :: internal
    END TYPE horsesPoints_t
    !
    TYPE foamMesh_t
        INTEGER                     :: nCells
        INTEGER                     :: Nout(3)
        INTEGER                     :: nPoints
		INTEGER                     :: nPointsMultiple
		INTEGER                     :: nPointsUnique
        INTEGER                     :: nFacesShared
        INTEGER                     :: nFacesUnshared
		INTEGER                     :: nMultipleFaces
        INTEGER                     :: nFacesBoundaries
		TYPE(horsesPoints_t)        ,ALLOCATABLE    :: horsesPoints(:) 
        TYPE(foamFaces_t)           ,ALLOCATABLE    :: faceUnshared(:)
        TYPE(foamFaces_t)           ,ALLOCATABLE    :: faceShared(:)
        TYPE(foamBoundaries_t)      ,ALLOCATABLE    :: boundaries(:)
        CONTAINS
            PROCEDURE   :: Construct     => constructMeshFaces
    END TYPE foamMesh_t
	!

!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine construct faces for foam Mesh format with derived type: foamMesh_t with input type Mesh_t
!
        SUBROUTINE constructMeshFaces(self, mesh)
            IMPLICIT NONE
            CLASS(foamMesh_t)                        :: self
            TYPE(Mesh_t),       INTENT(INOUT)        :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: eID, eID1, eIDf  ! Spectral (not finite) element counter
            INTEGER         :: i,j,k,l,i2,j2,k2
            INTEGER         :: pointID, pointID1, pointIDs          ! Finite Points counter
            INTEGER         :: N(3)
            INTEGER         :: faceUnsharedCount, faceSharedCount, faceBoundaryCount
            INTEGER         :: nFacesBoundary, nFacesBoundaries, nBoundaryEdgePoints, nElementCopy
			INTEGER         :: nCounter, nMultiplePoints, nMultiplePointsEdge
			LOGICAL         :: first
			REAL(KIND=RP)   :: x(3)
			REAL(KIND=RP)   :: tol
            INTEGER , ALLOCATABLE, DIMENSION(:,:,:)    :: CellPoints
            INTEGER , ALLOCATABLE, DIMENSION(:,:)      :: boundarySet
            INTEGER , ALLOCATABLE, DIMENSION(:)        :: boundaryFace, multiplePoints, multiplePointsEdge

			write(STD_OUT,'(15X,A,A)') ""
			write(STD_OUT,'(15X,A,A)') "Construct foamMesh Storage:"
			write(STD_OUT,'(15X,A,A)') "--------------------------"
			write(STD_OUT,'(20X,A,A)') "->  ", "Start Construct ..."
			
			tol =1.0e-10_RP

            self % Nout=mesh % elements(1) % Nout
            N = self % Nout
            self % nCells   = mesh % no_of_elements * ((N(1))*(N(2))*(N(3)))
            self % nPoints  = mesh % no_of_elements * ((N(1)+1)*(N(2)+1)*(N(3)+1))
			self % nPointsUnique = self % nPoints
            self % nFacesShared = mesh % no_of_elements * ((N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))

            faceSharedCount  =0
            faceUnsharedCount=0
            faceBoundaryCount=0
            nFacesBoundaries =0
!
!        Calculate the number of FE Faces on each Boundary -- Tested and Reliable
!        ------------------------------------------------------------------
           ASSOCIATE ( e => mesh % boundaries)
           ALLOCATE( self % boundaries(size(e)))
           DO i=1, size(e)
                nFacesBoundary=0
                DO j=1, e(i) % no_of_faces
                    SELECT CASE(e(i) % elementSides(j))
                        CASE(1,2)
                            nFacesBoundary = nFacesBoundary + N(1)*N(3)
                        CASE(3,5)
                            nFacesBoundary = nFacesBoundary + N(1)*N(2)
                        CASE(4,6)
                            nFacesBoundary = nFacesBoundary + N(2)*N(3)
                    END SELECT
					nElementCopy = 0
					DO k=1, size(e)
						IF (k .eq. i ) CYCLE ! ONLY LOOK AT DIFFERENT ZONE
						DO l=1, e(k) % no_of_faces 
							IF (e(i) % elements(j) .eq. e(k) % elements (l)) then
								nElementCopy = nElementCopy+1
								EXIT
							END IF 
						END DO 
						IF (nElementCopy .eq. 3) THEN
							EXIT
						END IF
					END DO
					e(i) % edgeDomain (j)   = .false.
					e(i) % cornerDomain (j) = .false.
					IF(nElementCopy.ge.2) e(i) % edgeDomain (j)   = .true.
					IF(nElementCopy.ge.3) e(i) % cornerDomain (j) = .true.
                 END DO
                 nFacesBoundaries=nFacesBoundaries+nFacesBoundary
                 ALLOCATE(self % boundaries(i) % faceBoundaries (nFacesBoundary))
            END DO
            END ASSOCIATE
            self % nFacesBoundaries = nFacesBoundaries
            self % boundaries(:) % nFace = 0
!
!        Store the number of Faces Unshared without boundary - Tested
!        --------------------------------------------
             self % nFacesUnshared   = (mesh % no_of_elements *2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1)))-nFacesBoundaries
!
!        Allocate point, faceUnshared and face shared
!        --------------------------------------------
            ALLOCATE (self % faceUnshared (1:self % nFacesUnshared),self % faceShared (1:self % nFacesShared))
            ALLOCATE(CellPoints(0:N(1),0:N(2),0:N(3)))
!
!        Construct Spectral Boundaries Array -- Tested and Reliable
!        ----------------------------------------------------------
            CALL constructSpectralBoundariesArray(mesh, boundarySet)
!
!        Assign FE faces for shared and unshared without boundary  to foamMesh; Store points into horsesPoints
!        -----------------------------------------------------------------------------------------------------
			write(STD_OUT,'(20X,A,A)') "->  ", "Construct all points storage ..."
			ALLOCATE (self % horsesPoints (self % nPoints))
			nMultiplePoints=0
			nMultiplePointsEdge=0
            DO eID = 1 , size(mesh % elements)

                pointID=((N(1)+1)*(N(2)+1)*(N(3)+1))*(eID-1)           ! FE Point ID after eID
!
!               Assign points at Spectral Element eID for FE faces extraction and store Point into horsesPoints
!               -----------------------------------------------------------------------------------------------
                DO k = 0, N(3) ; DO j = 0, N(2) ; DO i = 0, N(1)
                    CellPoints(i,j,k)=pointID
                    pointID=pointID+1
					ASSOCIATE(p =>self % horsesPoints (pointID))
					p % ID_OLD = pointID-1
					p % eID    = eID
					p % i      = i 
					p % j      = j
					p % k      = k 
					p % master = .true.            ! Update later if .false.
					p % copyID = pointID
					p % edgeline = .false.
					p % boundary = .false.
					p % internal = .true.
					if (((k.eq.0).or.(k.eq.N(3))).and.((j.eq.0).or.(j.eq.N(2))))  p % edgeline = .true.
					if (((k.eq.0).or.(k.eq.N(3))).and.((i.eq.0).or.(i.eq.N(1))))  p % edgeline = .true.
					if (((j.eq.0).or.(j.eq.N(2))).and.((i.eq.0).or.(i.eq.N(1))))  p % edgeline = .true.
					if ((k.eq.0).or.(k.eq.N(3)).or.(j.eq.0).or.(j.eq.N(2)).or.(i.eq.0).or.(i.eq.N(1)))  p % internal = .false.
					
					if ((.not. p % internal).and.(.not.p%edgeline)) nMultiplePoints=nMultiplePoints+1
					if (p % edgeline) nMultiplePointsEdge = nMultiplePointsEdge+1

					END ASSOCIATE
                END DO                ; END DO                ; END DO
!
!               Assign Boundary Faces from single Spectral Element eID into Finite Elements -- Tested
!               -------------------------------------------------------------------------------------------
                CALL finiteElementBoundariesFace (N, eID, size(boundarySet,2), boundarySet, boundaryFace)
!
!               Assign Faces from single Spectral Element eID into Finite Elements -- Tested
!               -------------------------------------------------------------------------------------------
                CALL finiteElementFacesSet(CellPoints, N, eID, boundaryFace, faceSharedCount, faceUnsharedCount, self)
            END DO
!
!           Assign boundary faces of all FE cell into foamMesh based on its boundary type -- Tested
!           ---------------------------------------------------------------------------------------
			write(STD_OUT,'(20X,A,A)') "->  ", "Construct boundaries ..."
            DO i=1,size(self % boundaries)
                IF (i.EQ.1) THEN
                    self % boundaries(i) % faceStart = size(self % faceShared) + size(self % faceUnshared)
                ELSE
                    self % boundaries(i) % faceStart = self % boundaries(i-1) % faceStart +self % boundaries(i-1) % nFace
                END IF
				DO j=1, self % boundaries(i) % nFace
					DO k=1,4
						self % horsesPoints (self % boundaries(i) % faceBoundaries(j) % facePoints (k) +1) % boundary = .true.
					END DO 
				END DO 
            END DO
!
!           Create a simple list of MultiplePoints
!           ------------------------------------------------	
			ALLOCATE(multiplePoints(nMultiplePoints))
			ALLOCATE(multiplePointsEdge(nMultiplePointsEdge))
			
			multiplePoints=0
			nMultiplePoints=0
			multiplePointsEdge=0
			nMultiplePointsEdge=0
			DO i=1, self % nPoints
				ASSOCIATE(p =>self % horsesPoints (i))
				if ((.not.p % internal).and.(.not.p % boundary).and.(.not.p % edgeline)) then
					nMultiplePoints=nMultiplePoints+1
					multiplePoints(nMultiplePoints)=i
				end if 
				if (p%edgeline) then
					nMultiplePointsEdge=nMultiplePointsEdge+1
					multiplePointsEdge(nMultiplePointsEdge)=i
				end if 
				END ASSOCIATE
			END DO 
!
!           Correct the number point with identical location
!           ------------------------------------------------
			write(STD_OUT,'(20X,A,A)') "->  ", "Looping to correct number point with identical location ..."
			l = 0
!$omp parallel shared(mesh, self, tol, l, multiplePoints, nMultiplePoints, multiplePointsEdge, nMultiplePointsEdge )
!$omp do schedule(runtime) private(k,i2,k2,x)	

			DO i=1, nMultiplePoints
!$omp critical
				l = l + 1
				IF (mod (l,int(nMultiplePoints/5)).eq.0) then
					write(STD_OUT,'(25X,A,A,I10,A,I10,A)') "->  ","Looping Multiple Points: ", l," of ", nMultiplePoints
				END IF
!$omp end critical
				
				k=multiplePoints(i)
				ASSOCIATE(p =>self % horsesPoints (k))
				if (.not. p%master) CYCLE
				x = mesh % elements (p % eID) % xOut (:,p % i,p % j,p % k)
				DO i2=1, nMultiplePoints-1
					k2 = i+INT((-1_RP)**(i2+1_RP)*CEILING(real(i2)/2_RP))
					if (k2.le.0) k2=k2+nMultiplePoints
					if (k2.gt.nMultiplePoints) k2=k2-nMultiplePoints
					k2=multiplePoints(k2)
					if (k2.eq.k) CYCLE
					ASSOCIATE(p2 =>self % horsesPoints (k2))
					if (.not. p2%master) CYCLE
					if (maxval(abs(mesh % elements (p2 % eID) % xOut (:,p2 % i,p2 % j,p2 % k)-x)).le.tol) then
						if (k2.lt.k) EXIT
						p2 % copyID = k
						p2 % master = .false.
						EXIT
					end if
					END ASSOCIATE
				END DO  
				END ASSOCIATE
			END DO 
!$omp end do
			
			l = 0
!$omp do schedule(runtime) private(k,i2,k2,x)	
			DO i=1, nMultiplePointsEdge
!$omp critical 
				l = l + 1
				IF (mod(l,int(nMultiplePointsEdge/5)).eq.0) then
					write(STD_OUT,'(25X,A,A,I10,A,I10,A)') "->  ","Looping Multiple Points Edge: ", l," of ", nMultiplePointsEdge
				END IF
!$omp end critical
				
				k=multiplePointsEdge(i)
				ASSOCIATE(p =>self % horsesPoints (k))
				if (.not. p%master) CYCLE
				x = mesh % elements (p % eID) % xOut (:,p % i,p % j,p % k)
				DO i2=1, int(nMultiplePointsEdge)
					k2=i-i2
					if ((k2.ge.1).and.(i2.gt.int(nMultiplePointsEdge/2))) k2=i+i2-int(nMultiplePointsEdge/2)
					if (k2.lt.1) k2=i-(i-i2)
					if (k2.gt.nMultiplePointsEdge) EXIT
					k2=multiplePointsEdge(k2)
					if (k2.eq.k) CYCLE
					ASSOCIATE(p2 =>self % horsesPoints (k2))
					if (.not. p2%master) CYCLE
					if (maxval(abs(mesh % elements (p2 % eID) % xOut (:,p2 % i,p2 % j,p2 % k)-x)).le.tol) then
						if (k2.lt.k) EXIT
						p2 % copyID = k
						p2 % master = .false.
					end if
					END ASSOCIATE
				END DO 
				END ASSOCIATE
			END DO 
!$omp end do
!$omp end parallel
			DEALLOCATE(multiplePoints, multiplePointsEdge)

			CALL updatePointID(self, mesh)
		
			CALL detectMultipleFaces(self)

			
			write(STD_OUT,'(20X,A,A)') "->  ", "Finish constructing foamMesh storage"
			write(STD_OUT,'(15X,A,A)') ""
			

        END SUBROUTINE constructMeshFaces
		
		subroutine updatePointID(self, mesh)
            IMPLICIT NONE
			CLASS(foamMesh_t)                                ,INTENT(INOUT):: self
			TYPE(Mesh_t),       INTENT(INOUT)       :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: eID, eID1, eIDf  ! Spectral (not finite) element counter
            INTEGER         :: i,j,k,l
            INTEGER         :: pointID       ! Finite Points counter
			INTEGER         :: nCounter
			REAL(KIND=RP)   :: x(3), tol
			
			tol = 1.0e-10_RP

!
!           Detect CopyID of multiple points, master/slave, and assign the associated new ID
!           --------------------------------------------------------------------------------
			write(STD_OUT,'(20X,A,A)') "->  ", "Assign new ID for master points ..."
			pointID = -1                                 ! New Point ID - start from 0
			i = 1
			self % nPointsUnique =0
			DO i=1, self % nPoints
				if (.not. self % horsesPoints (i) % master) cycle
				pointID=pointID+1											! Cannot parallel
				self % nPointsUnique = self % nPointsUnique+1				! Cannot parallel
				ASSOCIATE(p =>self % horsesPoints (i))
				p % ID_NEW = pointID                                           
				END ASSOCIATE
			END DO 
			write(STD_OUT,'(20X,A,A)') "->  ", "Assign new ID for copied points ..."
			
!$omp parallel shared(self)
!$omp do schedule(runtime) 			
			DO i=1, self % nPoints
				if (self % horsesPoints (i) % master) cycle
				ASSOCIATE(p =>self % horsesPoints (i))
				p % ID_NEW = self % horsesPoints (p % copyID) % ID_NEW
				END ASSOCIATE
			END DO 
!$omp end do
!$omp end parallel		
			

!
!           Update Point ID on Face Shared and Boudaries (straight as no double face)
!           -------------------------------------------------------------------------

			write(STD_OUT,'(20X,A,A)') "->  ", "Update PointID on faces ..."


!$omp parallel shared(self)
!$omp do schedule(runtime) private(j)
			DO i=1, self % nFacesShared
				ASSOCIATE(f =>self % faceShared (i))
				DO j=1,4
					f % facePoints(j) = self % horsesPoints (f % facePoints(j)+1) % ID_NEW
				END DO 
				END ASSOCIATE
			END DO 
!$omp end do
!$omp do schedule(runtime) private(j,k)
			DO i=1, size(self % boundaries)
				DO j=1, size(self % boundaries(i) % faceBoundaries)
				ASSOCIATE(f =>self % boundaries(i) % faceBoundaries(j))
					DO k=1,4
						f % facePoints(k) = self % horsesPoints (f % facePoints(k)+1) % ID_NEW
					END DO 
				END ASSOCIATE
				END DO 
			END DO 
!$omp end do
!$omp do schedule(runtime) private(j)
			DO i=1, self % nFacesUnshared
				ASSOCIATE(f =>self % faceUnshared (i))
				DO j=1,4
					f % facePoints(j) = self % horsesPoints (f % facePoints(j)+1) % ID_NEW
				END DO 
				f % master = .true.
				END ASSOCIATE
			END DO 		

!$omp end do
!$omp end parallel			
			
		end subroutine
		
		subroutine detectMultipleFaces(self)
            IMPLICIT NONE
			CLASS(foamMesh_t)                                ,INTENT(INOUT):: self
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,j2,k,l
			LOGICAL 		:: match 

!
!           Detect double faces on Face Unshared, assign neighbour
!           ------------------------------------------------------
			write(STD_OUT,'(20X,A,A)') "->  ", "Detect multiple faces ..."
			self % nMultipleFaces = self % nFacesUnshared/2                    ! Faces should either shared or on boundary
!$omp parallel shared(self)
!$omp do schedule(runtime) private(j,j2,k,l, match)
			DO i=1, self % nFacesUnshared
			
				IF (.not. self % faceUnshared (i) % master) CYCLE
				ASSOCIATE(f =>self % faceUnshared (i))

				DO j=1,	self % nFacesUnshared-1
					
					j2 = i+INT((-1_RP)**(j+1_RP)*CEILING(real(j)/2_RP))
					if (j2.le.0) j2=j2+self % nFacesUnshared
					if (j2.gt.self % nFacesUnshared) j2=j2-self % nFacesUnshared
					
					IF (.not. self % faceUnshared (j2) % master) CYCLE
					DO k=1,3, 2                                                 ! Only need to check 2 points opposite to each other
						match = .false.
						DO l=1,4
							if (f % facePoints(k).eq.self % faceUnshared(j2) % facePoints(l)) then
								match = .true.
								exit
							end if
						END DO
						if (.not.match) EXIT	
					END DO 
					
					IF( match) THEN
						if (j2.gt.i) EXIT 						  ! Master is Face with lower Face ID 
						self % faceUnshared(j2) % faceNeighbour = f % faceOwner
						f % faceNeighbour = self % faceUnshared(j2) % faceOwner
						f % master = .false.
						EXIT
					END IF
				END DO 

				END ASSOCIATE
			END DO
!$omp end do
!$omp end parallel	
			
		end subroutine
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine handle one spectral element´s nodes and construct the foamMesh faces
!
!
        SUBROUTINE finiteElementFacesSet(CellPoints,N,eID,boundaryFace, faceSharedCount, faceUnsharedCount, self)
            IMPLICIT NONE
            INTEGER                                         ,INTENT(IN)   :: N(3)
            INTEGER     ,DIMENSION(0:N(1),0:N(2),0:N(3))    ,INTENT(IN)   :: CellPoints
            INTEGER                                         ,INTENT(IN)   :: eID
            INTEGER                                         ,INTENT(INOUT):: faceUnsharedCount
            INTEGER                                         ,INTENT(INOUT):: faceSharedCount
            INTEGER     ,DIMENSION(1:6*N(1)*N(2)*N(3))      ,INTENT(IN)   :: boundaryFace
            CLASS(foamMesh_t)                                ,INTENT(INOUT):: self
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k
            INTEGER         :: feID,feIDstart,faceSharedStart, faceUnsharedStart, nUnsharedCount
            INTEGER         :: nStart, nEnd
            INTEGER         ,DIMENSION(1:2,1:2,1:2)                 ::  PointSet
            INTEGER         ,DIMENSION(1:7,1:6*N(1)*N(2)*N(3))      ::  FaceSet
            INTEGER         ,DIMENSION(1:7,1:(N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))      &
                                          ::  FaceSetShared
            INTEGER         ,DIMENSION(1:7,1:(2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1))))  &
                                          :: FaceSetUnshared


            feID              =(eID-1)*N(1)*N(2)*N(3)
            feIDstart         = feID
            faceSharedStart   = ((N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))*(eID-1)
            faceUnsharedStart = 2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1))*(eID-1)

            DO k=0,N(3)-1
                DO j=0,N(2)-1
                    DO i=0,N(1)-1
                        feID    =feID+1
                        nStart  =(feID-1)*6+1-feIDstart*6
                        nEnd    =(feID-1)*6+6-feIDstart*6
                        PointSet(:,:,:)=CellPoints(i:i+1,j:j+1,k:k+1)
                        FaceSet(1:4,nStart:nEnd)=facePointsOrder(PointSet)
                        FaceSet(5,nStart:nEnd)=feID-1
                        FaceSet(6,nStart:nEnd)=0
                        FaceSet(7,nStart:nEnd)=0
                    END DO
                END DO
            END DO

!
!        Assign Boundary ID to FE faces
!        --------------------------------------------
            FaceSet(7,:)=boundaryFace(:)
!
!        Sort array from min to max at 1st point
!        --------------------------------------------
            CALL sortItegerArrayMinMax(SHAPE(FaceSet), FaceSet, 1, 0)
!
!        Search for double faces, assign neighbour, and construct the FaceSetShared
!        --------------------------------------------
            CALL assignNeighbour(N,SHAPE(FaceSet) , FaceSet, FaceSetShared, FaceSetUnshared)
!
!        Sort array from min to max at boundary ID
!        --------------------------------------------
            CALL sortItegerArrayMinMax(SHAPE(FaceSetUnshared), FaceSetUnshared, 7, 0)
!
!
!        --------------------------------------------
            nUnsharedCount=0
            DO i=1, 2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1))
                IF (FaceSetUnshared(7,i).EQ.0) THEN
                    nUnsharedCount=nUnsharedCount+1
                ELSE IF (FaceSetUnshared(7,i).GT.0) THEN
                    ASSOCIATE(e => self % boundaries(FaceSetUnshared(7,i)))
                        e % nFace = e % nFace + 1
                        e % faceBoundaries (e % nFace) % facePoints     =  FaceSetUnshared(1:4,i)
                        e % faceBoundaries (e % nFace) % faceOwner      =  FaceSetUnshared(5,i)
                        e % faceBoundaries (e % nFace) % faceNeighbour  =  FaceSetUnshared(6,i)
                        e % faceBoundaries (e % nFace) % faceBoundary   =  FaceSetUnshared(7,i)
						e % faceBoundaries (e % nFace) % master         = .true.
                    END ASSOCIATE
                END IF
            END DO
!
!        Assign Shared Faces
!        --------------------------------------------
            DO i=1,((N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))
                faceSharedCount=faceSharedCount+1
                self % faceShared(faceSharedCount) % facePoints     = FaceSetShared(1:4,i)
                self % faceShared(faceSharedCount) % faceOwner      = FaceSetShared(5,i)
                self % faceShared(faceSharedCount) % faceNeighbour  = FaceSetShared(6,i)
                self % faceShared(faceSharedCount) % faceBoundary   = FaceSetShared(7,i)
				self % faceShared(faceSharedCount) % master         = .true.
            END DO
!
!        Assign Unshared Faces
!        --------------------------------------------
            DO i=1,nUnsharedCount
                faceUnsharedCount=faceUnsharedCount+1
                self % faceUnshared(faceUnsharedCount) % facePoints     = FaceSetUnshared(1:4,i)
                self % faceUnshared(faceUnsharedCount) % faceOwner      = FaceSetUnshared(5,i)
                self % faceUnshared(faceUnsharedCount) % faceNeighbour  = FaceSetUnshared(6,i)
                self % faceUnshared(faceUnsharedCount) % faceBoundary   = FaceSetUnshared(7,i)
				self % faceUnshared(faceUnsharedCount) % master         = .true.
            END DO

        END SUBROUTINE finiteElementFacesSet
!
!////////////////////////////////////////////////////////////////////////
!
! Given a set of points ID of a Hexa element, this function return the 6 set of faces construct from Point ID
! Point order pointing to the element center, modified from original horses2foam
!
        FUNCTION facePointsOrder(CellPoints) RESULT (FaceSet)
            IMPLICIT NONE
            INTEGER   ,DIMENSION(1:2,1:2,1:2)   ,INTENT(IN)     ::  CellPoints
            INTEGER   ,DIMENSION(1:4,1:6)                       ::  FaceSet
!
!        ---------------
!        Local variables
!        ---------------
!
!           Assign Face : Order follows horses .bmesh file

            FaceSet(1,1)=CellPoints(1,1,1)
            FaceSet(2,1)=CellPoints(2,1,1)
            FaceSet(3,1)=CellPoints(2,1,2)
            FaceSet(4,1)=CellPoints(1,1,2)

            FaceSet(1,2)=CellPoints(1,2,1)
            FaceSet(2,2)=CellPoints(1,2,2)
            FaceSet(3,2)=CellPoints(2,2,2)
            FaceSet(4,2)=CellPoints(2,2,1)

            FaceSet(1,3)=CellPoints(1,1,1)
            FaceSet(2,3)=CellPoints(1,2,1)
            FaceSet(3,3)=CellPoints(2,2,1)
            FaceSet(4,3)=CellPoints(2,1,1)

            FaceSet(1,4)=CellPoints(2,1,1)
            FaceSet(2,4)=CellPoints(2,2,1)
            FaceSet(3,4)=CellPoints(2,2,2)
            FaceSet(4,4)=CellPoints(2,1,2)

            FaceSet(1,5)=CellPoints(1,1,2)
            FaceSet(2,5)=CellPoints(2,1,2)
            FaceSet(3,5)=CellPoints(2,2,2)
            FaceSet(4,5)=CellPoints(1,2,2)

            FaceSet(1,6)=CellPoints(1,1,1)
            FaceSet(2,6)=CellPoints(1,1,2)
            FaceSet(3,6)=CellPoints(1,2,2)
            FaceSet(4,6)=CellPoints(1,2,1)

        END FUNCTION facePointsOrder
!
!////////////////////////////////////////////////////////////////////////
!
! This Recursive Subroutine sort Integer ArraySet at nColumn from its Minimum to Maximum Value
! When call for the first time k must be set to 0 -- WARNING DO NOT USE FOR HUGE ARRAY
!
        RECURSIVE SUBROUTINE sortItegerArrayMinMax(sizeArray, ArraySet, nColumn, k)
            IMPLICIT NONE
            INTEGER         ,DIMENSION(2)                               ,INTENT(IN)       :: sizeArray
            INTEGER         ,DIMENSION(1:sizeArray(1),1:sizeArray(2))   ,INTENT(INOUT)    :: ArraySet
            INTEGER                                                     ,INTENT(IN)       :: nColumn
            INTEGER                                                     ,INTENT(IN)       :: k
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         ,DIMENSION(1:sizeArray(1))  :: BufferSet
            INTEGER                                     :: i
!
!        ---------------
!        Perform Sorting Value
!        ---------------
!
            SELECT CASE(k)
            CASE(0)
                DO i=1,sizeArray(2)-1
                    IF (ArraySet(nColumn,i).GT.ArraySet(nColumn,i+1)) THEN
                        BufferSet(:)    = ArraySet(:,i)
                        ArraySet(:,i)   = ArraySet(:,i+1)
                        ArraySet(:,i+1) = BufferSet(:)
                        IF(i.GT.1) THEN
                            CALL sortItegerArrayMinMax(sizeArray, ArraySet, nColumn, i)
                        END IF
                    END IF
                END DO
            CASE DEFAULT ! For Recursive
                IF (ArraySet(nColumn,k-1).GT.ArraySet(nColumn,k)) THEN
                    BufferSet(:)    = ArraySet(:,k-1)
                    ArraySet(:,k-1) = ArraySet(:,k)
                    ArraySet(:,k)   = BufferSet(:)
                    IF(k.GT.2) THEN
                        CALL sortItegerArrayMinMax(sizeArray, ArraySet, nColumn, k-1)
                    END IF
                END IF
            END SELECT

        END SUBROUTINE sortItegerArrayMinMax
!
!////////////////////////////////////////////////////////////////////////
!
! This Subroutine assign Neighbour cell to the 1st face which are used by 2 cell
!
        SUBROUTINE assignNeighbour(N,sizeArray, FaceSet, FaceSetShared, FaceSetUnshared)
            IMPLICIT NONE
            INTEGER                                                     ,INTENT(IN)     :: N(3)
            INTEGER         ,DIMENSION(2)                               ,INTENT(IN)     :: sizeArray
            INTEGER         ,DIMENSION(1:sizeArray(1),1:sizeArray(2))   ,INTENT(INOUT)  :: FaceSet
            INTEGER         ,DIMENSION(1:sizeArray(1),1:(N(1)**2)*(N(1)-1)+(N(2)**2)*(N(2)-1)+(N(3)**2)*(N(3)-1))  &
                    ,INTENT(OUT)    :: FaceSetShared
            INTEGER         ,DIMENSION(1:sizeArray(1),1:2*(N(1)*N(2)+N(2)*N(3)+N(3)*N(1)))  &
                    ,INTENT(OUT)    :: FaceSetUnshared
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k
            k=0
!
!        ---------------
!        Search for double face (1st and 3rd Point Identical) -> Assign neighbour cell to both
!        ---------------
!
            DO i=1,sizeArray(2)-5
                DO j=1,5
                    IF ((FaceSet(1,i).EQ.FaceSet(1,i+j)).AND.(FaceSet(3,i).EQ.FaceSet(3,i+j))) THEN
                        k=k+1
                        FaceSet(6,i)=FaceSet(5,i+j)             ! Assign Neighbour both
                        FaceSet(6,i+j)=1                        ! Assign 1 to prevent error on unshared face filter
                        FaceSetShared(:,k)=FaceSet(:,i)         ! Store the first
                        EXIT
                    END IF
                END DO
            END DO

            k=0
            DO i=1,sizeArray(2)
                IF (FaceSet(6,i).EQ.0) THEN
                    k=k+1
                    FaceSetUnshared(:,k)=FaceSet(:,i)
                END IF
            END DO
        END SUBROUTINE assignNeighbour
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
! This Subroutine rearrange Boundary Set from spectra mesh boundary into array [CellID, elementSide, boundaryID]
! The array is sorted from low to high Cell ID
!
        SUBROUTINE constructSpectralBoundariesArray(mesh, boundarySet)
            IMPLICIT NONE
            CLASS(Mesh_t)                                        ,INTENT(IN)     :: mesh
            INTEGER   , ALLOCATABLE, DIMENSION(:,:)  ,INTENT(OUT)    :: boundarySet
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k=0
            INTEGER         :: nFacesTotal
!
!        ---------------
!        Construct Array for Spectral Boundary and sort the Cell from low to high value
!        ---------------
!
            nFacesTotal=SUM(mesh % boundaries (:) % no_of_faces)
            ALLOCATE(boundarySet(1:3,1:nFacesTotal))

            DO i=1, size(mesh % boundaries)
                DO j=1, mesh % boundaries(i) % no_of_faces
                    k=k+1
                    boundarySet(1,k)= (mesh % boundaries(i) % elements(j)) - 1  ! foam Cell start from 0
                    boundarySet(2,k)= mesh % boundaries(i) % elementSides(j)    ! face Side
                    boundarySet(3,k)= i                                         ! boundaryID
                END DO
            END DO

!        Sorting Cell Value
!        ---------------
            CALL sortItegerArrayMinMax(SHAPE(boundarySet), boundarySet, 1, 0)

        END SUBROUTINE constructSpectralBoundariesArray
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
! This Subroutine rearrange Boundary Set from spectra mesh boundary into array [CellID, elementSide, boundaryID]
! The array is sorted from low to high Cell ID
!
        SUBROUTINE finiteElementBoundariesFace (N, eID, sizeBoundary, boundarySet, boundaryFace)
            IMPLICIT NONE

            INTEGER                                  ,INTENT(IN)    :: N(3)
            INTEGER                                  ,INTENT(IN)    :: eID
            INTEGER                                  ,INTENT(IN)    :: sizeBoundary
            INTEGER   ,DIMENSION(3,sizeBoundary)     ,INTENT(IN)    :: boundarySet
            INTEGER   ,ALLOCATABLE  , DIMENSION(:)   ,INTENT(OUT)   :: boundaryFace
!
!        ---------------
!        Local variables
!        ---------------
!
            INTEGER         :: i,j,k
            INTEGER         :: nStart, nEnd

            k=0
            nStart=0
            nEnd=0
!
!           Allocate boundaryFace for each faces in FE (including double faces)
!           -----------------------------------
            ALLOCATE (boundaryFace(1:6*N(1)*N(2)*N(3)))
!
!        ---------------
!        Count Number of Spectra Faces for element eID
!        ---------------
!
            DO i=1,size(boundarySet,2)
                IF((boundarySet(1,i).EQ.(eID-1)).AND.(k.EQ.0)) THEN
                    k=1
                    nStart=i
                    nEnd  =i
                ELSE IF((boundarySet(1,i).EQ.(eID-1)).AND.(k.NE.0)) THEN
                    nEnd  =i
                ELSE IF (boundarySet(1,i).GT.(eID-1)) THEN
                    EXIT
                END IF
            END DO
            k=0

!
!        ---------------
!        Count Number of Spectra Faces for element eID
!        ---------------
!
            boundaryFace(:)=0               ! Fill 0
            DO i=1, nEnd-nStart+1
                SELECT CASE(boundarySet(2,nStart+i-1))
                    CASE(1)
                        DO j=1, N(3)
                            DO k=1,N(1)
                                boundaryFace((((j-1)*N(1)*N(2)+k)-1)*6+1)=boundarySet(3,nStart+i-1)
                            END DO
                        END DO
                    CASE(2)
                        DO j=1, N(3)
                            DO k=1,N(1)
                                boundaryFace((((j)*N(1)*N(2)-N(1)+k)-1)*6+2)=boundarySet(3,nStart+i-1)
                            END DO
                        END DO
                    CASE(3)
                        DO j=1, N(1)*N(2)
                            boundaryFace((j-1)*6+3)=boundarySet(3,nStart+i-1)
                        END DO
                    CASE(4)
                        DO j=N(1), (N(1)*N(2)*N(3)), N(1)
                            boundaryFace((j-1)*6+4)=boundarySet(3,nStart+i-1)
                        END DO
                    CASE(5)
                        DO j=(N(1)*N(2)*N(3))-(N(1)*N(2))+1, (N(1)*N(2)*N(3))
                            boundaryFace((j-1)*6+5)=boundarySet(3,nStart+i-1)
                        END DO
                    CASE(6)
                        DO j=1, (N(1)*N(2)*N(3)), N(1)
                            boundaryFace((j-1)*6+6)=boundarySet(3,nStart+i-1)
                        END DO
                 END SELECT
            END DO

        END SUBROUTINE finiteElementBoundariesFace
!
END MODULE foamMeshStorageConverter
!
!////////////////////////////////////////////// END OF FILE //////////////////////////////////////////////
!
