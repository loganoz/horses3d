!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!   HORSES3D to Foam Result - meshToFoam Module
!
!      This module convert Horses3D mesh storage to foamMesh Storage
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
MODULE meshToFoam
    USE SMConstants
    USE InterpolationMatrices
    USE SharedSpectralBasis
    USE createFoamMeshFile
    USE outputVariablesFoam
    USE storageFoam
    USE NodalStorageClass
    IMPLICIT NONE

!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE generateFoamMesh (meshFile, boundaryFile, resultFile, outputVariables, Nout, generateMesh)
			USE storageFoam
            IMPLICIT NONE
            CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: meshFile, boundaryFile, outputVariables
			CHARACTER(LEN=LINE_LENGTH), INTENT(IN)  :: resultFile
            INTEGER                   , INTENT(IN)     :: Nout(3)
			LOGICAL                   , INTENT(IN)     :: generateMesh
!
!        ---------------
!        Local variables
!        ---------------
!
            type(Mesh_t)                               :: mesh
            integer                                    :: eID, noOfSolutions
            real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
            integer                                    :: i,fid
			integer                       			   :: pos, pos2
			character(len=LINE_LENGTH) 				   :: dir, time
			character(len=LINE_LENGTH), allocatable    :: resultsFile(:)
			logical                                    :: generateMeshStatus
			
			generateMeshStatus = generateMesh
			noOfSolutions = getNoOfCommas(trim(resultFile)) + 1
			allocate( resultsFile(noOfSolutions) )
			if ( noOfSolutions .eq. 1 ) then
			   resultsFile(1) = ADJUSTL(TRIM(resultFile))
            else
               pos=0
               do i = 1, noOfSolutions-1
                  pos2 = index(trim(resultFile(pos+1:)),",") + pos
				  resultsFile(i) = ADJUSTL(TRIM(resultFile(pos+1:pos2-1)))
                  pos = pos2
               end do

               pos = index(trim(resultFile),",",BACK=.true.)
			   resultsFile(i) = ADJUSTL(TRIM(resultFile(pos+1:)))
            end if
!
!        Read the mesh and solution data
!        -------------------------------
            call mesh % ReadMesh(meshFile,boundaryFile)
			
		do i=1, noOfSolutions
			 call mesh % ReadSolution(resultsFile(i))
!
!        Set homogeneous nodes
!        ---------------------
         xi   = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(1),i=0,Nout(1)) /), (/ Nout(1)+1 /) )
         eta  = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(2),i=0,Nout(2)) /), (/ Nout(2)+1 /) )
         zeta = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(3),i=0,Nout(3)) /), (/ Nout(3)+1 /) )

!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )

            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)      ! TODO: check why it was Nmesh(1)
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)     ! TODO: check why it was Nmesh(1)

!
!           Construct interpolation matrices for the solution
!           -------------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), eta)        ! TODO: check why it was Nout(1)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), zeta)       ! TODO: check why it was Nout(1)
!
!           Perform interpolation
!           ---------------------

            call ProjectStorageHomogeneousPoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T, &
                                                     Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                     Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                     Tset(e % Nout(3), e % Nsol(3)) % T, &
                                                                    mesh % hasGradients)
            end associate
         end do
		 
!        Create the result directory
!        ---------------------------
		 CALL getcwd(dir)
		 CALL system('mkdir foamFiles')
		 CALL chdir('foamFiles')
		 
!        Generate Mesh
!        -------------
		 if (generateMeshStatus) then 
		    CALL system('mkdir constant')
			CALL chdir('constant')
			CALL system('mkdir polyMesh')
			CALL chdir('polyMesh')
			CALL createFoamMesh (mesh)
			generateMeshStatus=.FALSE.
		 end if 
		 
		 CALL chdir(trim(dir))
		 CALL chdir('foamFiles')
!
!------------------
!        File: outputVariable paraview fromat
!------------------
!        Create the Output file
!        ----------------------
			 write (time,'(F10.3)') mesh % time
			 CALL system('mkdir -p '// adjustl(trim( time ) ))
			 CALL chdir(adjustl(trim( time ) ))
			 CALL getOutputVariables(outputVariables, mesh)
			 CALL chdir((trim(dir)))
		 end do 
!        Create pointer for paraview
!        ---------------------------
		 CALL chdir('foamFiles')
		 open(fid, file="paraView.foam", status="unknown", action="write")
		 close(fid)

        END SUBROUTINE generateFoamMesh
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine convert horses3d mesh type variables into foamMesh using foamMeshStorage Module
!      Data are then written into files
!
        SUBROUTINE createFoamMesh (mesh)
            USE foamMeshStorage
            IMPLICIT NONE
            TYPE(Mesh_t)                 ,INTENT(INOUT)              :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
            CLASS(foamMesh_t)            ,POINTER      :: foamMesh   => NULL()
            INTEGER                                    :: fid
            INTEGER                                    :: nFace

            write(STD_OUT,'(10X,A,A)') "Creating Output Mesh File:"
            write(STD_OUT,'(10X,A,A)') "-------------------------"
            ALLOCATE (foamMesh)

!        Construct foamMesh Variable
            CALL foamMesh % Construct (mesh)
!
!------------------
!        File: points
!------------------
!           Create the file
            CALL createFilePointsHeader (fid)
!           Write elements
            CALL writeFoamPoints (fid,mesh,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "points"
!
!------------------
!        File: faces
!------------------
!           Create the file
            CALL createFileFacesHeader (fid)
!           Write faces
            CALL writeFoamFaces (fid,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "faces"
!
!------------------
!        File: neighbour
!------------------
!           Create the file
            CALL createFileNeighbourHeader (fid)
!           Write faces´s neighbour cell
            CALL writeFoamNeighbour (fid,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "neighbour"
!
!------------------
!        File: owner
!------------------
!           Create the file
            CALL createFileOwnerHeader (fid)
!           Write faces´s owner cell
            CALL writeFoamOwner (fid,foamMesh)
!           Close the file
            close(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "owner"
!
!------------------
!        File: boundary
!------------------
!           Create the file
            CALL createFileBoundaryHeader (fid)
!           Write boundary
            CALL writeFoamBoundary (fid,foamMesh, mesh)
!           Close the file
            CLOSE(fid)
			write(STD_OUT,'(30X,A,A30,A)') "->", "boundary"
!
!------------------
!        Foam Mesh Statistics
!------------------
            nFace=foamMesh % nFacesShared + foamMesh % nFacesUnshared + foamMesh % nFacesBoundaries
			write(STD_OUT,*)
            write(STD_OUT,'(10X,A,A)') "Output Mesh Statistics:"
            write(STD_OUT,'(10X,A,A)') "----------------------"
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Points:", foamMesh % nPoints
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Faces:", nFace
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Cells:", foamMesh % nCells
			write(STD_OUT,'(30X,A,A30,I10)') "->", "Number of Boundaries:", size(foamMesh % boundaries)
			write(STD_OUT,*)
!
!        Deallocate
            DEALLOCATE (foamMesh)

        END SUBROUTINE createFoamMesh
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine ProjectStorageHomogeneousPoints(e, TxMesh, TyMesh, TzMesh, TxSol, TySol, TzSol, hasGradients)
         use storageFoam
         use NodalStorageClass
         implicit none
         type(Element_t)     :: e
         real(kind=RP),       intent(in)  :: TxMesh(0:e % Nout(1), 0:e % Nmesh(1))
         real(kind=RP),       intent(in)  :: TyMesh(0:e % Nout(2), 0:e % Nmesh(2))
         real(kind=RP),       intent(in)  :: TzMesh(0:e % Nout(3), 0:e % Nmesh(3))
         real(kind=RP),       intent(in)  :: TxSol(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: TySol(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: TzSol(0:e % Nout(3), 0:e % Nsol(3))
         logical,             intent(in)  :: hasGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l, m, n
!
!        Project mesh
!        ------------
         allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % xOut = 0.0_RP

         do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

!
!        Project the solution
!        --------------------
         allocate( e % Qout(1:NVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % Qout = 0.0_RP

         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % Qout(:,i,j,k) = e % Qout(:,i,j,k) + e % Q(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

         if ( hasGradients ) then
            allocate( e % U_xout(1:NGRADVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)))
            e % U_xout = 0.0_RP

            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % U_xout(:,i,j,k) = e % U_xout(:,i,j,k) + e % U_x(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do

          allocate( e % U_yout(1:NGRADVARS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % U_yout = 0.0_RP

            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % U_yout(:,i,j,k) = e % U_yout(:,i,j,k) + e % U_y(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do

          allocate( e % U_zout(1:NGRADVARS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % U_zout = 0.0_RP

            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % U_zout(:,i,j,k) = e % U_zout(:,i,j,k) + e % U_z(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do

         end if

      end subroutine ProjectStorageHomogeneousPoints

END MODULE meshToFoam
