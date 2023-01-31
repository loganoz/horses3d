module ConstructMeshAndSpectralBasis_MOD
   use SMConstants
   use SolutionFile
   use PhysicsStorage
   use NodalStorageClass
   use ElementClass
   use mainKeywordsModule
   use FTValueDictionaryClass
   use HexMeshClass
   use SharedBCModule
   use FTValueDictionaryClass
   use FileReaders               , only: AssignControlFileName
   use InterpolationMatrices     , only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices

   private
   public   ConstructMeshAndSpectralBasis

   contains
      subroutine ConstructMeshAndSpectralBasis(meshFile, solutionFile, ControlFile, mesh, controlVariables)
         use ReadMeshFile
         use PartitionedMeshClass
         implicit none
         character(len=*),                intent(in)  :: meshFile
         character(len=*),                intent(in)  :: solutionFile
         character(len=*),                intent(in)  :: ControlFile
         type(HexMesh),                   intent(out) :: mesh
         type(FTValueDictionary)                      :: controlVariables
         integer                                      :: no_of_elements
         integer                                      :: nodeType, fileType, fid
         integer                                      :: dims(4), pos
         integer                                      :: eID, iter, padding
         integer, allocatable                         :: Nx(:), Ny(:), Nz(:)
         logical                                      :: success
         real(kind=RP)                                :: time
         integer                                      :: NDOF
         logical                                      :: useRelaxPeriodic
         
         call AssignControlFileName(controlFile)
         
!
!        Get number of elements, node types, and Euler/NS
!        ------------------------------------------------
         no_of_elements = GetSolutionFileNoOfElements(trim(solutionFile))
         nodeType       = GetSolutionFileNodeType(trim(solutionFile))
         fileType       = GetSolutionFileType(trim(solutionFile))
!
!        Preliminary read to gather polynomial orders from solution file
!        ---------------------------------------------------------------
         allocate(Nx(no_of_elements), Ny(no_of_elements), Nz(no_of_elements))
!
!        Set the file padding
!        --------------------
         if ( fileType .eq. SOLUTION_FILE ) then
            padding = NCONS
         else if ( fileType .eq. SOLUTION_AND_GRADIENTS_FILE ) then
            padding = NCONS + 3 * NGRAD
         end if
            
         fid = putSolutionFileInReadDataMode(trim(solutionFile))
!
!        Set the initial position
!        ------------------------
         pos = POS_INIT_DATA

         NDOF = 0
         do eID = 1, no_of_elements
            call getSolutionFileArrayDimensions(fid,dims,pos)         
            Nx(eID) = dims(2) - 1
            Ny(eID) = dims(3) - 1
            Nz(eID) = dims(4) - 1
!
!           Skip to the next element data
!           -----------------------------
            pos = pos + 5*SIZEOF_INT + &
                        padding*(Nx(eID)+1)*(Ny(eID)+1)*(Nz(eID)+1)*SIZEOF_RP 

            NDOF = NDOF + (Nx(eID)+1)*(Ny(eID)+1)*(Nz(eID)+1)
         end do
         close(fid)
!      
!        ----------------------
!        Allocate nodal storage and interpolation matrices
!        ----------------------
!
         call InitializeNodalStorage( nodeType, max(maxval(Nx), maxval(Ny), maxval(Nz)) )
      call NodalStorage(0) % Construct(CurrentNodes, 0)   ! Always construct orders 0
      call NodalStorage(1) % Construct(CurrentNodes, 1)   ! and 1
         
         do eID = 1, no_of_elements
            if( .not. NodalStorage(Nx(eID)) % Constructed) then   
               call NodalStorage(Nx(eID)) % Construct(nodeType, Nx(eID))
            end if

            if( .not. NodalStorage(Ny(eID)) % Constructed) then   
               call NodalStorage(Ny(eID)) % Construct(nodeType, Ny(eID))
            end if

            if( .not. NodalStorage(Nz(eID)) % Constructed) then   
               call NodalStorage(Nz(eID)) % Construct(nodeType, Nz(eID))
            end if
         end do
         
         call Initialize_InterpolationMatrices( max(maxval(Nx), maxval(Ny), maxval(Nz)))
!
!        --------------------------
!        Construct shared BC module
!        --------------------------
!
         call ConstructSharedBCModule
         CALL ConstructPhysicsStorage( controlVariables, success )
!      
!        --------------
!        Construct mesh
!        --------------
!      
         ! call Initialize_MPI_Partitions ( trim(controlVariables % stringValueForKey('partitioning', requestedLength = LINE_LENGTH)) )
         useRelaxPeriodic = controlVariables % logicalValueForKey("periodic relative tolerance")
         CALL constructMeshFromFile(mesh, trim(meshfile), nodeType, Nx, Ny, Nz, .true. , 0, useRelaxPeriodic, success )
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      call mesh % AllocateStorage(NDOF, controlVariables, .true.)
!
!     -------------
!     Read solution
!     -------------
!
      call mesh % LoadSolution(trim(solutionFile), iter, time)
      
      call Finalize_InterpolationMatrices

      end subroutine ConstructMeshAndSpectralBasis

end module ConstructMeshAndSpectralBasis_MOD