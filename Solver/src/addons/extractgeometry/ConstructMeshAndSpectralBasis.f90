!
!//////////////////////////////////////////////////////
!
!   @File:    ConstructMeshAndSpectralBasis.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Nov  1 19:56:53 2017
!   @Last revision date: Tue Feb 13 20:29:37 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 9cdffcbe5af1cc3ea1e17c83c91d73cc17fecde1
!
!//////////////////////////////////////////////////////
!
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

   private
   public   ConstructMeshAndSpectralBasis

   contains
      subroutine ConstructMeshAndSpectralBasis(meshFile, solutionFile, mesh)
         use ReadMeshFile
         implicit none
         character(len=*),                intent(in)  :: meshFile
         character(len=*),                intent(in)  :: solutionFile
         type(HexMesh),                   intent(out) :: mesh
         integer                                      :: no_of_elements
         integer                                      :: nodeType, fileType, fid
         integer                                      :: dims(4), pos
         integer                                      :: eID, iter, padding
         integer, allocatable                         :: Nx(:), Ny(:), Nz(:)
         logical                                      :: success
         real(kind=RP)                                :: time
!
!        Get number of elements, node types, and Euler/NS
!        ------------------------------------------------
         no_of_elements = GetSolutionFileNoOfElements(trim(solutionFile))
         nodeType       = GetSolutionFileNodeType(trim(solutionFile))
         fileType       = GetSolutionFileType(trim(solutionFile))
!
!        Preliminar read to gather polynomial orders from solution file
!        --------------------------------------------------------------          
         allocate(Nx(no_of_elements), Ny(no_of_elements), Nz(no_of_elements))
!
!        Set the file padding
!        --------------------
         if ( fileType .eq. SOLUTION_FILE ) then
            padding = NCONS
         else if ( fileType .eq. SOLUTION_AND_GRADIENTS_FILE ) then
            padding = NCONS + 3 * N_GRAD_EQN
         end if
            
         fid = putSolutionFileInReadDataMode(trim(solutionFile))
!
!        Set the initial position
!        ------------------------
         pos = POS_INIT_DATA

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
         end do
         close(fid)
!      
!        ----------------------
!        Allocate nodal storage      
!        ----------------------
!
         call InitializeNodalStorage( max(maxval(Nx), maxval(Ny), maxval(Nz)) )
         
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
!
!        --------------------------
!        Construct shared BC module
!        --------------------------
!
         call ConstructSharedBCModule
!      
!        --------------
!        Construct mesh
!        --------------
!      
         CALL constructMeshFromFile(mesh, trim(meshfile), nodeType, Nx, Ny, Nz, .true. , 0, success )
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      DO eID = 1, SIZE(mesh % elements) 
         CALL allocateElementStorage( mesh % elements(eID), Nx(eID),Ny(eID),Nz(eID), &
                                      NCONS, N_GRAD_EQN, .true. )
      END DO
!
!     -------------
!     Read solution
!     -------------
!
      call mesh % LoadSolution(trim(solutionFile), iter, time)

      end subroutine ConstructMeshAndSpectralBasis

end module ConstructMeshAndSpectralBasis_MOD
