!
!//////////////////////////////////////////////////////
!
!   @File:    ConstructMeshAndSpectralBasis.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Nov  1 19:56:53 2017
!   @Last revision date: Sun Nov  5 19:15:11 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 431d0b8be7da5b914d9e787fc6ac9e78aceca4ef
!
!//////////////////////////////////////////////////////
!
module ConstructMeshAndSpectralBasis_MOD
   use SolutionFile
   use NodalStorageClass
   use mainKeywordsModule
   use FTValueDictionaryClass
   use HexMeshClass
   use SharedBCModule

   private
   public   ConstructMeshAndSpectralBasis

   contains
      subroutine ConstructMeshAndSpectralBasis(meshFile, solutionFile, mesh, spA)
         use ReadMeshFile
         implicit none
         character(len=*),                intent(in)  :: meshFile
         character(len=*),                intent(in)  :: solutionFile
         type(HexMesh),                   intent(out) :: mesh
         type(NodalStorage), allocatable, intent(out) :: spA(:)
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
            padding = 4*NCONS
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
         allocate( spA(0:max(maxval(Nx), maxval(Ny), maxval(Nz))))
      
         do eID = 1, no_of_elements
            if( .not. spA(Nx(eID)) % Constructed) then   
               call spA(Nx(eID)) % Construct(nodeType, Nx(eID))
            end if

            if( .not. spA(Ny(eID)) % Constructed) then   
               call spA(Ny(eID)) % Construct(nodeType, Ny(eID))
            end if

            if( .not. spA(Nz(eID)) % Constructed) then   
               call spA(Nz(eID)) % Construct(nodeType, Nz(eID))
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
         CALL constructMeshFromFile(mesh, trim(meshfile), nodeType, spA, Nx, Ny, Nz, .true. , success )
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
