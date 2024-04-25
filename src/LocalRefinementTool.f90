!
!//////////////////////////////////////////////////////
!
! This module contains the pre-define local p-refinement of a mesh
!
!//////////////////////////////////////////////////////
!
Module LocalRefinementTool  !

    use SMConstants
    use FTValueDictionaryClass

    Implicit None

    private
    public LocalRef

!////////////////////////////////////////////////////////////////////////
    contains
!////////////////////////////////////////////////////////////////////////
!
    Subroutine LocalRef(controlVariables)

        use readHDF5
        use readSpecM
        use LocalRefinement
        use HexMeshClass
        use SMConstants
        use FTValueDictionaryClass
        use Headers
        use MPI_Process_Info

        Implicit None

        TYPE( FTValueDictionary)                :: controlVariables
!
!     ---------------
!     Local variables
!     ---------------
!
        CHARACTER(LEN=LINE_LENGTH)              :: fileName, meshFileName
        type(HexMesh)                           :: mesh
        logical                                 :: success
        type(LocalRef_t)                        :: locR
        character(len=LINE_LENGTH)              :: msg

        call CheckInputIntegrity(controlVariables, success)
        IF(.NOT. success)   error stop "Control file reading error"

!
!       ---------------------------
!       Set up the local refinement
!       ---------------------------
!
        call locR % Construct(controlVariables)
        meshFileName = controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH) 
        fileName = controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH) 

        print *, "start mesh reading"
        call ConstructSimpleMesh(mesh, meshFileName, locR)
!
!       -----------------
!       Describe the mesh
!       -----------------
!
        write(STD_OUT,'(/)')
        call Section_Header("Job description")

        write(msg,'(A,A,A)') 'Mesh file "',trim(meshFileName),'":'
        write(STD_OUT,'(/)')
        call SubSection_Header(trim(msg))
        write(STD_OUT,'(30X,A,A30,I0)') "->", "Number of elements: ", mesh % no_of_elements
        write(STD_OUT,'(/)')
!
!       ---------------------
!       Create the final file
!       ---------------------
!
        call mesh % ExportOrders(fileName)
!
!       --------------------------
!       Show simulation statistics
!       --------------------------
!
        call DisplaySimulationStatistics(mesh)

    End Subroutine LocalRef
!
!////////////////////////////////////////////////////////////////////////
!   Reading and displaying info routines
!////////////////////////////////////////////////////////////////////////
!
    Subroutine CheckInputIntegrity( controlVariables, success )  
       use SMConstants
       use Utilities, only: toLower
       USE FTValueDictionaryClass
       use FTValueClass
       IMPLICIT NONE
!
       TYPE(FTValueDictionary) :: controlVariables
       LOGICAL                 :: success
!
!      ---------------
!      Local variables
!      ---------------
!
       CLASS(FTObject), POINTER :: obj

       success = .TRUE.
!
!      Control variables with default value
!      ------------------------------------
       obj => controlVariables % objectForKey("z regions orders")
       if ( .not. associated(obj) ) then
          call controlVariables % addValueForKey("0","z regions orders")
       end if
       obj => controlVariables % objectForKey("z regions limits")
       if ( .not. associated(obj) ) then
          call controlVariables % addValueForKey("0","z regions limits")
       end if
!
!
!        Check the controlVariables created
!        ----------------------------------        
   ! DO i = 1, SIZE(mainKeywords)
   !    obj => controlVariables % objectForKey(mainKeywords(i))
   !    IF ( .NOT. ASSOCIATED(obj) )     THEN
   !       PRINT *, "Input file is missing entry for keyword: ",mainKeywords(i)
   !       success = .FALSE. 
   !    END IF  
   ! END DO  
   
    End Subroutine checkInputIntegrity
!
   Subroutine DisplaySimulationStatistics(mesh)
       use SMConstants
       use Headers
       use HexMeshClass

       implicit none
       type(HexMesh),   intent(in) :: mesh
!
!      ---------------
!      Local variables
!      ---------------
!
       integer                    :: eID
       integer                    :: NDOF, localNDOF, ierr
       real(kind=RP)              :: Naverage, localNaverage
     
       call Section_Header("Simulation statistics")
        write(STD_OUT,'(/)')
!
!      Get mesh-related quantities
!      ---------------------------
       NDOF = 0
       Naverage = 0
     
       do eID = 1, mesh % no_of_elements
          associate ( e => mesh % elements(eID) )
          NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1)      
          Naverage = Naverage + e % Nxyz(1) + e % Nxyz(2) + e % Nxyz(3)
          end associate
       end do

       Naverage = Naverage / (3.0_RP * mesh % no_of_elements)
!
!          Show preprocessing time
!          -----------------------
       call Subsection_Header("omesh file")

       write(STD_OUT,'(30X,A,A30,F5.2)')      "->   ", "Average polynomial order: ", Naverage
       write(STD_OUT,'(30X,A,A30,I0)')      "->   ", "Degrees of Freedom (NDOF): ", NDOF

!
    end subroutine DisplaySimulationStatistics
!
    Subroutine ConstructSimpleMesh(mesh, meshFileName, locR)

       use SMConstants
       use Headers
       use HexMeshClass
       use LocalRefinement
       use readHDF5
       use readSpecM
       use readGMSH
       use FileReadingUtilities, only: getFileExtension
       Implicit None

!          ---------------
!          Input variables
!          ---------------
!
           type(HexMesh)                    :: mesh
           CHARACTER(LEN=*)                 :: meshFileName
           type(LocalRef_t), intent(in)       :: locR
!
!          ---------------
!          Local variables
!          ---------------
!
           character(len=LINE_LENGTH) :: ext
           integer                    :: gmsh_version

           ext = getFileExtension(trim(meshFileName))
           if (trim(ext)=='h5') then
               call ConstructSimpleMesh_FromHDF5File_(mesh, meshFileName, locR=locR)
           elseif (trim(ext)=='mesh') then
               call ConstructSimpleMesh_FromSpecFile_(mesh, meshFileName, locR=locR)
           elseif (trim(ext)=='msh') then
               call CheckGMSHversion (meshFileName, gmsh_version)
               select case (gmsh_version)
                  case (4)
                     call ConstructSimpleMesh_FromGMSHFile_v4_( mesh, meshFileName, locR=locR )
                  case (2)
                     call ConstructSimpleMesh_FromGMSHFile_v2_( mesh, meshFileName, locR=locR )
                  case default
                     error stop "ReadMeshFile :: Unrecognized GMSH version."
               end select
           else
               error stop 'Mesh file extension not recognized.'
           end if

    End Subroutine ConstructSimpleMesh
!
End Module LocalRefinementTool
