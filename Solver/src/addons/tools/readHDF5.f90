Module  readHDF5
#ifdef HAS_HDF5
   use HDF5
#endif
   use HexMeshClass
   use ElementClass
   use SMConstants
   use ElementConnectivityDefinitions
   use LocalRefinement
Implicit None

!   
!  ----------------
!  Module variables
!  ----------------
!
#ifdef HAS_HDF5
   integer(HID_T) :: file_id       ! File identifier
#endif
   integer        :: iError        ! Error flag
   integer        :: idx = 0       ! Index of node to add to the list
   real(kind=RP)  :: Lref = 1.0_RP
!   
!  -----------------------------
!  Parameters defined in HOPR io
!  -----------------------------
!   
   INTEGER,PARAMETER              :: ELEM_FirstSideInd=3
   INTEGER,PARAMETER              :: ELEM_LastSideInd=4
   INTEGER,PARAMETER              :: ELEM_FirstNodeInd=5
   INTEGER,PARAMETER              :: ELEM_LastNodeInd=6


contains

! adapted from Horses core Read_HDF5Mesh_HOPR
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------
!  Get the number of elements of mesh
!  ----------------------------------
   function NumOfElems_HDF5( fileName ) result(nelem)
      !----------------------------------
      CHARACTER(LEN=*), intent(in) :: fileName
      integer                      :: nelem
      !----------------------------------
#ifdef HAS_HDF5

!     Initialize FORTRAN predefined datatypes
!     ---------------------------------------
      call h5open_f(iError)
      
!     Open the specified mesh file
!     ----------------------------
      call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, iError) ! instead of H5F_ACC_RDONLY_F one can also use  H5F_ACC_RDWR_F
      
!     Read the number of elements
!     ---------------------------
      CALL GetHDF5Attribute(File_ID,'nElems',1,IntegerScalar=nelem)
#else
      error stop ':: HDF5 is not linked correctly'
#endif
      
   end function NumOfElems_HDF5
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   subroutine ConstructSimpleMesh_FromHDF5File_( self, fileName, locR, AllNx, AllNy, AllNz )
      implicit none
      !-arguments--------------------------------------------------------------
      class(HexMesh)  ,           intent(inout) :: self
      ! class(SimpleHexMesh)  , intent(inout) :: self
      character(LEN=*), optional, intent(in)    :: fileName
      type(LocalRef_t), optional, intent(in)    :: locR
      integer         , optional, intent(in)    :: AllNx(:), AllNy(:), AllNz(:)     !<  Polynomial orders for all the elements
      !-local-variables---------------------------------------------------------
#ifdef HAS_HDF5
      ! Variables as called by Kopriva
      integer  :: numberOfElements  ! ...
      integer                         :: Nx, Ny, Nz     !<  Polynomial orders for each element
      integer  :: bFaceOrder        ! Polynomial order for aproximating curved faces
      INTEGER                          :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
      REAL(KIND=RP)                    :: corners(NDIM,NODES_PER_ELEMENT) ! Corners of element
      
      ! Variables as called in HOPR: For a description, see HOPR documentation
      integer                          :: nUniqueNodes
      integer                          :: nSides, nNodes
      integer         , allocatable    :: GlobalNodeIDs(:)
      double precision, allocatable    :: TempArray(:,:) !(kind=RP)
      real(kind=RP)   , allocatable    :: NodeCoords(:,:)
      integer         , allocatable    :: ElemInfo(:,:)
      integer         , allocatable    :: SideInfo(:,:)
      integer                          :: offset
      integer                          :: first, last
      INTEGER(HSIZE_T),POINTER         :: HSize(:)
      integer                          :: nBCs
      integer                          :: nDims
      CHARACTER(LEN=255), ALLOCATABLE  :: BCNames(:)
      
      ! Auxiliary variables
      integer :: i,j,k,l  ! Counters
      integer                    :: HOPRNodeID           ! Node ID in HOPR
      integer                    :: HCornerMap(8)        ! Map from the corner node index of an element to the local high-order node index used in HOPR
      integer                    :: HSideMap(6)          ! Map from the side index of an element in HORSES3D to the side index used in HOPR
      integer, allocatable       :: HNodeSideMap(:,:,:)  ! Map from the face-node-index of an element to the global node index of HOPR (for surface curvature)
      integer, allocatable       :: HOPRNodeMap(:)       ! Map from the global node index of HORSES3D to the global node index of HOPR
      real(kind=RP), allocatable :: TempNodes(:,:)       ! Nodes read from file to be exported to self % nodes
      logical                    :: CurveCondition
      integer, dimension(8)      :: falseNodeID          ! dummy variable needed to construct element, the actual nodes are not computed
      !---------------------------------------------------------------
      
!
!    ********************************
!    Check if a mesh partition exists
!    ********************************
!
      ! if ( mpi_partition % Constructed ) then
      !    call ConstructMeshPartition_FromHDF5File_( self, fileName, nodes, Nx, Ny, Nz, MeshInnerCurves, dir2D, success ) 
      !    return
      ! end if
      
!     Prepare to read file
!     ------------------------------------
      
      ! Initialize FORTRAN predefined datatypes
      call h5open_f(iError)
      
      ! Open the specified mesh file.
      call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, iError) ! instead of H5F_ACC_RDONLY_F one can also use  H5F_ACC_RDWR_F
        
!
!     Read important mesh attributes
!     ------------------------------
      CALL GetHDF5Attribute(File_ID,'nElems',1,IntegerScalar=numberOfElements)
      CALL GetHDF5Attribute(File_ID,'Ngeo',1,IntegerScalar=bFaceOrder)
      CALL GetHDF5Attribute(File_ID,'nSides',1,IntegerScalar=nSides)
      CALL GetHDF5Attribute(File_ID,'nNodes',1,IntegerScalar=nNodes)
      CALL GetHDF5Attribute(File_ID,'nUniqueNodes',1,IntegerScalar=nUniqueNodes)
      
      allocate(ElemInfo(6,1:numberOfElements))
      call ReadArrayFromHDF5(File_ID,'ElemInfo',2,(/6,numberOfElements/),0,IntegerArray=ElemInfo)
      
      offset=ElemInfo(ELEM_FirstNodeInd,1) ! hdf5 array starts at 0-> -1
      first=offset+1
      last =offset+nNodes
      
      ALLOCATE( GlobalNodeIDs(first:last) )
      CALL ReadArrayFromHDF5(File_ID,'GlobalNodeIDs',1,(/nNodes/),offset,IntegerArray=GlobalNodeIDs)
      
      allocate( NodeCoords(1:3,first:last),TempArray(1:3,first:last) )
      CALL ReadArrayFromHDF5(File_ID,'NodeCoords',2,(/3,nNodes/),offset,RealArray=TempArray)
      NodeCoords = TempArray
      deallocate (TempArray)
      
      offset=ElemInfo(ELEM_FirstSideInd,1) ! hdf5 array starts at 0-> -1  
      first=offset+1
      last =offset+nSides
      ALLOCATE(SideInfo(5,first:last))
      CALL ReadArrayFromHDF5(File_ID,'SideInfo',2,(/5,nSides/),offset,IntegerArray=SideInfo) ! There's a mistake in the documentation of HOPR regarding the SideInfo size!!
      
!!      
!     Some other initializations
!     ---------------------------------------
      self % no_of_elements    = numberOfElements
      self % no_of_allElements = numberOfElements
      ! numberOfBoundaryFaces = 0
      ! corners               = 0.0_RP
      
      ! HSideMap = HOPR2HORSESSideMap()
      HCornerMap = HOPR2HORSESCornerMap(bFaceOrder)
      ! call HOPR2HORSESNodeSideMap(bFaceOrder,HNodeSideMap)
      
      ALLOCATE( self % elements(numberOfelements) )
      
      ! allocate ( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
      ! self % Nx = Nx
      ! self % Ny = Ny
      ! self % Nz = Nz
     
      
!      
!     Now we construct the elements
!     ---------------------------------------

      do l = 1, numberOfElements

         DO k = 1, NODES_PER_ELEMENT
            HOPRNodeID     = ElemInfo(ELEM_FirstNodeInd,l) + HCornerMap(k)
            falseNodeID(k) = HOPRNodeID
            corners(:,k)   = NodeCoords(:,HOPRNodeID) / Lref
         END DO

         if( present(locR) ) then
            ! set dummy values
            falseNodeID = 0
            self % elements(l) % SurfInfo % corners = corners
            call locR% getOrderOfPosition(corners, Nx, Ny, Nz)
            call self % elements(l) % Construct (Nx, Ny, Nz, falseNodeID , l, l)
         else
            call self % elements(l) % Construct (AllNx(l), AllNy(l), AllNz(l), falseNodeID , l, l) 
         end if
         
         
      end do      ! l = 1, numberOfElements
      
!
!     Finish up
!     ---------
!      
      deallocate (ElemInfo)
      deallocate (SideInfo)
      deallocate (NodeCoords)
      deallocate (GlobalNodeIDs)
      
      ! if ( .not. MPI_Process % doMPIAction ) then
      !    deallocate (uNodes, vNodes, values)
      ! end if
!
#else
      error stop ':: HDF5 is not linked correctly'
#endif
   end subroutine ConstructSimpleMesh_FromHDF5File_
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------------------------------
   ! Copied from HOPR
   ! Copyright (C) 2015  Prof. Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
   ! This file is part of HOPR, a software for the generation of high-order meshes.
   !
   ! HOPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
   ! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  -----------------------------------------------------------------------------------------------------------------------
#ifdef HAS_HDF5
   SUBROUTINE GetHDF5Attribute(Loc_ID_in,AttribName,nVal,DatasetName,RealScalar,IntegerScalar,StrScalar,LogicalScalar,&
                                                                  RealArray,IntegerArray)
   !===================================================================================================================================
   ! Subroutine to read attributes from HDF5 file.
   !===================================================================================================================================
   ! MODULES
   ! IMPLICIT VARIABLE HANDLING
      IMPLICIT NONE
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! INPUT VARIABLES
      INTEGER(HID_T), INTENT(IN)           :: Loc_ID_in  ! ?
      INTEGER,INTENT(IN)                              :: nVal  ! ?
      CHARACTER(LEN=*), INTENT(IN)         :: AttribName  ! ?
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: DatasetName  ! ?
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! OUTPUT VARIABLES
      REAL              ,OPTIONAL,INTENT(OUT) :: RealArray(nVal)  ! ?
      INTEGER           ,OPTIONAL,INTENT(OUT) :: IntegerArray(nVal)  ! ?
      REAL              ,OPTIONAL,INTENT(OUT) :: RealScalar  ! ?
      INTEGER           ,OPTIONAL,INTENT(OUT) :: IntegerScalar  ! ?
      LOGICAL           ,OPTIONAL,INTENT(OUT) :: LogicalScalar  ! ?
      CHARACTER(LEN=255),OPTIONAL,INTENT(OUT) :: StrScalar  ! ?
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! LOCAL VARIABLES
      INTEGER(HID_T)                 :: Attr_ID, Type_ID,Loc_ID  ! ?
      INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf  ! ?
      INTEGER                        :: inttolog  ! ?
      !===================================================================================================================================

      Dimsf(1)=nVal
      IF(PRESENT(DatasetName))THEN
       ! Open dataset
        CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
      ELSE
        Loc_ID=Loc_ID_in
      END IF
      ! Create scalar data space for the attribute.
      ! Create the attribute for group Loc_ID.
      CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)
      ! Write the attribute data.
      IF(PRESENT(RealArray))THEN
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
      END IF
      IF(PRESENT(RealScalar))THEN
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
      END IF
      IF(PRESENT(IntegerArray))THEN
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER , IntegerArray, Dimsf, iError)
      END IF
      IF(PRESENT(IntegerScalar))THEN
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER , IntegerScalar, Dimsf, iError)
      END IF
      IF(PRESENT(LogicalScalar))THEN
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER , inttolog, Dimsf, iError)
        LogicalScalar=(inttolog.EQ.1)
      END IF
      IF(PRESENT(StrScalar))THEN
        CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)  ! Get HDF5 data type for character string
        CALL H5AREAD_F(Attr_ID, Type_ID, StrScalar, Dimsf, iError)
        CALL H5TCLOSE_F(Type_ID, iError)
        
      END IF
      ! Close the attribute.
      CALL H5ACLOSE_F(Attr_ID, iError)
      IF(PRESENT(DataSetName))THEN
        ! Close the dataset and property list.
        CALL H5DCLOSE_F(Loc_ID, iError)
      END IF

   END SUBROUTINE GetHDF5Attribute
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------------------------------
   ! Copied from HOPR.
   !  -> Corrected the fact that RealArray had to be defined as double precision!!
   ! Copyright (C) 2015  Prof. Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
   ! This file is part of HOPR, a software for the generation of high-order meshes.
   !
   ! HOPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
   ! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  -----------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ReadArrayFromHDF5(Loc_ID,ArrayName,Rank,nVal,Offset_in,RealArray,IntegerArray,StrArray)
   !===================================================================================================================================
   ! Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
   !===================================================================================================================================
   ! MODULES
   ! IMPLICIT VARIABLE HANDLING
      IMPLICIT NONE
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! INPUT VARIABLES
      INTEGER, INTENT(IN)                :: Rank ! ?
      INTEGER, INTENT(IN)                :: Offset_in  ! ?
      INTEGER, INTENT(IN)                            :: nVal(Rank)  ! ?
      INTEGER(HID_T), INTENT(IN)         :: Loc_ID  ! ?
      CHARACTER(LEN=*),INTENT(IN)        :: ArrayName  ! ?
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! OUTPUT VARIABLES
      double precision              ,DIMENSION(Rank),OPTIONAL,INTENT(OUT) :: RealArray  ! ?
      INTEGER           ,DIMENSION(Rank),OPTIONAL,INTENT(OUT) :: IntegerArray  ! ?
      CHARACTER(LEN=255),DIMENSION(Rank),OPTIONAL,INTENT(OUT) :: StrArray  ! ?
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! LOCAL VARIABLES
      INTEGER(HID_T)                 :: DSet_ID, Type_ID, MemSpace, FileSpace, PList_ID  ! ?
      INTEGER(HSIZE_T)               :: Offset(Rank),Dimsf(Rank)  ! ?
      !===================================================================================================================================
      ! Read array -----------------------------------------------------------------------------------------------------------------------
      Dimsf=nVal
      CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
      CALL H5DOPEN_F(Loc_ID, TRIM(ArrayName) , DSet_ID, iError)
      ! Define and select the hyperslab to use for reading.
      CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
      Offset(:)=0
      Offset(1)=Offset_in
      CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
      ! Create property list
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
      ! Read the data
      IF(PRESENT(RealArray))THEN
        CALL H5DREAD_F(DSet_ID,H5T_NATIVE_DOUBLE,&
                        RealArray    ,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
      END IF
      IF(PRESENT(IntegerArray))THEN
        CALL H5DREAD_F(DSet_ID,H5T_NATIVE_INTEGER, &
                        IntegerArray ,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
      END IF
      IF(PRESENT(StrArray))THEN
        ! Get datatype for the character string array
        CALL H5DGET_TYPE_F(DSet_ID, Type_ID, iError)
        CALL H5DREAD_F(DSet_ID,Type_ID,&
                        StrArray     ,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
        CALL H5TCLOSE_F(Type_ID, iError)
      END IF

      ! Close the property list
      CALL H5PCLOSE_F(PList_ID,iError)
      ! Close the file dataspace
      CALL H5SCLOSE_F(FileSpace,iError)
      ! Close the dataset
      CALL H5DCLOSE_F(DSet_ID, iError)
      ! Close the memory dataspace
      CALL H5SCLOSE_F(MemSpace,iError)
   
   END SUBROUTINE ReadArrayFromHDF5
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------------------------------
   ! Copied from HOPR.
   ! Copyright (C) 2015  Prof. Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
   ! This file is part of HOPR, a software for the generation of high-order meshes.
   !
   ! HOPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
   ! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  -----------------------------------------------------------------------------------------------------------------------
   SUBROUTINE GetHDF5DataSize(Loc_ID,DSetName,nDims,Size)
   !===================================================================================================================================
   ! Subroutine to get the size of an array stored in the hdf5 file
   !===================================================================================================================================
   ! MODULES
   ! IMPLICIT VARIABLE HANDLING
      IMPLICIT NONE
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! INPUT VARIABLES
      CHARACTER(LEN=*),INTENT(IN)               :: DSetName  ! ?
      INTEGER(HID_T),INTENT(IN)      :: Loc_ID  ! ?
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! OUTPUT VARIABLES
      INTEGER,INTENT(OUT)            :: nDims  ! ?
      INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)  ! ?
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! LOCAL VARIABLES
      INTEGER(HID_T)                 :: DSet_ID,FileSpace  ! ?
      INTEGER(HSIZE_T), POINTER      :: SizeMax(:)  ! ?
      !===================================================================================================================================
      !WRITE(UNIT_stdOut,'(A,A,A)')'GET SIZE OF "',TRIM(DSetName),'" IN HDF5 FILE... '
      ! Initialize FORTRAN predefined datatypes

      ! Get size of array ----------------------------------------------------------------------------------------------------------------
      ! Open the dataset with default properties.
      CALL H5DOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
      ! Get the data space of the dataset.
      CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
      ! Get number of dimensions of data space
      CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
      ! Get size and max size of data space
      ALLOCATE(Size(nDims),SizeMax(nDims))
      CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
      CALL H5SCLOSE_F(FileSpace, iError)
      CALL H5DCLOSE_F(DSet_ID, iError)

      !WRITE(UNIT_stdOut,*)'...DONE!'
      !WRITE(UNIT_stdOut,'(132("-"))')

   END SUBROUTINE GetHDF5DataSize

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------------------------------
!  Mapping between corner nodes needed in HORSES3D and the high-order nodes used by HOPR 
!  The mapping is the same as the one for the CNGS standard.
!  -----------------------------------------------------------------------------------------------------------------------
   pure function HOPR2HORSESCornerMap(N) result(CGNSCornerMap)
      implicit none
      !-----------------------------
      integer, intent(in)   :: N !<  Order of boundaries
      integer, dimension(8) :: CGNSCornerMap
      !-----------------------------
      CGNSCornerMap(1) =  1
      CGNSCornerMap(2) = (N+1)
      CGNSCornerMap(3) = (N+1)**2
      CGNSCornerMap(4) =  N*(N+1)+1
      CGNSCornerMap(5) =  N*(N+1)**2+1
      CGNSCornerMap(6) =  N*(N+1)**2+(N+1)
      CGNSCornerMap(7) = (N+1)**3
      CGNSCornerMap(8) =  N*(N+1)*(N+2)+1
      
   end function HOPR2HORSESCornerMap
!
#endif
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
End Module  readHDF5 