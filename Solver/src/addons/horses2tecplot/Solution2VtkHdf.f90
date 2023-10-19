#include "Includes.h"
module Solution2VtkHdfModule

   use SMConstants
   use SolutionFile
   use InterpolationMatrices
   use FileReadingUtilities, only: getFileName
#ifdef HAS_HDF5
   use iso_c_binding, only: c_char, c_loc, c_ptr
   use HDF5
#endif

   implicit none

   private
   public   Solution2VtkHdf
!
!  ========
   contains
!  ========
!
   subroutine Solution2VtkHdf(meshName, solutionName, Nout)
      use Storage
      use NodalStorageClass
      use SharedSpectralBasis
      use OutputVariables
      implicit none
      character(len=*), intent(in) :: meshName
      character(len=*), intent(in) :: solutionName
      integer,          intent(in) :: Nout(3)
#ifdef HAS_HDF5
!
!     ---------------
!     Local variables
!     ---------------
!
      type(Mesh_t)               :: mesh
      character(len=LINE_LENGTH) :: meshPltName
      character(len=LINE_LENGTH) :: solutionFile
      real(kind=RP), allocatable :: xi(:), eta(:), zeta(:)
      integer                    :: Nout__, Nout_(3)
      integer                    :: eID, vID, i, j, k

      character(len=16, kind=c_char), target :: filetype = "UnstructuredGrid"

      integer(HID_T)                                     :: fid, space, attr, dset
      integer(HID_T)                                     :: vtkhdf, fdata, pdata
      integer(HID_T)                                     :: mtype, ftype
      integer                                            :: error
      integer                                            :: nPoints
      integer                                            :: etype
      integer                                            :: ipt
      type(c_ptr),                   target              :: ptr
      integer,                       target              :: vec2(2)
      character(len=:, kind=c_char), target, allocatable :: title
      integer,                               allocatable :: localConnectivities(:)
      integer,                       target, allocatable :: connectivities(:)
      integer,                       target, allocatable :: offsets(:)
      integer,                       target, allocatable :: types(:)
      real(kind=RP),                 target, allocatable :: points(:, :)
      real(kind=RP),                 target, allocatable :: variable(:)
      character(len=STRING_CONSTANT_LENGTH), allocatable :: varNames(:)

!
!     Read the mesh and solution data
!     -------------------------------
      call mesh % ReadMesh(meshName)
      call mesh % ReadSolution(SolutionName)

      ! Use square domains
      if (maxval(Nout) > 0) then
         Nout_ = maxval(Nout)
      else
         Nout__ = 1  ! At least two nodes in each direction
         do eID = 1, mesh % no_of_elements
            Nout__ = max(Nout__, maxval(mesh % elements(eID) % Nsol))
         end do
         Nout_ = Nout__
      end if
      if (mesh % is2D) Nout_(3) = 0
!
!     Set homogeneous nodes
!     ---------------------
      allocate(xi(0:Nout_(1)))
      allocate(eta(0:Nout_(2)))
      allocate(zeta(0:Nout_(3)))

      xi   = [ (-1.0_RP + 2.0_RP * i / Nout_(1), i = 0, Nout_(1)) ]
      eta  = [ (-1.0_RP + 2.0_RP * i / Nout_(2), i = 0, Nout_(2)) ]
      if (mesh % is2D) then
         zeta = [ 0.0_RP ]
      else
         zeta = [ (-1.0_RP + 2.0_RP * i / Nout_(3), i = 0, Nout_(3)) ]
      end if
!
!     Write each element zone
!     -----------------------
      do eID = 1, mesh % no_of_elements
      associate ( e => mesh % elements(eID) )
         e % Nout = Nout_
!
!        Construct spectral basis for both mesh and solution
!        ---------------------------------------------------
         call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
         call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!        Construct interpolation matrices for the mesh
!        ---------------------------------------------
         call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
         call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)
         call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)
!
!        Construct interpolation matrices for the solution
!        -------------------------------------------------
         call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), xi)
         call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), eta)
         call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), zeta)
!
!        Perform interpolation
!        ---------------------
         call ProjectStorageHomogeneousPoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                 Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                 Tset(e % Nout(3), e % Nmesh(3)) % T, &
                                                 Tset(e % Nout(1), e % Nsol(1)) % T,  &
                                                 Tset(e % Nout(2), e % Nsol(2)) % T,  &
                                                 Tset(e % Nout(3), e % Nsol(3)) % T,  &
                                                 mesh % hasGradients, mesh % isStatistics )
      end associate
      end do
!
!     Write the solution file name
!     ----------------------------
      solutionFile = trim(getFileName(solutionName)) // ".hdf"
!
!     Create the file
!     ---------------
      call h5open_f(error)
      call h5fcreate_f(solutionFile, H5F_ACC_TRUNC_F, fid, error)

      call h5gcreate_f(fid, "VTKHDF", vtkhdf, error)

      ! Version
      vec2 = [1, 0]
      call h5screate_simple_f(1, [2_HSIZE_T], space, error)
      call h5acreate_f(vtkhdf, "Version", H5T_STD_I64LE, space, attr, error)

      call h5awrite_f(attr, H5T_NATIVE_INTEGER, c_loc(vec2), error)

      call h5aclose_f(attr, error)
      call h5sclose_f(space, error)

      ! File type
      call h5tcopy_f(H5T_FORTRAN_S1, mtype, error)
      call h5tset_size_f(mtype, len_trim(filetype, kind=8), error)
      call h5tcopy_f(H5T_C_S1, ftype, error)
      call h5tset_size_f(ftype, len_trim(filetype, kind=8), error)
      call h5tset_strpad_f(ftype, H5T_STR_NULLPAD_F, error)

      call h5screate_f(H5S_SCALAR_F, space, error)
      call h5acreate_f(vtkhdf, "Type", ftype, space, attr, error)

      call h5awrite_f(attr, mtype, c_loc(filetype(1:1)), error)

      call h5aclose_f(attr, error)
      call h5sclose_f(space, error)
      call h5tclose_f(mtype, error)
      call h5tclose_f(ftype, error)
!
!     Add the title
!     -------------
      title = "Generated from " // trim(meshName) // " and " // trim(solutionName)

      call h5tcopy_f(H5T_STRING, ftype, error)
      call h5tset_strpad_f(ftype, H5T_STR_NULLPAD_F, error)

      call h5gcreate_f(vtkhdf, "FieldData", fdata, error)
      call h5screate_simple_f(1, [1_HSIZE_T], space, error)
      call h5dcreate_f(fdata, "Title", ftype, space, dset, error)

      ptr = c_loc(title(1:1))
      call h5dwrite_f(dset, ftype, c_loc(ptr), error)

      call h5dclose_f(dset, error)
      call h5sclose_f(space, error)
      call h5gclose_f(fdata, error)
      call h5tclose_f(ftype, error)
!
!     Topology information required by VTKHDF
!     ---------------------------------------
      if (mesh % is2D) then
         etype = 70
      else
         etype = 72
      end if
      nPoints = mesh % no_of_elements * product(Nout_ + 1)

      allocate(points(3, nPoints))
      allocate(connectivities(nPoints))
      allocate(offsets(0:mesh % no_of_elements))
      allocate(types(mesh % no_of_elements))

      ipt = 0
      offsets(0) = 0
      localConnectivities = get_connectivities(etype, Nout_)

      do eID = 1, mesh % no_of_elements
      associate(e => mesh % elements(eID))

          ! Connectivities
          connectivities(ipt + 1:ipt + product(e % Nout + 1)) = ipt + localConnectivities

          ! Points
          do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
             ipt = ipt + 1
             points(:, ipt) = e % xOut(:, i, j, k)
          end do                ; end do                ; end do

          ! Offsets
          offsets(eID) = offsets(eID - 1) + product(e % Nout + 1)

          ! Element types
          types(eID) = etype

      end associate
      end do

      ! First, partition related variables
      call h5screate_simple_f(1, [1_HSIZE_T], space, error)

      vec2(1) = nPoints
      call h5dcreate_f(vtkhdf, "NumberOfPoints", H5T_STD_I64LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, c_loc(vec2), error)
      call h5dclose_f(dset, error)

      call h5dcreate_f(vtkhdf, "NumberOfConnectivityIds", H5T_STD_I64LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, c_loc(vec2), error)
      call h5dclose_f(dset, error)

      vec2(1) = mesh % no_of_elements
      call h5dcreate_f(vtkhdf, "NumberOfCells", H5T_STD_I64LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, c_loc(vec2), error)
      call h5dclose_f(dset, error)

      call h5sclose_f(space, error)

      ! Now, mesh information
      call h5screate_simple_f(2, [3_HSIZE_T, int(nPoints, kind=HSIZE_T)], space, error)
      call h5dcreate_f(vtkhdf, "Points", H5T_IEEE_F64LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, c_loc(points), error)
      call h5dclose_f(dset, error)
      call h5sclose_f(space, error)

      call h5screate_simple_f(1, [int(nPoints, kind=HSIZE_T)], space, error)
      call h5dcreate_f(vtkhdf, "Connectivity", H5T_STD_I64LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, c_loc(connectivities), error)
      call h5dclose_f(dset, error)
      call h5sclose_f(space, error)

      call h5screate_simple_f(1, [int(mesh % no_of_elements + 1, kind=HSIZE_T)], space, error)
      call h5dcreate_f(vtkhdf, "Offsets", H5T_STD_I64LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, c_loc(offsets), error)
      call h5dclose_f(dset, error)
      call h5sclose_f(space, error)

      call h5screate_simple_f(1, [int(mesh % no_of_elements, kind=HSIZE_T)], space, error)
      call h5dcreate_f(vtkhdf, "Types", H5T_STD_U8LE, space, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, c_loc(types), error)
      call h5dclose_f(dset, error)
      call h5sclose_f(space, error)

      deallocate(points)
      deallocate(connectivities)
      deallocate(offsets)
      deallocate(types)
!
!     Add the variables
!     -----------------
      call getOutputVariables()
      call getOutputVariablesList(varNames)

      allocate(variable(nPoints))

      call h5gcreate_f(vtkhdf, "PointData", pdata, error)
      call h5screate_simple_f(1, [int(nPoints, kind=HSIZE_T)], space, error)

      ! Compute output variables
      do eID = 1, mesh % no_of_elements
      associate(e => mesh % elements(eID))
         allocate(e % outputVars(1:no_of_outputVariables, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         call ComputeOutputVariables(no_of_outputVariables, outputVariableNames, e % Nout, e, e % outputVars, &
             mesh % refs, mesh % hasGradients, mesh % isStatistics, mesh % hasSensor)
      end associate
      end do

      ! Serialize and save
      do vID = 1, size(varNames)

         ipt = 0
         do eID = 1, mesh % no_of_elements
         associate(e => mesh % elements(eID))
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               ipt = ipt + 1
               variable(ipt) = e % outputVars(vID, i, j, k)
            end do                ; end do                ; end do
         end associate
         end do

         call h5dcreate_f(pdata, trim(varNames(vID)), H5T_IEEE_F64LE, space, dset, error)
         call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, c_loc(variable), error)
         call h5dclose_f(dset, error)

      end do

      call h5sclose_f(space, error)
      call h5gclose_f(pdata, error)
!
!     Close the file
!     --------------
      call h5gclose_f(vtkhdf, error)
      call h5fclose_f(fid, error)
      call h5close_f(error)

#endif
   end subroutine Solution2VtkHdf

   subroutine ProjectStorageHomogeneousPoints(e, TxMesh, TyMesh, TzMesh, TxSol, TySol, TzSol, hasGradients, hasStats)
      use Storage
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
      logical,             intent(in)  :: hasStats
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: i, j, k, l, m, n
!
!     Project mesh
!     ------------
      allocate( e % xOut(1:3, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
      e % xOut = 0.0_RP

      do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
         do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
            e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
         end do            ; end do            ; end do
      end do            ; end do            ; end do

!
!     Project the solution
!     --------------------
      allocate( e % Qout(1:NVARS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
      e % Qout = 0.0_RP

      do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
         do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
            e % Qout(:,i,j,k) = e % Qout(:,i,j,k) + e % Q(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
         end do            ; end do            ; end do
      end do            ; end do            ; end do

      if ( hasGradients ) then
         allocate( e % U_xout(1:NGRADVARS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)))
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

      if (hasStats) then
         allocate( e % statsout(1:NSTAT, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % statsout = 0.0_RP
         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % statsout(:,i,j,k) = e % statsout(:,i,j,k) + e % stats(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do
      end if

      if (hasUt_NS) then
         allocate( e % ut_NSout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % ut_NSout = 0.0_RP
         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % ut_NSout(:,i,j,k) = e % ut_NSout(:,i,j,k) + e % ut_NS(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do
      end if

      if (hasMu_NS) then
         allocate( e % mu_NSout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % mu_NSout = 0.0_RP
         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % mu_NSout(:,i,j,k) = e % mu_NSout(:,i,j,k) + e % mu_NS(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do
      end if

      if (hasWallY) then
         allocate( e % wallYout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % wallYout = 0.0_RP
         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % wallYout(:,i,j,k) = e % wallYout(:,i,j,k) + e % WallY(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do
      end if

      if (hasMu_sgs) then
         allocate( e % mu_sgsout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % mu_sgsout = 0.0_RP
         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % mu_sgsout(:,i,j,k) = e % mu_sgsout(:,i,j,k) + e % mu_sgs(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do
      end if

   end subroutine ProjectStorageHomogeneousPoints

   function get_connectivities(elemType, Nout) result(conn)
      implicit none
      integer, intent(in) :: elemType
      integer, intent(in) :: Nout(3)
      integer             :: conn(product(Nout + 1))

      integer              :: i
      integer              :: connectivities(0:Nout(1), 0:Nout(2), 0:Nout(3))
      integer, allocatable :: corners(:), edges(:), face(:), faces(:), volume(:)


      ! Array to vector indexing
      connectivities = reshape( [(i-1, i = 1, product(Nout + 1))], shape(connectivities) )

      ! Quadrangle
      if (elemType == 70) then
         corners = [ connectivities(0, 0, 0),             &
                     connectivities(Nout(1), 0, 0),       &
                     connectivities(Nout(1), Nout(2), 0), &
                     connectivities(0, Nout(2), 0) ]

         edges = [ connectivities(1:Nout(1)-1, 0, 0),       &
                   connectivities(Nout(1), 1:Nout(2)-1, 0), &
                   connectivities(1:Nout(1)-1, Nout(2), 0), &
                   connectivities(0, 1:Nout(2)-1, 0) ]

         face = pack(connectivities(1:Nout(1)-1, 1:Nout(2)-1, 0), .true.)

         conn = [ corners, edges, face ]

      ! Hexahedron
      else ! elemType == 72
          corners = [ connectivities(0, 0, 0),                   &
                      connectivities(Nout(1), 0, 0),             &
                      connectivities(Nout(2), Nout(2), 0),       &
                      connectivities(0, Nout(2), 0),             &
                      connectivities(0, 0, Nout(3)),             &
                      connectivities(Nout(1), 0, Nout(3)),       &
                      connectivities(Nout(1), Nout(2), Nout(3)), &
                      connectivities(0, Nout(2), Nout(3)) ]

          edges = [ connectivities(1:Nout(1)-1, 0, 0),             &
                    connectivities(Nout(1), 1:Nout(2)-1, 0),       &
                    connectivities(1:Nout(1)-1, Nout(2), 0),       &
                    connectivities(0, 1:Nout(2)-1, 0),             &
                    connectivities(1:Nout(1)-1, 0, Nout(3)),       &
                    connectivities(Nout(1), 1:Nout(2)-1, Nout(3)), &
                    connectivities(1:Nout(1)-1, Nout(2), Nout(3)), &
                    connectivities(0, 1:Nout(2)-1, Nout(3)),       &
                    connectivities(0, 0, 1:Nout(3)-1),             &
                    connectivities(Nout(1), 0, 1:Nout(3)-1),       &
                    connectivities(Nout(1), Nout(2), 1:Nout(3)-1), &
                    connectivities(0, Nout(2), 1:Nout(3)-1) ]

          faces = [ pack(connectivities(0, 1:Nout(2)-1, 1:Nout(3)-1), .true.),       &
                    pack(connectivities(Nout(1), 1:Nout(2)-1, 1:Nout(3)-1), .true.), &
                    pack(connectivities(1:Nout(1)-1, 0, 1:Nout(3)-1), .true.),       &
                    pack(connectivities(1:Nout(1)-1, Nout(2), 1:Nout(3)-1), .true.), &
                    pack(connectivities(1:Nout(1)-1, 1:Nout(2)-1, 0), .true.),       &
                    pack(connectivities(1:Nout(1)-1, 1:Nout(2)-1, Nout(3)), .true.) ]

          volume = pack(connectivities(1:Nout(1)-1, 1:Nout(2)-1, 1:Nout(3)-1), .true.)

          conn = [ corners, edges, faces, volume ]

      end if

   end function get_connectivities

end module Solution2VtkHdfModule