module GeometryClass
   use SMConstants
   use PhysicsStorage
   use NodalStorageClass
   use HexMeshClass
   use FTValueDictionaryClass
   use getInputData_MOD
   use FileReadingUtilities      , only: getFileName, getRealArrayFromString

   private
   public   Geometry_t, Construct

   type Geometry_t
      character(len=LINE_LENGTH) :: solutionName
      contains
         procedure ::  Construct => Geometry_Construct      
         procedure ::  Compute => Geometry_Compute      
         procedure ::  Describe => Geometry_Describe
   end type Geometry_t

   type, extends(Geometry_t)  :: Cube_t
      integer                    :: Npoints
      real(kind=RP)              :: L
      real(kind=RP)              :: xc(NDIM)
      real(kind=RP), allocatable :: x(:,:,:,:)
      real(kind=RP), allocatable :: U(:,:,:,:)
      contains
         procedure ::  Construct => Cube_Construct
         procedure ::  Compute => Cube_Compute
         procedure ::  Describe => Cube_Describe
   end type Cube_t

   type, extends(Geometry_t)  :: Cylinder_t

   end type Cylinder_t

   interface Construct
      module procedure ConstructGeometry
   end interface

   contains
      function ConstructGeometry(mesh, spA, controlVariables)
         implicit none
         class(HexMesh),           intent(in)    :: mesh
         class(NodalStorage_t),      intent(in)    :: spA(0:)
         class(FTValueDictionary), intent(in)    :: controlVariables
         class(Geometry_t),        pointer       :: ConstructGeometry
!
!        Allocate memory for the geometry
!        --------------------------------
         select case ( trim(controlVariables % StringValueForKey(GeometryKey,LINE_LENGTH)) )
         case ("Cube")
            allocate(Cube_t   :: ConstructGeometry)
         end select
!
!        Construct the class components
!        ------------------------------
         call ConstructGeometry % Construct(mesh, spA, controlVariables)
!
!        Add the solution file name
!        --------------------------
         ConstructGeometry % solutionName = trim(getFileName(trim(controlVariables % StringValueForKey(SolutionFileKey,LINE_LENGTH)))) 
!
!        Describe
!        --------
         call ConstructGeometry % Describe
   
      end function ConstructGeometry
!
!//////////////////////////////////////////////////////////
!
!     Geometry procedures
!     -------------------
!
!//////////////////////////////////////////////////////////
!
   subroutine Geometry_Construct(self, mesh, spA, controlVariables)
      implicit none
      class(Geometry_t),        intent(inout) :: self
      class(HexMesh),           intent(in)    :: mesh
      class(NodalStorage_t),      intent(in)    :: spA(0:)
      class(FTValueDictionary), intent(in)    :: controlVariables
!
!     -------------------------------
!     The template class does nothing
!     -------------------------------
!
   end subroutine Geometry_Construct

   subroutine Geometry_Compute(self, mesh, spA)
      implicit none
      class(Geometry_t),        intent(inout) :: self
      class(HexMesh),           intent(in)    :: mesh
      class(NodalStorage_t),      intent(in)    :: spA(0:)
!
!     -------------------------------
!     The template class does nothing
!     -------------------------------
!
   end subroutine Geometry_Compute

   subroutine Geometry_Describe(self)
      implicit none
      class(Geometry_t),   intent(in)  :: self
!
!     -------------------------------
!     The template class does nothing
!     -------------------------------
!
   end subroutine Geometry_Describe
!
!///////////////////////////////////////////////////////////////////////////
!
!     Cube procedures
!     ---------------
!
!///////////////////////////////////////////////////////////////////////////
!
   subroutine Cube_Construct(self, mesh, spA, controlVariables)
      implicit none
      class(Cube_t),            intent(inout) :: self
      class(HexMesh),           intent(in)    :: mesh
      class(NodalStorage_t),      intent(in)    :: spA(0:)
      class(FTValueDictionary), intent(in)    :: controlVariables
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: i, j, k, eID
      real(kind=RP)  :: NDOF
      real(kind=RP)  :: volume, xintegral(3)
      real(kind=RP), allocatable :: xc(:)

!
!     Compute domain volume
!     ---------------------
      volume = 0.0_RP
   
      do eID = 1, mesh % no_of_elements
         associate( e => mesh % elements(eID) ) 
         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )
         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
            volume = volume + spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) * e % geom % jacobian(i,j,k)
         end do                  ; end do                   ; end do
         end associate
         end associate
      end do
!
!     Get the number of points
!     ------------------------
      if ( controlVariables % ContainsKey(NPointsKey) ) then
         self % Npoints = controlVariables % IntegerValueForKey(NPointsKey)
         self % Npoints = self % Npoints + mod(self % Npoints,2)
      else
!
!        Perform an estimation based on the number of DOFs
!        -------------------------------------------------         
         NDOF = 0

         do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) ) 
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1)
            end associate
         end do

         self % Npoints = NDOF ** (1.0_RP / 3.0_RP)
         self % Npoints = self % Npoints + mod(self % Npoints,2)
      end if
!
!     Get the cube length
!     -------------------
      if ( controlVariables % ContainsKey(CubeLengthKey) ) then
         self % L = controlVariables % DoublePrecisionValueForKey(CubeLengthKey)
      else
!
!        Perform an estimation based on the domain volume
!        ------------------------------------------------
         self % L = (volume) ** (1.0_RP / 3.0_RP)

      end if
         
!
!     Get the cube center
!     -------------------
      if ( controlVariables % ContainsKey(CubeCenterKey) ) then
         xc = getRealArrayFromString( controlVariables % StringValueForKey(CubeCenterKey, LINE_LENGTH)) 
         self % xc = xc

      else
!
!        Perform an estimation based on the domain centroid
!        --------------------------------------------------
         xintegral = 0.0_RP
         do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) ) 
            associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                      spAeta  => NodalStorage(e % Nxyz(2)), &
                      spAzeta => NodalStorage(e % Nxyz(3)) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
               xintegral = xintegral + spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) * e % geom % jacobian(i,j,k) &
                                       * e % geom % x(:,i,j,k)
            end do                  ; end do                   ; end do
            end associate
            end associate
         end do
         self % xc = xintegral / volume

      end if

   end subroutine Cube_Construct

   subroutine Cube_Compute(self, mesh, spA)
         implicit none
         class(Cube_t),    intent(inout)  :: self
         class(HexMesh),   intent(in)  :: mesh
         class(NodalStorage_t), intent(in)  :: spA(0:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, allocatable       :: elems(:,:,:)
         real(kind=RP), allocatable :: xi(:,:,:,:)
         integer                    :: fid, i, j, k, n, neighs(7), lastelement
         real(kind=RP), allocatable :: x(:,:,:,:), U(:,:,:,:)
         real(kind=RP)              :: Q(NCONS)
         real(kind=RP)              :: dx
         logical                    :: success
         character(len=LINE_LENGTH) :: solutionFile
!
!        ************************
!        Construct auxiliary mesh
!        ************************
!
         allocate(elems(0:self % Npoints,0:self % Npoints,0:self % Npoints))
         allocate(xi(NDIM,0:self % Npoints,0:self % Npoints,0:self % Npoints))
         allocate(x(NDIM,0:self % Npoints,0:self % Npoints,0:self % Npoints))
         allocate(U(3,0:self % Npoints,0:self % Npoints,0:self % Npoints))

         lastelement = 1

!$omp parallel private(n,neighs,success,Q) firstprivate(lastelement) shared(x,xi,U,elems,mesh,spA)
!$omp do collapse(3) 
         do k = 0, self % Npoints ; do j = 0, self % Npoints  ; do i = 0, self % Npoints
!
!           Get coordinates
!           ---------------
            x(1,i,j,k) = self % xc(1) - 0.5_RP * self % L + self % L * i / self % Npoints
            x(2,i,j,k) = self % xc(2) - 0.5_RP * self % L + self % L * j / self % Npoints
            x(3,i,j,k) = self % xc(3) - 0.5_RP * self % L + self % L * k / self % Npoints
!
!           Get the neighbours. This achieves a (very) large speedup
!           --------------------------------------------------------
            neighs(1) = lastelement
            do n = 1, 6
               if ( mesh % elements(lastelement) % NumberOfConnections(n) .ne. 0 ) then
                  neighs(n+1) = mesh % elements(lastelement) % Connection(n) % globID     ! TODO: Careful if you are using MPI
               else
                  neighs(n+1) = -1
               end if
            end do
!
!           Get the local coordinates and check
!           -----------------------------------
            success = mesh % FindPointWithCoords(x(:,i,j,k),elems(i,j,k),xi(:,i,j,k),neighs)
            if ( .not. success ) then
               print*, "Not success in element ", i,j,k
               error stop
            end if

            lastelement = elems(i,j,k)
!
!           Get the solution in this point
!           ------------------------------
            Q = mesh % elements(elems(i,j,k)) % EvaluateSolutionAtPoint(NCONS,xi(:,i,j,k)) 
            U(1:3,i,j,k) = Q(IRHOU:IRHOW) / Q(IRHO)
         end do            ; end do             ; end do
!$omp end do
!$omp end parallel

!
!        Write a tecplot file
!        --------------------
         write(solutionFile,'(A,A)') trim(self % solutionName), ".cube.tec"
         open(newunit = fid, file = trim(solutionFile), status = "unknown", action = "write") 

         write(fid,'(A)') 'TITLE = "PRUEBA"'
         write(fid,'(A)') 'VARIABLES = "x" "y" "z" "U" "V" "W"'
         write(fid,'(A,I0,A,I0,A,I0,A)') 'ZONE I=',self % Npoints+1,", J=",self % Npoints+1,", K=",self % Npoints+1, ", F=POINT"
         do k = 0, self % Npoints ; do j = 0, self % Npoints  ; do i = 0, self % Npoints
            write(fid,'(6(ES24.16,1X))') x(1:3,i,j,k), U(1:3,i,j,k)
         end do            ; end do             ; end do 
         close(fid)

   end subroutine Cube_Compute

   subroutine Cube_Describe(self)
      use Headers
      implicit none
      class(Cube_t), intent(in)  :: self
      
      write(STD_OUT,'(/)')
      call Section_Header("Interpolation to new geometry")
      write(STD_OUT,'(/)')
      call Subsection_Header("Describing the cube geometry")

      write(STD_OUT,'(30X,A,A30,A,ES10.3,A,ES10.3,A,ES10.3,A)') "->", "Cube centroid: ","[",self % xc(1),", ",self % xc(2),", ",self % xc(3),"]."
      write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Cube side length: ", self % L, "."
      write(STD_OUT,'(30X,A,A30,I0,A)') "->", "Number of points per side: ", self % Npoints, "."
      
   end subroutine Cube_Describe

end module GeometryClass