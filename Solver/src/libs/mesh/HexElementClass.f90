!
!////////////////////////////////////////////////////////////////////////
!
!      ElementClass.f95
!      Created: 2008-06-04 15:34:44 -0400 
!      By: David Kopriva
!
!      Implements Algorithms:
!         Algorithm: 124: ElementClass (QuadElementClass)
!
!       The Quad Element class, Alg. 124. See Sec. 8.2.1.2. self has
!       been modified to add the association of a boundary name to an element
!       edge so that different boundary conditions can be applied to different
!       elements. The names of the boundaries (not necessarily the names of the
!       *boundary conditions* to be applied) are of length BC_STRING_LENGTH.
!       One will associate boundary conditions to boundaries in the routine
!       "ExternalState".
!
!       Modified 2D Code to move solution into element class. 5/14/15, 5:36 PM
!
!////////////////////////////////////////////////////////////////////////
!
      Module ElementClass
      USE SMConstants
      USE PolynomialInterpAndDerivsModule
      USE GaussQuadrature
      USE TransfiniteMapClass
      USE MappedGeometryClass
      USE MeshTypes
      USE ElementConnectivityDefinitions
      USE ConnectivityClass
      use StorageClass
      USE NodalStorageClass
      use PhysicsStorage
      IMPLICIT NONE
      
      
      TYPE Element
          integer                                        :: eID               ! ID of this element
          INTEGER                                        :: nodeIDs(8)
          INTEGER, DIMENSION(3)                          :: Nxyz              ! Polynomial orders in every direction (Nx,Ny,Nz)
          TYPE(MappedGeometry)                           :: geom
          CHARACTER(LEN=BC_STRING_LENGTH)                :: boundaryName(6)
          CHARACTER(LEN=BC_STRING_LENGTH)                :: boundaryType(6)
          INTEGER                                        :: NumberOfConnections(6)
          TYPE(Connectivity)                             :: Connection(6)
          type(Storage_t)                                :: storage
          type(NodalStorage), pointer                    :: spAxi
          type(NodalStorage), pointer                    :: spAeta
          type(NodalStorage), pointer                    :: spAzeta
          type(TransfiniteHexMap)                        :: hexMap            ! High-order mapper
          contains
            procedure   :: FindPointWithCoords => HexElement_FindPointWithCoords
            procedure   :: EvaluateSolutionAtPoint => HexElement_EvaluateSolutionAtPoint
      END TYPE Element 
      
!
!     -------------------------------------------------------------------------
!!    axisMap gives the element local coordinate number for the two directions
!!    on each face. The coordinate numbers are given by (xi,eta,zeta) = (1,2,3).
!!    For instance, the two coordinate directions on Face 1 are (xi,zeta).
!     -------------------------------------------------------------------------
!
      INTEGER, DIMENSION(2,6) :: axisMap =                        &
                                 RESHAPE( (/1, 3,                 & ! Face 1 (x,z)
                                            1, 3,                 & ! Face 2 (x,z)
                                            1, 2,                 & ! Face 3 (x,y)
                                            2, 3,                 & ! Face 4 (y,z)
                                            1, 2,                 & ! Face 5 (x,y)
                                            2, 3/)                & ! Face 6 (y,z)
                                 ,(/2,6/))
            
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructElementGeometry( self, spAxi, spAeta, spAzeta, nodeIDs, hexMap , eID)
         IMPLICIT NONE
         
         TYPE(Element)           :: self
         TYPE(NodalStorage), target :: spAxi
         TYPE(NodalStorage), target :: spAeta
         TYPE(NodalStorage), target :: spAzeta
         INTEGER                 :: nodeIDs(8)
         TYPE(TransfiniteHexMap) :: hexMap
         integer                 :: eID
         
         self % eID                   = eID
         self % nodeIDs               = nodeIDs
         self % Nxyz(1)               = spAxi   % N
         self % Nxyz(2)               = spAeta  % N
         self % Nxyz(3)               = spAzeta % N
         self % boundaryName          = emptyBCName
         self % boundaryType          = emptyBCName
         self % spAxi   => spAxi
         self % spAeta  => spAeta
         self % spAzeta => spAzeta
         self % hexMap = hexMap
!
!        --------
!        Geometry
!        --------
!
         CALL ConstructMappedGeometry( self % geom, spAxi, spAeta, spAzeta, hexMap )
!
!        ----------------------------------------
!        Solution Storage is allocated separately
!        ----------------------------------------
!
      END SUBROUTINE ConstructElementGeometry
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE allocateElementStorage(self, Nx, Ny, Nz, nEqn, nGradEqn, flowIsNavierStokes)  
         IMPLICIT NONE
         TYPE(Element)        :: self
         INTEGER, intent(in)  :: Nx, Ny, Nz, nEqn, nGradEqn
         LOGICAL, intent(in)  :: flowIsNavierStokes

         call self % Storage % Construct(Nx, Ny, Nz, nEqn, nGradEqn, flowIsNavierStokes)

      END SUBROUTINE allocateElementStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetElementBoundaryNames( self, names ) 
         IMPLICIT NONE
         TYPE(Element)                   :: self
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(6)
         INTEGER                         :: j
         
         DO j = 1, 6
            CALL toLower(names(j)) 
            self % boundaryName(j) = names(j)
         END DO  
      END SUBROUTINE SetElementBoundaryNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructElement( self )
         IMPLICIT NONE
         TYPE(Element) :: self
         
         CALL DestructMappedGeometry( self % geom )
         call self % Storage % Destruct   
         
         nullify( self % spAxi   )
         nullify( self % spAeta  )
         nullify( self % spAzeta )     

      END SUBROUTINE DestructElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintElement( self, id )
         IMPLICIT NONE 
         TYPE(Element) :: self
         INTEGER      :: id
         PRINT *, id, self % nodeIDs
         PRINT *, "   ",self % boundaryName
      END SUBROUTINE PrintElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionStorageToUnit( self, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: self
         INTEGER       :: fUnit
         
         WRITE(funit) self % storage % Q
      
      END SUBROUTINE SaveSolutionStorageToUnit
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionFromUnit( self, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: self
         INTEGER       :: fUnit
         
         READ(funit) self % storage % Q
      
      END SUBROUTINE LoadSolutionFromUnit
!
!////////////////////////////////////////////////////////////////////////
!
      logical function HexElement_FindPointWithCoords(self, spA, x, xi)
!
!        *************************+********************************
!          
!           This function finds whether a point is inside or not 
!           of the element. This is done solving
!           the mapping non-linear system
!
!        *************************+********************************
!          
!
         implicit none
         class(Element),      intent(in)  :: self
         class(NodalStorage), intent(in)  :: spA
         real(kind=RP),       intent(in)  :: x(NDIM)
         real(kind=RP),       intent(out) :: xi(NDIM)
!
!        ----------------------------------
!        Newton iterative solver parameters
!        ----------------------------------
!
         integer,       parameter   :: N_MAX_ITER = 50
         real(kind=RP), parameter   :: TOL = 1.0e-12_RP
         integer,       parameter   :: STEP = 1.0_RP
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j, k, iter
         real(kind=RP), parameter      :: INSIDE_TOL = 1.0e-08_RP
         real(kind=RP)                 :: lxi   (0:self % Nxyz(1)) 
         real(kind=RP)                 :: leta  (0:self % Nxyz(2)) 
         real(kind=RP)                 :: lzeta (0:self % Nxyz(3)) 
         real(kind=RP)                 :: dlxi   (0:self % Nxyz(1)) 
         real(kind=RP)                 :: dleta  (0:self % Nxyz(2)) 
         real(kind=RP)                 :: dlzeta (0:self % Nxyz(3)) 
         real(kind=RP)                 :: F(NDIM)
         real(kind=RP)                 :: Jac(NDIM,NDIM)
         real(kind=RP)                 :: dx(NDIM)
         interface
            function SolveThreeEquationLinearSystem(A,b)
               use SMConstants
               implicit none
               real(kind=RP), intent(in)  :: A(3,3)
               real(kind=RP), intent(in)  :: b(3)
               real(kind=RP)     :: SolveThreeEquationLinearSystem(3)
            end function SolveThreeEquationLinearSystem
         end interface
!
!        Initial seed
!        ------------      
         xi = 0.0_RP    

         do iter = 1 , N_MAX_ITER
!
!           Get Lagrange polynomials and derivatives
!           ----------------------------------------
            lxi     = spA % lxi   (xi(1))
            leta    = spA % leta  (xi(2))
            lzeta   = spA % lzeta (xi(3))
  
            F = 0.0_RP
            do k = 0, spA % Nz   ; do j = 0, spA % Ny ; do i = 0, spA % Nx
               F = F + self % geom % x(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do
   
            F = F - x
!
!           Stopping criteria: there are several
!           ------------------------------------
            if ( maxval(abs(F)) .lt. TOL ) exit
            if ( abs(xi(1)) .ge. 1.25_RP ) exit
            if ( abs(xi(2)) .ge. 1.25_RP ) exit
            if ( abs(xi(3)) .ge. 1.25_RP ) exit
!
!           Perform a step
!           --------------
            dlxi    = spA % dlxi  (xi(1))
            dleta   = spA % dleta (xi(2))
            dlzeta  = spA % dlzeta(xi(3))

            Jac = 0.0_RP
            do k = 0, spA % Nz   ; do j = 0, spA % Ny ; do i = 0, spA % Nx
               Jac(:,1) = Jac(:,1) + self % geom % x(:,i,j,k) * dlxi(i) * leta(j) * lzeta(k) 
               Jac(:,2) = Jac(:,2) + self % geom % x(:,i,j,k) * lxi(i) * dleta(j) * lzeta(k) 
               Jac(:,3) = Jac(:,3) + self % geom % x(:,i,j,k) * lxi(i) * leta(j) * dlzeta(k) 
            end do               ; end do             ; end do

            dx = solveThreeEquationLinearSystem( Jac , -F )
            xi = xi + STEP * dx
   
         end do

         if ( (abs(xi(1)) .lt. 1.0_RP + INSIDE_TOL) .and. &
              (abs(xi(2)) .lt. 1.0_RP + INSIDE_TOL) .and. &
              (abs(xi(3)) .lt. 1.0_RP + INSIDE_TOL)          ) then
!
!           Solution is valid
!           -----------------
            HexElement_FindPointWithCoords = .true.
   
         else
!
!           Solution is not valid
!           ---------------------
            HexElement_FindPointWithCoords = .false.
         
         end if

      end function HexElement_FindPointWithCoords

      function HexElement_EvaluateSolutionAtPoint(self, spA, xi)
         implicit none
         class(Element),   intent(in)    :: self
         class(NodalStorage), intent(in) :: spA
         real(kind=RP),    intent(in)    :: xi(NDIM)
         real(kind=RP)                   :: HexElement_EvaluateSolutionAtPoint(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: lxi(0:spA % Nx)
         real(kind=RP)  :: leta(0:spA % Ny)
         real(kind=RP)  :: lzeta(0:spA % Nz)
         real(kind=RP)  :: Q(NCONS)
!
!        Compute Lagrange basis
!        ----------------------
         lxi   = spA % lxi(xi(1))
         leta  = spA % leta(xi(2))
         lzeta = spA % lzeta(xi(3))
!
!        Compute the tensor product
!        --------------------------
         Q = 0.0_RP
      
         do k = 0, spA % Nz   ; do j = 0, spA % Ny ; do i = 0, spA % Nx
            Q = Q + self % storage % Q(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
         end do               ; end do             ; end do   

         HexElement_EvaluateSolutionAtPoint = Q

      end function HexElement_EvaluateSolutionAtPoint
      
      END Module ElementClass
