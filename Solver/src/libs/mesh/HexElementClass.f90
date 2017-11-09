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
      
      END Module ElementClass

