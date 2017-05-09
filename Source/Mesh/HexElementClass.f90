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
      IMPLICIT NONE
      
      
      TYPE Element
          INTEGER                                        :: nodeIDs(8)
          INTEGER                                        :: Nx, Ny, Nz
          TYPE(MappedGeometry)                           :: geom
          REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: Q, QDot, G
          REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: U_x, U_y, U_z
!
!         -------------------------------------------------------------
!         Boundary values of: The solution, the inviscid Riemann flux, 
!         the viscous riemann flux
!         -------------------------------------------------------------
!
          REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: Qb, Ub, U_xb, U_yb, U_zb, FStarb
          CHARACTER(LEN=BC_STRING_LENGTH)                :: boundaryName(6)
          CHARACTER(LEN=BC_STRING_LENGTH)                :: boundaryType(6)
          INTEGER                                        :: NumberOfConnections(6)
          TYPE(Connectivity)                             :: Connection(6)
      END TYPE Element 
            
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructElementGeometry( self, ng, nodeIDs, hexMap )
         USE NodalStorageClass
         IMPLICIT NONE
         
         TYPE(Element)           :: self
         TYPE(NodalStorage)      :: ng
         INTEGER                 :: nodeIDs(8)
         TYPE(TransfiniteHexMap) :: hexMap
         
         self % nodeIDs               = nodeIDs
         self % Nx                    = ng % Nx
         self % Ny                    = ng % Ny
         self % Nz                    = ng % Nz
         self % boundaryName          = emptyBCName
         self % boundaryType          = emptyBCName
!
!        --------
!        Geometry
!        --------
!
         CALL ConstructMappedGeometry( self % geom, ng, hexMap )
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
         TYPE(Element) :: self
         INTEGER       :: Nx, Ny, Nz, nEqn, nGradEqn
         LOGICAL       :: flowIsNavierStokes
!
!        ----------------
!        Volume variables
!        ----------------
!
         ALLOCATE( self % Q   (0:Nx,0:Ny,0:Nz,nEqn) )
         ALLOCATE( self % QDot(0:Nx,0:Ny,0:Nz,nEqn) )
         ALLOCATE( self % G   (0:Nx,0:Ny,0:Nz,nEqn) )
         
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( self % U_x(0:Nx,0:Ny,0:Nz,nGradEqn) )
            ALLOCATE( self % U_y(0:Nx,0:Ny,0:Nz,nGradEqn) )
            ALLOCATE( self % U_z(0:Nx,0:Ny,0:Nz,nGradEqn) )
         END IF
!
!        ---------------
!        Boundary values
!        ---------------
!
         ! Temporarily allocating with maximum (TODO: this is not very efficient and has to be changed) DGBoundaryStorage TYPE!!
         Nmax = MAX(Nx,Ny,Nz)
         ALLOCATE( self % Qb    (nEqn,0:Nmax,0:Nmax,6) )
         ALLOCATE( self % FStarb(nEqn,0:Nmax,0:Nmax,6) )
         
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( self % U_xb(nGradEqn,0:Nmax,0:Nmax,6) )
            ALLOCATE( self % U_yb(nGradEqn,0:Nmax,0:Nmax,6) )
            ALLOCATE( self % U_zb(nGradEqn,0:Nmax,0:Nmax,6) )
            ALLOCATE( self % Ub  (nGradEqn,0:Nmax,0:Nmax,6) )
         END IF
!
!        -----------------
!        Initialize memory
!        -----------------
!
         self % G           = 0.0_RP
         self % Q           = 0.0_RP
         self % QDot        = 0.0_RP
         self % Qb          = 0.0_RP
         self % FStarb      = 0.0_RP
      
         IF ( flowIsNavierStokes )     THEN
            self % Ub          = 0.0_RP
            self % U_x         = 0.0_RP
            self % U_y         = 0.0_RP
            self % U_z         = 0.0_RP
            self % U_xb        = 0.0_RP
            self % U_yb        = 0.0_RP
            self % U_zb        = 0.0_RP
         END IF

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
         
         DEALLOCATE( self % Q, self % QDot, self % G )
         DEALLOCATE( self % Qb, self % FStarb )
         
         IF ( ALLOCATED(self % Ub) )     THEN
            DEALLOCATE( self % Ub, self % U_x, self % U_y, self % U_xb, self % U_yb )
         END IF

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
         
         WRITE(funit) self % Q
      
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
         
         READ(funit) self % Q
      
      END SUBROUTINE LoadSolutionFromUnit
      
      END Module ElementClass

