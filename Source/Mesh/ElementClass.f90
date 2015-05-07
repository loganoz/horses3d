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
!       The Quad Element class, Alg. 124. See Sec. 8.2.1.2. This has
!       been modified to add the association of a boundary name to an element
!       edge so that different boundary conditions can be applied to different
!       elements. The names of the boundaries (not necessarily the names of the
!       *boundary conditions* to be applied) are of length BC_STRING_LENGTH.
!       One will associate boundary conditions to boundaries in the routine
!       "ExternalState".
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
      IMPLICIT NONE
      
      
      TYPE Element
          INTEGER                         :: nodeIDs(4)
          INTEGER                         :: N
          TYPE(MappedGeometry)            :: geom
          CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName(4)
      END TYPE Element
      
      INTEGER, PARAMETER :: edgeMap(2,4) = RESHAPE( (/1,2,2,3,4,3,1,4/), (/2,4/) )
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructElement( this, ng, nodeIDs, quadMap )
         USE Nodal2DStorageClass
         IMPLICIT NONE
         
         TYPE(Element)            :: this
         INTEGER                  :: nodeIDs(4)
!
!        --------------------------------------------------------------------------
!        Rem: To avoid interpolating the geometry later in the plotting stage, save
!        the quadMap as part of the element class and use it in the plotter.
!        --------------------------------------------------------------------------
!
         TYPE(TransfiniteQuadMap) :: quadMap
         TYPE(Nodal2DStorage)     :: ng
         
         this%nodeIDs      = nodeIDs
         this%N            = ng%N
         this%boundaryName = "---"
         CALL ConstructMappedGeometry( this%geom, ng, this%N, this%N, quadMap )
         
      END SUBROUTINE ConstructElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetElementBoundaryNames( this, names ) 
         IMPLICIT NONE
         TYPE(Element)                   :: this
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(4)
         this%boundaryName = names
      END SUBROUTINE SetElementBoundaryNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructElement( this )
         IMPLICIT NONE
         TYPE(Element) :: this
         CALL DestructMappedGeometry( this%geom )
      END SUBROUTINE DestructElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintElement( this, id )
         IMPLICIT NONE 
         TYPE(Element) :: this
         INTEGER      :: id
         PRINT *, id, this%nodeIDs, this%boundaryName
      END SUBROUTINE PrintElement
      
      END Module ElementClass
