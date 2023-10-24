!
!////////////////////////////////////////////////////////////////////////
!
!!     This module defines the connectivity in the master element.
!
!      PUBLIC DATA:
!        INTEGER, PARAMETER :: NODES_PER_FACE, NODES_PER_ELEMENT,
!                              FACES_PER_ELEMENT, EDGES_PER_ELEMENT
!        INTEGER, DIMENSION(4,6)  :: localFaceNode
!        INTEGER, DIMENSION(2,12) :: localEdgeNode
!
!!
! /////////////////////////////////////////////////////////////////////
!
      MODULE ElementConnectivityDefinitions
      USE SMConstants
      
      IMPLICIT NONE

      private
      public   NODES_PER_FACE, NODES_PER_ELEMENT, FACES_PER_ELEMENT, EDGES_PER_ELEMENT
      public   localFaceNode, localEdgeNode, sideMap, neighborFaces, axisMap, normalAxis

      INTEGER, PARAMETER :: NODES_PER_FACE    = 4
      INTEGER, PARAMETER :: NODES_PER_ELEMENT = 8
      INTEGER, PARAMETER :: FACES_PER_ELEMENT = 6
      INTEGER, PARAMETER :: EDGES_PER_ELEMENT = 12
 
!
!     -------------------------------------------------------------------------
!!    axisMap gives the element local coordinate number for the two directions
!!    on each face. The coordinate numbers are given by (xi,eta,zeta) = (1,2,3).
!!    For instance, the two coordinate directions on Face 1 are (xi,zeta).
!     -------------------------------------------------------------------------
!
      integer, dimension(2,6) :: axisMap =                        &
                                 RESHAPE( (/1, 3,                 & ! Face 1 (x,z)
                                            1, 3,                 & ! Face 2 (x,z)
                                            1, 2,                 & ! Face 3 (x,y)
                                            2, 3,                 & ! Face 4 (y,z)
                                            1, 2,                 & ! Face 5 (x,y)
                                            2, 3/)                & ! Face 6 (y,z)
                                 ,(/2,6/))
!
!     -------------------------------------------------------------
!     Definition of the normal axis to a face
!        Positive for axis going out and negative for axis going in
!     -------------------------------------------------------------
!
      integer, parameter :: normalAxis(6) = [ -2,  2, -3,  1,  3, -1]
!
!-------------------------------------------------------------------------
!  Definition of the neighbor faces to a given face
!-------------------------------------------------------------------------
!      
      integer, parameter :: neighborFaces(4,6) = reshape (  (/ 3, 4, 5, 6, &
                                                               3, 4, 5, 6, &
                                                               1, 2, 4, 6, &
                                                               1, 2, 3, 5, &
                                                               1, 2, 4, 6, &
                                                               1, 2, 3, 5 /) , (/4,6/) )

!
!-------------------------------------------------------------------------
!  Definition of the local node numbers for the faces of the master element
!--------------------------------------------------------------------------
!
!
      INTEGER, DIMENSION(4,6) :: localFaceNode = &
     &         RESHAPE( (/                    & 
     &         1, 2, 6, 5,                    & ! Node numbers, face 1 FRONT
     &         4, 3, 7, 8,                    & ! Node numbers, face 2 BACK
     &         1, 2, 3, 4,                    & ! Node numbers, face 3 BOTTOM
     &         2, 3, 7, 6,                    & ! Node numbers, face 4 RIGHT
     &         5, 6, 7, 8,                    & ! Node numbers, face 5 TOP
     &         1, 4, 8, 5                     & ! Node numbers, face 6 LEFT
     &         /),(/4,6/))
!
!-------------------------------------------------------------------------
!  Definition of the local node numbers for the edges of the master element
!--------------------------------------------------------------------------
!
      INTEGER, DIMENSION(2,12) :: localEdgeNode = &
     &         RESHAPE( (/                   &
     &         1, 2,                &! Node numbers, edge 1
     &         2, 6,                &! Node numbers, edge 2
     &         5, 6,                &! Node numbers, edge 3
     &         1, 5,                &! Node numbers, edge 4
     &         4, 3,                &! Node numbers, edge 5
     &         3, 7,                &! Node numbers, edge 6
     &         8, 7,                &! Node numbers, edge 7
     &         4, 8,                &! Node numbers, edge 8
     &         1, 4,                &! Node numbers, edge 9
     &         2, 3,                &! Node numbers, edge 10
     &         6, 7,                &! Node numbers, edge 11
     &         5, 8                 &! Node numbers, edge 12
     &         /), (/2,12/) )


!----------------------------------------------------------------------+
!                                                                      |
!  ELEMENT GEOMETRY, 8 NODES                                           |
!  ------------------------------------------------------              |
!                                                                      |
!                                                                      |
!                                   4 --------------- 3                |
!                               -                  -                   |
!                            -     |           -     |                 |
!                         -        |        -        |                 |
!                     -       (2)  |     -           |                 |
!                  -               |  -    (3)       |                 |
!                -                 -                 |                 |
!             8 --------------- 7                    |                 |
!                                  |    (4)          |                 |
!             |        (6)      |                                      |
!             |                 |  1 --------------- 2                 |
!             |                 |                 -       ETA (5)      |
!             |       (5)    -  |              -           |           |
!             |           -     |  (1)      -              |           |
!             |        -        |        -                 |   XI (4)  |
!             |     -           |     -                   /-------     |
!                -                 -                     /             |
!             5 --------------- 6                       / ZETA (6)     |
!----------------------------------------------------------------------|
!
! The following definitions are for the PATRAN defined HEX elements
!
      REAL(KIND=RP), PARAMETER :: oth = 1._RP/3._RP
      REAL(KIND=RP), PARAMETER :: tth = 2._RP/3._RP
!
!    ------------------------------------------------------------------
!     Mapping of element ordered nodes to the 6 faces of a HEX8
!     Mapping of the nodes on a face for the face interpolant of a HEX8
!    ------------------------------------------------------------------
!
      INTEGER, DIMENSION(4,6) :: faceMapHex8 = &
      RESHAPE((/1,2,6,5, &
                4,3,7,8, &
                1,2,3,4, &
                2,3,7,6, &
                5,6,7,8, &
                1,4,8,5/),(/4,6/))
!
      INTEGER, DIMENSION(2,2,6) :: intrpMapHex8 = &
      RESHAPE((/1,2,5,6, &
                4,3,8,7, &
                1,2,4,3, &
                2,3,6,7, &
                5,6,8,7, &
                1,4,5,8/),(/2,2,6/))
!
!    ------------------------------------------------------------------
!     Mapping of element ordered nodes to the 6 faces of a HEX26
!     Mapping of the nodes on a face for the face interpolant of a HEX26
!    ------------------------------------------------------------------
!
      INTEGER, DIMENSION(9,6) :: faceMapHex26 =   &
      RESHAPE((/1,2,6,5,9,14,17,13,25,            &
                4,3,7,8,11,15,19,16,26,           &
                1,2,3,4,9,10,11,12,21,            &
                2,3,7,6,10,15,18,14,24,           &
                5,6,7,8,17,18,19,20,22,           &
                1,4,8,5,12,16,20,13,23/),(/9,6/))

      INTEGER, DIMENSION(3,3,6) :: intrpMapHex26 = &
      RESHAPE((/1,9,2,13,25,14,5,17,6,             &
                4,11,3,16,26,15,8,19,7,              &
                1,9,2,12,21,10,4,11,3,               &
                2,10,3,14,24,15,6,18,7,              &
                5,17,6,20,22,18,8,19,7,              &
                1,12,4,13,23,16,5,20,8/),(/3,3,6/))
!
!    ------------------------------------------------------------------
!     Mapping of element ordered nodes to the 6 faces of a HEX64
!     Mapping of the nodes on a face for the face interpolant of a HEX64
!    ------------------------------------------------------------------
!
!
      INTEGER, DIMENSION(16,6) :: faceMapHex64 = &
      RESHAPE((/1,2,6,5,9,10,18,22,26,25,21,17,0,0,0,0,    &
                4,3,7,8,14,13,19,23,29,30,24,20,0,0,0,0,   &
                1,2,3,4,9,10,11,12,13,14,15,16,0,0,0,0,    &
                2,3,7,6,11,12,19,23,28,27,22,18,0,0,0,0,   &
                5,6,7,8,25,26,27,28,29,30,31,32,0,0,0,0,   &
                1,4,8,5,16,15,20,24,31,32,21,17,0,0,0,0/),(/16,6/))

      INTEGER, DIMENSION(4,4,6) :: intrpMapHex64 = &
      RESHAPE((/1,2,6,5,9,10,18,22,26,25,21,17,0,0,0,0,    &
                4,3,7,8,14,13,19,23,29,30,24,20,0,0,0,0,   &
                1,2,3,4,9,10,11,12,13,14,15,16,0,0,0,0,    &
                2,3,7,6,11,12,19,23,28,27,22,18,0,0,0,0,   &
                5,6,7,8,25,26,27,28,29,30,31,32,0,0,0,0,   &
                1,4,8,5,16,15,20,24,31,32,21,17,0,0,0,0/),(/4,4,6/))
!
      REAL(KIND=RP), DIMENSION(2) :: xiDistributionHex8  = (/0.0_RP, 1.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: xiDistributionHex26 = (/0.0_RP,0.5_RP, 1.0_RP/)
      REAL(KIND=RP), DIMENSION(4) :: xiDistributionHex64 = (/0.0_RP,oth,tth, 1.0_RP/)
!
!     -------------------------
!     list of defined constants
!     -------------------------
!
      INTEGER, PARAMETER :: ILLEG = -99, DFLT = 0, B1F2 = 1, B1B2 = 2, &
                            F1B2  = 3  , F2F1 = 4, B2F1 = 5, B2B1 = 6, &
                            F2B1  = 7
!
!     ------------------------------------------------------------------------
!     Orientation map: This is a 4x4 array that determines how a slave 
!     face is mapped onto a mortar. The value 1 or 2 indicates the index
!     while the sign indicates forward or backward. This array, then assigns
!     the constants defined above to the particular copy operation that must 
!     be done. For instance, (index1,index2) of the mortar mates with
!      (index2,-index1) of the slave, then the orientation is assigned the 
!     value F2B1 for "forward 2,backward 1".
!     -------------------------------------------------------------------------
!
      INTEGER, PARAMETER, DIMENSION(-2:2,-2:2) :: orientMap =               &
                          RESHAPE((/ILLEG,B1B2 ,ILLEG,F1B2 ,ILLEG,          &
                                    B2B1 ,ILLEG,ILLEG,ILLEG,F2B1,           &
                                    ILLEG,ILLEG,ILLEG,ILLEG,ILLEG,          &
                                    B2F1 ,ILLEG,ILLEG,ILLEG,F2F1,           &
                                    ILLEG,B1F2 ,ILLEG,DFLT ,ILLEG/),(/5,5/))
!
!     -------------------------------------------------------------------------
!     Side Map: Given two corner points, return the index in which an array 
!     runs between these two corners. This is shown in the figure below:
!
!               4 -----------3-------------3
!               |                          |
!          /\   |                          |
!          |    |                          |
!               |                          |
!          2    |                          |
!               4                          2
!          x    |                          |
!          e    |                          |
!          d    |                          |
!          n    |                          |
!          i    |                          |
!               |                          |
!               1 -----------1-------------2
!                        index 1 -->
!     -----------------------------------------------------------------------------
!
      INTEGER, PARAMETER, DIMENSION(4,4) :: sideMap =                           &
                                 RESHAPE((/ILLEG, -1   , ILLEG, -2   ,          &
                                           1    , ILLEG, -2   , ILLEG,          &
                                           ILLEG, 2    , ILLEG, 1    ,          &
                                           2    , ILLEG, -1   , ILLEG/),(/4,4/))

      END MODULE ElementConnectivityDefinitions