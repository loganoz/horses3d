!
!////////////////////////////////////////////////////////////////////////
!
!      FaceClass.f
!      Created: 2008-06-05 14:12:52 -0400 
!      By: David Kopriva  
!       
!      Modification history:
!           Modified to 3D             5/27/15, 11:13 AM: David A. Kopriva
!           Added isotropic mortars    4/26/17, 11:12 AM: Andrés Rueda
!           Added anisotropic mortars  5/16/17, 11:11 AM: Andrés Rueda
!
!      Implements Algorithms:
!         Algorithm 125: EdgeClass -> 3D
!
!      A face simply keeps track of which elements share a face and 
!      how they are oriented.
!
!      TODO: Remove projection matrices from each face (a global definition is more efficient)
!////////////////////////////////////////////////////////////////////////
!
      Module FaceClass
      USE SMConstants
      USE MeshTypes
      USE PolynomialInterpAndDerivsModule
      USE GaussQuadrature
      IMPLICIT NONE 
!
!     ----------------------------
!     Mortar Function definition
!     (used for storing the solution and flux information on the face)
!     ----------------------------
!
      TYPE MortarFunction
         REAL(KIND=RP), ALLOCATABLE   :: L(:,:,:)                    ! Solution on the left
         REAL(KIND=RP), ALLOCATABLE   :: R(:,:,:)                    ! Solution on the right
         REAL(KIND=RP), ALLOCATABLE   :: C(:,:,:)                    ! Solution on the mortar
         REAL(KIND=RP), ALLOCATABLE   :: Caux(:,:,:)                 ! ?
      END TYPE MortarFunction 
!
!     ---------------
!     Face definition
!     ---------------
!
      TYPE Face
         INTEGER                         :: FaceType                 ! Type of face: 0 = HMESH_BOUNDARY, 1 = HMESH_INTERIOR
         INTEGER                         :: nodeIDs(4)
         INTEGER                         :: elementIDs(2)            ! Convention is: 1 = Left, 2 = Right
         INTEGER                         :: elementSide(2)
         INTEGER                         :: rotation
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
         
         ! "Mortar" related variables
         TYPE(MortarFunction)            :: Phi                      ! Variable containing the solution information on the face
         REAL(KIND=RP), ALLOCATABLE      :: L2Phi(:,:)               ! Projection matrix from left element
         REAL(KIND=RP), ALLOCATABLE      :: R2Phi(:,:)               ! Projection matrix from right element
         REAL(KIND=RP), ALLOCATABLE      :: Phi2L(:,:)               ! Projection matrix to left element
         REAL(KIND=RP), ALLOCATABLE      :: Phi2R(:,:)               ! Projection matrix to left element
         INTEGER, DIMENSION(2)           :: NL, NR, NPhi             ! Orders of quadrature on left and right elements, and mortar
         INTEGER, DIMENSION(2)           :: NPhiR                    ! Order of quadrature on slave face of mortar (this is needed, since the solution is interpolated before rotating)
      END TYPE Face
!
!     ========
      CONTAINS
!     ========
!
      SUBROUTINE ConstructFace( self, nodeIDs, elementID, side )
         IMPLICIT NONE 
         TYPE(Face) :: self
         INTEGER    :: nodeIDs(4), elementID, side
         self % FaceType       = HMESH_UNDEFINED
         self % nodeIDS        = nodeIDs
         self % elementIDs     = HMESH_NONE
         self % elementSide    = HMESH_NONE
         self % elementIDs(1)  = elementID
         self % elementSide(1) = side
         self % boundaryName   = emptyBCName
         self % rotation       = 0
      END SUBROUTINE ConstructFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructFace( self )
         IMPLICIT NONE 
         TYPE(Face) :: self
         
         IF (self % FaceType == HMESH_INTERIOR) CALL DestructMortarStorage(self)
      END SUBROUTINE DestructFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintFace( self ) 
      IMPLICIT NONE
      TYPE(Face) :: self
      PRINT *, "Face TYPE = "   , self % FaceType
      PRINT *, "Element IDs: "  , self % elementIDs
      PRINT *, "Element Sides: ", self % elementSide
      IF ( self % FaceType == HMESH_INTERIOR )     THEN
         PRINT *, "Neighbor rotation: ", self  %  rotation
      ELSE
         PRINT *, "Boundary name = ", self % boundaryName
      END IF
      PRINT *, "-----------------------------------"
      END SUBROUTINE PrintFace
!
!////////////////////////////////////////////////////////////////////////
!
!  ROUTINE USED TO COMPUTE FACE ROTATION INDEXES
!     This routine takes indexes on the master Face of a mortar and
!     output the corresponding indexes on the slave Face 
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE iijjIndexes(i,j,Nx,Ny,rotation,ii,jj)
      IMPLICIT NONE
      
      INTEGER :: i,j       !<  Input indexes
      INTEGER :: Nx, Ny    !<  Polynomial orders
      INTEGER :: rotation  !<  Face rotation
      INTEGER :: ii,jj     !>  Output indexes
      
      SELECT CASE (rotation)
         CASE (0)
            ii = i
            jj = j
         CASE (1)
            ii = Ny - j
            jj = i
         CASE (2)
            ii = Nx - i
            jj = Ny - j
         CASE (3)
            ii = j
            jj = Nx - i
         CASE (4)
            ii = j
            jj = i
         CASE (5)
            ii = Nx - i
            jj = j
         CASE (6)
            ii = Ny - j
            jj = Nx - i
         CASE (7)
            ii = i
            jj = Ny - j
         CASE DEFAULT 
            PRINT *, "ERROR: Unknown rotation in element faces"
      END SELECT
      
   END SUBROUTINE iijjIndexes
!
!////////////////////////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ConstructMortarStorage( this, Neqn, NGradeqn, elems )
      USE ElementClass
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     Constructor of Mortar Storage:
!        Since the polynomial order is different in every local direction,
!        We will assume that the orientation of the "x" and "y" axes on 
!        the mortar are the same as on the left element's face
!     -------------------------------------------------------------------
!      
      TYPE(Face)   , INTENT(INOUT)      :: this       !<> Current face
      INTEGER      , INTENT(IN)         :: Neqn       !<  Number of equations
      INTEGER      , INTENT(IN)         :: NGradeqn   !<  Number of gradient equations
      TYPE(Element), INTENT(IN), TARGET :: elems(:)   !<  Elements in domain
!
!     --------------------
!     Internal variables  
!     --------------------
!
      
      INTEGER, DIMENSION(2)  :: NL, NR     ! Polynomial orders of left and right element
      INTEGER, DIMENSION(2)  :: NPhi       ! Polynomial orders of mortar 
      INTEGER, DIMENSION(2)  :: NPhiR      ! Polynomial orders of right face of mortar              
      TYPE(Element), POINTER :: eL         ! Element on the "left" of mortar
      TYPE(Element), POINTER :: eR         ! Element on the "right" of mortar
!
!     -----------------------------
!     Basic definitions
!     -----------------------------
!
      eL => elems(this % elementIDs(1)) 
      eR => elems(this % elementIDs(2))
!
!     ----------------------------------------------------------
!     The size of the mortar space is the maximum of each of the
!     contributing faces. Note that the order of operations
!     on a mortar will be to PROJECT and then to ROTATE (unlike the DSEM code).
!     Thus, the components of NPhi(!) will refer to different local face 
!     coordinate directions depending on the orientation.
!     ----------------------------------------------------------
!
      NL(1) = eL % Nxyz (axisMap(1,this % elementSide(1)))
      NL(2) = eL % Nxyz (axisMap(2,this % elementSide(1)))
      
      NR(1) = eR % Nxyz (axisMap(1,this % elementSide(2)))
      NR(2) = eR % Nxyz (axisMap(2,this % elementSide(2)))
      
      SELECT CASE ( this % rotation )
         CASE ( 0, 2, 5, 7 ) ! Local x and y axis are parallel or antiparallel
            NPhi(1)  = MAX(NL(1),NR(1))
            NPhi(2)  = MAX(NL(2),NR(2))
            NPhiR    = NPhi
         CASE ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
            NPhi(1)  = MAX(NL(1),NR(2))
            NPhi(2)  = MAX(NL(2),NR(1))
            NPhiR(1) = NPhi(2)
            NPhiR(2) = NPhi(1)
      END SELECT
      
      ! Save in mortar type
      this % NL    = NL
      this % NR    = NR
      this % NPhi  = NPhi
      this % NPhiR = NPhiR
!
!     -----------------------------
!     Allocate Mortar Storage
!     -----------------------------
!
      ALLOCATE( this%Phi%L    ( Neqn    , 0:NPhi (1), 0:NPhi (2) ) )  
      ALLOCATE( this%Phi%C    ( Neqn    , 0:NPhi (1), 0:NPhi (2) ) )
      ALLOCATE( this%Phi%Caux ( NGradeqn, 0:NPhi (1), 0:NPhi (2) ) )
      ALLOCATE( this%Phi%R    ( Neqn    , 0:NPhiR(1), 0:NPhiR(2) ) ) 
!
!     -----------------
!     Initialize memory
!     -----------------
!
      this % Phi % L    = 0._RP
      this % Phi % R    = 0._RP
      this % Phi % C    = 0._RP
      this % Phi % Caux = 0._RP
!
!     -----------------------------------------------------------------------
!     Construction of the projection matrices (simple Lagrange interpolation)
!        TODO: On Legendre-Gauss-Lobatto this is no longer exact...
!     -----------------------------------------------------------------------
!
      CALL ConstructMortarInterp(this % L2Phi, NL(1), NL(2), NPhi (1), NPhi (2))
      CALL ConstructMortarInterp(this % R2Phi, NR(1), NR(2), NPhiR(1), NPhiR(2))
      CALL ConstructMortar2ElInterp(this % Phi2L, NPhi (1), NPhi (2), NL(1), NL(2))
      CALL ConstructMortar2ElInterp(this % Phi2R, NPhiR(1), NPhiR(2), NR(1), NR(2))
      
   END SUBROUTINE ConstructMortarStorage
!
!////////////////////////////////////////////////////////////////////////
!   
   SUBROUTINE ConstructMortarInterp( Mat, N1x, N1y, N2x, N2y )
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     Constructor of Mortar Interpolation matrices using 0-indexing
!        Obtains a matrix that interpolates from (1) to (2)
!     -------------------------------------------------------------------
!     
      REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)  !> Interpolation matrix to be constructed
      INTEGER                     :: N1x, N1y  !< Polynomial orders on (1)
      INTEGER                     :: N2x, N2y  !< Polynomial orders on (2)
!
!     --------------------
!     Internal variables  
!     --------------------
!
      REAL(KIND=RP) :: x1(0:N1x), y1(0:N1y), w1x(0:N1x), w1y(0:N1y)   ! Nodes and weights
      REAL(KIND=RP) :: x2(0:N2x), y2(0:N2y), w2x(0:N2x), w2y(0:N2y)   ! Nodes and weights
      
      INTEGER       :: i, j, k, l, m, n
!
!     -----------
!     Allocations
!     -----------
!
      ALLOCATE(Mat(0:(N2x+1)*(N2y+1)-1,0:(N1x+1)*(N1y+1)-1))
!
!     ----------------------------------------------------
!     Obtain the quadrature nodes for both faces
!        Currently done only with Legendre-Gauss
!        TODO: Implement with Legendre-Gauss-Lobatto
!     ----------------------------------------------------
!
      CALL GaussLegendreNodesAndWeights(N1x, x1, w1x)
      CALL GaussLegendreNodesAndWeights(N1y, y1, w1y)
      CALL GaussLegendreNodesAndWeights(N2x, x2, w2x)
      CALL GaussLegendreNodesAndWeights(N2y, y2, w2y)
!
!     ----------------------------------------------------
!     Creation of the interpolation matrix
!     ----------------------------------------------------
!
      Mat = 0.0_RP
      
      DO l = 0, N1y
         DO k = 0, N1x
            m = k + l * (N1x+1) ! Column index
            DO j = 0, N2y
               DO i = 0, N2x
                  n = i + j * (N2x+1)       ! Row index
                  
                  Mat(n,m) = LagrangeInterpolationNoBar(x2(i),N1x,x1,k) * LagrangeInterpolationNoBar(y2(j),N1y,y1,l)
               END DO
            END DO
         END DO
      END DO
    
   END SUBROUTINE ConstructMortarInterp
!
!////////////////////////////////////////////////////////////////////////
!   
   SUBROUTINE ConstructMortar2ElInterp( Mat, N1x, N1y, N2x, N2y )
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     Constructor of Mortar Interpolation matrices using 0-indexing
!        Obtains a matrix that interpolates from (1) to (2)
!     -------------------------------------------------------------------
!     
      REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)  !> Interpolation matrix to be constructed
      INTEGER                     :: N1x, N1y  !< Polynomial orders on (1)
      INTEGER                     :: N2x, N2y  !< Polynomial orders on (2)
!
!     --------------------
!     Internal variables  
!     --------------------
!
      REAL(KIND=RP) :: x1(0:N1x), y1(0:N1y), w1x(0:N1x), w1y(0:N1y)   ! Nodes and weights
      REAL(KIND=RP) :: x2(0:N2x), y2(0:N2y), w2x(0:N2x), w2y(0:N2y)   ! Nodes and weights
      REAL(KIND=RP) :: MASSterm                                       ! Mass matrix term (this matrix is diagonal, so we only store one entry at a time)
      
      
      INTEGER       :: i, j, r, s, m, n
!
!     -----------
!     Allocations
!     -----------
!
      ALLOCATE(Mat(0:(N2x+1)*(N2y+1)-1,0:(N1x+1)*(N1y+1)-1))
!
!     ----------------------------------------------------
!     Obtain the quadrature nodes for both faces
!        Currently done only with Legendre-Gauss
!        TODO: Implement with Legendre-Gauss-Lobatto
!     ----------------------------------------------------
!
      CALL GaussLegendreNodesAndWeights(N1x, x1, w1x)
      CALL GaussLegendreNodesAndWeights(N1y, y1, w1y)
      CALL GaussLegendreNodesAndWeights(N2x, x2, w2x)
      CALL GaussLegendreNodesAndWeights(N2y, y2, w2y)
!
!     ----------------------------------------------------
!     Creation of the interpolation matrix
!     ----------------------------------------------------
!
      Mat = 0.0_RP
      
      ! Create S matrix and store it directly in "Mat"
      DO s = 0, N1y
         DO r = 0, N1x
            m = r + s * (N1x+1) ! Column index
            DO j = 0, N2y
               DO i = 0, N2x
                  n = i + j * (N2x+1)       ! Row index
                  
                  Mat(n,m) = LagrangeInterpolationNoBar(x1(r),N2x,x2,i) * LagrangeInterpolationNoBar(y1(s),N2y,y2,j) * &
                                                                                                              w1x(r) * w1y(s)
               END DO
            END DO
         END DO
      END DO
      
      ! Create Mass matrix and finish computing interpolation operator
      DO j = 0, N2y
         DO i = 0, N2x
            n = i + j * (N2x+1)       ! Row index
            
            MASSterm = w2x(i) * w2y(j)
            
            ! Matrix Multiplication I = M⁻¹S (taking advantage of the diagonal matrix)
            Mat(n,:) = Mat(n,:) / MASSterm
         END DO
      END DO
      
      
      
   END SUBROUTINE ConstructMortar2ElInterp
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ProjectToMortar(this, QL, QR, nEqn) 
      IMPLICIT NONE 
!
!     -------------------------------------------------------------------
!     This subroutine projects the value of the borders of the elements
!     to the mortars 
!     -------------------------------------------------------------------
!
!
!     ------
!     Inputs
!     ------
!
      TYPE(Face)   , INTENT(INOUT)     :: this                         !<> Face containing interface information
      REAL(KIND=RP), INTENT(IN)        :: QL (nEqn,0:this % NL(1), &
                                                   0:this % NL(2))     !<  Boundary solution of left element
      REAL(KIND=RP), INTENT(IN)        :: QR (nEqn,0:this % NR(1), &
                                                   0:this % NR(2))     !<  Boundary solution of right element
      INTEGER      , INTENT(IN)        :: nEqn                         !<  Number of equations  
!        
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER               :: iEQ         ! Equation counter
      INTEGER, DIMENSION(2) :: NL, NR      ! Polynomial orders of adjacent elements
      INTEGER, DIMENSION(2) :: NPhi, NPhiR ! Polynomial orders of mortar
      !--------------------------------------------------------------------------------------------
      
      NL   = this % NL
      NR   = this % NR
      NPhi = this % NPhi
      NPhiR= this % NPhiR
!
!     ------------------------------------------------------------------
!     Check if the polynomial orders are the same in both directions 
!     If so, we don't do the polynomial interpolation
!     to increase speed.         
!     ------------------------------------------------------------------
!
      ! Left element
      IF (ALL(NL == NPhi)) THEN
         this % Phi % L(1:nEqn,:,:) = QL(1:nEqn,:,:)
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = QL (iEQ,0:NL(1),0:NL(2))               , &
                                Q2     = this % Phi % L (iEQ,0:NPhi(1),0:NPhi(2)) , &
                                Interp = this % L2Phi                           , &
                                N1x    = NL(1)  , N1y = NL(2)                   , &
                                N2x    = NPhi(1), N2y = NPhi(2)    )
         ENDDO 
      END IF
      
      ! Right element
      IF (ALL(NR == NPhiR)) THEN
         this % Phi % R(1:nEqn,:,:) = QR(1:nEqn,:,:)
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = QR (iEQ,0:NR(1),0:NR(2))                   , &
                                Q2     = this % Phi % R (iEQ,0:NPhiR(1),0:NPhiR(2)) , &
                                Interp = this % R2Phi                               , &
                                N1x    = NR(1)   , N1y = NR(2)                      , &
                                N2x    = NPhiR(1), N2y = NPhiR(2)  )
         ENDDO
      END IF
      
   END SUBROUTINE ProjectToMortar
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ProjectFluxToElement(this,FStarbL,FStarbR,NEqn)
      IMPLICIT NONE
!
!     ---------------------------------------------------------------------
!     Performs the interpolation of the numerical flux from the Face to the
!     elements' FStarb
!     ---------------------------------------------------------------------
!
!     ------
!     Input
!     ------
!      
      TYPE(Face)   , INTENT(INOUT)     :: this                            !<> Face containing the computed flux
      REAL(KIND=RP), INTENT(OUT)       :: FStarbL(NEqn,0:this % NL(1), &
                                                       0:this % NL(2))    !>  Boundary solution of left element
      REAL(KIND=RP), INTENT(OUT)       :: FStarbR(NEqn,0:this % NR(1), &
                                                       0:this % NR(2))    !>  Boundary solution of right element
      INTEGER      , INTENT(IN)        :: nEqn                            !<  Number of equations 
!
!     -------------
!     Local variables
!     -------------
!
      INTEGER               :: iEQ         ! Equation counter
      INTEGER, DIMENSION(2) :: NL, NR      ! Polynomial orders of adjacent elements
      INTEGER, DIMENSION(2) :: NPhi, NPhiR ! Polynomial orders of mortar
      INTEGER               :: i,j,ii,jj   ! Counters
      INTEGER               :: rotation    ! Face rotation
!
!     ----------------------------------
!     Get polynomial orders and rotation
!     ----------------------------------
!
      NL   = this % NL
      NR   = this % NR
      NPhi = this % NPhi
      NPhiR= this % NPhiR
      
      rotation = this % rotation
!
!     ---------------------------------------------------------
!     Store flux in Phi%R taking into account rotation and sign
!        (both rotation and sign in Phi%L are the same as in Phi%C)
!     ---------------------------------------------------------
!
      DO j = 0, NPhi(2)
         DO i = 0, NPhi(1)
            CALL iijjIndexes(i,j,NPhi(1),NPhi(2),rotation,ii,jj)                              ! This turns according to the rotation of the elements
            this % Phi % R(:,ii,jj) = - this % Phi % C(:,i,j)
         END DO   
      END DO 
!
!     -----------------------------------------------------------------
!     Project flux back to element
!        Check if the polynomial orders are the same in both directions 
!        If so, we don't do the polynomial interpolation
!        to increase speed.
!     -----------------------------------------------------------------
!
      ! Left element
      IF (ALL(NL == NPhi)) THEN
         FStarbL(1:nEqn,:,:) = this % Phi % C(1:nEqn,:,:)               ! Phi%C is used instead on Phi%L to avoid the copying operation
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = this % Phi % C (iEQ,0:NPhi(1),0:NPhi(2)) , &   ! Phi%C is used instead on Phi%L to avoid the copying operation
                                Q2     = FStarbL(iEQ,0:NL(1),0:NL(2))             , &
                                Interp = this % Phi2L                             , &
                                N1x    = NPhi(1), N1y = NPhi(2)                   , &
                                N2x    = NL(1)  , N2y = NL(2)     )
         ENDDO 
      END IF
      
      ! Right element
      IF (ALL(NR == NPhiR)) THEN
         FStarbR(1:nEqn,:,:) = this % Phi % R(1:nEqn,:,:)
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = this % Phi % R (iEQ,0:NPhiR(1),0:NPhiR(2)) , &
                                Q2     = FStarbR(iEQ,0:NR(1),0:NR(2))               , &
                                Interp = this % Phi2R                               , &
                                N1x    = NPhiR(1), N1y = NPhiR(2)                   , &
                                N2x    = NR(1)   , N2y = NR(2)    )
         ENDDO
      END IF
      
   END SUBROUTINE ProjectFluxToElement
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ProjectToElement(this,C,UL,UR,NEqn)
      IMPLICIT NONE
!
!     ---------------------------------------------------------------------
!     Performs the interpolation of the numerical flux from the Face to the
!     elements' FStarb
!     ---------------------------------------------------------------------
!
!     ------
!     Input
!     ------
!      
      TYPE(Face)   , INTENT(INOUT)     :: this                             !<> Face 
      REAL(KIND=RP), INTENT(OUT)       :: C  (NEqn,0:this % NPhi(1), &
                                                   0:this % NPhi(2))       !>  Variable to be projected
      REAL(KIND=RP), INTENT(OUT)       :: UL (NEqn,0:this % NL(1)  , &
                                                   0:this % NL(2))         !>  Boundary solution of left element
      REAL(KIND=RP), INTENT(OUT)       :: UR (NEqn,0:this % NR(1)  , &
                                                   0:this % NR(2))         !>  Boundary solution of right element
      INTEGER      , INTENT(IN)        :: nEqn                             !<  Number of equations 
!
!     -------------
!     Local variables
!     -------------
!
      REAL(KIND=RP), ALLOCATABLE :: R (:,:,:)   ! Variable projected on the right element
      INTEGER                    :: iEQ         ! Equation counter
      INTEGER, DIMENSION(2)      :: NL, NR      ! Polynomial orders of adjacent elements
      INTEGER, DIMENSION(2)      :: NPhi, NPhiR ! Polynomial orders of mortar
      INTEGER                    :: i,j,ii,jj   ! Counters
      INTEGER                    :: rotation    ! Face rotation
      
!
!     ----------------------------------
!     Get polynomial orders and rotation
!     ----------------------------------
!
      NL   = this % NL
      NR   = this % NR
      NPhi = this % NPhi
      NPhiR= this % NPhiR
      
      rotation = this % rotation
      ALLOCATE(R(NEqn,0:NPhiR(1),0:NPhiR(2)))
!
!     ---------------------------------------------------------
!     Store flux in Phi%R taking into account rotation (I*C= UL)
!     ---------------------------------------------------------
!
      DO j = 0, NPhi(2)
         DO i = 0, NPhi(1)
            CALL iijjIndexes(i,j,NPhi(1),NPhi(2),rotation,ii,jj)                              ! This turns according to the rotation of the elements
            R(1:NEqn,ii,jj) = C(1:NEqn,i,j)
         END DO   
      END DO 
!
!     -----------------------------------------------------------------
!     Project flux back to element
!        Check if the polynomial orders are the same in both directions 
!        If so, we don't do the polynomial interpolation
!        to increase speed.
!     -----------------------------------------------------------------
!
      ! Left element
      IF (ALL(NL == NPhi)) THEN
         UL(1:NEqn,:,:) = C(1:NEqn,:,:)                                           ! C is used instead on L to avoid the copying operation
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = C (iEQ,0:NPhi(1),0:NPhi(2))   , &   ! C is used instead on L to avoid the copying operation
                                Q2     = UL(iEQ,0:NL(1),0:NL(2))       , &
                                Interp = this % Phi2L                  , &
                                N1x    = NPhi(1), N1y = NPhi(2)        , &
                                N2x    = NL(1)  , N2y = NL(2)  )
         ENDDO 
      END IF
      
      ! Right element
      IF (ALL(NR == NPhiR)) THEN
         UR(1:NEqn,:,:) = R(1:NEqn,:,:)
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = R (iEQ,0:NPhiR(1),0:NPhiR(2)), &
                                Q2     = UR(iEQ,0:NR   (1),0:NR   (2)), &
                                Interp = this % Phi2R                 , &
                                N1x    = NPhiR(1), N1y = NPhiR(2)     , &
                                N2x    = NR(1)   , N2y = NR(2)    )
         ENDDO
      END IF
      
      DEALLOCATE(R)
      
   END SUBROUTINE ProjectToElement
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Project1Eqn(Q1,Q2,Interp,N1x,N1y,N2x,N2y)
      IMPLICIT NONE
!
!     --------------------------------------------------------------------
!     This subroutine interpolates the boundary solution of an element
!     onto a face taking advantage of how fortran handles multidimensional
!     arrays 
!     --------------------------------------------------------------------
!
!
!     ------
!     Inputs
!     ------
!
      INTEGER        :: N1x, N1y, N2x, N2y             !<  Polynomial orders
      REAL(KIND=RP)  :: Q1(0:(N1x+1)*(N1y+1)-1)        !<  Solution to be interpolated (grid (1))
      REAL(KIND=RP)  :: Q2(0:(N2x+1)*(N2y+1)-1)        !>  Interpolated solution      (grid (2))
      REAL(KIND=RP)  :: Interp(0:(N2x+1)*(N2y+1)-1, & 
                               0:(N1x+1)*(N1y+1)-1)    !<  Interpolation matrix
      
      Q2 = MATMUL(Interp,Q1)
      
   END SUBROUTINE Project1Eqn
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE DestructMortarStorage( this )
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     Destructor of Mortar Storage
!     -------------------------------------------------------------------
!         
      TYPE(Face)       :: this
      !------------------------------------------
      CALL DestructMortarFunction(this%Phi)      
       
      DEALLOCATE(this%L2Phi)   
      DEALLOCATE(this%R2Phi)
      DEALLOCATE(this%Phi2R)
      DEALLOCATE(this%Phi2L)
    
   END SUBROUTINE DestructMortarStorage
!
!////////////////////////////////////////////////////////////////////////
!   
   SUBROUTINE DestructMortarFunction( this)
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     Destructor of Mortar Function
!     -------------------------------------------------------------------
!         
      TYPE(mortarFunction) :: this
      !-------------------------------------------
      DEALLOCATE(this%L)
      DEALLOCATE(this%R)
      DEALLOCATE(this%C)
      DEALLOCATE(this%Caux)
    
   END SUBROUTINE DestructMortarFunction   
!
!////////////////////////////////////////////////////////////////////////
!
END Module FaceClass
