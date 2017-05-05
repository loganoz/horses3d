!
!////////////////////////////////////////////////////////////////////////
!
!      FaceClass.f
!      Created: 2008-06-05 14:12:52 -0400 
!      By: David Kopriva  
!
!      Modified to 3D 5/27/15, 11:13 AM
!
!      Implements Algorithms:
!         Algorithm 125: EdgeClass -> 3D
!
!      A face simply keeps track of which elements share a face and 
!      how they are oriented.
!
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
         INTEGER                         :: NL, NR, N                ! Order of approximation
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
!     TODO: Check if this is enough or if one needs 8 indexes!!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE iijjIndexes(i,j,Nx,Ny,rotation,ii,jj)
      IMPLICIT NONE
      
      INTEGER :: i,j       !<  Input indexes
      INTEGER :: Nx,Ny     !<  Polynomial orders
      INTEGER :: rotation  !<  Face rotation
      INTEGER :: ii,jj     !>  Output indexes
      
      SELECT CASE (rotation)
         CASE (0)
            ii = i
            jj = j
         CASE (1) ! is this o.k?
            ii = Nx - i
            jj = j
         CASE (2)
            ii = Nx - i
            jj = Ny - j
         CASE (3)
            ii = j
            jj = Nx - i
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
!     Constructor of Mortar Storage
!     -------------------------------------------------------------------
!      
      TYPE(Face)   , INTENT(INOUT) :: this       !<> Current face
      INTEGER      , INTENT(IN)    :: Neqn       !<  Number of equations
      INTEGER      , INTENT(IN)    :: NGradeqn   !<  Number of gradient equations
      TYPE(Element), INTENT(IN)    :: elems(:)   !<  Elements in domain
!
!     --------------------
!     Internal variables  
!     
!     --------------------
!
      
      INTEGER  :: NL, NR     ! Polynomial orders of left and right element (TODO: NLx, NLy, NRx, NRy)
      INTEGER  :: NPhi       ! Polynomial order of mortar                  (TODO: Nphix, NPhiy)
!
!     -----------------------------
!     Allocate Mortar Storage
!     -----------------------------
!
      NL = elems(this % elementIDs(1)) % N
      NR = elems(this % elementIDs(2)) % N
      NPhi = MAX(NL,NR)
      
      this % NL = NL
      this % NR = NR
      this % N  = NPhi
      
      ALLOCATE( this%Phi%L    ( Neqn    , 0:NPhi, 0:NPhi ) )
      ALLOCATE( this%Phi%R    ( Neqn    , 0:NPhi, 0:NPhi ) )   
      ALLOCATE( this%Phi%C    ( Neqn    , 0:NPhi, 0:NPhi ) )
      ALLOCATE( this%Phi%Caux ( NGradeqn, 0:NPhi, 0:NPhi ) )
!
!     -----------------
!     Initialize memory
!     -----------------
!
      this % Phi % L = 0._RP
      this % Phi % R = 0._RP
      this % Phi % C = 0._RP
      this % Phi % Caux = 0._RP
!
!     -----------------------------------------------------------------------
!     Construction of the projection matrices (simple Lagrange interpolation)
!        TODO: a. Check if it changes a lot using finest grid for quadrature integral
!              b. On Legendre-Gauss-Lobatto this is no longer exact...
!     -----------------------------------------------------------------------
!
      CALL ConstructMortarInterp(this % L2Phi, NL, NL, NPhi, NPhi)
      CALL ConstructMortarInterp(this % R2Phi, NR, NR, NPhi, NPhi)
      CALL ConstructMortarInterp(this % Phi2L, NPhi, NPhi, NL, NL)
      CALL ConstructMortarInterp(this % Phi2R, NPhi, NPhi, NR, NR)
      
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
      TYPE(Face)   , INTENT(INOUT)     :: this         !<> Face containing interface information
      REAL(KIND=RP), INTENT(IN)        :: QL(:,0:,0:)  !<  Boundary solution of left element
      REAL(KIND=RP), INTENT(IN)        :: QR(:,0:,0:)  !<  Boundary solution of right element
      INTEGER      , INTENT(IN)        :: nEqn         !<  Number of equations  
!        
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER   :: iEQ                 ! Equation counter
      INTEGER   :: NLx, NLy, NRx, NRy  ! Polynomial orders of adjacent elements
      INTEGER   :: Nx, Ny              ! Polynomial orders of mortar
      !--------------------------------------------------------------------------------------------
      
      NLx  = this % NL        ! TODO: Change when anisotropic polynomials are implemented
      NLy  = this % NL
      NRx  = this % NR
      NRy  = this % NR
      Nx   = this % N
      Ny   = this % N
!
!     ------------------------------------------------------------------
!     Check if the polynomial orders are the same in both directions 
!     If so, we don't do the polynomial interpolation
!     to increase speed.         
!     ------------------------------------------------------------------
!
      ! Left element
      IF (NLx == Nx .AND. NLy == Ny) THEN
         this % Phi % L = QL
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = QL (iEQ,:,:)             , &
                                Q2     = this % Phi % L (iEQ,:,:) , &
                                Interp = this % L2Phi             , &
                                N1x    = NLx  , N1y = NLy         , &
                                N2x    = Nx   , N2y = Ny    )
         ENDDO 
      END IF
      
      ! Right element
      IF (NRx == Nx .AND. NRy == Ny) THEN
         this % Phi % R = QR
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = QR (iEQ,:,:)             , &
                                Q2     = this % Phi % R (iEQ,:,:) , &
                                Interp = this % R2Phi             , &
                                N1x    = NRx  , N1y = NRy         , &
                                N2x    = Nx   , N2y = Ny   )
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
      TYPE(Face)   , INTENT(INOUT)     :: this              !<> Face containing the computed flux
      REAL(KIND=RP), INTENT(OUT)       :: FStarbL(:,0:,0:)  !>  Boundary solution of left element
      REAL(KIND=RP), INTENT(OUT)       :: FStarbR(:,0:,0:)  !>  Boundary solution of right element
      INTEGER      , INTENT(IN)        :: nEqn              !<  Number of equations 
!
!     -------------
!     Local variables
!     -------------
!
      INTEGER   :: iEQ                 ! Equation counter
      INTEGER   :: NLx, NLy, NRx, NRy  ! Polynomial orders of adjacent elements
      INTEGER   :: Nx, Ny              ! Polynomial orders of mortar
      INTEGER   :: i,j,ii,jj           ! Counters
      INTEGER   :: rotation            ! Face rotation
!
!     ----------------------------------
!     Get polynomial orders and rotation
!     ----------------------------------
!
      NLx  = this % NL        ! TODO: Change when anisotropic polynomials are implemented
      NLy  = this % NL
      NRx  = this % NR
      NRy  = this % NR
      Nx   = this % N
      Ny   = this % N
      
      rotation = this % rotation
!
!     ---------------------------------------------------------
!     Store flux in Phi%R taking into account rotation and sign
!        (both rotation and sign in Phi%L are the same as in Phi%C)
!     ---------------------------------------------------------
!
      DO j = 0, Ny
         DO i = 0, Nx
            CALL iijjIndexes(i,j,Nx,Ny,rotation,ii,jj)                              ! This turns according to the rotation of the elements
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
      IF (NLx == Nx .AND. NLy == Ny) THEN
         FStarbL = this % Phi % C                                       ! Phi%C is used instead on Phi%L to avoid the copying operation
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = this % Phi % C (iEQ,:,:) , &   ! Phi%C is used instead on Phi%L to avoid the copying operation
                                Q2     = FStarbL(iEQ,:,:)         , &
                                Interp = this % Phi2L             , &
                                N1x    = Nx   , N1y = Ny          , &
                                N2x    = NLx  , N2y = NLy    )
         ENDDO 
      END IF
      
      ! Right element
      IF (NRx == Nx .AND. NRy == Ny) THEN
         FStarbR = this % Phi % R
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = this % Phi % R (iEQ,:,:) , &
                                Q2     = FStarbR(iEQ,:,:)         , &
                                Interp = this % Phi2R             , &
                                N1x    = Nx   , N1y = Ny          , &
                                N2x    = NRx  , N2y = NRy   )
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
      TYPE(Face)   , INTENT(INOUT)     :: this         !<> Face 
      REAL(KIND=RP), INTENT(OUT)       :: C (:,0:,0:)  !>  Variable to be projected
      REAL(KIND=RP), INTENT(OUT)       :: UL(:,0:,0:)  !>  Boundary solution of left element
      REAL(KIND=RP), INTENT(OUT)       :: UR(:,0:,0:)  !>  Boundary solution of right element
      INTEGER      , INTENT(IN)        :: nEqn              !<  Number of equations 
!
!     -------------
!     Local variables
!     -------------
!
      INTEGER   :: iEQ                 ! Equation counter
      INTEGER   :: NLx, NLy, NRx, NRy  ! Polynomial orders of adjacent elements
      INTEGER   :: Nx, Ny              ! Polynomial orders of mortar
      INTEGER   :: i,j,ii,jj           ! Counters
      INTEGER   :: rotation            ! Face rotation
!
!     ----------------------------------
!     Get polynomial orders and rotation
!     ----------------------------------
!
      NLx  = this % NL        ! TODO: Change when anisotropic polynomials are implemented
      NLy  = this % NL
      NRx  = this % NR
      NRy  = this % NR
      Nx   = this % N
      Ny   = this % N
      
      rotation = this % rotation
!
!     ---------------------------------------------------------
!     Store flux in Phi%R taking into account rotation (I*C= UL)
!     ---------------------------------------------------------
!
      DO j = 0, Ny
         DO i = 0, Nx
            CALL iijjIndexes(i,j,Nx,Ny,rotation,ii,jj)                              ! This turns according to the rotation of the elements
            this % Phi % R(:,ii,jj) = C(:,i,j)
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
      IF (NLx == Nx .AND. NLy == Ny) THEN
         UL = C                                           ! Phi%C is used instead on Phi%L to avoid the copying operation
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = C (iEQ,:,:)     , &   ! Phi%C is used instead on Phi%L to avoid the copying operation
                                Q2     = UL(iEQ,:,:)     , &
                                Interp = this % Phi2L    , &
                                N1x    = Nx   , N1y = Ny , &
                                N2x    = NLx  , N2y = NLy    )
         ENDDO 
      END IF
      
      ! Right element
      IF (NRx == Nx .AND. NRy == Ny) THEN
         UR = this % Phi % R
      ELSE
         DO iEQ = 1, nEqn
            CALL Project1Eqn  ( Q1     = this % Phi % R (iEQ,:,:) , &
                                Q2     = UR(iEQ,:,:)              , &
                                Interp = this % Phi2R             , &
                                N1x    = Nx   , N1y = Ny          , &
                                N2x    = NRx  , N2y = NRy   )
         ENDDO
      END IF
      
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
