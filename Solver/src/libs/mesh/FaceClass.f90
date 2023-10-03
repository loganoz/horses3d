!!
!////////////////////////////////////////////////////////////////////////
!
!      Modification history:
!           Modified to 3D             5/27/15, 11:13 AM: David A. Kopriva
!           Added isotropic mortars    4/26/17, 11:12 AM: Andrés Rueda
!           Added anisotropic mortars  5/16/17, 11:11 AM: Andrés Rueda
!           Embedded mortars          11/13/17, 05:37 PM: Juan Manzanero
!
!////////////////////////////////////////////////////////////////////////
!
      Module FaceClass
         use SMConstants
         use MeshTypes
         use PolynomialInterpAndDerivsModule
         use GaussQuadrature
         use MappedGeometryClass
         use StorageClass                    , only: FaceStorage_t
         use PhysicsStorage
         use NodalStorageClass
         use InterpolationMatrices           , only: Tset, TsetM
         IMPLICIT NONE 
   
         private
         public   Face
   !
   !     ************************************************************************************
   !
   !           Face derived type. The face connects two elements,
   !           =================  specified in "elementIDs". 
   !
   !        The face handles two problems:
   !
   !           1) Rotation: faces are not oriented in the same direction.
   !              --------  The direction adopted is always that of the left
   !                        element.
   !           2) Different polynomial order: elements at each sidecan have different
   !              --------------------------  polynomial orders.
   !
   !        We have defined the following quantities:
   !
   !           -> NelLeft/Right: The face original order in the element, without
   !              ~~~~~~~~~~~~~  rotation or projections.
   !
   !           -> NfLeft/Right: The face original order rotated to match each other,
   !              ~~~~~~~~~~~~  without projection.
   !                       ** Thus, always NfLeft = NelLeft, and
   !                          NfRight = rotation(NelRight)                           
   !
   !           -> Nf: The face polynomial order. This is the maximum along 
   !              ~~  each side of both directions.
   !
   !        The process from the element storage to face storage is:
   !
   !          (0:Nxyz)           (0:NelLeft)                (0:Nf)
   !        Left element -- Interpolation to Boundary -- Projection to face
   !
   !          (0:Nxyz)           (0:NelRight)             (0:NfRight)       (0:Nf)
   !        Right element -- Interpolation to Boundary --  Rotation -- Projection to face
   !
   !     ************************************************************************************
   !

         type Face
            type(Face),pointer              :: Mortar(:)                !Array of face pointers to slave mortar faces, (only for big master mortar)
            integer                         :: IsMortar                 !0 = conforming, 1 = big master mortar, 2 = small slave 
            logical                         :: flat
            integer                         :: ID                       ! face ID
            integer                         :: FaceType                 ! Type of face: 1 = HMESH_INTERIOR, 2 = HMESH_BOUNDARY, 3 = 2 = HMESH_MPI
            integer                         :: zone                     ! In the case of HMESH_BOUNDARY, which zone it belongs to
            integer                         :: rotation                 ! Relative orientation between faces
            integer                         :: NelLeft(2)               ! Left element face polynomial order
            integer                         :: NelRight(2)              ! Right element face polynomial order
            integer                         :: NfLeft(2)                ! Left face polynomial order
            integer                         :: NfRight(2)               ! Right face polynomial order
            integer                         :: Nf(2)                    ! Face polynomial order
            integer                         :: nodeIDs(4)               ! Face nodes ID
            integer                         :: elementIDs(2)            ! Convention is: 1 = Left, 2 = Right
            integer                         :: elementSide(2)     
            integer                         :: projectionType(2)  
            CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
            type(MappedGeometryFace)        :: geom
            type(FaceStorage_t)             :: storage(2)
            contains
               procedure   :: Construct                     => ConstructFace
               procedure   :: Destruct                      => DestructFace
               procedure   :: Print                         => PrintFace
               procedure   :: LinkWithElements              => Face_LinkWithElements
               procedure   :: AdaptSolutionToFace           => Face_AdaptSolutionToFace
               procedure   :: AdaptSolutionToMortarFace     => Face_AdaptSolutionToMortarFace
               procedure   :: AdaptGradientsToFace          => Face_AdaptGradientsToFace
               procedure   :: AdaptGradientsToMortarFace    => Face_AdaptGradientsToMortarFace
               procedure   :: AdaptAviscFluxToFace          => Face_AdaptAviscFluxToFace
               procedure   :: AdaptAviscFluxToMortarFace    => Face_AdaptAviscFluxToMortarFace

               procedure   :: ProjectFluxToElements                => Face_ProjectFluxToElements
               procedure   :: ProjectMortarFluxToElements          => Face_ProjectMortarFluxToElements
               procedure   :: ProjectGradientFluxToElements        => Face_ProjectGradientFluxToElements
               procedure   :: ProjectMortarGradientFluxToElements  => Face_ProjectMortarGradientFluxToElements
#if defined(NAVIERSTOKES)
               procedure   :: ProjectFluxJacobianToElements       => Face_ProjectFluxJacobianToElements
               procedure   :: ProjectMortarFluxJacobianToElements => Face_ProjectMortarFluxJacobianToElements
               procedure   :: ProjectGradJacobianToElements       => Face_ProjectGradJacobianToElements
               procedure   :: ProjectMortarGradJacobianToElements => Face_ProjectMortarGradJacobianToElements
               procedure   :: ProjectBCJacobianToElements         => Face_ProjectBCJacobianToElements
#endif
               procedure   :: copy           => Face_Assign
               generic     :: assignment(=)  => copy
         end type Face
   !
   !     ========
         CONTAINS
   !     ========
   !
         SUBROUTINE ConstructFace( self, ID, nodeIDs, elementID, side )
   !
   !        *******************************************************
   !
   !           This is not the full face construction, but a 
   !         variable initialization. Just nodes are introduced and
   !         the left element.
   !
   !        *******************************************************
   !
            IMPLICIT NONE 
            class(Face) :: self
            integer     :: ID, nodeIDs(4), elementID, side
   !
   !        --------------------------
   !        Set nodes and left element 
   !        --------------------------
   !
            self % ID             = ID
            self % nodeIDS        = nodeIDs
            self % elementIDs(1)  = elementID
            self % elementIDs(2)  = -1
            self % elementSide(1) = side
   !
   !        ------------
   !        Set defaults
   !        ------------
   !
            self % FaceType       = HMESH_UNDEFINED
            self % elementIDs(2)  = HMESH_NONE
            self % elementSide(2) = HMESH_NONE
            self % boundaryName   = emptyBCName
            self % rotation       = 0
            self % zone           = 0
            self % IsMortar       = 0
   
         end SUBROUTINE ConstructFace
   !
   !////////////////////////////////////////////////////////////////////////
   !
         elemental SUBROUTINE DestructFace( self )
            IMPLICIT NONE 
            class(Face), intent(inout) :: self
            
            self % ID = -1
            self % FaceType = HMESH_NONE
            self % rotation = 0
            self % NelLeft = -1
            self % NelRight = -1       
            self % NfLeft = -1       
            self % NfRight = -1       
            self % Nf = -1         
            self % nodeIDs = -1             
            self % elementIDs = -1
            self % elementSide = -1
            self % projectionType = -1
            self % boundaryName = ""
            
            call self % geom % Destruct
            call self % storage % Destruct
   
         end SUBROUTINE DestructFace
   !
   !////////////////////////////////////////////////////////////////////////
   !
         SUBROUTINE PrintFace( self ) 
         IMPLICIT NONE
         class(Face) :: self
         PRINT *, "Face type = "   , self % FaceType
         PRINT *, "Element IDs: "  , self % elementIDs
         PRINT *, "Element Sides: ", self % elementSide
         IF ( self % FaceType == HMESH_INTERIOR )     THEN
            PRINT *, "Neighbor rotation: ", self  %  rotation
         ELSE
            PRINT *, "Boundary name = ", self % boundaryName
         end IF
         PRINT *, "-----------------------------------"
         end SUBROUTINE PrintFace
   !
   !////////////////////////////////////////////////////////////////////////
   !
      SUBROUTINE Face_LinkWithElements( self, NelLeft, NelRight, nodeType, offset)
         IMPLICIT NONE
         class(Face)        ,     intent(INOUT) :: self        ! Current face
         integer,                 intent(in)    :: NelLeft(2)  ! Left element face polynomial order
         integer,                 intent(in)    :: NelRight(2) ! Right element face polynomial order
         integer,                 intent(in)    :: nodeType    ! Either Gauss or Gauss-Lobatto
         real(kind=RP), optional, intent(in)    :: offset

         real(kind=RP) :: test1(0:3)
         real(kind=RP) :: test2(0:3)
       real(kind=RP) :: test3(0:3)
       real(kind=RP) :: test4(1:3,0:3, 0:3)

       real(kind=RP) :: small(1:3, 0:3, 0:3,4)
       real(kind=RP) :: MInt(0 : 3, 0 : 3,1:2)
       real(kind=RP) :: tmp(1:3, 0:3,0:3)
       real(kind=RP)     :: Flux_tmp(1:3, 0:3, 0:3)
         real(kind=RP)     :: Flux_master(1:3, 0: 3, 0:3,4)
         !real(kind=RP)     :: small(1:3,0: self % NfLeft(1), 0:self % NfLeft(1),4)
         real(kind=RP)     :: Mout(0: 3, 0:3, 2)
         integer :: p,q,l,lm 
         !real(kind=RP) :: test4(self % NfLeft(1), self % Nf(1))
         !real(kind=RP) :: test5(self % NfLeft(1), self % Nf(1))
         integer :: i ,j
         real(kind=RP) ::error
         real(kind=RP) ::w(4)
#if (!defined(NAVIERSTOKES)) && (!defined(INCNS))
         logical  :: computeGradients = .true.
#endif
   !     i
   !     -------------------------------------------------------------
   !     First, get face elements polynomial orders (without rotation)
   !     -------------------------------------------------------------
   !
         self % NelLeft = NelLeft
         self % NelRight = NelRight
   !
   !     ---------------------------------------------------------
   !     Second, get face polynomial orders (considering rotation)
   !     ---------------------------------------------------------
   !
         self % NfLeft = self % NelLeft     ! Left elements are always oriented.
         
         SELECT CASE ( self % rotation )
         CASE ( 0, 2, 5, 7 ) ! Local x and y axis are parallel or antiparallel
            self % NfRight = self % NelRight
         CASE ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
            self % NfRight(2)  = self % NelRight(1)
            self % NfRight(1)  = self % NelRight(2)
         end SELECT
   !
   !     ------------------------------------------------------------------
   !     Third, the face polynomial order (the maximum for both directions)
   !     ------------------------------------------------------------------
   !
         self % Nf(1) = max(self % NfLeft(1), self % NfRight(1))
         self % Nf(2) = max(self % NfLeft(2), self % NfRight(2))
   !
   !     ------------------------------
   !     Construct needed nodal storage
   !     ------------------------------
   !
         if ( .not. NodalStorage(self % Nf(1)) % Constructed ) then
            call NodalStorage(self % Nf(1)) % Construct(nodeType, self % Nf(1))
         end if
   
         if ( .not. NodalStorage(self % Nf(2)) % Constructed ) then
            call NodalStorage(self % Nf(2)) % Construct(nodeType, self % Nf(2))
         end if
   !
   !     -----------------------------------------------------------------------
   !     Construction of the projection matrices (simple Lagrange interpolation)
   !     -----------------------------------------------------------------------
   !
         if (.not.present(offset)) then 
            call Tset(self % NfLeft(1), self % Nf(1)) % construct(self % NfLeft(1), self % Nf(1))
            call Tset(self % Nf(1), self % NfLeft(1)) % construct(self % Nf(1), self % NfLeft(1))
            
            call Tset(self % NfLeft(2), self % Nf(2)) % construct(self % NfLeft(2), self % Nf(2))
            call Tset(self % Nf(2), self % NfLeft(2)) % construct(self % Nf(2), self % NfLeft(2))

         else 
            call TsetM(self % NfLeft(1), self % Nf(1), 2, 1) % construct(self % NfLeft(1), self % Nf(1), 0.5_RP, 1)
            call TsetM(self % Nf(1), self % NfLeft(1), 2, 2) % construct(self % Nf(1), self % NfLeft(1), 0.5_RP, 2)
   
            call TsetM(self % NfLeft(1), self % Nf(1), 1, 1) % construct(self % NfLeft(1), self % Nf(1), -0.5_RP, 1)
            call TsetM(self % Nf(1), self % NfLeft(1), 1, 2) % construct(self % Nf(1), self % NfLeft(1), -0.5_RP, 2)
            
            call TsetM(self % NfLeft(2), self % Nf(2), 2, 1) % construct(self % NfLeft(2), self % Nf(2), 0.5_RP, 1)
            call TsetM(self % Nf(2), self % NfLeft(2), 2, 2) % construct(self % Nf(2), self % NfLeft(2), 0.5_RP, 2)
   
            call TsetM(self % NfLeft(2), self % Nf(2), 1, 1) % construct(self % NfLeft(2), self % Nf(2), -0.5_RP, 1)
            call TsetM(self % Nf(2), self % NfLeft(2), 1, 2) % construct(self % Nf(2), self % NfLeft(2), -0.5_RP, 2)

           ! MInt(:,:,1)=transpose(TsetM(self % NfLeft(1), self % Nf(1), 1, 1) % T)
           ! MInt(:,:,2)=transpose(TsetM(self % NfLeft(1), self % Nf(1), 2, 1) % T)
           ! CALL RANDOM_NUMBER(test4)
           ! write(*,*) 'test1 befor', test4
           ! do j=1,2
              ! tmp=0.0_RP
            !   do q=0, self%Nf(1) ; do p=0, self%Nf(1) ; do l=0, self%Nf(1)
                !     tmp(:,p,q)=tmp(:,p,q) +MInt(l,q,j)*test4(:,p,l)
                !  end do ; end do ; end do

              ! do i=1,2
              ! small(:,:,:,i+2*(j-1))=0.0_RP
               !   do q=0, self%Nf(1) ; do p=0, self%Nf(1) ; do l=0, self%Nf(1)
                !     small(:,p,q,i+2*(j-1))=small(:,p,q,i+2*(j-1)) +MInt(l,p,i)*tmp(:,l,q)
               !   end do ; end do ; end do 
              ! end do 
          !  end do 


          !  Mout(:,:,1)=transpose(TsetM(self % NfLeft(1), self % Nf(1), 1, 2) % T)
          !  Mout(:,:,2)=transpose(TsetM(self % NfLeft(1), self % Nf(1), 2, 2) % T)
          !  Flux_master=0.0_RP
           ! do lm=1,4
            !   do j=1,2
            !      Flux_tmp(:,:,:)=0.0_RP
             ! !    do i=1,2
               !      do q=0,self % NfLeft(1) ; do p=0,self % NfLeft(1) ; do l=0,self % NfLeft(1)
                !           Flux_tmp(:,p,q)=Flux_tmp(:,p,q) + Mout(l,p,i)*small(:,l,q,i+2*(j-1))
               !      end do ; end do ; end do 
             !     end do 
               
             !     do p=0,self % NfLeft(1) ; do q=0,self % NfLeft(1) ; do l=0,self % NfLeft(1)
             !           Flux_master(:,p,q,lm)=Flux_master(:,p,q,lm) + Mout(l,q,j)*Flux_tmp(:,p,l) 
              !    end do ; end do ; end do 
              ! end do 
            !end do 

           ! MInt(:,:,1)=(TsetM(self % NfLeft(1), self % Nf(1), 1, 1) % T)
           ! MInt(:,:,2)=(TsetM(self % NfLeft(1), self % Nf(1), 2, 1) % T)
           ! Mout(:,:,1)=(TsetM(self % NfLeft(1), self % Nf(1), 1, 2) % T)
            !Mout(:,:,2)=(TsetM(self % NfLeft(1), self % Nf(1), 2, 2) % T)
            
          !  write(*,*) 'matmul'
          !!  write(*,*)  MATMUL(Mout(:,:,1),MInt(:,:,1))+ MATMUL(Mout(:,:,2),MInt(:,:,2))
            !test1=MATMUL(TsetM(self % NfLeft(1), self % Nf(1), 2, 2) % T, TsetM(self % NfLeft(1), self % Nf(1), 2, 1) % T)  
            !test2=MATMUL(TsetM(self % NfLeft(1), self % Nf(1), 1, 2) % T, TsetM(self % NfLeft(1), self % Nf(1), 1, 1) % T)  
            !test3=test1+test2 
          !  write(*,*) 'test afterp', Flux_master(:,:,:,1)+Flux_master(:,:,:,2)+Flux_master(:,:,:,3)+Flux_master(:,:,:,4)
          !  test2 = 1.33_RP
          !  test1 = MATMUL(MInt(:,:,1),test2) !interpolate to mortar1
          !  error = SUM(ABS(test1-1.33_RP))/REAL(4)
         !   write(*,*)'Error of interpolate constant to mortar 1:',error
         !   test2 = MATMUL(MInt(:,:,2),test2) !interpolate to mortar2
         !   error = SUM(ABS(test2-1.33_RP))/REAL(4)
         !   write(*,*)'Error of interpolate constant to mortar 2:',error
           ! test2 = MATMUL(Mout(:,:,1),test1)+MATMUL(Mout(:,:,2),test2)
            !error = SUM(ABS(0.5_RP*test2-1.33_RP))/REAL(4)
           ! write(*,*)'Error of project constant back to big side:',error
          !  IF(error.GT. 100.*epsilon(1._RP)) THEN
             !  write(*,*) 'aie'
           ! else 
            !   write(*,*) 'good'
            ! END IF
            ! w=NodalStorage(3)%wb 


             !  test2 = 1.0d0
             !  test2(0) = -10.33d0
             !  error  = SUM(test2*w)  !save mean value of big side
             !  test1 = MATMUL(MInt(:,:,1),test2) !interpolate to mortar1
             !  test2 = MATMUL(MInt(:,:,2),test2) !interpolate to mortar2
             !  test2 = 0.5*(MATMUL(Mout(:,:,1),test1)+MATMUL(Mout(:,:,2),test2)) !project  back
             !  error  = error - SUM(test2*w)  !difference to initial mean value
             !  write(*,*)'Error of mean value of polynomial, projected back to big side:',error
             !  IF(error.GT. 100.*epsilon(1._RP)) THEN
              !    write(*,*) 'aie'
             !  else 
               !   write(*,*) 'Mortar operators build successfully.'
              !  END IF

         end if 
         
         call Tset(self % NfRight(1), self % Nf(1)) % construct(self % NfRight(1), self % Nf(1))
         call Tset(self % Nf(1), self % NfRight(1)) % construct(self % Nf(1), self % NfRight(1))
         
         call Tset(self % NfRight(2), self % Nf(2)) % construct(self % NfRight(2), self % Nf(2))
         call Tset(self % Nf(2), self % NfRight(2)) % construct(self % Nf(2), self % NfRight(2))


      
   !
   !     -----------------------  0- no projection
   !     Set the projection type: 1- x needs projection
   !     -----------------------  2- y needs projection
   !                              3- both x and y need projection
         if ( (self % NfLeft(1) .eq. self % Nf(1)) .and. (self % NfLeft(2) .eq. self % Nf(2)) ) then
            self % projectionType(1) = 0
         elseif ( (self % NfLeft(1) .ne. self % Nf(1)) .and. (self % NfLeft(2) .eq. self % Nf(2)) ) then
            self % projectionType(1) = 1
         elseif ( (self % NfLeft(1) .eq. self % Nf(1)) .and. (self % NfLeft(2) .ne. self % Nf(2)) ) then
            self % projectionType(1) = 2
         elseif ( (self % NfLeft(1) .ne. self % Nf(1)) .and. (self % NfLeft(2) .ne. self % Nf(2)) ) then
            self % projectionType(1) = 3
         end if
   
         if ( (self % NfRight(1) .eq. self % Nf(1)) .and. (self % NfRight(2) .eq. self % Nf(2)) ) then
            self % projectionType(2) = 0
         elseif ( (self % NfRight(1) .ne. self % Nf(1)) .and. (self % NfRight(2) .eq. self % Nf(2)) ) then
            self % projectionType(2) = 1
         elseif ( (self % NfRight(1) .eq. self % Nf(1)) .and. (self % NfRight(2) .ne. self % Nf(2)) ) then
            self % projectionType(2) = 2
         elseif ( (self % NfRight(1) .ne. self % Nf(1)) .and. (self % NfRight(2) .ne. self % Nf(2)) ) then
            self % projectionType(2) = 3
         end if
   
      end SUBROUTINE Face_LinkWithElements
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !
      subroutine Face_AdaptSolutionToFace(self, nEqn, Nelx, Nely, Qe, side, QdotE, computeQdot, fma, fmb, fmc, fmd)
         use MappedGeometryClass
         implicit none
         class(Face),   intent(inout)              :: self
         integer,       intent(in)                 :: nEqn
         integer,       intent(in)                 :: Nelx, Nely
         real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
         integer,       intent(in)                 :: side
         real(kind=RP), intent(in), optional       :: QdotE(1:nEqn, 0:Nelx, 0:Nely)
         logical,       intent(in), optional       :: computeQdot
         type(Face), optional, intent(inout)      ::fma
         type(Face), optional, intent(inout)      ::fmb
         type(Face), optional, intent(inout)      :: fmc
         type(Face), optional, intent(inout)      :: fmd 
         
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer       :: i, j, k, l, m, ii, jj, a, b, lm, inb, jnb , p, q
         real(kind=RP) :: Qe_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         ! real(kind=RP) :: QdotE_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         logical :: prolongQdot

         ! prolongQdot = present(QdotE)
         if (present(computeQdot)) then
             prolongQdot = computeQdot
         else
             prolongQdot = .FALSE.
         end if
   
         ! if (prolongQdot) then
         !     print *, "side: ", side
         !     ! print *, "projectionType 1: ",  self % projectionType(1)
         !     ! print *, "projectionType 2: ",  self % projectionType(2)
         !     print *, "projectionType side: ",  self % projectionType(side)
         ! end if
         !write(*,*) "proj sol toface", self % ID 
      if (self % isMortar==0 .OR. self % isMortar==2 ) then
         select case (side)
         case(1)
            if (self % IsMortar==2) then 
               write(*,*) 'slave mortar on wrong side'
            end if 
               associate(Qf => self % storage(1) % Q)
               select case ( self % projectionType(1) )
               case (0)
                  Qf = Qe
                  if (prolongQdot) self % storage(1) % Qdot = QdotE
   
               case (1)
                  Qf = 0.0_RP
                  do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
                     Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Qe(:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  Qf = 0.0_RP
                  do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
                     Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Qe(:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  Qf = 0.0_RP
                  do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
                     do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                        Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                                * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                                * Qe(:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
               end associate
         case(2)

            associate( Qf => self % storage(2) % Q )
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               Qe_rot(:,i,j) = Qe(:,ii,jj) 
            end do                        ; end do
   
            select case ( self % projectionType(2) )
            case (0)
               Qf = Qe_rot
            case (1)
               Qf = 0.0_RP
               do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Qe_rot(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               Qf = 0.0_RP
               do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Qe_rot(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               Qf = 0.0_RP
               do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
                  do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                     Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                               * Qe_rot(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate
         end select

      end if 

      end subroutine Face_AdaptSolutionToFace

      subroutine Face_AdaptSolutionToMortarFace(self, nEqn,Nelx,Nely, Qe, side, QdotE, computeQdot, fma, fmb, fmc, fmd)
         use MappedGeometryClass
         implicit none
         class(Face),   intent(inout)              :: self
         integer,       intent(in)                 :: nEqn
         integer,       intent(in)                 :: Nelx, Nely
         real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
         integer,       intent(in)                 :: side
         real(kind=RP), intent(in), optional       :: QdotE(1:nEqn, 0:Nelx, 0:Nely)
         logical,       intent(in), optional       :: computeQdot
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
         
         integer       :: i, j, l, lm , p, q


         logical :: prolongQdot
         real(kind=RP) :: small(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2),4)
         real(kind=RP) :: MInt(0:fma%Nf(1),0:fma%Nf(1),1:2)
         real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1),0:fma%Nf(1))


         if (self % isMortar==1)  then 

            if(side .ne. 1) then 
               write(*,*)'problem big mortar on wrong side'
            end if 
            if ((fma % IsMortar.ne.2) .or. (fmb % IsMortar.ne.2) .or. (fmc % IsMortar.ne.2) .or. (fmd % IsMortar.ne.2)) then 
               write(*,*) 'problem mortar slave projection sol to face'
            end if 

            MInt(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 1) % T)
            MInt(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 1) % T)

            do j=1,2
               tmp=0.0_RP
               do q=0, fma%Nf(1) ; do p=0, fma%Nf(1) ; do l=0, fma%Nf(1)
                     tmp(:,p,q)=tmp(:,p,q) +MInt(l,q,j)*Qe(:,p,l)
                  end do ; end do ; end do

               do i=1,2
               small(:,:,:,i+2*(j-1))=0.0_RP
                  do q=0, fma%Nf(1) ; do p=0, fma%Nf(1) ; do l=0, fma%Nf(1)
                     small(:,p,q,i+2*(j-1))=small(:,p,q,i+2*(j-1)) +MInt(l,p,i)*tmp(:,l,q)
                  end do ; end do ; end do 
               end do 
            end do 

            associate(Qf => fma % storage(1) % Q)
               Qf=small(:,:,:,1)
            end associate 
            associate(Qf => fmb % storage(1) % Q)
               Qf=small(:,:,:,2)
            end associate 
            associate(Qf => fmc % storage(1) % Q)
               Qf=small(:,:,:,3)
            end associate 
            associate(Qf => fmd % storage(1) % Q)
               Qf=small(:,:,:,4)
            end associate 

           ! associate(Qf => fma % storage(2) % Q)
           !    Qf=small(:,:,:,1)
           ! end associate 
           ! associate(Qf => fmb % storage(2) % Q)
               !Qf=small(:,:,:,2)
           !! end associate 
           ! associate(Qf => fmc % storage(2) % Q)
           !    Qf=small(:,:,:,4)
           ! end associate 
            !associate(Qf => fmd % storage(2) % Q)
            !   Qf=small(:,:,:,3)
           ! end associate 

         end if

      end subroutine Face_AdaptSolutionToMortarFace
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !
      subroutine Face_AdaptGradientsToFace(self, nEqn, Nelx, Nely, Uxe, Uye, Uze, side, fma, fmb, fmc ,fmd)
         use MappedGeometryClass
         implicit none
         class(Face),   intent(inout)  :: self
         integer,       intent(in)     :: nEqn 
         integer,       intent(in)     :: Nelx, Nely
         real(kind=RP), intent(in)     :: Uxe(nEqn, 0:Nelx, 0:Nely)
         real(kind=RP), intent(in)     :: Uye(nEqn, 0:Nelx, 0:Nely)
         real(kind=RP), intent(in)     :: Uze(nEqn, 0:Nelx, 0:Nely)
         integer,       intent(in)     :: side
         type(Face), optional, intent(inout)      ::fma
         type(Face), optional, intent(inout)      ::fmb
         type(Face), optional, intent(inout)      :: fmc
         type(Face), optional, intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer       :: i, j, k, l, m, ii, jj, a, b, lm
         real(kind=RP) :: Uxe_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         real(kind=RP) :: Uye_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         real(kind=RP) :: Uze_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      if (self % isMortar==0 .OR. self % isMortar==2 ) then
         select case (side)
         case(1)
            associate(Uxf => self % storage(1) % U_x, &
                      Uyf => self % storage(1) % U_y, &
                      Uzf => self % storage(1) % U_z   )
            select case ( self % projectionType(1) )
            case (0)
               Uxf = Uxe
               Uyf = Uye
               Uzf = Uze
            case (1)
               Uxf = 0.0_RP
               Uyf = 0.0_RP
               Uzf = 0.0_RP
               do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uxe(:,l,j)
                  Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uye(:,l,j)
                  Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uze(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               Uxf = 0.0_RP
               Uyf = 0.0_RP
               Uzf = 0.0_RP
               do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uxe(:,i,l)
                  Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uye(:,i,l)
                  Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uze(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               Uxf = 0.0_RP
               Uyf = 0.0_RP
               Uzf = 0.0_RP
               do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
                  do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                     Uxf(:,i,j) = Uxf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                               * Uxe(:,m,l)
                     Uyf(:,i,j) = Uyf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                               * Uye(:,m,l)
                     Uzf(:,i,j) = Uzf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                               * Uze(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate
         case(2)
            associate(Uxf => self % storage(2) % U_x, &
                      Uyf => self % storage(2) % U_y, &
                      Uzf => self % storage(2) % U_z   )
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               Uxe_rot(:,i,j) = Uxe(:,ii,jj) 
               Uye_rot(:,i,j) = Uye(:,ii,jj) 
               Uze_rot(:,i,j) = Uze(:,ii,jj) 
            end do                        ; end do
   
            select case ( self % projectionType(2) )
            case (0)
               Uxf = Uxe_rot
               Uyf = Uye_rot
               Uzf = Uze_rot
            case (1)
               Uxf = 0.0_RP
               Uyf = 0.0_RP
               Uzf = 0.0_RP
               do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uxe_rot(:,l,j)
                  Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uye_rot(:,l,j)
                  Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uze_rot(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               Uxf = 0.0_RP
               Uyf = 0.0_RP
               Uzf = 0.0_RP
               do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uxe_rot(:,i,l)
                  Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uye_rot(:,i,l)
                  Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uze_rot(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               Uxf = 0.0_RP
               Uyf = 0.0_RP
               Uzf = 0.0_RP
               do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
                  do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                     Uxf(:,i,j) = Uxf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                               * Uxe_rot(:,m,l)
                     Uyf(:,i,j) = Uyf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                               * Uye_rot(:,m,l)
                     Uzf(:,i,j) = Uzf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                               * Uze_rot(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate
         end select
     end if 
      
      end subroutine Face_AdaptGradientsToFace

      subroutine Face_AdaptGradientsToMortarFace(self, nEqn, Nelx, Nely, Uxe, Uye, Uze, side, fma, fmb, fmc ,fmd)
         use MappedGeometryClass
         implicit none
         class(Face),   intent(inout)  :: self
         integer,       intent(in)     :: nEqn 
         integer,       intent(in)     :: Nelx, Nely
         real(kind=RP), intent(in)     :: Uxe(nEqn, 0:Nelx, 0:Nely)
         real(kind=RP), intent(in)     :: Uye(nEqn, 0:Nelx, 0:Nely)
         real(kind=RP), intent(in)     :: Uze(nEqn, 0:Nelx, 0:Nely)
         integer,       intent(in)     :: side
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer       :: i, j, l, lm , p, q
         real(kind=RP) :: smallx(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2),4)
         real(kind=RP) :: smally(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2),4)
         real(kind=RP) :: smallz(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2),4)
         real(kind=RP) :: MInt(0:fma%Nf(1),0:fma%Nf(1),1:2)
         real(kind=RP) :: tmp_x(1:nEqn, 0:fma%Nf(1),0:fma%Nf(1))
         real(kind=RP) :: tmp_y(1:nEqn, 0:fma%Nf(1),0:fma%Nf(1))
         real(kind=RP) :: tmp_z(1:nEqn, 0:fma%Nf(1),0:fma%Nf(1))
   
         if (self % isMortar==1)  then 
   
            if(side .ne. 1) then 
               write(*,*)'problem big mortar on wrong side'
            end if 
            if ((fma % IsMortar.ne.2) .or. (fmb % IsMortar.ne.2) .or. (fmc % IsMortar.ne.2) .or. (fmd % IsMortar.ne.2)) then 
               write(*,*) 'problem mortar slave projection sol to face'
            end if 
   
            MInt(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 1) % T)
            MInt(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 1) % T)
            
            do j=1,2
               tmp_x=0.0_RP
               tmp_y=0.0_RP
               tmp_z=0.0_RP
               do q=0, fma%Nf(1) ; do p=0, fma%Nf(1) ; do l=0, fma%Nf(1)
                     tmp_x(:,p,q)=tmp_x(:,p,q) +MInt(l,q,j)*Uxe(:,p,l)
                     tmp_y(:,p,q)=tmp_y(:,p,q) +MInt(l,q,j)*Uye(:,p,l)
                     tmp_z(:,p,q)=tmp_z(:,p,q) +MInt(l,q,j)*Uze(:,p,l)
                  end do ; end do ; end do
   
               do i=1,2
               smallx(:,:,:,i+2*(j-1))=0.0_RP
               smally(:,:,:,i+2*(j-1))=0.0_RP
               smallz(:,:,:,i+2*(j-1))=0.0_RP
                  do q=0, fma%Nf(1) ; do p=0, fma%Nf(1) ; do l=0, fma%Nf(1)
                     smallx(:,p,q,i+2*(j-1))=smallx(:,p,q,i+2*(j-1)) +MInt(l,p,i)*tmp_x(:,l,q)
                     smally(:,p,q,i+2*(j-1))=smally(:,p,q,i+2*(j-1)) +MInt(l,p,i)*tmp_y(:,l,q)
                     smallz(:,p,q,i+2*(j-1))=smallz(:,p,q,i+2*(j-1)) +MInt(l,p,i)*tmp_z(:,l,q)
                  end do ; end do ; end do 
               end do 
            end do 
   
            associate(Uxf =>fma % storage(1) % U_x, &
               Uyf => fma % storage(1) % U_y, &
               Uzf => fma % storage(1) % U_z   )
               Uxf=smallx(:,:,:,1)
               Uyf=smally(:,:,:,1)
               Uzf=smallz(:,:,:,1)
            end associate 
            !associate(Uxf =>fma % storage(2) % U_x, &
            !   Uyf => fma % storage(1) % U_y, &
            !   Uzf => fma % storage(1) % U_z   )
            !   Uxf=smallx(:,:,:,1)
            !   Uyf=smally(:,:,:,1)
           !    Uzf=smallz(:,:,:,1)
           ! end associate 
   
            associate(Uxf =>fmb % storage(1) % U_x, &
               Uyf => fmb % storage(1) % U_y, &
               Uzf => fmb % storage(1) % U_z   )
               Uxf=smallx(:,:,:,2)
               Uyf=smally(:,:,:,2)
               Uzf=smallz(:,:,:,2)
            end associate 
           ! associate(Uxf =>fmb % storage(2) % U_x, &
            !   Uyf => fmb % storage(1) % U_y, &
            !   Uzf => fmb % storage(1) % U_z   )
            !  Uxf=smallx(:,:,:,2)
            !   Uyf=smally(:,:,:,2)
            !   Uzf=smallz(:,:,:,2)
           ! end associate 
   
   
            associate(Uxf =>fmc % storage(1) % U_x, &
               Uyf => fmc % storage(1) % U_y, &
               Uzf => fmc % storage(1) % U_z   )
               Uxf=smallx(:,:,:,3)
               Uyf=smally(:,:,:,3)
               Uzf=smallz(:,:,:,3)
            end associate
           ! associate(Uxf =>fmc % storage(2) % U_x, &
            !   Uyf => fmc % storage(1) % U_y, &
            !   Uzf => fmc % storage(1) % U_z   )
            !   Uxf=smallx(:,:,:,4)
             !  Uyf=smally(:,:,:,4)
             !  Uzf=smallz(:,:,:,4)
           ! end associate
   
            associate(Uxf =>fmd % storage(1) % U_x, &
               Uyf => fmd % storage(1) % U_y, &
               Uzf => fmd % storage(1) % U_z   )
               Uxf=smallx(:,:,:,4)
               Uyf=smally(:,:,:,4)
               Uzf=smallz(:,:,:,4)
            end associate
           ! associate(Uxf =>fmd % storage(2) % U_x, &
             !  Uyf => fmd % storage(1) % U_y, &
             !  Uzf => fmd % storage(1) % U_z   )
             !  Uxf=smallx(:,:,:,3)
             !  Uyf=smally(:,:,:,3)
            !   Uzf=smallz(:,:,:,3)
           ! end associate
         end if
   
   
      end subroutine Face_AdaptGradientsToMortarFace
   
   
      subroutine Face_AdaptAviscFluxToFace(self, nEqn, Nelx, Nely, AVn_e, side, fma, fmb, fmc, fmd)
         use MappedGeometryClass
         implicit none
         class(Face),   intent(inout)  :: self
         integer,       intent(in)     :: nEqn
         integer,       intent(in)     :: Nelx, Nely
         real(kind=RP), intent(in)     :: AVn_e(1:nEqn, 0:Nelx, 0:Nely)
         integer,       intent(in)     :: side
         type(Face), optional, intent(inout)      ::fma
         type(Face), optional, intent(inout)      ::fmb
         type(Face), optional, intent(inout)      :: fmc
         type(Face), optional, intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer       :: i, j, k, l, m, ii, jj, lm, a, b 
         real(kind=RP) :: AVn_e_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      if (self % isMortar==0 .OR. self % isMortar==2 ) then
         select case (side)
         case(1)
            associate(AVf => self % storage(1) % AviscFlux)
            select case ( self % projectionType(1) )
            case (0)
               AVf = AVn_e
            case (1)
               AVf = 0.0_RP
               do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * AVn_e(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               AVf = 0.0_RP
               do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * AVn_e(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               AVf = 0.0_RP
               do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
                  do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                     AVf(:,i,j) = AVf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                               * AVn_e(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            end associate
         case(2)
            associate( AVf => self % storage(2) % AviscFlux )
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               AVn_e_rot(:,i,j) = AVn_e(:,ii,jj) 
            end do                        ; end do
   
            select case ( self % projectionType(2) )
            case (0)
               AVf = AVn_e_rot
            case (1)
               AVf = 0.0_RP
               do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * AVn_e_rot(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               AVf = 0.0_RP
               do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * AVn_e_rot(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               AVf = 0.0_RP
               do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
                  do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                     AVf(:,i,j) = AVf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                               * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                               * AVn_e_rot(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
   
            AVf = -AVf
            end associate
         end select
      end if 
   
      end subroutine Face_AdaptAviscFluxToFace

      subroutine Face_AdaptAviscFluxToMortarFace(self, nEqn, Nelx, Nely, AVn_e, side, fma, fmb, fmc, fmd)
         use MappedGeometryClass
         implicit none
         class(Face),   intent(inout)  :: self
         integer,       intent(in)     :: nEqn
         integer,       intent(in)     :: Nelx, Nely
         real(kind=RP), intent(in)     :: AVn_e(1:nEqn, 0:Nelx, 0:Nely)
         integer,       intent(in)     :: side
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer       :: i, j, l, m, lm, p, q
         real(kind=RP) :: small(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2),4)
         real(kind=RP) :: MInt(0:fma%Nf(1),0:fma%Nf(1),1:2)
         real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1),0:fma%Nf(1))
   
         if (self % isMortar==1)  then 
   
            if(side .ne. 1) then 
               write(*,*)'problem big mortar on wrong side'
            end if 
            if ((fma % IsMortar.ne.2) .or. (fmb % IsMortar.ne.2) .or. (fmc % IsMortar.ne.2) .or. (fmd % IsMortar.ne.2)) then 
               write(*,*) 'problem mortar slave projection sol to face'
            end if 
   
            MInt(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 1) % T)
            MInt(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 1) % T)
   
            do j=1,2
               tmp=0.0_RP
               do q=0, fma%Nf(1) ; do p=0, fma%Nf(1) ; do l=0, fma%Nf(1)
                     tmp(:,p,q)=tmp(:,p,q) +MInt(l,q,j)* AVn_e(:,p,l)
                  end do ; end do ; end do
   
               do i=1,2
               small(:,:,:,i+2*(j-1))=0.0_RP
                  do q=0, fma%Nf(1) ; do p=0, fma%Nf(1) ; do l=0, fma%Nf(1)
                     small(:,p,q,i+2*(j-1))=small(:,p,q,i+2*(j-1)) +MInt(l,p,i)*tmp(:,l,q)
                  end do ; end do ; end do 
               end do 
            end do 
   
            associate(AVf => fma % storage(1) % AviscFlux)
               AVf=small(:,:,:,1)
            end associate 
            associate(AVf => fmb % storage(1) % AviscFlux)
               AVf=small(:,:,:,2)
            end associate 
            associate(AVf => fmc % storage(1) % AviscFlux)
               AVf=small(:,:,:,3)
            end associate 
            associate(AVf => fmd % storage(1) % AviscFlux)
               AVf=small(:,:,:,4)
            end associate 

           ! associate(AVf => fma % storage(2) % AviscFlux)
            !   AVf=small(:,:,:,1)
           ! end associate 
           ! associate(AVf => fmb % storage(2) % AviscFlux)
            !   AVf=small(:,:,:,2)
            !end associate 
           ! associate(AVf => fmc % storage(2) % AviscFlux)
            !   AVf=small(:,:,:,4)
           ! end associate 
           ! associate(AVf => fmd % storage(2) % AviscFlux)
           !    AVf=small(:,:,:,3)
           ! end associate 
   
         end if
   
      end subroutine Face_AdaptAviscFluxToMortarFace
   
      subroutine Face_ProjectFluxToElements(self, nEqn, flux, whichElements, flux_M1, flux_M2, flux_M3, flux_M4, fma, fmb, &
         fmc, fmd)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face)       :: self
         integer,       intent(in)  :: nEqn
         real(kind=RP), optional, intent(in)  :: flux(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         integer,       intent(in)  :: whichElements(2)
   
         real(kind=RP), optional, intent(in)  :: flux_M1(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         real(kind=RP), optional, intent(in)  :: flux_M2(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         real(kind=RP), optional, intent(in)  :: flux_M3(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         real(kind=RP), optional, intent(in)  :: flux_M4(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         type(Face), optional, intent(inout)      ::fma
         type(Face), optional, intent(inout)      ::fmb
         type(Face), optional, intent(inout)      :: fmc
         type(Face), optional, intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer           :: i, j, ii, jj, l, m, side, q, p, lm
         real(kind=RP)     :: fStarAux(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         real(kind=RP),allocatable     :: Flux_tmp(:, :, :)
         real(kind=RP),allocatable     :: Flux_master(:, :, :, :)
         real(kind=RP),allocatable     :: small(:, :, :, :)
         real(kind=RP),allocatable     :: Mout(:, :, :)

         do side = 1, 2
            select case ( whichElements(side) )
            case (1)    ! Prolong to left element
               if (self % IsMortar ==0) then
                  associate(fStar => self % storage(1) % Fstar)
                  select case ( self % projectionType(1) )
                  case (0)
                     fStar(1:nEqn,:,:) = flux
                  case (1)
                     fStar(1:nEqn,:,:) = 0.0
                     do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                        fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * flux(:,l,j)
                     end do                  ; end do                   ; end do
                     
                  case (2)
                     fStar(1:nEqn,:,:) = 0.0
                     do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                        fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * flux(:,i,l)
                     end do                  ; end do                   ; end do
            
                  case (3)
                     fStar(1:nEqn,:,:) = 0.0
                     do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                        do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                           fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) +   Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                  * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                  * flux(:,m,l)
                        end do                 ; end do
                     end do                  ; end do
                  end select
                  end associate
               end if 
 
   
            case (2)    ! Prolong to right element
   !      
   !           *********
   !           1st stage: Projection
   !           *********
   !           
                  if (self % IsMortar==0 .OR. self % IsMortar==2) then
                     select case ( self % projectionType(2) )
                     case (0)
                        fStarAux(1:nEqn,:,:) = flux
                     case (1)
                        fStarAux(1:nEqn,:,:) = 0.0
                        do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                           fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * flux(:,l,j)
                        end do                  ; end do                   ; end do
                        
                     case (2)
                        fStarAux(1:nEqn,:,:) = 0.0
                        do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                           fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * flux(:,i,l)
                        end do                  ; end do                   ; end do
               
                     case (3)
                        fStarAux(1:nEqn,:,:) = 0.0
                        do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                           do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                              fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                                     * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                                     * flux(:,m,l)
                           end do                 ; end do
                        end do                  ; end do
                     end select
         !      
         !           *********
         !           2nd stage: Rotation
         !           *********
         !      
                     associate(fStar => self % storage(2) % Fstar)
                     do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                        call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
                        fStar(1:nEqn,ii,jj) = fStarAux(1:nEqn,i,j) 
                     end do                        ; end do
         !
         !           *********
         !           3rd stage: Inversion
         !           *********
         !
                     
                     fStar = -fStar
                     
                     end associate
                  end if 
               end select
         end do
   
      end subroutine Face_ProjectFluxToElements

      subroutine Face_ProjectMortarFluxToElements(self, nEqn, flux, whichElements, flux_M1, flux_M2, flux_M3, flux_M4, fma, fmb, &
         fmc, fmd)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face)       :: self
         integer,       intent(in)  :: nEqn
         real(kind=RP), optional, intent(in)  :: flux(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         integer,       intent(in)  :: whichElements(2)
   
         real(kind=RP), intent(in)  :: flux_M1(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         real(kind=RP), intent(in)  :: flux_M2(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         real(kind=RP), intent(in)  :: flux_M3(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         real(kind=RP), intent(in)  :: flux_M4(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
   
         integer           :: i, j, ii, jj, l, m, side, q, p, lm
         real(kind=RP)     :: fStarAux(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         real(kind=RP)     :: Flux_tmp(1:nEqn, 0:fma % NfLeft(1), 0:fma % NfLeft(1))
         real(kind=RP)     :: Flux_master(1:nEqn, 0: fma % NfLeft(1), 0:fma % NfLeft(1),4)
         real(kind=RP)     :: small(1:nEqn,0: fma % NfLeft(1), 0:fma % NfLeft(1),4)
         real(kind=RP)     :: Mout(0: fma % NfLeft(1), 0:fma % NfLeft(1), 2)

             !  write(*,*) 'mortar prohection'
         if (self % IsMortar==1) then 

            small(:,:,:,1)=flux_M1
            small(:,:,:,2)=flux_M2
            small(:,:,:,3)=flux_M3
            small(:,:,:,4)=flux_M4

            Mout(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 2) % T)
            Mout(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 2) % T)
            Flux_master=0.0_RP
            do lm=1,4
               do j=1,2
                  Flux_tmp(:,:,:)=0.0_RP
                  do i=1,2
                     do q=0,fma % NfLeft(1) ; do p=0,fma % NfLeft(1) ; do l=0,fma % NfLeft(1)
                           Flux_tmp(:,p,q)=Flux_tmp(:,p,q) + Mout(l,p,i)*small(:,l,q,i+2*(j-1))
                     end do ; end do ; end do 
                  end do 
               
                  do p=0,fma % NfLeft(1) ; do q=0,fma % NfLeft(1) ; do l=0,fma % NfLeft(1)
                        Flux_master(:,p,q,lm)=Flux_master(:,p,q,lm) + Mout(l,q,j)*Flux_tmp(:,p,l) 
                  end do ; end do ; end do 
               end do 
            end do 

            associate(fStar => self % storage(1) % Fstar)

               fStar=0.25_RP*(Flux_master(:,:,:,1)+Flux_master(:,:,:,2)+Flux_master(:,:,:,4)+Flux_master(:,:,:,3))
      
            end associate 

            associate(fStar => self % storage(2) % Fstar)

               fStar=0.25_RP*(Flux_master(:,:,:,1)+Flux_master(:,:,:,2)+Flux_master(:,:,:,4)+Flux_master(:,:,:,3))
               fStar = -fStar
            end associate 
          end if 
      end subroutine Face_ProjectMortarFluxToElements
#if defined(NAVIERSTOKES) 
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !
   !
   !  Face_ProjectBCJacobianToElements:
   !  Routine to project the Jacobian of the boundary condition to the element
   !  --
      subroutine Face_ProjectBCJacobianToElements(self,nEqn,BCjacF)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face), target        :: self
         integer,       intent(in)  :: nEqn
         real(kind=RP), intent(in)  :: BCjacF(nEqn,nEqn,0:self % Nf(1), 0:self % Nf(2))
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer                :: i, j, l, m
         
   !     Always project to left element (this is only for boundaries)
   !     ------------------------------------------------------------
         associate(BCJac => self % storage(1) % BCJac)
         
         select case ( self % projectionType(1) )
            case (0)
               BCJac(1:nEqn,1:nEqn,:,:) = BCJacF
            case (1)
               BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                               + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * BCJacF(:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                               + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * BCJacF(:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                              +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                               * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                               * BCJacF(:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
         end select
         
         end associate
      end subroutine Face_ProjectBCJacobianToElements
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !            
      subroutine Face_ProjectFluxJacobianToElements(self, nEqn, whichElement,whichderiv)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face), target        :: self
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: whichElement
         integer,       intent(in)  :: whichderiv           !<  One can either transfer the derivative with respect to qL (LEFT) or to qR (RIGHT)
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer                :: i, j, ii, jj, l, m, side, lm, a, b
         real(kind=RP), pointer :: fluxDeriv(:,:,:,:)
         real(kind=RP)          :: fStarAux(nEqn,nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
         fluxDeriv(1:,1:,0:,0:) => self % storage(whichderiv) % dFStar_dqF
   
         select case ( whichElement )
            case (LEFT)    ! Prolong to left element
               if (self % IsMortar== 0) then 
                  associate(dFStar_dq => self % storage(1) % dFStar_dqEl)
                  select case ( self % projectionType(1) )
                  case (0)
                     dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = fluxDeriv
                  case (1)
                     dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
                     do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                        dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                                     + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,l,j)
                     end do                  ; end do                   ; end do
                     
                  case (2)
                     dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
                     do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                        dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                                     + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,i,l)
                     end do                  ; end do                   ; end do
            
                  case (3)
                     dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
                     do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                        do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                           dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                                  +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                     * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                     * fluxDeriv(:,:,m,l)
                        end do                 ; end do
                     end do                  ; end do
                  end select
                  end associate
               end if 
   
            case (RIGHT)    ! Prolong to right element
   !      
   !           *********
   !           1st stage: Projection
   !           *********
   !       
               if (self % IsMortar==0 .OR. self % IsMortar==2) then
                  select case ( self % projectionType(2) )
                  case (0)
                     fStarAux(1:nEqn,1:nEqn,:,:) = fluxDeriv
                  case (1)
                     fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
                     do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                        fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,l,j)
                     end do                  ; end do                   ; end do
                     
                  case (2)
                     fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
                     do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                        fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,i,l)
                     end do                  ; end do                   ; end do
            
                  case (3)
                     fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
                     do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                        do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                           fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                                  * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                                  * fluxDeriv(:,:,m,l)
                        end do                 ; end do
                     end do                  ; end do
                  end select
      !      
      !           *********
      !           2nd stage: Rotation
      !           *********
      !      
                  associate(dFStar_dq => self % storage(2) % dFStar_dqEl)
                  do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                     call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
                     dFStar_dq(1:nEqn,1:nEqn,ii,jj,whichderiv) = fStarAux(1:nEqn,1:nEqn,i,j) 
                  end do                        ; end do
      !
      !           *********
      !           3rd stage: Inversion
      !           *********
      !
                  dFStar_dq = -dFStar_dq
                  end associate
               end if 
         end select
         
   
      end subroutine Face_ProjectFluxJacobianToElements

      subroutine Face_ProjectMortarFluxJacobianToElements(self, nEqn, whichElement,whichderiv, fma, fmb, fmc, fmd) !!!ok
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face), target        :: self
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: whichElement
         integer,       intent(in)  :: whichderiv           !<  One can either transfer the derivative with respect to qL (LEFT) or to qR (RIGHT)
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer                :: i, j, l, m, side, lm, p, q
         real(kind=RP), pointer :: fluxDeriv(:,:,:,:)
         !write(*,*)'hum flux jacobian...'
         !luxDeriv(1:,1:,0:,0:) => self % storage(whichderiv) % dFStar_dqF
   
         real(kind=RP)     :: Flux_tmp(nEqn,nEqn,0: fma % NfLeft(1),0:fma % NfLeft(1))
         real(kind=RP)     :: Flux_master(nEqn,nEqn,0: fma % NfLeft(1),0:fma % NfLeft(1),2,4)
         real(kind=RP)     :: small(nEqn,nEqn,0: fma % NfLeft(1), 0:fma % NfLeft(1),4)
         real(kind=RP)     :: Mout(0: fma % NfLeft(1), 0:fma % NfLeft(1), 2)
   
   
          if (self % Ismortar==1) then 
            small(:,:,:,:,1) = fma  % storage(whichderiv) % dFStar_dqF
            small(:,:,:,:,2) = fmb % storage(whichderiv) % dFStar_dqF
            small(:,:,:,:,3) = fmc % storage(whichderiv) % dFStar_dqF
            small(:,:,:,:,4) = fmd % storage(whichderiv) % dFStar_dqF
            
             
            Mout(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 2) % T)
            Mout(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 2) % T)
            Flux_master=0.0_RP
            do lm=1,4
               do j=1,2
                  Flux_tmp(:,:,:,:)=0.0_RP
                  do i=1,2
                     do q=0,fma % NfLeft(1) ; do p=0,fma % NfLeft(1); do l=0,fma % NfLeft(1)
                        Flux_tmp(:,:,p,q)=Flux_tmp(:,:,p,q) + Mout(l,p,i)*small(:,:,l,q,i+2*(j-1))
                     end do ; end do ; end do 
                  end do 
                 
                  do p=0,fma % NfLeft(1) ; do q=0,fma % NfLeft(1) ; do l=0,fma % NfLeft(1)
                        Flux_master(1:nEqn,1:nEqn,p,q,whichderiv,lm)=Flux_master(1:nEqn,1:nEqn,p,q,whichderiv,lm) &
                         + Mout(l,q,j)*Flux_tmp(:,:,p,l) 
                  end do ; end do ; end do 
               end do 
            end do
   
            associate(dFStar_dq => self % storage(1) % dFStar_dqEl)
               dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv)=0.25_RP*(Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,1)+&
               Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,2)+ Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,4)+ &
               Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,3))
            end associate 
            associate(dFStar_dq => self % storage(2) % dFStar_dqEl)
               dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv)=0.25_RP*(Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,1)+&
               Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,2)+ Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,4)+ &
               Flux_master(1:nEqn,1:nEqn,:,:,whichderiv,3))
               dFStar_dq = -dFStar_dq
            end associate 
   
            end if 
         
      end subroutine Face_ProjectMortarFluxJacobianToElements
   
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !
#endif
      subroutine Face_ProjectGradientFluxToElements(self, nEqn, Hflux, whichElements, factor)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face)       :: self
         integer,       intent(in)  :: nEqn
         real(kind=RP), intent(in), optional  :: Hflux(nEqn, NDIM, 0:self % Nf(1), 0:self % Nf(2))
         integer,       intent(in)  :: whichElements(2)
         integer,       intent(in)  :: factor               ! A factor that relates LEFT and RIGHT fluxes
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer           :: i, j, ii, jj, l, m, side,p,q,lm
         real(kind=RP)     :: hStarAux(nEqn, NDIM, 0:self % NfRight(1), 0:self % NfRight(2))
         real(kind=RP), allocatable     :: Flux_tmp(:,:,:,:)
         real(kind=RP), allocatable     :: Flux_master(:,:,:,:,:)
         real(kind=RP), allocatable     ::small(:,:,:,:,:)
         real(kind=RP), allocatable     ::Mout( :,:,:)

         do side = 1, 2
            select case ( whichElements(side) )
            case (1)    ! project from left element
            if (self % IsMortar==0) then 
               associate(unStar => self % storage(1) % unStar)
               select case ( self % projectionType(1) )
               case (0)
                  unStar(:,:,:,:) = Hflux
               case (1)
                  unStar(:,:,:,:) = 0.0
                  do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                     unStar(:,:,i,j) = unStar(:,:,i,j) + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * Hflux(:,:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  unStar(:,:,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                     unStar(:,:,i,j) = unStar(:,:,i,j) + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * Hflux(:,:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  unStar(:,:,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                        unStar(:,:,i,j) = unStar(:,:,i,j) +   Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                * Hflux(:,:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
               end associate
            end if 
   
            case (2)    ! Prolong from right element
   !      
   !           *********
   !           1st stage: Projection
   !           *********
   !      
            if (self % IsMortar==0 .OR. self % IsMortar==2) then 
               select case ( self % projectionType(2) )
               case (0)
                  HstarAux = Hflux
               case (1)
                  HstarAux = 0.0
                  do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                     HstarAux(:,:,i,j) = HstarAux(:,:,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * Hflux(:,:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  HstarAux = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                     HstarAux(:,:,i,j) = HstarAux(:,:,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * Hflux(:,:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  HstarAux = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                        HstarAux(:,:,i,j) = HstarAux(:,:,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                                * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                                * Hflux(:,:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
   !      
   !           *********
   !           2nd stage: Rotation
   !           *********
   !      
               associate(unStar => self % storage(2) % unStar)
               do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
                  unStar(:,:,ii,jj) = HstarAux(:,:,i,j) 
               end do                        ; end do
   !
   !           *********
   !           3rd stage: Multiplication by a factor (inversion usually)
   !           *********
   !
               unStar = factor * unStar
   
               end associate
            end if 
            end select
         end do
   
      end subroutine Face_ProjectGradientFluxToElements

      subroutine Face_ProjectMortarGradientFluxToElements(self, nEqn, Hflux, whichElements, factor,fma,fmb,fmc,fmd, HfluxM1, HfluxM2, HfluxM3, HfluxM4)   !!!ok
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         class(Face)       :: self
         integer,       intent(in)  :: nEqn
         real(kind=RP), intent(in), optional  :: Hflux(nEqn, NDIM, 0:self % Nf(1), 0:self % Nf(2))
         integer,       intent(in)  :: whichElements(2)
         integer,       intent(in)  :: factor               ! A factor that relates LEFT and RIGHT fluxes
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
         real(kind=RP), intent(in)  :: HfluxM1(nEqn, NDIM, 0:fma % Nf(1), 0:fma % Nf(2))
         real(kind=RP), intent(in)  :: HfluxM2(nEqn, NDIM, 0:fma % Nf(1), 0:fma % Nf(2))
         real(kind=RP), intent(in)  :: HfluxM3(nEqn, NDIM, 0:fma % Nf(1), 0:fma % Nf(2))
         real(kind=RP), intent(in)  :: HfluxM4(nEqn, NDIM, 0:fma % Nf(1), 0:fma % Nf(2))
   
   !
   !     ---------------
   !     Local variables
   !     --------------
   !
         integer           :: i, j, l, m,p,q,lm
         real(kind=RP)     :: Flux_tmp(1:nEqn,NDIM, 0:fma % NfLeft(1), 0:fma % NfLeft(1))
         real(kind=RP)     :: Flux_master(1:nEqn,NDIM, 0: fma % NfLeft(1), 0:fma % NfLeft(1),4)
         real(kind=RP)     :: small(1:nEqn,NDIM,0: fma % NfLeft(1), 0:fma % NfLeft(1),4)
         real(kind=RP)     :: Mout(0: fma % NfLeft(1), 0:fma % NfLeft(1), 2)
   
         if (self % IsMortar==1) then 
            
            small(:,:,:,:,1)=HfluxM1
            small(:,:,:,:,2)=HfluxM2
            small(:,:,:,:,3)=HfluxM3
            small(:,:,:,:,4)=HfluxM4
   
            Mout(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 2) % T)
            Mout(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 2) % T)
            Flux_master=0.0_RP
            do lm=1,4
               do j=1,2
                  Flux_tmp(:,:,:,:)=0.0_RP
                  do i=1,2
                     do q=0,fma % NfLeft(1) ; do p=0,fma % NfLeft(1); do l=0,fma % NfLeft(1)
                        Flux_tmp(:,:,p,q)=Flux_tmp(:,:,p,q) + Mout(l,p,i)*small(:,:,l,q,i+2*(j-1))
                     end do ; end do ; end do 
                  end do 
                 
                  do p=0,fma % NfLeft(1) ; do q=0,fma % NfLeft(1) ; do l=0,fma % NfLeft(1)
                        Flux_master(:,:,p,q,lm)=Flux_master(:,:,p,q,lm) + Mout(l,q,j)*Flux_tmp(:,:,p,l) 
                  end do ; end do ; end do 
               end do 
            end do
   
            associate(unStar => self % storage(1) % unStar)
               unStar=0.25_RP*(Flux_master(:,:,:,:,1)+Flux_master(:,:,:,:,2)+Flux_master(:,:,:,:,4)+Flux_master(:,:,:,:,3))
            end associate 
            associate(unStar => self % storage(2) % unStar)
               unStar=0.25_RP*(Flux_master(:,:,:,:,1)+Flux_master(:,:,:,:,2)+Flux_master(:,:,:,:,4)+Flux_master(:,:,:,:,3))
               unStar = factor * unStar
            end associate 
         end if 
   
      end subroutine Face_ProjectMortarGradientFluxToElements
     
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !
#if defined(NAVIERSTOKES)
      subroutine Face_ProjectGradJacobianToElements(self, whichElement, whichderiv)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         !---------------------------------------------------------
         class(Face), target        :: self
         integer,       intent(in)  :: whichElement
         integer,       intent(in)  :: whichderiv
         !---------------------------------------------------------
         integer                :: i, j, ii, jj, l, m, side, lm, a, b
         real(kind=RP), pointer :: fluxDeriv(:,:,:,:,:)
         real(kind=RP)          :: fStarAux(NCONS,NCONS,1:NDIM, 0:self % NfRight(1), 0:self % NfRight(2))
   
         
         fluxDeriv(1:,1:,1:,0:,0:) => self % storage(whichderiv) % dFv_dGradQF
         
         select case ( whichElement )
            case (LEFT)    ! Prolong to left element
            if (self % IsMortar==0) then 
               associate(dFv_dGradQEl => self % storage(1) % dFv_dGradQEl)
               
               select case ( self % projectionType(1) )
               case (0)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = fluxDeriv
               case (1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
                  do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                     dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                             + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
                  do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                     dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                             + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                        dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                                      +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                                       * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                                       * fluxDeriv(:,:,:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
               
               end associate
            end if 
   
            case (RIGHT)    ! Prolong to right element
            if (self % IsMortar==0 .OR. self % IsMortar==2) then
   !      
   !           *********
   !           1st stage: Projection
   !           *********
   !      
               select case ( self % projectionType(2) )
               case (0)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = fluxDeriv
               case (1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
                  do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                     fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                     fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                        fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                                * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                                * fluxDeriv(:,:,:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
   !      
   !           ********* 
   !           2nd stage: Rotation
   !           *********
   !      
               associate(dFv_dGradQEl => self % storage(2) % dFv_dGradQEl)
               
               do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,ii,jj,whichderiv) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) 
               end do                        ; end do
   !
   !           *********
   !           3rd stage: Inversion
   !           *********
   !
               dFv_dGradQEl = -dFv_dGradQEl
               
               end associate
            end if 
         end select
         
         
      end subroutine Face_ProjectGradJacobianToElements
      
      subroutine Face_ProjectMortarGradJacobianToElements(self, whichElement, whichderiv, fma, fmb, fmc, fmd)
         use MappedGeometryClass
         use PhysicsStorage
         implicit none
         !---------------------------------------------------------
         class(Face), target        :: self
         integer,       intent(in)  :: whichElement
         integer,       intent(in)  :: whichderiv
         type(Face), intent(inout)      ::fma
         type(Face), intent(inout)      ::fmb
         type(Face), intent(inout)      :: fmc
         type(Face), intent(inout)      :: fmd 
         !---------------------------------------------------------
   
         real(kind=RP), pointer :: fluxDeriv(:,:,:,:,:)
   
         integer           :: i, j, l, m, p, q, lm
         real(kind=RP)     :: Flux_tmp(1:NCONS,1:NCONS,1:NDIM,0: fma % NfLeft(1),0:fma % NfLeft(1))
         real(kind=RP)     :: Flux_master(1:NCONS,1:NCONS,1:NDIM,0: fma % NfLeft(1),0:fma % NfLeft(1),2,4)
         real(kind=RP)     :: Mout(0: fma % NfLeft(1), 0:fma % NfLeft(1), 2)
         real(kind=RP)     :: small(1:NCONS,1:NCONS,NDIM,0: fma % NfLeft(1), 0:fma % NfLeft(1),4)
         
          if (self % Ismortar==1) then 
            small(:,:,:,:,:,1)=fma % storage(whichderiv) % dFv_dGradQF
            
            small(:,:,:,:,:,2) = fmb % storage(whichderiv) % dFv_dGradQF

            small(:,:,:,:,:,3) = fmc % storage(whichderiv) % dFv_dGradQF

            small(:,:,:,:,:,4) = fmd % storage(whichderiv) % dFv_dGradQF


            Mout(:,:,1)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 1, 2) % T)
            Mout(:,:,2)=transpose(TsetM(fma % NfLeft(1), fma % Nf(1), 2, 2) % T)
            Flux_master=0.0_RP
            do lm=1,4
               do j=1,2
                  Flux_tmp(:,:,:,:,:)=0.0_RP
                  do i=1,2
                     do q=0,fma % NfLeft(1) ; do p=0,fma % NfLeft(1); do l=0,fma % NfLeft(1)
                        Flux_tmp(:,:,:,p,q)=Flux_tmp(:,:,:,p,q) + Mout(l,p,i)*small(:,:,:,l,q,i+2*(j-1))
                     end do ; end do ; end do 
                  end do 
                 
                  do p=0,fma % NfLeft(1) ; do q=0,fma % NfLeft(1) ; do l=0,fma % NfLeft(1)
                        Flux_master(1:NCONS,1:NCONS,1:NDIM,p,q,whichderiv,lm)=Flux_master(1:NCONS,1:NCONS,1:NDIM,p,q,whichderiv,lm) &
                         + Mout(l,q,j)*Flux_tmp(:,:,:,p,l) 
                  end do ; end do ; end do 
               end do 
            end do
   
            associate(dFv_dGradQEl => self % storage(1) % dFv_dGradQEl)
               dFv_dGradQEl=0.25_RP*(Flux_master(:,:,:,:,:,:,1)+Flux_master(:,:,:,:,:,:,2)+Flux_master(:,:,:,:,:,:,4)+Flux_master(:,:,:,:,:,:,3))
            end associate 
            associate(dFv_dGradQEl => self % storage(2) % dFv_dGradQEl)
               dFv_dGradQEl=0.25_RP*(Flux_master(:,:,:,:,:,:,1)+Flux_master(:,:,:,:,:,:,2)+Flux_master(:,:,:,:,:,:,4)+Flux_master(:,:,:,:,:,:,3))
            dFv_dGradQEl = -dFv_dGradQEl
            end associate 
         end if 
         
      end subroutine Face_ProjectMortarGradJacobianToElements

#endif
   !
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !
      SUBROUTINE InterpolateToBoundary( u, v, Nx, Ny, Nz, which_dim , bValue , NEQ)
         use SMConstants
         use PhysicsStorage
         IMPLICIT NONE
         integer                      , INTENT(IN)    :: Nx, Ny, Nz
         real(kind=RP)                , intent(in)    :: u(1:NEQ,0:Nx,0:Ny,0:Nz) , v(0:)
         integer                      , intent(in)    :: which_dim
         REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
         integer                      , intent(in)    :: NEQ
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
         integer                                    :: i , j , k , eq
   
         select case (which_dim)
         case (IX)
   
            do j = 0 , Nz ; do i = 0 , Ny ; do k = 0 , Nx
               bValue(:,i,j) = bValue(:,i,j) + u(:,k,i,j) * v(k)
            end do        ; end do        ; end do
   
         case (IY)
   
            do j = 0 , Nz ; do k = 0 , Ny ; do i = 0 , Nx
               bValue(:,i,j) = bValue(:,i,j) + u(:,i,k,j) * v(k)
            end do        ; end do        ; end do
   
         case (IZ)
   
            do k = 0 , Nz ; do j = 0 , Ny ; do i = 0 , Nx
               bValue(:,i,j) = bValue(:,i,j) + u(:,i,j,k) * v(k)
            end do        ; end do        ; end do
   
         end select
         
         end SUBROUTINE InterpolateToBoundary
   !
   !////////////////////////////////////////////////////////////////////////
   !
         elemental subroutine Face_Assign(to,from)
            implicit none
            class(Face), intent(inout) :: to
            class(Face), intent(in)    :: from
            
            
            to % flat = from % flat
            to % ID = from % ID
            to % FaceType = from % FaceType
            to % zone = from % zone
            to % rotation = from % rotation
            to % NelLeft = from % NelLeft
            to % NelRight = from % NelRight
            to % NfLeft = from % NfLeft
            to % NfRight = from % NfRight
            to % Nf = from % Nf
            to % nodeIDs = from % nodeIDs
            to % elementIDs = from % elementIDs
            to % elementSide = from % elementSide
            to % projectionType = from % projectionType
            to % boundaryName = from % boundaryName
            to % geom = from % geom
            to % storage = from % storage
            to % IsMortar = from % IsMortar
         end subroutine Face_Assign

   end Module FaceClass
   