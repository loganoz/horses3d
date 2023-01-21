#include "Includes.h"
MODULE NodalStorageClass
   USE SMConstants
   USE PolynomialInterpAndDerivsModule
   USE GaussQuadrature
   use FTValueDictionaryClass          , only: FTValueDictionary
   use mainKeywordsModule              , only: discretizationNodesKey
   IMPLICIT NONE

   private
   public GAUSS, GAUSSLOBATTO                                              ! parameters
   public NodalStorage_t                                                   ! type
   public NodalStorage, NodalStorage_Gauss, NodalStorage_GaussLobatto      ! Nodal storage variables
   public InitializeNodalStorage, DestructGlobalNodalStorage, CurrentNodes ! Main nodal storage used in the simulation

   integer, parameter      :: GAUSS = 1
   integer, parameter      :: GAUSSLOBATTO = 2

!
!  -------------------
!  Nodal storage class
!  -------------------
!
   type NodalStorage_t
      logical                                    :: Constructed = .FALSE.     ! Constructed flag
      integer                                    :: nodes                     ! Either GAUSS or GAUSSLOBATTO
      integer                                    :: N                         ! Polynomial order
      real(kind=RP), dimension(:)  , allocatable :: x                         ! Node position
      real(kind=RP), dimension(:)  , allocatable :: w                         ! Weights
      real(kind=RP), dimension(:)  , allocatable :: wb                        ! Barycentric weights
      real(kind=RP), dimension(:,:), allocatable :: v                         ! Boundary interpolation vector
      real(kind=RP), dimension(:,:), allocatable :: b                         ! Boundary interpolation vector scaled with weight
      real(kind=RP), dimension(:,:), allocatable :: vd                        ! Boundary derivative vector
      real(kind=RP), dimension(:,:), allocatable :: bd                        ! Boundary derivative vector scaled with weight
      real(kind=RP), dimension(:,:), allocatable :: D                         ! DG derivative matrix
      real(kind=RP), dimension(:,:), allocatable :: DT                        ! Transposed DG derivative matrix
      real(kind=RP), dimension(:,:), allocatable :: hatD                      ! Weak form derivative matrix
      real(kind=RP), dimension(:,:), allocatable :: hatG                      ! Weak form Laplacian derivative matrix hatG := W⁻¹DᵀWD (where W is the diagonal matrix with the quadrature weights)
      real(kind=RP), dimension(:,:), allocatable :: sharpD                    ! (Two times) the strong form derivative matrix
      real(kind=RP), dimension(:,:), allocatable :: Fwd                       ! Projection matrix from Lagrange to Legendre
      real(kind=RP), dimension(:,:), allocatable :: Bwd                       ! Projection matrix from Legendre to Lagrange
      real(kind=RP), dimension(:)  , allocatable :: Lw                        ! Norm of the Legendre polynomials, ||L_i||^2
      real(kind=RP), dimension(:)  , allocatable :: xCGL, wbCGL
      real(kind=RP), dimension(:,:), allocatable :: DCGL, TCheb2Gauss
      contains
         procedure :: construct => ConstructNodalStorage
         procedure :: destruct  => DestructNodalStorage
         procedure :: lj        => NodalStorage_getlj
         procedure :: dlj       => NodalStorage_getdlj

   END TYPE NodalStorage_t

!  ------------------------------------------------
!  NodalStorage contains the nodal storage information
!  for every possible polynomial order of the mesh
!  ------------------------------------------------
   type(NodalStorage_t), target, allocatable :: NodalStorage_Gauss(:)
   type(NodalStorage_t), target, allocatable :: NodalStorage_GaussLobatto(:)


   type(NodalStorage_t), pointer :: NodalStorage(:)   ! Default nodal storage
   integer  :: CurrentNodes

   interface InitializeNodalStorage
      module procedure InitializeNodalStorage_controlVars, InitializeNodalStorage_nodeType
   end interface InitializeNodalStorage

!
!     ========
      CONTAINS
!     ========
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine InitializeNodalStorage_controlVars(controlVariables, Nmax)
      implicit none
      !---------------------------------------
      type(FTValueDictionary), intent(in) :: controlVariables
      integer                , intent(in) :: Nmax
      !---------------------------------------

      select case ( trim(controlVariables % stringValueForKey(trim(discretizationNodesKey), requestedLength = LINE_LENGTH)) )
         case("Gauss")
            safedeallocate(NodalStorage_Gauss)
            allocate ( NodalStorage_Gauss (0:Nmax) )

            NodalStorage => NodalStorage_Gauss
            CurrentNodes = GAUSS

         case("Gauss-Lobatto")
            safedeallocate(NodalStorage_GaussLobatto)
            allocate ( NodalStorage_GaussLobatto (0:Nmax) )

            NodalStorage => NodalStorage_GaussLobatto
            CurrentNodes = GAUSSLOBATTO
         case default
            print*, "Unknown discretization nodes."
            print*, "Options available are:"
            print*, "   * Gauss"
            print*, "   * Gauss-Lobatto"
            errorMessage(STD_OUT)
            stop
      end select

   end subroutine InitializeNodalStorage_controlVars

   subroutine InitializeNodalStorage_nodeType(nodeType, Nmax)
      implicit none
      !---------------------------------------
      integer, intent(in) :: nodeType
      integer, intent(in) :: Nmax
      !---------------------------------------

      CurrentNodes = nodeType

      select case ( nodeType )
         case(GAUSS)
            safedeallocate(NodalStorage_Gauss)
            allocate ( NodalStorage_Gauss (0:Nmax) )

            NodalStorage => NodalStorage_Gauss


         case(GAUSSLOBATTO)
            safedeallocate(NodalStorage_GaussLobatto)
            allocate ( NodalStorage_GaussLobatto (0:Nmax) )

            NodalStorage => NodalStorage_GaussLobatto
         case default
            print*, "Unknown discretization nodes."
            print*, "Options available are:"
            print*, "   * Gauss"
            print*, "   * Gauss-Lobatto"
            errorMessage(STD_OUT)
            stop
      end select

   end subroutine InitializeNodalStorage_nodeType
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine DestructGlobalNodalStorage
      implicit none
      !---------------------
      integer :: k
      !---------------------

      if ( allocated(NodalStorage_Gauss) ) then
         do k=lbound(NodalStorage_Gauss,1), ubound(NodalStorage_Gauss,1)
            IF (.NOT. NodalStorage_Gauss(k) % Constructed) cycle
            call NodalStorage_Gauss(k) % destruct()
         end do
         deallocate (NodalStorage_Gauss)
      end if

      if ( allocated(NodalStorage_GaussLobatto) ) then
         do k=lbound(NodalStorage_GaussLobatto,1), ubound(NodalStorage_GaussLobatto,1)
            IF (.NOT. NodalStorage_GaussLobatto(k) % Constructed) cycle
            call NodalStorage_GaussLobatto(k) % destruct()
         end do
         deallocate (NodalStorage_GaussLobatto)
      end if

      nullify (NodalStorage)
   end subroutine DestructGlobalNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructNodalStorage( this, nodes, N)
      implicit none
      class(NodalStorage_t)    :: this       !<> Nodal storage being constructed
      integer, intent(in)      :: nodes
      integer, intent(in)      :: N          !<  Polynomial order
      !--------------------------------------
      integer            :: i,j,k
      integer, PARAMETER :: LEFT = 1, RIGHT = 2
      real(kind=RP)      :: wb(0:N)
      real(kind=RP)      :: Lkj(0:N,0:N), dLk_dummy
      !--------------------------------------

      if (this % Constructed) return

      this % nodes = nodes
      this % N = N

      ALLOCATE( this % x     (0:N) )
      ALLOCATE( this % w     (0:N) )
      ALLOCATE( this % wb    (0:N) )
      ALLOCATE( this % v     (0:N,2) )
      ALLOCATE( this % b     (0:N,2) )
      ALLOCATE( this % vd    (0:N,2) )
      ALLOCATE( this % bd    (0:N,2) )
      ALLOCATE( this % D     (0:N,0:N) )
      ALLOCATE( this % DT    (0:N,0:N) )
      ALLOCATE( this % hatD  (0:N,0:N) )
      ALLOCATE( this % hatG  (0:N,0:N) )
      ALLOCATE( this % Fwd   (0:N,0:N) )
      ALLOCATE( this % Bwd   (0:N,0:N) )
      ALLOCATE( this % Lw    (0:N) )
      ALLOCATE( this % xCGL  (0:N) )
      ALLOCATE( this % wbCGL (0:N) )
      ALLOCATE( this % DCGL  (0:N,0:N) )
      ALLOCATE( this % TCheb2Gauss (0:N,0:N) )
!
!     -----------------
!     Nodes and weights
!     -----------------
!
      select case (this % nodes)
         case (GAUSS)
            CALL GaussLegendreNodesAndWeights  ( N, this % x  , this % w )
         case (GAUSSLOBATTO)
            if ( N .ne. 0 ) then
               CALL LegendreLobattoNodesAndWeights( N, this % x  , this % w )

            else
               this % x = 0.0_RP
               this % w = 2.0_RP

            end if

            allocate( this % sharpD(0:N,0:N) )
         case default
            print*, "Undefined nodes choice"
            errorMessage(STD_OUT)
            stop
      end select
!
!     -----------------
!     Derivative Matrix
!     -----------------
!
      CALL PolynomialDerivativeMatrix( N, this % x, this % D )
      DO j = 0, N
         DO i = 0, N
            this % DT  (i,j) = this % D (j,i)
            this % hatD(i,j) = this % DT(i,j) * this % w(j) / this % w(i)
         END DO
      END DO

      this % hatG = matmul(this % hatD, this % D)
!
!     --------------------------------------------------------------
!     Construct the strong form derivative matrices (Skew-Symmetric)
!     --------------------------------------------------------------
!
      if ( this % nodes .eq. GAUSSLOBATTO ) then
         if ( N .ne. 0 ) then
            this % sharpD = 2.0_RP * this % D
            this % sharpD(0,0) = 2.0_RP * this % D(0,0) + 1.0_RP / this % w(0)
            this % sharpD(N,N) = 2.0_RP * this % D(N,N) - 1.0_RP / this % w(N)

         else
            this % sharpD = 0.0_RP

         end if
      end if
!
!     ---------------------
!     Interpolation vectors
!     ---------------------
!
      CALL BarycentricWeights( N, this % x, wb )

      CALL InterpolatingPolynomialVector(  1.0_RP, N, this % x, wb, this % b(:,RIGHT) )
      CALL InterpolatingPolynomialVector( -1.0_RP, N, this % x, wb, this % b(:,LEFT)  )

      this % wb = wb
      this % v  = this % b

      this % b(0:N,LEFT)  = this % b(0:N,LEFT) /this % w
      this % b(0:N,RIGHT) = this % b(0:N,RIGHT)/this % w

!
!     ------------------
!     Derivative vectors
!     ------------------
!
      CALL PolyDerivativeVector(  1.0_RP, N, this % x, this % bd(:,RIGHT) )
      CALL PolyDerivativeVector( -1.0_RP, N, this % x, this % bd(:,LEFT)  )

      this % vd = this % bd

      this % bd(0:N,LEFT)  = this % bd(0:N,LEFT) /this % w
      this % bd(0:N,RIGHT) = this % bd(0:N,RIGHT)/this % w
!
!     -------------------------------------------
!     Construct Chebyshev-Gauss-Lobatto framework
!     -------------------------------------------
!
      this % DCGL = 0.0_RP
      this % TCheb2Gauss = 0.0_RP

      if ( N .ne. 0 ) then
         this % xCGL = (/ (-cos(1.0_RP*i*PI/this % N), i = 0, this % N) /)

      else
         this % xCGL = 0.0_RP

      end if

      call BarycentricWeights(N, this % xCGL, this % wbCGL)
      call PolynomialDerivativeMatrix( this % N, this % xCGL, this % DCGL)
      call PolynomialInterpolationMatrix(this % N, this % N, this % xCGL, this % wbCGL, this % x,&
                                         this % TCheb2Gauss)
!
!     ---------------------------------------------------------------------
!     Construct projection matrices from/to Lagrange to/from Legendre basis
!     ---------------------------------------------------------------------
!
!     Get the evaluation of Legendre polynomials at the interpolation nodes
      do j = 0 , N ;    do k = 0 , N
         call LegendrePolyAndDerivative(k, this % x(j), Lkj(k,j), dLk_dummy)
      end do       ;    end do
!
!     Get the norm of Legendre polynomials
      do k = 0 , N
         this % Lw(k) = 0.0_RP
         do j = 0 , N
            this % Lw(k) = this % Lw(k) + this % w(j) * Lkj(k,j) * Lkj(k,j)
         end do
      end do
!
!     Get the transformation from Nodal to Modal and vice-versa matrices
      do k = 0 , N   ; do i = 0 , N
         this % Fwd(k,i) = this % w(i) * Lkj(k,i) / this % Lw(k)
         this % Bwd(i,k) = Lkj(k,i)
      end do         ; end do
!
!     Set the nodal storage as constructed
!     ------------------------------------
      this % Constructed = .TRUE.

   END SUBROUTINE ConstructNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
   elemental subroutine DestructNodalStorage( this )
      IMPLICIT NONE
      class(NodalStorage_t), intent(inout) :: this
!
!     Attempting to destruct a non-constructed nodal storage
!     ------------------------------------------------------
      if (.not. this % constructed ) return
!
!     Destruct otherwise
!     ------------------
      this % constructed = .FALSE.

      DEALLOCATE( this % x )
      DEALLOCATE( this % w  )
      DEALLOCATE( this % wb  )
      DEALLOCATE( this % D  )
      DEALLOCATE( this % DT )
      DEALLOCATE( this % hatD)
      DEALLOCATE( this % hatG)
      DEALLOCATE( this % xCGL )
      DEALLOCATE( this % wbCGL )
      DEALLOCATE( this % DCGL )
      DEALLOCATE( this % TCheb2Gauss )
      DEALLOCATE( this % v )
      DEALLOCATE( this % b )
      DEALLOCATE( this % vd )
      DEALLOCATE( this % bd )
      DEALLOCATE( this % Fwd )
      DEALLOCATE( this % Bwd )
      DEALLOCATE( this % Lw )
      safedeallocate( this % sharpD )    !  This matrices are just generated for Gauss-Lobatto discretizations.

      end subroutine DestructNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
      function NodalStorage_getlj(self, xi)
         implicit none
         class(NodalStorage_t), intent(in)  :: self
         real(kind=RP),       intent(in)  :: xi
         real(kind=RP)                    :: NodalStorage_getlj(0:self % N)

         call InterpolatingPolynomialVector(xi , self % N , self % x , self % wb , NodalStorage_getlj)

      end function NodalStorage_getlj

      function NodalStorage_getdlj(self, xi)
         implicit none
         class(NodalStorage_t), intent(in)  :: self
         real(kind=RP),       intent(in)  :: xi
         real(kind=RP)                    :: NodalStorage_getdlj(0:self % N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i

         do i = 0 , self % N
            NodalStorage_getdlj(i) = EvaluateLagrangePolyDerivative( i, xi, self % N , self % x)
         end do

      end function NodalStorage_getdlj

END Module NodalStorageClass
