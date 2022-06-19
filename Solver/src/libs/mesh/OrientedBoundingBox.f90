!
!//////////////////////////////////////////////////////
!
!   @File:    OrientedBoundingBox.f90
!   @Author:  Stefano Colombo (stefano.colombo@upm.es) 
!   @Created: Mon Aug  9 11:44:25 2018
!   @Last revision date: Mon Aug 09 17:24:24 2021
!   @Last revision author: Stefano Colombo (stefano.colombo@upm.es)
!   @Last revision commit: 
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module OrientedBoundingBox

   use SMConstants
   use Utilities
   use TessellationTypes

   implicit none   
   
   integer,       parameter :: TASK_THRESHOLD = 10000
   real(kind=RP), parameter :: SAFETY_FACTOR = 0.001_RP

   public OBB

!
!  **************************************************
!  Main type for the Convex Hull computations
!  **************************************************
   
   type Hull_type
    
      type(point_type), dimension(:), allocatable :: Points 
      integer                                     :: NumOfPoints

   end type
!
!  **************************************************
!  Main type for a rectangle
!  **************************************************
   type rectangle
      
      real(kind=rp), dimension(2,4)  :: Vertices
      real(kind=rp), dimension(2)    :: Center
      real(kind=rp), dimension(NDIM) :: normal, t1, t2
      real(kind=rp)                  :: Length, Width, Angle
   
      contains
         procedure :: ComputeVertices
   
   end type   
!
!  **************************************************
!  Main type for the Oriented Bounding Box computations
!  **************************************************
   type OBB_type
      
      type(point_type), dimension(:), allocatable :: Points, HullPoints
      type(rectangle)                             :: MBR
      real(kind=rp),    dimension(NDIM,8)         :: vertices, LocVertices
      real(kind=rp),    dimension(NDIM,NDIM)      :: R, invR
      real(kind=rp),    dimension(NDIM)           :: CloudCenter, LocFrameCenter
      real(kind=rp)                               :: nMin, nMax
      integer                                     :: NumOfPoints, center, left, right
      character(len=LINE_LENGTH)                  :: filename
      logical                                     :: verbose
      
      contains
         procedure :: construct             => OBB_construct
         procedure :: ReadStorePoints       => OBB_ReadStorePoints
         procedure :: ComputeAngle          => OBB_ComputeAngle
         procedure :: SortingNodes          => OBB_SortingNodes
         procedure :: isPointInside         => OBB_isPointInside
         procedure :: ChangeObjsRefFrame    => OBB_ChangeObjsRefFrame
         procedure :: STL_rotate            => OBB_STL_rotate
         procedure :: STL_translate         => OBB_STL_translate
         procedure :: ChangeRefFrame
         procedure :: ComputeRotationMatrix 
         procedure :: plot                  => OBB_plot
         procedure :: isInsidePolygon       => OBB_isInsidePolygon

   end type
   
   type(OBB_type), allocatable :: OBB(:)

contains
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the determinant of a 3x3 Matrix. 
!  -------------------------------------------------
   
   function Determinant( Mat ) result( det )

      implicit none
      !-arguments-----------------------------
      real(kind=rp), dimension(:,:), intent(in) :: Mat
      !-local-variables-------------------------
      real(kind=rp) :: det

      det = Mat(1,1)*( Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2) ) - &
            Mat(1,2)*( Mat(2,1)*Mat(3,3) - Mat(2,3)*Mat(3,1) ) + &
            Mat(1,3)*( Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1) )

   end  function Determinant
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the eigenvalues & eigenvectors of a 3x3 matrix. 
! For details, see https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf .
!  -------------------------------------------------
   
   subroutine ComputeEigenStructure( A, EigenVal, EigenVec1, EigenVec2, EigenVec3 )
      use MappedGeometryClass
      implicit none
      !-arguments----------------------------------------------
      real(kind=rp), dimension(:,:),   intent(in)  :: A
      real(kind=rp),                   intent(out) :: EigenVal
      real(kind=rp), dimension(NDIM), intent(out) :: EigenVec1, EigenVec2, EigenVec3
      !-local-variables------------------------------------------
      real(kind=rp),  dimension(NDIM,NDIM) :: Ident, B, EigenVecs
       real(kind=rp), dimension(NDIM)      :: EigenVals
      real(kind=rp)                         :: p1, p2, p, q, r, phi
      integer                               :: i
      integer, dimension(1)                 :: min_loc, max_loc

      p1 = POW2( A(1,2) ) + POW2( A(1,3) ) + POW2( A(2,3) )    
      
!      
!     Diagonal Matrix
!     --------------
      if( AlmostEqual(p1,0.0_RP) ) then
         EigenVals(1) = A(1,1)
         EigenVals(2) = A(2,2)
         EigenVals(3) = A(3,3)
         min_loc   = minloc(EigenVals)
         max_loc   = maxloc(EigenVals)
         EigenVal  = EigenVals(min_loc(1))
         EigenVec1 = 0.0_RP; EigenVec2 = 0.0_RP; EigenVec3 = 0.0_RP
         EigenVec3(min_loc(1)) = 1.0_RP
         EigenVec1(max_loc(1)) = 1.0_RP
         EigenVals(min_loc(1)) = huge(1.0_RP)
         min_loc = minloc(EigenVals)
         EigenVec2(min_loc(1)) = 1.0_RP
         return
      end if
     
      q = 0._RP; Ident = 0._RP
     
      do i = 1, NDIM
         q         = q + A(i,i)
         Ident(i,i) = 1.0_RP
      end do
      
      q = q/3.0_RP
      
      p2 = POW2( A(1,1) - q ) + POW2( A(2,2) - q ) + POW2( A(3,3) - q ) + 2.0_RP*p1
          
      p = sqrt(p2/6.0_RP) !sqrt{ [tr(A-qI)^2]/6 }
!
!     Matrix B has the same eigenstraucture of A since it's translated by a factor = q
!     ----------------------------------------------------------------------
      B = ( 1.0_RP/p ) * ( A - q * Ident )
          
      r = Determinant(B)/2.0_RP
      
      if( r <= -1.0_RP ) then
         phi = PI/3.0_RP
      elseif( r >= 1 ) then
         phi = 0.0_RP
      else
         phi = acos(r)/3.0_RP
      end if
!
!     Ordered eigen values \lambda_1 > \lambda_2 > \lambda_3
!
      EigenVals(1)  = q + 2.0_RP * p * cos(phi)
      EigenVals(3) = q + 2.0_RP * p * cos(phi + 2.0_RP/3.0_RP*PI )
      EigenVals(2) = 3.0_RP * q - EigenVals(1) - Eigenvals(3)
      
      ! (A - \lambda_i I)· v_i = 0 ==> v_i*· (A - \lambda_i I)e_j = 0 with e_j unit basis vector
      ! v_i = (A^j - \lambda_i e_j) x (A^k - \lambda_i e_k)
      ! v_j = (A^i - \lambda_j e_i) x (A^k - \lambda_j e_k)
      ! v_k = v_i x v_j      
      
      call vcross( A(:,2)-EigenVals(1)*(/ 0.0_RP,1.0_RP,0.0_RP /), &
                   A(:,3)-EigenVals(1)*(/ 0.0_RP,0.0_RP,1.0_RP /), &
                   EigenVecs(:,1)                                  )
      call vcross( A(:,1)-EigenVals(2)*(/ 1.0_RP,0.0_RP,0.0_RP /), &
                   A(:,3)-EigenVals(2)*(/ 0.0_RP,0.0_RP,1.0_RP /), &
                   EigenVecs(:,2)                                  )
      call vcross( EigenVecs(:,1),EigenVecs(:,2), EigenVecs(:,3))

      
      EigenVecs(:,1) = EigenVecs(:,1)/norm2(EigenVecs(:,1))
      EigenVecs(:,2) = EigenVecs(:,2)/norm2(EigenVecs(:,2))
      EigenVecs(:,3) = EigenVecs(:,3)/norm2(EigenVecs(:,3))

!
!     Saving the last eigenvalue
!     ------------------------
      EigenVal = EigenVals(3)
      
!
!     Saving the Eigenvectors
!     ------------------------
     EigenVec1 = EigenVecs(:,1); EigenVec2 = EigenVecs(:,2); EigenVec3 = EigenVecs(:,3)   
       
   end  subroutine ComputeEigenStructure
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!   
!  -------------------------------------------------------------------------------------------
!  Subroutine that computes the vertices of a rectangle
!  -------------------------------------------------------------------------------------------
   subroutine ComputeVertices( this )
      
      implicit none
      !-arguments----------------------
      class(rectangle), intent(inout) :: this
      
      this% vertices(:,1) = 0.5_RP*(/ -this% Length, -this% Width /)
      this% vertices(:,2) = 0.5_RP*(/ this% Length, -this% Width /)
      this% vertices(:,3) = 0.5_RP*(/ this% Length, this% Width /)
      this% vertices(:,4) = 0.5_RP*(/ -this% Length, this% Width /)
      
   end subroutine ComputeVertices
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!   
!  -------------------------------------------------------------------------------------------
!  Subroutine for reading the nodes coming from a .dat file. In the future it will a .stl/.obj file
!  -------------------------------------------------------------------------------------------
   subroutine OBB_ReadStorePoints( this, stl )

      implicit none
      !-arguments----------------------
      class(OBB_type), intent(inout) :: this
      type(STLfile),   intent(in)    :: stl
      !-local.variables-------------------
      integer                        :: i, j, n
      
      this% NumOfPoints = 3*stl% NumOfObjs
      
      allocate(this% Points(this% NumOfPoints))

      n = 0

      associate( Objs => stl% ObjectsList )

      do i = 1, stl% NumOfObjs
         do j = 1, Objs(i)% NumOfVertices
            n = n + 1
            this% Points(n)% coords = Objs(i)% vertices(j)% coords
            this% Points(n)% index = n
         end do 
      end do

      end associate

      this% filename = stl% filename

   end subroutine OBB_ReadStorePoints
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!   
!  -------------------------------------------------------------------------------------------
!  Subroutine for plotting the OBB
!  -------------------------------------------------------------------------------------------
   subroutine OBB_plot( this )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------
      class(OBB_type), intent(inout) :: this
      !-local-variables----------------------
      integer                        :: i, funit

      if( .not. MPI_Process% isRoot ) return

      funit = UnusedUnit()

      open(funit,file='IBM/OrientedBoundingBox_'//trim(this% filename)//'.tec', status='unknown')

      write(funit,"(a28)") 'TITLE = "OrientedBoudingBox"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'
    
      do i = 1, 8
         write(funit,'(3E13.5)') Lref*this% vertices(1,i), Lref*this% vertices(2,i), Lref*this% vertices(3,i)
      end do 

      write(funit,'(4i2)') 1, 2, 3, 4
      write(funit,'(4i2)') 1, 5, 8, 4
      write(funit,'(4i2)') 5, 6, 7, 8
      write(funit,'(4i2)') 2, 3, 7, 6
      write(funit,'(4i2)') 4, 8, 7, 3
      write(funit,'(4i2)') 1, 2, 6, 5

      close(unit=funit)

   end subroutine OBB_plot

   subroutine ChangeRefFrame(this, v, FRAME, vNew)

      implicit none
      !-arguments----------------------------------------------------
      class(OBB_type),                 intent(inout) :: this
      real(kind=rp), dimension(:),     intent(in)    :: v
      character(len=*),                intent(in)    :: FRAME
      real(kind=rp), dimension(NDIM), intent(out)   :: vNew
      !-local-variables----------------------------------------------
      real(kind=rp), dimension(NDIM)      :: b
      real(kind=rp), dimension(NDIM,NDIM) :: T, invT

      T = 0.0_RP
      T(NDIM,NDIM) = 1.0_RP
      T(1,1) = cos(this% MBR% Angle); T(2,2)  = T(1,1)
      T(1,2) = sin(this% MBR% Angle); T(2,1) = -T(1,2)
            
      invT(:,1) = T(1,:)
      invT(:,2) = T(2,:)
      invT(:,3) = T(3,:)
      
      select case( trim(FRAME) )
         
         case('local')
         
!~             b = matmul( this% invR,(v - this% CloudCenter))
!~             b(1:2) = b(1:2) - this% MBR% Center
!~             vNew = matmul(invR,b)
            b = matmul( this% R,(v - this% CloudCenter))
            b(1:2) = b(1:2) - this% MBR% Center
            vNew = matmul(T,b)

         case('global')
         
!~             b = matmul(R,v)
!~             b(1:2) = b(1:2) + this% MBR% center
!~             vNew = this% CloudCenter + matmul(this% R,b)
            b = matmul(invT,v)
            b(1:2) = b(1:2) + this% MBR% center
            vNew = this% CloudCenter + matmul(this% invR,b)
            
      end select

   end subroutine ChangeRefFrame
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes a rotation Matrix assuming a starting orthonormal 
! base (1,0,0);(0,1,0);(0,0,1) and the final orthonormal one u, v, w
!  -------------------------------------------------   
   subroutine ComputeRotationMatrix( this, u, v, w )

      implicit none
      !-arguments------------------------------------------------
      class(OBB_type),                 intent(inout) :: this
      real(kind=rp),  dimension(NDIM), intent(in)    :: u, v, w
      
      this% R(1,:) = (/ u(1), u(2), u(3) /)
      this% R(2,:) = (/ v(1), v(2), v(3) /)
      this% R(3,:) = (/ w(1), w(2), w(3) /)
      
      this% invR(:,1) = this% R(1,:)
      this% invR(:,2) = this% R(2,:)
      this% invR(:,3) = this% R(3,:)

   end subroutine ComputeRotationMatrix
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the Oriented Bounding Box. 
!  -------------------------------------------------
   subroutine OBB_construct( this, stl, isPlot )
   
      implicit none
      !-arguments-------------------------------
      class(OBB_type), intent(inout) :: this
      type(STLfile),   intent(in)    :: stl
      logical,         intent(in)    :: isPlot
      !-local-variables-------------------------
      type(Hull_type) :: Hull
      real(kind=rp)   :: EigenVal
      integer         :: i
! 
!     Reading the data
!     ---------------
      call this% ReadStorePoints( stl )
!
!     Computing center of the points cloud
!     ---------------------------------
      this% CloudCenter = ComputeCentroid( this )

!
!     Diagonalization of the cloud
!     -------------------------
      call PointsCloudDiagonalization( this, EigenVal, this% MBR% t1, this% MBR% t2, this% MBR% normal ) 

!
!     Computing rotation matrix
!     ------------------------
      call this% ComputeRotationMatrix( this% MBR% t1, this% MBR% t2, this% MBR% normal ) 

!
!     Point projection
!     ---------------
      call ProjectPointsOnPlane( this )

!
!     Computing convex hull
!     ---------------------
      call ConvexHull( Hull, this )

!
!     Minimum Bounding Rectangle
!     ---------------------------
      call RotatingCalipers( Hull, this% MBR% Width, this% MBR% Length, this% MBR% Angle, this% MBR% Center )

!
!     Setting vertices of the MBR
!    -----------------------------
     call this% MBR% ComputeVertices()

!
!     Extrusion
!     ---------
      call ExtrudeMBR(  this )

      if( isPlot ) call this% plot()

      if(allocated(this% HullPoints)) deallocate(this% HullPoints)
      allocate( this% HullPoints(Hull% NumOfPoints) )

      do i = 1, Hull% NumOfPoints
         this% HullPoints(i) = Hull% Points(i)
      end do

      deallocate(this% Points, Hull% Points)

   end subroutine OBB_construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function computes the orientation of the rotation of 2D points. 
! Rotation := 1 -> Clockwise
! Rotation := 2 -> Counterclockwise
! Rotation := 0 -> Aligned
!  -------------------------------------------------
   
   integer function RotationType( Point1, Point2, Point3 ) result( Rotation )

      implicit none
      !-arguments-----------------------------------
      type(point_type), intent(in) :: Point1, Point2, Point3
      !-local-variables-------------------------------
      real(kind=rp), dimension(NDIM-1) :: b, c
      real(kind=rp)                    :: angle

      b = Point3% coords(1:2) - Point2% coords(1:2); c = Point2% coords(1:2) - Point1% coords(1:2)
      
      angle = b(1)*c(2) - b(2)*c(1)
      
      if( AlmostEqual(angle,0.0_RP) ) then
         Rotation = 0
      elseif( angle .gt. 0.0_RP ) then
         Rotation = 1
      else
         Rotation = 2
      end if
      
   end function RotationType
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------------------------
! This function returns the index of the left most point of a cloud. 
!  -------------------------------------------------------------------
   
   subroutine OBB_ComputeAngle( this, LowestIndex ) 

      implicit none
      !-arguments------------------------
      class(OBB_type), intent(inout) :: this
      integer,         intent(in)    :: LowestIndex
      !-local-variables--------------------
      real(kind=rp), dimension(NDIM) :: v
      integer                        :: i    

!$omp parallel shared(this, LowestIndex,i)
!$omp do schedule(runtime) private(v)
      do i = 1, this% NumOfPoints
         if( i .eq. LowestIndex ) then
            this% Points(i)% theta = -2.0_RP
            cycle
         end if
         v = this% Points(i)% coords-this% Points(LowestIndex)% coords
         if( almostEqual(norm2(v(1:2)),0.0_RP)  ) then 
            this% Points(i)% theta = -1.0_RP 
            cycle
         end if
         this% Points(i)% theta = acos(v(1)/norm2(v(1:2)))
      end do
!$omp end do
!$omp end parallel

   end  subroutine OBB_ComputeAngle
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ------------------------------------------------------------------
! This function returns the index of the left most point of a cloud. 
!  ------------------------------------------------------------------

   subroutine ComputeHullExtremePoint( OBB, LowestIndex ) 

      implicit none
      !-arguments------------------------
      class(OBB_type),  intent(inout) :: OBB
      integer,          intent(out)   :: LowestIndex
      !-local-variables--------------------
      type(point_type) :: LowestPoint
      integer          :: i
      
      LowestPoint = OBB% Points(1)
      
      do i = 2, OBB% NumOfPoints
         if( OBB% Points(i)% coords(2) .lt. LowestPoint% coords(2) ) then
            LowestPoint = OBB% Points(i)
         elseif( AlmostEqual( OBB% Points(i)% coords(2), LowestPoint% coords(2)) ) then
            if( OBB% Points(i)% coords(1) .lt. LowestPoint% coords(1) ) then
               LowestPoint = OBB% Points(i)
            end if
         end if
      end do
       
      LowestIndex = LowestPoint% index
      
   end subroutine ComputeHullExtremePoint
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine orders the nodes according to the angle theta
!  ------------------------------------------------    
   
   subroutine OBB_SortingNodes( this, left, right )
      use StopwatchClass 
      implicit none
      !-arguments--------------------------------------------------------------
      class(OBB_type), intent(inout) :: this
      integer,         intent(in)    :: left, right
      !-local-variables--------------------------------------------------------
      integer :: i
      
!$omp parallel shared(this,left, right)
!$omp single
      call sort( this% Points(:)% theta, this% Points(:)% index, this% Points(:)% coords(1), &
                 this% Points(:)% coords(2), this% Points(:)% coords(3), left, right )
!$omp end single
!$omp end parallel
      
      do i = 1, this% NumOfPoints-1
         if( almostEqual(this% Points(i+1)% theta, this% Points(i)% theta) ) this% Points(i+1)% delete = .true.
      end do

   end subroutine OBB_SortingNodes
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine sort the elements of a from the lowest to the grates. b saves the old indeces. 
!  ------------------------------------------------    
   recursive subroutine sort( a, b, coordx, coordy, coordz, left, right )
   
      implicit none
      !-arguments-------------------------------------------------
      real(kind=rp), dimension(:),   intent(inout) :: a
      integer,       dimension(:),   intent(inout) :: b
      real(kind=rp), dimension(:),   intent(inout) :: coordx, coordy, coordz
      integer,                       intent(in)    :: left, right
      !-local-variables-------------------------------------------
      real(kind=rp) :: x, t
      integer       :: i, j, ind, left_cycle, &
                       right_cycle
                       
      if( left .ge. right ) return
   
      x = a( (left + right)/2 )
    
      i = left; j = right

      do
         do while( a(i) .lt. x ) 
            i = i+1
         end do
         do while( x .lt. a(j) )
            j = j-1
         end do
         if( i .ge. j ) exit
         t   = a(i); a(i) = a(j); a(j) = t
         ind = b(i); b(i) = b(j); b(j) = ind
         t   = coordx(i); coordx(i) = coordx(j); coordx(j) = t
         t   = coordy(i); coordy(i) = coordy(j); coordy(j) = t
         t   = coordz(i); coordz(i) = coordz(j); coordz(j) = t
         
         i = i+1
         j = j-1
      end do
      
      left_cycle = (i-1)-left
      right_cycle = right-(j+1)
      
!$omp task shared(a,b,coordx,coordy,coordz,left,i) if(left_cycle > TASK_THRESHOLD)
         call sort( a, b, coordx, coordy, coordz, left, i-1 )
!$omp end task

!$omp task shared(a,b,coordx,coordy,coordz,j,right) if(right_cycle > TASK_THRESHOLD)
      call sort( a, b, coordx, coordy, coordz, j+1, right )
!$omp end task
!$omp taskwait

   end subroutine sort
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function computes the Centroid of the points cloud. 
!  -------------------------------------------------
   
   function ComputeCentroid( OBB ) result( CloudCenter )

      implicit none
      !-arguments-----------------------------------
      type(OBB_type),  intent(in) :: OBB
      !-local-variables-------------------------------
      real(kind=rp), dimension(NDIM) :: CloudCenter
      integer                        :: i

      CloudCenter = 0.0_RP

      do i = 1, OBB% NumOfPoints
         CloudCenter = CloudCenter + OBB% Points(i)% coords 
      end do

      CloudCenter = CloudCenter/OBB% NumOfPoints 
      
   end  function ComputeCentroid
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the covariance matrix of a clouds of points given its center. 
!  -------------------------------------------------
   
   subroutine CovarianceMatrix( OBB, CovMat )

      implicit none
      !-arguments----------------------------------------------
     type(OBB_type),                       intent(in)   :: OBB
      real(kind=rp), dimension(NDIM,NDIM), intent(out) :: CovMat
      !-local-variables----------------------------------------------
      real(kind=rp), dimension(NDIM) :: d
      integer                        :: i, j, k

      CovMat = 0.0_RP

      do i = 1, OBB% NumOfPoints
         d = OBB% Points(i)% coords - OBB% CloudCenter
         do k = 1, NDIM; do j = 1,  NDIM
            CovMat(k,j) = CovMat(k,j) + d(k) * d(j) 
         end do; end do
      end do
      
      CovMat = CovMat/OBB% NumOfPoints
          
   end  subroutine CovarianceMatrix
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the projection of the points cloud on a plane defined by a normal vector. 
!  -------------------------------------------------
   
   subroutine ProjectPointsOnPlane( OBB  )
      use MappedGeometryClass
      implicit none
      !-arguments----------------------------------------------------------
      type(OBB_type), intent(inout) :: OBB
      !-local-variables------------------------------------------------------
      real(kind=rp), dimension(NDIM) :: d
      real(kind=rp)                  :: N_point
      integer :: i 
!
!     Initialization of nMin and nMax, i.e. of the min and max normal distance of the cloud w.r.t 
!     the projection plane 
!     -------------------------------------------------------------------------------
      OBB% nMax = -huge(1.0_RP); OBB% nMin = huge(1.0_RP)

      do i = 1, OBB% NumOfPoints
         d = OBB% Points(i)% coords - OBB% CloudCenter
         N_Point = vdot( d, OBB% MBR% normal )   
         OBB% Points(i)% coords = (/ vdot(d, OBB% MBR% t1), vdot(d, OBB% MBR% t2), 0.0_RP /)   
         OBB% nMax = max(OBB% nMax,N_point)
         OBB% nMin = min(OBB% nMin,N_point)
      end do

      OBB% nMax = (1.0_RP+SAFETY_FACTOR)*OBB% nMax
      OBB% nMin = (1.0_RP+SAFETY_FACTOR)*OBB% nMin 

   end  subroutine ProjectPointsOnPlane
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the eigenstructure of the OBB. 
!  -------------------------------------------------
   
   subroutine PointsCloudDiagonalization( OBB, EigenVal, EigenVec1, EigenVec2, EigenVec3  ) 

      implicit none
      !-arguments---------------------------------------------
      type(OBB_type),                  intent(in)  :: OBB
      real(kind=rp),                   intent(out) :: EigenVal
      real(kind=rp), dimension(NDIM), intent(out) :: EigenVec1, EigenVec2, EigenVec3    
      !-local-variables----------------------------
      real(kind=rp), dimension(NDIM,NDIM) :: CovMat

      call CovarianceMatrix( OBB, CovMat )      

      call ComputeEigenStructure( CovMat, EigenVal, EigenVec1, EigenVec2, EigenVec3 )
   
   end subroutine PointsCloudDiagonalization
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the Convex Hull of a set of points. 
! The points are ordered according to the value of the x-coord.
!  -------------------------------------------------

   subroutine ConvexHull( Hull, OBB )

      implicit none
      !-arguments------------------------------------------------
      type(Hull_type), intent(inout) :: Hull
      type(OBB_type),  intent(inout) :: OBB
      !-local-variables--------------------------------------------
      type(point_type) ,pointer :: p, p1
      integer                   :: i, LowestIndex, start
      type(PointLinkedList)     :: convPoints
!     
!     Compute Left most point & set Hull's head
!     -----------------------------------------
      call ComputeHullExtremePoint( OBB, LowestIndex )

!     
!     Compute the angles
!     ------------------
      call OBB% ComputeAngle( LowestIndex ) 

!
!     Sorting points according to their angle  
!     ---------------------------------------
      call OBB% SortingNodes(1, OBB% NumOfPoints )  

      convPoints = PointLinkedList()
      call convPoints% add( OBB% Points(1) )
      
      !Check duplicate
      if( almostEqual(OBB% Points(1)% coords(1), OBB% Points(2)% coords(1)) .and. &
          almostEqual(OBB% Points(1)% coords(2), OBB% Points(2)% coords(2)) ) OBB% Points(2)% delete = .true.

      do i = 2, OBB% NumOfPoints     
         if( .not. OBB% Points(i)% delete ) then
            call convPoints% add( OBB% Points(i) )
            start = i+1
            exit
         end if
      end do   
        
      Hull% NumOfPoints = 2

      p  => convPoints% head
      p1 => convPoints% head% next

      do i = start, OBB% NumOfPoints
         if( OBB% Points(i)% delete ) cycle
         do while( RotationType(p, p1, OBB% Points(i) ) .eq. 1  )
            p => p% prev; p1 => p1% prev
            call convPoints% RemoveLast()
            Hull% NumOfPoints = Hull% NumOfPoints - 1
         end do
         call convPoints% add( OBB% Points(i) )
         p1 => convPoints% head% prev; p => p1% prev
         Hull% NumOfPoints = Hull% NumOfPoints + 1 
      end do      

      allocate(Hull% Points(Hull% NumOfPoints))
      
      p => convPoints% head

      do i = 1, Hull% NumOfPoints
         Hull% Points(i) = p
         p => p% next
      end do
         
      call convPoints% destruct()

   end subroutine ConvexHull
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the 2-pi normalized-angle between 2 vectors. 
!  -------------------------------------------------
   real(kind=rp) function ComputingAngle( a, b ) result( Angle )
   
      implicit none
      !-arguments------------------------------------
      real(kind=rp), dimension(:), intent(in) :: a, b
      !-local-variables------------------------------ 
      real(kind=rp) :: theta1, theta2
   
      theta1 = atan2( a(2), a(1) )
      theta2 = atan2( b(2), b(1) )
      
      Angle = theta2 - theta1
   
      Angle = mod(Angle+2.0_RP*PI, 2.0_RP*PI) 
   
   end function ComputingAngle
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine counterclockwise-rotates 2D vector. 
!  -------------------------------------------------   
   subroutine RotateVector( RotVec, theta )

      implicit none
      !-arguments------------------------------------------
      real(kind=rp), dimension(:), intent(inout) :: RotVec
      real(kind=rp),               intent(in)    :: theta
      !-local-variables------------------------------------
      real(kind=rp), dimension(size(RotVec)) :: TempVec
!
!     Temporary vector
!     ----------------
      TempVec = RotVec  
      
      RotVec(1) = TempVec(1)*cos(theta) - TempVec(2)*sin(theta)
      RotVec(2) = TempVec(1)*sin(theta) + TempVec(2)*cos(theta)
      
   end subroutine RotateVector
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine the parameter t_s use for defining the position of the point S on a line.
! x_s = x_p + t_s(x_q-x_p)
! y_s = y_p + t_s(y_q-y_p)
!  -------------------------------------------------
   function CurvAbscissa( p1, p2, p3 ) result( t )
   
      implicit none
      !-arguments--------------------------------
      type(point_type), intent(in) :: p1, p2, p3
      !-local-variables--------------------------
      real(kind=rp) :: t

     t = ( p3% coords(1) - p1% coords(1) )*( p2% coords(1) - p1% coords(1) ) + &
         ( p3% coords(2) - p1% coords(2) )*( p2% coords(2) - p1% coords(2) )

     t = t/(  POW2(p2% coords(1) - p1% coords(1)) + POW2(p2% coords(2) - p1% coords(2))  )

   end function CurvAbscissa
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the distance of a point (p3) from a line defined by
! the points p1 and p2 
!  -------------------------------------------------
    real(kind=rp) function ComputeWidth( p1, p2, p3 ) result( width )
    
       implicit none
       !-arguments--------------------------------
       type(point_type), intent(in) :: p1, p2, p3
       !-local-variables--------------------------
       real(kind=rp), dimension(NDIM-1) :: vec
       real(kind=rp)                     :: t
!
!     Computing t-values for the point p3
!     --------------------------------
       t = CurvAbscissa( p1, p2, p3 )
       
!
!     Position of the projection of p3 point on the line p1-p2
!     ------------------------------------------------
       vec = p1% coords(1:2) + t *( p2% coords(1:2)- p1% coords(1:2) )
       
!
!     Computing the distance
!     ---------------------
       width = norm2( p3% coords(1:2) - vec(1:2) )
 
   end function ComputeWidth 
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the rotating calipers method. 
!  -------------------------------------------------
   subroutine RotatingCalipers( Hull, rectWidth, rectLength, rectAngle, rectCenter )

      implicit none
      !-arguments-----------------------------------------------------------------
      type(Hull_type),             intent(inout) :: Hull
      real(kind=rp), dimension(:), intent(out)   :: rectCenter
      real(kind=rp),               intent(out)   :: rectWidth, rectLength, rectAngle
      !-local-variables-------------------------------------------------------------
      real(kind=rp), dimension(NDIM-1)            :: Caliper1, Caliper2, v1, v2,  PointMin, PointMax
      integer,       dimension(1)                 :: loc, MinIndex, MaxIndex
      real(kind=rp)                               :: RotAngle, width, minwidth, theta1, theta2, dtheta
      real(kind=rp), dimension(Hull% NumOfPoints) :: x, y
      integer                                     :: indexMin1, indexMin2, indexMax1, indexMax2
!
!     Initialization
!     ------------
      minwidth = huge(1.0_RP); rectAngle = 0.0_RP; RotAngle = 0.0_RP
! 
!     Initializing caliper
!     ----------------
      Caliper1 = (/ 1.0_RP, 0.0_RP/)
      Caliper2 = (/ -1.0_RP, 0.0_RP/)
      
      MinIndex = minloc(Hull% Points(:)% coords(2))
      MaxIndex = maxloc(Hull% Points(:)% coords(2))
      
      indexMin1 = MinIndex(1); indexMax1 = MaxIndex(1)   
      
      do while( RotAngle .lt. PI )
      
         indexMin2 = mod(indexMin1, Hull% NumOfPoints) + 1
         indexMax2 = mod(indexMax1, Hull% NumOfPoints) + 1
      
!
!        Looping Counterclockwise from y-min point (p1) 
!        and from y-max point (p2). 
!        -------------------------------------------
         v1 = Hull% Points(indexMin2)% coords(1:2) - Hull% Points(indexMin1)% coords(1:2) 
         v2 = Hull% Points(indexMax2)% coords(1:2) - Hull% Points(indexMax1)% coords(1:2)
         
!
!       Computing the angles between calipers and vectors v1,v2
!       ---------------------------------------------------    
         theta1 = ComputingAngle( Caliper1, v1 ) 
         theta2 = ComputingAngle( Caliper2, v2 )
      
         dtheta = min(theta1,theta2)
         loc    = minloc((/ theta1,theta2 /))
         call RotateVector(Caliper1, dtheta)
         call RotateVector(Caliper2, dtheta)
         
         RotAngle = RotAngle + dtheta
         
         if( loc(1) .eq. 1 ) then 
            width = ComputeWidth( Hull% Points(indexMin1), Hull% Points(indexMin2), Hull% Points(indexMax1) )
            indexMin1 = mod(indexMin1, Hull% NumOfPoints) + 1
         else
            width = ComputeWidth( Hull% Points(indexMax1), Hull% Points(indexMax2), Hull% Points(indexMin1) )
            indexMax1 = mod(indexMax1, Hull% NumOfPoints) + 1
         end if
      
         if( width .lt. minwidth ) then
            minwidth  = width
            rectAngle = RotAngle
         end if
         
      end do 
      
      rectCenter(1) = sum(Hull% Points(:)% coords(1))/Hull% NumOfPoints
      rectCenter(2) = sum(Hull% Points(:)% coords(2))/Hull% NumOfPoints
      
!
!     Initialization of the extreme points
!     ---------------------------------------
      PointMax = -huge(1.0_RP); PointMin = huge(1.0_RP)
!
!     ClockWise rotation looking for extreme points
!     ---------------------------------------------
    
      x = (Hull% Points(:)% coords(1)-rectCenter(1))*cos(rectAngle) + &
          (Hull% Points(:)% coords(2)-rectCenter(2))*sin(rectAngle)
      y = -(Hull% Points(:)% coords(1)-rectCenter(1))*sin(rectAngle) + &
           (Hull% Points(:)% coords(2)-rectCenter(2))*cos(rectAngle)
    
      PointMax(1) = maxval(x)
      PointMax(2) = maxval(y)
      PointMin(1) = minval(x)
      PointMin(2) = minval(y)
      
!
!     Minimum Bounding Rectangle properties
!     -------------------------------------
      rectLength = PointMax(1) - PointMin(1);
      rectWidth  = PointMax(2) - PointMin(2);
      
      !Safety factor
      rectLength = (1.0_RP+SAFETY_FACTOR)*rectLength
      rectWidth  = (1.0_RP+SAFETY_FACTOR)*rectWidth
      
      rectCenter = 0.5_RP* (/ PointMax(1) + PointMin(1),PointMax(2) + PointMin(2) /)

      call RotateVector(rectCenter,rectAngle)
      
      rectCenter(1) = rectCenter(1) + sum(Hull% Points(:)% coords(1))/Hull% NumOfPoints
      rectCenter(2) = rectCenter(2) + sum(Hull% Points(:)% coords(2))/Hull% NumOfPoints
      
   end subroutine RotatingCalipers
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine extrude the minimum bounding rectangle to get the full OBB.
! MBR -> Minimum Bounding Rectangle
! MRB := ( rectWidth, rectLength, rectAngle, rectCenter )
!  -------------------------------------------------
   subroutine ExtrudeMBR( OBB )
   
      implicit none
      !-arguments-------------------------------------------------
      type(OBB_type), intent(inout) :: OBB
      !-local-variables-------------------------------------------
      integer                       :: i
!
!     Extrusion
!     ---------

      do i = 1, 4 
         call OBB% ChangeRefFrame((/ OBB% MBR% vertices(1,i), OBB% MBR% vertices(2,i), OBB% nMin /), 'global', OBB% vertices(:,i) )
         call OBB% ChangeRefFrame((/ OBB% MBR% vertices(1,i), OBB% MBR% vertices(2,i), OBB% nMax /), 'global', OBB% vertices(:,i+4) )
         OBB% LocVertices(:,i)   = (/ OBB% MBR% vertices(:,i), OBB% nMin /)
         OBB% LocVertices(:,i+4) = (/ OBB% MBR% vertices(:,i), OBB% nMax /)
      end do

      OBB% LocFrameCenter = OBB% CloudCenter + matmul(OBB% R, (/ OBB% MBR% Center,0.0_RP/))
      
   end subroutine ExtrudeMBR   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function check if a point is inside the OBB or not
!  -------------------------------------------------   
   function OBB_isPointInside( this, coords, coeff ) result( isInsideOBB )
   
      implicit none
      !-arguments------------------------------------
      class(OBB_type),             intent(inout) :: this
      real(kind=rp), dimension(:), intent(in)    :: coords
      real(kind=rp),               intent(in)    :: coeff
      logical                                    :: isInsideOBB
      !-local-variables--------------------------------
      real(kind=rp), dimension(NDIM) :: Point
      real(kind=rp)                  :: Multcoeff
      
      optional :: coeff
      
      if( present(coeff) ) then
         Multcoeff = coeff
      else
         Multcoeff = 1.0_RP
      end if
      
      call this% ChangeRefFrame(coords, 'local', Point)
      
      isInsideOBB = .false.
      
      !check x y z     
       if( isInsideBox(Point, Multcoeff*this% LocVertices) ) isInsideOBB = .true.

   end function OBB_isPointInside
   
   
   subroutine OBB_ChangeObjsRefFrame( this, objs )
   
      implicit none
      !-arguments---------------------------------------------
      class(OBB_type),                 intent(inout) :: this
      type(Object_type), dimension(:), intent(inout) :: objs
      !-local-variables---------------------------------------
      integer :: i, j
!$omp parallel shared(this,objs,i)
!$omp do schedule(runtime) private(j)
      do i = 1, size(objs)
         do j = 1, size(objs(i)% vertices)
            call this% ChangeRefFrame( objs(i)% vertices(j)% coords, 'local', objs(i)% vertices(j)% coords ) 
         end do
      end do
!$omp end do
!$omp end parallel

   end subroutine OBB_ChangeObjsRefFrame
   
      
   subroutine OBB_STL_rotate( this, stl )
   
      implicit none
      !-arguments-----------------------------
      class(OBB_type), intent(inout):: this
      type(STLfile),   intent(inout):: stl
      !-local-variables-----------------------
      real(kind=RP) :: meanCoords(NDIM), meanCoordsNew(NDIM), shiftVec(NDIM)
      integer :: i, j, NumOfPoints
      
      meanCoords = 0.0_RP; meanCoordsNew = 0.0_RP; NumOfPoints = 0
      
!$omp parallel shared(this,stl,meanCoords,meanCoordsNew,NumOfPoints,i)
!$omp do schedule(runtime) private(j)
      do i = 1, stl% NumOfObjs
         do j = 1, size(stl% ObjectsList(i)% vertices)
            call this% ChangeRefFrame( stl% ObjectsList(i)% vertices(j)% coords, 'global', &
                                       stl% ObjectsList(i)% vertices(j)% coords             ) 
!$omp critical
            meanCoords = meanCoords + stl% ObjectsList(i)% vertices(j)% coords 
            NumOfPoints = NumOfPoints + 1
!$omp end critical
            stl% ObjectsList(i)% vertices(j)% coords = matmul( stl% rotationMatrix, stl% ObjectsList(i)% vertices(j)% coords )   
!$omp critical
            meanCoordsNew = meanCoordsNew + stl% ObjectsList(i)% vertices(j)% coords
!$omp end critical
         end do
      end do
!$omp end do 
!$omp end parallel   
      
      meanCoords = meanCoords/NumOfPoints; meanCoordsNew = meanCoordsNew/NumOfPoints
      
      shiftVec = meanCoordsNew - meanCoords
!$omp parallel shared(i,shiftVec)
!$omp do schedule(runtime) private(j)
      do i = 1, stl% NumOfObjs
         do j = 1, size(stl% ObjectsList(i)% vertices)
            stl% ObjectsList(i)% vertices(j)% coords = stl% ObjectsList(i)% vertices(j)% coords - shiftVec
         end do
      end do
!$omp end do
!$omp end parallel
   end subroutine OBB_STL_rotate
   
   subroutine OBB_STL_translate( this, stl )
   
      implicit none
      !-arguments-----------------------------
      class(OBB_type), intent(inout):: this
      type(STLfile),   intent(inout):: stl
      !-local-variables-----------------------
      real(kind=RP) :: meanCoords(NDIM)
      integer :: i, j, NumOfPoints
      
      meanCoords = 0.0_RP; NumOfPoints = 0
      
!$omp parallel shared(this,stl,meanCoords,NumOfPoints,i)
!$omp do schedule(runtime) private(j)
      do i = 1, stl% NumOfObjs
         do j = 1, size(stl% ObjectsList(i)% vertices)
            call this% ChangeRefFrame( stl% ObjectsList(i)% vertices(j)% coords, 'global', &
                                       stl% ObjectsList(i)% vertices(j)% coords             ) 
            stl% ObjectsList(i)% vertices(j)% coords(stl% motionAxis) = stl% ObjectsList(i)% vertices(j)% coords(stl% motionAxis) + stl% ds
!$omp critical
            meanCoords = meanCoords + stl% ObjectsList(i)% vertices(j)% coords 
            NumOfPoints = NumOfPoints + 1
!$omp end critical
         end do
      end do
!$omp end do 
!$omp end parallel   

   end subroutine OBB_STL_translate

!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function check if a point is inside a generic box
!  -------------------------------------------------   

   logical function isInsideBox( Point, vertices, equal, coeff ) result( isInside )
   
      implicit none

      real(kind=rp), dimension(:),   intent(in) :: Point
      real(kind=rp), dimension(:,:), intent(in) :: vertices
      logical,                       intent(in) :: equal
      real(kind=rp), optional,       intent(in) :: coeff
      
      optional :: equal

      isInside = .false.

      if( present(equal) ) then
         if( .not. equal ) then
            if( (Point(1) > vertices(1,1) .and. Point(1) < vertices(1,7)) .and. &
                (Point(2) > vertices(2,1) .and. Point(2) < vertices(2,7)) .and. &
                (Point(3) > vertices(3,1) .and. Point(3) < vertices(3,7)) ) then
                isInside = .true.
            end if
            return
         end if
      end if

      if( (Point(1) >= vertices(1,1) .and. Point(1) <= vertices(1,7)) .and. &
          (Point(2) >= vertices(2,1) .and. Point(2) <= vertices(2,7)) .and. &
          (Point(3) >= vertices(3,1) .and. Point(3) <= vertices(3,7)) ) isInside = .true.

   end function isInsideBox
   
   
   
   logical function OBB_isInsidePolygon( OBB, Point, coeff ) result( isInside )
      use MappedGeometryClass
      implicit none
      
      class(OBB_type),               intent(inout) :: OBB
      real(kind=rp),   dimension(:), intent(in)    :: Point
      real(kind=rp),                 intent(in)    :: coeff
      !-local-variables------------------------------------------
      real(kind=rp) :: d(NDIM), coords(NDIM), v1(NDIM-1), v2(NDIM-1), &
                       v3(NDIM-1), v4(NDIM-1), length, N_Point
      logical :: Intersection
      integer :: NumOfIntersections, i
      
      d       = Point - OBB% CloudCenter
      N_Point = vdot( d, OBB% MBR% normal ) 
      coords  = (/ vdot(d, OBB% MBR% t1), vdot(d, OBB% MBR% t2), 0.0_RP /) 
 
      length = maxval(abs(OBB% HullPoints(:)% coords(1))) 
      
      v1 = coords(1:2)
      v2 = coords(1:2); v2(1) = v2(1) + 1000_RP*length
      
      NumOfIntersections = 0

      do i = 1, size(OBB% HullPoints)-1
         
         v3 = OBB% HullPoints(i+1)% coords(1:NDIM-1)
         v4 = OBB% HullPoints(i)% coords(1:NDIM-1)
      
         Intersection = TwoD_RayTracing( v1, v2, coeff*v3, coeff*v4 ) 
      
         if( Intersection ) NumOfIntersections = NumOfIntersections + 1 
      
      end do      
 
      if( mod(NumOfIntersections,2) .eq. 0 ) then
         isInside = .false.
      else
         isInside = .true.
      end if
      
      if( N_Point .gt. coeff*OBB% nMax .or. N_Point .lt. coeff*OBB% nMin ) then
         isInside = .false.
      end if 
      
   end function OBB_isInsidePolygon
   
   
   logical function TwoD_RayTracing( v1, v2, v3, v4 ) result( Intersection )
   
      implicit none
      !-arguments--------------------------------------------------------
      real(kind=rp), dimension(NDIM-1), intent(in) :: v1, v2, v3, v4
      !-local-varirables------------------------------------------------
      real(kind=rp) :: a1, a2, b1, b2, c1, c2, d1, d2
      
      Intersection = .false.
      
      a1 = v2(2) - v1(2)
      b1 = v1(1) - v2(1)
      c1 = (v2(1)*v1(2)) - (v1(1)*v2(2))
      
      d1 = a1*v3(1) + b1*v3(2) + c1
      d2 = a1*v4(1) + b1*v4(2) + c1
   
      if( d1 * d2 > 0.0_RP ) then
         Intersection = .false.
         return
      end if

      a2 = v4(2) - v3(2)
      b2 = v3(1) - v4(1)
      c2 = (v4(1)*v3(2)) - (v3(1)*v4(2))

      d1 = a2*v1(1) + b2*v1(2) + c2
      d2 = a2*v2(1) + b2*v2(2) + c2
      
      if( d1 * d2 > 0.0_RP ) then
         Intersection = .false.
         return
      end if
      
      if( AlmostEqual( ((a1*b2) - (a2*b1)) ,0.0_RP) )then
         Intersection = .false.
         return
      end if
      
      Intersection = .true.
      
   end function TwoD_RayTracing

end module OrientedBoundingBox
