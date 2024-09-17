#include "Includes.h"
module MPI_IBMUtilities

   use SMConstants
   use Utilities
   use TessellationTypes
   use OrientedBoundingBox
   use MPI_Process_info
   use FaceClass
#ifdef _HAS_MPI_
   use mpi
#endif
   
   implicit none
   
   private
   public :: IBMpoints
   public :: SendSTL2Partitions, receiveSTLpartitions
   public :: recvScalarPlotRoot, sendScalarPlotRoot
   public :: recvVectorPlotRoot, sendVectorPlotRoot
   public :: MPIProcedures_IBM_HO_faces, IBM_HO_findElements
   public :: Set_IBM_HO_faces, GatherHOfacesState
   public :: plotSTL, SendOBB, RecvOBB, SendAxis, RecvAxis
   public :: castMaskNumOfObjs, castMask, castMaskPlot, castIsInsideBody
   public :: gatherMaskGeom
   public :: castStateBandRegion, castGradientsBandRegion
   public :: IBM_HO_GetState
   public :: FixingmpiFaces
   public :: MPIProcedures_IBM_HOIntegrationPoints
   public :: IBM_HO_GetGradient, GatherHOIntegrationPointsState
   public :: Mask2Root, Mask2Partitions, MaskLogical2Root
   public :: StateMask2Root, StateMask2Partitions

   type IBMpoints

      type(point_type), allocatable :: x(:)
      real(kind=RP),    allocatable :: coords(:,:), Q(:,:), U_x(:,:), U_y(:,:), U_z(:,:), dist(:), normal(:,:), xi(:,:), V(:,:), coordsNEW(:,:), lj(:)
      integer,          allocatable :: element_index(:), local_position(:,:), NumOfIntersections(:), domain(:), N(:), fIDs(:), STLNum(:)
      logical,          allocatable :: isInsideBody(:)
      integer                       :: LocNumOfObjs, NumOfObjs
      logical                       :: computeV = .true. 

      contains 
         procedure :: build           => IBMpoints_build 
         procedure :: buildPoints     => IBMpoints_buildPoints 
         procedure :: buildstencil    => IBMpoints_buildstencil 
         procedure :: copy            => IBMpoints_copy 
         procedure :: copyPoints      => IBMpoints_copyPoints 
         procedure :: destroy         => IBMpoints_destroy 
         procedure :: destroyPoints   => IBMpoints_destroyPoints 
         procedure :: buildState      => IBMpoints_buildState
         procedure :: buildBandRegion => IBMpoints_buildBandRegion
   end type

contains

   subroutine IBMpoints_build( this, NumOfObjs )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 
      integer,           intent(in)    :: NumOfObjs
      
      allocate( this% coords(NumOfObjs,NDIM),         &
                this% element_index(NumOfObjs),       &
                this% local_position(NumOfObjs,NDIM), &
                this% NumOfIntersections(NumOfObjs),  & 
                this% isInsideBody(NumOfObjs),        & 
                this% STLNum(NumOfObjs)               )

      this% NumOfObjs          = NumOfObjs 
      this% NumOfIntersections = 0
      this% isInsideBody       = .false.
      this% STLNum             = 0

   end subroutine IBMpoints_build

   subroutine IBMpoints_buildstencil( this, NumOfObjs )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 
      integer,           intent(in)    :: NumOfObjs
      
      allocate( this% coords(NumOfObjs,NDIM),         &
                this% element_index(NumOfObjs),       &
                this% local_position(NumOfObjs,NDIM), &
                this% xi(NumOfObjs,NDIM),             &
                this% N(NumOfObjs),                   & 
                this% fIDs(NumOfObjs),                & 
                this% domain(NumOfObjs),              & 
                this% STLNum(NumOfObjs)               ) 

      this% NumOfObjs     = NumOfObjs 
      this% N             = 0 
      this% domain        = 0 
      this% element_index = 0
      this% fIDs          = 0
      this% STLNum        = 0

   end subroutine IBMpoints_buildstencil 

   subroutine IBMpoints_buildBandRegion( this, NumOfObjs )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 
      integer,           intent(in)    :: NumOfObjs
      
      allocate( this% xi(NumOfObjs,NDIM), &
                this% domain(NumOfObjs),  &
                this% N(NumOfObjs)        )
      
      this% N      = 0 
      this% domain = 0

   end subroutine IBMpoints_buildBandRegion

   subroutine IBMpoints_buildState( this, NumOfObjs, nEqn )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 
      integer,           intent(in)    :: NumOfObjs, nEqn

      allocate( this% Q(NumOfObjs,nEqn)  , &
                this% U_x(NumOfObjs,nEqn), &
                this% U_y(NumOfObjs,nEqn), &
                this% U_z(NumOfObjs,nEqn)  )

      this% Q   = 0.0_RP
      this% U_x = 0.0_RP
      this% U_y = 0.0_RP
      this% U_z = 0.0_RP

   end subroutine IBMpoints_buildState

   subroutine IBMpoints_copy( this, tocopy )

      implicit none 

      class(IBMPoints), intent(inout) :: this
      type(IBMPoints),  intent(in)    :: tocopy 

      integer :: i 

      do i = 1, tocopy% NumOfObjs
         this% coords(i,:)         = tocopy% coords(i,:)
         this% element_index(i)    = tocopy% element_index(i)
         this% local_position(i,:) = tocopy% local_position(i,:)
         this% STLNum(i)           = tocopy% STLNum(i)
      end do

   end subroutine IBMpoints_copy

   subroutine IBMpoints_buildPoints( this, NumOfObjs )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 
      integer,           intent(in)    :: NumOfObjs
      
      allocate( this% x(NumOfObjs) )
      this% NumOfObjs = NumOfObjs
      
   end subroutine IBMpoints_buildPoints

   subroutine IBMpoints_copyPoints( this, tocopy )

      implicit none 

      class(IBMPoints), intent(inout) :: this
      type(IBMPoints),  intent(in)    :: tocopy 

      integer :: i 

      do i = 1, tocopy% NumOfObjs
         this% x(i)% coords         = tocopy% x(i)% coords
         this% x(i)% element_index  = tocopy% x(i)% element_index
         this% x(i)% local_position = tocopy% x(i)% local_position
         this% x(i)% STLNum         = tocopy% x(i)% STLNum
      end do 

   end subroutine IBMpoints_copyPoints

   subroutine IBMpoints_destroy( this )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 

      if( allocated(this% coords            ) ) deallocate( this% coords             )
      if( allocated(this% element_index     ) ) deallocate( this% element_index      )
      if( allocated(this% local_position    ) ) deallocate( this% local_position     )
      if( allocated(this% NumOfIntersections) ) deallocate( this% NumOfIntersections )
      if( allocated(this% isInsideBody      ) ) deallocate( this% isInsideBody       )
      if( allocated(this% STLNum            ) ) deallocate( this% STLNum             )   
      if( allocated(this% Q                 ) ) deallocate( this% Q                  )
      if( allocated(this% U_x               ) ) deallocate( this% U_x                )
      if( allocated(this% U_y               ) ) deallocate( this% U_y                )
      if( allocated(this% U_z               ) ) deallocate( this% U_z                )

      this% NumOfObjs = 0

   end subroutine IBMpoints_destroy

   subroutine IBMpoints_destroyPoints( this )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 

      if( allocated(this% x) ) deallocate( this% x )

      this% NumOfObjs = 0

   end subroutine IBMpoints_destroyPoints

   subroutine StateMask2Root( IBMmask, domain, nEqn )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain, nEqn

      if( .not. MPI_Process% isRoot ) then 
         call sendStateMask2Root( IBMmask, domain, nEqn )
      end if 
      
      if( MPI_Process% isRoot ) then 
         call recvStateMaskRoot( IBMmask, domain, nEqn )
      end if 
      
   end subroutine StateMask2Root

   subroutine StateMask2Partitions( IBMmask, domain, nEqn )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain, nEqn 

      integer :: domains 
     
      if( MPI_Process% isRoot ) then 
         call sendStateMask2partitions( IBMmask, domain, nEqn )
      end if

      if( .not. MPI_Process% isRoot ) then 
         call recvStateMaskPartitions( IBMmask, domain, nEqn )
      end if

   end subroutine StateMask2Partitions

   subroutine sendStateMask2Root( IBMmask, domain, nEqn )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain, nEqn
#ifdef _HAS_MPI_
      integer   :: i, sendQ_req(nEqn), sendUx_req(nEqn), sendUy_req(nEqn), sendUz_req(nEqn), ierr

      do i = 1, nEqn
         call mpi_isend( IBMmask(domain)% Q(:,i)  , IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, sendQ_req(i) , ierr ) 
         call mpi_isend( IBMmask(domain)% U_x(:,i), IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, sendUx_req(i), ierr ) 
         call mpi_isend( IBMmask(domain)% U_y(:,i), IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, sendUy_req(i), ierr ) 
         call mpi_isend( IBMmask(domain)% U_z(:,i), IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, sendUz_req(i), ierr ) 
      end do 

      call mpi_waitall( nEqn, sendQ_req , MPI_STATUSES_IGNORE, ierr )
      call mpi_waitall( nEqn, sendUx_req, MPI_STATUSES_IGNORE, ierr )
      call mpi_waitall( nEqn, sendUy_req, MPI_STATUSES_IGNORE, ierr )
      call mpi_waitall( nEqn, sendUz_req, MPI_STATUSES_IGNORE, ierr )
#endif
   end subroutine sendStateMask2Root

   subroutine recvStateMaskRoot(IBMmask, domain, nEqn )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain, nEqn
#ifdef _HAS_MPI_
      integer :: domains, i, ierr
      integer   :: recvQ_req(MPI_Process% nProcs-1,nEqn)
      integer   :: recvUx_req(MPI_Process% nProcs-1,nEqn)
      integer   :: recvUy_req(MPI_Process% nProcs-1,nEqn)
      integer   :: recvUz_req(MPI_Process% nProcs-1,nEqn)

      do domains = 2, MPI_Process% nProcs
         do i = 1, nEqn
            call mpi_irecv( IBMmask(domains)% Q(:,i)  , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recvQ_req(domains-1,i) , ierr)
            call mpi_irecv( IBMmask(domains)% U_x(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recvUx_req(domains-1,i), ierr)
            call mpi_irecv( IBMmask(domains)% U_y(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recvUy_req(domains-1,i), ierr)
            call mpi_irecv( IBMmask(domains)% U_z(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recvUz_req(domains-1,i), ierr)
         end do 
      end do 

      do i = 1, nEqn
         call mpi_waitall( MPI_Process% nProcs-1, recvQ_req(:,i) , MPI_STATUSES_IGNORE, ierr )
         call mpi_waitall( MPI_Process% nProcs-1, recvUx_req(:,i), MPI_STATUSES_IGNORE, ierr )
         call mpi_waitall( MPI_Process% nProcs-1, recvUy_req(:,i), MPI_STATUSES_IGNORE, ierr )
         call mpi_waitall( MPI_Process% nProcs-1, recvUz_req(:,i), MPI_STATUSES_IGNORE, ierr )
      end do 
#endif
   end subroutine recvStateMaskRoot

   subroutine sendStateMask2Partitions( IBMmask, domain, nEqn )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain, nEqn
#ifdef _HAS_MPI_
      integer   :: i, domains, domains_, msg, msg1, ierr 
      integer   :: sendQ_req(MPI_Process% nProcs-1,MPI_Process% nProcs,nEqn)
      integer   :: sendUx_req(MPI_Process% nProcs-1,MPI_Process% nProcs,nEqn)
      integer   :: sendUy_req(MPI_Process% nProcs-1,MPI_Process% nProcs,nEqn)
      integer   :: sendUz_req(MPI_Process% nProcs-1,MPI_Process% nProcs,nEqn)

      do domains = 2, MPI_Process% nProcs 
         do domains_ = 1, MPI_Process% nProcs 
            do i = 1, nEqn
               call mpi_isend( IBMmask(domains_)% Q(:,i)  , IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, sendQ_req(domains-1,domains_,i) , ierr ) 
               call mpi_isend( IBMmask(domains_)% U_x(:,i), IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, sendUx_req(domains-1,domains_,i), ierr ) 
               call mpi_isend( IBMmask(domains_)% U_y(:,i), IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, sendUy_req(domains-1,domains_,i), ierr ) 
               call mpi_isend( IBMmask(domains_)% U_z(:,i), IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, sendUz_req(domains-1,domains_,i), ierr ) 
            end do 
         end do
      end do 

      do msg = 1, MPI_Process% nProcs 
         do msg1 = 1, nEqn
            call mpi_waitall( MPI_Process% nProcs-1, sendQ_req(:,msg,msg1) , MPI_STATUSES_IGNORE, ierr )
            call mpi_waitall( MPI_Process% nProcs-1, sendUx_req(:,msg,msg1), MPI_STATUSES_IGNORE, ierr )
            call mpi_waitall( MPI_Process% nProcs-1, sendUy_req(:,msg,msg1), MPI_STATUSES_IGNORE, ierr )
            call mpi_waitall( MPI_Process% nProcs-1, sendUz_req(:,msg,msg1), MPI_STATUSES_IGNORE, ierr )
         end do 
      end do 
#endif
   end subroutine sendStateMask2Partitions

   subroutine recvStateMaskPartitions( IBMmask, domain, nEqn )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain, nEqn
#ifdef _HAS_MPI_
      integer   :: i, domains, msg, ierr 
      integer   :: recvQ_req(MPI_Process% nProcs,nEqn)
      integer   :: recvUx_req(MPI_Process% nProcs,nEqn)
      integer   :: recvUy_req(MPI_Process% nProcs,nEqn)
      integer   :: recvUz_req(MPI_Process% nProcs,nEqn)

      do domains = 1, MPI_Process% nProcs
         do i = 1, nEqn 
            call mpi_irecv( IBMmask(domains)% Q(:,i)  , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recvQ_req(domains,i) , ierr)
            call mpi_irecv( IBMmask(domains)% U_x(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recvUx_req(domains,i), ierr)
            call mpi_irecv( IBMmask(domains)% U_y(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recvUy_req(domains,i), ierr)
            call mpi_irecv( IBMmask(domains)% U_z(:,i), IBMmask(domains)% NumOfObjs, MPI_INT   , 0, MPI_ANY_TAG, MPI_COMM_WORLD, recvUz_req(domains,i), ierr)
         end do
      end do 

      do msg = 1, nEqn
         call mpi_waitall( MPI_Process% nProcs, recvQ_req(:,msg) , MPI_STATUSES_IGNORE, ierr )
         call mpi_waitall( MPI_Process% nProcs, recvUx_req(:,msg), MPI_STATUSES_IGNORE, ierr )
         call mpi_waitall( MPI_Process% nProcs, recvUy_req(:,msg), MPI_STATUSES_IGNORE, ierr )
         call mpi_waitall( MPI_Process% nProcs, recvUz_req(:,msg), MPI_STATUSES_IGNORE, ierr )
      end do
#endif
   end subroutine recvStateMaskPartitions

   subroutine castStateBandRegion( IBMmask, nEqn )
      
      implicit none 
      
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: nEqn
#ifdef _HAS_MPI_
      integer :: domains, i, ierr 

      do domains = 1, MPI_Process% nProcs 
         do i = 1, nEqn
            call mpi_bcast( IBMmask(domains)% Q(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         end do 
      end do
#endif 
   end subroutine castStateBandRegion

   subroutine castGradientsBandRegion( IBMmask, nEqn )
      
      implicit none 
      
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: nEqn
#ifdef _HAS_MPI_
      integer :: domains, i, ierr 
      
      do domains = 1, MPI_Process% nProcs 
         do i = 1, nEqn
            call mpi_bcast( IBMmask(domains)% U_x(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
            call mpi_bcast( IBMmask(domains)% U_y(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
            call mpi_bcast( IBMmask(domains)% U_z(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         end do 
      end do
#endif 
   end subroutine castGradientsBandRegion
!__________________________________________________

   subroutine castMaskPlot( IBMMask )
      
      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
! #ifdef _HAS_MPI_
!       integer :: domains, ierr

!       if ( .not. MPI_Process % doMPIAction ) return

!       do domains = 1, MPI_Process% nProcs 
!          call mpi_bcast( IBMmask(domains)% x(:)% isInsideBody, IBMmask(domains)% NumOfObjs, MPI_LOGICAL, domains-1, MPI_COMM_WORLD, ierr )     
!       end do 
! #endif 
   end subroutine castMaskPlot


   subroutine gatherMaskGeom( IBMmask )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: dist(:), normal(:,:) 
      integer                    :: i, domain, domains, ierr, index, msg, startIdx, endIdx
      integer                    :: gather_req(MPI_Process% nProcs, 4)

      if ( .not. MPI_Process % doMPIAction ) return

      domain = MPI_Process% rank + 1
      
      allocate( dist(IBMmask(domain)% NumOfObjs*MPI_Process% nProcs),       &
                normal(IBMmask(domain)% NumOfObjs*MPI_Process% nProcs,NDIM) )

      do domains = 1, MPI_Process% nProcs 
         call mpi_igather( IBMmask(domains)% dist        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE,  dist         ,                             &
                           IBMmask(domains)% NumOfObjs   , MPI_DOUBLE                 , domains-1 , MPI_COMM_WORLD, gather_req(domains,1), ierr )
         call mpi_igather( IBMmask(domains)% normal(:,IX), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, normal(:,IX)  ,                             &
                           IBMmask(domains)% NumOfObjs   , MPI_DOUBLE                 , domains-1 , MPI_COMM_WORLD, gather_req(domains,2), ierr )
         call mpi_igather( IBMmask(domains)% normal(:,IY), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, normal(:,IY)  ,                             &
                           IBMmask(domains)% NumOfObjs   , MPI_DOUBLE                 , domains-1 , MPI_COMM_WORLD, gather_req(domains,3), ierr )
         call mpi_igather( IBMmask(domains)% normal(:,IZ), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, normal(:,IZ)  ,                             &
                           IBMmask(domains)% NumOfObjs   , MPI_DOUBLE                 , domains-1 , MPI_COMM_WORLD, gather_req(domains,4), ierr )
      end do 

      do msg = 1, 4 
         call mpi_waitall( MPI_Process% nProcs, gather_req(:,msg), MPI_STATUSES_IGNORE, ierr )
      end do

      do i = 1, IBMmask(domain)% NumOfObjs 
         do domains = 1, MPI_Process% nProcs
            index = (domains-1)*IBMmask(domain)% NumOfObjs + i   
            if( IBMmask(domain)% dist(i) .gt. dist(index) ) then 
               IBMmask(domain)% dist(i)   = dist(index) 
               IBMmask(domain)% normal(i,:) = normal(index,:)
            end if 
         end do
      end do 

      deallocate( dist, normal )
#endif
   end subroutine gatherMaskGeom

   subroutine Mask2Root( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 

      if( .not. MPI_Process% isRoot ) then 
         call sendMask2Root( IBMmask, domain )
      end if 
      
      if( MPI_Process% isRoot ) then 
         call recvMask2Root( IBMmask, domain )
      end if 
      
   end subroutine Mask2Root

   subroutine Mask2Partitions( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 

      integer :: domains 

      if( MPI_Process% isRoot ) then 
         call sendMask2partitions( IBMmask, domain )
      end if

      if( .not. MPI_Process% isRoot ) then 
         call recvMask2Partitions( IBMmask, domain )
      end if

   end subroutine Mask2Partitions

   subroutine sendMask2Root( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: coords_x(:), coords_y(:), coords_z(:)
      integer,       allocatable :: local_position_x(:), local_position_y(:), local_position_z(:), element_index(:)
      integer                    :: NumOfObjs, send_req(8), ierr

      call mpi_isend( IBMmask(domain)% NumOfObjs           ,                          1, MPI_INT,    0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )
      call mpi_isend( IBMmask(domain)% coords(:,IX)        , IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr ) 
      call mpi_isend( IBMmask(domain)% coords(:,IY)        , IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr ) 
      call mpi_isend( IBMmask(domain)% coords(:,IZ)        , IBMmask(domain)% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr ) 
      call mpi_isend( IBMmask(domain)% element_index       , IBMmask(domain)% NumOfObjs, MPI_INT   , 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(5), ierr )
      call mpi_isend( IBMmask(domain)% local_position(:,IX), IBMmask(domain)% NumOfObjs, MPI_INT   , 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(6), ierr )
      call mpi_isend( IBMmask(domain)% local_position(:,IY), IBMmask(domain)% NumOfObjs, MPI_INT   , 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(7), ierr )
      call mpi_isend( IBMmask(domain)% local_position(:,IZ), IBMmask(domain)% NumOfObjs, MPI_INT   , 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(8), ierr )

      call mpi_waitall( 8, send_req, MPI_STATUSES_IGNORE, ierr )
#endif
   end subroutine sendMask2Root

   subroutine recvMask2Root( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer                    :: domains, NumOfObjs, msg, ierr
      integer                    :: start_index, final_index, biggerdomains, elems_per_domain(MPI_Process% nProcs)
      integer                    :: recvFirst_req(MPI_Process% nProcs-1), recv_req(MPI_Process% nProcs-1,7)


      integer  :: recv_req_(8)
      real(kind=RP), allocatable :: coords_x(:), coords_y(:), coords_z(:)
      integer,       allocatable :: local_position_x(:), local_position_y(:), local_position_z(:), element_index(:)

      do domains = 2, MPI_Process% nProcs 
         call mpi_irecv(IBMmask(domains)% NumOfObjs, 1, MPI_INT, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recvFirst_req(domains-1), ierr)
      end do

      call mpi_waitall(  MPI_Process% nProcs-1, recvFirst_req, MPI_STATUSES_IGNORE, ierr )

      do domains = 2, MPI_Process% nProcs 
         call IBMmask(domains)% build(IBMmask(domains)% NumOfObjs)
      end do

      do domains = 2, MPI_Process% nProcs
         call mpi_irecv( IBMmask(domains)% coords(:,IX)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,1), ierr)
         call mpi_irecv( IBMmask(domains)% coords(:,IY)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,2), ierr)
         call mpi_irecv( IBMmask(domains)% coords(:,IZ)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,3), ierr)
         call mpi_irecv( IBMmask(domains)% element_index       , IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,4), ierr)
         call mpi_irecv( IBMmask(domains)% local_position(:,IX), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,5), ierr)
         call mpi_irecv( IBMmask(domains)% local_position(:,IY), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,6), ierr)
         call mpi_irecv( IBMmask(domains)% local_position(:,IZ), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1,7), ierr)
      end do 

      do msg = 1, 7 
         call mpi_waitall( MPI_Process% nProcs-1, recv_req(:,msg), MPI_STATUSES_IGNORE, ierr )
      end do 
#endif
   end subroutine recvMask2Root

   subroutine sendMask2partitions( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer  :: send_req(MPI_Process% nProcs-1,MPI_Process% nProcs,7), sendFirst_req(MPI_Process% nProcs-1,MPI_Process% nProcs)
      integer  :: domains, domains_, msg, msg1, ierr

      do domains = 2, MPI_Process% nProcs 
         do domains_ = 1, MPI_Process% nProcs
            call mpi_isend( IBMmask(domains_)% NumOfObjs, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, sendFirst_req(domains-1,domains_), ierr )
         end do 
      end do

      do domains = 2, MPI_Process% nProcs 
         do domains_ = 1, MPI_Process% nProcs
            call mpi_isend( IBMmask(domains_)% coords(:,IX)        , IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,1), ierr ) 
            call mpi_isend( IBMmask(domains_)% coords(:,IY)        , IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,2), ierr ) 
            call mpi_isend( IBMmask(domains_)% coords(:,IZ)        , IBMmask(domains_)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,3), ierr ) 
            call mpi_isend( IBMmask(domains_)% element_index       , IBMmask(domains_)% NumOfObjs, MPI_INT   , domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,4), ierr )
            call mpi_isend( IBMmask(domains_)% local_position(:,IX), IBMmask(domains_)% NumOfObjs, MPI_INT   , domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,5), ierr )
            call mpi_isend( IBMmask(domains_)% local_position(:,IY), IBMmask(domains_)% NumOfObjs, MPI_INT   , domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,6), ierr )
            call mpi_isend( IBMmask(domains_)% local_position(:,IZ), IBMmask(domains_)% NumOfObjs, MPI_INT   , domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,domains_,7), ierr )
         end do
      end do 

      do msg = 1, MPI_Process% nProcs 
         call mpi_waitall(  MPI_Process% nProcs-1, sendFirst_req(:,msg), MPI_STATUSES_IGNORE, ierr )
         do msg1 = 1, 7
            call mpi_waitall(  MPI_Process% nProcs-1, send_req(:,msg,msg1), MPI_STATUSES_IGNORE, ierr )
         end do 
      end do 
#endif
   end subroutine sendMask2partitions

   subroutine recvMask2Partitions( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer                    :: recvFirst_req(MPI_Process% nProcs), ierr
      integer                    :: recv_req(MPI_Process% nProcs,7), domains, msg

      integer  :: recv_req_(7), NumOfObjs
      real(kind=RP), allocatable :: coords_x(:), coords_y(:), coords_z(:)
      integer,       allocatable :: local_position_x(:), local_position_y(:), local_position_z(:), element_index(:)

      do domains = 1, MPI_Process% nProcs 
         call mpi_irecv(IBMmask(domains)% NumOfObjs, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recvFirst_req(domains), ierr)
      end do

      call mpi_waitall( MPI_Process% nProcs, recvFirst_req, MPI_STATUSES_IGNORE, ierr )

      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         call IBMmask(domains)% build(IBMmask(domains)% NumOfObjs)
      end do

      do domains = 1, MPI_Process% nProcs
         call mpi_irecv( IBMmask(domains)% coords(:,IX)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,1), ierr)
         call mpi_irecv( IBMmask(domains)% coords(:,IY)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,2), ierr)
         call mpi_irecv( IBMmask(domains)% coords(:,IZ)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,3), ierr)
         call mpi_irecv( IBMmask(domains)% element_index       , IBMmask(domains)% NumOfObjs, MPI_INT   , 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,4), ierr)
         call mpi_irecv( IBMmask(domains)% local_position(:,IX), IBMmask(domains)% NumOfObjs, MPI_INT   , 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,5), ierr)
         call mpi_irecv( IBMmask(domains)% local_position(:,IY), IBMmask(domains)% NumOfObjs, MPI_INT   , 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,6), ierr)
         call mpi_irecv( IBMmask(domains)% local_position(:,IZ), IBMmask(domains)% NumOfObjs, MPI_INT   , 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,7), ierr)
      end do 

      do msg = 1, 7
         call mpi_waitall( MPI_Process% nProcs, recv_req(:,msg), MPI_STATUSES_IGNORE, ierr )
      end do
#endif
   end subroutine recvMask2Partitions

   subroutine MaskLogical2Root( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain

      if( .not. MPI_Process% isRoot ) then 
         call sendLogical2Root( IBMmask, domain )
      end if 

      if( MPI_Process% isRoot ) then 
         call recvLogical2Root( IBMmask, domain )
      end if  

   end subroutine MaskLogical2Root

   subroutine sendLogical2Root( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer  :: sendFirst_req, ierr

      call mpi_isend( IBMmask(domain)% isInsideBody, IBMmask(domain)% NumOfObjs, MPI_LOGICAL, 0, DEFAULT_TAG, MPI_COMM_WORLD, sendFirst_req, ierr )
      
      call mpi_wait(sendFirst_req, MPI_STATUS_IGNORE, ierr)
#endif
   end subroutine sendLogical2Root
   
   subroutine recvLogical2Root( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer  :: domains, recv_req(MPI_Process% nProcs-1), ierr

      do domains = 2, MPI_Process% nProcs 
         call mpi_irecv( IBMmask(domains)% isInsideBody, IBMmask(domains)% NumOfObjs, MPI_LOGICAL, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains-1), ierr)
      end do 

      call mpi_waitall( MPI_Process% nProcs-1, recv_req, MPI_STATUSES_IGNORE, ierr )
#endif
   end subroutine recvLogical2Root

   subroutine castMaskNumOfObjs( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer :: domains, ierr

      if ( .not. MPI_Process % doMPIAction ) return

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMmask(domains)% NumOfObjs, 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )    
      end do 

      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         call IBMmask(domains)% build(IBMmask(domains)% NumOfObjs)
      end do
#endif
   end subroutine castMaskNumOfObjs

   subroutine castMask( IBMmask, domain )
      
      implicit none 
      
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain
#ifdef _HAS_MPI_
      integer :: domains, ierr 

      if ( .not. MPI_Process % doMPIAction ) return

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMmask(domains)% coords(:,IX)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% coords(:,IY)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% coords(:,IZ)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% element_index       , IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% local_position(:,IX), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% local_position(:,IY), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% local_position(:,IZ), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
      end do 

#endif
   end subroutine castMask

   subroutine castIsInsideBody( IBMmask, domain )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain
#ifdef _HAS_MPI_
      integer :: domains, ierr 

      if ( .not. MPI_Process % doMPIAction ) return

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMmask(domains)% isInsideBody, IBMmask(domains)% NumOfObjs, MPI_LOGICAL, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% STLNum      , IBMmask(domains)% NumOfObjs, MPI_INT    , domains-1, MPI_COMM_WORLD, ierr )
      end do 
#endif
   end subroutine castIsInsideBody

   subroutine SendAxis( stl ) 

      implicit none 

      type(STLfile), intent(in) :: stl 
#ifdef _HAS_MPI_
      integer :: domains, send_req(MPI_Process% nProcs-1), &
                 msg, ierr

      if( .not. MPI_Process% isRoot ) return 
      
      do domains = 2, MPI_Process% nProcs
         call mpi_isend(stl% maxAxis, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1), ierr )
      end do

      call mpi_waitall( MPI_Process% nProcs-1, send_req, MPI_STATUSES_IGNORE, ierr )
#endif
   end subroutine SendAxis

   subroutine RecvAxis( stl )

      implicit none 

      type(STLfile), intent(inout) :: stl
#ifdef _HAS_MPI_
      integer :: recv_req, ierr
      
      if( MPI_Process% isRoot ) return 

      call mpi_irecv(stl% maxAxis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req, ierr)

      call mpi_wait(recv_req, MPI_STATUS_IGNORE, ierr)
#endif
   end subroutine RecvAxis


   subroutine SendOBB( STLNum )

      implicit none 

      integer, intent(in) :: STLNum 
#ifdef _HAS_MPI_
      integer              :: domains, i, ierr, domain, msg, &
                              send_req(MPI_Process% nProcs-1,17)

      if ( .not. MPI_Process% isRoot ) return

      do domains = 2, MPI_Process% nProcs

         call mpi_isend(OBB(STLNum)% MBR% angle, 1, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,1), ierr )

         call mpi_isend(OBB(STLNum)% MBR% center, NDIM-1, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,2), ierr )
         
         call mpi_isend(OBB(STLNum)% CloudCenter, NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,3), ierr )
     
         do i = 1, NDIM 
            call mpi_isend(OBB(STLNum)% R(:,i), NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,3+i), ierr )
         end do 

         do i = 1, NDIM 
            call mpi_isend(OBB(STLNum)% invR(:,i), NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,6+i), ierr )
         end do 

         do i = 1, BOXVERTICES
            call mpi_isend(OBB(STLNum)% LocVertices(:,i), NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,9+i), ierr )
         end do 
      
      end do 

      do msg = 1, 17 
         call mpi_waitall(MPI_Process% nProcs-1, send_req(:,msg), MPI_STATUSES_IGNORE, ierr)
      end do
#endif
   end subroutine SendOBB
   
   subroutine recvOBB( STLNum )

      implicit none 

      integer, intent(in) :: STLNum

      integer :: domain, i
#ifdef _HAS_MPI_
      integer :: domains, ierr, recv_req(17)

      if( MPI_Process% isRoot ) return 

      call mpi_irecv(OBB(STLNum) % MBR% angle,  1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr)
      
      call mpi_irecv(OBB(STLNum)%  MBR% center, NDIM-1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr)

      call mpi_irecv(OBB(STLNum)% CloudCenter, NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr)

      do i = 1, NDIM 
         call mpi_irecv(OBB(STLNum)% R(:,i), NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3+i), ierr)
      end do 

      do i = 1, NDIM 
         call mpi_irecv(OBB(STLNum)% invR(:,i), NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6+i), ierr)
      end do 

      do i = 1, BOXVERTICES 
         call mpi_irecv(OBB(STLNum)% LocVertices(:,i), NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9+i), ierr)
      end do 

      call mpi_waitall(17, recv_req, MPI_STATUSES_IGNORE, ierr)
#endif
   end subroutine recvOBB

   subroutine SendSTL2Partitions( STL, STLNum, axis )
      implicit none
      !-arguments---------------------------------------------------------------------
      type(STLfile), intent(inout) :: STL
      integer,       intent(in)    :: STLNum, axis
      !-local-variables---------------------------------------------------------------
#ifdef _HAS_MPI_
      real(kind=RP),     allocatable :: Bar(:), coord(:),                         &
                                        normals_x(:), normals_y(:), normals_z(:), &
                                        vertices_x(:,:), vertices_y(:,:),         &
                                        vertices_z(:,:)
      integer                        :: NumOfObjs,                                &
                                        start_index, final_index, i, j,           &
                                        nProcs, ierr, biggerdomains,              &
                                        elems_per_domain(MPI_Process% nProcs),    &
                                        msg, send_req(MPI_Process% nProcs-1,12),  &
                                        NumOfObjsPartion(MPI_Process% nProcs)
      integer, allocatable           :: SortedIndex(:)

      STL% partition = MPI_Process% rank

      NumOfObjs = STL% NumOfObjs

      allocate( Bar(NumOfObjs),                      &
                coord(NumOfObjs),                    &
                SortedIndex(NumOfObjs),              &
                normals_x(NumOfObjs),                &
                normals_y(NumOfObjs),                &
                normals_z(NumOfObjs),                &
                vertices_x(NumOfObjs,NumOfVertices), &
                vertices_y(NumOfObjs,NumOfVertices), &
                vertices_z(NumOfObjs,NumOfVertices)  )

      do i = 1, NumOfObjs
         SortedIndex(i) = STL% ObjectsList(i)% index
         Bar(i)         = maxval(STL% ObjectsList(i)% vertices(1:NumOfVertices)% coords(axis))
      end do

      call sort( Bar, SortedIndex, coord, coord, coord, 1, NumOfObjs )

      deallocate(Bar,coord)

      elems_per_domain = NumOfObjs/MPI_Process% nProcs
      biggerdomains    = mod(NumOfObjs,MPI_Process% nProcs)
      elems_per_domain(1:biggerdomains) = elems_per_domain(1:biggerdomains) + 1

      start_index = 1
      do nProcs = 1, MPI_Process% nProcs
         final_index = start_index + elems_per_domain(nProcs) - 1

         do i = start_index, final_index
            do j = 1, NumOfVertices
               vertices_x(i,j) = STL% ObjectsList(SortedIndex(i))% vertices(j)% coords(1)
               vertices_y(i,j) = STL% ObjectsList(SortedIndex(i))% vertices(j)% coords(2)
               vertices_z(i,j) = STL% ObjectsList(SortedIndex(i))% vertices(j)% coords(3)
            end do
            normals_x(i) = STL% ObjectsList(SortedIndex(i))% normal(1)
            normals_y(i) = STL% ObjectsList(SortedIndex(i))% normal(2)
            normals_z(i) = STL% ObjectsList(SortedIndex(i))% normal(3)
         end do  

         if( STL% show ) call DescribeSTLPartitions(nProcs-1,(final_index-start_index+1))

         NumOfObjsPartion(nProcs) = final_index-start_index+1

         start_index = final_index + 1

      end do
      
      do nProcs = 2, MPI_Process% nProcs
         call mpi_send(NumOfObjsPartion(nProcs), 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, ierr )
      end do
      
      start_index = elems_per_domain(1) + 1

      do nProcs = 2, MPI_Process% nProcs

         final_index = start_index + elems_per_domain(nProcs) - 1

         call mpi_isend(normals_x(start_index:final_index), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )

         call mpi_isend(normals_y(start_index:final_index), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )

         call mpi_isend(normals_z(start_index:final_index), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )

         call mpi_isend(vertices_x(start_index:final_index,1), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )

         call mpi_isend(vertices_y(start_index:final_index,1), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,5), ierr )

         call mpi_isend(vertices_z(start_index:final_index,1), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,6), ierr )

         call mpi_isend(vertices_x(start_index:final_index,2), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,7), ierr )

         call mpi_isend(vertices_y(start_index:final_index,2), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,8), ierr )

         call mpi_isend(vertices_z(start_index:final_index,2), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,9), ierr )

         call mpi_isend(vertices_x(start_index:final_index,3), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,10), ierr )

         call mpi_isend(vertices_y(start_index:final_index,3), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,11), ierr )

         call mpi_isend(vertices_z(start_index:final_index,3), NumOfObjsPartion(nProcs), MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,12), ierr )

         start_index = final_index + 1

      end do   

      do msg = 1, 12
         call mpi_waitall(MPI_Process% nProcs-1, send_req(:,msg), MPI_STATUSES_IGNORE, ierr)
      end do

      call STL% destroy()

      allocate(STL% ObjectsList(elems_per_domain(1)))

      do i = 1, elems_per_domain(1)
         allocate(STL% ObjectsList(i)% vertices(NumOfVertices))
         STL% ObjectsList(i)% normal(1) = normals_x(i)
         STL% ObjectsList(i)% normal(2) = normals_y(i)
         STL% ObjectsList(i)% normal(3) = normals_z(i)
         do j = 1, NumOfVertices
            STL% ObjectsList(i)% vertices(j)% coords(1) = vertices_x(i,j)
            STL% ObjectsList(i)% vertices(j)% coords(2) = vertices_y(i,j)
            STL% ObjectsList(i)% vertices(j)% coords(3) = vertices_z(i,j)
         end do 
         STL% ObjectsList(i)% index = i
      end do 

      STL% NumOfObjs = elems_per_domain(1)

      deallocate(vertices_x, vertices_y, vertices_z, normals_x, normals_y, normals_z, SortedIndex)
#endif
   end subroutine SendSTL2Partitions

   subroutine receiveSTLpartitions( STL, STLNum )
      implicit none
      !-arguments-----------------------------------------------------
      type(STLfile), intent(inout) :: STL
      integer,       intent(in)    :: STLNum
#ifdef _HAS_MPI_
      !-local-variables-----------------------------------------------
      real(kind=RP), allocatable :: normals_x(:), normals_y(:), normals_z(:), &
                                    vertices_x(:,:), vertices_y(:,:),         &
                                    vertices_z(:,:)
      integer                    :: NumOfObjs, ierr, recv_req(12), i, j

      if( MPI_Process% isRoot ) return

      call mpi_recv( STL% NumOfObjs, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

      NumOfObjs = STL% NumOfObjs

      allocate( normals_x(NumOfObjs),       &
                normals_y(NumOfObjs),       &
                normals_z(NumOfObjs),       &
                vertices_x(NumOfObjs,NDIM), &
                vertices_y(NumOfObjs,NDIM), &
                vertices_z(NumOfObjs,NDIM)  )

      call mpi_irecv( normals_x      , NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )  
      call mpi_irecv( normals_y      , NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr )  
      call mpi_irecv( normals_z      , NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr )  
      call mpi_irecv( vertices_x(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr )  
      call mpi_irecv( vertices_y(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr )
      call mpi_irecv( vertices_z(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr )  
      call mpi_irecv( vertices_x(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr )
      call mpi_irecv( vertices_y(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr )
      call mpi_irecv( vertices_z(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9), ierr )
      call mpi_irecv( vertices_x(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(10), ierr )
      call mpi_irecv( vertices_y(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(11), ierr )
      call mpi_irecv( vertices_z(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(12), ierr )

      call mpi_waitall(12, recv_req, MPI_STATUSES_IGNORE, ierr)

      allocate(STL% ObjectsList(NumOfObjs))

      do i = 1, NumOfObjs
         STL% ObjectsList(i)% NumOfVertices = NumOfVertices
         allocate(STL% ObjectsList(i)% vertices(NumOfVertices))
         do j = 1, NumOfVertices
            STL% ObjectsList(i)% vertices(j)% coords(1) = vertices_x(i,j)
            STL% ObjectsList(i)% vertices(j)% coords(2) = vertices_y(i,j)
            STL% ObjectsList(i)% vertices(j)% coords(3) = vertices_z(i,j)
         end do
         STL% ObjectsList(i)% normal(1) = normals_x(i)
         STL% ObjectsList(i)% normal(2) = normals_y(i)
         STL% ObjectsList(i)% normal(3) = normals_z(i)
         STL% ObjectsList(i)% index = i
      end do
      
      deallocate( normals_x,  &
                  normals_y,  &
                  normals_z,  &
                  vertices_x, &
                  vertices_y, &
                  vertices_z  )
#endif
   end subroutine receiveSTLpartitions

   subroutine sendSTLRoot( STL )

      implicit none 

      type(STLfile), intent(in) :: STL 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: normals_x(:), normals_y(:), normals_z(:),   &
                                    vertices_x(:,:), vertices_y(:,:), vertices_z(:,:)
      integer                    :: NumOfObjs, i, j, ierr, send_req(13), &
                                    array_of_statuses(MPI_STATUS_SIZE,13)

      if( MPI_Process% isRoot ) return

      NumOfObjs = stl% NumOfObjs

      allocate( normals_x(NumOfObjs),                & 
                normals_y(NumOfObjs),                &
                normals_z(NumOfObjs),                &
                vertices_x(NumOfObjs,NumOfVertices), &
                vertices_y(NumOfObjs,NumOfVertices), &
                vertices_z(NumOfObjs,NumOfVertices)  ) 
      !MUST BE IN THE GLOBAL REF FRAME
      do i = 1, NumOfObjs
         normals_x(i) = STL% ObjectsList(i)% normal(IX) 
         normals_y(i) = STL% ObjectsList(i)% normal(IY)  
         normals_z(i) = STL% ObjectsList(i)% normal(IZ)  
         do j = 1, NumOfVertices
            vertices_x(i,j) = STL% ObjectsList(i)% vertices(j)% coords(IX) 
            vertices_y(i,j) = STL% ObjectsList(i)% vertices(j)% coords(IY) 
            vertices_z(i,j) = STL% ObjectsList(i)% vertices(j)% coords(IZ) 
         end do 
      end do 
      
      call mpi_isend(NumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )

      call mpi_wait(send_req(1),MPI_STATUS_IGNORE,ierr)

      call mpi_isend(normals_x, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )

      call mpi_isend(normals_y, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )

      call mpi_isend(normals_z, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )

      call mpi_isend(vertices_x(:,1), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(5), ierr )

      call mpi_isend(vertices_y(:,1), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(6), ierr )

      call mpi_isend(vertices_z(:,1), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(7), ierr )

      call mpi_isend(vertices_x(:,2), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(8), ierr )

      call mpi_isend(vertices_y(:,2), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(9), ierr )

      call mpi_isend(vertices_z(:,2), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(10), ierr )

      call mpi_isend(vertices_x(:,3), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(11), ierr )

      call mpi_isend(vertices_y(:,3), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(12), ierr )

      call mpi_isend(vertices_z(:,3), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(13), ierr )

      call mpi_waitall(13, send_req, array_of_statuses, ierr)
#endif
   end subroutine sendSTLRoot

   subroutine recvSTLRootandPlot( STL, iter )

      implicit none 

      type(STLfile), intent(inout) :: STL 
      integer,       intent(in)    :: iter
#ifdef _HAS_MPI_
      type(STLfile)              :: STLplot 
      real(kind=RP), allocatable :: normals_x(:), normals_y(:), normals_z(:),          & 
                                    vertices_x(:,:), vertices_y(:,:), vertices_z(:,:)
      integer,       allocatable :: NumOfPartitionObjs(:), recv_req(:,:)
      integer                    :: nProcs, ierr, array_of_statuses(MPI_STATUS_SIZE,13), &
                                    i, j, start 

      if( .not. MPI_Process% isRoot ) return

      allocate( NumOfPartitionObjs(MPI_Process% nProcs-1), &
                recv_req(MPI_Process% nProcs-1,13)         )

      do nProcs = 2, MPI_Process% nProcs 
         call mpi_irecv( NumOfPartitionObjs(nProcs-1), 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )  
         call mpi_wait(recv_req(nProcs-1,1), MPI_STATUS_IGNORE, ierr)
      end do 

      STLplot% NumOfObjs = sum(NumOfPartitionObjs) + STL% NumOfObjs
      STLplot% filename  = STL% filename
      
      allocate(STLplot% ObjectsList(STLplot% NumOfObjs))

      do i = 1, STL% NumOfObjs
         STLplot% ObjectsList(i)% normal = STL% ObjectsList(i)% normal 
         allocate(STLplot% ObjectsList(i)% vertices(NumOfVertices))
         do j = 1, NumOfVertices
            STLplot% ObjectsList(i)% vertices(j)% coords = STL% ObjectsList(i)% vertices(j)% coords 
         end do 
      end do 
         
      start = STL% NumOfObjs 

      do nProcs =  2, MPI_Process% nProcs 

         allocate( normals_x(NumOfPartitionObjs(nProcs-1)),                &
                   normals_y(NumOfPartitionObjs(nProcs-1)),                &
                   normals_z(NumOfPartitionObjs(nProcs-1)),                &
                   vertices_x(NumOfPartitionObjs(nProcs-1),NumOfVertices), &
                   vertices_y(NumOfPartitionObjs(nProcs-1),NumOfVertices), &
                   vertices_z(NumOfPartitionObjs(nProcs-1),NumOfVertices)  )

         call mpi_irecv( normals_x, NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )  

         call mpi_irecv( normals_y, NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )  
                      
         call mpi_irecv( normals_z, NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )  

         call mpi_irecv( vertices_x(:,1), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,5), ierr )  

         call mpi_irecv( vertices_y(:,1), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,6), ierr )  
                      
         call mpi_irecv( vertices_z(:,1), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,7), ierr ) 
         
         call mpi_irecv( vertices_x(:,2), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,8), ierr )  

         call mpi_irecv( vertices_y(:,2), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,9), ierr )  

         call mpi_irecv( vertices_z(:,2), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,10), ierr ) 
         
         call mpi_irecv( vertices_x(:,3), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,11), ierr )  

         call mpi_irecv( vertices_y(:,3), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,12), ierr )  
                      
         call mpi_irecv( vertices_z(:,3), NumOfPartitionObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,13), ierr ) 

         call mpi_waitall(13, recv_req(nProcs-1,:), array_of_statuses, ierr)

         do i = 1, NumOfPartitionObjs(nProcs-1)
            start = start + 1
            STLplot% ObjectsList(start)% normal(1) = normals_x(i)
            STLplot% ObjectsList(start)% normal(2) = normals_y(i)
            STLplot% ObjectsList(start)% normal(3) = normals_z(i)
            allocate(STLplot% ObjectsList(start)% vertices(NumOfVertices))
            do j = 1, NumOfVertices
               STLplot% ObjectsList(start)% vertices(j)% coords(1) = vertices_x(i,j) 
               STLplot% ObjectsList(start)% vertices(j)% coords(2) = vertices_y(i,j) 
               STLplot% ObjectsList(start)% vertices(j)% coords(3) = vertices_z(i,j) 
            end do 
         end do 

         deallocate( normals_x,  &
                     normals_y,  & 
                     normals_z,  &  
                     vertices_x, & 
                     vertices_y, & 
                     vertices_z  )

      end do 

      deallocate( NumOfPartitionObjs, &
                  recv_req            )

      call STLplot% plot( iter )
      call STLplot% destroy()
#else 
      call STL% plot( iter )
#endif
   end subroutine recvSTLRootandPlot

   subroutine plotSTL( STL, iter )

      implicit none 

      type(STLfile), intent(inout) :: STL 
      integer,       intent(in)    :: iter

      if( MPI_Process% doMPIAction ) then 
         call sendSTLRoot( STL )
      end if 

      if( MPI_Process% isRoot ) then 
         call recvSTLRootandPlot( STL, iter )
      end if 

   end subroutine plotSTL
!____________________________________________________________________________________________________________________________

   subroutine MPIProcedures_IBM_HO_faces( IBMStencilPoints, nEqn )

      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      integer :: nEqn 
     
      call CastHOfacesNumOfObjs( IBMStencilPoints, nEqn )
      call CastHOfaces( IBMStencilPoints )

   end subroutine MPIProcedures_IBM_HO_faces

   subroutine Set_IBM_HO_faces( IBMStencilPoints, faces, nEqn )

      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)     
      type(face),      intent(inout) :: faces(:)
      integer,         intent(in)    :: nEqn 

      integer :: i, j, k, fID, N, domain, counter, counter1

      domain = MPI_Process% rank + 1 

      IBMStencilPoints(domain)% NumOfObjs = 0 

      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then 
            IBMStencilPoints(domain)% NumOfObjs = IBMStencilPoints(domain)% NumOfObjs + max(f% Nf(1)+1,f% Nf(2)+1) * (f% Nf(1)+1) * (f% Nf(2)+1)
         end if 
         end associate 
      end do

      call IBMStencilPoints(domain)% buildstencil( IBMStencilPoints(domain)% NumOfObjs ) 
      call IBMStencilPoints(domain)% buildstate(   IBMStencilPoints(domain)% NumOfObjs, nEqn ) 

      counter = 0 

      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then 
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               N = f% stencil(i,j)% N
               allocate( f% stencil(i,j)% Qsb(nEqn),  &
                         f% stencil(i,j)% Q(nEqn,0:N) )
               do k = 0, N
                  counter = counter + 1
                  IBMStencilPoints(domain)% coords(counter,:)         = f% stencil(i,j)% x_s(:,k)
                  IBMStencilPoints(domain)% N(counter)                = f% stencil(i,j)% N
                  IBMStencilPoints(domain)% local_position(counter,:) = (/i,j,k/) 
                  IBMStencilPoints(domain)% fIDs(counter)             = fID
                  IBMStencilPoints(domain)% STLNum(counter)           = f% STLNum 
               end do
            end do; end do 
         end if 
         end associate 
      end do
      
   end subroutine Set_IBM_HO_faces

   subroutine CastHOfacesNumOfObjs( IBMStencilPoints, nEqn )

      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      integer,         intent(in)    :: nEqn
#ifdef _HAS_MPI_
      integer :: domain, domains, ierr

      if( .not. MPI_Process% doMPIAction ) return

      domain = MPI_Process% rank + 1 
      
      do domains = 1, MPI_Process% nProcs   
         call mpi_bcast( IBMStencilPoints(domains)% NumOfObjs, 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )     
      end do 
      
      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         call IBMStencilPoints(domains)% buildstencil( IBMStencilPoints(domains)% NumOfObjs ) 
         call IBMStencilPoints(domains)% buildstate(   IBMStencilPoints(domains)% NumOfObjs, nEqn ) 
      end do 
#endif
   end subroutine CastHOfacesNumOfObjs

   subroutine CastHOfaces( IBMStencilPoints )

      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
#ifdef _HAS_MPI_
      integer :: m, i, j, k, domain, domains, counter, counter1, ierr

      if( .not. MPI_Process% doMPIAction ) return

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMStencilPoints(domains)% coords(:,IX)        , IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )       
         call mpi_bcast( IBMStencilPoints(domains)% coords(:,IY)        , IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )       
         call mpi_bcast( IBMStencilPoints(domains)% coords(:,IZ)        , IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )        
         call mpi_bcast( IBMStencilPoints(domains)% N                   , IBMStencilPoints(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )       
         call mpi_bcast( IBMStencilPoints(domains)% STLNum              , IBMStencilPoints(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )       
         call mpi_bcast( IBMStencilPoints(domains)% local_position(:,IZ), IBMStencilPoints(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )       
      end do  
#endif 
   end subroutine CastHOfaces

!================================================== HO faces state ==================================================
   subroutine IBM_HO_findElements( IBMStencilPoints, elements, NumOfSTL, clipAxis )
      use ElementClass
      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: NumOfSTL, clipAxis

      real(kind=RP) :: xi(NDIM), coords(NDIM)
      logical       :: FOUND
      integer       :: domain, domains, fID, i, j, k, eID, counter, counter1, N, STLNum  

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs
         IBMStencilPoints(domains)% element_index = 0  
         do counter = 1, IBMStencilPoints(domains)% NumOfObjs 
            do eID = 1, size(elements)
               associate( e => elements(eID) )
               if( e% HO_IBM ) cycle
               do STLNum = 1, NumOfSTL
                  coords = IBMStencilPoints(domains)% coords(counter,:)
                  !if( OBB(STLNum)% isPointInside(coords, 2.0_RP) ) then 
                     if( e% FindPointWithCoords(coords, clipAxis, xi) ) then
                        IBMStencilPoints(domains)% domain(counter)        = domain
                        IBMStencilPoints(domains)% element_index(counter) = eID 
                        IBMStencilPoints(domains)% xi(counter,:)          = xi 
                        exit 
                     end if   
                  !end if
               end do 
               end associate
            end do 
         end do 
      end do
      
   end subroutine IBM_HO_findElements

   subroutine IBM_HO_GetState( IBMStencilPoints, elements, nEqn )
      use ElementClass
      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: nEqn

      real(kind=RP) :: xi(NDIM)
      integer       :: domain, domains, i, eID 

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs 
         do i = 1, IBMStencilPoints(domains)% NumOfObjs
            if( IBMStencilPoints(domains)% domain(i) .eq. domain ) then 
               eID = IBMStencilPoints(domains)% element_index(i) 
               xi  = IBMStencilPoints(domains)% xi(i,:)
               IBMStencilPoints(domains)% Q(i,:) = elements(eID)% EvaluateSolutionAtPoint(nEqn, xi)
            end if 
         end do
      end do

   end subroutine IBM_HO_GetState

   subroutine GatherHOfacesState( IBMStencilPoints, nEqn, faces )
      use PhysicsStorage
      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      integer,         intent(in)    :: nEqn
      type(face),      intent(inout) :: faces(:)  
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: Q(:,:)
      integer,       allocatable :: eIDsdomain(:)
      integer                    :: domain, domains, i, j, k, fID, n, index
      integer                    :: msg, ierr, gather_req(MPI_Process% nProcs,nEqn+1)

      domain = MPI_Process% rank + 1

      if( MPI_Process% doMPIAction ) then 
      
         allocate( Q(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),    &
                   eIDsdomain(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs) )

         do domains = 1, MPI_Process% nProcs 
            do i = 1, nEqn
               call mpi_igather( IBMStencilPoints(domains)% Q(:,i), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, Q(:,i),            &
                                IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, gather_req(domains,i), ierr )
            end do 
            call mpi_igather( IBMStencilPoints(domains)% domain, IBMStencilPoints(domains)% NumOfObjs, MPI_INT, eIDsdomain,             &
                             IBMStencilPoints(domains)% NumOfObjs, MPI_INT, domains-1, MPI_COMM_WORLD, gather_req(domains,nEqn+1), ierr )
         end do
         
         do msg = 1, nEqn+1
            call mpi_waitall( MPI_Process% nProcs, gather_req(:,msg), MPI_STATUS_IGNORE, ierr )
         end do

         do n = 1, IBMStencilPoints(domain)% NumOfObjs
            do domains = 1, MPI_Process% nProcs 
               index = (domains-1)*IBMStencilPoints(domain)% NumOfObjs + n 
               if( eIDsdomain(index) .ne. 0 ) then 
                  fID = IBMStencilPoints(domain)% fIDs(n)
                  i   = IBMStencilPoints(domain)% local_position(n,IX)
                  j   = IBMStencilPoints(domain)% local_position(n,IY)
                  k   = IBMStencilPoints(domain)% local_position(n,IZ) 

                  faces(fID)% stencil(i,j)% Q(:,k) = Q(index,:)  
               end if 
            end do
         end do
         
         deallocate( Q, eIDsdomain )

      else 

         do n = 1, IBMStencilPoints(domain)% NumOfObjs
            fID = IBMStencilPoints(domain)% fIDs(n)
            i   = IBMStencilPoints(domain)% local_position(n,IX)
            j   = IBMStencilPoints(domain)% local_position(n,IY)
            k   = IBMStencilPoints(domain)% local_position(n,IZ)

            faces(fID)% stencil(i,j)% Q(:,k) = IBMStencilPoints(domain)% Q(n,:)
         end do 

      end if
#endif
   end subroutine GatherHOfacesState

   subroutine sendHO_IBMfacesMPI( faces, domain, MPIfaces )
      use MPI_Face_Class
      implicit none 

      type(face),           intent(in) :: faces(:) 
      integer,              intent(in) :: domain 
      type(MPI_FacesSet_t), intent(in) :: MPIfaces
   
      logical :: HO_IBM(MPIfaces% faces(domain)% no_of_faces)
      integer :: STLNum(MPIfaces% faces(domain)% no_of_faces)
      integer :: mpifID, fID 
#ifdef _HAS_MPI_
      integer :: send_req(2), ierr

      HO_IBM = .false. 
   
      if( MPIfaces% faces(domain)% no_of_faces .eq. 0 ) return 

      do mpifID = 1, MPIfaces% faces(domain)% no_of_faces
         fID = MPIfaces% faces(domain)% faceIDs(mpifID)
         if( faces(fID)% HO_IBM ) then 
            HO_IBM(mpifID) = .true. 
            STLNum(mpifID) = faces(fID)% STLNum
         else 
            HO_IBM(mpifID)= .false.
            STLNum(mpifID) = 0
         end if 
      end do

      call mpi_isend( HO_IBM, MPIfaces% faces(domain)% no_of_faces, MPI_LOGICAL, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )
      call mpi_isend( STLNum, MPIfaces% faces(domain)% no_of_faces, MPI_INT    , domain-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )

      call mpi_waitall(2, send_req, MPI_STATUS_IGNORE, ierr)
#endif
   end subroutine sendHO_IBMfacesMPI

   subroutine recvHO_IBMfacesMPI( faces, domain, MPIfaces, additional_NumOfMaskObjs )
      use MPI_Face_Class
      implicit none 

      type(face),           intent(inout) :: faces(:) 
      integer,              intent(in)    :: domain 
      type(MPI_FacesSet_t), intent(in)    :: MPIfaces
      integer,              intent(inout) :: additional_NumOfMaskObjs

      logical :: HO_IBM(MPIfaces% faces(domain)% no_of_faces)
      integer :: STLNum(MPIfaces% faces(domain)% no_of_faces)
      integer :: mpifID, fID, i, j
#ifdef _HAS_MPI_
      integer :: recv_req(2), ierr

      if( MPIfaces% faces(domain)% no_of_faces .eq. 0 ) return

      call mpi_irecv( HO_IBM, MPIfaces% faces(domain)% no_of_faces, MPI_LOGICAL, domain-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )
      call mpi_irecv( STLNum, MPIfaces% faces(domain)% no_of_faces, MPI_INT    , domain-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr )
   
      call mpi_waitall(2, recv_req, MPI_STATUS_IGNORE, ierr)

      do mpifID = 1, MPIfaces% faces(domain)% no_of_faces
         fID = MPIfaces% faces(domain) % faceIDs(mpifID)
         if( HO_IBM(mpifID) ) then
            faces(fID)% HO_IBM = HO_IBM(mpifID)
            faces(fID)% HOSIDE = maxloc(faces(fID)% elementIDs, dim=1) 
            faces(fID)% STLNum = STLNum(mpifID) 
            allocate(faces(fID)% stencil(0:faces(fID)% Nf(1),0:faces(fID)% Nf(2)))
            do j = 0, faces(fID)% Nf(2); do i = 0, faces(fID)% Nf(1)
               faces(fID)% stencil(i,j)% x = faces(fID)% geom% x(:,i,j)
               additional_NumOfMaskObjs    = additional_NumOfMaskObjs + 1
            end do; end do
         end if 
      end do
#endif
   end subroutine recvHO_IBMfacesMPI

   subroutine FixingmpiFaces( faces, MPIfaces, NumOfMaskObjs, MPIfixed )
      use MPI_Face_Class
      use MeshTypes, only: HMESH_MPI
      implicit none  
 
      type(face),           intent(inout) :: faces(:) 
      type(MPI_FacesSet_t), intent(inout) :: MPIfaces
      integer,              intent(inout) :: NumOfMaskObjs
      logical,              intent(inout) :: MPIfixed
#ifdef _HAS_MPI_
      integer :: domains, additional_NumOfMaskObjs, ierr 
      
      if( .not. MPI_Process% doMPIAction ) return 

      if( MPIfixed ) return 
      
      do domains = 1, MPI_Process% nProcs 
         call sendHO_IBMfacesMPI( faces, domains, MPIfaces )
      end do
      additional_NumOfMaskObjs = 0 
      do domains = 1, MPI_Process% nProcs 
         call recvHO_IBMfacesMPI( faces, domains, MPIfaces, additional_NumOfMaskObjs )
      end do

      NumOfMaskObjs = NumOfMaskObjs + additional_NumOfMaskObjs

      MPIfixed = .true. 
#endif
   end subroutine FixingmpiFaces

!__________________________________________________Procedures for HO integration____________________________________________________

   subroutine MPIProcedures_IBM_HOIntegrationPoints( IBM_HOIntegrationPoints )

      implicit none 

      type(IBMpoints), intent(inout) :: IBM_HOIntegrationPoints(:) 
     
      call CastHOIntegrationPointsNumOfObjs( IBM_HOIntegrationPoints )
      call CastHOIntegrationPoints( IBM_HOIntegrationPoints )

   end subroutine MPIProcedures_IBM_HOIntegrationPoints
   
   subroutine CastHOIntegrationPointsNumOfObjs( IBM_HOIntegrationPoints )

      implicit none 

      type(IBMpoints), intent(inout) :: IBM_HOIntegrationPoints(:)
#ifdef _HAS_MPI_      
      integer :: domain, domains, ierr

      if( .not. MPI_Process% doMPIAction ) return 

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBM_HOIntegrationPoints(domains)% NumOfObjs, 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )  
      end do 

      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         call IBM_HOIntegrationPoints(domains)% build(IBM_HOIntegrationPoints(domains)% NumOfObjs)
      end do
#endif 
   end subroutine CastHOIntegrationPointsNumOfObjs

   subroutine CastHOIntegrationPoints( IBM_HOIntegrationPoints )

      implicit none 

      type(IBMpoints), intent(inout) :: IBM_HOIntegrationPoints(:)
#ifdef _HAS_MPI_      
      integer :: domains, ierr

      if( .not. MPI_Process% doMPIAction ) return

      ! do domains = 1, MPI_Process% nProcs 
      !    call mpi_bcast( IBM_HOIntegrationPoints(domains)% x(:)% coords(IX), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )  
      !    call mpi_bcast( IBM_HOIntegrationPoints(domains)% x(:)% coords(IY), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )  
      !    call mpi_bcast( IBM_HOIntegrationPoints(domains)% x(:)% coords(IZ), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )  
      ! end do 
#endif
   end subroutine CastHOIntegrationPoints

   subroutine IBM_HO_GetGradient( IBMStencilPoints, elements, nEqn )
      use ElementClass
      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: nEqn

      real(kind=RP) :: xi(NDIM)
      integer       :: domain, domains, i, eID 

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs 
         do i = 1, IBMStencilPoints(domains)% NumOfObjs
            if( IBMStencilPoints(domains)% domain(i) .eq. domain ) then 
               eID = IBMStencilPoints(domains)% element_index(i) 
               xi  = IBMStencilPoints(domains)% xi(i,:)
               IBMStencilPoints(domains)% U_x(i,:) = elements(eID)% EvaluateGradientAtPoint(nEqn, xi, IX)
               IBMStencilPoints(domains)% U_y(i,:) = elements(eID)% EvaluateGradientAtPoint(nEqn, xi, IY)
               IBMStencilPoints(domains)% U_z(i,:) = elements(eID)% EvaluateGradientAtPoint(nEqn, xi, IZ)
            end if 
         end do
      end do

   end subroutine IBM_HO_GetGradient

   subroutine GatherHOIntegrationPointsState( IBMStencilPoints, ObjectsList, zoneMask, nEqn )
      use PhysicsStorage
      implicit none 

      type(IBMpoints),   intent(inout) :: IBMStencilPoints(:)
      type(Object_type), intent(inout) :: ObjectsList(:)
      logical,           intent(in)    :: zoneMask
      integer,           intent(in)    :: nEqn
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: Q(:,:), U_x(:,:), U_y(:,:), U_z(:,:)
      integer,       allocatable :: eIDsdomain(:)
      integer                    :: domain, domains, i, j, k, n, index
      integer                    :: msg, ierr, gather_req(MPI_Process% nProcs,4*nEqn+1)

      domain = MPI_Process% rank + 1

         allocate( Q(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),    &
                   U_x(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),  &
                   U_y(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),  &
                   U_z(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),  &
                   eIDsdomain(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs) )

         do domains = 1, MPI_Process% nProcs 
            do i = 1, nEqn
               call mpi_igather( IBMStencilPoints(domains)% Q(:,i), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, Q(:,i),     &
                                 IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD,                     &
                                 gather_req(domains,i), ierr                                                                      )
               call mpi_igather( IBMStencilPoints(domains)% U_x(:,i), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, U_x(:,i), &
                                 IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD,                     &
                                 gather_req(domains,i+nEqn), ierr                                                                 )
               call mpi_igather( IBMStencilPoints(domains)% U_y(:,i), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, U_y(:,i), &
                                 IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD,                     &
                                 gather_req(domains,i+2*nEqn), ierr                                                               )
               call mpi_igather( IBMStencilPoints(domains)% U_z(:,i), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, U_z(:,i), &
                                 IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD,                     &
                                 gather_req(domains,i+3*nEqn), ierr                                                               )
            end do 
            call mpi_igather( IBMStencilPoints(domains)% domain, IBMStencilPoints(domains)% NumOfObjs, MPI_INT, eIDsdomain, &
                             IBMStencilPoints(domains)% NumOfObjs, MPI_INT, domains-1, MPI_COMM_WORLD,                      &
                             gather_req(domains,4*nEqn+1), ierr                                                             )
         end do

         do msg = 1, 4*nEqn+1
            call mpi_waitall( MPI_Process% nProcs, gather_req(:,msg), MPI_STATUS_IGNORE, ierr )
         end do
 
         if( zoneMask ) then 
            do n = 1, IBMStencilPoints(domain)% NumOfObjs
               do domains = 1, MPI_Process% nProcs 
                  index = (domains-1)*IBMStencilPoints(domain)% NumOfObjs + n 
                  if( eIDsdomain(index) .ne. 0 ) then 
                     i = IBMStencilPoints(domain)% local_position(n,IX)
                     j = IBMStencilPoints(domain)% local_position(n,IY)
                     k = IBMStencilPoints(domain)% local_position(n,IZ)

                     
                     ObjectsList(i)% IntegrationVertices(j)% Q   = ObjectsList(i)% IntegrationVertices(j)% Q   + IBMStencilPoints(domain)% lj(k) * Q(index,:)
                     ObjectsList(i)% IntegrationVertices(j)% U_x = ObjectsList(i)% IntegrationVertices(j)% U_x + IBMStencilPoints(domain)% lj(k) * U_x(index,:)
                     ObjectsList(i)% IntegrationVertices(j)% U_y = ObjectsList(i)% IntegrationVertices(j)% U_y + IBMStencilPoints(domain)% lj(k) * U_y(index,:)
                     ObjectsList(i)% IntegrationVertices(j)% U_z = ObjectsList(i)% IntegrationVertices(j)% U_y + IBMStencilPoints(domain)% lj(k) * U_z(index,:)
                  end if 
               end do
            end do
         else
            do n = 1, IBMStencilPoints(domain)% NumOfObjs
               do domains = 1, MPI_Process% nProcs 
                  index = (domains-1)*IBMStencilPoints(domain)% NumOfObjs + n 
                  if( eIDsdomain(index) .ne. 0 ) then 
                     i = IBMStencilPoints(domain)% local_position(n,IX)
                     j = IBMStencilPoints(domain)% local_position(n,IY)

                     ObjectsList(i)% IntegrationVertices(j)% Q   = Q(index,:)
                     ObjectsList(i)% IntegrationVertices(j)% U_x = U_x(index,:)
                     ObjectsList(i)% IntegrationVertices(j)% U_y = U_y(index,:)
                     ObjectsList(i)% IntegrationVertices(j)% U_z = U_z(index,:)
                  end if 
               end do
            end do
         end if 

         deallocate( Q, U_x, U_y, U_z, eIDsdomain )
#endif
      end subroutine GatherHOIntegrationPointsState
!________________________________________________________________________________________________________________________________

   subroutine recvScalarPlotRoot( stl, STLNum, x, y, z, scalar )

      implicit none

      type(STLfile),              intent(in)    :: stl
      integer,                    intent(in)    :: STLNum
      real(kind=RP), allocatable, intent(inout) :: x(:), y(:), z(:), scalar(:)
      !-local-variables-------------------------------------------------------------------
      real(kind=RP)              :: coords(NDIM)
      integer                    :: i, j, k, NumOfObjs       
#ifdef _HAS_MPI_      
      integer                    :: ierr, nProcs, start_index, final_index, &
                                    array_of_statuses(MPI_STATUS_SIZE,4)
      integer,       allocatable :: recv_req(:,:), recv_Firstreq(:), NumOfObjspartition(:)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state(:)
#endif
      if( .not. MPI_Process% isRoot ) return

      NumOfObjs = NumOfVertices * stl% NumOfObjs
#ifdef _HAS_MPI_ 
      allocate( NumOfObjspartition(MPI_Process% nProcs-1), &
                recv_Firstreq(MPI_Process% nProcs-1)       )

      do nProcs = 2, MPI_Process% nProcs
         
         call mpi_irecv( NumOfObjspartition(nProcs-1), 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_Firstreq(nProcs-1), ierr ) 
         call mpi_wait( recv_Firstreq(nProcs-1), MPI_STATUS_IGNORE, ierr )

      end do 

      deallocate(recv_Firstreq)

      allocate( x(NumOfObjs+sum(NumOfObjspartition)),      &
                y(NumOfObjs+sum(NumOfObjspartition)),      &
                z(NumOfObjs+sum(NumOfObjspartition)),      &
                scalar(NumOfObjs+sum(NumOfObjspartition)), &
                recv_req(MPI_Process% nProcs-1,4)          )

      start_index = NumOfObjs
 
      do nProcs = 2, MPI_Process% nProcs  

         allocate( COORD_x(NumOfObjspartition(nProcs-1)), &
                   COORD_y(NumOfObjspartition(nProcs-1)), & 
                   COORD_z(NumOfObjspartition(nProcs-1)), & 
                   state(NumOfObjspartition(nProcs-1))    )

         call mpi_irecv( COORD_x, NumOfObjspartition(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )

         call mpi_irecv( COORD_y, NumOfObjspartition(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )
        
         call mpi_irecv( COORD_z, NumOfObjspartition(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )
        
         call mpi_irecv( state,   NumOfObjspartition(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )

         call mpi_waitall(4, recv_req(nProcs-1,:), array_of_statuses, ierr)  

         do i = 1, NumOfObjspartition(nProcs-1)
            x(start_index+i)      = COORD_x(i)
            y(start_index+i)      = COORD_y(i)
            z(start_index+i)      = COORD_z(i)
            scalar(start_index+i) = state(i)
         end do

         start_index = start_index + NumOfObjspartition(nProcs-1)

         deallocate(COORD_x, COORD_y, COORD_z, state)

      end do

      deallocate( NumOfObjspartition, recv_req )
#else
      allocate( x(NumOfObjs),     &
                y(NumOfObjs),     &
                z(NumOfObjs),     &
                scalar(NumOfObjs) )
#endif  
      k = 0
      do i = 1, stl% NumOfObjs
         associate(obj => stl% ObjectsList(i))
         do j = 1, NumOfVertices
            call OBB(STLNum)% ChangeRefFrame(obj% vertices(j)% coords, GLOBAL, coords)
            k         = k + 1
            x(k)      = coords(1)
            y(k)      = coords(2)
            z(k)      = coords(3)
            scalar(k) = obj% IntegrationVertices(j)% ScalarValue
         end do
         end associate
      end do

   end subroutine recvScalarPlotRoot

   subroutine sendScalarPlotRoot( stl, STLNum )

      implicit none 
      !-arguments-------------------------------------------------------------------------
      type(STLfile), intent(in) :: stl
      integer,       intent(in) :: STLNum    
#ifdef _HAS_MPI_      
      !-local-variables-------------------------------------------------------------------
      integer                    :: ierr, i, j, k, NumOfObjs, status(MPI_STATUS_SIZE), &
                                    array_of_statuses(MPI_STATUS_SIZE,4),              &
                                    send_Firstreq, send_req(4)
      real(kind=RP)              :: coords(NDIM)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state(:)

      if( MPI_Process% isRoot ) return

      NumOfObjs = NumOfVertices * stl% NumOfObjs

      call mpi_isend( NumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_Firstreq, ierr )

      call mpi_wait(send_Firstreq, status, ierr)

      allocate( COORD_x(NumOfObjs), &
                COORD_y(NumOfObjs), &
                COORD_z(NumOfObjs), &
                state(NumOfObjs)    )
                  
      k = 0
      do i = 1, stl% NumOfObjs
         associate(obj => stl% ObjectsList(i))
         do j = 1, NumOfVertices
            call OBB(STLNum)% ChangeRefFrame(obj% vertices(j)% coords, GLOBAL, coords)
            k          = k + 1
            COORD_x(k) = coords(1)
            COORD_y(k) = coords(2)
            COORD_z(k) = coords(3)
            state(k)   = obj% IntegrationVertices(j)% ScalarValue
         end do
         end associate
      end do

      call mpi_isend( COORD_x, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )
      
      call mpi_isend( COORD_y, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )
      
      call mpi_isend( COORD_z, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )
                         
      call mpi_isend( state, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )

      call mpi_waitall( 4, send_req, array_of_statuses, ierr )                            

      deallocate( COORD_x, COORD_y, COORD_z, state )     
#endif
   end subroutine sendScalarPlotRoot

   subroutine recvVectorPlotRoot( stl, STLNum, x, y, z, vector_x, vector_y, vector_z )

      implicit none
      type(STLfile),              intent(in)    :: stl
      integer,                    intent(in)    :: STLNum
      real(kind=RP), allocatable, intent(inout) :: x(:), y(:), z(:), vector_x(:), &
                                                   vector_y(:), vector_z(:)
      !-local-variables-------------------------------------------------------------------
      real(kind=RP)              :: coords(NDIM)
      integer                    :: i, j, k, NumOfObjs       
#ifdef _HAS_MPI_      
      integer                    :: ierr, nProcs, start_index, final_index, &
                                    array_of_statuses(MPI_STATUS_SIZE,6)
      integer,       allocatable :: recv_req(:,:), recv_Firstreq(:), NumOfObjsPartion(:)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state_x(:),   &
                                    state_y(:), state_z(:)
#endif
      if( .not. MPI_Process% isRoot ) return

      NumOfObjs = NumOfVertices * stl% NumOfObjs
#ifdef _HAS_MPI_ 
      allocate( NumOfObjsPartion(MPI_Process% nProcs-1), &
                recv_Firstreq(MPI_Process% nProcs-1)     )

      NumOfObjsPartion = 0

      do nProcs = 2, MPI_Process% nProcs
         
         call mpi_irecv( NumOfObjsPartion(nProcs-1), 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_Firstreq(nProcs-1), ierr ) 
         call mpi_wait( recv_Firstreq(nProcs-1), MPI_STATUS_IGNORE, ierr )

      end do 

      deallocate(recv_Firstreq)
  
      allocate( x(NumOfObjs+sum(NumOfObjsPartion)),        &
                y(NumOfObjs+sum(NumOfObjsPartion)),        &
                z(NumOfObjs+sum(NumOfObjsPartion)),        &
                vector_x(NumOfObjs+sum(NumOfObjsPartion)), &
                vector_y(NumOfObjs+sum(NumOfObjsPartion)), &
                vector_z(NumOfObjs+sum(NumOfObjsPartion)), &
                recv_req(MPI_Process% nProcs-1,6)          )

      start_index = NumOfObjs

      do nProcs = 2, MPI_Process% nProcs  

         allocate( COORD_x(NumOfObjsPartion(nProcs-1)), &
                   COORD_y(NumOfObjsPartion(nProcs-1)), &
                   COORD_z(NumOfObjsPartion(nProcs-1)), &
                   state_x(NumOfObjsPartion(nProcs-1)), &
                   state_y(NumOfObjsPartion(nProcs-1)), &
                   state_z(NumOfObjsPartion(nProcs-1))  )

         call mpi_irecv( COORD_x, NumOfObjsPartion(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )
         
         call mpi_irecv( COORD_y, NumOfObjsPartion(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )
        
         call mpi_irecv( COORD_z, NumOfObjsPartion(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )
        
         call mpi_irecv( state_x, NumOfObjsPartion(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )
        
         call mpi_irecv( state_y, NumOfObjsPartion(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,5), ierr )

         call mpi_irecv( state_z, NumOfObjsPartion(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,6), ierr )

         call mpi_waitall(6, recv_req(nProcs-1,:), array_of_statuses, ierr)    

         do i = 1, NumOfObjsPartion(nProcs-1)
            x(start_index+i)        = COORD_x(i)
            y(start_index+i)        = COORD_y(i)
            z(start_index+i)        = COORD_z(i)
            vector_x(start_index+i) = state_x(i)
            vector_y(start_index+i) = state_y(i)
            vector_z(start_index+i) = state_z(i)
         end do
         start_index = start_index + NumOfObjsPartion(nProcs-1)

         deallocate(COORD_x, COORD_y, COORD_z, state_x, state_y, state_z)
      end do

      deallocate( NumOfObjsPartion, recv_req )
#else
      allocate( x(NumOfObjs),        &
                y(NumOfObjs),        &
                z(NumOfObjs),        &
                vector_x(NumOfObjs), &
                vector_y(NumOfObjs), &
                vector_z(NumOfObjs)  )
#endif  
      k = 0
      do i = 1, stl% NumOfObjs
         associate(obj => stl% ObjectsList(i) )
         do j = 1, NumOfVertices
            call OBB(STLNum)% ChangeRefFrame(obj% vertices(j)% coords,GLOBAL,coords)
            k           = k + 1
            x(k)        = coords(1)
            y(k)        = coords(2)
            z(k)        = coords(3)
            vector_x(k) = obj% IntegrationVertices(j)% VectorValue(IX) 
            vector_y(k) = obj% IntegrationVertices(j)% VectorValue(IY)
            vector_z(k) = obj% IntegrationVertices(j)% VectorValue(IZ)
         end do
         end associate 
      end do

   end subroutine recvVectorPlotRoot

   subroutine sendVectorPlotRoot( stl, STLNum )
      implicit none 

       !-arguments-------------------------------------------------------------------------
      type(STLfile), intent(in) :: stl
      integer,       intent(in) :: STLNum    
#ifdef _HAS_MPI_      
      !-local-variables-------------------------------------------------------------------
      integer                    :: ierr, i, j, k, NumOfObjs, status(MPI_STATUS_SIZE), &
                                    array_of_statuses(MPI_STATUS_SIZE,6),              &
                                    send_Firstreq, send_req(6)
      real(kind=RP)              :: coords(NDIM)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state_x(:),    &
                                    state_y(:), state_z(:)

      if( MPI_Process% isRoot ) return

      NumOfObjs = NumOfVertices * stl% NumOfObjs

      call mpi_isend( NumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_Firstreq, ierr )

      call mpi_wait(send_Firstreq, status, ierr)

      allocate( COORD_x(NumOfObjs), &
                COORD_y(NumOfObjs), &
                COORD_z(NumOfObjs), &
                state_x(NumOfObjs), &
                state_y(NumOfObjs), &
                state_z(NumOfObjs)  )
                  
      k = 0
      do i = 1, stl% NumOfObjs
         associate( obj => stl% ObjectsList(i))
         do j = 1, NumOfVertices
            call OBB(STLNum)% ChangeRefFrame(obj% vertices(j)% coords,GLOBAL,coords)
            k          = k + 1
            COORD_x(k) = coords(1)
            COORD_y(k) = coords(2)
            COORD_z(k) = coords(3)
            state_x(k) = obj% IntegrationVertices(j)% VectorValue(IX)
            state_y(k) = obj% IntegrationVertices(j)% VectorValue(IY)
            state_z(k) = obj% IntegrationVertices(j)% VectorValue(IZ)
         end do
         end associate
      end do

      call mpi_isend( COORD_x, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )
      
      call mpi_isend( COORD_y, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )
      
      call mpi_isend( COORD_z, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )
                         
      call mpi_isend( state_x, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )
                         
      call mpi_isend( state_y, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(5), ierr )

      call mpi_isend( state_z, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(6), ierr )

      call mpi_waitall( 6, send_req, array_of_statuses, ierr )                            

      deallocate( COORD_x, COORD_y, COORD_z, state_x, state_y, state_z )     
#endif
   end subroutine sendVectorPlotRoot

end module MPI_IBMUtilities
