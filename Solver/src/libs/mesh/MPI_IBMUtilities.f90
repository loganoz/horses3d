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
   public :: RecvPointsListRoot, SendPointsList2Root
   public :: RecvPointsListPartitions, SendPointsList2partitions 
   public :: recvScalarPlotRoot, sendScalarPlotRoot
   public :: recvVectorPlotRoot, sendVectorPlotRoot
   public :: MPIProcedures_IBM_HO_faces, IBM_HO_findElements
   public :: MPIStateProcedures_IBM_HO_faces, Set_IBM_HO_faces
   public :: plotSTL, SendOBB, RecvOBB, sendMask, recvMask, SendAxis, RecvAxis
   public :: sendMaskInsideBody, recvMaskInsideBody
   public :: sendMaskGeom, recvMaskGeom
   public :: sendStateBandRegion, recvStateBandRegion
   public :: sendGradientsBandRegion, recvGradientsBandRegion
   
   type IBMpoints

      type(point_type), allocatable :: x(:)
      integer                       :: LocNumOfObjs, NumOfObjs

   end type

   type Intersections_t

      integer, allocatable :: intersections(:)

   end type Intersections_t

   type(IBMpoints), public :: Mask

contains

   subroutine sendStateBandRegion( IBMmask, domain, nEqn )

      implicit none 

      type(IBMpoints), intent(in) :: IBMmask
      integer,         intent(in) :: domain, nEqn
#ifdef _HAS_MPI_ 
      real(kind=RP), allocatable :: Q(:,:)
      integer                    :: domains, i, send_req(nEqn), ierr, & 
                                    array_of_statuses(MPI_STATUS_SIZE,nEqn)

      if( .not. MPI_Process% doMPIAction ) return 

      do domains = 1, MPI_Process% nProcs
         
         if( domains .eq. domain ) cycle 

         if( IBMmask% NumOfObjs .eq. 0 ) cycle 

         allocate( Q(IBMmask% NumOfObjs,nEqn) )

         do i = 1, IBMmask% NumOfObjs
            Q(i,1:nEqn) = IBMmask% x(i)% Q(1:nEqn) 
         end do 

         do i = 1, nEqn 
            call mpi_isend( Q(:,i), IBMmask% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(i), ierr ) 
         end do

         call mpi_waitall( nEqn, send_req, array_of_statuses, ierr)

         deallocate( Q )

      end do 
#endif 
   end subroutine sendStateBandRegion 

   subroutine recvStateBandRegion( IBMmask, nEqn )

      implicit none 

      type(IBMpoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: nEqn 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: Q(:,:)
      integer                    :: domains, domain, i, ierr, recv_req(nEqn), & 
                                    array_of_statuses(MPI_STATUS_SIZE,nEqn)

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle 

         allocate(Q(IBMmask(domains)% NumOfObjs,nEqn))

         do i = 1, nEqn
            call mpi_irecv(Q(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(i), ierr)
         end do
         
         call mpi_waitall( nEqn, recv_req(:), array_of_statuses, ierr )

         do i = 1, IBMmask(domains)% NumOfObjs
            IBMmask(domains)% x(i)% Q(1:nEqn) = Q(i,1:nEqn)
         end do 

         deallocate(Q)

      end do 
#endif
   end subroutine recvStateBandRegion


   subroutine sendGradientsBandRegion( IBMmask, domain, nEqn )

      implicit none 

      type(IBMpoints), intent(in) :: IBMmask
      integer,         intent(in) :: domain, nEqn
#ifdef _HAS_MPI_ 
      real(kind=RP), allocatable :: U_x(:,:), U_y(:,:), U_z(:,:)
      integer                    :: domains, i, send_req(3*nEqn), ierr,       & 
                                    array_of_statuses(MPI_STATUS_SIZE,3*nEqn)

      if( .not. MPI_Process% doMPIAction ) return 

      do domains = 1, MPI_Process% nProcs
         
         if( domains .eq. domain ) cycle 

         if( IBMmask% NumOfObjs .eq. 0 ) cycle 

         allocate( U_x(IBMmask% NumOfObjs,nEqn), &
                   U_y(IBMmask% NumOfObjs,nEqn), &
                   U_z(IBMmask% NumOfObjs,nEqn)  )

         do i = 1, IBMmask% NumOfObjs
            U_x(i,:) = IBMmask% x(i)% U_x
            U_y(i,:) = IBMmask% x(i)% U_y
            U_z(i,:) = IBMmask% x(i)% U_z
         end do 

         do i = 1, nEqn 
            call mpi_isend( U_x(:,i), IBMmask% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(i), ierr ) 
         end do

         do i = 1, nEqn 
            call mpi_isend( U_y(:,i), IBMmask% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nEqn+i), ierr ) 
         end do

         do i = 1, nEqn 
            call mpi_isend( U_z(:,i), IBMmask% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2*nEqn+i), ierr ) 
         end do

         call mpi_waitall( 3*nEqn, send_req, array_of_statuses, ierr)

         deallocate( U_x, U_y, U_z )

      end do 
#endif 
   end subroutine sendGradientsBandRegion 

   subroutine recvGradientsBandRegion( IBMmask, nEqn )

      implicit none 

      type(IBMpoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: nEqn 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: U_x(:,:), U_y(:,:), U_z(:,:)
      integer                    :: domains, domain, i, ierr, recv_req(3*nEqn), & 
                                    array_of_statuses(MPI_STATUS_SIZE,3*nEqn)

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle 

         allocate( U_x(IBMmask(domains)% NumOfObjs,nEqn), &
                   U_y(IBMmask(domains)% NumOfObjs,nEqn), &
                   U_z(IBMmask(domains)% NumOfObjs,nEqn)  )

         do i = 1, nEqn
            call mpi_irecv(U_x(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(i), ierr)
         end do

         do i = 1, nEqn
            call mpi_irecv(U_y(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nEqn+i), ierr)
         end do

         do i = 1, nEqn
            call mpi_irecv(U_z(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2*nEqn+i), ierr)
         end do
         
         call mpi_waitall(3*nEqn, recv_req(:), array_of_statuses, ierr )

         do i = 1, IBMmask(domains)% NumOfObjs
            IBMmask(domains)% x(i)% U_x = U_x(i,:)
            IBMmask(domains)% x(i)% U_y = U_y(i,:)
            IBMmask(domains)% x(i)% U_z = U_z(i,:)
         end do 

         deallocate( U_x, U_y, U_z )

      end do 
#endif
   end subroutine recvGradientsBandRegion
!__________________________________________________
   subroutine sendMaskGeom( IBMmask, domain )

      implicit none 

      type(IBMPoints), intent(in) :: IBMmask(:)
      integer,         intent(in) :: domain 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: dist(:), normal(:,:)
      integer                    :: i, j, domains, send_req(4),           &
                                    array_of_statuses(MPI_STATUS_SIZE,4), &
                                    ierr

      if ( .not. MPI_Process % doMPIAction ) return

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle

         if( IBMmask(domains)% NumOfObjs .eq. 0 ) cycle 

         allocate( dist(IBMmask(domains)% NumOfObjs),       & 
                   normal(IBMmask(domains)% NumOfObjs,NDIM) )

         do i = 1, IBMmask(domains)% NumOfObjs
            dist(i)     = IBMmask(domains)% x(i)% dist
            normal(i,:) = IBMmask(domains)% x(i)% normal
         end do

         call mpi_isend( dist, IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr ) 

         do i = 1, NDIM 
            call mpi_isend( normal(:,i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1+i), ierr ) 
         end do

         call mpi_waitall( 4, send_req, array_of_statuses, ierr) 

         deallocate(dist,normal)

      end do 
#endif
   end subroutine sendMaskGeom

   subroutine recvMaskGeom( IBMmask )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: dist(:), normal(:,:) 
      integer                    :: i, j, domains, domain, recv_req(4),   &
                                    array_of_statuses(MPI_STATUS_SIZE,4), &
                                    ierr
      
      if ( .not. MPI_Process % doMPIAction ) return

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle

         if( IBMmask(domain)% NumOfObjs .eq. 0 ) cycle

         allocate( dist(IBMmask(domain)% NumOfObjs),       &
                   normal(IBMmask(domain)% NumOfObjs,NDIM) )

         call mpi_irecv(dist, IBMmask(domain)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr)

         do i = 1, NDIM
            call mpi_irecv(normal(:,i), IBMmask(domain)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1+i), ierr)
         end do
         
         call mpi_waitall( 4, recv_req(:), array_of_statuses, ierr )

         do i = 1, IBMmask(domain)% NumOfObjs
            if( dist(i) .lt. IBMmask(domain)% x(i)% dist  ) then 
               IBMmask(domain)% x(i)% dist   = dist(i)
               IBMmask(domain)% x(i)% normal = normal(i,:) 
            end if 
         end do

         deallocate(dist,normal)

      end do 
#endif
   end subroutine recvMaskGeom

   subroutine sendMaskInsideBody( IBMmask, domain )

      implicit none 

      type(IBMPoints), intent(in) :: IBMmask(:)
      integer,         intent(in) :: domain 
#ifdef _HAS_MPI_
      integer, allocatable :: intersections(:)
      integer              :: i, domains, domainsSend, ierr, send_req(1)

      if ( .not. MPI_Process % doMPIAction ) return

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle

         if( IBMmask(domains)% NumOfObjs .eq. 0 ) cycle 

         allocate( intersections(IBMmask(domains)% NumOfObjs) )
   
         do i = 1, IBMmask(domains)% NumOfObjs
            intersections(i) = IBMmask(domains)% x(i)% NumOfIntersections
         end do

         call mpi_isend( intersections, IBMmask(domains)% NumOfObjs, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr ) 

         call mpi_wait( send_req(1), MPI_STATUS_IGNORE, ierr )

         deallocate(intersections)

      end do 
#endif
   end subroutine sendMaskInsideBody

   subroutine recvMaskInsideBody( IBMmask )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      integer,               allocatable :: intersections(:)
      integer                            :: i, domains, domain, domainsRecv, ierr, recv_req(1)
      
      if ( .not. MPI_Process % doMPIAction ) return

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle

         if( IBMmask(domain)% NumOfObjs .eq. 0 ) cycle

         allocate( intersections(IBMmask(domain)% NumOfObjs) ) 

         call mpi_irecv(intersections, IBMmask(domain)% NumOfObjs, MPI_INT, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr)

         call mpi_wait( recv_req(1), MPI_STATUS_IGNORE, ierr ) 

         do i = 1, IBMmask(domain)% NumOfObjs
            IBMmask(domain)% x(i)% NumOfIntersections = IBMmask(domain)% x(i)% NumOfIntersections + intersections(i)
         end do
         
         deallocate(intersections)

      end do
#endif
   end subroutine recvMaskInsideBody 

   subroutine sendMask( IBMmask, domain )

      implicit none 

      type(IBMPoints), intent(in) :: IBMmask
      integer,         intent(in) :: domain 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: coords(:,:)
      integer,       allocatable :: eIDs(:), local_Position(:,:) , send_req(:,:)
      integer                    :: i, domains, ierr, array_of_statuses(MPI_STATUS_SIZE,8)

      if ( .not. MPI_Process % doMPIAction ) return

      allocate( coords(IBMmask% NumOfObjs,NDIM),         &
                eIDs(IBMmask% NumOfObjs),                &
                local_Position(IBMmask% NumOfObjs,NDIM), &
                send_req(MPI_Process% nProcs,8)          )

      do i = 1, IBMmask% NumOfObjs
         coords(i,:)         = IBMmask% x(i)% coords
         eIDs(i)             = IBMmask% x(i)% element_index
         local_Position(i,:) = IBMmask% x(i)% local_Position
      end do
 
      do domains = 1, MPI_Process% nProcs 

         if( domains .eq. domain ) cycle

         call mpi_isend( IBMmask% NumOfObjs, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,1), ierr ) 

         call mpi_wait( send_req(domains,1), MPI_STATUS_IGNORE, ierr )

         if( IBMmask% NumOfObjs .eq. 0 ) cycle 

         do i = 1, NDIM 
            call mpi_isend(coords(:,i), IBMmask% NumOfObjs, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,1+i), ierr )
         end do

         call mpi_isend(eIDs, IBMmask% NumOfObjs, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,5), ierr )

         do i = 1, NDIM 
            call mpi_isend(local_Position(:,i), IBMmask% NumOfObjs, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,5+i), ierr )
         end do

         call mpi_waitall(8, send_req(domains,:), array_of_statuses, ierr)

      end do 

      deallocate(coords,eIDs,local_Position,send_req)
#endif
   end subroutine sendMask

   subroutine recvMask( IBMmask )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: coords(:,:)
      integer,       allocatable :: eIDs(:), local_Position(:,:), recv_req(:,:)
      integer                    :: domain, i, domains, NumOfObjs, ierr, array_of_statuses(MPI_STATUS_SIZE,8)

      domain = MPI_Process% rank + 1

      allocate(recv_req(MPI_Process% nProcs,8))

      do domains = 1, MPI_Process% nProcs

         if( domains .eq. domain ) cycle

         call mpi_irecv(NumOfObjs, 1, MPI_INT, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,1), ierr)

         call mpi_wait( recv_req(domains,1), MPI_STATUS_IGNORE, ierr )

         IBMmask(domains)% NumOfObjs = NumOfObjs

         if( NumOfObjs .eq. 0 ) cycle 

         allocate( coords(NumOfObjs,NDIM),        &
                   eIDs(NumOfObjs),               &
                   local_Position(NumOfObjs,NDIM) )

         do i = 1, NDIM 
            call mpi_irecv(coords(:,i), NumOfObjs, MPI_DOUBLE, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,1+i), ierr)   
         end do 
         
         call mpi_irecv(eIDs, NumOfObjs, MPI_INT, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,5), ierr)

         do i = 1, NDIM 
            call mpi_irecv(local_position(:,i), NumOfObjs, MPI_INT, domains-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(domains,5+i), ierr)
         end do 

         call mpi_waitall(8, recv_req(domains,:), array_of_statuses, ierr)

         allocate( IBMmask(domains)% x(NumOfObjs) )

         do i = 1, NumOfObjs
            IBMmask(domains)% x(i)% coords         = coords(i,:) 
            IBMmask(domains)% x(i)% element_index  = eIDs(i) 
            IBMmask(domains)% x(i)% local_position = local_position(i,:) 
         end do

         deallocate(coords, eIDs, local_position)
      end do

      deallocate(recv_req)
#endif
   end subroutine recvMask

   subroutine SendAxis( STLNum ) 

      implicit none 

      integer, intent(in) :: STLNum 
#ifdef _HAS_MPI_
      integer, allocatable :: send_req(:,:)
      integer              :: domains, ierr, array_of_statuses(MPI_STATUS_SIZE,2)

      if( .not. MPI_Process% isRoot ) return 

      allocate(send_req(MPI_Process% nProcs,2))
      
      do domains = 2, MPI_Process% nProcs

         call mpi_isend(OBB(STLNum)% maxAxis, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,1), ierr )

         call mpi_isend(OBB(STLNum)% minAxis, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,2), ierr )

         call mpi_waitall( 2, send_req(domains,:), array_of_statuses, ierr )
      end do

      deallocate(send_req)
#endif
   end subroutine SendAxis

   subroutine RecvAxis( STLNum )

      implicit none 

      integer, intent(in) :: STLNum
#ifdef _HAS_MPI_
      integer :: recv_req(2), ierr, array_of_statuses(MPI_STATUS_SIZE,2)
      
      if( MPI_Process% isRoot ) return 

      call mpi_irecv(OBB(STLNum)% MaxAxis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr)

      call mpi_irecv(OBB(STLNum)% minAxis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr)

      call mpi_waitall( 2, recv_req, array_of_statuses, ierr )
#endif
   end subroutine RecvAxis


   subroutine SendOBB( STLNum )

      implicit none 

      integer, intent(in) :: STLNum 
#ifdef _HAS_MPI_
      integer, allocatable :: send_req(:,:)
      integer              :: domains, i, ierr, domain, &
                              array_of_statuses(MPI_STATUS_SIZE,17)

      if ( .not. MPI_Process% isRoot ) return

      domain = MPI_Process% rank+1

      allocate(send_req(MPI_Process% nProcs,17))

      do domains = 2, MPI_Process% nProcs

         if( domains .eq. domain ) cycle

         call mpi_isend(OBB(STLNum)% MBR% angle, 1, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,1), ierr )

         call mpi_isend(OBB(STLNum)% MBR% center, NDIM-1, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,2), ierr )
         
         call mpi_isend(OBB(STLNum)% CloudCenter, NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,3), ierr )
     
         do i = 1, NDIM 
            call mpi_isend(OBB(STLNum)% R(:,i), NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,3+i), ierr )
         end do 

         do i = 1, NDIM 
            call mpi_isend(OBB(STLNum)% invR(:,i), NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,6+i), ierr )
         end do 

         do i = 1, BOXVERTICES
            call mpi_isend(OBB(STLNum)% LocVertices(:,i), NDIM, MPI_DOUBLE, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains,9+i), ierr )
         end do 
 
         call mpi_waitall(17, send_req(domains,:), array_of_statuses, ierr)
      
      end do 

      deallocate(send_req)
#endif
   end subroutine SendOBB
   
   subroutine recvOBB( STLNum )

      implicit none 

      integer, intent(in) :: STLNum

      integer :: domain, i
#ifdef _HAS_MPI_
      integer :: domains, ierr, array_of_statuses(MPI_STATUS_SIZE,17), &
                 recv_req(17)
#endif
      if( MPI_Process% isRoot ) return 
#ifdef _HAS_MPI_
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

      call mpi_waitall(17, recv_req, array_of_statuses, ierr)
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
      integer                        :: NumOfObjs, NumOfObjsPartion,           &
                                        start_index, final_index, i, j,        &
                                        nProcs, ierr, biggerdomains,           &
                                        elems_per_domain(MPI_Process% nProcs), &
                                        array_of_statuses(MPI_STATUS_SIZE,13)
      integer, allocatable           :: SortedIndex(:), send_req(:,:)

      STL% partition = MPI_Process% rank

      NumOfObjs = STL% NumOfObjs

      allocate( Bar(NumOfObjs),                          &
                coord(NumOfObjs),                        &
                SortedIndex(NumOfObjs),                  &
                normals_x(NumOfObjs),                    &
                normals_y(NumOfObjs),                    &
                normals_z(NumOfObjs),                    &
                vertices_x(NumOfObjs,NumOfVertices),     &
                vertices_y(NumOfObjs,NumOfVertices),     &
                vertices_z(NumOfObjs,NumOfVertices),     &
                send_req(MPI_Process% nProcs-1,13)       )

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

         start_index = final_index + 1

      end do

      start_index = elems_per_domain(1)+1
      do nProcs = 2, MPI_Process% nProcs

         final_index = start_index + elems_per_domain(nProcs) - 1
            
         NumOfObjsPartion = final_index-start_index+1

         call mpi_isend(NumOfObjsPartion, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )

         call mpi_wait(send_req(nProcs-1,1),MPI_STATUS_IGNORE,ierr)

         call mpi_isend(normals_x(start_index:final_index), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )

         call mpi_isend(normals_y(start_index:final_index), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )

         call mpi_isend(normals_z(start_index:final_index), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )

         call mpi_isend(vertices_x(start_index:final_index,1), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,5), ierr )

         call mpi_isend(vertices_y(start_index:final_index,1), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,6), ierr )

         call mpi_isend(vertices_z(start_index:final_index,1), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,7), ierr )

         call mpi_isend(vertices_x(start_index:final_index,2), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,8), ierr )

         call mpi_isend(vertices_y(start_index:final_index,2), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,9), ierr )

         call mpi_isend(vertices_z(start_index:final_index,2), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,10), ierr )

         call mpi_isend(vertices_x(start_index:final_index,3), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,11), ierr )

         call mpi_isend(vertices_y(start_index:final_index,3), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,12), ierr )

         call mpi_isend(vertices_z(start_index:final_index,3), NumOfObjsPartion, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,13), ierr )

         call mpi_waitall(13, send_req(nProcs-1,:), array_of_statuses, ierr)

         start_index = final_index + 1

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

      deallocate(send_req, vertices_x, vertices_y, vertices_z, normals_x, normals_y, normals_z)
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
      integer                    :: NumOfObjs, ierr, recv_req(13), i, j,      &
                                    array_of_statuses(MPI_STATUS_SIZE,13)

      if( MPI_Process% isRoot ) return

      STL% partition = MPI_Process% rank 

      call mpi_irecv( STL% NumOfObjs, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )

      call mpi_wait(recv_req(1), MPI_STATUS_IGNORE, ierr)

      NumOfObjs = STL% NumOfObjs

      allocate( normals_x(NumOfObjs),       &
                normals_y(NumOfObjs),       &
                normals_z(NumOfObjs),       &
                vertices_x(NumOfObjs,NDIM), &
                vertices_y(NumOfObjs,NDIM), &
                vertices_z(NumOfObjs,NDIM)  )

      call mpi_irecv( normals_x, NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr )  

      call mpi_irecv( normals_y, NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr )  
                      
      call mpi_irecv( normals_z, NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr )  
                      
      call mpi_irecv( vertices_x(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr )  

      call mpi_irecv( vertices_y(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr )
                      
      call mpi_irecv( vertices_z(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr )  

      call mpi_irecv( vertices_x(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr )

      call mpi_irecv( vertices_y(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9), ierr )

      call mpi_irecv( vertices_z(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(10), ierr )

      call mpi_irecv( vertices_x(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(11), ierr )
      
      call mpi_irecv( vertices_y(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(12), ierr )
      
      call mpi_irecv( vertices_z(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(13), ierr )

      call mpi_waitall(13, recv_req, array_of_statuses, ierr)

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
         normals_x(i) = STL% ObjectsList(i)% normal(1) 
         normals_y(i) = STL% ObjectsList(i)% normal(2)  
         normals_z(i) = STL% ObjectsList(i)% normal(3)  
         do j = 1, NumOfVertices
            vertices_x(i,j) = STL% ObjectsList(i)% vertices(j)% coords(1) 
            vertices_y(i,j) = STL% ObjectsList(i)% vertices(j)% coords(2) 
            vertices_z(i,j) = STL% ObjectsList(i)% vertices(j)% coords(3) 
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

      type(STLfile), intent(in) :: STL 
      integer,       intent(in) :: iter
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

      type(STLfile), intent(in) :: STL 
      integer,       intent(in) :: iter

      if( MPI_Process% doMPIAction ) then 
         call sendSTLRoot( STL )
      end if 

      if( MPI_Process% isRoot ) then 
         call recvSTLRootandPlot( STL, iter )
      end if 

   end subroutine plotSTL


   subroutine RecvPointsListpartitions( PointsList )
      implicit none 
      !-arguments-----------------------------------------------------------------
      type(IBMpoints), intent(inout) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables-----------------------------------------------------------
      real(kind=RP), allocatable :: coords(:,:)
      integer,       allocatable :: local_position(:,:), element_index(:),   &
                                    partition(:), indeces(:)
      integer                    :: i, nProcs, ierr, recv_req(9),            &
                                    array_of_statuses(MPI_STATUS_SIZE,9) 

      if( MPI_Process% isRoot ) return

      allocate( coords(PointsList% NumOfObjs,NDIM),         &
                local_position(PointsList% NumOfObjs,NDIM), &
                element_index(PointsList% NumOfObjs),       &
                partition(PointsList% NumOfObjs),           &
                indeces(PointsList% NumOfObjs)              )

      call mpi_irecv( coords(:,1), PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )

      call mpi_irecv( coords(:,2), PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr )

      call mpi_irecv( coords(:,3), PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr )

      call mpi_irecv( local_position(:,1), PointsList% NumOfObjs, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr )

      call mpi_irecv( local_position(:,2), PointsList% NumOfObjs, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr )

      call mpi_irecv( local_position(:,3), PointsList% NumOfObjs, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr )

      call mpi_irecv( element_index, PointsList% NumOfObjs, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr )

      call mpi_irecv( partition, PointsList% NumOfObjs, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr )

      call mpi_irecv( indeces, PointsList% NumOfObjs, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9), ierr )

      call mpi_waitall(9, recv_req, array_of_statuses, ierr)

      do i = 1, PointsList% NumOfObjs
         PointsList% x(i)% coords         = coords(i,:)
         PointsList% x(i)% local_position = local_position(i,:)
         PointsList% x(i)% element_index  = element_index(i)
         PointsList% x(i)% partition      = partition(i)
         PointsList% x(i)% index          = indeces(i)
      end do

      deallocate( coords,         &
                  local_position, &
                  element_index,  &
                  partition,      &
                  indeces         )
#endif
   end subroutine RecvPointsListpartitions


   subroutine SendPointsList2partitions( PointsList )
      implicit none 
      !-arguments----------------------------------------------------------------
      type(IBMPoints), intent(inout) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables----------------------------------------------------------
      real(kind=RP), allocatable :: coords(:,:), normals(:,:), Dist(:)
      integer,       allocatable :: local_position(:,:), element_index(:),   &
                                    partition(:), indeces(:), send_req(:,:)
      integer                    :: i, nProcs, ierr,                         &
                                    array_of_statuses(MPI_STATUS_SIZE,9)

      allocate( coords(PointsList% NumOfObjs,NDIM),         &
                local_position(PointsList% NumOfObjs,NDIM), &
                element_index(PointsList% NumOfObjs),       &
                partition(PointsList% NumOfObjs),           &
                indeces(PointsList% NumOfObjs),             &
                send_req(MPI_Process% nProcs-1,9)          )

      do i = 1, PointsList% NumOfObjs
         coords(i,:)         = PointsList% x(i)% coords
         local_position(i,:) = PointsList% x(i)% local_position
         element_index(i)    = PointsList% x(i)% element_index
         partition(i)        = PointsList% x(i)% partition
         indeces(i)          = PointsList% x(i)% index
      end do 

      do nProcs = 2, MPI_Process% nProcs

         call mpi_isend( coords(:,1), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )

         call mpi_isend( coords(:,2), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )

         call mpi_isend( coords(:,3), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )

         call mpi_isend( local_position(:,1), PointsList% NumOfObjs, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )

         call mpi_isend( local_position(:,2), PointsList% NumOfObjs, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,5), ierr )

         call mpi_isend( local_position(:,3), PointsList% NumOfObjs, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,6), ierr )

         call mpi_isend( element_index, PointsList% NumOfObjs, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,7), ierr )
      
         call mpi_isend( partition, PointsList% NumOfObjs, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,8), ierr )

         call mpi_isend( indeces, PointsList% NumOfObjs, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,9), ierr )

         call mpi_waitall(9, send_req(nProcs-1,:), array_of_statuses, ierr)

      end do

      deallocate( coords,         &
                  local_position, &
                  element_index,  &
                  partition,      &
                  indeces,        &
                  send_req        )
#endif
   end subroutine SendPointsList2partitions

   subroutine RecvPointsListRoot( PointsList )
      implicit none 
      !-arguments-------------------------------------------------------------------------
      type(IBMPoints),           intent(inout) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables-------------------------------------------------------------------
      real(kind=RP), allocatable :: coords(:,:), normals_x(:,:), normals_y(:,:),        &
                                    normals_z(:,:)
      integer,       allocatable :: local_position(:,:),  element_index(:),             &
                                    partition(:), recv_req(:,:)
      integer                    :: i, LocNumOfObjs, start_index, final_index, ierr,    &
                                    rank, nProcs, array_of_statuses(MPI_STATUS_SIZE,9)

      allocate( coords(PointsList% NumOfObjs,NDIM),                     &
                local_position(PointsList% NumOfObjs,NDIM),             &
                element_index(PointsList% NumOfObjs),                   &
                partition(PointsList% NumOfObjs),                       &
                recv_req(MPI_Process% nProcs-1,9)                      ) 

      start_index = PointsList% LocNumOfObjs; final_index = PointsList% LocNumOfObjs

      do nProcs = 2, MPI_Process% nProcs 

         call mpi_irecv( LocNumOfObjs, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )

         call mpi_wait(recv_req(nProcs-1,1), MPI_STATUS_IGNORE, ierr) 

         if( LocNumOfObjs .eq. 0 ) cycle

         call mpi_irecv( coords(1:LocNumOfObjs,1), LocNumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )

         call mpi_irecv( coords(1:LocNumOfObjs,2), LocNumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )

         call mpi_irecv( coords(1:LocNumOfObjs,3), LocNumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )

         call mpi_irecv( local_position(1:LocNumOfObjs,1), LocNumOfObjs, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,5), ierr )

         call mpi_irecv( local_position(1:LocNumOfObjs,2), LocNumOfObjs, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,6), ierr )

         call mpi_irecv( local_position(1:LocNumOfObjs,3), LocNumOfObjs, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,7), ierr )

         call mpi_irecv( element_index(1:LocNumOfObjs), LocNumOfObjs, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,8), ierr )

         call mpi_irecv( partition(1:LocNumOfObjs), LocNumOfObjs, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,9), ierr )

         call mpi_waitall(9, recv_req(nProcs-1,:), array_of_statuses, ierr)

         start_index = final_index
         final_index = start_index + LocNumOfObjs

         do i = 1, LocNumOfObjs
            PointsList% x(start_index+i)% coords         = coords(i,:)
            PointsList% x(start_index+i)% local_position = local_position(i,:)
            PointsList% x(start_index+i)% element_index  = element_index(i)
            PointsList% x(start_index+i)% partition      = partition(i)
            PointsList% x(start_index+i)% index          = start_index + i 
         end do

      end do

      deallocate( coords,         &
                  local_position, &
                  element_index,  &
                  partition,      &
                  recv_req        )
#endif
   end subroutine RecvPointsListRoot


   subroutine sendPointsList2Root( PointsList )
      implicit none 
      !-arguments-----------------------------
      type(IBMPoints), intent(inout) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables--------------
      real(kind=RP), allocatable :: coords(:,:), normals_x(:), normals_y(:), &
                                    normals_z(:)
      integer,       allocatable :: local_position(:,:),  element_index(:),  &
                                    partition(:)
      integer                    :: i, send_req(9), ierr,                   &
                                    array_of_statuses(MPI_STATUS_SIZE,9)

      if( MPI_Process% isRoot ) return

      call mpi_isend( PointsList% LocNumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )

      call mpi_wait(send_req(1), MPI_STATUS_IGNORE, ierr)

      if( PointsList% LocNumOfObjs .eq. 0 ) return

      allocate( coords(PointsList% LocNumOfObjs,NDIM),         &
                local_position(PointsList% LocNumOfObjs,NDIM), &
                element_index(PointsList% LocNumOfObjs),       &
                partition(PointsList% LocNumOfObjs)            )

      do i = 1, PointsList% LocNumOfObjs
         coords(i,:)         = PointsList% x(i)% coords
         local_position(i,:) = PointsList% x(i)% local_position
         element_index(i)    = PointsList% x(i)% element_index
         partition(i)        = PointsList% x(i)% partition
      end do 

      call mpi_isend( coords(:,1), PointsList% LocNumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )

      call mpi_isend( coords(:,2), PointsList% LocNumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )

      call mpi_isend( coords(:,3), PointsList% LocNumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )

      call mpi_isend( local_position(:,1), PointsList% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(5), ierr )

      call mpi_isend( local_position(:,2), PointsList% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(6), ierr )

      call mpi_isend( local_position(:,3), PointsList% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(7), ierr )

      call mpi_isend( element_index, PointsList% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(8), ierr )
      
      call mpi_isend( partition, PointsList% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(9), ierr )

      call mpi_waitall(9, send_req, array_of_statuses, ierr)

      deallocate( coords,         &
                  local_position, &
                  element_index,  &
                  partition       )
#endif
   end subroutine sendPointsList2Root

!____________________________________________________________________________________________________________________________

   subroutine RecvFacesPartitions( )
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables-----------------------------------------------------------
      real(kind=RP), allocatable :: coords(:,:), xiB(:), xiI(:)
      integer,       allocatable :: face_order(:,:), face_index(:)
      integer                    :: i, j, k, fID, m, l, nProcs, ierr, domain,       &
                                    recv_req(11), NumOfFaces, NumOfObjs, NumOfDoFs, &
                                    array_of_statuses(MPI_STATUS_SIZE,11) 

      if( MPI_Process% isRoot ) return
      
      do domain = 0, MPI_Process% nProcs-1

         call mpi_irecv( NumOfFaces, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )

         call mpi_wait(recv_req(1), MPI_STATUS_IGNORE, ierr) 

         IBM_HO_faces(domain)% NumOfFaces = NumOfFaces

         if( NumOfFaces .eq. 0 ) cycle

         call mpi_irecv( NumOfObjs, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr )

         call mpi_wait(recv_req(2), MPI_STATUS_IGNORE, ierr)

         IBM_HO_faces(domain)% NumOfObjs = NumOfObjs

         call mpi_irecv( NumOfDoFs, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr )

         call mpi_wait(recv_req(3), MPI_STATUS_IGNORE, ierr)

         IBM_HO_faces(domain)% NumOfDoFs = NumOfDoFs

         allocate( coords(NumOfObjs,NDIM),        &
                   xiB(NumOfDoFs),                &
                   xiI(NumOfDoFs),                &
                   face_order(NumOfFaces,NDIM-1), &
                   face_index(NumOfFaces)         )

         call mpi_irecv( coords(:,1), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr )

         call mpi_irecv( coords(:,2), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr )

         call mpi_irecv( coords(:,3), NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr )

         call mpi_irecv( xiB, NumOfDoFs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr )

         call mpi_irecv( xiI, NumOfDoFs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr )

         call mpi_irecv( face_order(:,1), NumOfFaces, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9), ierr )

         call mpi_irecv( face_order(:,2), NumOfFaces, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(10), ierr )

         call mpi_irecv( face_index, NumOfFaces, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(11), ierr )

         call mpi_waitall(11, recv_req, array_of_statuses, ierr)

         if( domain .ne. MPI_Process% rank ) then 
            allocate(IBM_HO_faces(domain)% faces(NumOfFaces))
            m = 0; l = 0
            do fID = 1, NumOfFaces
               associate( f => IBM_HO_faces(domain)% faces(fID) )
               f% Nf(1) = face_order(fID,1)
               f% Nf(2) = face_order(fID,2)
               f% ID   = face_index(fID)
               call f% StencilConstruct()
               do i = 0, face_order(fID,1); do j = 0, face_order(fID,2)
                  l = l + 1
                  f% stencil(i,j)% xiB = xiB(l)
                  f% stencil(i,j)% xiI = xiI(l)
                  do k = 0, f% stencil(i,j)% N
                     m = m + 1 
                     f% stencil(i,j)% x_s(:,k) = coords(m,:)
                  end do
               end do; end do 
               end associate
            end do
         end if 

         deallocate( coords,     &
                     xiB,        &
                     xiI,        &
                     face_order, &
                     face_index  )

      end do
#endif
   end subroutine RecvFacesPartitions

   subroutine SendFaces2Partitions( )
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables----------------------------------------------------------
      real(kind=RP), allocatable :: coords(:,:), xiB(:), xiI(:)
      integer,       allocatable :: face_order(:,:), face_index(:),           &
                                    send_req(:,:)
      integer                    :: i, j, k, m, l, nProcs, ierr, domain, fID, &
                                    NumOfObjs, NumOfFaces, NumOfDoFs,         &
                                    array_of_statuses(MPI_STATUS_SIZE,11)

      allocate(send_req(MPI_Process% nProcs-1,11))

      do nProcs = 2, MPI_Process% nProcs
         do domain = 0, MPI_Process% nProcs-1
            NumOfFaces = IBM_HO_faces(domain)% NumOfFaces
            NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs
            NumOfDoFs  = IBM_HO_faces(domain)% NumOfDoFs

            call mpi_isend( NumOfFaces, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )

            call mpi_wait(send_req(nProcs-1,1), MPI_STATUS_IGNORE, ierr)
      
            if( NumOfFaces .eq. 0 ) cycle
      
            call mpi_isend( NumOfObjs, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )
      
            call mpi_wait(send_req(nProcs-1,2), MPI_STATUS_IGNORE, ierr)

            call mpi_isend( NumOfDoFs, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )
      
            call mpi_wait(send_req(nProcs-1,3), MPI_STATUS_IGNORE, ierr)
 
            allocate( coords(NumOfObjs,NDIM),        &
                      xiB(NumOfDoFs),                &
                      xiI(NumOfDoFs),                &
                      face_order(NumOfFaces,NDIM-1), &
                      face_index(NumOfFaces)         )
            
            m = 0; l = 0
            do fID = 1, NumOfFaces
               associate( f => IBM_HO_faces(domain)% faces(fID) )
               face_index(fID)   = f% ID 
               face_order(fID,1) = f% Nf(1) 
               face_order(fID,2) = f% Nf(2) 
               do i = 0, f% Nf(1); do j = 0, f% Nf(2)
                  l = l + 1
                  xiB(l) = f% stencil(i,j)% xiB 
                  xiI(l) = f% stencil(i,j)% xiI 
                  do k = 0, f% stencil(i,j)% N
                     m = m + 1 
                     coords(m,1) = f% stencil(i,j)% x_s(1,k)
                     coords(m,2) = f% stencil(i,j)% x_s(2,k)
                     coords(m,3) = f% stencil(i,j)% x_s(3,k)
                  end do 
               end do; end do 
               end associate
            end do   

            call mpi_isend( coords(:,1), NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )

            call mpi_isend( coords(:,2), NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,5), ierr )

            call mpi_isend( coords(:,3), NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,6), ierr )

            call mpi_isend( xiB, NumOfDoFs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,7), ierr )

            call mpi_isend( xiI, NumOfDoFs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,8), ierr )

            call mpi_isend( face_order(:,1), NumOfFaces, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,9), ierr )

            call mpi_isend( face_order(:,2), NumOfFaces, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,10), ierr )

            call mpi_isend( face_index, NumOfFaces, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,11), ierr )
      
            call mpi_waitall(11, send_req(nProcs-1,:), array_of_statuses, ierr)

            deallocate( coords,     &
                        xiB,        &
                        xiI,        &
                        face_order, &
                        face_index  )
         end do 
      end do
      
      deallocate( send_req )
#endif
   end subroutine SendFaces2Partitions

   subroutine RecvFacesRoot( )
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables-------------------------------------------------------------------
      real(kind=RP), allocatable :: coords(:,:), xiB(:), xiI(:)
      integer,       allocatable :: face_order(:,:), face_index(:), recv_req(:,:)
      integer                    :: i, j, k, m, l, N, fID, ierr,                         &
                                    rank, nProcs, array_of_statuses(MPI_STATUS_SIZE,11), &
                                    NumOfObjs, NumOfFaces, NumOfDoFs

      allocate( recv_req(MPI_Process% nProcs-1,11) )

      do nProcs = 2, MPI_Process% nProcs 

         call mpi_irecv( NumOfFaces, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )

         call mpi_wait(recv_req(nProcs-1,1), MPI_STATUS_IGNORE, ierr) 

         IBM_HO_faces(nProcs-1)% NumOfFaces = NumOfFaces

         if( NumOfFaces .eq. 0 ) cycle

         call mpi_irecv( NumOfObjs, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )

         call mpi_wait(recv_req(nProcs-1,2), MPI_STATUS_IGNORE, ierr)

         IBM_HO_faces(nProcs-1)% NumOfObjs = NumOfObjs

         call mpi_irecv( NumOfDoFs, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )

         call mpi_wait(recv_req(nProcs-1,3), MPI_STATUS_IGNORE, ierr)

         IBM_HO_faces(nProcs-1)% NumOfDoFs = NumOfDoFs

         allocate( coords(NumOfObjs,NDIM),        &
                   xiB(NumOfDoFs),                &
                   xiI(NumOfDoFs),                &
                   face_order(NumOfFaces,NDIM-1), &
                   face_index(NumOfFaces)         )

         call mpi_irecv( coords(:,1), NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )

         call mpi_irecv( coords(:,2), NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,5), ierr )

         call mpi_irecv( coords(:,3), NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,6), ierr )

         call mpi_irecv( xiB, NumOfDoFs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,7), ierr )

         call mpi_irecv( xiI, NumOfDoFs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,8), ierr )
         
         call mpi_irecv( face_order(:,1), NumOfFaces, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,9), ierr )
         
         call mpi_irecv( face_order(:,2), NumOfFaces, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,10), ierr )

         call mpi_irecv( face_index, NumOfFaces, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,11), ierr )

         call mpi_waitall(11, recv_req(nProcs-1,:), array_of_statuses, ierr)         
         
         allocate( IBM_HO_faces(nProcs-1)% faces(NumOfFaces) )
         m = 0; l = 0
         do fID = 1, NumOfFaces
            associate( f => IBM_HO_faces(nProcs-1)% faces(fID) )
            f% Nf(1) = face_order(fID,1)
            f% Nf(2) = face_order(fID,2)
            f% ID    = face_index(fID)
            call f% StencilConstruct()
            do i = 0, face_order(fID,1); do j = 0, face_order(fID,2) 
               l = l + 1
               f% stencil(i,j)% xiB = xiB(l)
               f% stencil(i,j)% xiI = xiI(l)
               do k = 0, f% stencil(i,j)% N
                  m = m + 1 
                  f% stencil(i,j)% x_s(:,k) = coords(m,:) 
               end do 
            end do; end do
            end associate
         end do 

         deallocate( coords,     &
                     xiB,        &
                     xiI,        &
                     face_order, &
                     face_index  )

      end do
      
      deallocate( recv_req )
#endif
   end subroutine RecvFacesRoot

   subroutine SendFaces2Root()
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables--------------
      real(kind=RP), allocatable :: coords(:,:), xiB(:), xiI(:)
      integer,       allocatable :: face_order(:,:), face_index(:)
      integer                    :: fID, i, j, k, m, l, send_req(11), ierr,  &
                                    array_of_statuses(MPI_STATUS_SIZE,11),   &
                                    domain, NumOfFaces, NumOfObjs, NumOfDoFs

      if( MPI_Process% isRoot ) return

      domain     = MPI_Process% rank   
      NumOfFaces = IBM_HO_faces(domain)% NumOfFaces
      NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs 
      NumOfDoFs  = IBM_HO_faces(domain)% NumOfDoFs 

      call mpi_isend( NumOfFaces, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )
      call mpi_wait( send_req(1), MPI_STATUS_IGNORE, ierr )

      if( NumOfFaces .eq. 0 ) return

      call mpi_isend( NumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )
      call mpi_wait( send_req(2), MPI_STATUS_IGNORE, ierr )

      call mpi_isend( NumOfDoFs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )
      call mpi_wait( send_req(3), MPI_STATUS_IGNORE, ierr )

      allocate( coords(NumOfObjs,NDIM),        &
                xiB(NumOfDoFs),                &
                xiI(NumOfDoFs),                &
                face_order(NumOfFaces,NDIM-1), &
                face_index(NumOfFaces)         )

      m = 0; l = 0
      do fID = 1, NumOfFaces
         associate(f => IBM_HO_faces(domain)% faces(fID))
         face_order(fID,1) = f% Nf(1)
         face_order(fID,2) = f% Nf(2)
         face_index(fID)   = f% ID 
         do i = 0, f% Nf(1); do j = 0, f% Nf(2)
            l = l + 1 
            xiB(l) = f% stencil(i,j)% xiB
            xiI(l) = f% stencil(i,j)% xiI
            do k = 0, f% stencil(i,j)% N
               m = m + 1 
               coords(m,1) = f% stencil(i,j)% x_s(1,k)
               coords(m,2) = f% stencil(i,j)% x_s(2,k)
               coords(m,3) = f% stencil(i,j)% x_s(3,k)
            end do 
         end do; end do 
         end associate
      end do

      call mpi_isend( coords(:,1), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )

      call mpi_isend( coords(:,2), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(5), ierr )

      call mpi_isend( coords(:,3), NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(6), ierr )

      call mpi_isend( xiB, NumOfDoFs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(7), ierr )

      call mpi_isend( xiI, NumOfDoFs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(8), ierr )

      call mpi_isend( face_order(:,1), NumOfFaces, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(9), ierr )

      call mpi_isend( face_order(:,2), NumOfFaces, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(10), ierr )

      call mpi_isend( face_index, NumOfFaces, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(11), ierr )
      
      call mpi_waitall(11, send_req, array_of_statuses, ierr)

      deallocate( coords,     &
                  xiB,        &
                  xiI,        &
                  face_order, &
                  face_index  )
#endif
   end subroutine SendFaces2Root


   subroutine RecvFacesStatePartitions( )
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables-----------------------------------------------------------
      real(kind=RP), allocatable :: u(:), v(:), w(:)
      integer                    :: i, j, k, fID, m, nProcs, ierr, domain, &
                                    recv_req(3), NumOfFaces, NumOfObjs,    &
                                    array_of_statuses(MPI_STATUS_SIZE,3) 

      if( MPI_Process% isRoot ) return
      
      do domain = 0, MPI_Process% nProcs-1

         NumOfFaces = IBM_HO_faces(domain)% NumOfFaces
         NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs

         if( NumOfFaces .eq. 0 ) cycle

         allocate( u(NumOfObjs), &
                   v(NumOfObjs), &
                   w(NumOfObjs)  )

         call mpi_irecv( u, NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )

         call mpi_irecv( v, NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr )

         call mpi_irecv( w, NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr )

         call mpi_waitall(3, recv_req, array_of_statuses, ierr)

         m = 0 
         do fID = 1, NumOfFaces
            associate( f => IBM_HO_faces(domain)% faces(fID) )
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               do k = 0, f% stencil(i,j)% N
                  m = m + 1 
                  f% stencil(i,j)% UU(IX,k) = u(m)
                  f% stencil(i,j)% UU(IY,k) = v(m)
                  f% stencil(i,j)% UU(IZ,k) = w(m)
               end do
            end do; end do 
            end associate
         end do

         deallocate( u, &
                     v, &
                     w  )

      end do
#endif
   end subroutine RecvFacesStatePartitions

   subroutine SendFacesState2Partitions( )
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables----------------------------------------------------------
      real(kind=RP), allocatable :: u(:), v(:), w(:)
      integer,       allocatable :: send_req(:,:)
      integer                    :: i, j, k, m, nProcs, ierr, domain, fID, &
                                    NumOfObjs, NumOfFaces,                 &
                                    array_of_statuses(MPI_STATUS_SIZE,3)

      allocate(send_req(MPI_Process% nProcs-1,3))

      do nProcs = 2, MPI_Process% nProcs
         do domain = 0, MPI_Process% nProcs-1

            NumOfFaces = IBM_HO_faces(domain)% NumOfFaces
            NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs

            if( NumOfFaces .eq. 0 ) cycle
 
            allocate( u(NumOfObjs), &
                      v(NumOfObjs), &
                      w(NumOfObjs)  )
            
            m = 0 
            do fID = 1, NumOfFaces
               associate( f => IBM_HO_faces(domain)% faces(fID) )
               do i = 0, f% Nf(1); do j = 0, f% Nf(2)
                  do k = 0, f% stencil(i,j)% N
                     m = m + 1 
                     u(m) = f% stencil(i,j)% UU(IX,k)
                     v(m) = f% stencil(i,j)% UU(IY,k)
                     w(m) = f% stencil(i,j)% UU(IZ,k)
                  end do 
               end do; end do 
               end associate
            end do   

            call mpi_isend( u, NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )

            call mpi_isend( v, NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )
 
            call mpi_isend( w, NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )
      
            call mpi_waitall(3, send_req(nProcs-1,:), array_of_statuses, ierr)

            deallocate( u, &
                        v, &
                        w  )
         end do 
      end do
      
      deallocate( send_req )
#endif
   end subroutine SendFacesState2Partitions

   subroutine RecvFacesStateRoot()
      implicit none 
#ifdef _HAS_MPI_
      !-local-variables-------------------------------------------------------------------
      real(kind=RP), allocatable :: u(:), v(:), w(:)
      logical,       allocatable :: state(:)
      integer,       allocatable :: recv_req(:,:)
      integer                    :: i, j, k, m, N, fID, ierr, l,                  &
                                    nProcs, array_of_statuses(MPI_STATUS_SIZE,4), &
                                    NumOfObjs, NumOfFaces, NumOfDoFs, domain

      if( .not. MPI_Process% isRoot ) return

      allocate( recv_req(MPI_Process% nProcs-1,4) )

      do nProcs = 2, MPI_Process% nProcs 
         do domain = 0, MPI_Process% nProcs-1

            NumOfFaces = IBM_HO_faces(domain)% NumOfFaces
            NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs
            NumOfDoFs  = IBM_HO_faces(domain)% NumOfDoFs

            if( NumOfFaces .eq. 0 ) cycle 

            allocate( u(NumOfObjs),    &
                      v(NumOfObjs),    &
                      w(NumOfObjs),    &
                      state(NumOfDoFs) )

            call mpi_irecv( u, NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )

            call mpi_irecv( v, NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )

            call mpi_irecv( w, NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )
            
            call mpi_irecv( state, NumOfDoFs, MPI_LOGICAL, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )

            call mpi_waitall(4, recv_req(nProcs-1,:), array_of_statuses, ierr)         
         
            m = 0; l = 0
            do fID = 1, NumOfFaces
               associate( f => IBM_HO_faces(domain)% faces(fID) )
               do i = 0, f% Nf(1); do j = 0, f% Nf(2) 
                  l = l + 1 
                  if( .not. state(l) ) then 
                     m = m + (f% stencil(i,j)% N + 1)
                  else
                     do k = 0, f% stencil(i,j)% N
                        m = m + 1 
                        f% stencil(i,j)% UU(IX,k) = u(m) 
                        f% stencil(i,j)% UU(IY,k) = v(m) 
                        f% stencil(i,j)% UU(IZ,k) = w(m) 
                     end do 
                  end if 
               end do; end do
               end associate
            end do 

            deallocate( u,    &
                        v,    &
                        w,    &
                        state )

         end do
      end do
      
      deallocate( recv_req )
#endif
   end subroutine RecvFacesStateRoot

   subroutine SendFacesState2Root()
      
      implicit none 
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: u(:), v(:), w(:)
      logical,       allocatable :: state(:)
      integer                    :: NumOfObjs, domain, fID, m, i, j, k, l, ierr, &
                                    array_of_statuses(MPI_STATUS_SIZE,4),        &
                                    send_req(4), NumOfFaces, NumOfDoFs

      if( MPI_Process% isRoot ) return

      do domain = 0, MPI_Process% nProcs-1

         NumOfFaces = IBM_HO_faces(domain)% NumOfFaces
         NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs  
         NumOfDoFs  = IBM_HO_faces(domain)% NumOfDoFs  

         if( NumOfFaces .eq. 0 ) cycle 

         allocate( u(NumOfObjs),    &
                   v(NumOfObjs),    &
                   w(NumOfObjs),    &
                   state(NumOfDoFs) )

         m = 0; l = 0
         do fID = 1, NumOfFaces
            associate( f => IBM_HO_faces(domain)% faces(fID) )
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               l = l + 1
               state(l) = f% stencil(i,j)% state
               do k = 0, f% stencil(i,j)% N 
                  m    = m + 1
                  u(m) = f% stencil(i,j)% UU(IX,k)
                  v(m) = f% stencil(i,j)% UU(IY,k)
                  w(m) = f% stencil(i,j)% UU(IZ,k)
               end do 
            end do; end do 
            end associate 
         end do 

         call mpi_isend( u, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )

         call mpi_isend( v, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )

         call mpi_isend( w, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )

         call mpi_isend( state, NumOfDoFs, MPI_LOGICAL, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )
         
         call mpi_waitall(4, send_req, array_of_statuses, ierr)

         deallocate( u,    &
                     v,    & 
                     w,    &
                     state )

      end do
#endif
   end subroutine SendFacesState2Root

   subroutine MPIProcedures_IBM_HO_faces()

      implicit none 

      if( MPI_Process% doMPIAction ) then
         call SendFaces2Root()
      end if
      
      if( MPI_Process% doMPIRootAction ) then 
         call RecvFacesRoot()
      end if

      if( MPI_Process% doMPIRootAction ) then 
         call SendFaces2partitions()
      end if
 
      if( MPI_Process% doMPIAction ) then
         call RecvFacespartitions()
      end if

   end subroutine MPIProcedures_IBM_HO_faces

   subroutine MPIStateProcedures_IBM_HO_faces()

      implicit none 

      if( MPI_Process% doMPIAction ) then
         call SendFacesState2Root()
      end if
      
      if( MPI_Process% doMPIRootAction ) then 
         call RecvFacesStateRoot()
      end if

      if( MPI_Process% doMPIRootAction ) then 
         call SendFacesState2partitions()
      end if
 
      if( MPI_Process% doMPIAction ) then
         call RecvFacesStatepartitions()
      end if

   end subroutine MPIStateProcedures_IBM_HO_faces

   subroutine IBM_HO_findElements( elements )
      use ElementClass
      implicit none 

      type(element), intent(inout) :: elements(:)

      real(kind=RP) :: xi(NDIM)
      logical       :: FOUND
      integer       :: domain, fID, i, j, k, eID
      
      do domain = 0, MPI_Process% nProcs-1
         do fID = 1, IBM_HO_faces(domain)% NumOfFaces
            associate( f => IBM_HO_faces(domain)% faces(fID) )
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               f% stencil(i,j)% state = .false.
               do eID = 1, size(elements)
                  associate(e => elements(eID))
                  if( e% FindPointWithCoords(f% stencil(i,j)% x_s(:,0), 0, xi) ) then
                     f% stencil(i,j)% state     = .true.
                     f% stencil(i,j)% partition = MPI_Process% rank 
                     f% stencil(i,j)% xi_s(:,0) = xi 
                     f% stencil(i,j)% eIDs(0)   = eID
                     do k = 1, f% stencil(i,j)% N 
                        FOUND = e% FindPointWithCoords(f% stencil(i,j)% x_s(:,k), 0, f% stencil(i,j)% xi_s(:,k))
                        f% stencil(i,j)% eIDs(k) = eID
                     end do
                     exit 
                  end if 
                  end associate
               end do 
            end do; end do 
            end associate 
         end do
      end do
  
   end subroutine IBM_HO_findElements

   subroutine Set_IBM_HO_faces( faces )

      implicit none 

      type(face), intent(inout) :: faces(:)

      integer :: i, j, k, fID, m, domain

      allocate( IBM_HO_faces(0:MPI_Process% nProcs-1) )
 
      domain = MPI_Process% rank

      IBM_HO_faces(domain)% NumOfObjs  = 0
      IBM_HO_faces(domain)% NumOfFaces = 0 
      IBM_HO_faces(domain)% NumOfDoFs  = 0 
      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then 
            IBM_HO_faces(domain)% NumOfFaces = IBM_HO_faces(domain)% NumOfFaces + 1
            IBM_HO_faces(domain)% NumOfObjs  = IBM_HO_faces(domain)% NumOfObjs + (f% stencil(i,j)% N+1)*(f% Nf(1)+1)*(f% Nf(2)+1) 
            IBM_HO_faces(domain)% NumOfDoFs  = IBM_HO_faces(domain)% NumOfDoFs + (f% Nf(1)+1)*(f% Nf(2)+1) 
         end if 
         end associate 
      end do
      
      allocate(IBM_HO_faces(domain)% faces(IBM_HO_faces(domain)% NumOfFaces))

      m = 0 
      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then 
            m = m + 1 
            IBM_HO_faces(domain)% faces(m)% Nf(1) = f% Nf(1)
            IBM_HO_faces(domain)% faces(m)% Nf(2) = f% Nf(2)
            IBM_HO_faces(domain)% faces(m)% ID    = f% ID 
            f% HO_ID                              = m
            f% domain                             = MPI_Process% rank
            call IBM_HO_faces(domain)% faces(m)% StencilConstruct()
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               IBM_HO_faces(domain)% faces(m)% stencil(i,j)% xiB = f% stencil(i,j)% xiB
               IBM_HO_faces(domain)% faces(m)% stencil(i,j)% xiI = f% stencil(i,j)% xiI
               do k = 0, f% stencil(i,j)% N 
                  IBM_HO_faces(domain)% faces(m)% stencil(i,j)% x_s(:,k) = f% stencil(i,j)% x_s(:,k)
               end do 
            end do; end do 
         end if 
         end associate 
      end do

   end subroutine Set_IBM_HO_faces

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
