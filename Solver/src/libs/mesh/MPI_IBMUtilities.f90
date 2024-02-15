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
   public :: castMaskNumOfObjs, castMask, castMaskPlot
   public :: gatherMaskGeom
   public :: castStateBandRegion, castGradientsBandRegion
   public :: IBM_HO_GetState
   public :: FixingmpiFaces
   public :: IBM_HOintegration_findElements, MPIProcedures_IBM_HOIntegrationPoints
   public :: IBM_HO_GetGradient, GatherHOIntegrationPointsState

   type IBMpoints

      type(point_type), allocatable :: x(:)
      integer                       :: LocNumOfObjs, NumOfObjs

      contains 
         procedure :: build    => IBMpoints_build 

   end type

contains

   subroutine IBMpoints_build( this, NumOfObjs )

      implicit none 

      class(IBMPoints),  intent(inout) :: this 
      integer,           intent(in)    :: NumOfObjs

      integer :: i
      
      allocate(this% x(NumOfObjs))
      this% NumOfObjs = NumOfObjs 

   end subroutine IBMpoints_build

   subroutine castStateBandRegion( IBMmask, nEqn )
      
      implicit none 
      
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: nEqn
#ifdef _HAS_MPI_
      integer :: domains, i, ierr 

      do domains = 1, MPI_Process% nProcs 
         do i = 1, nEqn
            call mpi_bcast( IBMmask(domains)% x(:)% Q(i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
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
            call mpi_bcast( IBMmask(domains)% x(:)% U_x(i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
            call mpi_bcast( IBMmask(domains)% x(:)% U_y(i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
            call mpi_bcast( IBMmask(domains)% x(:)% U_z(i), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         end do 
      end do
#endif 
   end subroutine castGradientsBandRegion
!__________________________________________________

   subroutine castMaskPlot( IBMMask )
      
      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      integer :: domains, ierr

      if ( .not. MPI_Process % doMPIAction ) return

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMmask(domains)% x(:)% isInsideBody, IBMmask(domains)% NumOfObjs, MPI_LOGICAL, domains-1, MPI_COMM_WORLD, ierr )     
      end do 
#endif 
   end subroutine castMaskPlot


   subroutine gatherMaskGeom( IBMmask )

      implicit none 

      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: dist(:), normal(:,:) 
      integer                    :: i, domain, domains, ierr, index

      if ( .not. MPI_Process % doMPIAction ) return

      domain = MPI_Process% rank + 1
      
      allocate( dist(IBMmask(domain)% NumOfObjs*MPI_Process% nProcs),       &
                normal(IBMmask(domain)% NumOfObjs*MPI_Process% nProcs,NDIM) )

      do domains = 1, MPI_Process% nProcs 
         call mpi_gather( IBMmask(domains)% x(:)% dist      , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, dist,          &
                          IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                   )
         call mpi_gather( IBMmask(domains)% x(:)% normal(IX), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, normal(:,IX), &
                          IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                   )
         call mpi_gather( IBMmask(domains)% x(:)% normal(IY), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, normal(:,IY), &
                          IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                   )
         call mpi_gather( IBMmask(domains)% x(:)% normal(IZ), IBMmask(domains)% NumOfObjs, MPI_DOUBLE, normal(:,IZ), &
                          IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                   )
      end do 

      do i = 1, IBMmask(domain)% NumOfObjs 
         do domains = 1, MPI_Process% nProcs
            index = (domains-1)*IBMmask(domain)% NumOfObjs + i   
            if( IBMmask(domain)% x(i)% dist .gt. dist(index) ) then 
               IBMmask(domain)% x(i)% dist   = dist(index) 
               IBMmask(domain)% x(i)% normal = normal(index,:)
            end if 
         end do
      end do 

      deallocate( dist, normal )
#endif
   end subroutine gatherMaskGeom

   subroutine castMaskNumOfObjs( IBMmask, domain )

      implicit none 
      type(IBMPoints), intent(inout) :: IBMmask(:)
      integer,         intent(in)    :: domain 
#ifdef _HAS_MPI_
      integer :: domains, ierr

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMmask(domains)% NumOfObjs, 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )     
      end do 

      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         allocate( IBMmask(domains)% x(IBMmask(domains)% NumOfObjs) )
      end do
#endif 
   end subroutine castMaskNumOfObjs

   subroutine castMask( IBMmask )
      
      implicit none 
      
      type(IBMPoints), intent(inout) :: IBMmask(:)
#ifdef _HAS_MPI_
      integer :: domains, ierr 

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBMmask(domains)% x(:)% coords(IX)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% x(:)% coords(IY)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% x(:)% coords(IZ)        , IBMmask(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% x(:)% element_index     , IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% x(:)% local_position(IX), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% x(:)% local_position(IY), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBMmask(domains)% x(:)% local_position(IZ), IBMmask(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )
      end do
#endif 
   end subroutine castMask

   subroutine SendAxis( STLNum ) 

      implicit none 

      integer, intent(in) :: STLNum 
#ifdef _HAS_MPI_
      integer             :: domains, send_req(MPI_Process% nProcs-1,2), &
                             msg, ierr

      if( .not. MPI_Process% isRoot ) return 
      
      do domains = 2, MPI_Process% nProcs

         call mpi_isend(OBB(STLNum)% maxAxis, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,1), ierr )
         call mpi_isend(OBB(STLNum)% minAxis, 1, MPI_INT, domains-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(domains-1,2), ierr )

      end do

      do msg = 1, 2 
         call mpi_waitall( MPI_Process% nProcs-1, send_req(:,msg), MPI_STATUSES_IGNORE, ierr )
      end do 
#endif
   end subroutine SendAxis

   subroutine RecvAxis( STLNum )

      implicit none 

      integer, intent(in) :: STLNum
#ifdef _HAS_MPI_
      integer :: recv_req(2), ierr
      
      if( MPI_Process% isRoot ) return 

      call mpi_irecv(OBB(STLNum)% MaxAxis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr)
      call mpi_irecv(OBB(STLNum)% minAxis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr)

      call mpi_waitall( 2, recv_req, MPI_STATUSES_IGNORE, ierr )
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

      deallocate(vertices_x, vertices_y, vertices_z, normals_x, normals_y, normals_z)
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

      integer :: i, j, k, fID, m, domain, counter, counter1

      allocate( IBM_HO_faces(MPI_Process% nProcs) )
 
      domain = MPI_Process% rank+1

      IBM_HO_faces(domain)% NumOfObjs     = 0
      IBMStencilPoints(domain)% NumOfObjs = 0 

      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then 
            IBM_HO_faces(domain)% NumOfFaces    = IBM_HO_faces(domain)% NumOfFaces    + 1
            IBMStencilPoints(domain)% NumOfObjs = IBMStencilPoints(domain)% NumOfObjs + max(f% Nf(1)+1,f% Nf(2)+1) * (f% Nf(1)+1) * (f% Nf(2)+1) 
         end if 
         end associate 
      end do
      
      allocate( IBM_HO_faces(domain)% faces(IBM_HO_faces(domain)% NumOfFaces) )

      call IBMStencilPoints(domain)% build(IBMStencilPoints(domain)% NumOfObjs)
      
      m = 0; counter = 1; counter1 = 1
      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then 
            m = m + 1 
            IBM_HO_faces(domain)% faces(m)% Nf(1) = f% Nf(1)
            IBM_HO_faces(domain)% faces(m)% Nf(2) = f% Nf(2)
            IBM_HO_faces(domain)% faces(m)% ID    = f% ID 
            f% HO_ID                              = m
            f% domain                             = domain
            call IBM_HO_faces(domain)% faces(m)% StencilConstruct()
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               IBM_HO_faces(domain)% faces(m)% stencil(i,j)% xiB = f% stencil(i,j)% xiB 
               IBM_HO_faces(domain)% faces(m)% stencil(i,j)% xiI = f% stencil(i,j)% xiI 
               do k = 0, f% stencil(i,j)% N
                  IBMStencilPoints(domain)% x(counter1)% coords         = f% stencil(i,j)% x_s(:,k)
                  IBMStencilPoints(domain)% x(counter1)% N              = f% stencil(i,j)% N
                  IBMStencilPoints(domain)% x(counter1)% local_position =(/i,j,k/)
                  IBMStencilPoints(domain)% x(counter1)% faceID         = m
                  IBMStencilPoints(domain)% x(counter1)% element_index  = f% HOeID
                  IBMStencilPoints(domain)% x(counter1)% domain         = f% HOdomain
                  counter1 = counter1 + 1
               end do 
               counter = counter + 1
            end do; end do 
         end if 
         end associate 
      end do
      
      counter1 = 1
      do m = 1, IBM_HO_faces(domain)% NumOfFaces
         associate( f => IBM_HO_faces(domain)% faces(m) )
         do i = 0, f% Nf(1); do j = 0, f% Nf(2)
            do k = 0, f% stencil(i,j)% N
               f% stencil(i,j)% x_s(:,k) = IBMStencilPoints(domain)% x(counter1)% coords
               counter1 = counter1 + 1
            end do
         end do; end do 
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
         call mpi_bcast( IBM_HO_faces(domains)% NumOfFaces   , 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )     
         call mpi_bcast( IBMStencilPoints(domains)% NumOfObjs, 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )     
      end do 
      
      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         allocate( IBM_HO_faces(domains)% faces(IBM_HO_faces(domains)% NumOfFaces) )
         call IBMStencilPoints(domains)% build(IBMStencilPoints(domains)% NumOfObjs)
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
         call mpi_bcast( IBM_HO_faces(domains)% faces(:)% Nf(1)     , IBM_HO_faces(domains)% NumOfFaces   , MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )     
         call mpi_bcast( IBM_HO_faces(domains)% faces(:)% Nf(2)     , IBM_HO_faces(domains)% NumOfFaces   , MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )        
         call mpi_bcast( IBMStencilPoints(domains)% x(:)% coords(IX), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )     
         call mpi_bcast( IBMStencilPoints(domains)% x(:)% coords(IY), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )     
         call mpi_bcast( IBMStencilPoints(domains)% x(:)% coords(IZ), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )     
         call mpi_bcast( IBMStencilPoints(domains)% x(:)% N         , IBMStencilPoints(domains)% NumOfObjs, MPI_INT   , domains-1, MPI_COMM_WORLD, ierr )       
      end do 

      do domains = 1, MPI_Process% nProcs
         if( domains .eq. domain ) cycle 
         counter1 = 1
         do m = 1, IBM_HO_faces(domains)% NumOfFaces
            call IBM_HO_faces(domains)% faces(m)% StencilConstruct()
            associate( f => IBM_HO_faces(domains)% faces(m) )
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               do k = 0, f% stencil(i,j)% N
                  f% stencil(i,j)% x_s(:,k) = IBMStencilPoints(domains)% x(counter1)% coords
                  counter1 = counter1 + 1
               end do
            end do; end do 
            end associate
         end do
      end do 
#endif 
   end subroutine CastHOfaces

!================================================== HO faces state ==================================================
   subroutine IBM_HO_findElements( IBMStencilPoints, elements )
      use ElementClass
      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      type(element),   intent(inout) :: elements(:)

      real(kind=RP) :: xi(NDIM)
      logical       :: FOUND
      integer       :: domain, domains, fID, i, j, k, eID, counter, counter1 

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs  
         !counter = 1
         counter = 0
         !do while( counter .le. IBMStencilPoints(domains)% NumOfObjs )
         do while( counter .lt. IBMStencilPoints(domains)% NumOfObjs )
            counter = counter + IBMStencilPoints(domains)% x(counter+1)% N + 1 
            FOUND = .false.
            do eID = 1, size(elements)
               associate( e => elements(eID) )
               if( ALL(e% MaskCorners) ) cycle
               if( e% FindPointWithCoords(IBMStencilPoints(domains)% x(counter)% coords, 0, xi) ) then
                  IBMStencilPoints(domains)% x(counter)% domain        = domain
                  IBMStencilPoints(domains)% x(counter)% element_index = eID 
                  IBMStencilPoints(domains)% x(counter)% xi            = xi 
                  counter1 = counter 
                  do k = 1, IBMStencilPoints(domains)% x(counter1)% N 
                     counter1 = counter1 - 1 
                     !counter1 = counter1 + 1 
                     FOUND = e% FindPointWithCoords(IBMStencilPoints(domains)% x(counter1)% coords, 0, xi)
                     IBMStencilPoints(domains)% x(counter1)% domain        = domain
                     IBMStencilPoints(domains)% x(counter1)% element_index = eID 
                     IBMStencilPoints(domains)% x(counter1)% xi            = xi 
                  end do
                  FOUND = .true. 
                  exit 
               end if
               end associate
            end do 
            !counter = counter + (IBMStencilPoints(domains)% x(counter)% N + 1) 
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
            if( IBMStencilPoints(domains)% x(i)% domain .eq. domain ) then 
               eID = IBMStencilPoints(domains)% x(i)% element_index 
               xi  = IBMStencilPoints(domains)% x(i)% xi
               IBMStencilPoints(domains)% x(i)% Q = elements(eID)% EvaluateSolutionAtPoint(nEqn, xi)
            end if 
         end do
      end do

   end subroutine IBM_HO_GetState

   subroutine GatherHOfacesState( IBMStencilPoints, nEqn )
      use PhysicsStorage
      implicit none 

      type(IBMpoints), intent(inout) :: IBMStencilPoints(:)
      integer,         intent(in)    :: nEqn
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: Q(:,:)
      integer,       allocatable :: eIDsdomain(:)
      integer                    :: domain, domains, i, j, k, m, n, ierr, index

      domain = MPI_Process% rank + 1

      if( MPI_Process% doMPIAction ) then 
      
         allocate( Q(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),    &
                   eIDsdomain(IBMStencilPoints(domain)% NumOfObjs*MPI_Process% nProcs) )

         do domains = 1, MPI_Process% nProcs 
            do i = 1, nEqn
               call mpi_gather( IBMStencilPoints(domains)% x(:)% Q(i), IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, Q(:,i), &
                                IBMStencilPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                )
            end do 
            call mpi_gather( IBMStencilPoints(domains)% x(:)% domain, IBMStencilPoints(domains)% NumOfObjs, MPI_INT, eIDsdomain, &
                             IBMStencilPoints(domains)% NumOfObjs, MPI_INT, domains-1, MPI_COMM_WORLD, ierr                      )
         end do
         
         do n = 1, IBMStencilPoints(domain)% NumOfObjs
            do domains = 1, MPI_Process% nProcs 
               index = (domains-1)*IBMStencilPoints(domain)% NumOfObjs + n 
               if( eIDsdomain(index) .ne. 0 ) then 
                  m = IBMStencilPoints(domain)% x(n)% faceID
                  i = IBMStencilPoints(domain)% x(n)% local_position(IX)
                  j = IBMStencilPoints(domain)% x(n)% local_position(IY)
                  k = IBMStencilPoints(domain)% x(n)% local_position(IZ) 

                  IBM_HO_faces(domain)% faces(m)% stencil(i,j)% Q(:,k) = Q(index,:)
               end if 
            end do
         end do
         
         deallocate( Q, eIDsdomain )

      else 

         do n = 1, IBMStencilPoints(domain)% NumOfObjs
            m = IBMStencilPoints(domain)% x(n)% faceID
            i = IBMStencilPoints(domain)% x(n)% local_position(IX)
            j = IBMStencilPoints(domain)% x(n)% local_position(IY)
            k = IBMStencilPoints(domain)% x(n)% local_position(IZ)

            IBM_HO_faces(domain)% faces(m)% stencil(i,j)% Q(:,k) = IBMStencilPoints(domain)% x(n)% Q
         end do 

      end if
#endif
   end subroutine GatherHOfacesState

   subroutine FixingmpiFaces( faces, MPIfaces, NumOfMaskObjs )
      use MPI_Face_Class
      use MeshTypes, only: HMESH_MPI
      implicit none 
 
      type(face),           intent(inout) :: faces(:) 
      type(MPI_FacesSet_t), intent(inout) :: MPIfaces
      integer,              intent(inout) :: NumOfMaskObjs
#ifdef _HAS_MPI_
      integer              :: domain, domains, m, mpifID, fID, i, j, ierr

      if( .not. MPI_Process% doMPIAction ) return 

      domain = MPI_Process% rank + 1
      
      allocate(IBM_HO_mpifaces(MPI_Process% nProcs))

      IBM_HO_mpifaces(domain)% NumOfFaces = 0

      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         do mpifID = 1, MPIfaces% faces(domains)% no_of_faces
            fID = MPIfaces% faces(domains)% faceIDs(mpifID)
            if( faces(fID)% HO_IBM ) then 
               IBM_HO_mpifaces(domain)% NumOfFaces = IBM_HO_mpifaces(domain)% NumOfFaces + 1
            end if 
         end do
      end do 

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBM_HO_mpifaces(domains)% NumOfFaces, 1, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )  
      end do 
      
      do domains = 1, MPI_Process% nProcs 
         allocate(IBM_HO_mpifaces(domains)% faces(IBM_HO_mpifaces(domains)% NumOfFaces) )
      end do
      
      m = 0

      do domains = 1, MPI_Process% nProcs 
         if( domains .eq. domain ) cycle 
         do mpifID = 1, MPIfaces% faces(domains)% no_of_faces
            fID = MPIfaces% faces(domains)% faceIDs(mpifID)
            if( faces(fID)% HO_IBM ) then 
               m = m + 1
               IBM_HO_mpifaces(domain)% faces(m)% ID     = mpifID
               IBM_HO_mpifaces(domain)% faces(m)% domain = domains
               faces(fID)% HO_IBM   = .false.
               faces(fID)% corrGrad = .true. 
               deallocate(faces(fID)% stencil)
               NumOfMaskObjs = NumOfMaskObjs - (faces(fID)% Nf(1)+1)*(faces(fID)% Nf(2)+1)
            end if 
         end do
      end do
      
      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBM_HO_mpifaces(domains)% faces(:)% ID    , IBM_HO_mpifaces(domains)% NumOfFaces, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )
         call mpi_bcast( IBM_HO_mpifaces(domains)% faces(:)% domain, IBM_HO_mpifaces(domains)% NumOfFaces, MPI_INT, domains-1, MPI_COMM_WORLD, ierr )
      end do 

      do domains = 1, MPI_Process% nProcs
         do m = 1, IBM_HO_mpifaces(domains)% NumOfFaces
            if( IBM_HO_mpifaces(domains)% faces(m)% domain .eq. domain ) then 
               mpifID = IBM_HO_mpifaces(domains)% faces(m)% ID
               fID    = MPIfaces% faces(domains)% faceIDs(mpifID)
               faces(fID)% HO_IBM = .true. 
               faces(fID)% HOSIDE = maxloc(faces(fID)% elementIDs, dim=1)
               allocate(faces(fID)% stencil(0:faces(fID)% Nf(1),0:faces(fID)% Nf(2)))
               do j = 0, faces(fID)% Nf(2); do i = 0, faces(fID)% Nf(1)
                  faces(fID)% stencil(i,j)% x = faces(fID)% geom% x(:,i,j)
               end do; end do
               NumOfMaskObjs = NumOfMaskObjs + (faces(fID)% Nf(1)+1)*(faces(fID)% Nf(2)+1)
            end if
         end do
      end do

      do domains = 1, MPI_Process% nProcs
         deallocate(IBM_HO_mpifaces(domains)% faces)
      end do
      
      deallocate(IBM_HO_mpifaces)
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

      do domains = 1, MPI_Process% nProcs 
         call mpi_bcast( IBM_HOIntegrationPoints(domains)% x(:)% coords(IX), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )  
         call mpi_bcast( IBM_HOIntegrationPoints(domains)% x(:)% coords(IY), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )  
         call mpi_bcast( IBM_HOIntegrationPoints(domains)% x(:)% coords(IZ), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr )  
      end do 
#endif
   end subroutine CastHOIntegrationPoints

   subroutine IBM_HOintegration_findElements( IBM_HOIntegrationPoints, clipAxis, elements )
      use ElementClass
      implicit none 

      type(IBMpoints), intent(inout) :: IBM_HOIntegrationPoints(:)
      integer,         intent(in)    :: clipAxis
      type(element),   intent(inout) :: elements(:)

      real(kind=RP) :: xi(NDIM), x(NDIM)
      integer       :: domain, domains, i, eID, ierr
      logical       :: found = .false.

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs  
         do i = 1, IBM_HOIntegrationPoints(domains)% NumOfObjs 
            do eID = 1, size(elements)
               associate( e => elements(eID) ) 
               if( any(e% MaskCorners) ) then 
                   if( e% FindPointWithCoords(IBM_HOIntegrationPoints(domains)% x(i)% coords, clipAxis, xi) ) then
                     found = e% FindPointWithCoords(IBM_HOIntegrationPoints(domains)% x(i)% coords, 0, xi)
                     IBM_HOIntegrationPoints(domains)% x(i)% domain        = domain
                     IBM_HOIntegrationPoints(domains)% x(i)% element_index = eID 
                     IBM_HOIntegrationPoints(domains)% x(i)% xi            = xi 
                     exit 
                  end if
               end if 
               end associate 
            end do
         end do
      end do
      
   end subroutine IBM_HOintegration_findElements

   subroutine IBM_HO_GetGradient( IBM_HOIntegrationPoints, elements, nEqn )
      use ElementClass
      implicit none 

      type(IBMpoints), intent(inout) :: IBM_HOIntegrationPoints(:)
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: nEqn

      real(kind=RP) :: xi(NDIM)
      integer       :: domain, domains, i, eID 

      domain = MPI_Process% rank + 1

      do domains = 1, MPI_Process% nProcs 
         do i = 1, IBM_HOIntegrationPoints(domains)% NumOfObjs
            if( IBM_HOIntegrationPoints(domains)% x(i)% domain .eq. domain ) then 
               eID = IBM_HOIntegrationPoints(domains)% x(i)% element_index 
               xi  = IBM_HOIntegrationPoints(domains)% x(i)% xi
               IBM_HOIntegrationPoints(domains)% x(i)% U_x = elements(eID)% EvaluateGradientAtPoint(nEqn, xi, IX)
               IBM_HOIntegrationPoints(domains)% x(i)% U_y = elements(eID)% EvaluateGradientAtPoint(nEqn, xi, IY)
               IBM_HOIntegrationPoints(domains)% x(i)% U_z = elements(eID)% EvaluateGradientAtPoint(nEqn, xi, IZ)
            end if 
         end do
      end do

   end subroutine IBM_HO_GetGradient

   subroutine GatherHOIntegrationPointsState( IBM_HOIntegrationPoints, ObjectsList, nEqn )
      use PhysicsStorage
      implicit none 

      type(IBMpoints),   intent(inout) :: IBM_HOIntegrationPoints(:)
      type(Object_type), intent(inout) :: ObjectsList(:)
      integer,           intent(in)    :: nEqn
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: Q(:,:), U_x(:,:), U_y(:,:), U_z(:,:)
      integer,       allocatable :: eIDsdomain(:)
      integer                    :: domain, domains, i, j, k, n, ierr, index

      domain = MPI_Process% rank + 1
      
      if( MPI_Process% doMPIAction ) then 

         allocate( Q(IBM_HOIntegrationPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),    &
                   U_x(IBM_HOIntegrationPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),  &
                   U_y(IBM_HOIntegrationPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),  &
                   U_z(IBM_HOIntegrationPoints(domain)% NumOfObjs*MPI_Process% nProcs,nEqn),  &
                   eIDsdomain(IBM_HOIntegrationPoints(domain)% NumOfObjs*MPI_Process% nProcs) )

         do domains = 1, MPI_Process% nProcs 
            do i = 1, nEqn
               call mpi_gather( IBM_HOIntegrationPoints(domains)% x(:)% Q(i), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, Q(:,i)    , &
                                IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                           )
               call mpi_gather( IBM_HOIntegrationPoints(domains)% x(:)% U_x(i), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, U_x(:,i), &
                                IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                           )
               call mpi_gather( IBM_HOIntegrationPoints(domains)% x(:)% U_y(i), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, U_y(:,i), &
                                IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                           )
               call mpi_gather( IBM_HOIntegrationPoints(domains)% x(:)% U_z(i), IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, U_z(:,i), &
                                IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_DOUBLE, domains-1, MPI_COMM_WORLD, ierr                           )
            end do 
            call mpi_gather( IBM_HOIntegrationPoints(domains)% x(:)% domain, IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_INT, eIDsdomain, &
                             IBM_HOIntegrationPoints(domains)% NumOfObjs, MPI_INT, domains-1, MPI_COMM_WORLD, ierr                             )
         end do
         
         do n = 1, IBM_HOIntegrationPoints(domain)% NumOfObjs
            do domains = 1, MPI_Process% nProcs 
               index = (domains-1)*IBM_HOIntegrationPoints(domain)% NumOfObjs + n 
               if( eIDsdomain(index) .ne. 0 ) then 
                  i = IBM_HOIntegrationPoints(domain)% x(n)% local_position(IX)
                  j = IBM_HOIntegrationPoints(domain)% x(n)% local_position(IY)

                  ObjectsList(i)% IntegrationVertices(j)% Q   = Q(index,:)
                  ObjectsList(i)% IntegrationVertices(j)% U_x = U_x(index,:)
                  ObjectsList(i)% IntegrationVertices(j)% U_y = U_y(index,:)
                  ObjectsList(i)% IntegrationVertices(j)% U_z = U_z(index,:)
               end if 
            end do
         end do

         deallocate( Q, U_x, U_y, U_z, eIDsdomain )

      else 

         do n = 1, IBM_HOIntegrationPoints(domain)% NumOfObjs
            i = IBM_HOIntegrationPoints(domain)% x(n)% local_position(IX)
            j = IBM_HOIntegrationPoints(domain)% x(n)% local_position(IY)

            ObjectsList(i)% IntegrationVertices(j)% Q   = IBM_HOIntegrationPoints(domain)% x(n)% Q
            ObjectsList(i)% IntegrationVertices(j)% U_x = IBM_HOIntegrationPoints(domain)% x(n)% U_x
            ObjectsList(i)% IntegrationVertices(j)% U_y = IBM_HOIntegrationPoints(domain)% x(n)% U_y
            ObjectsList(i)% IntegrationVertices(j)% U_z = IBM_HOIntegrationPoints(domain)% x(n)% U_z
         end do 

      end if 
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
