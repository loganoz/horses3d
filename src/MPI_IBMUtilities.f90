#include "Includes.h"
module MPI_IBMUtilities

   use SMConstants
   use Utilities
   use TessellationTypes
   use OrientedBoundingBox
   use MPI_Process_info
#ifdef _HAS_MPI_
   use mpi
#endif
   
   implicit none
   
   private
   public :: IBMpoints
   public :: recvOBB, sendOBB 
   public :: SendSTL2Partitions, receiveSTLpartitions
   public :: GetBRvertices
   public :: getmaskcandidates
   public :: RecvPointsListRoot, SendPointsList2Root
   public :: RecvPointsListPartitions, SendPointsList2partitions
   public :: recvPointsMaskRoot, sendPointsMask2Root
   public :: recvPointsMaskPartitions, sendPointsMask2Partitions
   public :: recvNormalsRoot, sendNormals2Root
   public :: recvDistancesANDNormalspartitions, sendDistanceANDNormals2partitions   
   public :: recvScalarPlotRoot, sendScalarPlotRoot
   public :: recvVectorPlotRoot, sendVectorPlotRoot

   type IBMpoints

      type(point_type), allocatable :: x(:)
      real(kind=RP),    allocatable :: Q(:,:), U_x(:,:), U_y(:,:), U_z(:,:)
      integer                       :: LocNumOfObjs, NumOfObjs

   end type
   
   type(IBMpoints), public :: Mask

contains

 subroutine recvOBB( STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------------
      integer, intent(in) :: STLNum
#ifdef _HAS_MPI_
      !-local-variables-----------------------------------------------------------
      real(kind=RP) :: vertex_x(8), vertex_y(8), vertex_z(8), &
                        angle, center(2), CloudCenter(NDIM),  &
                        R1(NDIM), R2(NDIM), R3(NDIM), Length,  &
                        Width, nMax, nMin
      integer       :: ierr, recv_req(13), array_of_statuses(MPI_STATUS_SIZE,13)
      
      if( MPI_Process% isRoot ) return
      
      call mpi_irecv( vertex_x, 8, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr ) 
      
      call mpi_irecv( vertex_y, 8, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr ) 
      
      call mpi_irecv( vertex_z, 8, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr ) 
      
      call mpi_irecv( angle, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr ) 
      
      call mpi_irecv( center, 2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr ) 
      
      call mpi_irecv( CloudCenter, NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr ) 
      
      call mpi_irecv( R1, NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr ) 
      
      call mpi_irecv( R2, NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr ) 
      
      call mpi_irecv( R3, NDIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9), ierr ) 
      
      call mpi_irecv( Length, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(10), ierr ) 
      
      call mpi_irecv( Width, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(11), ierr ) 
      
      call mpi_irecv( nMax, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(12), ierr ) 
      
      call mpi_irecv( nMin, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(13), ierr ) 
    
      call mpi_waitall(13, recv_req, array_of_statuses, ierr)
   
      OBB(STLNum)% LocVertices(1,:) = vertex_x
      OBB(STLNum)% LocVertices(2,:) = vertex_y
      OBB(STLNum)% LocVertices(3,:) = vertex_z
      OBB(STLNum)% MBR% angle       = angle 
      OBB(STLNum)% MBR% center      = center 
      OBB(STLNum)% CloudCenter      = CloudCenter
      OBB(STLNum)% R(:,1)           = R1
      OBB(STLNum)% R(:,2)           = R2
      OBB(STLNum)% R(:,3)           = R3
      OBB(STLNum)% MBR% Length      = Length
      OBB(STLNum)% MBR% Width       = Width
      OBB(STLNum)% nMax             = nMax
      OBB(STLNum)% nMin             = nMin

      OBB(STLNum)% invR(:,1) = OBB(STLNum)% R(1,:)
      OBB(STLNum)% invR(:,2) = OBB(STLNum)% R(2,:)
      OBB(STLNum)% invR(:,3) = OBB(STLNum)% R(3,:)
      
      OBB(STLNum)% MBR% t1     = OBB(STLNum)% R(1,:)
      OBB(STLNum)% MBR% t2     = OBB(STLNum)% R(2,:)
      OBB(STLNum)% MBR% normal = OBB(STLNum)% R(3,:)
#endif
   end subroutine recvOBB
   
   subroutine sendOBB( STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------------------------
      integer, intent(in) :: STLNum
#ifdef _HAS_MPI_
      !-local-variables----------------------------------------------------
      real(kind=RP)        :: vertex_x(8), vertex_y(8), vertex_z(8), &
                              angle, center(2), CloudCenter(NDIM),  &
                              R1(NDIM), R2(NDIM), R3(NDIM), Length,  &
                              Width, nMax, nMin
      integer, allocatable :: send_req(:,:)
      integer              :: array_of_statuses(MPI_STATUS_SIZE,13),  &
                              nProcs, ierr
                              
      allocate(send_req(MPI_Process% nProcs-1,13))
      
      vertex_x    = OBB(STLNum)% LocVertices(1,:)
      vertex_y    = OBB(STLNum)% LocVertices(2,:)
      vertex_z    = OBB(STLNum)% LocVertices(3,:)
      angle       = OBB(STLNum)% MBR% angle
      center      = OBB(STLNum)% MBR% center
      CloudCenter = OBB(STLNum)% CloudCenter
      R1          = OBB(STLNum)% R(:,1)
      R2          = OBB(STLNum)% R(:,2)
      R3          = OBB(STLNum)% R(:,3)
      Length      = OBB(STLNum)% MBR% Length
      Width       = OBB(STLNum)% MBR% Width
      nMax        = OBB(STLNum)% nMax
      nMin        = OBB(STLNum)% nMin
      

      do nProcs = 2, MPI_Process% nProcs
   
         call mpi_isend(vertex_x, 8, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )
         
         call mpi_isend(vertex_y, 8, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )
         
         call mpi_isend(vertex_z, 8, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )
         
         call mpi_isend(angle, 1, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )
         
         call mpi_isend(center, 2, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,5), ierr )
         
         call mpi_isend(CloudCenter, NDIM, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,6), ierr )
         
         call mpi_isend(R1, NDIM, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,7), ierr )
         
         call mpi_isend(R2, NDIM, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,8), ierr )
         
         call mpi_isend(R3, NDIM, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,9), ierr )
         
         call mpi_isend(Length, 1, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,10), ierr )
         
         call mpi_isend(Width, 1, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,11), ierr )
         
         call mpi_isend(nMax, 1, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,12), ierr )
         
         call mpi_isend(nMin, 1, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,13), ierr )

         call mpi_waitall(13, send_req(nProcs-1,:), array_of_statuses, ierr)
      
      end do
      
      deallocate(send_req)
#endif
   end subroutine sendOBB

   subroutine GetVertices( vertices, axis, SplittingPlanes )
   
      implicit none
      !-arguments-----------------------------------------------------------
      real(kind=RP), intent(inout) :: vertices(:,:)
      real(kind=RP), intent(in)    :: SplittingPlanes(2)
      integer,       intent(in)    :: axis 
      !-local-variables-----------------------------------------------------
      integer :: v_indeces(4,2)
   
      if( axis .eq. 1 ) then
         v_indeces(:,1) = (/1,4,8,5/)
         v_indeces(:,2) = (/2,3,7,6/)
      elseif( axis .eq. 2 ) then
         v_indeces(:,1) = (/1,5,6,2/)
         v_indeces(:,2) = (/4,8,7,3/)
      else
         v_indeces(:,1) = (/1,2,3,4/)
         v_indeces(:,2) = (/5,6,7,8/)      
      end if

      vertices(axis,v_indeces(:,1)) = SplittingPlanes(1)
      vertices(axis,v_indeces(:,2)) = SplittingPlanes(2)
   
   end subroutine GetVertices

   subroutine SendSTL2Partitions( rootSTL, STLNum, rootVertices, rootAxis )
      implicit none
      !-arguments---------------------------------------------------------------------
      type(STLfile), intent(inout) :: rootSTL
      integer,       intent(in)    :: STLNum
      real(kind=RP), intent(inout) :: rootVertices(:,:)
      integer,       intent(inout) :: rootAxis
      !-local-variables---------------------------------------------------------------
      integer                        :: maxvec(1)
#ifdef _HAS_MPI_
      real(kind=RP),     allocatable :: locVertices(:,:,:) , Bar(:), coord(:),    &
                                        normals_x(:), normals_y(:), normals_z(:), &
                                        vertices_x(:,:), vertices_y(:,:),         &
                                        vertices_z(:,:)
      real(kind=RP)                  :: SplittingPlanes(2), kdtreevertices(NDIM,8)
      integer                        :: NumOfObjs, NumOfObjsPP, NumOfObjsPartion, &
                                        start_index, final_index, i, j,           &
                                        nProcs, ierr,                             &
                                        array_of_statuses(MPI_STATUS_SIZE,17)
      integer, allocatable           :: SortedIndex(:), send_req(:,:)
#endif
      maxvec   = maxloc((/ OBB(STLNum)% MBR% Length,OBB(STLNum)% MBR% Width,abs(OBB(STLNum)% nMax) + abs(OBB(STLNum)% nMin) /))  
      rootAxis = maxvec(1)
#ifdef _HAS_MPI_
      rootSTL% partition = MPI_Process% rank

      NumOfObjs   = size(rootSTL% ObjectsList)

      allocate( Bar(NumOfObjs),                          &
                coord(NumOfObjs),                        &
                SortedIndex(NumOfObjs),                  &
                normals_x(NumOfObjs),                    &
                normals_y(NumOfObjs),                    &
                normals_z(NumOfObjs),                    &
                vertices_x(NumOfObjs,3),                 &
                vertices_y(NumOfObjs,3),                 &
                vertices_z(NumOfObjs,3),                 &
                locVertices(8,NDIM,MPI_Process% nProcs), &
                send_req(MPI_Process% nProcs-1,17)       )

      do i = 1, NumOfObjs
         Bar(i) = 0.0_RP
         SortedIndex(i) = rootSTL% ObjectsList(i)% index
         do j = 1, size(rootSTL% ObjectsList(i)% vertices)
            Bar(i) = Bar(i) + rootSTL% ObjectsList(i)% vertices(j)% coords(rootAxis)
         end do
         Bar(i) = Bar(i)/size(rootSTL% ObjectsList(i)% vertices)
      end do
!$omp parallel 
!$omp single
      call sort( Bar, SortedIndex, coord, coord, coord, 1, NumOfObjs )
!$omp end single
!$omp end parallel
      deallocate(Bar,coord)

      NumOfObjsPP = NumOfObjs/MPI_Process% nProcs

      start_index = 0; final_index = 0

      do nProcs = 1, MPI_Process% nProcs
         start_index = final_index 
         final_index = start_index + NumOfObjsPP

         if( nProcs .eq. MPI_Process% nProcs ) final_index = NumOfObjs

         SplittingPlanes(1) =  huge(1.0_RP)
         SplittingPlanes(2) = -huge(1.0_RP)

         do i = 1, (final_index-start_index) 
            do j = 1, 3
               SplittingPlanes(1) = min(SplittingPlanes(1),rootSTL% ObjectsList(SortedIndex(start_index+i))% vertices(j)% coords(rootAxis))
               SplittingPlanes(2) = max(SplittingPlanes(2),rootSTL% ObjectsList(SortedIndex(start_index+i))% vertices(j)% coords(rootAxis))
               vertices_x(start_index+i,j) = rootSTL% ObjectsList(SortedIndex(start_index+i))% vertices(j)% coords(1)
               vertices_y(start_index+i,j) = rootSTL% ObjectsList(SortedIndex(start_index+i))% vertices(j)% coords(2)
               vertices_z(start_index+i,j) = rootSTL% ObjectsList(SortedIndex(start_index+i))% vertices(j)% coords(3)
            end do
            normals_x(start_index+i) = rootSTL% ObjectsList(SortedIndex(start_index+i))% normal(1)
            normals_y(start_index+i) = rootSTL% ObjectsList(SortedIndex(start_index+i))% normal(2)
            normals_z(start_index+i) = rootSTL% ObjectsList(SortedIndex(start_index+i))% normal(3)
         end do  

         kdtreevertices = OBB(STLNum)% LocVertices
   
         call GetVertices( kdtreevertices, rootAxis, SplittingPlanes )

         do j = 1, 8
            locVertices(j,:,nProcs) = kdtreevertices(:,j)
         end do

         if( rootSTL% show ) call DescribeSTLPartitions(nProcs-1,(final_index-start_index))

      end do

      start_index = NumOfObjsPP; final_index = NumOfObjsPP

      do nProcs = 2, MPI_Process% nProcs

         start_index  = final_index + 1
         final_index  = (start_index-1) + NumOfObjsPP

         if( nProcs .eq. MPI_Process% nProcs ) final_index = NumOfObjs
            
         NumOfObjsPartion = (final_index-start_index) + 1

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
 
         call mpi_isend(locVertices(:,1,nProcs), 8, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,14), ierr )
 
         call mpi_isend(locVertices(:,2,nProcs), 8, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,15), ierr )
  
         call mpi_isend(locVertices(:,3,nProcs), 8, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,16), ierr )

         call mpi_isend(rootAxis, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,17), ierr )  

         call mpi_waitall(17, send_req(nProcs-1,:), array_of_statuses, ierr)

      end do  

      call rootSTL% destroy()

      allocate(rootSTL% ObjectsList(NumOfObjsPP))

      do i = 1, NumOfObjsPP
         allocate(rootSTL% ObjectsList(i)% vertices(NDIM))
         rootSTL% ObjectsList(i)% normal(1) = normals_x(i)
         rootSTL% ObjectsList(i)% normal(2) = normals_y(i)
         rootSTL% ObjectsList(i)% normal(3) = normals_z(i)
         do j = 1, NDIM
            rootSTL% ObjectsList(i)% vertices(j)% coords(1) = vertices_x(i,j)
            rootSTL% ObjectsList(i)% vertices(j)% coords(2) = vertices_y(i,j)
            rootSTL% ObjectsList(i)% vertices(j)% coords(3) = vertices_z(i,j)
         end do 
         rootSTL% ObjectsList(i)% index = i
      end do 

      rootSTL% NumOfObjs = NumOfObjsPP

      do j = 1, 8
         rootVertices(:,j) = locVertices(j,:,1)
      end do

      deallocate(send_req, locVertices, vertices_x, vertices_y, vertices_z, normals_x, normals_y, normals_z)
#else
      if( rootSTL% show ) call rootSTL% Describe( rootSTL% filename )
#endif

   end subroutine SendSTL2Partitions


   subroutine receiveSTLpartitions( partitionSTL, STLNum, partitionVertex, partitionAxis )
      implicit none
      !-arguments-----------------------------------------------------
      type(STLfile), intent(inout) :: partitionSTL
      integer,       intent(in)    :: STLNum
      real(kind=RP), intent(inout) :: partitionVertex(:,:)
      integer,       intent(inout) :: partitionAxis
#ifdef _HAS_MPI_
      !-local-variables-----------------------------------------------
      real(kind=RP), allocatable :: normals_x(:), normals_y(:), normals_z(:), &
                                    vertices_x(:,:), vertices_y(:,:),         &
                                    vertices_z(:,:), locVertices(:,:)
      integer                    :: NumOfObjs, ierr, recv_req(17), i, j,      &
                                    array_of_statuses(MPI_STATUS_SIZE,17)
 
      if( MPI_Process% isRoot ) return

      partitionSTL% partition = MPI_Process% rank 

      call mpi_irecv( partitionSTL% NumOfObjs, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr )

      call mpi_wait(recv_req(1), MPI_STATUS_IGNORE, ierr)

      NumOfObjs = partitionSTL% NumOfObjs

      allocate( normals_x(NumOfObjs),       &
                normals_y(NumOfObjs),       &
                normals_z(NumOfObjs),       &
                vertices_x(NumOfObjs,NDIM), &
                vertices_y(NumOfObjs,NDIM), &
                vertices_z(NumOfObjs,NDIM), &
                locVertices(8,NDIM)         )

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

      call mpi_irecv( locVertices(:,1), 8, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(14), ierr )    
                   
      call mpi_irecv( locVertices(:,2), 8, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(15), ierr )  
                      
      call mpi_irecv( locVertices(:,3), 8, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(16), ierr ) 

      call mpi_irecv( partitionAxis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(17), ierr  )

      call mpi_waitall(17, recv_req, array_of_statuses, ierr)

      allocate(partitionSTL% ObjectsList(NumOfObjs))

      do i = 1, NumOfObjs
         partitionSTL% ObjectsList(i)% NumOfVertices = 3
         allocate(partitionSTL% ObjectsList(i)% vertices(3))
         do j = 1, partitionSTL% ObjectsList(i)% NumOfVertices
            partitionSTL% ObjectsList(i)% vertices(j)% coords(1) = vertices_x(i,j)
            partitionSTL% ObjectsList(i)% vertices(j)% coords(2) = vertices_y(i,j)
            partitionSTL% ObjectsList(i)% vertices(j)% coords(3) = vertices_z(i,j)
         end do
         partitionSTL% ObjectsList(i)% normal(1) = normals_x(i)
         partitionSTL% ObjectsList(i)% normal(2) = normals_y(i)
         partitionSTL% ObjectsList(i)% normal(3) = normals_z(i)
         partitionSTL% ObjectsList(i)% index = i
      end do

      partitionVertex(1,:) = locVertices(:,1)
      partitionVertex(2,:) = locVertices(:,2)
      partitionVertex(3,:) = locVertices(:,3)

      deallocate( normals_x,  &
                  normals_y,  &
                  normals_z,  &
                  vertices_x, &
                  vertices_y, &
                  vertices_z, &
                  locVertices )
#endif
   end subroutine receiveSTLpartitions

   subroutine GetBRvertices( vertices, BandRegionCoeff, axis, STLNum, BRvertices )
   
      implicit none
      !-arguments---------------------------------------------
      real(kind=RP), intent(in)  :: vertices(:,:),      &
                                    BandRegionCoeff 
      integer,       intent(in)  :: axis, STLNum
      real(kind=RP), intent(out) :: BRvertices(NDIM,8)
      !-local-variables---------------------------------------
      integer :: v_indeces(4,2)

      if( axis .eq. 1 ) then
         v_indeces(:,1) = (/1,4,8,5/)
         v_indeces(:,2) = (/2,3,7,6/)
      elseif( axis .eq. 2 ) then
         v_indeces(:,1) = (/1,5,6,2/)
         v_indeces(:,2) = (/4,8,7,3/)
      else
         v_indeces(:,1) = (/1,2,3,4/)
         v_indeces(:,2) = (/5,6,7,8/)      
      end if

      BRvertices = BandRegionCoeff*OBB(STLNum)% LocVertices

      if( MPI_Process% rank .ne. 0 ) then
         BRvertices(axis,v_indeces(:,1)) = vertices(axis,v_indeces(:,1))  
      end if
      if( MPI_Process% rank .ne. (MPI_Process% nProcs-1) ) then
         BRvertices(axis,v_indeces(:,2)) = vertices(axis,v_indeces(:,2))
      end if
      
   end subroutine GetBRvertices


   subroutine GetMaskCandidates( elements, no_of_elements, no_of_DoFs, STLNum, NumOfSTL ) 
      use ElementClass
      implicit none
      !-arguments----------------------------------------------------------------
      type(element),  intent(inout) :: elements(:)
      integer,        intent(in)    :: no_of_elements, no_of_DoFs, &
                                       STLNum, NumOfSTL
      !-local-variables-----------------------------------------------------------
      type(point_type), allocatable :: x(:)
      integer                       :: eID, i, j, k, ierr

      allocate(x(no_of_DoFs))
      Mask% LocNumOfObjs = 0

      do eID = 1, no_of_elements
         if( .not. allocated(elements(eID)% isInsideBody) ) then
            call elements(eID)% ConstructIBM(elements(eID)% Nxyz(1), elements(eID)% Nxyz(2), elements(eID)% Nxyz(3), NumOfSTL )
         end if
         
         do k = 0, elements(eID)% Nxyz(3); do j = 0, elements(eID)% Nxyz(2) ; do i = 0, elements(eID)% Nxyz(1)

            if( elements(eID)% isInsideBody(i,j,k) ) cycle

            elements(eID)% isInsideBody(i,j,k) = OBB(STLNum)% isPointInside( coords = elements(eID)% geom% x(:,i,j,k) )

            if( elements(eID)% isInsideBody(i,j,k) ) then
               Mask% LocNumOfObjs                    = Mask% LocNumOfObjs + 1
               x(Mask% LocNumOfObjs)% coords         = elements(eID)% geom% x(:,i,j,k)
               x(Mask% LocNumOfObjs)% local_Position = (/i,j,k/)
               x(Mask% LocNumOfObjs)% element_index  = eID
               x(Mask% LocNumOfObjs)% partition      = MPI_Process% rank
            end if
         end do; end do; end do
          
      end do
#ifdef _HAS_MPI_
      call mpi_allreduce(Mask% LocNumOfObjs, Mask% NumOfObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      Mask% NumOfObjs = Mask% LocNumOfObjs
#endif
      allocate(Mask% x(Mask% NumOfObjs))

      do i = 1, Mask% LocNumOfObjs
         Mask% x(i)% coords            = x(i)% coords
         Mask% x(i)% local_Position    = x(i)% local_Position
         Mask% x(i)% element_index     = x(i)% element_index
         Mask% x(i)% partition         = x(i)% partition
         Mask% x(i)% index             = i     
      end do

      deallocate(x)

   end subroutine GetMaskCandidates


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



   subroutine recvPointsMaskRoot()
      implicit none
      !-local-variables--------------------------------------------------
      integer              :: i 
#ifdef _HAS_MPI_      
      integer              :: nProcs, ierr, recv_req,   &
                              status(MPI_STATUS_SIZE) 
      integer, allocatable :: NumOfIntersectionsLoc(:), &
                              NumOfIntersections(:)  

      allocate( NumOfIntersectionsLoc(Mask% NumOfObjs), &
                NumOfIntersections(Mask% NumOfObjs)     )
            
      NumOfIntersections = 0

      do nProcs = 2, MPI_Process% nProcs
      
         call mpi_irecv( NumOfIntersectionsLoc, Mask% NumOfObjs, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req, ierr )  
            
         call mpi_wait(recv_req, status, ierr)                      
                      
         NumOfIntersections = NumOfIntersections + NumOfIntersectionsLoc

      end do

      do i = 1, Mask% NumOfObjs
         Mask% x(i)% NumOfIntersections = Mask% x(i)% NumOfIntersections + NumOfIntersections(i)
         Mask% x(i)% isInsideBody = .false.
         if(mod(Mask% x(i)% NumOfIntersections,2) .ne. 0 ) Mask% x(i)% isInsideBody = .true.
      end do

      deallocate(NumOfIntersections,NumOfIntersectionsLoc)
#else
      do i = 1, Mask% NumOfObjs
         Mask% x(i)% isInsideBody = .false.
         if(mod(Mask% x(i)% NumOfIntersections,2) .ne. 0 ) Mask% x(i)% isInsideBody = .true.
      end do
#endif          
   end subroutine recvPointsMaskRoot
   
   subroutine sendPointsMask2Root()
      implicit none
#ifdef _HAS_MPI_  
      !-local-variables-----------------------------------------------------------
      integer              :: ierr, send_req, i,      &
                              status(MPI_STATUS_SIZE)
      integer, allocatable :: NumOfIntersectionsLoc(:)
      
      if( MPI_Process% isRoot ) return
     
      allocate( NumOfIntersectionsLoc(Mask% NumOfObjs) )

      do i = 1, Mask% NumOfObjs
         NumOfIntersectionsLoc(i) = Mask% x(i)% NumOfIntersections
      end do

      call mpi_isend( NumOfIntersectionsLoc, Mask% NumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req, ierr )
      
      call mpi_wait(send_req, status, ierr)              
  
      deallocate( NumOfIntersectionsLoc )
#endif  
   end subroutine sendPointsMask2Root
   
 
   subroutine recvPointsMaskPartitions()
      implicit none
#ifdef _HAS_MPI_    
      !-local-variables-----------------------------------------------------------  
      integer              :: ierr, i, recv_req,      &
                              status(MPI_STATUS_SIZE)
      logical, allocatable :: isInsideBody(:)
      
      if( MPI_Process% isRoot ) return 

      allocate( isInsideBody(Mask% NumOfObjs) )      
 
      call mpi_irecv( isInsideBody, Mask% NumOfObjs, MPI_LOGICAL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req, ierr )  

      call mpi_wait(recv_req, status, ierr)     

      do i = 1, Mask% NumOfObjs
         Mask% x(i)% isInsideBody = isInsideBody(i)  
      end do 
 
      deallocate( isInsideBody )
#endif
   end subroutine recvPointsMaskPartitions

   subroutine sendPointsMask2Partitions()
      implicit none
#ifdef _HAS_MPI_  
      !-local-variables-----------------------------------------------------------------------------
      integer              :: nProcs, ierr, i, &
                              statuses(MPI_STATUS_SIZE)
      integer              :: send_req
      logical, allocatable :: isInsideBody(:)
      
      allocate( isInsideBody(Mask% NumOfObjs) )
        
      do i = 1, Mask% NumOfObjs
         isInsideBody(i) = Mask% x(i)% isInsideBody
      end do
      
      do nProcs = 2, MPI_Process% nProcs          
      
         call mpi_isend( isInsideBody, Mask% NumOfObjs, MPI_LOGICAL, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req, ierr )
         
         call mpi_wait(send_req, statuses, ierr )

      end do
      
      deallocate( isInsideBody )
#endif   
   end subroutine sendPointsMask2Partitions

   subroutine recvNormalsRoot( PointsList, ranks )
      implicit none 
      !-arguments-------------------------------------------------------------------------
      type(IBMPoints), intent(inout) :: PointsList
      real(kind=RP),   intent(in)    :: ranks(:)
      !-local-variables-------------------------------------------------------------------
#ifdef _HAS_MPI_
      real(kind=RP), allocatable :: normals_x(:,:), normals_y(:,:), normals_z(:,:)
      integer,       allocatable :: recv_req(:,:)
      integer                    :: i, ierr, rank, nProcs,                          &
                                    array_of_statuses(MPI_STATUS_SIZE,3)

      allocate( normals_x(PointsList% NumOfObjs,MPI_Process% nProcs-1), &
                normals_y(PointsList% NumOfObjs,MPI_Process% nProcs-1), &
                normals_z(PointsList% NumOfObjs,MPI_Process% nProcs-1), &
                recv_req(MPI_Process% nProcs-1,3)                       )
 
      do nProcs = 2, MPI_Process% nProcs 
 
         call mpi_irecv( normals_x(:,nProcs-1), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )

         call mpi_irecv( normals_y(:,nProcs-1), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )

         call mpi_irecv( normals_z(:,nProcs-1), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )

         call mpi_waitall(3, recv_req(nProcs-1,:), array_of_statuses, ierr ) 

      end do

      do i = 1, PointsList% NumOfObjs
         rank = ranks(i)
         if( rank .eq. 0 ) cycle
         PointsList% x(i)% normal(1) = normals_x(i,rank)
         PointsList% x(i)% normal(2) = normals_y(i,rank)
         PointsList% x(i)% normal(3) = normals_z(i,rank)
      end do 

      deallocate( normals_x, &
                  normals_y, &
                  normals_z, &
                  recv_req   )
#endif 
   end subroutine recvNormalsRoot

   subroutine sendNormals2Root( PointsList )
      implicit none 
      !-arguments-----------------------------------------------------------------
      type(IBMPoints), intent(in) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables----------------------------------------------------------
      real(kind=RP), allocatable :: normals(:,:)
      integer                    :: send_req(3), ierr, i,                &
                                    array_of_statuses(MPI_STATUS_SIZE,3)
          
      if( MPI_Process% isRoot ) return

      allocate(normals(PointsList% NumOfObjs,NDIM))

      do i = 1, PointsList% NumOfObjs
         normals(i,:) = PointsList% x(i)% normal
      end do

      call mpi_isend( normals(:,1), PointsList% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )

      call mpi_isend( normals(:,2), PointsList% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )

      call mpi_isend( normals(:,3), PointsList% NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )

      call mpi_waitall(3, send_req, array_of_statuses, ierr)

      deallocate(normals)
#endif
   end subroutine sendNormals2Root

   subroutine recvDistancesANDNormalspartitions( PointsList )
      implicit none
      !-arguments------------------------------------------------------------------------------
      type(IBMPoints), intent(inout) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables------------------------------------------------------------------------
      real(kind=RP), allocatable :: normals(:,:), Dist(:)
      integer                    :: i, ierr, recv_req(4),                &
                                    array_of_statuses(MPI_STATUS_SIZE,4)

      if( MPI_Process% isRoot ) return

      allocate( normals(PointsList% NumOfObjs,NDIM), &
                Dist(PointsList% NumOfObjs)          )

      call mpi_irecv( normals(:,1), PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr ) 

      call mpi_irecv( normals(:,2), PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr ) 
      
      call mpi_irecv( normals(:,3), PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr ) 

      call mpi_irecv( Dist, PointsList% NumOfObjs, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr ) 

      call mpi_waitall(4, recv_req, array_of_statuses, ierr)

      do i = 1, PointsList% NumOfObjs
         PointsList% x(i)% normal = normals(i,:)
         PointsList% x(i)% Dist   = Dist(i)
      end do

      deallocate( normals, &
                  Dist     )
#endif
   end subroutine recvDistancesANDNormalspartitions

   subroutine sendDistanceANDNormals2partitions( PointsList )
      implicit none
      !-arguments------------------------------------------------------------------------------
      type(IBMPoints), intent(in) :: PointsList
#ifdef _HAS_MPI_
      !-local-variables------------------------------------------------------------------------
      real(kind=RP), allocatable :: normals(:,:), Dist(:)
      integer,       allocatable :: send_req(:,:)
      integer                    :: i, ierr, nProcs, array_of_statuses(MPI_STATUS_SIZE,4)

      allocate( normals(PointsList% NumOfObjs,NDIM), &
                Dist(PointsList% NumOfObjs),         &
                send_req(MPI_Process% nProcs-1,4)    )

      do i = 1, PointsList% NumOfObjs
         normals(i,:) = PointsList% x(i)% normal
         Dist(i)      = PointsList% x(i)% Dist
      end do

      do nProcs = 2, MPI_Process% nProcs

         call mpi_isend( normals(:,1), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,1), ierr )

         call mpi_isend( normals(:,2), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,2), ierr )

         call mpi_isend( normals(:,3), PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,3), ierr )
         
         call mpi_isend( Dist, PointsList% NumOfObjs, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )

         call mpi_waitall(4, send_req(nProcs-1,:), array_of_statuses, ierr)

      end do

      deallocate( normals, &
                  Dist,    &
                  send_req )
#endif
   end subroutine sendDistanceANDNormals2partitions


   subroutine recvScalarPlotRoot( ObjectsList, STLNum, rootScalar, x, y, z, scalar )

      implicit none
      type(Object_type),          intent(in)    :: ObjectsList(:)
      integer,                    intent(in)    :: STLNum
      real(kind=RP),              intent(in)    :: rootScalar(:,:)
      real(kind=RP), allocatable, intent(inout) :: x(:), y(:), z(:), scalar(:)
      !-local-variables-------------------------------------------------------------------
      real(kind=RP)              :: coords(NDIM)
      integer                    :: i, j, k, rootNumOfObjs       
#ifdef _HAS_MPI_      
      integer                    :: ierr, nProcs, rank, ObjsSize, start_index,        &
                                    final_index, array_of_statuses(MPI_STATUS_SIZE,4)
      integer,       allocatable :: recv_req(:,:), recv_Firstreq(:), NumOfObjs(:)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state(:)
#endif
      if( .not. MPI_Process% isRoot ) return

      rootNumOfObjs = 3*size(ObjectsList)
#ifdef _HAS_MPI_ 

      allocate( NumOfObjs(MPI_Process% nProcs-1),    &
                recv_Firstreq(MPI_Process% nProcs-1) )

      do nProcs = 2, MPI_Process% nProcs
         
         call mpi_irecv( ObjsSize, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_Firstreq(nProcs-1), ierr ) 
         call mpi_wait( recv_Firstreq(nProcs-1), MPI_STATUS_IGNORE, ierr )

         NumOfObjs(nProcs-1) = ObjsSize

      end do 

      deallocate(recv_Firstreq)

      allocate( x(rootNumOfObjs+sum(NumOfObjs)),      &
                y(rootNumOfObjs+sum(NumOfObjs)),      &
                z(rootNumOfObjs+sum(NumOfObjs)),      &
                scalar(rootNumOfObjs+sum(NumOfObjs)), &
                COORD_x(sum(NumOfObjs)),              &
                COORD_y(sum(NumOfObjs)),              &
                COORD_z(sum(NumOfObjs)),              &
                state(sum(NumOfObjs)),                &
                recv_req(MPI_Process% nProcs-1,4)     )
 
      do nProcs = 2, MPI_Process% nProcs  

         start_index = sum(NumOfObjs(1:nProcs-2)) + 1 
         final_index = (start_index-1) + NumOfObjs(nProcs-1)

         call mpi_irecv( COORD_x(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )

         call mpi_irecv( COORD_y(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )
        
         call mpi_irecv( COORD_z(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )
        
         call mpi_irecv( state(start_index:final_index),   NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )

         call mpi_waitall(4, recv_req(nProcs-1,:), array_of_statuses, ierr)  

      end do
      do i = 1, sum(NumOfObjs)
         x(rootNumOfObjs+i)      = COORD_x(i)
         y(rootNumOfObjs+i)      = COORD_y(i)
         z(rootNumOfObjs+i)      = COORD_z(i)
         scalar(rootNumOfObjs+i) = state(i)
      end do  

      deallocate( NumOfObjs, COORD_x, COORD_y, COORD_z, state, recv_req )
#else
      allocate( x(rootNumOfObjs),      &
                y(rootNumOfObjs),      &
                z(rootNumOfObjs),      &
                scalar(rootNumOfObjs)  )
#endif  
      k = 0
      do i = 1, size(ObjectsList)
         do j = 1, size(ObjectsList(i)% vertices)
            call OBB(STLNum)% ChangeRefFrame(ObjectsList(i)% vertices(j)% coords,GLOBAL,coords)
            k         = k + 1
            x(k)      = coords(1)
            y(k)      = coords(2)
            z(k)      = coords(3)
            scalar(k) = Rootscalar(i,j)
         end do
      end do

   end subroutine recvScalarPlotRoot

   subroutine sendScalarPlotRoot( ObjectsList, STLNum, partitionScalar )
      implicit none 

      !-arguments-------------------------------------------------------------------------
      type(Object_type), intent(in) :: ObjectsList(:)
      integer,           intent(in) :: STLNum
      real(kind=RP),     intent(in) :: partitionScalar(:,:)      
#ifdef _HAS_MPI_      
      !-local-variables-------------------------------------------------------------------
      integer                    :: ierr, i, j, k, NumOfObjs, status(MPI_STATUS_SIZE), &
                                    array_of_statuses(MPI_STATUS_SIZE,4),              &
                                    send_Firstreq, send_req(4)
      real(kind=RP)              :: coords(NDIM)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state(:)

      if( MPI_Process% isRoot ) return

      NumOfObjs = 3*size(ObjectsList)

      call mpi_isend( NumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_Firstreq, ierr )

      call mpi_wait(send_Firstreq, status, ierr)

      allocate( COORD_x(NumOfObjs), &
                COORD_y(NumOfObjs), &
                COORD_z(NumOfObjs), &
                state(NumOfObjs)    )
                  
      k = 0
      do i = 1, size(ObjectsList)
         do j = 1, size(ObjectsList(i)% vertices)
            call OBB(STLNum)% ChangeRefFrame(ObjectsList(i)% vertices(j)% coords,GLOBAL,coords)
            k          = k + 1
            COORD_x(k) = coords(1)
            COORD_y(k) = coords(2)
            COORD_z(k) = coords(3)
            state(k)   = partitionScalar(i,j)
         end do
      end do

      call mpi_isend( COORD_x, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(1), ierr )
      
      call mpi_isend( COORD_y, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(2), ierr )
      
      call mpi_isend( COORD_z, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(3), ierr )
                         
      call mpi_isend( state, NumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_req(4), ierr )

      call mpi_waitall( 4, send_req, array_of_statuses, ierr )                            

      deallocate( COORD_x, COORD_y, COORD_z, state )     
#endif
   end subroutine sendScalarPlotRoot


   subroutine recvVectorPlotRoot( ObjectsList, STLNum, rootVector, x, y, z, vector_x, vector_y, vector_z )

      implicit none
      type(Object_type),          intent(in)    :: ObjectsList(:)
      integer,                    intent(in)    :: STLNum
      real(kind=RP),              intent(in)    :: rootVector(:,:,:)
      real(kind=RP), allocatable, intent(inout) :: x(:), y(:), z(:), vector_x(:), &
                                                   vector_y(:), vector_z(:)
      !-local-variables-------------------------------------------------------------------
      real(kind=RP)              :: coords(NDIM)
      integer                    :: i, j, k, rootNumOfObjs       
#ifdef _HAS_MPI_      
      integer                    :: ierr, nProcs, rank, ObjsSize, start_index,        &
                                    final_index, array_of_statuses(MPI_STATUS_SIZE,6)
      integer,       allocatable :: recv_req(:,:), recv_Firstreq(:), NumOfObjs(:)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state_x(:),   &
                                    state_y(:), state_z(:)
#endif
      if( .not. MPI_Process% isRoot ) return

      rootNumOfObjs = 3*size(ObjectsList)
#ifdef _HAS_MPI_ 

      allocate( NumOfObjs(MPI_Process% nProcs-1),    &
                recv_Firstreq(MPI_Process% nProcs-1) )

      NumOfObjs = 0

      do nProcs = 2, MPI_Process% nProcs
         
         call mpi_irecv( ObjsSize, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_Firstreq(nProcs-1), ierr ) 
         call mpi_wait( recv_Firstreq(nProcs-1), MPI_STATUS_IGNORE, ierr )

         NumOfObjs(nProcs-1) = ObjsSize

      end do 

      deallocate(recv_Firstreq)

      allocate( x(rootNumOfObjs+sum(NumOfObjs)),        &
                y(rootNumOfObjs+sum(NumOfObjs)),        &
                z(rootNumOfObjs+sum(NumOfObjs)),        &
                vector_x(rootNumOfObjs+sum(NumOfObjs)), &
                vector_y(rootNumOfObjs+sum(NumOfObjs)), &
                vector_z(rootNumOfObjs+sum(NumOfObjs)), &
                COORD_x(sum(NumOfObjs)),                &
                COORD_y(sum(NumOfObjs)),                &
                COORD_z(sum(NumOfObjs)),                &
                state_x(sum(NumOfObjs)),                &
                state_y(sum(NumOfObjs)),                &
                state_z(sum(NumOfObjs)),                & 
                recv_req(MPI_Process% nProcs-1,6)       )

      do nProcs = 2, MPI_Process% nProcs  

         start_index = sum(NumOfObjs(1:nProcs-2)) + 1 
         final_index = (start_index-1) + NumOfObjs(nProcs-1)

         call mpi_irecv( COORD_x(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,1), ierr )
         
         call mpi_irecv( COORD_y(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,2), ierr )
        
         call mpi_irecv( COORD_z(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,3), ierr )
        
         call mpi_irecv( state_x(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,4), ierr )
        
         call mpi_irecv( state_y(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,5), ierr )

         call mpi_irecv( state_z(start_index:final_index), NumOfObjs(nProcs-1), MPI_DOUBLE, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(nProcs-1,6), ierr )

         call mpi_waitall(6, recv_req(nProcs-1,:), array_of_statuses, ierr)    
           
      end do

      do i = 1, sum(NumOfObjs)
         x(rootNumOfObjs+i)        = COORD_x(i)
         y(rootNumOfObjs+i)        = COORD_y(i)
         z(rootNumOfObjs+i)        = COORD_z(i)
         vector_x(rootNumOfObjs+i) = state_x(i)
         vector_y(rootNumOfObjs+i) = state_y(i)
         vector_z(rootNumOfObjs+i) = state_z(i)
      end do

      deallocate( NumOfObjs, COORD_x, COORD_y, COORD_z, state_x, state_y, state_z, recv_req )
#else
      allocate( x(rootNumOfObjs),        &
                y(rootNumOfObjs),        &
                z(rootNumOfObjs),        &
                vector_x(rootNumOfObjs), &
                vector_y(rootNumOfObjs), &
                vector_z(rootNumOfObjs)  )
#endif  
      k = 0
      do i = 1, size(ObjectsList)
         do j = 1, size(ObjectsList(i)% vertices)
            call OBB(STLNum)% ChangeRefFrame(ObjectsList(i)% vertices(j)% coords,GLOBAL,coords)
            k           = k + 1
            x(k)        = coords(1)
            y(k)        = coords(2)
            z(k)        = coords(3)
            vector_x(k) = rootVector(1,i,j)
            vector_y(k) = rootVector(2,i,j)
            vector_z(k) = rootVector(3,i,j)
         end do
      end do

   end subroutine recvVectorPlotRoot

   subroutine sendVectorPlotRoot( ObjectsList, STLNum, partitionVector )
      implicit none 

       !-arguments-------------------------------------------------------------------------
      type(Object_type), intent(in) :: ObjectsList(:)
      integer,           intent(in) :: STLNum
      real(kind=RP),     intent(in) :: partitionVector(:,:,:)      
#ifdef _HAS_MPI_      
      !-local-variables-------------------------------------------------------------------
      integer                    :: ierr, i, j, k, NumOfObjs, status(MPI_STATUS_SIZE), &
                                    array_of_statuses(MPI_STATUS_SIZE,6),              &
                                    send_Firstreq, send_req(6)
      real(kind=RP)              :: coords(NDIM)
      real(kind=RP), allocatable :: COORD_x(:), COORD_y(:), COORD_z(:), state_x(:),    &
                                    state_y(:), state_z(:)

      if( MPI_Process% isRoot ) return

      NumOfObjs = 3*size(ObjectsList)

      call mpi_isend( NumOfObjs, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, send_Firstreq, ierr )

      call mpi_wait(send_Firstreq, status, ierr)

      allocate( COORD_x(NumOfObjs), &
                COORD_y(NumOfObjs), &
                COORD_z(NumOfObjs), &
                state_x(NumOfObjs), &
                state_y(NumOfObjs), &
                state_z(NumOfObjs)  )
                  
      k = 0
      do i = 1, size(ObjectsList)
         do j = 1, size(ObjectsList(i)% vertices)
            call OBB(STLNum)% ChangeRefFrame(ObjectsList(i)% vertices(j)% coords,GLOBAL,coords)
            k          = k + 1
            COORD_x(k) = coords(1)
            COORD_y(k) = coords(2)
            COORD_z(k) = coords(3)
            state_x(k) = partitionVector(1,i,j)
            state_y(k) = partitionVector(2,i,j)
            state_z(k) = partitionVector(3,i,j)
         end do
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
