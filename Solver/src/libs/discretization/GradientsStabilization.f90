#include "Includes.h"
module GradientsStabilization
#if defined(CAHNHILLIARD)
   use SMConstants
   use MeshTypes
   use HexMeshClass
   use PhysicsStorage
   use VariableConversion
   use BoundaryConditions, only: BCs
   implicit none

   private
   public   LambdaEstimator_f

   public   StabilizeGradients

   abstract interface
      subroutine LambdaEstimator_f(mesh, fID, i, j, lambda)
         use SMConstants
         use HexMeshClass
         implicit none
         class(HexMesh),   intent(in)  :: mesh
         integer,          intent(in)  :: fID, i, j
         real(kind=RP),    intent(out) :: lambda
      end subroutine LambdaEstimator_f
   end interface
!
!  ========
   contains
!  ========
!
      subroutine StabilizeGradients(mesh, time)
         implicit none
         class(HexMesh)                :: mesh
         real(kind=RP),    intent(in)  :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, fID, i, j, k
!
!        **************************************************
!        Compute the face stabilization flux in local faces
!        **************************************************
!
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               if (f % IsMortar==0) then 
               call GradientsStabilization_InteriorFace(f) 
               elseif (f % Ismortar==1) then 
               call GradientsStabilization_InteriorFace(f,mesh % faces(fID+1),mesh % faces(fID+2),mesh % faces(fID+3),mesh % faces(fID+4)) 
               end if 
            
            case (HMESH_BOUNDARY) 
               call GradientsStabilization_BoundaryFace(f, time) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
#ifdef _HAS_MPI_
            if ( e % hasSharedFaces ) cycle
#endif
!
!           Revert the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) * e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
!
!           Add the surface integrals
!           -------------------------
            call PerformSurfaceIntegrals( e, mesh)
!
!           Perform the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do

            end associate
         end do
!$omp end do

#ifdef _HAS_MPI_
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
               call GradientsStabilization_MPIFace(f) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            if ( .not. e % hasSharedFaces ) cycle
!
!           Revert the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) * e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do

!
!           Add the surface integrals
!           -------------------------
            call PerformSurfaceIntegrals( e, mesh)
!
!           Perform the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do

            end associate
         end do
!$omp end do
#endif
      end subroutine StabilizeGradients

      subroutine GradientsStabilization_InteriorFace(f, fma, fmb, fmc, fmd)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)             :: f
         type(Face),optional    :: fma
         type(Face),optional    :: fmb
         type(Face),optional    :: fmc
         type(Face),optional    :: fmd
         procedure(LambdaEstimator_f)  :: LambdaEstimator
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i,j, lm
         real(kind=RP) :: cL(NCOMP), cR(NCOMP), uHat(NCOMP)
         real(kind=RP) :: Hflux(NCOMP,NDIM,0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP) :: vAver(NDIM)

         real(kind=RP), allocatable :: HfluxM1(:,:,:,:)
         real(kind=RP), allocatable :: HfluxM2(:,:,:,:)
         real(kind=RP), allocatable :: HfluxM3(:,:,:,:)
         real(kind=RP), allocatable :: HfluxM4(:,:,:,:)
         integer :: Nfm(4,2) 

         if (f % IsMortar ==0 ) then 
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            cL = f % storage(1) % c(:,i,j)
            cR = f % storage(2) % c(:,i,j)
  
            vAver = AVERAGE( f % storage(1) % v(:,i,j) , f % storage(2) % v(:,i,j))

            Uhat = 0.5_RP * sign2(dot_product(vAver, f % geom % normal(:,i,j)))*(cL - cR) * f % geom % jacobian(i,j) 
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         call f % ProjectGradientFluxToElements(NCOMP, HFlux,(/1,2/),-1)
      end if 
      if (f % IsMortar ==1) then 
         do lm=1,4
            Nfm(1,:)=fma % Nf
            Nfm(2,:)=fmb % Nf
            Nfm(3,:)=fmc % Nf
            Nfm(4,:)=fmd % Nf
            if (lm==1) allocate(HfluxM1(NCOMP,NDIM,0:fma% Nf(1), 0:fma % Nf(2)))
            if (lm==2) allocate(HfluxM2(NCOMP,NDIM,0:fmb% Nf(1), 0:fmb % Nf(2)))
            if (lm==3) allocate(HfluxM2(NCOMP,NDIM,0:fmc% Nf(1), 0:fmc % Nf(2)))
            if (lm==4) allocate(HfluxM2(NCOMP,NDIM,0:fmd% Nf(1), 0:fmd % Nf(2)))

            do j = 0, Nfm(lm,2)  ; do i = 0, Nfm(lm,1)
               select case (lm)
               case(1)
               cL = fma % storage(1) % c(:,i,j)
               cR = fma % storage(2) % c(:,i,j)
   
               vAver = AVERAGE( fma % storage(1) % v(:,i,j) , fma % storage(2) % v(:,i,j))
               Uhat = 0.5_RP * sign2(dot_product(vAver, fma % geom % normal(:,i,j)))*(cL - cR) * fma % geom % jacobian(i,j) 
               case(2)
               cL = fmb % storage(1) % c(:,i,j)
               cR = fmb % storage(2) % c(:,i,j)
   
               vAver = AVERAGE( fmb % storage(1) % v(:,i,j) , fmb % storage(2) % v(:,i,j))
               Uhat = 0.5_RP * sign2(dot_product(vAver, fmb % geom % normal(:,i,j)))*(cL - cR) * fmb % geom % jacobian(i,j)
               case(3)
               cL = fmc % storage(1) % c(:,i,j)
               cR = fmc % storage(2) % c(:,i,j)
   
               vAver = AVERAGE( fmc % storage(1) % v(:,i,j) , fmc % storage(2) % v(:,i,j))
               Uhat = 0.5_RP * sign2(dot_product(vAver, fmc % geom % normal(:,i,j)))*(cL - cR) * fmc % geom % jacobian(i,j)
               case(4)
               cL = fmd % storage(1) % c(:,i,j)
               cR = fmd % storage(2) % c(:,i,j)
   
               vAver = AVERAGE( fmd % storage(1) % v(:,i,j) , fmd % storage(2) % v(:,i,j))
               Uhat = 0.5_RP * sign2(dot_product(vAver, fmd % geom % normal(:,i,j)))*(cL - cR) * fmd % geom % jacobian(i,j)
               end select 

               select case (lm)
               case (1)
                  HfluxM1(:,IX,i,j) = Uhat * fma % geom % normal(IX,i,j)
                  HfluxM1(:,IY,i,j) = Uhat * fma % geom % normal(IY,i,j)
                  HfluxM1(:,IZ,i,j) = Uhat * fma % geom % normal(IZ,i,j)
               case (2)
                  HfluxM2(:,IX,i,j) = Uhat * fmb % geom % normal(IX,i,j)
                  HfluxM2(:,IY,i,j) = Uhat * fmb % geom % normal(IY,i,j)
                  HfluxM2(:,IZ,i,j) = Uhat * fmb % geom % normal(IZ,i,j)
               case(3)
                  HfluxM3(:,IX,i,j) = Uhat * fmc % geom % normal(IX,i,j)
                  HfluxM3(:,IY,i,j) = Uhat * fmc % geom % normal(IY,i,j)
                  HfluxM3(:,IZ,i,j) = Uhat * fmc % geom % normal(IZ,i,j)
               case(4)
                  HfluxM4(:,IX,i,j) = Uhat * fmd % geom % normal(IX,i,j)
                  HfluxM4(:,IY,i,j) = Uhat * fmd % geom % normal(IY,i,j)
                  HfluxM4(:,IZ,i,j) = Uhat * fmd % geom % normal(IZ,i,j)
               end select 
            end do               ; end do
         end do 
         call f % ProjectMortarGradientFluxToElements(NCOMP, HFlux,(/1,0/),-1, fma, fmb, fmc, fmd, HfluxM1, HfluxM2, HfluxM3, HfluxM4)
         call fma % ProjectGradientFluxToElements(NCOMP, HfluxM1,(/0,2/),-1)
         call fmb % ProjectGradientFluxToElements(NCOMP, HfluxM2,(/0,2/),-1)
         call fmc % ProjectGradientFluxToElements(NCOMP, HfluxM3,(/0,2/),-1)
         call fmd % ProjectGradientFluxToElements(NCOMP, HfluxM4,(/0,2/),-1)
         deallocate(HfluxM1)
         deallocate(HfluxM2)
         deallocate(HfluxM2)
         deallocate(HfluxM2)
      end if 
         
      end subroutine GradientsStabilization_InteriorFace   

      subroutine GradientsStabilization_MPIFace(f)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)    :: f
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i,j, thisSide
         real(kind=RP) :: cL(NCOMP), cR(NCOMP), Uhat(NCOMP)
         real(kind=RP) :: Hflux(NCOMP,NDIM,0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP) :: vAver(NDIM)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            cL = f % storage(1) % c(:,i,j)
            cR = f % storage(2) % c(:,i,j)
   
            vAver = AVERAGE( f % storage(1) % v(:,i,j) , f % storage(2) % v(:,i,j))

            Uhat = 0.5_RP * sign2(dot_product(vAver, f % geom % normal(:,i,j)))*(cL - cR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(NCOMP, HFlux, (/thisSide, HMESH_NONE/), -1)
         
      end subroutine GradientsStabilization_MPIFace   

      subroutine GradientsStabilization_BoundaryFace(f, time)
         use Physics
         use FaceClass
         implicit none
         type(Face)             :: f
         real(kind=RP)          :: time
         integer                :: i, j
         real(kind=RP)          :: Uhat(NCOMP), cL(NCOMP), cR(NCOMP)
         real(kind=RP)          :: bvExt(NCOMP)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

            bvExt =  f % storage(1) % c(:,i,j)
   
            call BCs(f % zone) % bc % PhaseFieldState( f % geom % x(:,i,j), & 
                                time               , &
                                f % geom % normal(:,i,j)      , &
                                bvExt              )
!   
!           -------------------
!           u, v, w, T averages
!           -------------------
!   
            cL = f % storage(1) % Q(:,i,j)
            cR = bvExt
   
            Uhat = 0.5_RP * dot_product(f % storage(1) % v(:,i,j), f % geom % normal(:,i,j))*(cL - cR) * f % geom % jacobian(i,j) / (abs(dot_product(f % storage(1) % v(:,i,j),f % geom % normal(:,i,j))) + epsilon(1.0_RP))
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine GradientsStabilization_BoundaryFace

      subroutine PerformSurfaceIntegrals( e, mesh )
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         class(Element)                      :: e
         class(HexMesh)                      :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: faceInt_x(NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_y(NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_z(NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )

         call VectorWeakIntegrals % StdFace(e, NCOMP, &
               mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % unStar, &
               mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % unStar, &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar, &
               mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % unStar, &
               mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % unStar, &
               mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % unStar, &
               faceInt_x, faceInt_y, faceInt_z )

         e % storage % c_x = e % storage % c_x + faceInt_x
         e % storage % c_y = e % storage % c_y + faceInt_y
         e % storage % c_z = e % storage % c_z + faceInt_z

      end subroutine PerformSurfaceIntegrals
#endif
end module GradientsStabilization