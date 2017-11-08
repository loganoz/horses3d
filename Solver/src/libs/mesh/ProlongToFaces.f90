module ProlongToFacesProcedures
   implicit none
   
   CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongToFaces( e )
!
!     -----------------------------------------------------------
!     For Gauss point approximations, we interpolate to each face
!     of the element and store the result in the face solution 
!     array, Qb
!     -----------------------------------------------------------
!
         USE PhysicsStorage
         USE ElementClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(Element)      :: e
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: Nx, Ny, Nz, i, j, k, nv
         
         Nx = e % Nxyz(1)
         Ny = e % Nxyz(2)
         Nz = e % Nxyz(3)
!
!        --------------
!        Initialization
!        --------------
!
         e % storage % Qb = 0.0_RP
!
!        --------------
!        Left and right
!        --------------
!
         CALL InterpolateToBoundary(e % storage % Q, e % spAxi % v(:,LEFT) , Nx,Ny,Nz, IX, e % storage % Qb(:,:,:,ELEFT) , N_EQN)
         CALL InterpolateToBoundary(e % storage % Q, e % spAxi % v(:,RIGHT), Nx,Ny,Nz, IX, e % storage % Qb(:,:,:,ERIGHT), N_EQN)
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary(e % storage % Q, e % spAeta % v(:,FRONT), Nx,Ny,Nz, IY, e % storage % Qb(:,:,:,EFRONT), N_EQN)
         CALL InterpolateToBoundary(e % storage % Q, e % spAeta % v(:,BACK) , Nx,Ny,Nz, IY, e % storage % Qb(:,:,:,EBACK) , N_EQN)
!
!        --------------
!        Bottom and Top
!        --------------
!
        CALL InterpolateToBoundary(e % storage % Q, e % spAzeta % v(:,BOTTOM), Nx,Ny,Nz, IZ, e % storage % Qb(:,:,:,EBOTTOM), N_EQN)
        CALL InterpolateToBoundary(e % storage % Q, e % spAzeta % v(:,TOP)   , Nx,Ny,Nz, IZ, e % storage % Qb(:,:,:,ETOP)   , N_EQN)

      END SUBROUTINE ProlongToFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongGradientToFaces( e )
!
!     -----------------------------------------------------------
!     For Gauss point approximations, we interpolate to each face
!     of the element and store the result in the face solution 
!     array, Qb
!     -----------------------------------------------------------
!
         USE PhysicsStorage
         USE ElementClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(Element)      :: e
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: Nx, Ny, Nz
         INTEGER :: i, j, k, nv
         
         Nx = e % Nxyz(1)
         Ny = e % Nxyz(2)
         Nz = e % Nxyz(3)
!
!        --------------
!        Initialization
!        --------------
!  
         e % storage % U_xb = 0.0_RP
         e % storage % U_yb = 0.0_RP
         e % storage % U_zb = 0.0_RP

!
!        --------------
!        Left and Right
!        --------------
!
         call InterpolateToBoundary( e % storage % U_x , e % spAxi % v(:,LEFT ) , Nx, Ny, Nz , IX , e % storage % U_xb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % storage % U_x , e % spAxi % v(:,RIGHT) , Nx, Ny, Nz , IX , e % storage % U_xb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % storage % U_y , e % spAxi % v(:,LEFT ) , Nx, Ny, Nz , IX , e % storage % U_yb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % storage % U_y , e % spAxi % v(:,RIGHT) , Nx, Ny, Nz , IX , e % storage % U_yb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % storage % U_z , e % spAxi % v(:,LEFT ) , Nx, Ny, Nz , IX , e % storage % U_zb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % storage % U_z , e % spAxi % v(:,RIGHT) , Nx, Ny, Nz , IX , e % storage % U_zb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % storage % U_x , e % spAeta % v(:,FRONT) , Nx, Ny, Nz , IY , e % storage % U_xb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_x , e % spAeta % v(:,BACK)  , Nx, Ny, Nz , IY , e % storage % U_xb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_y , e % spAeta % v(:,FRONT) , Nx, Ny, Nz , IY , e % storage % U_yb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_y , e % spAeta % v(:,BACK)  , Nx, Ny, Nz , IY , e % storage % U_yb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_z , e % spAeta % v(:,FRONT) , Nx, Ny, Nz , IY , e % storage % U_zb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_z , e % spAeta % v(:,BACK)  , Nx, Ny, Nz , IY , e % storage % U_zb(:,:,:,EBACK  ) , N_GRAD_EQN )
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % storage % U_x, e % spAzeta % v(:,BOTTOM), Nx, Ny, Nz , IZ , e % storage % U_xb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_x, e % spAzeta % v(:,TOP)   , Nx, Ny, Nz , IZ , e % storage % U_xb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_y, e % spAzeta % v(:,BOTTOM), Nx, Ny, Nz , IZ , e % storage % U_yb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_y, e % spAzeta % v(:,TOP)   , Nx, Ny, Nz , IZ , e % storage % U_yb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_z, e % spAzeta % v(:,BOTTOM), Nx, Ny, Nz , IZ , e % storage % U_zb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % storage % U_z, e % spAzeta % v(:,TOP)   , Nx, Ny, Nz , IZ , e % storage % U_zb(:,:,:,ETOP)    , N_GRAD_EQN )                  

      END SUBROUTINE ProlongGradientToFaces 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE InterpolateToBoundary( u, v, Nx, Ny, Nz, which_dim , bValue , NEQ)
!
!     -------------------------------------------------------------
!     Interpolation to the boundary is a dot product for each row or
!     column. Using here the intrinsic Fortran function, without
!     having tested that it is faster or slower than a direct
!     computation for the values of N that we used in the DGSEM.
!     -------------------------------------------------------------
!
         USE SMConstants
         USE Physics
         IMPLICIT NONE
         INTEGER                      , INTENT(IN)    :: Nx, Ny, Nz
         real(kind=RP)                , intent(in)    :: u(1:NEQ,0:Nx,0:Ny,0:Nz) , v(0:)
         integer                      , intent(in)    :: which_dim
         REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
         integer                      , intent(in)    :: NEQ
!
!        ---------------
!        Local variables
!        ---------------
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
         
      END SUBROUTINE InterpolateToBoundary
      
end module ProlongToFacesProcedures
