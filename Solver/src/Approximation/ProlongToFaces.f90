module ProlongToFacesProcedures
   implicit none
   
   CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongToFaces( e, spA )
!
!     -----------------------------------------------------------
!     For Gauss point approximations, we interpolate to each face
!     of the element and store the result in the face solution 
!     array, Qb
!     -----------------------------------------------------------
!
         USE PhysicsStorage
         USE NodalStorageClass
         USE ElementClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(NodalStorage) :: spA
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
         e % Qb = 0.0_RP
!
!        --------------
!        Left and right
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % vx(:,LEFT) , Nx, Ny, Nz, IX, e % Qb(:,:,:,ELEFT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % vx(:,RIGHT), Nx, Ny, Nz, IX, e % Qb(:,:,:,ERIGHT), N_EQN)
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % vy(:,FRONT), Nx, Ny, Nz, IY, e % Qb(:,:,:,EFRONT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % vy(:,BACK) , Nx, Ny, Nz, IY, e % Qb(:,:,:,EBACK)  , N_EQN)
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % vz(:,BOTTOM), Nx, Ny, Nz, IZ, e % Qb(:,:,:,EBOTTOM) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % vz(:,TOP)   , Nx, Ny, Nz, IZ, e % Qb(:,:,:,ETOP)    , N_EQN)

      END SUBROUTINE ProlongToFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongGradientToFaces( e, spA )
!
!     -----------------------------------------------------------
!     For Gauss point approximations, we interpolate to each face
!     of the element and store the result in the face solution 
!     array, Qb
!     -----------------------------------------------------------
!
         USE PhysicsStorage
         USE NodalStorageClass
         USE ElementClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(NodalStorage) :: spA
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
         e % U_xb = 0.0_RP
         e % U_yb = 0.0_RP
         e % U_zb = 0.0_RP

!
!        --------------
!        Left and Right
!        --------------
!
         call InterpolateToBoundary( e % U_x , spA % vx(:,LEFT ) , Nx, Ny, Nz , IX , e % U_xb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_x , spA % vx(:,RIGHT) , Nx, Ny, Nz , IX , e % U_xb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % vx(:,LEFT ) , Nx, Ny, Nz , IX , e % U_yb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % vx(:,RIGHT) , Nx, Ny, Nz , IX , e % U_yb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % vx(:,LEFT ) , Nx, Ny, Nz , IX , e % U_zb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % vx(:,RIGHT) , Nx, Ny, Nz , IX , e % U_zb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % U_x , spA % vy(:,FRONT) , Nx, Ny, Nz , IY , e % U_xb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x , spA % vy(:,BACK)  , Nx, Ny, Nz , IY , e % U_xb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % vy(:,FRONT) , Nx, Ny, Nz , IY , e % U_yb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % vy(:,BACK)  , Nx, Ny, Nz , IY , e % U_yb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % vy(:,FRONT) , Nx, Ny, Nz , IY , e % U_zb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % vy(:,BACK)  , Nx, Ny, Nz , IY , e % U_zb(:,:,:,EBACK  ) , N_GRAD_EQN )
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % U_x, spA % vz(:,BOTTOM), Nx, Ny, Nz , IZ , e % U_xb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x, spA % vz(:,TOP)   , Nx, Ny, Nz , IZ , e % U_xb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % vz(:,BOTTOM), Nx, Ny, Nz , IZ , e % U_yb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % vz(:,TOP)   , Nx, Ny, Nz , IZ , e % U_yb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % vz(:,BOTTOM), Nx, Ny, Nz , IZ , e % U_zb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % vz(:,TOP)   , Nx, Ny, Nz , IZ , e % U_zb(:,:,:,ETOP)    , N_GRAD_EQN )                  

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
         real(kind=RP)                , intent(in)    :: u(0:Nx,0:Ny,0:Nz,1:NEQ) , v(0:)
         integer                      , intent(in)    :: which_dim
         REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
         integer                      , intent(in)    :: NEQ
!        --------------------------------------------------------------------------------
         integer                                    :: i , j , k , eq

         select case (which_dim)
            case (IX)

               do eq = 1 , NEQ
                  do j = 0 , Nz
                     do i = 0 , Ny
                        do k = 0 , Nx

                           bValue(eq,i,j) = bValue(eq,i,j) + u(k,i,j,eq) * v(k)
   
                        end do
                     end do
                  end do
               end do

            case (IY)

                do eq = 1 , NEQ
                  do j = 0 , Nz
                     do k = 0 , Ny
                        do i = 0 , Nx

                           bValue(eq,i,j) = bValue(eq,i,j) + u(i,k,j,eq) * v(k)
   
                        end do
                     end do
                  end do
               end do

            case (IZ)

               do eq = 1 , NEQ
                  do k = 0 , Nz
                     do j = 0 , Ny
                        do i = 0 , Nx

                           bValue(eq,i,j) = bValue(eq,i,j) + u(i,j,k,eq) * v(k)
   
                        end do
                     end do
                  end do
               end do

         end select
         
      END SUBROUTINE InterpolateToBoundary
      
end module ProlongToFacesProcedures
