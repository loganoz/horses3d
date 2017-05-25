module ProlongToFacesProcedures
   implicit none
   
   contains

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
         INTEGER :: N, i, j, k, nv
         
         N = e % N
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
         CALL InterpolateToBoundary( e % Q, spA % v(:,LEFT) , N, IX, e % Qb(:,:,:,ELEFT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % v(:,RIGHT), N, IX, e % Qb(:,:,:,ERIGHT), N_EQN)
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % v(:,FRONT), N, IY, e % Qb(:,:,:,EFRONT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % v(:,BACK) , N, IY, e % Qb(:,:,:,EBACK) , N_EQN)
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % v(:,BOTTOM), N, IZ, e % Qb(:,:,:,EBOTTOM) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % v(:,TOP)   , N, IZ, e % Qb(:,:,:,ETOP)    , N_EQN)

      END SUBROUTINE ProlongToFaces


   

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
         INTEGER :: N, i, j, k, nv
         
         N = e % N
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
         call InterpolateToBoundary( e % U_x , spA % v(:,LEFT ) , N , IX , e % U_xb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_x , spA % v(:,RIGHT) , N , IX , e % U_xb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % v(:,LEFT ) , N , IX , e % U_yb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % v(:,RIGHT) , N , IX , e % U_yb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % v(:,LEFT ) , N , IX , e % U_zb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % v(:,RIGHT) , N , IX , e % U_zb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % U_x , spA % v(:,FRONT) , N , IY , e % U_xb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x , spA % v(:,BACK)  , N , IY , e % U_xb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % v(:,FRONT) , N , IY , e % U_yb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % v(:,BACK)  , N , IY , e % U_yb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % v(:,FRONT) , N , IY , e % U_zb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % v(:,BACK)  , N , IY , e % U_zb(:,:,:,EBACK  ) , N_GRAD_EQN )
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % U_x, spA % v(:,BOTTOM), N, IZ , e % U_xb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x, spA % v(:,TOP)   , N, IZ , e % U_xb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % v(:,BOTTOM), N, IZ , e % U_yb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % v(:,TOP)   , N, IZ , e % U_yb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % v(:,BOTTOM), N, IZ , e % U_zb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % v(:,TOP)   , N, IZ , e % U_zb(:,:,:,ETOP)    , N_GRAD_EQN )                  

      END SUBROUTINE ProlongGradientToFaces   

      SUBROUTINE InterpolateToBoundary( u, v, N, which_dim , bValue , NEQ)
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
         INTEGER                      , INTENT(IN)  :: N
         real(kind=RP)                , intent(in)  :: u(0:,0:,0:,1:) , v(0:)
         integer                      , intent(in)  :: which_dim
         REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
         integer                      , intent(in)  :: NEQ
!        --------------------------------------------------------------------------------
         integer                                    :: i , j , k , eq

         select case (which_dim)
            case (IX)

               do eq=1,NEQ ; do j = 0,N ; do i = 0,N ; do k=0,N
                  bValue(eq,i,j) = bValue(eq,i,j) + u(k,i,j,eq) * v(k)
               end do ;      end do     ; end do     ; end do

            case (IY)

               do eq=1,NEQ ; do j = 0,N ; do k = 0,N ; do i=0,N
                  bValue(eq,i,j) = bValue(eq,i,j) + u(i,k,j,eq) * v(k)
               end do ;      end do     ; end do     ; end do

            case (IZ)

               do eq=1,NEQ ; do k = 0,N ; do j = 0,N ; do i=0,N
                  bValue(eq,i,j) = bValue(eq,i,j) + u(i,j,k,eq) * v(k)
               end do ;      end do     ; end do     ; end do

         end select
         
      END SUBROUTINE InterpolateToBoundary

end module ProlongToFacesProcedures
