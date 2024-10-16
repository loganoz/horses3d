#include "Includes.h"
!  **************
   module Physics_SLR
!  **************
!
      use SMConstants
      use PhysicsStorage_SLR
      use VariableConversion_SLR
      use FluidData_SLR
      implicit none
      
      private
      public  slrViscousFlux, ViscousJacobian

!
!     ========
      CONTAINS 
!     ========
!
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine slrViscousFlux(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn     )
         real(kind=RP), intent(in)  :: U_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)

         ! Change the indecies to slr
         F(1,IX) = mu * U_x(1)

         F(1,IY) = mu * U_y(1)

         F(1,IZ) = mu * U_z(1)


      end subroutine slrViscousFlux
 
!     -------------------------------------------------------------------------------
!     Subroutine for computing the Jacobian of the inviscid flux when it has the form 
!
!        F = f*iHat + g*jHat + h*kHat
!
!     First index indicates the flux term and second index indicates the conserved 
!     variable term. For example:
!           dfdq     := df/dq
!                       d f(2) |
!           dfdq(2,4) = ------ |
!                       d q(4) |q
!     ***** This routine is necessary for computing the analytical Jacobian. *****
!     -------------------------------------------------------------------------------
      pure subroutine InviscidJacobian(q,dfdq,dgdq,dhdq)
         implicit none
         !-------------------------------------------------
         real(kind=RP), intent (in)  :: q(NCONS)
         real(kind=RP), intent (out) :: dfdq(NCONS,NCONS)
         real(kind=RP), intent (out) :: dgdq(NCONS,NCONS)
         real(kind=RP), intent (out) :: dhdq(NCONS,NCONS)
         !-------------------------------------------------
      !    real(kind=RP)  :: u,v,w ! Velocity components
      !    real(kind=RP)  :: V2    ! Total velocity squared
      !    real(kind=RP)  :: p     ! Pressure
      !    real(kind=RP)  :: H     ! Total enthalpy
      !    !-------------------------------------------------
         
      !    associate( gammaMinus1 => thermodynamics % gammaMinus1, & 
      !               gamma => thermodynamics % gamma )
         
      !    u  = q(IRHOU) / q(IRHO)
      !    v  = q(IRHOV) / q(IRHO)
      !    w  = q(IRHOW) / q(IRHO)
      !    V2 = u*u + v*v + w*w
      !    p  = Pressure(q)
      !    H  = (q(IRHOE) + p) / q(IRHO)
!
!        Flux in the x direction (f)
!        ---------------------------

         dfdq(1,1) = 0._RP
         
!
!        Flux in the y direction (g)
!        ---------------------------
         
         dgdq(1,1) = 0._RP
!
!        Flux in the z direction (h)
!        ---------------------------
         
         dhdq(1,1) = 0._RP
         
         !end associate
         
      end subroutine InviscidJacobian

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!     ------------------------------------------------------------------------------------------
!     Subroutine for computing the Jacobians of the viscous fluxes 
!
!     1. Jacobian with respect to the gradients of the conserved variables: df/d(∇q)
!        Every direction of the viscous flux can be expressed as
!
!        f_i = \sum_j df_dgradq(:,:,j,i) dq/dx_j
!
!        In Navier-Stokes, this dependence is linear. Therefore:
!
!              df_dgradq         := df/d(∇q)
!                                   d f_i(2)   |
!              df_dgradq(2,4,j,i) = ---------- |
!                                   d(∇q)_j(4) |q=cons,
!        where (∇q)_j = dq/dx_j
!        
!        Following Hartmann's notation, G_{ij} = df_dgradq(:,:,j,i). --> R. Hartmann. "Discontinuous Galerkin methods for compressible flows: higher order accuracy, error estimation and adaptivity". 2005.
!
!     2. Jacobian with respect to the conserved variables: df/dq
!
!              df_dq       := df/d(∇q)
!                             d f_i(2) |
!              df_dq(2,4,i) = -------- |
!                             dq(4)    |∇q=cons
!
!     NOTE 1: Here the thermal conductivity and the viscosity are computed using Sutherland's law!     
!     NOTE 2: The dependence of the temperature on q is not considered in Sutherland's law
!
!     ***** This routine is necessary for computing the analytical Jacobian. *****
!     ------------------------------------------------------------------------------------------
      !pure subroutine ViscousJacobian(q, Q_x, Q_y, Q_z, df_dgradq, df_dq)
       subroutine ViscousJacobian(q, Q_x, Q_y, Q_z, df_dgradq, df_dq)
         implicit none
         !-------------------------------------------------
         real(kind=RP), intent(in)  :: q(NCONS)                      !< Conserved variables state
         real(kind=RP), intent(in)  :: Q_x (1:NGRAD)
         real(kind=RP), intent(in)  :: Q_y (1:NGRAD)
         real(kind=RP), intent(in)  :: Q_z (1:NGRAD) ! , intent(in)
         real(kind=RP), intent(out) :: df_dgradq(NCONS,NCONS,NDIM,NDIM)
         real(kind=RP), intent(out) :: df_dq    (NCONS,NCONS,NDIM)
         !-------------------------------------------------
!          real(kind=RP)            :: T , sutherLaw
!          real(kind=RP)            :: u , v , w, E, u2, v2, w2, Vel2
!          real(kind=RP)            :: vv_x, vv_y, vv_z
!          real(kind=RP)            :: V_gradU, V_gradV, V_gradW, gradE(3)
!          real(kind=RP)            :: gamma_Pr
!          real(kind=RP)            :: rho_DivV, V_gradRho
!          real(kind=RP)            :: invRho, invRho2, uDivRho(NDIM), U_x(NDIM), U_y(NDIM), U_z(NDIM)
!          real(kind=RP)            :: F(NCONS,NDIM)
!          real(kind=RP)            :: dMu_dQ(NCONS)
!          real(kind=RP), parameter :: lambda = 1._RP/3._RP
!          real(kind=RP), parameter :: f4_3 = 4._RP/3._RP
!          real(kind=RP), parameter :: f2_3 = 2._RP/3._RP
!          !-------------------------------------------------
         
!          invRho  = 1._RP / Q(IRHO)
!          invRho2 = invRho * invRho
         
!          uDivRho = [Q(IRHOU) , Q(IRHOV) , Q(IRHOW) ] * invRho2
         
!          u_x = invRho * Q_x(IRHOU:IRHOW) - uDivRho * Q_x(IRHO)
!          u_y = invRho * Q_y(IRHOU:IRHOW) - uDivRho * Q_y(IRHO)
!          u_z = invRho * Q_z(IRHOU:IRHOW) - uDivRho * Q_z(IRHO)
         
!          u  = Q(IRHOU) * invRho
!          v  = Q(IRHOV) * invRho
!          w  = Q(IRHOW) * invRho
         
!          E  = Q(IRHOE) * invRho
!          u2 = u*u
!          v2 = v*v
!          w2 = w*w
!          Vel2 = u2 + v2 + w2

!          T     = Temperature(q)
!          sutherLaw = SutherlandsLaw(T)
         
!          associate ( gamma => thermodynamics % gamma, & 
!                      gammaM2 => dimensionless % gammaM2, &
!                      gammaMinus1 => thermodynamics % gammaMinus1, &
!                      Re    => dimensionless % Re    , &
!                      Pr    => dimensionless % Pr ) 
         
!          gamma_Pr = gamma/Pr
         
!
!        *****************************
!        Derivative with respect to ∇q
!        *****************************
!
         
!
!        Flux in the x direction: f = G_{1:} · ∇q
!        ----------------------------------------
           
            
!          ! G_{11}
           df_dgradq(:,1,1,1) = (/ 1.0_RP /)
         
!          ! G_{12}
           df_dgradq(:,1,2,1) = ( 0.0_RP )
         
!          ! G_{13}
           df_dgradq(:,1,3,1) = ( 0.0_RP)

         
!
!        Flux in the y direction: g = G_{2:} · ∇q
!        ----------------------------------------
         
         ! G_{21}
         df_dgradq(:,1,1,2) = ( 0.0_RP)

         
!          ! G_{22}
          df_dgradq(:,1,2,2) = (1.0_RP)

         
!          ! G_{23}
          df_dgradq(:,1,3,2) = ( 0.0_RP)
         
!
!        Flux in the z direction: h = G_{3:} · ∇q
!        ----------------------------------------
         
!          ! G_{31}
          df_dgradq(:,1,1,3) = ( 0.0_RP)
         
!          ! G_{32}
          df_dgradq(:,1,2,3) = ( 0.0_RP )
         
!          ! G_{33}
          df_dgradq(:,1,3,3) = ( 1.0_RP )

         
! !
! !        Scale with mu/(rho*Re) .or. kappa/(rho*Re)
! !        ------------------------------------------
         
!          df_dgradq = df_dgradq * sutherLaw / ( Q(IRHO) * Re )
         
! !
! !        ****************************
! !        Derivative with respect to q
! !        ****************************
! !
! !        Auxiliary variables
! !        ------------------
         
!          rho_DivV      = Q(IRHO) * ( U_x(IX) + U_y(IY) + U_z(IZ) )       ! rho ∇ · v
!          V_gradRho     = u * Q_x(IRHO) + v * Q_y(IRHO) + w * Q_z(IRHO)   ! v · ∇rho
!          V_gradU       = u * U_x(IX) + v * U_y(IX) + w * U_z(IX)
!          V_gradV       = u * U_x(IY) + v * U_y(IY) + w * U_z(IY)
!          V_gradW       = u * U_x(IZ) + v * U_y(IZ) + w * U_z(IZ)
         
!          vv_x = 2._RP * Q(IRHO) * ( u * U_x(IX) + v * U_x(IY) + w * U_x(IZ) )
!          vv_y = 2._RP * Q(IRHO) * ( u * U_y(IX) + v * U_y(IY) + w * U_y(IZ) )
!          vv_z = 2._RP * Q(IRHO) * ( u * U_z(IX) + v * U_z(IY) + w * U_z(IZ) )
         
!          gradE(1) = Q_x(IRHOE) * invRho - E * invRho * Q_x(IRHO)
!          gradE(2) = Q_y(IRHOE) * invRho - E * invRho * Q_y(IRHO)
!          gradE(3) = Q_z(IRHOE) * invRho - E * invRho * Q_z(IRHO)
         
! !        Jacobian entries
! !        ----------------
         
!          ! A_1
         df_dq(:,1,1) = (0.0_RP)
         
!          ! A_2
          df_dq(:,1,2) = (0.0_RP)
         
         
!          ! A_3
         df_dq(:,1,3) = (0.0_RP)
         
         
! !
! !        Scale with mu/(rho² Re) .or. kappa/(rho² Re)
! !        --------------------------------------------
         
!          df_dq = df_dq * sutherLaw / ( Q(IRHO)**2 * Re ) 
         
! !
! !        Correct with the derivative of the Sutherland's law
! !        --------------------------------------------------
!          dMu_dQ    = SutherlandsLawDeriv(Q,T)
         
!          call ViscousFlux_STATE(NCONS, NCONS, Q, Q_x, Q_y, Q_z, dimensionless % mu, 0._RP, dimensionless % kappa, F)
!          F = F / sutherLaw
         
!          df_dq(:,:,1) = df_dq(:,:,1) + outer_product(F(:,1),dMu_dQ)
!          df_dq(:,:,2) = df_dq(:,:,2) + outer_product(F(:,2),dMu_dQ)
!          df_dq(:,:,3) = df_dq(:,:,3) + outer_product(F(:,3),dMu_dQ)
         
!          end associate
      end subroutine ViscousJacobian
   END Module Physics_SLR
