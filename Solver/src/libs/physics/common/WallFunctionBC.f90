#if defined(NAVIERSTOKES)
MODULE WallFunctionBC

   USE SMConstants
   USE PhysicsStorage
   IMPLICIT NONE

!   
!  *****************************
!  Default everything to private
!  *****************************
!
   PRIVATE
!
!  ****************
!  Public variables
!  ****************
!
!
!  ******************
!  Public definitions
!  ******************
!
   PUBLIC WallViscousFlux, wall_shear, u_tau_f, u_tau_f_ABL, y_plus_f, u_plus_f, WallFUnctionBC_FlowNeumann_HOIBM
!

   CONTAINS 
!   
!------------------------------------------------------------------------------------------------------------------------
!
      ! Wall shear stress is computed here using the Reichardt model
      ! following Frere et al 2017.
      ! Wall shear stress (tau_w) is then used to compute the viscous flux.
      ! The viscous flux is set in the same direction as the parallel velocity. 
      ! If the reference velocity parallel to the wall is less than a 
      ! minimum value, flux is set to zero to avoid numerical issues. 

   SUBROUTINE WallViscousFlux (U_inst, dWall, nHat, rho, mu, U_avg, visc_flux, u_tau)

      use WallFunctionDefinitions, only: useAverageV

      IMPLICIT NONE

      REAL(kind=RP) , INTENT(IN)     :: U_inst(NDIM)       ! Instantaneous velocity from LES solver
      REAL(kind=RP) , INTENT(IN)     :: dWall              ! Normal wall distance
      REAL(kind=RP),  INTENT(IN)     :: nHat(NDIM)         ! Unitary vector normal to wall
      REAL(kind=RP) , INTENT(IN)     :: rho                ! Density
      REAL(kind=RP) , INTENT(IN)     :: mu                 ! Dynamic viscosity
      REAL(kind=RP) , INTENT(IN)     :: U_avg(NDIM)        ! Mean (in time) velocity from LES solver
      REAL(kind=RP) , INTENT(INOUT)  :: visc_flux(NCONS)   ! Viscous Flux
      REAL(kind=RP) , INTENT(INOUT)  :: u_tau              ! friction velocity

      !local variables
      REAL(kind=RP), DIMENSION(NDIM) :: x_II               ! Unitary vector parallel to wall
      REAL(kind=RP), DIMENSION(NDIM) :: U_ref              ! Reference Velocity
      REAL(kind=RP), DIMENSION(NDIM) :: u_parallel         ! Velocity parallel to wall
      REAL(kind=RP), DIMENSION(NDIM) :: u_parallel_aux     ! Velocity parallel to wall auxiliary, is the instantaneous when the average is used
      REAL(kind=RP)                  :: u_II               ! Velocity magnitude parallel to wall, used in Wall Function
      REAL(kind=RP)                  :: tau_w              ! Wall shear stress
      REAL(kind=RP)                  :: beta               ! damping factor from Thomas et. al

      if (useAverageV) then
          U_ref = U_avg
      else
          U_ref = U_inst
      end if 

      u_parallel = U_ref - (dot_product(U_ref, nHat) * nHat)
      x_II = u_parallel / norm2(u_parallel)
      u_II = dot_product(U_ref, x_II)

      ! Wall model only modifies momentum viscous fluxes. 
      call wall_shear(u_II, dWall, rho, mu, tau_w, u_tau)
      ! visc_flux(IRHOU:IRHOW) = - tau_w_f (u_II,dWall,rho,mu) * x_II 
      if (useAverageV) then
          u_parallel_aux = U_inst - (dot_product(U_inst, nHat) * nHat)
          x_II = u_parallel_aux / u_II !schuman, the direction scales with the instantaneous values
          ! beta = 0.3_RP ! thomas arbitrary value. If set to 0.0_RP, the schuman eq is recovered
          ! x_II = (beta * u_parallel_aux + (1-beta) * u_parallel) / u_II ! thomas, the direction scales with both the instantaneous and the mean
      end if 
      visc_flux(IRHOU:IRHOW) = - tau_w * x_II 
      ! print *, "visc_flux: ", visc_flux

   END SUBROUTINE 
!   
!------------------------------------------------------------------------------------------------------------------------
!
   ! FUNCTION tau_w_f (u_II, y, rho, mu)
   SUBROUTINE wall_shear(u_II, y, rho, mu, tau_w, u_tau)
      
      USE WallFunctionDefinitions, ONLY: wallFuncIndex, STD_WALL, ABL_WALL
      IMPLICIT NONE

      REAL(kind=RP), INTENT(IN)     :: u_II    ! Mean streamwise velocity parallel to the wall
      REAL(kind=RP), INTENT(IN)     :: y       ! Normal wall distance
      REAL(kind=RP), INTENT(IN)     :: rho     ! Density
      REAL(kind=RP), INTENT(IN)     :: mu      ! Dynamic viscosity

      REAL(kind=RP), INTENT(OUT)    :: tau_w   ! (OUT) Wall shear stress
      REAL(kind=RP), INTENT(INOUT)  :: u_tau   ! Friction velocity
      ! REAL(kind=RP)              :: tau_w_f ! (OUT) Wall shear stress

      ! REAL(kind=RP)              :: u_tau   ! Friction velocity
      REAL(kind=RP)                 :: nu      ! Kinematic viscosity

      nu = mu / rho 

      select case (wallFuncIndex)
          case (STD_WALL)
              ! u_tau is computed by solving Eq. (3) in Frere et al 2017
              ! along with the definitions of u+ and y+.
              ! use previous solution as the new starting point
              u_tau = u_tau_f( u_II, y, nu, u_tau )
          case (ABL_WALL)
              ! print *, "u_II: ", u_II, " z: ", y
              u_tau = u_tau_f_ABL( u_II, y, nu )
              ! print *, "u_tau: ", u_tau
      end select

      ! then the definition of the wall shear stress is used
      tau_w = rho * u_tau * u_tau

   END SUBROUTINE
   ! END FUNCTION  
!   
!------------------------------------------------------------------------------------------------------------------------
!
   FUNCTION u_tau_f (u_II,y,nu, u_tau0)

      USE WallFunctionDefinitions, ONLY: newtonTol, newtonAlpha, newtonMaxIter
      ! USE WallFunctionDefinitions, ONLY: newtonTol, newtonAlpha, newtonMaxIter, u_tau0
      IMPLICIT NONE

      REAL(kind=RP), INTENT(IN)        :: u_II         ! Mean streamwise velocity parallel to the wall
      REAL(kind=RP), INTENT(IN)        :: y            ! Normal wall distance
      REAL(kind=RP), INTENT(IN)        :: nu           ! Kinematic viscosity   
      REAL(kind=RP), INTENT(IN)        :: u_tau0       ! 

      REAL(kind=RP)                    :: u_tau_f      ! Friction velocity
      

      REAL(kind=RP)                    :: u_tau            ! Previous value for Newton method
      REAL(kind=RP)                    :: u_tau_next       ! Next value for Newton method
      INTEGER                          :: i                ! Counter
      
      REAL(kind=RP)                    :: Aux_x0           ! Objective function evaluated at x0
      REAL(kind=RP)                    :: JAC              ! Derivate of objective function evaluated at x0
      REAL(kind=RP)                    :: eps              ! Size of the perturbation to compute numerical der.
      REAL(kind=RP)                    :: alpha            ! Parameter for the damped Newton method

      ! The value of u_tau is found by solving a non linear equation.
      ! The damped Newton method is used. 

      ! Initial seed for Newton method
      u_tau = u_tau0

      ! Iterate in Newton's method until convergence criteria is met

      DO i = 1, newtonMaxIter
         
         ! Evaluate auxiliary function at u_tau
         Aux_x0  =   Aux_f ( u_tau      , u_II, y, nu )

         ! Compute numerical derivative of auxiliary function at u_tau
         eps     = ABS(u_tau) * 1.0E-8_RP
         JAC     = ( Aux_f ( u_tau + eps, u_II, y, nu ) - Aux_x0 ) / eps

         ! Default value for alpha (Newton method)
         alpha = newtonAlpha
         ! Damped alpha parameter for the Damped Newton method 
         DO WHILE ( ABS( Aux_x0 ) < ABS( Aux_f(u_tau - Aux_x0 / JAC * alpha, u_II, y, nu ) ) )
            alpha = alpha / 2
         END DO  

         ! Damped Newton step
         u_tau_next = u_tau - Aux_x0 / JAC * alpha 

         ! Convergence criteria for Newton's method
         IF ( ( ABS ( ( u_tau_next - u_tau ) / u_tau ) < newtonTol ) &
                           .AND. &
              ( ABS ( Aux_x0 )                         < newtonTol ) ) THEN

            ! Assign output value to u_tau
            u_tau_f = u_tau_next  
            RETURN

         END IF

         ! Set value for u_tau for next iteration
         u_tau = u_tau_next

      END DO

      error stop "DAMPED NEWTON METHOD IN WALL FUNCTION DOES NOT CONVERGE."

   END FUNCTION 
!   
!------------------------------------------------------------------------------------------------------------------------
!
   FUNCTION Aux_f (u_tau, u_II, y, nu)

      IMPLICIT NONE

      REAL(kind=RP), INTENT(IN)  :: u_tau
      REAL(kind=RP), INTENT(IN)  :: u_II
      REAL(kind=RP), INTENT(IN)  :: y
      REAL(kind=RP), INTENT(IN)  :: nu

      REAL(kind=RP)              :: Aux_f ! (OUT)
        
      ! Auxiliary function is evaluated at x0
      ! When Aux_f = 0 The definition of the 
      ! dimensionless mean streamwise velocity 
      ! parallel to the wall is recovered and 
      ! a valid value for u_tau is found. 
      ! Aux_f = 0 -> u_II / u_tau = u_plus

      Aux_f = u_II - u_plus_f ( y_plus_f ( y, u_tau, nu ) ) * u_tau 

   END FUNCTION
!   
!------------------------------------------------------------------------------------------------------------------------
!
   PURE FUNCTION u_plus_f (y_plus)
       USE WallFunctionDefinitions, ONLY: kappa, WallC
      IMPLICIT NONE
      ! Definition of u_plus
      ! Reichardt law-of-the-wall (taken from Frere et al 2017 Eq. (3))
      REAL(kind=RP), INTENT(IN)  :: y_plus
      
      REAL(kind=RP)              :: u_plus_f ! (OUT)
      
      u_plus_f = 1.0_RP / kappa * log( 1.0_RP + kappa * y_plus ) &
                 + ( WallC - 1.0_RP / kappa * log(kappa) ) * & 
                   ( 1.0_RP - exp( -y_plus / 11.0_RP )-( y_plus / 11.0_RP ) * exp( -y_plus / 3.0_RP ) )

   END FUNCTION
!   
!------------------------------------------------------------------------------------------------------------------------
!
   PURE FUNCTION y_plus_f (y, u_tau, nu)
      ! Definition of y_plus (taken from Frere et al 2017)

      REAL(kind=RP), INTENT(IN)  :: y
      REAL(kind=RP), INTENT(IN)  :: u_tau
      REAL(kind=RP), INTENT(IN)  :: nu

      REAL(kind=RP)              :: y_plus_f ! (OUT)

      y_plus_f = y * u_tau / nu 

   END FUNCTION
!   
!------------------------------------------------------------------------------------------------------------------------
!
   FUNCTION u_tau_f_ABL (u_II,y,nu)
      USE WallFunctionDefinitions, ONLY: y0, d, kappa
      IMPLICIT NONE

      REAL(kind=RP), INTENT(IN)        :: u_II         ! Mean streamwise velocity parallel to the wall
      REAL(kind=RP), INTENT(IN)        :: y            ! Normal wall distance
      REAL(kind=RP), INTENT(IN)        :: nu           ! Kinematic viscosity   
      REAL(kind=RP)                    :: u_tau_f_ABL  ! (OUT) Friction velocity

      u_tau_f_ABL = kappa * u_II / log( (y-d) / y0 )

   END FUNCTION u_tau_f_ABL
!   
!------------------------------------------------------------------------------------------------------------------------
!

   REAL(KIND=RP) FUNCTION tau_w_IBM( u_tau, y, rho, nu )
      USE WallFunctionDefinitions, ONLY: kappa, WallC
      IMPLICIT NONE 

      real(kind=RP), intent(in) :: u_tau, y, rho, nu 

      real(kind=RP) :: y_plus, gradU  
      ! 
      ! y_plus = y * u_tau/nu
      ! \partial y_plus/ \partial y = u_tau/nu 
      ! _______________________________________________________________________________________________
      ! u = u_tau * ( 1/kappa * log( 1 + kappa * y_plus ) ) &
      !   + u_tau * ( WallC - 1 / kappa * log(kappa) ) * & 
      !   ( 1 - exp( -y_plus / 11 )-( y_plus / 11 ) * exp( -y_plus / 3 ) )
      ! _______________________________________________________________________________________________
      ! tau(y) = mu * (\partial u/\partial y)|_y
      ! _______________________________________________________________________________________________
      ! (\partial u/\partial y)|_y = u_tau * (u_tau/nu)/( 1 + kappa * y_plus )    &
      !                            + u_tau * ( WallC - 1 / kappa * log(kappa) ) * &
      !                             ( exp( -y_plus / 11 )*(u_tau/nu)/11  - (u_tau/nu)/11 * exp( -y_plus / 3 ) + y_plus/11 * exp(-y_plus/3) * (u_tau/nu)/3 )


      y_plus = y_plus_f (y, u_tau, nu)
      gradU  = u_tau*u_tau/( 1.0_RP + kappa * y_plus ) + ( WallC - 1.0_RP / kappa * log(kappa) ) * u_tau * u_tau /(11.0_RP) * &
               ( exp(-y_plus/11.0_RP) - exp(-y_plus/3.0_RP) + y_plus/3.0_RP * exp(-y_plus/3.0_RP) )

      tau_w_IBM = rho * gradU 

   END FUNCTION tau_w_IBM

   SUBROUTINE WallFunctionBC_FlowNeumann_HOIBM( N, Q, dWall, nHat, xsb, nodes, u_tau, visc_fluxsb )   
      use PolynomialInterpAndDerivsModule
      use VariableConversion
      IMPLICIT NONE 

      real(kind=RP), intent(in)    :: Q(NCONS,0:N)
      integer,       intent(in)    :: N 
      real(kind=RP), intent(in)    :: dWall(0:N), xsb, nodes(0:N)
      real(kind=RP), intent(in)    :: nHat(NDIM)
      real(kind=RP), intent(inout) :: u_tau
      real(kind=RP), intent(inout) :: visc_fluxsb(NCONS)

      real(kind=RP) :: U_ref(NDIM), u_parallel(NDIM), x_II(NDIM), u_II, Q_ref(NCONS)
      real(kind=RP) :: lj(0:N), w(0:N), den, visc_flux(NDIM,0:N)
      real(kind=RP) :: mu_ref, kappa_ref, nu_ref, mu, nu
      integer       :: i

      visc_fluxsb = 0.0_RP 
      den         = 0.0_RP
      Q_ref       = Q(:,N)

      U_ref      = Q_ref(IRHOU:IRHOW)/Q_ref(IRHO)
      u_parallel = U_ref - (dot_product(U_ref, nHat) * nHat)
      x_II       = u_parallel / norm2(u_parallel)
      u_II       = dot_product(U_ref, x_II)

      call get_laminar_mu_kappa(Q_ref, mu_ref, kappa_ref)

      nu_ref = mu_ref/Q_ref(IRHO)
      u_tau  = u_tau_f (u_II, dWall(N), nu_ref, u_tau)

      do i = 0, N 
         call get_laminar_mu_kappa(Q(:,i), mu, kappa_ref)
         nu = mu/Q(IRHO,i)
         visc_flux(:,i) = tau_w_IBM( u_tau, dWall(i), Q(IRHO,i), nu ) * x_II 
      end do 

      do i = 0, N
         lj(i) = LagrangeInterpolatingPolynomial( i, xsb, N, nodes )
         w(i)  = 1.0_RP/(abs(-1.0_RP - nodes(i))**2 + 1.0e-10)
         den   = den + w(i) * lj(i)
      end do

      do i = 0, N
         visc_fluxsb(IRHOU) = visc_fluxsb(IRHOU) + visc_flux(IX,i) * lj(i) * w(i)
         visc_fluxsb(IRHOV) = visc_fluxsb(IRHOV) + visc_flux(IY,i) * lj(i) * w(i)
         visc_fluxsb(IRHOW) = visc_fluxsb(IRHOW) + visc_flux(IZ,i) * lj(i) * w(i)
      end do 

      visc_fluxsb = -visc_fluxsb/den 

   END SUBROUTINE WallFunctionBC_FlowNeumann_HOIBM

END MODULE 
#endif