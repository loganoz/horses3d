#include "Includes.h"

module WallDistance
#if defined(NAVIERSTOKES)
   use SMConstants
   use NodalStorageClass
   use Utilities, only: almostEqual

   private
   public   ComputeWallDistance
!
!  ************
!  Return codes
!  ************
!
   integer, parameter       :: SUCCESS_RC = 0
   integer, parameter       :: EXCEED_ITERS_RC = -1
!
!  ***********************
!  Quasi-Newton parameters
!  ***********************
!
   integer,       parameter :: N_MAX_ITER_NEWTON  =  50
   real(kind=RP), parameter :: MIN_GRAD_TOL       =  1.0e-11_RP
   real(kind=RP), parameter :: ALPHA_LIMITER      =  1.0_RP
!
!  *****************************************
!  Global quantities to store previous steps
!  *****************************************
!
   real(kind=RP)     :: x0(2)
   real(kind=RP)     :: phi0
   real(kind=RP)     :: dphi0(2)
   real(kind=RP)     :: B0(2,2)
!
!  *************************
!  Line search HZ parameters
!  *************************
!
   integer,       parameter  :: MAX_ITERS_LS_HZ           = 100
   integer,       parameter  :: MAX_ITERS_LS_HZ_bracket   = 50
   integer,       parameter  :: MAX_ITERS_LS_HZ_update    = 100
   real(kind=RP), parameter  :: GAMMA_LS_HZ               = 0.66_RP
   real(kind=RP), parameter  :: PSI0_LS_HZ                = 0.01_RP
   real(kind=RP), parameter  :: PSI1_LS_HZ                = 0.1_RP
   real(kind=RP), parameter  :: PSI2_LS_HZ                = 2._RP
   real(kind=RP), parameter  :: TOLA_LS_HZ_init           = 1e-5_RP
   real(kind=RP), parameter  :: RHO_LS_HZ                 = 5._RP
   real(kind=RP), parameter  :: EPSILON_LS_HZ             = 1e-6_RP
   real(kind=RP), parameter  :: THETA_LS_HZ               = 0.5_RP
   real(kind=RP), parameter  :: ONEMTHETA_LS_HZ           = 1._RP-THETA_LS_HZ
   real(kind=RP), parameter  :: DELTA_LS_HZ               = 0.1_RP
   real(kind=RP), parameter  :: DELTAAP_LS_HZ             = 2._RP*DELTA_LS_HZ -1._RP
   real(kind=RP), parameter  :: SIGMA_LS_HZ               = 0.9_RP
!
!  ****************************
!  Hessian Estimation functions 
!  ****************************
!
   abstract interface
      subroutine HE_F(x, dphi, B)
         use SMConstants
         implicit none
         real(kind=RP), intent(in)  :: x(2)
         real(kind=RP), intent(in)  :: dphi(2)
         real(kind=RP), intent(out) :: B(2,2)
      end subroutine HE_F
   end interface

   contains
      subroutine ComputeWallDistance(xP, Nf, xf, spAf, xSeed, d, xi_Wall, RC) 
!
!        *****************************************************************
!
!              This method performs an optimization of the distance to
!           the wall given by its nodal coordinates "xf".
!           The algorithm consists in a Quasi-Newton optimization core
!           with a gradient projection method to consider the domain
!           end points.
!
!        *****************************************************************
!
         implicit none
         real(kind=RP),       intent(in)  :: xP(NDIM)
         integer,             intent(in)  :: Nf(2)
         real(kind=RP),       intent(in)  :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t), intent(in)  :: spAf(2)
         real(kind=RP),       intent(in)  :: xSeed(2)
         real(kind=RP),       intent(out) :: d
         real(kind=RP),       intent(out) :: xi_Wall(2)
         integer,             intent(out) :: RC
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iter
         real(kind=RP)  :: x(2), phi, dphi(2), dir_corr(2)
!
!        Initialize fitness function and gradient
!        ----------------------------------------
         x0 = xSeed
         call FitnessFunctionAndGradient(x0, xP, Nf, xf, spAf, phi0, dphi0)
!
!        Initialize the inverse hessian matrix
!        -------------------------------------         
         B0 = 0.0_RP
         B0(1,1) = 1.0_RP
         B0(2,2) = 1.0_RP

         x = x0

         do iter = 0, N_MAX_ITER_NEWTON
!
!           Get new fitness function values
!           -------------------------------
            call FitnessFunctionAndGradient(x, xP, Nf, xf, spAf, phi, dphi)
!
!           Check the exit criteria
!           -----------------------
            if ( exitCriteria(x, dphi) ) then
               xi_Wall = x
               d       = sqrt(phi)
               RC      = SUCCESS_RC
               return
            
            end if

            call OptimizationStep(xP, Nf, xf, spAf, phi, dphi, x)

         end do

         xi_Wall = x
         d       = sqrt(phi)
         RC      = EXCEED_ITERS_RC

      end subroutine ComputeWallDistance

      logical function exitCriteria(x, dphi)
!
!        ****************************************
!        The only exit criteria considered is the
!        size of the projected gradient
!        ****************************************
!
         implicit none
         real(kind=RP), intent(in)  :: x(2)
         real(kind=RP), intent(in)  :: dphi(2)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: dir(2)

         dir = -dphi

         call CorrectDirection(x,dir)
      
         if ( norm2(dir) .le. MIN_GRAD_TOL ) then
            exitCriteria = .true.

         else
            exitCriteria = .false.

         end if

      end function exitCriteria

      subroutine OptimizationStep(xP, Nf, xf, spAf, phi, dphi, x)
!
!        *************************************************************
!
!              An optimization step consists in:
!           1/ Estimate the Hessian inverse
!           2/ Compute the Newton step direction
!           3/ Correct the direction according to the domain limits
!           4/ Compute the maximum Quasi-Newton alpha coefficient
!           5/ Perform a Line Search along the direction
!
!        *************************************************************
!
         implicit none
         real(kind=RP), intent(in)         :: xP(NDIM)
         integer,       intent(in)         :: Nf(2)
         real(kind=RP), intent(in)         :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t), intent(in) :: spAf(2)
         real(kind=RP), intent(in)         :: phi
         real(kind=RP), intent(in)         :: dphi(2)
         real(kind=RP), intent(inout)      :: x(2)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: B(2,2), dir(2), alphaMax, alpha
         integer        :: RC
         real(kind=RP)  :: y(2), f_change
         real(kind=RP)  :: phiOut, dPhiOut(2)
!
!        Update Hessian matrix inverse
!        -----------------------------
         call HessianEstimation_DFP(x,dphi,B)
!
!        Update values for the next iteration
!        ------------------------------------
         x0    = x
         phi0  = phi
         dphi0 = dphi
         B0    = B
!
!        Compute the Newton method direction 
!        -----------------------------------
         dir = -matmul(B,dphi)
!
!        Correct the direction to account for domain limits
!        --------------------------------------------------
         call correctDirection(x, dir)
!
!        Get an estimation of the maximum alpha coefficient
!        --------------------------------------------------
         call getAlphaMax(x, dir, alphaMax)
!
!        Perform a line search along the direction
!        -----------------------------------------
         phiOut  = phi
         dphiOut = dphi
         alpha   = alphaMax
         RC = LS_HZ(xP, Nf, xf, spAf, alpha, x, dir, phiOut, dphiOut, y, f_change)
!
!        Get a new alpha estimation with the LineSearch resulting direction
!        ------------------------------------------------------------------
         call getAlphaMax(x0, dir, alphaMax)
!
!        Perform a limiting if alpha has been exceeded
!        ---------------------------------------------
         if ( alpha .gt. alphaMax ) then
            x = x0 + alphaMax * dir
         end if

      end subroutine OptimizationStep

      subroutine HessianEstimation_DFP_DFP(x, dphi, B)
!
!        ******************************************************
!              This (inverse) Hessian estimation is performed
!           using the Davidon-Fletcher-Powell formula
!        ******************************************************
!
         implicit none
         real(kind=RP), intent(in)  :: x(2)
         real(kind=RP), intent(in)  :: dphi(2)
         real(kind=RP), intent(out) :: B(2,2)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j
         real(kind=RP) :: u(2), g(2)
         real(kind=RP) :: uu(2,2), gg(2,2)

         u = x-x0
         g = dphi - dphi0

         do j = 1, 2 ; do i = 1, 2
           uu(i,j) = u(i) * u(j)
           gg(i,j) = g(i) * g(j) 
         end do      ; end do

         if ( almostEqual(dot_product(u,g), 0.0_RP) ) then
            B = B0
         else 
            B = B0 + uu / dot_product(u,g) - matmul(B0,matmul(gg,B0))/ (dot_product(g,matmul(B0,g)))
         end if

      end subroutine HessianEstimation_DFP_DFP

      subroutine correctDirection(x, dir)
!
!        ****************************************************************
!
!              Correct the search direction so that it remains
!           parallel to the domain bounds if the point lays in those
!
!        ****************************************************************
!
         implicit none
         real(kind=RP), intent(in)    :: x(2)
         real(kind=RP), intent(inout) :: dir(2)
      
         if ( almostEqual(x(1), -1.0_RP) .and. (dir(1) .lt. 0.0_RP) ) then
            dir(1) = 0.0_RP
   
         elseif ( almostEqual(x(1), 1.0_RP) .and. (dir(1) .gt. 0.0_RP) ) then
            dir(1) = 0.0_RP

         end if

         if ( almostEqual(x(2), -1.0_RP) .and. (dir(2) .lt. 0.0_RP) ) then
            dir(2) = 0.0_RP
   
         elseif ( almostEqual(x(2), 1.0_RP) .and. (dir(2) .gt. 0.0_RP) ) then
            dir(2) = 0.0_RP

         end if

      end subroutine correctDirection

      subroutine getAlphaMax(x, dir, alphaMax)
!
!        *******************************************************
!
!              The maximum alpha is estimated by computing the
!           minimum distance to the domain limits following
!           the search direction
!
!        *******************************************************
!
         implicit none
         real(kind=RP), intent(in)  :: x(2)
         real(kind=RP), intent(in)  :: dir(2)
         real(kind=RP), intent(out) :: alphaMax
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: alpha(4)

         if ( (.not. almostEqual(dir(1), 0.0_RP)) .and. (.not. almostEqual(x(1), -1.0_RP)) ) then
            alpha(1) = (-1.0_RP - x(1)) / dir(1)

         else
            alpha(1) = -1.0_RP

         end if

         if ( (.not. almostEqual(dir(1), 0.0_RP)) .and. (.not. almostEqual(x(1), 1.0_RP)) ) then
            alpha(2) = (1.0_RP - x(1)) / dir(1)

         else
            alpha(2) = -1.0_RP

         end if

         if ( (.not. almostEqual(dir(2), 0.0_RP)) .and. (.not. almostEqual(x(2), -1.0_RP)) ) then
            alpha(3) = (-1.0_RP - x(2)) / dir(2)

         else
            alpha(3) = -1.0_RP

         end if

         if ( (.not. almostEqual(dir(2), 0.0_RP)) .and. (.not. almostEqual(x(2), 1.0_RP)) ) then
            alpha(4) = (1.0_RP - x(2)) / dir(2)

         else
            alpha(4) = -1.0_RP

         end if

         alphaMax = minval(alpha,(alpha .gt. 0.0_RP))

         alphaMax = min(abs(alphaMax), ALPHA_LIMITER)
   
      end subroutine getAlphaMax

      subroutine FitnessFunction(x, xP, Nf, xf, spAf, phi)
!
!        ************************************************************
!              This subroutine evaluates the fitness function in
!           the point x. The fitness function is the squared
!           distance to the wall
!        ************************************************************
!
         implicit none
         real(kind=RP),       intent(in)  :: x(2)
         real(kind=RP),       intent(in)  :: xP(NDIM)
         integer,             intent(in)  :: Nf(2)
         real(kind=RP),       intent(in)  :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t), intent(in)  :: spAf(2)
         real(kind=RP),       intent(out) :: phi
!
!        ---------------
!        Local variables
!        ---------------
!
         integer           :: i, j
         real(kind=RP)     :: xWall(NDIM)
         real(kind=RP)     :: lxi(0:Nf(1)), leta(0:Nf(2))
!
!        Get the interpolating polynomials
!        ---------------------------------
         lxi  = spAf(1) % lj(x(1))
         leta = spAf(2) % lj(x(2))
!
!        Get the surface point
!        ---------------------
         xWall = 0.0_RP
         do j = 0, Nf(2)   ; do i = 0, Nf(1)
            xWall = xWall + xf(:,i,j) * lxi(i) * leta(j)
         end do            ; end do
!
!        Compute the fitness function
!        ----------------------------
         phi = sum( POW2(xWall-xP) )

      end subroutine FitnessFunction

      subroutine FitnessFunctionAndGradient(x, xP, Nf, xf, spAf, phi, dphi)
!
!        ************************************************************
!              This subroutine evaluates the fitness function and
!           its gradient in the point x.
!        ************************************************************
!
         implicit none
         real(kind=RP),       intent(in)  :: x(2)
         real(kind=RP),       intent(in)  :: xP(NDIM)
         integer,             intent(in)  :: Nf(2)
         real(kind=RP),       intent(in)  :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t), intent(in)  :: spAf(2)
         real(kind=RP),       intent(out) :: phi
         real(kind=RP),       intent(out) :: dphi(2)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer           :: i, j
         real(kind=RP)     :: xWall(NDIM), xWall_xi(NDIM), xWall_eta(NDIM)
         real(kind=RP)     :: lxi(0:Nf(1))  , leta(0:Nf(2))
         real(kind=RP)     :: dlxi(0:Nf(1)) , dleta(0:Nf(2))
!
!        Get the interpolating polynomials and derivatives
!        -------------------------------------------------
         lxi  = spAf(1) % lj(x(1))
         leta = spAf(2) % lj(x(2))

         dlxi  = spAf(1) % dlj(x(1))
         dleta = spAf(2) % dlj(x(2))
!
!        Get the surface point and derivatives
!        -------------------------------------
         xWall     = 0.0_RP
         xWall_xi  = 0.0_RP
         xWall_eta = 0.0_RP
         do j = 0, Nf(2)   ; do i = 0, Nf(1)
            xWall     = xWall     + xf(:,i,j) * lxi(i)  * leta(j)
            xWall_xi  = xWall_xi  + xf(:,i,j) * dlxi(i) * leta(j)
            xWall_eta = xWall_eta + xf(:,i,j) * lxi(i)  * dleta(j)
         end do            ; end do
!
!        Compute the fitness function
!        ----------------------------
         phi = sum( POW2(xWall-xP) )
!
!        Compute the fitness function gradient
!        -------------------------------------
         dphi(1) = 2.0_RP * sum((xWall - xP)* xWall_xi )
         dphi(2) = 2.0_RP * sum((xWall - xP)* xWall_eta)

      end subroutine FitnessFunctionAndGradient
!
!///////////////////////////////////////////////////////////////////////////////
!
!        Line Search algorithms
!        ----------------------
!
!///////////////////////////////////////////////////////////////////////////////
!
      function LS_HZ(xP, Nf, xf, spAf, alpha, x, dir, fval, g, y, f_change) result(exit_flag)
!      
!        ***************************************************************************
!           This is an implementation of the Hager-Zhang line search algorithm
!        ***************************************************************************
         implicit none
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         real(kind=RP),               intent(inout) :: alpha
         real(kind=RP), dimension(2), intent(inout) :: x
         real(kind=RP), dimension(2), intent(inout) :: dir
         real(kind=RP),               intent(inout) :: fval
         real(kind=RP), dimension(2), intent(inout) :: g
         real(kind=RP), dimension(2), intent(out)   :: y
         real(kind=RP),               intent(out)   :: f_change
         integer                                    :: exit_flag
!      
!        ---------------
!        Local variables
!        ---------------
!      
         real (kind=RP)                :: a, b, w, wohi, wolo, awohi, fpert
         integer                       :: i
         real (kind=RP)                :: phi, dphi
         real (kind=RP)                :: phia, dphia
         real (kind=RP), dimension (2) :: g1, ga
         real (kind=RP), dimension (2) :: d1, da
         real (kind=RP), dimension (2) :: x1, xa
!      
!        Compute Wolfe condition bounds
!        ------------------------------
         awohi = dot_product(g, dir)
         wohi  = DELTA_LS_HZ * awohi      ! Wolfe upper bound
         wolo  = SIGMA_LS_HZ * awohi      ! Wolfe lower bound
         awohi = DELTAAP_LS_HZ * awohi    ! Approximated Wolfe upper bound
         fpert = fval + EPSILON_LS_HZ     ! Perturbation for the AP Wolfe: phi(alpha) = phi(0) + eps_k
         
         if (LS_HZ_init(alpha, x, dir, fval, g,  wohi, wolo, awohi, fpert, y, f_change, xP, Nf, xf, spAf)) then
             exit_flag = 0
             return
         elseif (LS_HZ_bracket(a, b, alpha, x, dir, fval, g, y , f_change , fpert, wohi, wolo, awohi, exit_flag, xP, Nf, xf, spAf).or. &
             (exit_flag.ne.0)) then
         return
         end if
         
         w=b-a
         
         da = a*dir
         xa = x+ da
         call FitnessFunctionAndGradient(xa, xP, Nf, xf, spAf, phia, ga)
         dphia = dot_product(ga, dir)
         
         d1 = b*dir
         x1 = x+ d1
         call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
         dphi = dot_product(g1, dir)
         
         do i=1, MAX_ITERS_LS_HZ
         
             if ( LS_HZ_secant2(a,b,alpha, x, dir, g, y, fval, f_change, phi, dphi, phia, dphia, g1, d1, x1, &
                                ga, da, xa, wohi, wolo, awohi,fpert, exit_flag, xP, Nf, xf, spAf).or.&
                 (exit_flag.ne.0)) return
         
             if((b-a)>GAMMA_LS_HZ*w) then
                 alpha=(a+b)*0.5_RP
                 d1 = alpha* dir
                 x1 = x+ d1
                 call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
                 dphi = dot_product(g1, dir)
         
                 if (LS_HZ_check(alpha, fval, phi, dphi, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
                     exit_flag = 0
                     f_change  = phi-fval
                     fval      = phi
                     y         = g1-g
                     g         = g1
                     dir       = d1
                     x         = x + alpha * dir
                     return
                 end if
                 exit_flag= LS_HZ_update(a,b,alpha,x, dir, phi, dphi, fpert, g1, x1, d1,ga, xa, da, phia, dphia, xP, Nf, xf, spAf)
                 if(exit_flag.ne.0) return
             end if
             w=b-a
         end do
         
         exit_flag=-2
         
      end function
      
      function LSinit_HZ(alpha, x, dir, fval, g, y, f_change, xnorm, gnorm, xP, Nf, xf, spAf) result(exit_flag)
         implicit none
         real(kind=RP),               intent(inout) :: alpha
         real(kind=RP), dimension(2), intent(inout) :: x
         real(kind=RP), dimension(2), intent(inout) :: dir
         real(kind=RP),               intent(inout) :: fval
         real(kind=RP), dimension(2), intent(inout) :: g
         real(kind=RP), dimension(2), intent(out)   :: y
         real(kind=RP),               intent(out)   :: f_change
         real(kind=RP),               intent(in)    :: xnorm
         real(kind=RP),               intent(in)    :: gnorm
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         integer                                    :: exit_flag
!      
!        ---------------
!        Local variables   
!        ---------------
!      
         real (kind=RP)                :: a, b, w,wohi, wolo, awohi, fpert
         integer                       :: i
         real (kind=RP)                :: phi, dphi
         real (kind=RP)                :: phia, dphia
         real (kind=RP), dimension (2) :: g1, ga
         real (kind=RP), dimension (2) :: d1, da
         real (kind=RP), dimension (2) :: x1, xa
         
         
         awohi = dot_product(g, dir)
         wohi  = DELTA_LS_HZ * awohi
         wolo  = SIGMA_LS_HZ * awohi
         awohi = DELTAAP_LS_HZ * awohi
         fpert = fval + EPSILON_LS_HZ
         
         if (LS_HZ_init1(alpha, x, dir, fval, g,  wohi, wolo, awohi, fpert, y, f_change, xnorm, gnorm, xP, Nf, xf, spAf)) then
             exit_flag =0
             return
         elseif (LS_HZ_bracket(a, b, alpha, x, dir, fval, g, y , f_change , fpert, wohi, wolo, awohi, exit_flag, xP, Nf, xf, spAf).or. &
             (exit_flag.ne.0))  then
         return
         end if
         
         w=b-a
         
         da = a*dir
         xa = x+ da
!         call Compute_func_gradient (phia, ga, xa)
         call FitnessFunctionAndGradient(xa, xP, Nf, xf, spAf, phia, ga)
         dphia = dot_product(ga, dir)
         
         d1 = b*dir
         x1 = x+ d1
!         call Compute_func_gradient (phi, g1, x1)
         call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
         dphi = dot_product(g1, dir)
         
         do i=1, MAX_ITERS_LS_HZ
         
             if ( LS_HZ_secant2(a,b,alpha, x, dir, g, y, fval, f_change, phi, dphi, phia, dphia, g1, d1, x1, ga, &
                               da, xa, wohi, wolo, awohi,fpert, exit_flag, xP, Nf, xf, spAf).or.&
                 (exit_flag.ne.0)) return
         
             if((b-a)>GAMMA_LS_HZ*w) then
                 alpha=(a+b)*0.5_RP
                 d1 = alpha* dir
                 x1 = x+ d1
!                 call compute_func_gradient (phi, g1, x1)
                 call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
                 dphi = dot_product(g1, dir)
         
                 if (LS_HZ_check(alpha, fval, phi, dphi, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
                     exit_flag = 0
                     f_change  = phi-fval
                     fval      = phi
                     y         = g1-g
                     g         = g1
                     dir       = d1
                     x         = x1
                     return
                 end if
                 exit_flag= LS_HZ_update(a,b,alpha,x, dir, phi, dphi, fpert, g1, x1, d1,ga, xa, da, phia, dphia, xP, Nf, xf, spAf)
                 if(exit_flag.ne.0) return
             end if
             w=b-a
         end do
         
         exit_flag=-2
         
      end function
      
      function LS_HZ_secant2(a,b,alpha, x, dir, g, y, fval, f_change, phi, dphi, phia, dphia, &
                             g1, d1, x1, ga, da, xa, wohi, wolo, awohi,fpert, flag, xP, Nf, xf, spAf) result (found)
         implicit none
         real(kind=RP),               intent(inout) :: a
         real(kind=RP),               intent(inout) :: b
         real(kind=RP),               intent(inout) :: alpha
         real(kind=RP), dimension(2), intent(inout) :: x
         real(kind=RP), dimension(2), intent(inout) :: dir
         real(kind=RP), dimension(2), intent(inout) :: g
         real(kind=RP), dimension(2), intent(inout) :: y
         real(kind=RP),               intent(inout) :: fval
         real(kind=RP),               intent(inout) :: f_change
         real(kind=RP),               intent(inout) :: phi, dphi
         real(kind=RP),               intent(inout) :: phia, dphia
         real(kind=RP), dimension(2), intent(inout) :: g1, ga
         real(kind=RP), dimension(2), intent(inout) :: d1, da
         real(kind=RP), dimension(2), intent(inout) :: x1, xa
         real(kind=RP),               intent(in)    :: wohi
         real(kind=RP),               intent(in)    :: wolo
         real(kind=RP),               intent(in)    :: awohi
         real(kind=RP),               intent(in)    :: fpert
         integer,                     intent(out)   :: flag
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         logical                                    :: found
!      
!        ---------------
!        Local variables   
!        ---------------
!      
         real (kind=RP)   :: c, aa, bb, dphiaa, dphibb
         
         c = (a*dphi-b*dphia)/(dphi-dphia)
         
         aa=a
         bb=b
         flag=LS_HZ_update(aa,bb,c,x, dir, phi, dphibb, fpert, g1, x1, d1, &
             ga, xa, da, phia, dphiaa, xP, Nf, xf, spAf)
         
         if(flag.ne.0) then
             found = .false.
             return
         elseif (c.eq.bb) then
             c = (b*dphibb-bb*dphi)/(dphibb-dphi)
             goto 100
         elseif (c.eq.aa) then
             c = (a*dphiaa-aa*dphia)/(dphiaa-dphia)
             goto 100
         else
             dphia=dphiaa
             dphi=dphibb
             goto 200
         end if
         
         100 dphia=dphiaa;dphi=dphibb
         flag=LS_HZ_update(aa,bb,c,x, dir, phi, dphi, fpert, g1, x1, d1, &
             ga, xa, da, phia, dphia, xP, Nf, xf, spAf)
         
         
         if(flag.ne.0) then
             found = .false.
             return
         end if
         
         
         200 a= aa; b=bb
         flag = 0
         if (LS_HZ_check(a, fval, phia, dphia, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
             alpha     = a
             f_change  = phia-fval
             fval      = phia
             y         = ga-g
             g         = ga
             dir       = da
             x         = xa
             found     = .true.
         elseif (LS_HZ_check(b, fval, phi, dphi, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
             alpha     = b
             f_change  = phi-fval
             fval      = phi
             y         = g1-g
             g         = g1
             dir       = d1
             x         = x1
             found     = .true.
         else
             found     = .false.
         end if
         return
         
         
      end function
      
      function LS_HZ_update(a,b,alpha,x, dir, phi, dphi, fpert, g1, x1, d1, &
          ga, xa, da, phia, dphia, xP, Nf, xf, spAf) result(flag)
         implicit none
         real(kind=RP),               intent(inout) :: a
         real(kind=RP),               intent(inout) :: b
         real(kind=RP),               intent(in)    :: alpha
         real(kind=RP), dimension(2), intent(in)    :: x
         real(kind=RP), dimension(2), intent(in)    :: dir
         real(kind=RP),               intent(inout) :: phi
         real(kind=RP),               intent(inout) :: dphi
         real(kind=RP),               intent(out)   :: phia
         real(kind=RP),               intent(out)   :: dphia
         real(kind=RP),               intent(in)    :: fpert
         real(kind=RP), dimension(2), intent(out)   :: g1, ga
         real(kind=RP), dimension(2), intent(out)   :: x1, xa
         real(kind=RP), dimension(2), intent(out)   :: d1, da
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         integer                                    :: flag
!      
!        ---------------
!        Local variables   
!        ---------------
!      
         real(kind=RP)               :: dphi_alpha, phi_alpha
         real(kind=RP), dimension(2) :: d_alpha, x_alpha, g_alpha
         integer                     :: i
         
         if ((alpha<a).or.(alpha>b)) then
             flag = 0
             return
         end if
         
         d_alpha = alpha*dir
         x_alpha = x+ d_alpha
!         call Compute_func_gradient (phi_alpha, g_alpha, x_alpha)
         call FitnessFunctionAndGradient(x_alpha, xP, Nf, xf, spAf, phi_alpha, g_alpha)
         dphi_alpha = dot_product(g_alpha, dir)
         
         
         if (dphi_alpha>=0._RP) then
             b=alpha
             phi=phi_alpha
             dphi=dphi_alpha
             g1=g_alpha
             d1=d_alpha
             x1=x_alpha
             flag = 0
             return
         elseif (phi_alpha<=fpert) then
             a= alpha
             phia=phi_alpha
             dphia=dphi_alpha
             ga=g_alpha
             da=d_alpha
             xa=x_alpha
             flag = 0
             return
         end if
         
         b = alpha
         
         do i = 1, MAX_ITERS_LS_HZ_update
             phi_alpha = ONEMTHETA_LS_HZ*a+THETA_LS_HZ*b
             d1 = phi_alpha * dir
             x1 = x + d1
!             call Compute_func_gradient ( phi, g1, x1)
             call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
             dphi = dot_product(g1, dir)
             if (dphi>=0._RP) then
                 b = phi_alpha
                 flag = 0
                 return
             elseif (phi<=fpert) then
                 a = phi_alpha
                 ga = g1
                 xa = x1
                 da = d1
                 phia = phi
                 dphia = dphi
                 cycle
             end if
             b = phi_alpha
         end do
         flag = -5
      end function
      
!DEC$ ATTRIBUTES FORCEINLINE :: LS_HZ_check
      function LS_HZ_check(alpha, fval, phi, dphi, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf) result(ok)
         implicit none
         real(kind=RP),         intent(inout) :: alpha
         real(kind=RP),         intent(in)    :: fval
         real(kind=RP),         intent(in)    :: phi
         real(kind=RP),         intent(in)    :: dphi
         real(kind=RP),         intent(in)    :: wohi
         real(kind=RP),         intent(in)    :: wolo
         real(kind=RP),         intent(in)    :: awohi
         real(kind=RP),         intent(in)    :: fpert
         real(kind=RP),         intent(in)    :: xP(3)
         integer,               intent(in)    :: Nf(2)
         real(kind=RP),         intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t), intent(in)    :: spAf(2)
         logical                              :: ok
         
         ok = ((dphi>=wolo).and.(((phi-fval)<=(alpha*wohi)).or.((phi<=fpert).and.(dphi<=awohi))))
         
      end function
      
      function LS_HZ_init(alpha, x, dir, f, g,  wohi, wolo, awohi, fpert, y, f_change, xP, Nf, xf, spAf) result(found)
         implicit none
         
         real(kind=RP),               intent(inout) :: alpha
         real(kind=RP), dimension(2), intent(inout) :: x
         real(kind=RP), dimension(2), intent(inout) :: dir
         real(kind=RP),               intent(inout) :: f
         real(kind=RP), dimension(2), intent(inout) :: g
         real(kind=RP),               intent(in)    :: wohi
         real(kind=RP),               intent(in)    :: wolo
         real(kind=RP),               intent(in)    :: awohi
         real(kind=RP),               intent(in)    :: fpert
         real(kind=RP), dimension(2), intent(out)   :: y
         real(kind=RP),               intent(out)   :: f_change
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         logical                                    :: found
!      
!        ---------------
!        Local variables   
!        ---------------
!      
         real (kind=RP), dimension(2) :: g1, x1, d1
         real (kind=RP)               :: alpha_q, a, phi, dphi
      
         alpha_q=PSI1_LS_HZ*alpha
         d1 = alpha_q*dir
         x1 = x + d1
!         call Compute_func_gradient(phi, g1, x1)
         call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
         if (phi.le.f) then
             dphi=alpha_q*alpha_q
             a=dot_product(g,dir)
         
             alpha_q=0.5_RP*dphi*a/(f-phi+alpha_q*a)
         
             a=(-f+phi*(1._RP-a))/dphi
             if ((alpha_q.ge.0._RP).and.(a>TOLA_LS_HZ_init)) then
                 alpha=alpha_q
                 d1 = alpha*dir
                 x1 = x + d1
!                 call Compute_func_gradient(phi, g1, x1)
                 call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
                 dphi = dot_product(g1, dir)
                 goto 100
             end if
         end if
         
         
         alpha = PSI2_LS_HZ*alpha
         d1 = alpha*dir
         x1 = x + d1
!         call Compute_func_gradient(phi, g1, x1)
         call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
         dphi = dot_product(g1,dir)
         
         
         
         100 if (LS_HZ_check(alpha, f, phi, dphi, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
             f_change  = phi-f
             f         = phi
             y         = g1-g
             g         = g1
             dir       = d1
             x         = x1
             found = .true.
         else
             found = .false.
         end if
      
      end function
      
      function LS_HZ_init1(alpha, x, dir, f, g,  wohi, wolo, awohi, fpert, y, f_change, xnorm, gnorm, xP, Nf, xf, spAf) result(found)
         implicit none
         
         real(kind=RP),               intent(inout) :: alpha
         real(kind=RP), dimension(2), intent(inout) :: x
         real(kind=RP), dimension(2), intent(inout) :: dir
         real(kind=RP),               intent(inout) :: f
         real(kind=RP), dimension(2), intent(inout) :: g
         real(kind=RP),               intent(in)    :: wohi
         real(kind=RP),               intent(in)    :: wolo
         real(kind=RP),               intent(in)    :: awohi
         real(kind=RP),               intent(in)    :: fpert
         real(kind=RP), dimension(2), intent(out)   :: y
         real(kind=RP),               intent(out)   :: f_change
         real(kind=RP),               intent(in)    :: xnorm
         real(kind=RP),               intent(in)    :: gnorm
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         logical                                    :: found
!      
!        ---------------
!        Local variables   
!        ---------------
!      
         real(kind=RP), dimension(2) :: g1, x1, d1
         real(kind=RP)               :: phi, dphi
         
         if (xnorm>0._RP) then
             alpha=PSI0_LS_HZ*xnorm/gnorm
         elseif (f.ne.0._RP) then
             alpha=PSI0_LS_HZ*abs(f)/dot_product(g, g)
         else
             alpha = 1._RP
         end if
         
         d1 = alpha*dir
         x1 = x + d1
!         call Compute_func_gradient(phi, g1, x1)
         call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi, g1)
         dphi = dot_product(g1,dir)
         if (LS_HZ_check(alpha, f, phi, dphi, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
             f_change  = phi-f
             f         = phi
             y         = g1-g
             g         = g1
             dir       = d1
             x         = x1
             found = .true.
         else
             found = .false.
         end if
         
      end function
      
      function LS_HZ_bracket (a, b, alpha, x, dir, fval, g, y , f_change , &
          fpert, wohi, wolo, awohi, flag, xP, Nf, xf, spAf) result(found)
         implicit none
         real(kind=RP),               intent(out)   :: a
         real(kind=RP),               intent(out)   :: b
         real(kind=RP),               intent(inout) :: alpha
         real(kind=RP), dimension(2), intent(inout) :: x
         real(kind=RP), dimension(2), intent(inout) :: dir
         real(kind=RP),               intent(inout) :: fval
         real(kind=RP), dimension(2), intent(inout) :: g
         real(kind=RP), dimension(2), intent(out)   :: y
         real(kind=RP),               intent(out)   :: f_change
         real(kind=RP),               intent(in)    :: fpert
         real(kind=RP),               intent(in)    :: wohi
         real(kind=RP),               intent(in)    :: wolo
         real(kind=RP),               intent(in)    :: awohi
         integer,                     intent(out)   :: flag
         real(kind=RP),               intent(in)    :: xP(3)
         integer,                     intent(in)    :: Nf(2)
         real(kind=RP),               intent(in)    :: xf(NDIM,0:Nf(1),0:Nf(2))
         class(NodalStorage_t),       intent(in)    :: spAf(2)
         logical                                    :: found
!      
!        ---------------
!        Local variables   
!        ---------------
!      
         real(kind=RP)               :: cj, a_bar, b_bar, phi_der_cj, phi_cj, phia, dphia
         real(kind=RP), dimension(2) :: g1, ga, d1, x1, da, xa
         integer                     :: i, j
         
         cj=alpha
         a=0._RP
         do j=1,  MAX_ITERS_LS_HZ_bracket
             d1 = cj*dir
             x1 = x + d1
!             call compute_func_gradient(phi_cj,g1, x1)
             call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi_cj, g1)
             phi_der_cj=dot_product(g1, dir)
         
             if (phi_der_cj>=0._RP) then
                 b=cj
                 goto 100!return
             else if ((phi_der_cj<0._RP).and.(phi_cj>fpert)) then
                 a_bar=0._RP
                 b_bar=cj
                 do i=1,  MAX_ITERS_LS_HZ_bracket
                     cj=(THETA_LS_HZ)*b_bar+(ONEMTHETA_LS_HZ)*a_bar
                     d1 = cj*dir
                     x1 = x + d1
!                     call compute_func_gradient(phi_cj,g1, x1)
                     call FitnessFunctionAndGradient(x1, xP, Nf, xf, spAf, phi_cj, g1)
                     phi_der_cj=dot_product(g1, dir)
                     if(phi_der_cj>=0._RP) then
                         a=a_bar
                         b=cj
                         goto 100 !return
                     elseif(phi_cj<=fpert) then
                         a_bar=cj
                         phia = phi_cj
                         dphia = phi_der_cj
                         ga = g1
                         da = d1
                         xa = x1
                         cycle
                     end if
                     b_bar=cj
                 end do
                 goto 200 !err
             else
                 if (phi_cj<=fpert) then
                     a=cj
                     phia = phi_cj
                     dphia = phi_der_cj
                     ga = g1
                     da = d1
                     xa = x1
                 end if
                 cj=RHO_LS_HZ*cj
             end if
         end do
         goto 200 !err
         
         
         100 flag = 0
         if (LS_HZ_check(a, fval, phia, dphia, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
             alpha     = a
             f_change  = phia-fval
             fval      = phia
             y         = ga-g
             g         = ga
             dir       = da
             x         = xa
             found     = .true.
         elseif (LS_HZ_check(b, fval, phi_cj, phi_der_cj, wohi, wolo, awohi, fpert, xP, Nf, xf, spAf)) then
             alpha     = b
             f_change  = phi_cj-fval
             fval      = phi_cj
             y         = g1-g
             g         = g1
             dir       = d1
             x         = x1
             found     = .true.
         else
             found     = .false.
         end if
         return
         
         200 flag=-3
         
      end function
#endif
end module WallDistance
