!
!//////////////////////////////////////////////////////
!
!   @File:    ManufacturedSolutions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 13:23:11 2018
!   @Last revision date: Wed Apr 18 20:19:09 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 0d746cd20d04ebda97f349d7f3b0b0fe00b5d7ca
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      ManufacturedSolutions.f90
!      Created: 2017-05-18 17:00:00 +0100 
!      By: Andrés Rueda
!
!      Manufactured solutions definitions for 2D and 3D cases (Euler and Navier-Stokes)
!
!
!////////////////////////////////////////////////////////////////////////////////////////
MODULE ManufacturedSolutionsNSSA
   USE SMConstants
   USE PhysicsStorage_NSSA
   USE Physics_NSSA
   use FluidData_NSSA
   IMPLICIT NONE

   private
   public   InitializeManufacturedSol, ManufacturedSolP
   public   ManufacturedSolutionState
   public   ManufacturedSolutionSourceNSSA

   real(kind=RP) :: dwall, cv2, sbar, somega
   real(kind=RP) :: rhoinit, rhox, rhoy, rhoxy, arhox, arhoy, arhoxy 
   real(kind=RP) :: uinit, ux, uy, uxy, aux, auy, auxy 
   real(kind=RP) :: vinit, vx, vy, vxy, avx, avy, avxy 
   real(kind=RP) :: pinit, px, py, pxy, apx, apy, apxy 
   real(kind=RP) :: thetainit, thetax, thetay, thetaxy, athetax, athetay, athetaxy 
!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE InitializeManufacturedSol(ProblemDIM)
      IMPLICIT NONE
!
!     --------------------------------------------------------------------------------------
!     Constant definitions for the manufactured solutions. 
!
!        The primitive variable j(rho, u, v, w, p) is defined as:
!           j = jC(1) + jC(2)*Sin(pi*jC(5)*x) + jC(3)*Sin(pi*jC(6)*y) + jC(4)*Sin(pi*jC(7)*z) 
!     See:
!        |> Roy, Chris, Curt Ober, and Tom Smith. "Verification of a compressible CFD code 
!           using the method of manufactured solutions." 32nd AIAA Fluid Dynamics Conference 
!           and Exhibit. 2002.
!     --------------------------------------------------------------------------------------
!
      CHARACTER(LEN=*) :: ProblemDIM
            
            rhoinit = 1.0_RP
            rhox = 0.1_RP
            rhoy = -0.2_RP
            rhoxy = 0.1_RP
            arhox = 1.0_RP
            arhoy =1.0_RP
            arhoxy = 1.0_RP

            uinit = 2.0_RP
            ux = 0.3_RP
            uy = 0.3_RP
            uxy = 0.3_RP
            aux = 3.0_RP
            auy =1.0_RP
            auxy =1.0_RP

            vinit = 2.0_RP
            vx = 0.3_RP
            vy = 0.3_RP
            vxy = 0.3_RP
            avx = 1.0_RP
            avy =1.0_RP
            avxy = 1.0_RP

            pinit = 10.0_RP
            px = 1.0_RP
            py = 1.0_RP
            pxy = 0.5_RP
            apx = 2.0_RP
            apy = 1.0_RP
            apxy = 1.0_RP


            thetainit = 0.6_RP
            thetax = -0.03_RP
            thetay = -0.02_RP
            thetaxy = 0.02_RP
            athetax = 2.0_RP
            athetay = 1.0_RP
            athetaxy = 3.0_RP

      
   END SUBROUTINE InitializeManufacturedSol
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   FUNCTION ManufacturedSolP(x) RESULT(p)
      IMPLICIT NONE
      REAL(KIND=RP) :: x(3), p
            
   END FUNCTION ManufacturedSolP
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ManufacturedSolutionState( xvec, Q )
      IMPLICIT NONE
!
!     ----------------------------------------------------
!     Manufactured solution state (valid for Euler and NS)
!     ----------------------------------------------------
!      
      
      REAL(KIND=RP), intent(in)    :: xvec(3)
      REAL(KIND=RP), intent(inout) :: Q(NCONS)
      
      REAL(KIND=RP) :: p

      associate ( x => xvec(1), &
                  y => xvec(2), &
                  z => xvec(3), &
                  gamma => thermodynamics % gamma)


      Q(1) = rhoinit+rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y)+rhoy*Cos(arhoy*Pi*y)+rhox*Sin(arhox*Pi*x)
      Q(2) = uinit+uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)+uy*Cos(auy*Pi*y)+ux*Sin(aux*Pi*x)
      Q(3) = vinit+vx*Cos(avx*Pi*x)+vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)+vy*Sin(avy*Pi*y)
      Q(4) = 0.0_RP
      Q(5) = (pinit+px*Cos(apx*Pi*x)+pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y)+py*Sin(apy*Pi*y))/&
             ((-1+gamma)*(rhoinit+rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y)+&
              rhoy*Cos(arhoy*Pi*y)+rhox*Sin(arhox*Pi*x)))+&
             ((uinit+uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)+uy*Cos(auy*Pi*y)+&
              ux*Sin(aux*Pi*x))**2+(vinit+vx*Cos(avx*Pi*x)+vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)+&
              vy*Sin(avy*Pi*y))**2)/2
      Q(6) = thetainit+thetax*Cos(athetax*Pi*x)+thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y)+thetay*Cos(athetay*Pi*y)
      
      end associate

   END SUBROUTINE ManufacturedSolutionState 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
subroutine ManufacturedSolutionSourceNSSA(xvec, dwall, time, S)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
            use SMConstants
            IMPLICIT NONE
            real(kind=RP),             intent(in)       :: xvec(NDIM)
            real(kind=RP),             intent(in)       :: dwall
            real(kind=RP),             intent(in)       :: time
            real(kind=RP),             intent(inout)    :: S(NCONS)

            ! --------- my ----------
            real(kind=RP)             :: Q(NCONS)
            real(kind=RP) :: cv2, sbar, somega, S51, S52, cwone, tRatio
            real(kind=RP) :: rhoinit, rhox, rhoy, rhoxy, arhox, arhoy, arhoxy 
            real(kind=RP) :: uinit, ux, uy, uxy, aux, auy, auxy 
            real(kind=RP) :: vinit, vx, vy, vxy, avx, avy, avxy 
            real(kind=RP) :: pinit, px, py, pxy, apx, apy, apxy 
            real(kind=RP) :: thetainit, thetax, thetay, thetaxy, athetax, athetay, athetaxy 
            !real(kind=RP), pointer :: x, y, z, gamma, Mach, Re, Pr, Prt, gammaM2


!            rhoinit = 1.0_RP
!            rhox = 0.1_RP
!            rhoy = -0.2_RP
!            rhoxy = 0.1_RP
!            arhox = 1.0_RP
!            arhoy =1.0_RP
!            arhoxy = 1.0_RP!

!            uinit = 2.0_RP
!            ux = 0.3_RP
!            uy = 0.3_RP
!            uxy = 0.3_RP
!            aux = 3.0_RP
!            auy =1.0_RP
!            auxy =1.0_RP!

!            vinit = 2.0_RP
!            vx = 0.3_RP
!            vy = 0.3_RP
!            vxy = 0.3_RP
!            avx = 1.0_RP
!            avy =1.0_RP
!            avxy = 1.0_RP!

!            pinit = 10.0_RP
!            px = 1.0_RP
!            py = 1.0_RP
!            pxy = 0.5_RP
!            apx = 2.0_RP
!            apy = 1.0_RP
!            apxy = 1.0_RP!
!

!            thetainit = 0.6_RP
!            thetax = -0.03_RP
!            thetay = -0.02_RP
!            thetaxy = 0.02_RP
!            athetax = 2.0_RP
!            athetay = 1.0_RP
!            athetaxy = 3.0_RP!

!            cv2 = 0.7_RP
!            cwone = 0.1355_RP / (0.41_RP*0.41_RP) + (1.0_RP + 0.622_RP)/ (2.0_RP/3.0_RP)!

!               
!!      associate ( x => xvec(1), &
!!                  y => xvec(2), &
!!                  z => xvec(3), &
!!                  gamma => thermodynamics % gamma, &
!!                  Mach => dimensionless % Mach, &
!!                  Re => dimensionless % Re, &
!!                  Pr => dimensionless % Pr, &
!!                  Prt => dimensionless % Prt, &
!!                  gammaM2 => dimensionless % gammaM2, &
!!                  tRatio => S_div_Tref_Sutherland )
!!           -------------
!  
!            somega =  Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                              avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                              auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + auy*Pi*uy*Sin(auy*Pi*y))**2) !

!            sbar   =  (2.4390243902439024*(thetainit + thetax*Cos(athetax*Pi*x) + &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))&
!                      *(1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y))*&
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                        (gammaM2*(1 + tRatio)*&
!                          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                            py*Sin(apy*Pi*y))*&
!                          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                               (tRatio + (gammaM2*&
!                                     (pinit + px*Cos(apx*Pi*x) + &
!                                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                             (gammaM2**6*(1 + tRatio)**4*&
!                               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                               (357.91099999999994 + &
!                                 ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                       thetay*Cos(athetay*Pi*y))**3*&
!                                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                    (tRatio + &
!                                       (gammaM2*&
!                                          (pinit + px*Cos(apx*Pi*x) + &
!                                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)&
!                                            ))/&
!                                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                    (pinit + px*Cos(apx*Pi*x) + &
!                                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3*&
!                                    ((gammaM2*&
!                                         (pinit + px*Cos(apx*Pi*x) + &
!                                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))&
!                                         )/&
!                                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))))))/dwall**2 
!                       
!      call ManufacturedSolutionState(xvec, Q)

      !if (Q(6) >= 0.0_RP) then!

!            S(1)=         -((arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!          (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) + &
!            ux*Sin(aux*Pi*x))) - (rhoinit + &
!          rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(aux*Pi*ux*Cos(aux*Pi*x) - &
!          auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x)) - &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(avy*Pi*vy*Cos(avy*Pi*y) - &
!          avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) - &
!       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y))*&
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) + &
!          vy*Sin(avy*Pi*y)) !

!        S(2)=          apx*Pi*px*Sin(apx*Pi*x) + apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x) -  &
!       (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!           ux*Sin(aux*Pi*x))**2 -  &
!       2*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(aux*Pi*ux*Cos(aux*Pi*x) -  &
!          auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x)) +  &
!       (((gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**4* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!            (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!            (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!            (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!            (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                 athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**4* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                 (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                 (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                    ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                           apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                      (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                         (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                           arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                 (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                    ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                           apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                      (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                         (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                           arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                  (2.*gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!                 (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**2* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                    (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                      athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2))* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/Re &
!         - (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(avy*Pi*vy*Cos(avy*Pi*y) -  &
!          avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) +  &
!       (((gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!          (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!            auy**2*Pi**2*uy*Cos(auy*Pi*y) +  &
!            avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/Re +  &
!       (2*(-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!             aux**2*Pi**2*ux*Sin(aux*Pi*x))* &
!           ((gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))) +  &
!          2*(aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!           ((gammaM2*(1 + tRatio)* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**4* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**3* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                  athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                     (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**2* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                       athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                   ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                         thetay*Cos(athetay*Pi*y))**3* &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                      (tRatio + (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                    (gammaM2**3*(1 + tRatio)**3* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                      ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2)) &
!            - (2*((gammaM2*(1 + tRatio)* &
!                  (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!               (gammaM2*(1 + tRatio)* &
!                  (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!               (gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                  ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                         arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!               (gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                         arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!               (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**4* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!               (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                  (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!               (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**3* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                    athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!               (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                  ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                         arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!               (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                  ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                         arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                (2.*gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!               ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                  ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                    (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                       (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                         arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                    (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**2* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                         athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                    (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                       ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                              apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                         (gammaM2* &
!                            (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                              arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                    (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                       ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                              apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                         (gammaM2* &
!                            (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                              arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                     (2.*gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2 &
!                  ))*(aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/3. -  &
!          (2*((gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!               ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!             (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!               aux**2*Pi**2*ux*Sin(aux*Pi*x) +  &
!               avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/3.)/Re -  &
!       (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y)) - (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!        (-(auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y)) - auy*Pi*uy*Sin(auy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y))        !
!
!
!
!

!       S(3) =               -(apy*Pi*py*Cos(apy*Pi*y)) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!          uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!        (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x)) +  &
!       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y) +  &
!       (((gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!          (-(avx**2*Pi**2*vx*Cos(avx*Pi*x)) -  &
!            avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!            auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y)))/Re +  &
!       (((gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            (gammaM2*(1 + tRatio)* &
!               (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!            (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**4* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!            (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!            (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                 athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!               (357.91099999999994 +  &
!                 ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**4* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**4* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                 (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                    (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                 (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**2* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                      athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                 (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                    ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                           apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                      (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                           arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                  (gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                 (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                       thetay*Cos(athetay*Pi*y))**3* &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                    (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                    ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                           apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                      (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                           arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                  (2.*gammaM2**3*(1 + tRatio)**3* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                       py*Sin(apy*Pi*y))**3* &
!                    ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!               (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2))* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/Re &
!         - (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(vinit + vx*Cos(avx*Pi*x) +  &
!          vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) + vy*Sin(avy*Pi*y)) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(aux*Pi*ux*Cos(aux*Pi*x) -  &
!          auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y)) - 2*(rhoinit +  &
!          rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(avy*Pi*vy*Cos(avy*Pi*y) -  &
!          avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y)) - (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)* &
!             Sin(arhoxy*Pi*y)) - arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!           vy*Sin(avy*Pi*y))**2 +  &
!       (2*((gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**4* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**3* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                  athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**2* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                          Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) &
!                     )/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                   ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                         thetay*Cos(athetay*Pi*y))**3* &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                      (tRatio + (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                    (gammaM2**3*(1 + tRatio)**3* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                      ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2)) &
!            *(avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) &
!           - (2*((gammaM2*(1 + tRatio)* &
!                  (apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!               (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**4* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!               (gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!               (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!               (gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                  ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                         arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!               (gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                         arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!               (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                  ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                         arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!               (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                  ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                    (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                         arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                (2.*gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!               (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**3* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                  (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                    athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!               ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                  ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (apy*Pi*py*Cos(apy*Pi*y) -  &
!                         apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                    (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                         arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                    (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                       ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                              apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                         (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                            (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                              arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                    (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                       ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                              apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                         (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                            (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                              arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                     (2.*gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!                    (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**2* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                       (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                            Sin(athetaxy*Pi*y)) -  &
!                         athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2 &
!                  ))*(aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/3. +  &
!          2*((gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!           (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) -  &
!             avy**2*Pi**2*vy*Sin(avy*Pi*y)) -  &
!          (2*((gammaM2*(1 + tRatio)* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))* &
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                  (tRatio + (gammaM2* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!               ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**4* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                (gammaM2**3*(1 + tRatio)**3* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3* &
!                  ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                  (357.91099999999994 +  &
!                    ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                          thetay*Cos(athetay*Pi*y))**3* &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                       (tRatio + (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                     (gammaM2**3*(1 + tRatio)**3* &
!                       (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                       ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!             (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) +  &
!               auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y) -  &
!               avy**2*Pi**2*vy*Sin(avy*Pi*y)))/3.)/Re !
!
!
!
!
!

!       S(4) = 0!
!
!
!
!
!
!

!        ! 
!        S(5)=            (((gammaM2*(-(apx**2*Pi**2*px*Cos(apx*Pi*x)) -  &
!                  apxy**2*Pi**2*pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!             (2*gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2 -  &
!             (gammaM2*(-(arhoxy**2*Pi**2*rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y)) -  &
!                  arhox**2*Pi**2*rhox*Sin(arhox*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2 +  &
!             (2*gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                   arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))**2* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**3)* &
!           ((gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))) +  &
!          ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!             (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!           ((gammaM2*(1 + tRatio)* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**4* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**3* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                  athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                     (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**2* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                       athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                   ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                         thetay*Cos(athetay*Pi*y))**3* &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                      (tRatio + (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                    (gammaM2**3*(1 + tRatio)**3* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                      ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2)) &
!            + (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x))* &
!           ((gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!           (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!             auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)) +  &
!          (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!           (2*(aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!              ((gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))) &
!              - (2*((gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                  )*(aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/3.) +  &
!          (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!             ux*Sin(aux*Pi*x))*(2* &
!              (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!                aux**2*Pi**2*ux*Sin(aux*Pi*x))* &
!              ((gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))) &
!              + 2*(aux*Pi*ux*Cos(aux*Pi*x) -  &
!                auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!              ((gammaM2*(1 + tRatio)* &
!                   (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                (gammaM2*(1 + tRatio)* &
!                   (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                     arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                (gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                   ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!                (gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**4* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 + (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                   (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                     arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 + (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**3* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                     athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 + (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                   ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 - (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                   ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 (2.*gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 - ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                   ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                     (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                        (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                     (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**2* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                          athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x)) &
!                         *(tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                     (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                        ((gammaM2* &
!                             (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                               apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                          (gammaM2* &
!                             (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                               arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                     (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                        ((gammaM2* &
!                             (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                               apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                          (gammaM2* &
!                             (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                               arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                      (2.*gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                      ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                            thetay*Cos(athetay*Pi*y))**3* &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                         (tRatio +  &
!                            (gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                       (gammaM2**3*(1 + tRatio)**3* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3 &
!                          *((gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))** &
!                    2)) - (2*((gammaM2*(1 + tRatio)* &
!                     (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                  (gammaM2*(1 + tRatio)* &
!                     (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                  (gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!                  (gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    + (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    + (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                       athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    + (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    - (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    - ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           4*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        + (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                          (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        + (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**2* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)* &
!                             Sin(athetaxy*Pi*x))* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        + (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                          ((gammaM2* &
!                               (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                            (gammaM2* &
!                               (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        - (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                          ((gammaM2* &
!                               (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                            (gammaM2* &
!                               (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                        (2.*gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5))) &
!                    /(gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                        ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                              thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                              thetay*Cos(athetay*Pi*y))**3* &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                           (tRatio +  &
!                              (gammaM2* &
!                                 (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                         (gammaM2**3*(1 + tRatio)**3* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             **3*((gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)) &
!                       **2))*(aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/3. -  &
!             (2*((gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                  )*(-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!                  aux**2*Pi**2*ux*Sin(aux*Pi*x) +  &
!                  avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/3.) +  &
!          ((gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!           (-(avx**2*Pi**2*vx*Cos(avx*Pi*x)) -  &
!             avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!             auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y))* &
!           (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!             vy*Sin(avy*Pi*y)) + ((gammaM2*(1 + tRatio)* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) + &
!                (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**4* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**3* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                  athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                     (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**2* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                       athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                   ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                         thetay*Cos(athetay*Pi*y))**3* &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                      (tRatio + (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                    (gammaM2**3*(1 + tRatio)**3* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                      ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2)) &
!            *(-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!             auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))* &
!           (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!             vy*Sin(avy*Pi*y)))/Re +  &
!              (((gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!           (-((gammaM2*(-(arhoxy**2*Pi**2*rhoxy*Cos(arhoxy*Pi*x)* &
!                       Cos(arhoxy*Pi*y)) - arhoy**2*Pi**2*rhoy*Cos(arhoy*Pi*y))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2) +  &
!             (gammaM2*(-(apxy**2*Pi**2*pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y)) -  &
!                  apy**2*Pi**2*py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!             (2*gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2 +  &
!             (2*gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                   arhoy*Pi*rhoy*Sin(arhoy*Pi*y))**2)/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**3) +  &
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!             (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!           ((gammaM2*(1 + tRatio)* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**4* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(-1 + gamma)*Mach**2*Pr* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**3* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                  athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**2* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                          Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) &
!                     )/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!              ((-1 + gamma)*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                   ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                         thetay*Cos(athetay*Pi*y))**3* &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                      (tRatio + (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                    (gammaM2**3*(1 + tRatio)**3* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                      ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2)) &
!            + (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!             ux*Sin(aux*Pi*x))*((gammaM2*(1 + tRatio)* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**4* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              (2.*gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) +  &
!             (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**3* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                  athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) -  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**2* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                          Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) &
!                     )/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                   ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                         thetay*Cos(athetay*Pi*y))**3* &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                      (tRatio + (gammaM2* &
!                            (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))) &
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                    (gammaM2**3*(1 + tRatio)**3* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                      ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))**2)) &
!            *(-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!             auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)) +  &
!          ((gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!           (-(auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y)) - auy*Pi*uy*Sin(auy*Pi*y))* &
!           (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!             auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)) +  &
!          (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!             ux*Sin(aux*Pi*x))*((gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))**4* &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!              (gammaM2**3*(1 + tRatio)**3* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y))**3* &
!                ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                (357.91099999999994 +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!           (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!             auy**2*Pi**2*uy*Cos(auy*Pi*y) +  &
!             avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)) +  &
!          (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!           (2*((gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!              (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) &
!              - (2*((gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                  )*(aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/3.) +  &
!          (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!             vy*Sin(avy*Pi*y))*(2* &
!              ((gammaM2*(1 + tRatio)* &
!                   (apy*Pi*py*Cos(apy*Pi*y) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (apy*Pi*py*Cos(apy*Pi*y) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**4* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 - (gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                   (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                     arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                   (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                     arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 - (gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                   ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!                (gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                   ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 - (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                   ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                     (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                 (2.*gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 + (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**3* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                   (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                     athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                 - ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                   ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (apy*Pi*py*Cos(apy*Pi*y) -  &
!                          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                     (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!                     (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                        ((gammaM2* &
!                             (apy*Pi*py*Cos(apy*Pi*y) -  &
!                               apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                          (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                              *(-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)* &
!                                 Sin(arhoxy*Pi*y)) - arhoy*Pi*rhoy*Sin(arhoy*Pi*y))) &
!                            /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!                     (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                        ((gammaM2* &
!                             (apy*Pi*py*Cos(apy*Pi*y) -  &
!                               apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                          (gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                              *(-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)* &
!                                 Sin(arhoxy*Pi*y)) - arhoy*Pi*rhoy*Sin(arhoy*Pi*y))) &
!                            /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                      (2.*gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!                     (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**2* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                        (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                             Sin(athetaxy*Pi*y)) -  &
!                          athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                      ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                            thetay*Cos(athetay*Pi*y))**3* &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                         (tRatio +  &
!                            (gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                       (gammaM2**3*(1 + tRatio)**3* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3 &
!                          *((gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))** &
!                    2))*(avy*Pi*vy*Cos(avy*Pi*y) -  &
!                avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) -  &
!             (2*((gammaM2*(1 + tRatio)* &
!                     (apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**4* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    - (gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                  (7*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    - (gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!                  (gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                  (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    - (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                       (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                   (2.*gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    + (4*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**3* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                          Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) &
!                     )/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    - ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                     ((-3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (apy*Pi*py*Cos(apy*Pi*y) -  &
!                            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           4*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        + (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                            arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        + (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!                          ((gammaM2* &
!                               (apy*Pi*py*Cos(apy*Pi*y) -  &
!                                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                            (gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y))* &
!                               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)* &
!                                 Sin(arhoxy*Pi*y)) - arhoy*Pi*rhoy*Sin(arhoy*Pi*y))) &
!                              /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) &
!                        - (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                          ((gammaM2* &
!                               (apy*Pi*py*Cos(apy*Pi*y) -  &
!                                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                            (gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y))* &
!                               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)* &
!                                 Sin(arhoxy*Pi*y)) - arhoy*Pi*rhoy*Sin(arhoy*Pi*y))) &
!                              /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!                        (2.*gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) &
!                        + (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**2* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!                          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!                               Sin(athetaxy*Pi*y)) -  &
!                            athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                    /(gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                        ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                              thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                              thetay*Cos(athetay*Pi*y))**3* &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                           (tRatio +  &
!                              (gammaM2* &
!                                 (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                         (gammaM2**3*(1 + tRatio)**3* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             **3*((gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)) &
!                       **2))*(aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/3. +  &
!             2*((gammaM2*(1 + tRatio)* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))* &
!                   Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                 ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                   (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                      thetay*Cos(athetay*Pi*y))**4* &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                   (tRatio + (gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                 (gammaM2**3*(1 + tRatio)**3* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y))**3* &
!                   ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                   (357.91099999999994 +  &
!                     ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                           thetay*Cos(athetay*Pi*y))**3* &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                        (tRatio +  &
!                           (gammaM2* &
!                              (pinit + px*Cos(apx*Pi*x) +  &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y) &
!                                ))/ &
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                      (gammaM2**3*(1 + tRatio)**3* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                        ((gammaM2* &
!                             (pinit + px*Cos(apx*Pi*x) +  &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)) &
!                             )/ &
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))))* &
!              (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) -  &
!                avy**2*Pi**2*vy*Sin(avy*Pi*y)) -  &
!             (2*((gammaM2*(1 + tRatio)* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     Sqrt((gammaM2* &
!                         (pinit + px*Cos(apx*Pi*x) +  &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!                   ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                     (tRatio + (gammaM2* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!                  ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                        thetay*Cos(athetay*Pi*y))**4* &
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**7* &
!                     (tRatio + (gammaM2* &
!                           (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                   (gammaM2**3*(1 + tRatio)**3* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**3* &
!                     ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5* &
!                     (357.91099999999994 +  &
!                       ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                             thetay*Cos(athetay*Pi*y))**3* &
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!                          (tRatio +  &
!                             (gammaM2* &
!                                (pinit + px*Cos(apx*Pi*x) +  &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                  py*Sin(apy*Pi*y)))/ &
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!                        (gammaM2**3*(1 + tRatio)**3* &
!                          (pinit + px*Cos(apx*Pi*x) +  &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))** &
!                           3*((gammaM2* &
!                               (pinit + px*Cos(apx*Pi*x) +  &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                                 py*Sin(apy*Pi*y)))/ &
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))) &
!                  )*(-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) +  &
!                  auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y) -  &
!                  avy**2*Pi**2*vy*Sin(avy*Pi*y)))/3.))/Re -  &
!       (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!        (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!          py*Sin(apy*Pi*y) + (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.)) -  &
!       (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!        (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!          py*Sin(apy*Pi*y) + (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.)) -  &
!       (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x) +  &
!          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) -  &
!             ((arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y)))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2) +  &
!             (2*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))* &
!                 (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                   auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x)) +  &
!                2*(-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!                   avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x))* &
!                 (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y)))/2.) +  &
!          (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!             arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.)) -  &
!       (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y))*(apy*Pi*py*Cos(apy*Pi*y) -  &
!          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y) +  &
!          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((apy*Pi*py*Cos(apy*Pi*y) - apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) -  &
!             ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2) +  &
!             (2*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))* &
!                 (-(auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y)) -  &
!                   auy*Pi*uy*Sin(auy*Pi*y)) +  &
!                2*(avy*Pi*vy*Cos(avy*Pi*y) -  &
!                   avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!                 (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y)))/2.) +  &
!          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!             arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.))       !
!
!
!
!
!
!

!!END S5 THETA>0
!              if (sbar .GT. -cv2 * somega * Re) then!

!               S(6) =         (-6.4948969840282045*(thetainit + thetax*Cos(athetax*Pi*x) + &
!             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!             thetay*Cos(athetay*Pi*y))**2*&
!          (Min(2.0_RP,(5.948839976204641*&
!                (thetainit + thetax*Cos(athetax*Pi*x) + &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                  thetay*Cos(athetay*Pi*y)))/&
!              (dwall**2*Re*((2.4390243902439024*&
!                     (thetainit + thetax*Cos(athetax*Pi*x) + &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                       thetay*Cos(athetay*Pi*y))*&
!                     (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y))*&
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                          (tRatio + &
!                            (gammaM2*&
!                               (pinit + px*Cos(apx*Pi*x) + &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                 py*Sin(apy*Pi*y)))/&
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                        (gammaM2*(1 + tRatio)*&
!                          (pinit + px*Cos(apx*Pi*x) + &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                          Sqrt((gammaM2*&
!                              (pinit + px*Cos(apx*Pi*x) + &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)&
!                                ))/&
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                               (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                             (gammaM2**6*(1 + tRatio)**4*&
!                               (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                               (357.91099999999994 + &
!                                 ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                  Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                      avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                      auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                      auy*Pi*uy*Sin(auy*Pi*y))**2)))) + &
!            0.3*(-Min(2.0_RP,(5.948839976204641*&
!                    (thetainit + thetax*Cos(athetax*Pi*x) + &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                      thetay*Cos(athetay*Pi*y)))/&
!                  (dwall**2*Re*((2.4390243902439024*&
!                         (thetainit + thetax*Cos(athetax*Pi*x) + &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                           thetay*Cos(athetay*Pi*y))*&
!                         (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                thetay*Cos(athetay*Pi*y))*&
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                              (tRatio + &
!                                (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                 (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                            (gammaM2*(1 + tRatio)*&
!                              (pinit + px*Cos(apx*Pi*x) + &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)&
!                                )*&
!                              Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                              (1 + &
!                                ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                 (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                      Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                          avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                          auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                          auy*Pi*uy*Sin(auy*Pi*y))**2)))) + &
!               Min(2.0_RP,(5.948839976204641*&
!                    (thetainit + thetax*Cos(athetax*Pi*x) + &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                      thetay*Cos(athetay*Pi*y)))/&
!                  (dwall**2*Re*((2.4390243902439024*&
!                         (thetainit + thetax*Cos(athetax*Pi*x) + &
!                           thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                           thetay*Cos(athetay*Pi*y))*&
!                         (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                thetay*Cos(athetay*Pi*y))*&
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                              (tRatio + &
!                                (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                 (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                            (gammaM2*(1 + tRatio)*&
!                              (pinit + px*Cos(apx*Pi*x) + &
!                                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)&
!                                )*&
!                              Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                              (1 + &
!                                ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                 (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                      Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                          avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                          auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                          auy*Pi*uy*Sin(auy*Pi*y))**2))))**6))*&
!          (1/(64 + (Min(2.0_RP,(5.948839976204641*&
!                      (thetainit + thetax*Cos(athetax*Pi*x) + &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                        thetay*Cos(athetay*Pi*y)))/&
!                    (dwall**2*Re*((2.4390243902439024*&
!                           (thetainit + thetax*Cos(athetax*Pi*x) + &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                             thetay*Cos(athetay*Pi*y))*&
!                           (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                              (gammaM2*(1 + tRatio)*&
!                                (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                        Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                            avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                            auy*Pi*uy*Sin(auy*Pi*y))**2)))) + &
!                  0.3*(-Min(2.0_RP,(5.948839976204641*&
!                          (thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y)))/&
!                        (dwall**2*Re*&
!                          ((2.4390243902439024*&
!                               (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                 thetay*Cos(athetay*Pi*y))*&
!                               (1 - &
!                                 ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                            Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2)))) + &
!                     Min(2.0_RP,(5.948839976204641*&
!                          (thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y)))/&
!                        (dwall**2*Re*&
!                          ((2.4390243902439024*&
!                               (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                 thetay*Cos(athetay*Pi*y))*&
!                               (1 - &
!                                 ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                            Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2))))**6))**6))**&
!           0.16666666666666666*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))/(dwall**2*Re) - &
!       (thetainit + thetax*Cos(athetax*Pi*x) + &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))*&
!        (arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x)) - &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!          athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x)) + &
!       (3*gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))*&
!          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetax**2*Pi**2*thetax*Cos(athetax*Pi*x)) - &
!            athetaxy**2*Pi**2*thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetaxy**2*Pi**2*thetaxy*Cos(athetaxy*Pi*x)*&
!               Cos(athetaxy*Pi*y)) - athetay**2*Pi**2*thetay*Cos(athetay*Pi*y))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*((gammaM2*&
!               (-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (4.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (-(((thetainit + thetax*Cos(athetax*Pi*x) + &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                   thetay*Cos(athetay*Pi*y))*&
!                 (-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                   apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))*&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                 (tRatio + (gammaM2*&
!                      (pinit + px*Cos(apx*Pi*x) + &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!               (gammaM2*(1 + tRatio)*&
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                    py*Sin(apy*Pi*y))**2*&
!                 Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) + &
!            (2*(thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!               (arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!            ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!                 athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (2.*gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (thetainit + thetax*Cos(athetax*Pi*x) + &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))*&
!        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)) - &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*&
!             Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) + &
!       (3*gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) - &
!            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!            arhoy*Pi*rhoy*Sin(arhoy*Pi*y))*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) + &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (4.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y))*&
!          (-(((thetainit + thetax*Cos(athetax*Pi*x) + &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                   thetay*Cos(athetay*Pi*y))*&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                 (apy*Pi*py*Cos(apy*Pi*y) - &
!                   apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))*&
!                 (tRatio + (gammaM2*&
!                      (pinit + px*Cos(apx*Pi*x) + &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!               (gammaM2*(1 + tRatio)*&
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                    py*Sin(apy*Pi*y))**2*&
!                 Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) + &
!            (2*(thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (2.*gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) + &
!            ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!                 athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (0.933*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          ((-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!               athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))**2 + &
!            (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!               athetay*Pi*thetay*Sin(athetay*Pi*y))**2))/Re + &
!       0.1355*(thetainit + thetax*Cos(athetax*Pi*x) + &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))*&
!        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*((2.4390243902439024*&
!             (thetainit + thetax*Cos(athetax*Pi*x) + &
!               thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!               thetay*Cos(athetay*Pi*y))*&
!             (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                    thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                    thetay*Cos(athetay*Pi*y))*&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                  (tRatio + (gammaM2*&
!                       (pinit + px*Cos(apx*Pi*x) + &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                       rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                (gammaM2*(1 + tRatio)*&
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                    py*Sin(apy*Pi*y))*&
!                  Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                  (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                          thetay*Cos(athetay*Pi*y))**4*&
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                       (tRatio + (gammaM2*&
!                             (pinit + px*Cos(apx*Pi*x) + &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))&
!                             )/&
!                           (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                     (gammaM2**6*(1 + tRatio)**4*&
!                       (pinit + px*Cos(apx*Pi*x) + &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**6*&
!                       (357.91099999999994 + &
!                         ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                               thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                               thetay*Cos(athetay*Pi*y))**3*&
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                            (tRatio + &
!                               (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                          (gammaM2**3*(1 + tRatio)**3*&
!                            (pinit + px*Cos(apx*Pi*x) + &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))&
!                              **3*&
!                            ((gammaM2*&
!                                 (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)&
!                         ))))))/(dwall**2*Re) + &
!          Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!              avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!              auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + auy*Pi*uy*Sin(auy*Pi*y))**2)) 
!!END S6 THETA>0
!    else!

!    S(6)=         (-6.4948969840282045*(thetainit + thetax*Cos(athetax*Pi*x) + &
!             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!             thetay*Cos(athetay*Pi*y))**2*&
!          (Min(2.0_RP,(5.948839976204641*&
!                (thetainit + thetax*Cos(athetax*Pi*x) + &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                  thetay*Cos(athetay*Pi*y)))/&
!              (dwall**2*Re*(Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                      avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                      auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                      auy*Pi*uy*Sin(auy*Pi*y))**2) + &
!                  (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                         avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                         auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                         auy*Pi*uy*Sin(auy*Pi*y))**2)*&
!                     ((2.1951219512195124*&
!                          (thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y))*&
!                          (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                 thetay*Cos(athetay*Pi*y))*&
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                               (tRatio + &
!                                 (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                             (gammaM2*(1 + tRatio)*&
!                               (pinit + px*Cos(apx*Pi*x) + &
!                                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                 py*Sin(apy*Pi*y))*&
!                               Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                 (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                               (1 + &
!                                 ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                       0.48999999999999994*&
!                        Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                            avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                            auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!                   ((-2.4390243902439024*&
!                        (thetainit + thetax*Cos(athetax*Pi*x) + &
!                          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                          thetay*Cos(athetay*Pi*y))*&
!                        (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                               thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                               thetay*Cos(athetay*Pi*y))*&
!                             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                             (tRatio + &
!                               (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                           (gammaM2*(1 + tRatio)*&
!                             (pinit + px*Cos(apx*Pi*x) + &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))&
!                              *Sqrt((gammaM2*&
!                                 (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                             (1 + &
!                               ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) - &
!                     0.4999999999999999*&
!                      Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                          avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                          auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                          auy*Pi*uy*Sin(auy*Pi*y))**2))))) + &
!            0.3*(-Min(2.0_RP,(5.948839976204641*&
!                    (thetainit + thetax*Cos(athetax*Pi*x) + &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                      thetay*Cos(athetay*Pi*y)))/&
!                  (dwall**2*Re*(Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                          avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                          auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                          auy*Pi*uy*Sin(auy*Pi*y))**2) + &
!                      (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                             auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                             auy*Pi*uy*Sin(auy*Pi*y))**2)*&
!                         ((2.1951219512195124*&
!                              (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                thetay*Cos(athetay*Pi*y))*&
!                              (1 - &
!                                ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                 (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                           0.48999999999999994*&
!                            Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!                       ((-2.4390243902439024*&
!                            (thetainit + thetax*Cos(athetax*Pi*x) + &
!                              thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                              thetay*Cos(athetay*Pi*y))*&
!                            (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                 (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                 (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                               (gammaM2*(1 + tRatio)*&
!                                 (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                 Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                 (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) - &
!                         0.4999999999999999*&
!                          Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                              avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                              auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                              auy*Pi*uy*Sin(auy*Pi*y))**2))))) + &
!               Min(2.0_RP,(5.948839976204641*&
!                    (thetainit + thetax*Cos(athetax*Pi*x) + &
!                      thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                      thetay*Cos(athetay*Pi*y)))/&
!                  (dwall**2*Re*(Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                          avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                          auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                          auy*Pi*uy*Sin(auy*Pi*y))**2) + &
!                      (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                             avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                             auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                             auy*Pi*uy*Sin(auy*Pi*y))**2)*&
!                         ((2.1951219512195124*&
!                              (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                thetay*Cos(athetay*Pi*y))*&
!                              (1 - &
!                                ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                 (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                           0.48999999999999994*&
!                            Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!                       ((-2.4390243902439024*&
!                            (thetainit + thetax*Cos(athetax*Pi*x) + &
!                              thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                              thetay*Cos(athetay*Pi*y))*&
!                            (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                 (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                 (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                               (gammaM2*(1 + tRatio)*&
!                                 (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                 Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                 (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) - &
!                         0.4999999999999999*&
!                          Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                              avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                              auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                              auy*Pi*uy*Sin(auy*Pi*y))**2)))))**6))*&
!          (1/(64 + (Min(2.0_RP,(5.948839976204641*&
!                      (thetainit + thetax*Cos(athetax*Pi*x) + &
!                        thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                        thetay*Cos(athetay*Pi*y)))/&
!                    (dwall**2*Re*(Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                            avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                            auy*Pi*uy*Sin(auy*Pi*y))**2) + &
!                        (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                               auy*Pi*uy*Sin(auy*Pi*y))**2)*&
!                           ((2.1951219512195124*&
!                                (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                (1 - &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                             0.48999999999999994*&
!                              Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!                         ((-2.4390243902439024*&
!                              (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                thetay*Cos(athetay*Pi*y))*&
!                              (1 - &
!                                ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                 (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) - &
!                           0.4999999999999999*&
!                            Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2))))) + &
!                  0.3*(-Min(2.0_RP,(5.948839976204641*&
!                          (thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y)))/&
!                        (dwall**2*Re*&
!                          (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2) + &
!                            (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2)*&
!                               ((2.1951219512195124*&
!                                  (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (1 - &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                                 0.48999999999999994*&
!                                  Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!                             ((-2.4390243902439024*&
!                                  (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (1 - &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) - &
!                               0.4999999999999999*&
!                                Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2))))) + &
!                     Min(2.0_RP,(5.948839976204641*&
!                          (thetainit + thetax*Cos(athetax*Pi*x) + &
!                            thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                            thetay*Cos(athetay*Pi*y)))/&
!                        (dwall**2*Re*&
!                          (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                auy*Pi*uy*Sin(auy*Pi*y))**2) + &
!                            (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2)*&
!                               ((2.1951219512195124*&
!                                  (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (1 - &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) + &
!                                 0.48999999999999994*&
!                                  Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!                             ((-2.4390243902439024*&
!                                  (thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (1 - &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                                  (gammaM2*(1 + tRatio)*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))*&
!                                  Sqrt((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                                  (1 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**4*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                                  (gammaM2**6*(1 + tRatio)**4*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**6*&
!                                  (357.91099999999994 + &
!                                  ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                  (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                                  (gammaM2**3*(1 + tRatio)**3*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                  ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                  )))))))/(dwall**2*Re) - &
!                               0.4999999999999999*&
!                                Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                                  auy*Pi*uy*Sin(auy*Pi*y))**2)))))**6))**6))**&
!           0.16666666666666666*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))/(dwall**2*Re) - &
!       (thetainit + thetax*Cos(athetax*Pi*x) + &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))*&
!        (arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x)) - &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!          athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x)) + &
!       (3*gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))*&
!          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetax**2*Pi**2*thetax*Cos(athetax*Pi*x)) - &
!            athetaxy**2*Pi**2*thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetaxy**2*Pi**2*thetaxy*Cos(athetaxy*Pi*x)*&
!               Cos(athetaxy*Pi*y)) - athetay**2*Pi**2*thetay*Cos(athetay*Pi*y))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*((gammaM2*&
!               (-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (4.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2*&
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (-(((thetainit + thetax*Cos(athetax*Pi*x) + &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                   thetay*Cos(athetay*Pi*y))*&
!                 (-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                   apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))*&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                 (tRatio + (gammaM2*&
!                      (pinit + px*Cos(apx*Pi*x) + &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!               (gammaM2*(1 + tRatio)*&
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                    py*Sin(apy*Pi*y))**2*&
!                 Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) + &
!            (2*(thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!               (arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!            ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!                 athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (2.*gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (thetainit + thetax*Cos(athetax*Pi*x) + &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))*&
!        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)) - &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*&
!             Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) + &
!       (3*gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) - &
!            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!            arhoy*Pi*rhoy*Sin(arhoy*Pi*y))*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) + &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))*&
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!        (4.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) + &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!            athetay*Pi*thetay*Sin(athetay*Pi*y))*&
!          (-(((thetainit + thetax*Cos(athetax*Pi*x) + &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                   thetay*Cos(athetay*Pi*y))*&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                 (apy*Pi*py*Cos(apy*Pi*y) - &
!                   apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))*&
!                 (tRatio + (gammaM2*&
!                      (pinit + px*Cos(apx*Pi*x) + &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!               (gammaM2*(1 + tRatio)*&
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                    py*Sin(apy*Pi*y))**2*&
!                 Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) + &
!            (2*(thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) - &
!            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                 thetay*Cos(athetay*Pi*y))*&
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) - &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) - &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/&
!             (2.*gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) + &
!            ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!               (tRatio + (gammaM2*&
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                      py*Sin(apy*Pi*y)))/&
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!               (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!                 athetay*Pi*thetay*Sin(athetay*Pi*y)))/&
!             (gammaM2*(1 + tRatio)*&
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                 py*Sin(apy*Pi*y))*&
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/&
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) + &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) + &
!       (0.933*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))*&
!          ((-(athetax*Pi*thetax*Sin(athetax*Pi*x)) - &
!               athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))**2 + &
!            (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) - &
!               athetay*Pi*thetay*Sin(athetay*Pi*y))**2))/Re + &
!       0.1355*(thetainit + thetax*Cos(athetax*Pi*x) + &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))*&
!        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!              avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!              auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + auy*Pi*uy*Sin(auy*Pi*y))**&
!            2) + (Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                 avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                 auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + auy*Pi*uy*Sin(auy*Pi*y)&
!                 )**2)*((2.1951219512195124*&
!                  (thetainit + thetax*Cos(athetax*Pi*x) + &
!                    thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                    thetay*Cos(athetay*Pi*y))*&
!                  (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                         thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                         thetay*Cos(athetay*Pi*y))*&
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                       (tRatio + (gammaM2*&
!                            (pinit + px*Cos(apx*Pi*x) + &
!                              pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))&
!                           /(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                     (gammaM2*(1 + tRatio)*&
!                       (pinit + px*Cos(apx*Pi*x) + &
!                         pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                       Sqrt((gammaM2*&
!                           (pinit + px*Cos(apx*Pi*x) + &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                         (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                           rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                       (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                               thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                               thetay*Cos(athetay*Pi*y))**4*&
!                            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                            (tRatio + &
!                               (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                          (gammaM2**6*(1 + tRatio)**4*&
!                            (pinit + px*Cos(apx*Pi*x) + &
!                               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))&
!                              **6*&
!                            (357.91099999999994 + &
!                              ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                                 (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                                 (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                               (gammaM2**3*(1 + tRatio)**3*&
!                                 (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                                 ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                                 )))))))/(dwall**2*Re) + &
!               0.48999999999999994*&
!                Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                    avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                    auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                    auy*Pi*uy*Sin(auy*Pi*y))**2)))/&
!           ((-2.4390243902439024*(thetainit + thetax*Cos(athetax*Pi*x) + &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                  thetay*Cos(athetay*Pi*y))*&
!                (1 - ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                       thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                       thetay*Cos(athetay*Pi*y))*&
!                     (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2*&
!                     (tRatio + (gammaM2*&
!                          (pinit + px*Cos(apx*Pi*x) + &
!                            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/&
!                   (gammaM2*(1 + tRatio)*&
!                     (pinit + px*Cos(apx*Pi*x) + &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))*&
!                     Sqrt((gammaM2*&
!                         (pinit + px*Cos(apx*Pi*x) + &
!                           pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/&
!                       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                         rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))*&
!                     (1 + ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                             thetay*Cos(athetay*Pi*y))**4*&
!                          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**10*&
!                          (tRatio + &
!                             (gammaM2*&
!                                (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**4)/&
!                        (gammaM2**6*(1 + tRatio)**4*&
!                          (pinit + px*Cos(apx*Pi*x) + &
!                             pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))**&
!                           6*(357.91099999999994 + &
!                            ((thetainit + thetax*Cos(athetax*Pi*x) + &
!                                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + &
!                                  thetay*Cos(athetay*Pi*y))**3*&
!                               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6*&
!                               (tRatio + &
!                                  (gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/&
!                             (gammaM2**3*(1 + tRatio)**3*&
!                               (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y))**3*&
!                               ((gammaM2*&
!                                  (pinit + px*Cos(apx*Pi*x) + &
!                                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + &
!                                  py*Sin(apy*Pi*y)))/&
!                                  (rhoinit + &
!                                  rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + &
!                                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5&
!                               )))))))/(dwall**2*Re) - &
!             0.4999999999999999*Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) - &
!                  avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) + &
!                  auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + &
!                  auy*Pi*uy*Sin(auy*Pi*y))**2))) 
!             end if!

!      else!

!      S(1)=         -((arhox*Pi*rhox*Cos(arhox*Pi*x) - &
!            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))*&
!          (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) + &
!            ux*Sin(aux*Pi*x))) - (rhoinit + &
!          rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(aux*Pi*ux*Cos(aux*Pi*x) - &
!          auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x)) - &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) + &
!          rhox*Sin(arhox*Pi*x))*(avy*Pi*vy*Cos(avy*Pi*y) - &
!          avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) - &
!       (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) - &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y))*&
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) + &
!          vy*Sin(avy*Pi*y)) !
!
!
!
!
!
!
!
!

!           S(2) =                apx*Pi*px*Sin(apx*Pi*x) + apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x) -  &
!       (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!           ux*Sin(aux*Pi*x))**2 -  &
!       2*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(aux*Pi*ux*Cos(aux*Pi*x) -  &
!          auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x)) +  &
!       (gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) -  &
!            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!            arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!          uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!        (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)) +  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!            auy**2*Pi**2*uy*Cos(auy*Pi*y) +  &
!            avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       ((2*gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!               apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!             (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (2*gammaM2*(1 + tRatio)* &
!             (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!               aux**2*Pi**2*ux*Sin(aux*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!               arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!             (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!          (gammaM2*(1 + tRatio)*(aux*Pi*ux*Cos(aux*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*((gammaM2* &
!                  (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!               apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (2*gammaM2*(1 + tRatio)* &
!             (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!               arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) -  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!               aux**2*Pi**2*ux*Sin(aux*Pi*x) +  &
!               avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))/Re -  &
!       (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y)) - (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!          rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!        (-(auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y)) - auy*Pi*uy*Sin(auy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y))  !
!
!
!
!
!
!
!
!
!
!

!        S(3)=                -(apy*Pi*py*Cos(apy*Pi*y)) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!          uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!        (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x)) +  &
!       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y) +  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(avx**2*Pi**2*vx*Cos(avx*Pi*x)) -  &
!            avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!            auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       (gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (gammaM2*(1 + tRatio)*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!       (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (-(avx*Pi*vx*Sin(avx*Pi*x)) - avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!        (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(vinit + vx*Cos(avx*Pi*x) +  &
!          vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) + vy*Sin(avy*Pi*y)) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(aux*Pi*ux*Cos(aux*Pi*x) -  &
!          auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y)) - 2*(rhoinit +  &
!          rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(avy*Pi*vy*Cos(avy*Pi*y) -  &
!          avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y)) - (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)* &
!             Sin(arhoxy*Pi*y)) - arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!        (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!           vy*Sin(avy*Pi*y))**2 +  &
!       ((2*gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) -  &
!               apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!               arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!             (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (apy*Pi*py*Cos(apy*Pi*y) - apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!               arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) -  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) -  &
!               avy**2*Pi**2*vy*Sin(avy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (2*gammaM2*(1 + tRatio)* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) +  &
!               auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y) -  &
!               avy**2*Pi**2*vy*Sin(avy*Pi*y)))/ &
!           (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))/Re  !
!
!
!
!
!
!
!
!

!        S(4) = 0.0_RP!
!

!        S(5)=            ((gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!               apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!               arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                     arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)**2)/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                     arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                   (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)**2)/ &
!           (2.*(-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx**2*Pi**2*px*Cos(apx*Pi*x)) -  &
!                    apxy**2*Pi**2*pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (2*gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                  (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2 -  &
!               (gammaM2*(-(arhoxy**2*Pi**2*rhoxy*Cos(arhoxy*Pi*x)* &
!                       Cos(arhoxy*Pi*y)) - arhox**2*Pi**2*rhox*Sin(arhox*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2 +  &
!               (2*gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                     arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))**2* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**3))/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))) &
!            /((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!           ((2*gammaM2*(1 + tRatio)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) +  &
!          (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!             ux*Sin(aux*Pi*x))*((2*gammaM2*(1 + tRatio)* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (2*gammaM2*(1 + tRatio)* &
!                (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!                  aux**2*Pi**2*ux*Sin(aux*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (2*gammaM2*(1 + tRatio)* &
!                (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                       arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!                  aux**2*Pi**2*ux*Sin(aux*Pi*x) +  &
!                  avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avx**2*Pi**2*vx*Cos(avx*Pi*x)) -  &
!               avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!               auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y))* &
!             (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!               vy*Sin(avy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!               apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))* &
!             (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!               vy*Sin(avy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!               arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))* &
!             (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!               vy*Sin(avy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))* &
!             (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!               vy*Sin(avy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                    arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))* &
!             (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!               vy*Sin(avy*Pi*y)))/ &
!           (2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))))/Re +  &
!       ((gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) -  &
!               apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!               arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                   (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                     arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)**2)/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                     apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                   (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                     arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)**2)/ &
!           (2.*(-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-((gammaM2*(-(arhoxy**2*Pi**2*rhoxy*Cos(arhoxy*Pi*x)* &
!                         Cos(arhoxy*Pi*y)) - arhoy**2*Pi**2*rhoy*Cos(arhoy*Pi*y))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2) +  &
!               (gammaM2*(-(apxy**2*Pi**2*pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y)) -  &
!                    apy**2*Pi**2*py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (2*gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2 +  &
!               (2*gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                     arhoy*Pi*rhoy*Sin(arhoy*Pi*y))**2)/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**3))/ &
!           ((-1 + gamma)*Mach**2*Pr* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!               uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!             (apy*Pi*py*Cos(apy*Pi*y) - apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))) &
!            /((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!               uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!               arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))) &
!            /((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!          (gammaM2*(1 + tRatio)*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!               uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))) &
!            /((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!          (gammaM2*(1 + tRatio)*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!               uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*((gammaM2* &
!                  (apy*Pi*py*Cos(apy*Pi*y) -  &
!                    apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!               (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                  (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                    arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))) &
!            /(2.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!               pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!             Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                   pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y)) -  &
!               auy*Pi*uy*Sin(auy*Pi*y))* &
!             (-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!               avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) -  &
!               auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) - auy*Pi*uy*Sin(auy*Pi*y))) &
!            /((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (gammaM2*(1 + tRatio)*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) +  &
!               uy*Cos(auy*Pi*y) + ux*Sin(aux*Pi*x))* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!               py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                   py*Sin(apy*Pi*y)))/ &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!             (-(auxy**2*Pi**2*uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y)) -  &
!               auy**2*Pi**2*uy*Cos(auy*Pi*y) +  &
!               avxy**2*Pi**2*vxy*Sin(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!           ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                    pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!          (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!           ((2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) +  &
!          (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!             vy*Sin(avy*Pi*y))*((2*gammaM2*(1 + tRatio)* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!                (avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!                (avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!                (avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (apy*Pi*py*Cos(apy*Pi*y) -  &
!                  apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) -  &
!             (gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                       apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                  (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                     (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                       arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!                (aux*Pi*ux*Cos(aux*Pi*x) + avy*Pi*vy*Cos(avy*Pi*y) -  &
!                  auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x) -  &
!                  avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) -  &
!                  avy**2*Pi**2*vy*Sin(avy*Pi*y)))/ &
!              ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!             (2*gammaM2*(1 + tRatio)* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!                (-(avxy**2*Pi**2*vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y)) +  &
!                  auxy**2*Pi**2*uxy*Sin(auxy*Pi*x)*Sin(auxy*Pi*y) -  &
!                  avy**2*Pi**2*vy*Sin(avy*Pi*y)))/ &
!              (3.*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!                (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))))/Re -  &
!       (aux*Pi*ux*Cos(aux*Pi*x) - auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x))* &
!        (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!          py*Sin(apy*Pi*y) + (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.)) -  &
!       (avy*Pi*vy*Cos(avy*Pi*y) - avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!        (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!          py*Sin(apy*Pi*y) + (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.)) -  &
!       (uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!          ux*Sin(aux*Pi*x))*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!          apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x) +  &
!          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) -  &
!             ((arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                  arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y)))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2) +  &
!             (2*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))* &
!                 (aux*Pi*ux*Cos(aux*Pi*x) -  &
!                   auxy*Pi*uxy*Cos(auxy*Pi*y)*Sin(auxy*Pi*x)) +  &
!                2*(-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!                   avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x))* &
!                 (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y)))/2.) +  &
!          (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!             arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.)) -  &
!       (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!          vy*Sin(avy*Pi*y))*(apy*Pi*py*Cos(apy*Pi*y) -  &
!          apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y) +  &
!          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!           ((apy*Pi*py*Cos(apy*Pi*y) - apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) -  &
!             ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))* &
!                (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                  arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2) +  &
!             (2*(uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))* &
!                 (-(auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y)) -  &
!                   auy*Pi*uy*Sin(auy*Pi*y)) +  &
!                2*(avy*Pi*vy*Cos(avy*Pi*y) -  &
!                   avxy*Pi*vxy*Cos(avxy*Pi*x)*Sin(avxy*Pi*y))* &
!                 (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y)))/2.) +  &
!          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!             arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!           ((pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))/ &
!              ((-1 + gamma)*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))) +  &
!             ((uinit + uxy*Cos(auxy*Pi*x)*Cos(auxy*Pi*y) + uy*Cos(auy*Pi*y) +  &
!                   ux*Sin(aux*Pi*x))**2 +  &
!                (vinit + vx*Cos(avx*Pi*x) + vxy*Cos(avxy*Pi*x)*Cos(avxy*Pi*y) +  &
!                   vy*Sin(avy*Pi*y))**2)/2.))         !
!
!
!
!

!      S(6)=  (cwone*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!             thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!             thetay*Cos(athetay*Pi*y))**2* &
!          (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))/(dwall**2*Re) -  &
!       (thetainit + thetax*Cos(athetax*Pi*x) +  &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))* &
!        (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!          arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x)) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!          athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x)) +  &
!       (3*gammaM2*(1 + tRatio)*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!            apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       (3*gammaM2*(1 + tRatio)*(-(athetax**2*Pi**2*thetax*Cos(athetax*Pi*x)) -  &
!            athetaxy**2*Pi**2*thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       (3*gammaM2*(1 + tRatio)*(-(athetaxy**2*Pi**2*thetaxy*Cos(athetaxy*Pi*x)* &
!               Cos(athetaxy*Pi*y)) - athetay**2*Pi**2*thetay*Cos(athetay*Pi*y))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (3*gammaM2*(1 + tRatio)*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!            arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!          (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*((gammaM2* &
!               (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (4.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       (3*gammaM2*(1 + tRatio)*(-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!            athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!          (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!            py*Sin(apy*Pi*y))*Sqrt((gammaM2* &
!              (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))* &
!                 (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                   apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                 (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!               (gammaM2*(1 + tRatio)* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))**2* &
!                 Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) +  &
!            (2*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                 athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**4* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!               (arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                 arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**2* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!                 athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               ((gammaM2*(-(apx*Pi*px*Sin(apx*Pi*x)) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*y)*Sin(apxy*Pi*x)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(arhox*Pi*rhox*Cos(arhox*Pi*x) -  &
!                      arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*y)*Sin(arhoxy*Pi*x))* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (thetainit + thetax*Cos(athetax*Pi*x) +  &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))* &
!        (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!          arhoy*Pi*rhoy*Sin(arhoy*Pi*y)) -  &
!       (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)* &
!             Sin(athetaxy*Pi*y)) - athetay*Pi*thetay*Sin(athetay*Pi*y)) +  &
!       (3*gammaM2*(1 + tRatio)*(apy*Pi*py*Cos(apy*Pi*y) -  &
!            apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))* &
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))* &
!          (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!            arhoy*Pi*rhoy*Sin(arhoy*Pi*y))* &
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!             rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))* &
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                  pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!              (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2) +  &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5))* &
!          ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!            (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2)* &
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!            athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!        (4.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       (3*gammaM2*(1 + tRatio)*(pinit + px*Cos(apx*Pi*x) +  &
!            pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!          Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!            (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!              rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!          (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!            athetay*Pi*thetay*Sin(athetay*Pi*y))* &
!          (-(((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                   thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                   thetay*Cos(athetay*Pi*y))* &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!                 (apy*Pi*py*Cos(apy*Pi*y) -  &
!                   apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!                 (tRatio + (gammaM2* &
!                      (pinit + px*Cos(apx*Pi*x) +  &
!                        pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                    (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                      rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))/ &
!               (gammaM2*(1 + tRatio)* &
!                 (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                    py*Sin(apy*Pi*y))**2* &
!                 Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))))) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (apy*Pi*py*Cos(apy*Pi*y) -  &
!                 apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y))* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3)/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**4* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            (2*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                 rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            (6*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                 arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) -  &
!            ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                 thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                 thetay*Cos(athetay*Pi*y))* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) +  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5) -  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**3* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               ((gammaM2*(apy*Pi*py*Cos(apy*Pi*y) -  &
!                      apxy*Pi*pxy*Cos(apxy*Pi*x)*Sin(apxy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)) -  &
!                 (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y))* &
!                    (-(arhoxy*Pi*rhoxy*Cos(arhoxy*Pi*x)*Sin(arhoxy*Pi*y)) -  &
!                      arhoy*Pi*rhoy*Sin(arhoy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2))/ &
!             (2.*gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2.5) +  &
!            ((rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**2* &
!               (tRatio + (gammaM2* &
!                    (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                      py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))* &
!               (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                 athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!             (gammaM2*(1 + tRatio)* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                 py*Sin(apy*Pi*y))* &
!               Sqrt((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!            (3*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                  thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                  thetay*Cos(athetay*Pi*y))**2* &
!               (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                  rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**6* &
!               (tRatio + (gammaM2* &
!                     (pinit + px*Cos(apx*Pi*x) +  &
!                       pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                   (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**3* &
!               (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!                 athetay*Pi*thetay*Sin(athetay*Pi*y)))/ &
!             (gammaM2**3*(1 + tRatio)**3* &
!               (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                  py*Sin(apy*Pi*y))**3* &
!               ((gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                      pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                    rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**1.5)))/ &
!        (2.*Re*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                 pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!               rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))) +  &
!       (0.933*(rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!            rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))* &
!          ((-(athetax*Pi*thetax*Sin(athetax*Pi*x)) -  &
!               athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*y)*Sin(athetaxy*Pi*x))**2 +  &
!            (-(athetaxy*Pi*thetaxy*Cos(athetaxy*Pi*x)*Sin(athetaxy*Pi*y)) -  &
!               athetay*Pi*thetay*Sin(athetay*Pi*y))**2))/Re +  &
!       0.1355*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!          thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) + thetay*Cos(athetay*Pi*y))* &
!        (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) + rhoy*Cos(arhoy*Pi*y) +  &
!          rhox*Sin(arhox*Pi*x))*(1 -  &
!          (1000*(thetainit + thetax*Cos(athetax*Pi*x) +  &
!                thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                thetay*Cos(athetay*Pi*y))**2* &
!             (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!             (tRatio + (gammaM2*(pinit + px*Cos(apx*Pi*x) +  &
!                     pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                 (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                   rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2)/ &
!           (gammaM2**3*(1 + tRatio)**2* &
!             (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                py*Sin(apy*Pi*y))**3* &
!             (1 + ((thetainit + thetax*Cos(athetax*Pi*x) +  &
!                     thetaxy*Cos(athetaxy*Pi*x)*Cos(athetaxy*Pi*y) +  &
!                     thetay*Cos(athetay*Pi*y))**2* &
!                  (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                     rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x))**5* &
!                  (tRatio + (gammaM2* &
!                        (pinit + px*Cos(apx*Pi*x) +  &
!                          pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) + py*Sin(apy*Pi*y)))/ &
!                      (rhoinit + rhoxy*Cos(arhoxy*Pi*x)*Cos(arhoxy*Pi*y) +  &
!                        rhoy*Cos(arhoy*Pi*y) + rhox*Sin(arhox*Pi*x)))**2)/ &
!                (gammaM2**3*(1 + tRatio)**2* &
!                  (pinit + px*Cos(apx*Pi*x) + pxy*Cos(apxy*Pi*x)*Cos(apxy*Pi*y) +  &
!                     py*Sin(apy*Pi*y))**3))))* &
!        Sqrt((-(avx*Pi*vx*Sin(avx*Pi*x)) -  &
!            avxy*Pi*vxy*Cos(avxy*Pi*y)*Sin(avxy*Pi*x) +  &
!            auxy*Pi*uxy*Cos(auxy*Pi*x)*Sin(auxy*Pi*y) + auy*Pi*uy*Sin(auy*Pi*y))**2) 
!
!      end if 

!      end associate
   END SUBROUTINE ManufacturedSolutionSourceNSSA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   

END MODULE ManufacturedSolutionsNSSA
