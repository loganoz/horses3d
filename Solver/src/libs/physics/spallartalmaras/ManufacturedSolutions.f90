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
   use VariableConversion_NSSA
   use SpallartAlmarasTurbulence
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
   SUBROUTINE ManufacturedSolutionState( xvec, Q, thetaeddy )
      IMPLICIT NONE
!
!     ----------------------------------------------------
!     Manufactured solution state (valid for Euler and NS)
!     ----------------------------------------------------
!      
      
      REAL(KIND=RP), intent(in)    :: xvec(3)
      REAL(KIND=RP), intent(inout) :: Q(NCONS)
      REAL(KIND=RP), intent(inout) :: thetaeddy
      REAL(KIND=RP)                :: rho , u, v, w, en

                

      associate ( x => xvec(1), &
                  y => xvec(2), &
                  z => xvec(3), &
                  gamma => thermodynamics % gamma)

                   rho = 2.0 + 0.1*Sin(Pi*(x + y + z)) 
                   u = 1.0 + 0.1*Sin(Pi*(x + y + z))
                   v = -4.0 + 2.0*(2.0 + 0.1*Sin(Pi*(x + y + z)))
                   w = 5.0 - 3.0*(2.0 + 0.1*Sin(Pi*(x + y + z)))
                   en = 2.0 + 0.1*Sin(Pi*(x + y + z))
                   thetaeddy = Cos(Pi*(x + y + z)) 


                   Q(1) = rho
                   Q(2) = rho * u
                   Q(3) = rho * v
                   Q(4) = rho * w
                   Q(5) = rho * en
                   Q(6) = rho * thetaeddy
      
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
            real(kind=RP),             intent(inout)       :: dwall
            real(kind=RP),             intent(in)       :: time
            real(kind=RP),             intent(inout)    :: S(NCONS)

            ! --------- my ----------
            real(kind=RP)             :: Q(NCONS)
            real(kind=RP) :: cv2, sbar, somega, S51, S52, S531, S532, S533, S534, S535, S536 ,S541, S542,S55,S56, S57, cwone, theta, eta
            theta = 0.0_RP
            cv2 = 0.7_RP
            cwone = 0.1355_RP / (0.41_RP*0.41_RP) + (1.0_RP + 0.622_RP)/ (2.0_RP/3.0_RP)!
            dwall  = xvec(2) + 1.0_RP    
            call ManufacturedSolutionState( xvec, Q, theta)

            associate ( x => xvec(1), &
                        y => xvec(2), &
                        z => xvec(3), &
                        gamma => thermodynamics % gamma, &
                        Mach => dimensionless % Mach, &
                        Re => dimensionless % Re, &
                        Pr => dimensionless % Pr, &
                        Prt => dimensionless % Prt, &
                        gammaM2 => dimensionless % gammaM2, &
                        tRatio => S_div_Tref_Sutherland )
!           -------------

            somega = 2.0359847367938517*Sqrt(Cos(Pi*(x + y + z))**2)
            
            sbar   =           (5.948839976204641*Cos(Pi*(x + y + z))*  &
        (1 - (Cos(Pi*(x + y + z))*  &
             (tRatio + (-1 + gamma)*gammaM2*  &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +   &
                  0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/  &
           ((-1 + gamma)*gammaM2*(1 + tRatio)*  &
             (1 + (Cos(Pi*(x + y + z))**4*  &
                  (tRatio + (-1 + gamma)*gammaM2*  &
                      (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +   &
                        0.1*Sin(Pi*(x + y + z))))**4*  &
                  (2. + 0.1*Sin(Pi*(x + y + z)))**4)/  &
                ((-1 + gamma)**6*gammaM2**6*(1 + tRatio)**4*  &
                  (357.91099999999994 +   &
                    (Cos(Pi*(x + y + z))**3*  &
                       (tRatio +   &
                          (-1 + gamma)*gammaM2*  &
                           (2. +   &
                             (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +   &
                             0.1*Sin(Pi*(x + y + z))))**3*  &
                       (2. + 0.1*Sin(Pi*(x + y + z)))**3)/  &
                     ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3*  &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +   &
                          0.1*Sin(Pi*(x + y + z)))**3*  &
                       ((-1 + gamma)*gammaM2*  &
                          (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +   &
                            0.1*Sin(Pi*(x + y + z))))**1.5))*  &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +   &
                     0.1*Sin(Pi*(x + y + z)))**6))*  &
             (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)))*  &
             Sqrt((-1 + gamma)*gammaM2*  &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -   &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))  &
                 )))))/dwall**2             

                     

            S(1)=           0. - 0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))      





      if (theta .GE. 0.0_RP) then

        S(2)=                 -0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))**2 -  &
       (-1. + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             1.2566370614359172*Cos(Pi*(x + y + z))* &
              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))))) &
         *(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*(-1. + gamma)*Cos(Pi*(x + y + z))* &
        (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) -  &
       (0.9869604401089358*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)))/Re +  &
       (0.3141592653589793*Cos(Pi*(x + y + z))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2))/Re +  &
       (-1.9739208802178716*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) +  &
          0.6283185307179586*Cos(Pi*(x + y + z))* &
           (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
              (2.*(tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )**2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) -  &
             (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) +  &
             (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**2.* &
                (-0.5*(Cos(Pi*(x + y + z))* &
                      (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                        0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                            (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           1.2566370614359172*Cos(Pi*(x + y + z))* &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           0.6283185307179586*Cos(Pi*(x + y + z))* &
                            (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      ((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))**2* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) -  &
                  (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                     Sin(Pi*(x + y + z)))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))))/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) -  &
             (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**5.* &
                (-0.5*(Cos(Pi*(x + y + z))* &
                      (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                        0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                            (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           1.2566370614359172*Cos(Pi*(x + y + z))* &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           0.6283185307179586*Cos(Pi*(x + y + z))* &
                            (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      ((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))**2* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) -  &
                  (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                     Sin(Pi*(x + y + z)))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))))/ &
              (357.91099999999994 +  &
                 ((Cos(Pi*(x + y + z))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((-1. + gamma)*gammaM2*(1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      Sqrt((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))))**3.)**2))/Re             
!

!
!

       S(3) =                -0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       (-1. + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             1.2566370614359172*Cos(Pi*(x + y + z))* &
              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))))) &
         *(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.6283185307179586*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.6283185307179586*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*(-1. + gamma)*Cos(Pi*(x + y + z))* &
        (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) -  &
       (1.9739208802178716*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)))/Re +  &
       (0.6283185307179586*Cos(Pi*(x + y + z))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2))/Re +  &
       (-3.947841760435743*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) +  &
          1.2566370614359172*Cos(Pi*(x + y + z))* &
           (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
              (2.*(tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )**2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) -  &
             (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) +  &
             (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**2.* &
                (-0.5*(Cos(Pi*(x + y + z))* &
                      (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                        0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                            (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           1.2566370614359172*Cos(Pi*(x + y + z))* &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           0.6283185307179586*Cos(Pi*(x + y + z))* &
                            (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      ((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))**2* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) -  &
                  (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                     Sin(Pi*(x + y + z)))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))))/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) -  &
             (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**5.* &
                (-0.5*(Cos(Pi*(x + y + z))* &
                      (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                        0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                            (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           1.2566370614359172*Cos(Pi*(x + y + z))* &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           0.6283185307179586*Cos(Pi*(x + y + z))* &
                            (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      ((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))**2* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) -  &
                  (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                     Sin(Pi*(x + y + z)))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))))/ &
              (357.91099999999994 +  &
                 ((Cos(Pi*(x + y + z))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((-1. + gamma)*gammaM2*(1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      Sqrt((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))))**3.)**2))/Re 




!!!
!

       S(4) =                     -0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       (-1. + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             1.2566370614359172*Cos(Pi*(x + y + z))* &
              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))))) &
         *(2. + 0.1*Sin(Pi*(x + y + z))) +  &
       0.9424777960769379*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
       0.9424777960769379*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
       0.9424777960769379*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*(-1. + gamma)*Cos(Pi*(x + y + z))* &
        (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) +  &
       (2.9608813203268074*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)))/Re -  &
       (0.9424777960769379*Cos(Pi*(x + y + z))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2))/Re +  &
       (5.921762640653615*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) -  &
          1.8849555921538759*Cos(Pi*(x + y + z))* &
           (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
              (2.*(tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )**2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) -  &
             (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) +  &
             (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**2.* &
                (-0.5*(Cos(Pi*(x + y + z))* &
                      (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                        0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                            (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           1.2566370614359172*Cos(Pi*(x + y + z))* &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           0.6283185307179586*Cos(Pi*(x + y + z))* &
                            (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      ((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))**2* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) -  &
                  (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                     Sin(Pi*(x + y + z)))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))))/ &
              (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.) -  &
             (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**5.* &
                (-0.5*(Cos(Pi*(x + y + z))* &
                      (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                        0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                            (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           1.2566370614359172*Cos(Pi*(x + y + z))* &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                           0.6283185307179586*Cos(Pi*(x + y + z))* &
                            (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      ((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))**2* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) +  &
                  (Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))) -  &
                  (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                     Sin(Pi*(x + y + z)))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))))/ &
              (357.91099999999994 +  &
                 ((Cos(Pi*(x + y + z))* &
                      (tRatio + (-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))) &
                     /((-1. + gamma)*gammaM2*(1 + tRatio)* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))* &
                      Sqrt((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))))**3.)**2))/Re       
!


        S(5) =         0. - (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (0.6283185307179586*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)*(2. + 0.1*Sin(Pi*(x + y + z))) &
           + 0.3141592653589793*(-1 + gamma)*Cos(Pi*(x + y + z))* &
           (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (0.6283185307179586*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)*(2. + 0.1*Sin(Pi*(x + y + z))) &
           + 0.3141592653589793*(-1 + gamma)*Cos(Pi*(x + y + z))* &
           (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       (0.6283185307179586*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)*(2. + 0.1*Sin(Pi*(x + y + z))) &
           + 0.3141592653589793*(-1 + gamma)*Cos(Pi*(x + y + z))* &
           (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))))* &
        (1. + 0.1*Sin(Pi*(x + y + z)))  






       S(5) =     S(5) +               (0.3141592653589793*Cos(Pi*(x + y + z))* &
          (0. + 0.6283185307179586*Cos(Pi*(x + y + z))* &
             (((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.))) +  &
         1.184352528130723*Cos(Pi*(x + y + z))**2* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)) +  &
         1.9739208802178716*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) -  &
         2.9608813203268074*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) -  &
         0.6283185307179586*Cos(Pi*(x + y + z))* &
          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2) +  &
         0.9424777960769379*Cos(Pi*(x + y + z))* &
          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2) +  &
         (1. + 0.1*Sin(Pi*(x + y + z)))* &
          (-1.9739208802178716*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) +  &
            0.6283185307179586*Cos(Pi*(x + y + z))* &
             (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )/ &
                (2.*(tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))* &
                  Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))) +  &
               ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))**2 +  &
               (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) -  &
               (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) +  &
               (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**2.* &
                  (-0.5*(Cos(Pi*(x + y + z))* &
                        (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             1.2566370614359172*Cos(Pi*(x + y + z))* &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             0.6283185307179586*Cos(Pi*(x + y + z))* &
                              (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        ((-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))**2* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z)))))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                     ((1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) -  &
                    (Pi*(tRatio +  &
                         (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z)))))))/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) -  &
               (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**5.* &
                  (-0.5*(Cos(Pi*(x + y + z))* &
                        (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             1.2566370614359172*Cos(Pi*(x + y + z))* &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             0.6283185307179586*Cos(Pi*(x + y + z))* &
                              (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        ((-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))**2* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z)))))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                     ((1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) -  &
                    (Pi*(tRatio +  &
                         (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z)))))))/ &
                (357.91099999999994 +  &
                   ((Cos(Pi*(x + y + z))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        Sqrt((-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))))**3.)**2)))/Re







                S(5) = S(5) +            ((-1 + gamma)*gammaM2*((Cos(Pi*(x + y + z))**4* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (gammaM2*(1 + tRatio)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
                *Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))))* &
          (-0.9869604401089358*Sin(Pi*(x + y + z)) +  &
            (-2.76348923230502*Cos(Pi*(x + y + z))**2 -  &
               5.921762640653615*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                Sin(Pi*(x + y + z)) +  &
               3.947841760435743*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                Sin(Pi*(x + y + z)) +  &
               1.9739208802178716*(1. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) &
               )/2.))/Re + ((-1 + gamma)*gammaM2* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
            (1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               )/2.)*((-3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /(2.*(-1 + gamma)**3*gammaM2**2*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2.5) -  &
            (3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**4*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (1.2566370614359172*Cos(Pi*(x + y + z))**5* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**3) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**2*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**3*gammaM2**2*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
               )/ &
             (2.*Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))* &
               Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z))))) +  &
            (gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) -  &
            ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
                *Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**2) -  &
            (4*Pi*Cos(Pi*(x + y + z))**3* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4* &
               Sin(Pi*(x + y + z)))/ &
             ((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) -  &
            (Cos(Pi*(x + y + z))**4* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4* &
               ((-3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  (2.*(-1 + gamma)**2*gammaM2**2*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**2.5) -  &
                 (3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**4* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) +  &
                 (0.9424777960769379*Cos(Pi*(x + y + z))**4* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**2)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) +  &
                 (3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**2* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (3*Pi*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3*Sin(Pi*(x + y + z)))/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5)))/ &
             ((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                  (Cos(Pi*(x + y + z))**3* &
                     (tRatio + (-1 + gamma)*gammaM2* &
                         (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                           0.1*Sin(Pi*(x + y + z))))**3* &
                     (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                   ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                        0.1*Sin(Pi*(x + y + z)))**3* &
                     ((-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5))**2* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5)))/Re     
           


                S(5) = S(5) +           (0.6283185307179586*Cos(Pi*(x + y + z))* &
          (0. + 1.2566370614359172*Cos(Pi*(x + y + z))* &
             (((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.))) +  &
         0.5921762640653615*Cos(Pi*(x + y + z))**2* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)) +  &
         0.9869604401089358*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) -  &
         2.9608813203268074*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))* &
          (1. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) -  &
         0.3141592653589793*Cos(Pi*(x + y + z))* &
          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2) +  &
         0.9424777960769379*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2) +  &
         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (-3.947841760435743*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) +  &
            1.2566370614359172*Cos(Pi*(x + y + z))* &
             (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )/ &
                (2.*(tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))* &
                  Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))) +  &
               ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))**2 +  &
               (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) -  &
               (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) +  &
               (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**2.* &
                  (-0.5*(Cos(Pi*(x + y + z))* &
                        (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             1.2566370614359172*Cos(Pi*(x + y + z))* &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             0.6283185307179586*Cos(Pi*(x + y + z))* &
                              (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        ((-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))**2* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z)))))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                     ((1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) -  &
                    (Pi*(tRatio +  &
                         (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z)))))))/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) -  &
               (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**5.* &
                  (-0.5*(Cos(Pi*(x + y + z))* &
                        (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             1.2566370614359172*Cos(Pi*(x + y + z))* &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             0.6283185307179586*Cos(Pi*(x + y + z))* &
                              (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        ((-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))**2* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z)))))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                     ((1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) -  &
                    (Pi*(tRatio +  &
                         (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z)))))))/ &
                (357.91099999999994 +  &
                   ((Cos(Pi*(x + y + z))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        Sqrt((-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))))**3.)**2)))/Re 






                S(5) = S(5) +         ((-1 + gamma)*gammaM2*((Cos(Pi*(x + y + z))**4* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (gammaM2*(1 + tRatio)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
                *Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))))* &
          (-0.9869604401089358*Sin(Pi*(x + y + z)) +  &
            (-2.76348923230502*Cos(Pi*(x + y + z))**2 -  &
               5.921762640653615*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                Sin(Pi*(x + y + z)) +  &
               3.947841760435743*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                Sin(Pi*(x + y + z)) +  &
               1.9739208802178716*(1. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) &
               )/2.))/Re + ((-1 + gamma)*gammaM2* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
            (1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               )/2.)*((-3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /(2.*(-1 + gamma)**3*gammaM2**2*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2.5) -  &
            (3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**4*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (1.2566370614359172*Cos(Pi*(x + y + z))**5* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**3) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**2*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**3*gammaM2**2*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
               )/ &
             (2.*Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))* &
               Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z))))) +  &
            (gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) -  &
            ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
                *Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**2) -  &
            (4*Pi*Cos(Pi*(x + y + z))**3* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4* &
               Sin(Pi*(x + y + z)))/ &
             ((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) -  &
            (Cos(Pi*(x + y + z))**4* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4* &
               ((-3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  (2.*(-1 + gamma)**2*gammaM2**2*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**2.5) -  &
                 (3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**4* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) +  &
                 (0.9424777960769379*Cos(Pi*(x + y + z))**4* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**2)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) +  &
                 (3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**2* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (3*Pi*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3*Sin(Pi*(x + y + z)))/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5)))/ &
             ((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                  (Cos(Pi*(x + y + z))**3* &
                     (tRatio + (-1 + gamma)*gammaM2* &
                         (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                           0.1*Sin(Pi*(x + y + z))))**3* &
                     (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                   ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                        0.1*Sin(Pi*(x + y + z)))**3* &
                     ((-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5))**2* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5)))/Re  


            S(5) = S(5) +             (-0.9424777960769379*Cos(Pi*(x + y + z))* &
          (0. - 1.8849555921538759*Cos(Pi*(x + y + z))* &
             (((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.))) -  &
         0.3947841760435743*Cos(Pi*(x + y + z))**2* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)) +  &
         0.9869604401089358*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) +  &
         1.9739208802178716*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.))* &
          (1. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) -  &
         0.3141592653589793*Cos(Pi*(x + y + z))* &
          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2) -  &
         0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2 + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.)/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) +  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**2.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**3.) -  &
            (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
               ((Cos(Pi*(x + y + z))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))))**5.* &
               (-0.5*(Cos(Pi*(x + y + z))* &
                     (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                       0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          1.2566370614359172*Cos(Pi*(x + y + z))* &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                          0.6283185307179586*Cos(Pi*(x + y + z))* &
                           (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((1 + tRatio)*(2. -  &
                       0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     ((-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))**2* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) +  &
                 (Cos(Pi*(x + y + z))* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                      0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                    (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                  ((1 + tRatio)*(2. -  &
                      0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))) -  &
                 (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                    Sin(Pi*(x + y + z)))/ &
                  ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))* &
                    Sqrt((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z)))))))/ &
             (357.91099999999994 +  &
                ((Cos(Pi*(x + y + z))* &
                     (tRatio + (-1. + gamma)*gammaM2* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))* &
                     Sqrt((-1. + gamma)*gammaM2* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z))))))**3.)**2) +  &
         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
          (5.921762640653615*(((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) + (Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.))*Sin(Pi*(x + y + z)) -  &
            1.8849555921538759*Cos(Pi*(x + y + z))* &
             (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )/ &
                (2.*(tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))* &
                  Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))) +  &
               ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )) - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z)))))/ &
                (tRatio + (-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))**2 +  &
               (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) -  &
               (Pi*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.*Sin(Pi*(x + y + z)))/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) +  &
               (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**2.* &
                  (-0.5*(Cos(Pi*(x + y + z))* &
                        (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             1.2566370614359172*Cos(Pi*(x + y + z))* &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             0.6283185307179586*Cos(Pi*(x + y + z))* &
                              (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        ((-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))**2* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z)))))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                     ((1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) -  &
                    (Pi*(tRatio +  &
                         (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z)))))))/ &
                (357.91099999999994 +  &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**3.) -  &
               (3.*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  ((Cos(Pi*(x + y + z))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))))**5.* &
                  (-0.5*(Cos(Pi*(x + y + z))* &
                        (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             1.2566370614359172*Cos(Pi*(x + y + z))* &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                             0.6283185307179586*Cos(Pi*(x + y + z))* &
                              (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        ((-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                       )/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))**2* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                       (tRatio + (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z)))))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) +  &
                    (Cos(Pi*(x + y + z))* &
                       (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                         0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                             (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            1.2566370614359172*Cos(Pi*(x + y + z))* &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                            0.6283185307179586*Cos(Pi*(x + y + z))* &
                             (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                     ((1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z))))) -  &
                    (Pi*(tRatio +  &
                         (-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))* &
                       (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
                     ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                       (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                         0.1*Sin(Pi*(x + y + z)))* &
                       Sqrt((-1. + gamma)*gammaM2* &
                         (2. - 0.5* &
                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                           0.1*Sin(Pi*(x + y + z)))))))/ &
                (357.91099999999994 +  &
                   ((Cos(Pi*(x + y + z))* &
                        (tRatio +  &
                          (-1. + gamma)*gammaM2* &
                           (2. - 0.5* &
                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                             0.1*Sin(Pi*(x + y + z))))* &
                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                      ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                        (2. - 0.5* &
                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                          0.1*Sin(Pi*(x + y + z)))* &
                        Sqrt((-1. + gamma)*gammaM2* &
                          (2. - 0.5* &
                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                            0.1*Sin(Pi*(x + y + z))))))**3.)**2)))/Re  


   
                S(5) = S(5) +         ((-1 + gamma)*gammaM2*((Cos(Pi*(x + y + z))**4* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (gammaM2*(1 + tRatio)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
                *Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))))* &
          (-0.9869604401089358*Sin(Pi*(x + y + z)) +  &
            (-2.76348923230502*Cos(Pi*(x + y + z))**2 -  &
               5.921762640653615*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                Sin(Pi*(x + y + z)) +  &
               3.947841760435743*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                Sin(Pi*(x + y + z)) +  &
               1.9739208802178716*(1. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) &
               )/2.))/Re + ((-1 + gamma)*gammaM2* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
            (1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               )/2.)*((-3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /(2.*(-1 + gamma)**3*gammaM2**2*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2.5) -  &
            (3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**4*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (1.2566370614359172*Cos(Pi*(x + y + z))**5* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**3) &
              /((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            (3*Cos(Pi*(x + y + z))**4* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**2*(2. + 0.1*Sin(Pi*(x + y + z)))**4) &
              /((-1 + gamma)**3*gammaM2**2*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) +  &
            ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
               )/ &
             (2.*Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))* &
               Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z))))) +  &
            (gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                 (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) -  &
            ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                 (1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))) &
                *Sqrt((-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))/ &
             (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**2) -  &
            (4*Pi*Cos(Pi*(x + y + z))**3* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4* &
               Sin(Pi*(x + y + z)))/ &
             ((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                 (Cos(Pi*(x + y + z))**3* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5))* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5) -  &
            (Cos(Pi*(x + y + z))**4* &
               (tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))**3*(2. + 0.1*Sin(Pi*(x + y + z)))**4* &
               ((-3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  (2.*(-1 + gamma)**2*gammaM2**2*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**2.5) -  &
                 (3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**4* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) +  &
                 (0.9424777960769379*Cos(Pi*(x + y + z))**4* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**2)/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) +  &
                 (3*Cos(Pi*(x + y + z))**3* &
                    (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                      (1.8849555921538759*Cos(Pi*(x + y + z))* &
                          (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         1.2566370614359172*Cos(Pi*(x + y + z))* &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                         0.6283185307179586*Cos(Pi*(x + y + z))* &
                          (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**2* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                  ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5) -  &
                 (3*Pi*Cos(Pi*(x + y + z))**2* &
                    (tRatio + (-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**3* &
                    (2. + 0.1*Sin(Pi*(x + y + z)))**3*Sin(Pi*(x + y + z)))/ &
                  ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))**3* &
                    ((-1 + gamma)*gammaM2* &
                       (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                            (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                         0.1*Sin(Pi*(x + y + z))))**1.5)))/ &
             ((-1 + gamma)**4*gammaM2**3*Mach**2*Prt*(1 + tRatio)**3* &
               (357.91099999999994 +  &
                  (Cos(Pi*(x + y + z))**3* &
                     (tRatio + (-1 + gamma)*gammaM2* &
                         (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                              (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                           0.1*Sin(Pi*(x + y + z))))**3* &
                     (2. + 0.1*Sin(Pi*(x + y + z)))**3)/ &
                   ((-1 + gamma)**3*gammaM2**3*(1 + tRatio)**3* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                        0.1*Sin(Pi*(x + y + z)))**3* &
                     ((-1 + gamma)*gammaM2* &
                        (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                             (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                          0.1*Sin(Pi*(x + y + z))))**1.5))**2* &
               (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )**3*((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**1.5)))/Re  


! CORRECT FLUXES THETA > 0

                     S(6) =            0. - 0.3141592653589793*Cos(Pi*(x + y + z))**2* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))**2* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))**2*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       (44.41321980490211*(-1. + gamma)*gammaM2*(1 + tRatio)*Cos(Pi*(x + y + z))* &
          (1. + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ))*(2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) +  &
       Pi*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
        Sin(Pi*(x + y + z)) + Pi*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) +  &
       Pi*(1. + 0.1*Sin(Pi*(x + y + z)))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
        Sin(Pi*(x + y + z)) - (7.0685834705770345*(-1. + gamma)**2*gammaM2**2* &
          (1 + tRatio)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
            0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               ))*(1. + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ))*(2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sin(Pi*(x + y + z)))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) -  &
       (14.137166941154069*(-1. + gamma)*gammaM2*(1 + tRatio)* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
            0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               ))*(1. + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ))*Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sin(Pi*(x + y + z)))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) +  &
       (14.137166941154069*(-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
            0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               ))*(1. + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ))*(2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sin(Pi*(x + y + z)))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
              (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                   (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                   (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))**2 &
          ) - (14.137166941154069*(-1. + gamma)*gammaM2*(1 + tRatio)* &
          (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sin(Pi*(x + y + z))*(-0.5* &
             (Cos(Pi*(x + y + z))* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
              ((1 + tRatio)*(2. -  &
                  0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                ((-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))**1.5) - (Cos(Pi*(x + y + z))* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))** &
                2*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 ))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + (Cos(Pi*(x + y + z))* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
               )/ &
             ((1 + tRatio)*(2. - 0.5* &
                  ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) - (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &!
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))))       
  







!CORRECT SOURCE TERMS THETA >= 0
        if (sbar .GE. -cv2 * somega * Re  ) then 

                    S(6) = S(6) +         (-2.0051747451504216*cwone*Cos(Pi*(x + y + z))**2.*  &
                       (Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/  &
                           (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                                 0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                               (5.948839976204641*Cos(Pi*(x + y + z))*  &
                                  (1. - (Cos(Pi*(x + y + z))*  &
                                       (tRatio +   &
                                         (-1. + gamma)*gammaM2*  &
                                          (2. -   &
                                            0.5*  &
                                             ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                            0.1*Sin(Pi*(x + y + z))))*  &
                                       (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                     ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                       (1. + (Cos(Pi*(x + y + z))*  &
                                            (tRatio +   &
                                              (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                            (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                            ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                          ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                            (357.91099999999994 +   &
                                              ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                            (2. -   &
                                              0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                              0.1*Sin(Pi*(x + y + z)))*  &
                                            Sqrt((-1. + gamma)*gammaM2*  &
                                              (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))*  &
                                       (2. - 0.5*  &
                                          ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                         0.1*Sin(Pi*(x + y + z)))*  &
                                       Sqrt((-1. + gamma)*gammaM2*  &
                                         (2. -   &
                                           0.5*  &
                                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))) +   &
                         0.3*(-Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/  &
                               (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                                     0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                                   (5.948839976204641*Cos(Pi*(x + y + z))*  &
                                      (1. - (Cos(Pi*(x + y + z))*  &
                                           (tRatio +   &
                                             (-1. + gamma)*gammaM2*  &
                                              (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                         ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                           (1. +   &
                                             (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                              ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (357.91099999999994 +   &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))*  &
                                           (2. -   &
                                             0.5*  &
                                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                             0.1*Sin(Pi*(x + y + z)))*  &
                                           Sqrt((-1. + gamma)*gammaM2*  &
                                             (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))) +   &
                            Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/  &
                               (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                                     0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                                   (5.948839976204641*Cos(Pi*(x + y + z))*  &
                                      (1. - (Cos(Pi*(x + y + z))*  &
                                           (tRatio +   &
                                             (-1. + gamma)*gammaM2*  &
                                              (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                         ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                           (1. +   &
                                             (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                              ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (357.91099999999994 +   &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))*  &
                                           (2. -   &
                                             0.5*  &
                                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                             0.1*Sin(Pi*(x + y + z)))*  &
                                           Sqrt((-1. + gamma)*gammaM2*  &
                                             (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))))**6.)  &
                         )*(1/  &
                          (64. + (Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/  &
                                 (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                                       0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                                     (5.948839976204641*Cos(Pi*(x + y + z))*  &
                                        (1. - (Cos(Pi*(x + y + z))*  &
                                             (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                             (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                             (1. +   &
                                               (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (357.91099999999994 +   &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))*  &
                                             (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                             Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))) +   &
                               0.3*(-Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/  &
                                     (dwall**2*Re*  &
                                       (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                                         (5.948839976204641*Cos(Pi*(x + y + z))*  &
                                            (1. -   &
                                              (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (1. +   &
                                               (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (357.91099999999994 +   &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))) +   &
                                  Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/  &
                                     (dwall**2*Re*  &
                                       (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                                         (5.948839976204641*Cos(Pi*(x + y + z))*  &
                                            (1. -   &
                                              (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (1. +   &
                                               (Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (357.91099999999994 +   &
                                               ((Cos(Pi*(x + y + z))*  &
                                               (tRatio +   &
                                               (-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                               (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                               ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))*  &
                                               Sqrt((-1. + gamma)*gammaM2*  &
                                               (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))))**6.)  &
                               )**6.))**0.16666666666666666*(2. + 0.1*Sin(Pi*(x + y + z))))/  &
                     (dwall**2.*Re) + 0.1355*Cos(Pi*(x + y + z))*  &
                     (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +   &
                         0.09869604401089357*Cos(Pi*(x + y + z))**2.) +   &
                       (5.948839976204641*Cos(Pi*(x + y + z))*  &
                          (1. - (Cos(Pi*(x + y + z))*  &
                               (tRatio + (-1. + gamma)*gammaM2*  &
                                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                    0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/  &
                             ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                               (1. + (Cos(Pi*(x + y + z))*  &
                                    (tRatio + (-1. + gamma)*gammaM2*  &
                                       (2. - 0.5*  &
                                          ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                            (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                            (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                         0.1*Sin(Pi*(x + y + z))))*  &
                                    (2. + 0.1*Sin(Pi*(x + y + z)))*  &
                                    ((Cos(Pi*(x + y + z))*  &
                                         (tRatio +   &
                                           (-1. + gamma)*gammaM2*  &
                                            (2. -   &
                                              0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                              0.1*Sin(Pi*(x + y + z))))*  &
                                         (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                       ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                         (2. -   &
                                           0.5*  &
                                            ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                              (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                           0.1*Sin(Pi*(x + y + z)))*  &
                                         Sqrt((-1. + gamma)*gammaM2*  &
                                           (2. -   &
                                             0.5*  &
                                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                             0.1*Sin(Pi*(x + y + z))))))**3.)/  &
                                  ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                    (357.91099999999994 +   &
                                      ((Cos(Pi*(x + y + z))*  &
                                           (tRatio +   &
                                             (-1. + gamma)*gammaM2*  &
                                              (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))*  &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/  &
                                         ((-1. + gamma)*gammaM2*(1 + tRatio)*  &
                                           (2. -   &
                                             0.5*  &
                                              ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                             0.1*Sin(Pi*(x + y + z)))*  &
                                           Sqrt((-1. + gamma)*gammaM2*  &
                                             (2. -   &
                                               0.5*  &
                                               ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                               (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                               0.1*Sin(Pi*(x + y + z))))))**3.)*  &
                                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                      0.1*Sin(Pi*(x + y + z)))*  &
                                    Sqrt((-1. + gamma)*gammaM2*  &
                                      (2. - 0.5*  &
                                         ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                        0.1*Sin(Pi*(x + y + z))))))*  &
                               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))  &
                                *Sqrt((-1. + gamma)*gammaM2*  &
                                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +   &
                                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) +   &
                                   0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))*  &
                     (2. + 0.1*Sin(Pi*(x + y + z))) +   &
                    (27.625022718649113*(2. + 0.1*Sin(Pi*(x + y + z)))*  &
                       (-Sin(Pi*(x + y + z)))**2.)/Re  

                    else 

                    S(6) = S(6) +         (-2.0051747451504216*cwone*Cos(Pi*(x + y + z))**2.* &
                   (Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/ &
                       (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                             0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                           (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                              (0.48999999999999994* &
                                 Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                   0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                (5.353955978584177*Cos(Pi*(x + y + z))* &
                                   (1. - (Cos(Pi*(x + y + z))* &
                                        (tRatio +  &
                                          (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                      ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                        (1. +  &
                                          (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                        (2. -  &
                                          0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                          0.1*Sin(Pi*(x + y + z)))* &
                                        Sqrt((-1. + gamma)*gammaM2* &
                                          (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                            (-0.4999999999999999* &
                               Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                 0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                              (5.948839976204641*Cos(Pi*(x + y + z))* &
                                 (1. - (Cos(Pi*(x + y + z))* &
                                      (tRatio +  &
                                        (-1. + gamma)*gammaM2* &
                                         (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                      (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                    ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                      (1. +  &
                                        (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                         ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                      (2. -  &
                                        0.5* &
                                         ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                        0.1*Sin(Pi*(x + y + z)))* &
                                      Sqrt((-1. + gamma)*gammaM2* &
                                        (2. -  &
                                          0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                          0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))))) +  &
                     0.3*(-Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/ &
                           (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                 0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                               (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                    0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                                  (0.48999999999999994* &
                                     Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                       0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                    (5.353955978584177*Cos(Pi*(x + y + z))* &
                                       (1. -  &
                                         (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                          ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                                (-0.4999999999999999* &
                                   Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                     0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                                  (5.948839976204641*Cos(Pi*(x + y + z))* &
                                     (1. -  &
                                       (Cos(Pi*(x + y + z))* &
                                          (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                          (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                        ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                          (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                          (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                          Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))))) +  &
                        Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/ &
                           (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                 0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                               (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                    0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                                  (0.48999999999999994* &
                                     Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                       0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                    (5.353955978584177*Cos(Pi*(x + y + z))* &
                                       (1. -  &
                                         (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                          ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                                (-0.4999999999999999* &
                                   Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                     0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                                  (5.948839976204641*Cos(Pi*(x + y + z))* &
                                     (1. -  &
                                       (Cos(Pi*(x + y + z))* &
                                          (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                          (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                        ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                          (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                          (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                          Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))))**6. &
                        ))*(1/ &
                      (64. + (Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/ &
                             (dwall**2*Re*(Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                   0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                 (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                      0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                                    (0.48999999999999994* &
                                       Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                         0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                      (5.353955978584177*Cos(Pi*(x + y + z))* &
                                         (1. -  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                                  (-0.4999999999999999* &
                                     Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                       0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                                    (5.948839976204641*Cos(Pi*(x + y + z))* &
                                       (1. -  &
                                         (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                          ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))))) +  &
                           0.3*(-Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/ &
                                 (dwall**2*Re* &
                                   (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                       0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                     (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                          0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                                        (0.48999999999999994* &
                                           Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                          (5.353955978584177*Cos(Pi*(x + y + z))* &
                                           (1. -  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                                      (-0.4999999999999999* &
                                         Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                                        (5.948839976204641*Cos(Pi*(x + y + z))* &
                                           (1. -  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re))))) +  &
                              Min(2.,(5.948839976204641*Cos(Pi*(x + y + z)))/ &
                                 (dwall**2*Re* &
                                   (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                       0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                     (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                          0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                                        (0.48999999999999994* &
                                           Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                                          (5.353955978584177*Cos(Pi*(x + y + z))* &
                                           (1. -  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                                      (-0.4999999999999999* &
                                         Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                                        (5.948839976204641*Cos(Pi*(x + y + z))* &
                                           (1. -  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (1. +  &
                                           (Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (357.91099999999994 +  &
                                           ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                           ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))))**6. &
                              ))**6.))**0.16666666666666666*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                 (dwall**2.*Re) + 0.1355*Cos(Pi*(x + y + z))* &
                 (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                     0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                   (Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                        0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
                      (0.48999999999999994* &
                         Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                           0.09869604401089357*Cos(Pi*(x + y + z))**2.) +  &
                        (5.353955978584177*Cos(Pi*(x + y + z))* &
                           (1. - (Cos(Pi*(x + y + z))* &
                                (tRatio + (-1. + gamma)*gammaM2* &
                                   (2. - 0.5* &
                                      ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                     0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
                                )/ &
                              ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                (1. + (Cos(Pi*(x + y + z))* &
                                     (tRatio +  &
                                       (-1. + gamma)*gammaM2* &
                                        (2. -  &
                                          0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                          0.1*Sin(Pi*(x + y + z))))* &
                                     (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                     ((Cos(Pi*(x + y + z))* &
                                          (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                          (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                        ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                          (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                          Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                   ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                     (357.91099999999994 +  &
                                       ((Cos(Pi*(x + y + z))* &
                                           (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                           (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                          ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                           Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                     (2. -  &
                                       0.5* &
                                        ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                       0.1*Sin(Pi*(x + y + z)))* &
                                     Sqrt((-1. + gamma)*gammaM2* &
                                       (2. -  &
                                         0.5* &
                                          ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                         0.1*Sin(Pi*(x + y + z))))))* &
                                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                  0.1*Sin(Pi*(x + y + z)))* &
                                Sqrt((-1. + gamma)*gammaM2* &
                                  (2. - 0.5* &
                                     ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                    0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))/ &
                    (-0.4999999999999999*Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
                         0.09869604401089357*Cos(Pi*(x + y + z))**2.) -  &
                      (5.948839976204641*Cos(Pi*(x + y + z))* &
                         (1. - (Cos(Pi*(x + y + z))* &
                              (tRatio + (-1. + gamma)*gammaM2* &
                                 (2. - 0.5* &
                                    ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                   0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                            ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                              (1. + (Cos(Pi*(x + y + z))* &
                                   (tRatio +  &
                                     (-1. + gamma)*gammaM2* &
                                      (2. -  &
                                        0.5* &
                                         ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                        0.1*Sin(Pi*(x + y + z))))* &
                                   (2. + 0.1*Sin(Pi*(x + y + z)))* &
                                   ((Cos(Pi*(x + y + z))* &
                                        (tRatio +  &
                                          (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                        (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                      ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                        (2. -  &
                                          0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                          0.1*Sin(Pi*(x + y + z)))* &
                                        Sqrt((-1. + gamma)*gammaM2* &
                                          (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)/ &
                                 ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                   (357.91099999999994 +  &
                                     ((Cos(Pi*(x + y + z))* &
                                          (tRatio +  &
                                           (-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))* &
                                          (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                                        ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                                          (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z)))* &
                                          Sqrt((-1. + gamma)*gammaM2* &
                                           (2. -  &
                                           0.5* &
                                           ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                           0.1*Sin(Pi*(x + y + z))))))**3.)* &
                                   (2. - 0.5* &
                                      ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                     0.1*Sin(Pi*(x + y + z)))* &
                                   Sqrt((-1. + gamma)*gammaM2* &
                                     (2. -  &
                                       0.5* &
                                        ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                       0.1*Sin(Pi*(x + y + z))))))* &
                              (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                   (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                   (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                0.1*Sin(Pi*(x + y + z)))* &
                              Sqrt((-1. + gamma)*gammaM2* &
                                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                                  0.1*Sin(Pi*(x + y + z)))))))/(dwall**2*Re)))* &
                 (2. + 0.1*Sin(Pi*(x + y + z))) +  &
                (27.625022718649113*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                   (-Sin(Pi*(x + y + z)))**2.)/Re 



                   end if





      


         else









          S(2)=                (0.3141592653589793*Cos(Pi*(x + y + z))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2))/Re - 0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))**2 -  &
       (-1. + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             1.2566370614359172*Cos(Pi*(x + y + z))* &
              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))))) &
         *(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*(-1. + gamma)*Cos(Pi*(x + y + z))* &
        (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) -  &
       (0.9869604401089358*(0. + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) &
           *Sin(Pi*(x + y + z)))/Re +  &
       (0.6283185307179586*Cos(Pi*(x + y + z))* &
           (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
              (2.*(tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )**2) - 1.9739208802178716* &
           (0. + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             )*Sin(Pi*(x + y + z)))/Re    


           
            S(3) =               -0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
       (0.6283185307179586*Cos(Pi*(x + y + z))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2))/Re - 0.3141592653589793*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       (-1. + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             1.2566370614359172*Cos(Pi*(x + y + z))* &
              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))))) &
         *(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.6283185307179586*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.6283185307179586*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*(-1. + gamma)*Cos(Pi*(x + y + z))* &
        (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) -  &
       (1.9739208802178716*(0. + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) &
           *Sin(Pi*(x + y + z)))/Re +  &
       (1.2566370614359172*Cos(Pi*(x + y + z))* &
           (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
              (2.*(tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )**2) - 3.947841760435743* &
           (0. + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             )*Sin(Pi*(x + y + z)))/Re  



      S(4) =               -0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
       0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       (0.9424777960769379*Cos(Pi*(x + y + z))* &
          (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
             (2.*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               **2))/Re - 0.3141592653589793*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       (-1. + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
          0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             1.2566370614359172*Cos(Pi*(x + y + z))* &
              (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
             0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))))) &
         *(2. + 0.1*Sin(Pi*(x + y + z))) +  &
       0.9424777960769379*Cos(Pi*(x + y + z))* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
       0.9424777960769379*Cos(Pi*(x + y + z))* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
       0.9424777960769379*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) -  &
       0.3141592653589793*(-1. + gamma)*Cos(Pi*(x + y + z))* &
        (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
             (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) +  &
       (2.9608813203268074*(0. + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               )/ &
             (tRatio + (-1. + gamma)*gammaM2* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) &
           *Sin(Pi*(x + y + z)))/Re +  &
       (-1.8849555921538759*Cos(Pi*(x + y + z))* &
           (((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))/ &
              (2.*(tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  )) + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               - ((-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )**2) + 5.921762640653615* &
           (0. + ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                Sqrt((-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                  ))/ &
              (tRatio + (-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
             )*Sin(Pi*(x + y + z)))/Re    





       S(5)=                0. - (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (0.6283185307179586*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)*(2. + 0.1*Sin(Pi*(x + y + z))) &
           + 0.3141592653589793*(-1 + gamma)*Cos(Pi*(x + y + z))* &
           (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (0.6283185307179586*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)*(2. + 0.1*Sin(Pi*(x + y + z))) &
           + 0.3141592653589793*(-1 + gamma)*Cos(Pi*(x + y + z))* &
           (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       (0.6283185307179586*Cos(Pi*(x + y + z))*(2. + 0.1*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)*(2. + 0.1*Sin(Pi*(x + y + z))) &
           + 0.3141592653589793*(-1 + gamma)*Cos(Pi*(x + y + z))* &
           (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z))))* &
        (1. + 0.1*Sin(Pi*(x + y + z))) +  &
       (0.6283185307179586*Cos(Pi*(x + y + z))* &
           (0. + 1.2566370614359172*Cos(Pi*(x + y + z))* &
              (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))))) -  &
          0.3141592653589793*Cos(Pi*(x + y + z))* &
           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*(tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))) -  &
             ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2) +  &
          (-1 + gamma)*gammaM2*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
           (((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             (gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))) -  &
             ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                   (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z))))**2)) +  &
          0.5921762640653615*Cos(Pi*(x + y + z))**2* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z))))) +  &
          0.9424777960769379*Cos(Pi*(x + y + z))* &
           (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*(tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))) -  &
             ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2)*(1. + 0.1*Sin(Pi*(x + y + z))) +  &
          0.9869604401089358*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z)) -  &
          2.9608813203268074*(0. +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
           Sin(Pi*(x + y + z)) + (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (1.2566370614359172*Cos(Pi*(x + y + z))* &
              (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))/ &
                 (2.*(tRatio + (-1 + gamma)*gammaM2* &
                      (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                        0.1*Sin(Pi*(x + y + z))))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z))))) +  &
                ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))) -  &
                ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z))))**2) -  &
             3.947841760435743*(0. +  &
                ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*gammaM2*(0. +  &
             (gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))))* &
           (-0.9869604401089358*Sin(Pi*(x + y + z)) +  &
             (-2.76348923230502*Cos(Pi*(x + y + z))**2 -  &
                5.921762640653615*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                 Sin(Pi*(x + y + z)) +  &
                3.947841760435743*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                 Sin(Pi*(x + y + z)) +  &
                1.9739208802178716*(1. + 0.1*Sin(Pi*(x + y + z)))* &
                 Sin(Pi*(x + y + z)))/2.))/Re +  &
       (0.3141592653589793*Cos(Pi*(x + y + z))* &
           (0. + 0.6283185307179586*Cos(Pi*(x + y + z))* &
              (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))))) -  &
          0.6283185307179586*Cos(Pi*(x + y + z))* &
           (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*(tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))) -  &
             ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2) +  &
          0.9424777960769379*Cos(Pi*(x + y + z))* &
           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*(tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))) -  &
             ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2) +  &
          (-1 + gamma)*gammaM2*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
           (((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             (gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))) -  &
             ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                   (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z))))**2)) +  &
          1.184352528130723*Cos(Pi*(x + y + z))**2* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z))))) +  &
          1.9739208802178716*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z)) -  &
          2.9608813203268074*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z)) +  &
          (1. + 0.1*Sin(Pi*(x + y + z)))* &
           (0.6283185307179586*Cos(Pi*(x + y + z))* &
              (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))/ &
                 (2.*(tRatio + (-1 + gamma)*gammaM2* &
                      (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                        0.1*Sin(Pi*(x + y + z))))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z))))) +  &
                ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))) -  &
                ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z))))**2) -  &
             1.9739208802178716*(0. +  &
                ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*gammaM2*(0. +  &
             (gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))))* &
           (-0.9869604401089358*Sin(Pi*(x + y + z)) +  &
             (-2.76348923230502*Cos(Pi*(x + y + z))**2 -  &
                5.921762640653615*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                 Sin(Pi*(x + y + z)) +  &
                3.947841760435743*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                 Sin(Pi*(x + y + z)) +  &
                1.9739208802178716*(1. + 0.1*Sin(Pi*(x + y + z)))* &
                 Sin(Pi*(x + y + z)))/2.))/Re +  &
       (-0.9424777960769379*Cos(Pi*(x + y + z))* &
           (0. - 1.8849555921538759*Cos(Pi*(x + y + z))* &
              (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))))) -  &
          0.3141592653589793*Cos(Pi*(x + y + z))* &
           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*(tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))) -  &
             ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2) +  &
          (-1 + gamma)*gammaM2*(0.3141592653589793*Cos(Pi*(x + y + z)) +  &
             (1.8849555921538759*Cos(Pi*(x + y + z))* &
                 (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                1.2566370614359172*Cos(Pi*(x + y + z))* &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                0.6283185307179586*Cos(Pi*(x + y + z))* &
                 (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
           (((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             (gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))) -  &
             ((-1 + gamma)*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                   (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z))))**2)) -  &
          0.3947841760435743*Cos(Pi*(x + y + z))**2* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z))))) -  &
          0.6283185307179586*Cos(Pi*(x + y + z))* &
           (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  ))/ &
              (2.*(tRatio + (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))) +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))) -  &
             ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                  (1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z))))**2)*(1. + 0.1*Sin(Pi*(x + y + z))) +  &
          0.9869604401089358*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (0. + ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z)) +  &
          1.9739208802178716*(0. +  &
             ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (tRatio + (-1 + gamma)*gammaM2* &
                 (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                   0.1*Sin(Pi*(x + y + z)))))*(1. + 0.1*Sin(Pi*(x + y + z)))* &
           Sin(Pi*(x + y + z)) + (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
           (-1.8849555921538759*Cos(Pi*(x + y + z))* &
              (((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))/ &
                 (2.*(tRatio + (-1 + gamma)*gammaM2* &
                      (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                        0.1*Sin(Pi*(x + y + z))))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z))))) +  &
                ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))) -  &
                ((-1 + gamma)**2*gammaM2**2*(1 + tRatio)* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) +  &
                     (1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z))))/2.)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z))))**2) +  &
             5.921762640653615*(0. +  &
                ((-1 + gamma)*gammaM2*(1 + tRatio)* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z)))* &
                   Sqrt((-1 + gamma)*gammaM2* &
                     (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                 (tRatio + (-1 + gamma)*gammaM2* &
                    (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                      0.1*Sin(Pi*(x + y + z)))))*Sin(Pi*(x + y + z))) +  &
          (-1 + gamma)*gammaM2*(0. +  &
             (gammaM2*(1 + tRatio)* &
                (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. + 0.1*Sin(Pi*(x + y + z)) &
                  )*Sqrt((-1 + gamma)*gammaM2* &
                  (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                    0.1*Sin(Pi*(x + y + z)))))/ &
              (Mach**2*Pr*(tRatio +  &
                  (-1 + gamma)*gammaM2* &
                   (2. + (-(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 -  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2)/2. +  &
                     0.1*Sin(Pi*(x + y + z))))))* &
           (-0.9869604401089358*Sin(Pi*(x + y + z)) +  &
             (-2.76348923230502*Cos(Pi*(x + y + z))**2 -  &
                5.921762640653615*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                 Sin(Pi*(x + y + z)) +  &
                3.947841760435743*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
                 Sin(Pi*(x + y + z)) +  &
                1.9739208802178716*(1. + 0.1*Sin(Pi*(x + y + z)))* &
                 Sin(Pi*(x + y + z)))/2.))/Re                        



!   THEETA < = 0 CORRECT FLUXES 

                   S(6) =          0. - 0.3141592653589793*Cos(Pi*(x + y + z))**2* &
        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))**2* &
        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) -  &
       0.3141592653589793*Cos(Pi*(x + y + z))**2*(1. + 0.1*Sin(Pi*(x + y + z))) -  &
       (44.41321980490211*(-1. + gamma)*gammaM2*(1 + tRatio)*Cos(Pi*(x + y + z))* &
          (1 + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + 0.5*((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**2.)* &
          (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) +  &
       Pi*(5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
        Sin(Pi*(x + y + z)) + Pi*(-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))* &
        (2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)) +  &
       Pi*(1. + 0.1*Sin(Pi*(x + y + z)))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
        Sin(Pi*(x + y + z)) - (7.0685834705770345*(-1. + gamma)**2*gammaM2**2* &
          (1 + tRatio)*(0.3141592653589793*Cos(Pi*(x + y + z)) -  &
            0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               ))*(1 + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + 0.5*((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**2.)* &
          (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sin(Pi*(x + y + z)))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) -  &
       (14.137166941154069*(-1. + gamma)*gammaM2*(1 + tRatio)* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
            0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               ))*(1 + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + 0.5*((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**2.)* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sin(Pi*(x + y + z)))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))) +  &
       (14.137166941154069*(-1. + gamma)**2*gammaM2**2*(1 + tRatio)* &
          (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
            0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               1.2566370614359172*Cos(Pi*(x + y + z))* &
                (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
               0.6283185307179586*Cos(Pi*(x + y + z))*(1. + 0.1*Sin(Pi*(x + y + z))) &
               ))*(1 + (Cos(Pi*(x + y + z))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + 0.5*((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**2.)* &
          (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sin(Pi*(x + y + z)))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
              (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                   (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                   (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))**2 &
          ) - (14.137166941154069*(-1. + gamma)*gammaM2*(1 + tRatio)* &
          (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
               (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
          Sqrt((-1. + gamma)*gammaM2* &
            (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                 (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))))* &
          Sin(Pi*(x + y + z))*(-0.5* &
             (Cos(Pi*(x + y + z))* &
                (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                  0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                      (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     1.2566370614359172*Cos(Pi*(x + y + z))* &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                     0.6283185307179586*Cos(Pi*(x + y + z))* &
                      (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                (tRatio + (-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
              ((1 + tRatio)*(2. -  &
                  0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
                ((-1. + gamma)*gammaM2* &
                   (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     ))**1.5) - (Cos(Pi*(x + y + z))* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                     (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))** &
                2*Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
               (tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 ))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + (Cos(Pi*(x + y + z))* &
               (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                 0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                     (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    1.2566370614359172*Cos(Pi*(x + y + z))* &
                     (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                    0.6283185307179586*Cos(Pi*(x + y + z))* &
                     (1. + 0.1*Sin(Pi*(x + y + z)))))*(2. + 0.1*Sin(Pi*(x + y + z))) &
               )/ &
             ((1 + tRatio)*(2. - 0.5* &
                  ((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) - (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                 )*(2. + 0.1*Sin(Pi*(x + y + z)))*Sin(Pi*(x + y + z)))/ &
             ((-1. + gamma)*gammaM2*(1 + tRatio)* &
               (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                    (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))* &
               Sqrt((-1. + gamma)*gammaM2* &
                 (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                      (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))) &
               ) + 1.*((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**1.* &
             (-0.5*(Cos(Pi*(x + y + z))* &
                   (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                     0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                         (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                        1.2566370614359172*Cos(Pi*(x + y + z))* &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                        0.6283185307179586*Cos(Pi*(x + y + z))* &
                         (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                   (tRatio + (-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                 ((1 + tRatio)*(2. -  &
                     0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )*((-1. + gamma)*gammaM2* &
                      (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                           (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                        0.1*Sin(Pi*(x + y + z))))**1.5) -  &
               (Cos(Pi*(x + y + z))* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                        (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)) &
                     )**2*Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))) +  &
               (0.3141592653589793*Cos(Pi*(x + y + z))**2* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z)))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))) +  &
               (Cos(Pi*(x + y + z))* &
                  (0.3141592653589793*Cos(Pi*(x + y + z)) -  &
                    0.5*(-1.8849555921538759*Cos(Pi*(x + y + z))* &
                        (5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       1.2566370614359172*Cos(Pi*(x + y + z))* &
                        (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z)))) +  &
                       0.6283185307179586*Cos(Pi*(x + y + z))* &
                        (1. + 0.1*Sin(Pi*(x + y + z)))))* &
                  (2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((1 + tRatio)*(2. -  &
                    0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))) -  &
               (Pi*(tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z)))* &
                  Sin(Pi*(x + y + z)))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))))/ &
        (Re*(tRatio + (-1. + gamma)*gammaM2* &
             (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                  (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z)))))      
                  

!CORRECT SOURCE TERMS THETA < 0

                S(6) = S(6) +         (cwone*Cos(Pi*(x + y + z))**2*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
        (dwall**2*Re) + 0.1355*Cos(Pi*(x + y + z))* &
        Sqrt(4.046537804446636*Cos(Pi*(x + y + z))**2 +  &
          0.09869604401089357*Cos(Pi*(x + y + z))**2.)* &
        (1 - (1000.*((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**2.)/ &
           (1 + ((Cos(Pi*(x + y + z))* &
                  (tRatio + (-1. + gamma)*gammaM2* &
                     (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                          (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                       0.1*Sin(Pi*(x + y + z))))*(2. + 0.1*Sin(Pi*(x + y + z))))/ &
                ((-1. + gamma)*gammaM2*(1 + tRatio)* &
                  (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                       (1. + 0.1*Sin(Pi*(x + y + z)))**2) + 0.1*Sin(Pi*(x + y + z))) &
                   *Sqrt((-1. + gamma)*gammaM2* &
                    (2. - 0.5*((5. - 3.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (-4. + 2.*(2. + 0.1*Sin(Pi*(x + y + z))))**2 +  &
                         (1. + 0.1*Sin(Pi*(x + y + z)))**2) +  &
                      0.1*Sin(Pi*(x + y + z))))))**2.))* &
        (2. + 0.1*Sin(Pi*(x + y + z))) +  &
       (27.625022718649113*(2. + 0.1*Sin(Pi*(x + y + z)))* &
          (-Sin(Pi*(x + y + z)))**2.)/Re 

                  end if 



      end associate
   END SUBROUTINE ManufacturedSolutionSourceNSSA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   

END MODULE ManufacturedSolutionsNSSA
