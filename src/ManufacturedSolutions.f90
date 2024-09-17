!
!////////////////////////////////////////////////////////////////////////
!
!      Manufactured solutions definitions for 3D case (Navier-Stokes/Spalart-Almaras)
!
!////////////////////////////////////////////////////////////////////////////////////////
MODULE ManufacturedSolutionsNSSA
   USE SMConstants
   USE PhysicsStorage_NSSA
   USE Physics_NSSA
   use FluidData_NSSA
   use VariableConversion_NSSA
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

                   rho = exp(x)*1.0_RP/10.0_RP + exp(y)*1.0_RP/10.0_RP + exp(z)*1.0_RP/10.0_RP + 1.0_RP
                   u = COS(PI*(x + y +z )) 
                   v = SIN(PI*(x + y +z ))
                   w = 0.5_RP*COS(PI*(x + y +z )) + 0.5_RP*SIN(PI*(x + y +z ))
                   en = 1.0_RP/(thermodynamics % gamma - 1._RP) + 3.0_RP/2.0_RP
                   thetaeddy = COS(PI*(x + y +z ))


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
            real(kind=RP),             intent(inout)    :: dwall
            real(kind=RP),             intent(in)       :: time
            real(kind=RP),             intent(inout)    :: S(NCONS)

           
   END SUBROUTINE ManufacturedSolutionSourceNSSA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   

END MODULE ManufacturedSolutionsNSSA
