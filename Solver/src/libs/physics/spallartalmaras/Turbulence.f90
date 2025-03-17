#include "Includes.h"
module SpallartAlmarasTurbulence
   use SMConstants
   use VariableConversion_NSSA
   use FluidData_NSSA
   use PhysicsStorage_NSSA
   use FTValueDictionaryClass
   use Utilities
   use mainKeywordsModule
   use Headers
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
   public SAmodel, InitializeTurbulenceModel, Spalart_Almaras_t, SA_ComputeViscosity

!  ****************
!  Class definition
!  ****************
!
   type Spalart_Almaras_t
      logical                    :: constructed = .false.
      real(kind=RP)              :: cv1 = 7.1_RP
      real(kind=RP)              :: cv2 = 0.7_RP
      real(kind=RP)              :: cv3 = 0.9_RP
      real(kind=RP)              :: cb1 = 0.1355_RP
      real(kind=RP)              :: cb2 = 0.622_RP
      real(kind=RP)              :: cw2 = 0.3_RP
      real(kind=RP)              :: cw3 = 2.0_RP
      real(kind=RP)              :: kappa = 0.41_RP
      real(kind=RP)              :: sigma_sa = 2.0_RP/3.0_RP
      real(kind=RP)              :: rmax  = 2.0_RP
      real(kind=RP)              :: cw1

      contains
!        procedure         :: Destruct          => Class_Destruct
         procedure         :: Initialize        => SAmodel_Initialize
         
         procedure, private ::      Compute_chi
         procedure, private ::      Compute_fv1
         procedure, private ::      Compute_fv2
         procedure, private ::      Compute_sbar
         procedure, private ::      Compute_modifiedvorticity
         procedure, private ::      Compute_gn
         procedure, private ::      Compute_g
         procedure, private ::      Compute_fw
         procedure, private ::      Compute_ProductionTerm
         procedure, private ::      Compute_DestructionTerm
         procedure, private ::      Compute_AdditionalSourceTermKappa

         procedure, public  :: ComputeViscosity    => SA_ComputeViscosity
         procedure          :: ComputeSourceTerms  => SA_Compute_SourceTerms          
   end type Spalart_Almaras_t
!

!
   class(Spalart_Almaras_t), allocatable   :: SAmodel

!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
      subroutine InitializeTurbulenceModel(model, controlVariables)
         implicit none
         class(Spalart_Almaras_t), allocatable        :: model
         class(FTValueDictionary),  intent(in) :: controlVariables
!
!        ---------------
!        Local variables         
!        ---------------

         if (.not. allocated(model)) allocate( Spalart_Almaras_t :: model)

         call model % Initialize(controlVariables)
            
      end subroutine InitializeTurbulenceModel


      subroutine SAmodel_Initialize(self, controlVariables)
         implicit none
         class(Spalart_Almaras_t)              :: self
         class(FTValueDictionary),  intent(in) :: controlVariables

         self % cw1  = self % cb1 / POW2(self % kappa)  + (1.0_RP + self % cb2)/ self % sigma_sa

      end subroutine SAmodel_Initialize
!

!/////////////////////////////////////////////////////////
!  
!        Suitable subroutines for Variable_procedure
!        ---------------------------------------------
!
!/////////////////////////////////////////////////////////
!

      subroutine SA_ComputeViscosity(self, rhotheta, kinematic_viscocity, rho, mu, mu_t, eta, xvec)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: rhotheta
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: kinematic_viscocity
         real(kind=RP), intent(in)            :: mu
         real(kind=RP), intent(inout)         :: mu_t
         real(kind=RP), intent(inout)         :: eta
         real(kind=RP), intent(in) , optional :: xvec(3)

         real(kind=RP)            :: fv1
         real(kind=RP)            :: chi
         real(kind=RP)            :: theta, etaexact

         theta = rhotheta / rho
         call self % Compute_chi(theta, kinematic_viscocity, chi)
         call self % Compute_fv1(chi, fv1)

         if (theta .GE. 0.0_RP ) then
            mu_t = rho * theta * fv1 
            eta  = dimensionless % mu * mu * (1.0_RP + chi) / self % sigma_sa
        else
             mu_t = 0.0_RP       
             eta  = dimensionless % mu * mu * (1.0_RP + chi + (POW2(chi))/2.0_RP) / self % sigma_sa
        endif 

      end subroutine SA_ComputeViscosity

      !/////////////////////////////////////////////////////////
! Compute Spallart Almaras Source  terms

       subroutine SA_Compute_SourceTerms(self, rhotheta, kinematic_viscocity, rho, dwall,&
                                         Q, Q_x, Q_y, Q_z, S_SA, xvec)
         implicit none
         !-----------------------------------------------------------
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: rhotheta
         real(kind=RP), intent(in)            :: kinematic_viscocity
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(in)            :: Q(NCONS)
         real(kind=RP), intent(in)            :: Q_x(NCONS)
         real(kind=RP), intent(in)            :: Q_y(NCONS)
         real(kind=RP), intent(in)            :: Q_z(NCONS)
         real(kind=RP), intent(out)           :: S_SA(NCONS)
         real(kind=RP), intent(in) , optional :: xvec(NDIM)

         !-----------------------------------------------------------
         real(kind=RP)           :: chi , sbar, Production, thetaexact, rhoexact, DESTRUCTION2, DESTRUCTION1, DESTRUCTION3
         real(kind=RP)           :: theta
         real(kind=RP)           :: vort, somega
         real(kind=RP)           :: fv1
         real(kind=RP)           :: production_G
         real(kind=RP)           :: destruction_Y
         real(kind=RP)           :: source_Kappa
         real(kind=RP)           :: U_x(NDIM)
         real(kind=RP)           :: U_y(NDIM)
         real(kind=RP)           :: U_z(NDIM)         
         real(kind=RP)           :: Theta_x, Theta_y, Theta_z
         theta = rhotheta / rho
         call getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call ComputeVorticity(U_x, U_y, U_z, vort)
         call geteddyviscositygradients(Q, Q_x, Q_y, Q_z, Theta_x, Theta_y, Theta_z)

         call self % Compute_chi(theta, kinematic_viscocity, chi)
         call self % Compute_fv1(chi, fv1)

         if (dwall .GT. 0.0_RP) then 
            call self % Compute_ProductionTerm(chi, fv1, vort, rho, theta, dwall,  production_G)
            call self % Compute_DestructionTerm(chi, fv1, vort, rho, theta, dwall,  destruction_Y)       
         else
            production_G = 0.0_RP
            destruction_Y = 0.0_RP
         endif

         call self % Compute_AdditionalSourceTermKappa(rho, Theta_x, Theta_y, Theta_z, source_Kappa)   

         S_SA(1:5) = 0.0_RP
         S_SA(6)   = production_G - destruction_Y + source_Kappa 
            
     end subroutine SA_Compute_SourceTerms



!/////////////////////////////////////////////////////////
      subroutine Compute_chi(self, theta, kinematic_viscocity, chi)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: kinematic_viscocity
         real(kind=RP), intent(out)           :: chi

         chi = theta / kinematic_viscocity

      end subroutine Compute_chi

      subroutine Compute_fv1(self, chi, fv1)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(out)            :: fv1

         fv1 = (chi**3.0_RP)/((chi**3.0_RP) + (self % cv1)**3.0_RP)

      end subroutine Compute_fv1

!/////////////////////////////////////////////////////////
! compute production terms 

      subroutine Compute_ProductionTerm(self, chi, fv1, vort, rho, theta, dwall, production_G)
         implicit none
         !-----------------------------------------------------------
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(in)            :: fv1
         real(kind=RP), intent(in)            :: vort
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(inout)         :: production_G
         !-----------------------------------------------------------
         real(kind=RP)                        :: fv2
         real(kind=RP)                        :: sbar, stilda 
         real(kind=RP)                        :: gn, PRODUCTION1, PRODUCTION2, PRODUCTION3

      if (theta .GE. 0.0_RP ) then

            fv2    =  Compute_fv2(self, chi, fv1)
            sbar   =  Compute_sbar(self, theta, dwall, fv2)
            stilda =  Compute_modifiedvorticity(self, vort, sbar)
            production_G = self % cb1 * stilda * rho * theta
            
         else
            gn = Compute_gn(self, chi)
            production_G = self % cb1 * vort * rho * theta * gn 
            
         end if 

      end subroutine Compute_ProductionTerm

      function Compute_fv2(self, chi, fv1) result(fv2)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(in)            :: fv1
         real(kind=RP)                        :: fv2

         fv2 = 1.0_RP - chi / (1.0_RP + chi * fv1) 
     
      end function Compute_fv2

      function Compute_sbar(self, theta, dwall, fv2) result(sbar)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: fv2
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP)                        :: sbar

         sbar =  dimensionless % mu * theta * fv2/ (POW2(self % kappa) * POW2(dwall) ) 
     
      end function Compute_sbar

      function Compute_modifiedvorticity(self, vort, sbar) result(stilda)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: vort
         real(kind=RP), intent(in)            :: sbar
         real(kind=RP)                        :: stilda

         IF (sbar .GE. - self % cv2 * vort) then
            stilda = vort + sbar 
         ELSE
            stilda =  vort +  (vort*(POW2(self % cv2)*vort + self % cv3 * sbar))/((self % cv3 - 2.0_RP* self % cv2)*vort - sbar)
         END IF     
      end function Compute_modifiedvorticity

      function Compute_gn(self, chi) result (gn)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP)                        :: gn

         gn = 1.0_RP - (1000.0_RP*POW2(chi))/(1.0_RP + POW2(chi) )

     end function Compute_gn

!/////////////////////////////////////////////////////////
! compute destruction terms 
      
      subroutine Compute_DestructionTerm(self, chi, fv1, vort, rho, theta, dwall, destruction_Y)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(in)            :: fv1
         real(kind=RP), intent(in)            :: vort
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(inout)         :: destruction_Y
         
         real(kind=RP)            :: g
         real(kind=RP)            :: fw
         real(kind=RP)            :: fv2
         real(kind=RP)            :: sbar, stilda 

         if (theta .GE. 0.0_RP ) then
            fv2    =  Compute_fv2(self, chi, fv1)
            sbar   =  Compute_sbar(self, theta, dwall, fv2)
            stilda =  Compute_modifiedvorticity(self, vort, sbar)
            g  = Compute_g(self, theta, dwall, stilda)
            fw = Compute_fw(self, g)
            destruction_Y =    dimensionless % mu * self % cw1 * fw * (rho*POW2(theta))/(POW2(dwall))          
         else
            destruction_Y =  - dimensionless % mu * self % cw1 * (rho*POW2(theta))/(POW2(dwall)) 
         end if 

      end subroutine Compute_DestructionTerm

      function Compute_g(self, theta, dwall, stilda) result(g)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(in)            :: stilda
         real(kind=RP)                        :: g
         real(kind=RP)                        :: r
         
         if (stilda .EQ. 0.0_RP) then
            r = self % rmax
         else 
            r = min(dimensionless % mu * theta/(stilda * POW2(self % kappa) * POW2(dwall)), self % rmax )
         endif

         g = r + self % cw2 * ( (r**6.0_RP) - r)

      end function Compute_g

      function Compute_fw(self, g) result(fw)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: g
         real(kind=RP)            :: fw

         fw = g * ((1.0_RP + self % cw3**6.0_RP) / (g**6.0_RP + self % cw3**6.0_RP))**(1.0_RP/6.0_RP)

      end function Compute_fw

!/////////////////////////////////////////////////////////
! compute additional source  terms

       subroutine Compute_AdditionalSourceTermKappa(self, rho, Theta_x, Theta_y, Theta_z, source_Kappa)
         implicit none
         class(Spalart_Almaras_t), intent(inout) :: self
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: Theta_x
         real(kind=RP), intent(in)            :: Theta_y
         real(kind=RP), intent(in)            :: Theta_z
         real(kind=RP), intent(inout)           :: source_Kappa

         source_Kappa  = dimensionless % mu * self % cb2 * rho * (POW2(Theta_x) + POW2(Theta_y) + POW2(Theta_z)) / self % sigma_sa
               
      end subroutine Compute_AdditionalSourceTermKappa

   
end module SpallartAlmarasTurbulence