!
!//////////////////////////////////////////////////////
!
!   @File:    template.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jun 27 11:29:03 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module SpallartAlmarasTurbulence
   use SMConstants, only: RP
   use VariableConversion
   use FluidData
   use PhysicsStorage
   use Physics
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
   public SAmodel, InitializeTurbulenceModel, Spalart_Almaras_t

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
      real(kind=RP)              :: sigma = 2.0_RP/3.0_RP
      real(kind=RP)              :: rmax  = 2.0_RP
      real(kind=RP)              :: cw1
      real(kind=RP)              :: stilda

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

         procedure            :: ComputeViscosity    => SA_ComputeViscosity
         procedure            :: ComputeSourceTerms  => SA_Compute_SourceTerms          
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

            write(STD_OUT,'(A,A,A)') "RANS Model: ",trim(modelName), "."

            call model % Initialize(controlVariables)
         
        
      end subroutine InitializeTurbulenceModel


      subroutine SAmodel_Initialize(self, controlVariables)
         implicit none
         class(Spalart_Almaras_t)              :: self
         class(FTValueDictionary),  intent(in) :: controlVariables

         self % cw1  = self % cb1 / POW2(self % kappa)  + (1.0_RP + self % cb2)/ self % sigma

      end subroutine SAmodel_Initialize
!

!/////////////////////////////////////////////////////////
!  
!        Suitable subroutines for Variable_procedure
!        ---------------------------------------------
!
!/////////////////////////////////////////////////////////
!

      subroutine SA_ComputeViscosity(self, rhotheta, kinematic_viscocity, rho, mu_t, eta)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: rhotheta
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: kinematic_viscocity
         real(kind=RP), intent(out)           :: mu_t
         real(kind=RP), intent(out)           :: eta

         real(kind=RP)            :: fv1
         real(kind=RP)            :: chi
         real(kind=RP)            :: theta

         theta = rhotheta / rho

         call Compute_chi(theta, kinematic_viscocity, chi)
         call Compute_fv1(chi, fv1)

         IF (theta .GT. 0.0_RP ) then
            mu_t = rho * theta * fv1
            eta = mu * (1.0_RP + chi) / self % sigma

         ELSE
            mu_t = 0.0_RP       
            eta = mu * (1.0_RP + chi + (chi**2)/2.0_RP) / self % sigma
         END IF 


      end subroutine SA_ComputeViscosity

      !/////////////////////////////////////////////////////////
! Compute Spallart Almaras Source  terms

       subroutine SA_Compute_SourceTerms(self, theta, kinematic_viscocity, rho, dwall,&
                                         Qx, Qy, Qz, source_Kappa)
         implicit none
         !-----------------------------------------------------------
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: kinematic_viscocity
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(in)            :: Qx(NCONS)
         real(kind=RP), intent(in)            :: Qy(NCONS)
         real(kind=RP), intent(in)            :: Qz(NCONS)
         real(kind=RP), intent(out)           :: S_SA
         !-----------------------------------------------------------
         real(kind=RP)           :: chi
         real(kind=RP)           :: vort
         real(kind=RP)           :: fv1
         real(kind=RP)           :: production_G
         real(kind=RP)           :: destruciton_Y
         real(kind=RP)           :: source_Kappa
         real(kind=RP)           :: U_x(NDIM)
         real(kind=RP)           :: U_y(NDIM)
         real(kind=RP)           :: U_z(NDIM)         

         call getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call ComputeVorticity(U_x, U_y, U_z, vort)

         call Compute_chi(theta, kinematic_viscocity, chi)
         call Compute_fv1(chi, fv1)

         call Compute_ProductionTerm(chi, fv1, vort, rho, theta, production_G)

         call Compute_DestructionTerm(theta, rho, dwall, destruciton_Y)

         call Compute_AdditionalSourceTermKappa(rho, Theta_x, Theta_y, Theta_z, source_Kappa)   

         S_SA = production_G - destruciton_Y + source_Kappa    
      
      end subroutine SA_Compute_SourceTerms



!/////////////////////////////////////////////////////////
      subroutine Compute_chi(self, theta, kinematic_viscocity, chi)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: kinematic_viscocity
         real(kind=RP), intent(in)            :: chi

         chi = theta / kinematic_viscocity

      end subroutine Compute_chi

      subroutine Compute_fv1(self, chi, fv1)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(in)            :: fv1

         fv1 = chi**3/(chi**3 + self % cv1 **3)

      end subroutine Compute_fv1

!/////////////////////////////////////////////////////////
! compute production terms 

      subroutine Compute_ProductionTerm(self, chi, fv1, vort, rho, theta, production_G)
         implicit none
         !-----------------------------------------------------------
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(in)            :: fv1
         real(kind=RP), intent(in)            :: vort
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(out)           :: production_G
         !-----------------------------------------------------------
         real(kind=RP)                        :: fv2
         real(kind=RP)                        :: sbar
         real(kind=RP)                        :: gn

         IF (theta .GT. 0.0_RP ) then
            
            fv2    =  Compute_fv2(self, chi, fv1)
            sbar   =  Compute_sbar(theta, fv2)
            self % stilda =  Compute_modifiedvorticity(self, vort, sbar)
            
            production_G = self % cb1 * self % stilda * rho * theta
         
         ELSE
            
            gn = Compute_gn(chi)
            
            production_G = self % cb1 * vort * rho * theta * gn 
         
         END IF
      
      end subroutine Compute_ProductionTerm

      pure function Compute_fv2(self, chi, fv1) result(fv2)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP), intent(in)            :: fv1
         real(kind=RP)                        :: fv2

         fv2 = 1.0_RP - chi / (1.0_RP + chi * fv1) 
     
      end function Compute_fv2

      pure function Compute_sbar(self, theta, fv2) result(sbar)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: fv2
         real(kind=RP)                        :: sbar

         sbar = theta * fv2/ (POW2(self % kappa) + POW2(self % dwall) )
     
      end function Compute_sbar

      pure function Compute_modifiedvorticity(self, vort, sbar) result(stilda)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: fv2
         real(kind=RP), intent(in)            :: sbar
         real(kind=RP)                        :: stilda

         IF (sbar .GT. - self % cv2 * vort ) then
            stilda = vort + sbar 
         ELSE
            stilda =  vort +  (vort*(POW2(self % cv2)*vort + self % cv3 * sbar))/((self % cv3 - 2.0_RP* self % cv2)*vort - sbar)
         END IF     
      end function Compute_modifiedvorticity

      function Compute_gn(self, chi) result (gn)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: chi
         real(kind=RP)                        :: gn

         gn = 1.0_RP - (1000.0_RP*POW2(chi))/(1.0_RP + POW2(chi) )

      end subroutine Compute_gn

!/////////////////////////////////////////////////////////
! compute destruciton terms 
      
      subroutine Compute_DestructionTerm(self, theta, rho, dwall, destruciton_Y)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: fw
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(in)            :: y_destr

         IF (theta .GT. 0.0_RP ) then
         
         g  = Compute_g(theta, dwall)
         fw = Compute_fw(g)
         
         destruciton_Y =    self % cw1 * fw * (rho*POW2(theta))/(POW2(dwall)) 
         
         ELSE
         
         destruciton_Y =  - self % cw1 * (rho*POW2(theta))/(POW2(dwall)) 
      
         END IF
      
      end subroutine Compute_ProductionTerm

      function Compute_g(self, theta, dwall) result(g)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(out)           :: g

         real(kind=RP)                        :: r

         r = min(theta/( self % stilda * POW2(self % kappa) * POW2(dwall)), self % rmax )

         g = r + self % cw2 * ( r**6 - r)

      end function Compute_g

      function Compute_fw(self, g) result(fw)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: theta
         real(kind=RP), intent(in)            :: stilda
         real(kind=RP), intent(in)            :: dwall
         real(kind=RP), intent(out)            :: g

         fw = g * (1.0_RP + self % cw3**6 / g**6 + self % cw3**6)**(1.0_RP/6.0_RP)

      end function Compute_fw

!/////////////////////////////////////////////////////////
! compute additional source  terms

       subroutine Compute_AdditionalSourceTermKappa(self, rho, Theta_x, Theta_y, Theta_z, source_Kappa)
         implicit none
         class(Spalart_Almaras_t), intent(in) :: self
         real(kind=RP), intent(in)            :: rho
         real(kind=RP), intent(in)            :: Theta_x
         real(kind=RP), intent(in)            :: Theta_y
         real(kind=RP), intent(in)            :: Theta_z
         real(kind=RP), intent(in)            :: Kappa

         source_Kappa  = self % cb2 * rho * sum (POW2(Theta_x) + POW2(Theta_y) + POW2(Theta_z))
               
      end subroutine Compute_AdditionalSourceTermKappa

   
end module SpallartAlmarasTurbulences
