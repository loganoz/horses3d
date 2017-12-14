module ViscousMethodClass
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   implicit none
!
!
   private
   public   ViscousMethod_t 

   type ViscousMethod_t
      contains
         procedure      :: Initialize         => BaseClass_Initialize
         procedure      :: ComputeGradient    => BaseClass_ComputeGradient
         procedure      :: ComputeInnerFluxes => BaseClass_ComputeInnerFluxes
         procedure      :: RiemannSolver      => BaseClass_RiemannSolver
   end type ViscousMethod_t
!
!  ========
   contains
!  ========
!
      subroutine BaseClass_Initialize(self, controlVariables)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t)                :: self
         class(FTValueDictionary),  intent(in) :: controlVariables
         interface
            subroutine toLower(str)
               character(*), intent(in out) :: str
            end subroutine toLower
         end interface

      end subroutine BaseClass_Initialize
      subroutine BaseClass_ComputeGradient( self , mesh , time , externalStateProcedure , externalGradientsProcedure)
!
!        *****************************************************
!           BaseClass computes Local Gradients by default
!        *****************************************************
!           
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(ViscousMethod_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         external                         :: externalStateProcedure
         external                         :: externalGradientsProcedure
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID

!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % ComputeLocalGradient
         end do
!$omp end do

      end subroutine BaseClass_ComputeGradient

      subroutine BaseClass_ComputeInnerFluxes( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t) ,  intent (in)   :: self
         type(Element)                           :: e
         real(kind=RP)           ,  intent (out) :: contravariantFlux(1:N_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxes

      subroutine BaseClass_RiemannSolver ( self, f, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                           nHat, flux )
         use SMConstants
         use PhysicsStorage
         use FaceClass
         implicit none
         class(ViscousMethod_t)               :: self
         class(Face),   intent(in)            :: f
         real(kind=RP), dimension(N_EQN)      :: QLeft
         real(kind=RP), dimension(N_EQN)      :: QRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP), dimension(N_EQN)      :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolver
end module ViscousMethodClass
