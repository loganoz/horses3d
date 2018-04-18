module EllipticDiscretizationClass
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use DGSEMClass, only: BCState_FCN
   implicit none
!
   private
   public   EllipticDiscretization_t, EllipticFlux0D_f, EllipticFlux2D_f, EllipticFlux3D_f

   public BaseClass_ComputeGradient

   type EllipticDiscretization_t
      contains
         procedure      :: Initialize                => BaseClass_Initialize
         procedure      :: ComputeGradient           => BaseClass_ComputeGradient
         procedure      :: ComputeInnerFluxes        => BaseClass_ComputeInnerFluxes
         procedure      :: ComputeInnerFluxesWithSGS => BaseClass_ComputeInnerFluxesWithSGS
         procedure      :: RiemannSolver             => BaseClass_RiemannSolver
         procedure      :: RiemannSolverWithSGS      => BaseClass_RiemannSolverWithSGS
         procedure      :: Describe                  => BaseClass_Describe
   end type EllipticDiscretization_t

   abstract interface
      pure subroutine EllipticFlux0D_f(Q, U_x, U_y, U_z, mu, kappa, F)
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)  :: Q   (1:NCONS     )
         real(kind=RP), intent(in)  :: U_x (1:NGRAD)
         real(kind=RP), intent(in)  :: U_y (1:NGRAD)
         real(kind=RP), intent(in)  :: U_z (1:NGRAD)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:NCONS, 1:NDIM)
      end subroutine EllipticFlux0D_f

      pure subroutine EllipticFlux2D_f( N, Q, U_x, U_y, U_z, mu, kappa, F)
         use SMConstants
         use PhysicsStorage
         implicit none
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:NGRAD, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:NGRAD, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:NGRAD, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 1:NDIM, 0:N(1), 0:N(2))
      end subroutine EllipticFlux2D_f

      pure subroutine EllipticFlux3D_f( N, Q, U_x, U_y, U_z, mu, kappa, F)
         use SMConstants
         use PhysicsStorage
         implicit none
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:NGRAD, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:NGRAD, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:NGRAD, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 0:N(1), 0:N(2), 0:N(3), 1:NDIM )
      end subroutine EllipticFlux3D_f

   end interface

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
         class(EllipticDiscretization_t)                :: self
         class(FTValueDictionary),  intent(in) :: controlVariables

      end subroutine BaseClass_Initialize

      subroutine BaseClass_Describe(self)
         implicit none
         class(EllipticDiscretization_t), intent(in)  :: self

      end subroutine BaseClass_Describe

      subroutine BaseClass_ComputeGradient( self , mesh , time , externalStateProcedure)
!
!        *****************************************************
!           BaseClass computes Local Gradients by default
!        *****************************************************
!           
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(EllipticDiscretization_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(BCState_FCN)           :: externalStateProcedure
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

      subroutine BaseClass_ComputeInnerFluxes( self , e , EllipticFlux, contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(EllipticDiscretization_t) ,  intent (in) :: self
         type(Element)                                  :: e
         procedure(EllipticFlux3D_f)                    :: EllipticFlux
         real(kind=RP)           ,  intent (out)        :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxes

      subroutine BaseClass_ComputeInnerFluxesWithSGS( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(EllipticDiscretization_t) ,  intent (in)   :: self
         type(Element)                           :: e
         real(kind=RP)           ,  intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxesWithSGS

      subroutine BaseClass_RiemannSolver ( self, f, EllipticFlux, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                           nHat, dWall, flux )
         use SMConstants
         use PhysicsStorage
         use FaceClass
         implicit none
         class(EllipticDiscretization_t)         :: self
         class(Face),   intent(in)               :: f
         procedure(EllipticFlux0D_f)             :: EllipticFlux
         real(kind=RP), dimension(NCONS)         :: QLeft
         real(kind=RP), dimension(NCONS)         :: QRight
         real(kind=RP), dimension(NGRAD)    :: U_xLeft
         real(kind=RP), dimension(NGRAD)    :: U_yLeft
         real(kind=RP), dimension(NGRAD)    :: U_zLeft
         real(kind=RP), dimension(NGRAD)    :: U_xRight
         real(kind=RP), dimension(NGRAD)    :: U_yRight
         real(kind=RP), dimension(NGRAD)    :: U_zRight
         real(kind=RP), dimension(NDIM)          :: nHat
         real(kind=RP)                           :: dWall
         real(kind=RP), dimension(NCONS)         :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolver

      subroutine BaseClass_RiemannSolverWithSGS ( self, f, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                                  nHat, dWall, flux )
         use SMConstants
         use PhysicsStorage
         use FaceClass
         implicit none
         class(EllipticDiscretization_t)               :: self
         class(Face),   intent(in)            :: f
         real(kind=RP), dimension(NCONS)      :: QLeft
         real(kind=RP), dimension(NCONS)      :: QRight
         real(kind=RP), dimension(NGRAD) :: U_xLeft
         real(kind=RP), dimension(NGRAD) :: U_yLeft
         real(kind=RP), dimension(NGRAD) :: U_zLeft
         real(kind=RP), dimension(NGRAD) :: U_xRight
         real(kind=RP), dimension(NGRAD) :: U_yRight
         real(kind=RP), dimension(NGRAD) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP)                        :: dWall
         real(kind=RP), dimension(NCONS)      :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolverWithSGS
end module EllipticDiscretizationClass
