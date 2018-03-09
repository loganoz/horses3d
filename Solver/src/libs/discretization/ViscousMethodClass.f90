module ViscousMethodClass
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use DGSEMClass, only: BCState_FCN
   implicit none
!
   private
   public   ViscousMethod_t, BaseClass_ComputeGradient

   type ViscousMethod_t
      contains
         procedure      :: Initialize                => BaseClass_Initialize
         procedure      :: ComputeGradient           => BaseClass_ComputeGradient
         procedure      :: ComputeInnerFluxes        => BaseClass_ComputeInnerFluxes
         procedure      :: ComputeInnerFluxJacobian  => BaseClass_ComputeInnerFluxJacobian
         procedure      :: ComputeInnerFluxesWithSGS => BaseClass_ComputeInnerFluxesWithSGS
         procedure      :: RiemannSolver             => BaseClass_RiemannSolver
         procedure      :: RiemannSolverWithSGS      => BaseClass_RiemannSolverWithSGS
         procedure      :: Describe                  => BaseClass_Describe
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

      end subroutine BaseClass_Initialize

      subroutine BaseClass_Describe(self)
         implicit none
         class(ViscousMethod_t), intent(in)  :: self

      end subroutine BaseClass_Describe

      subroutine BaseClass_ComputeGradient( self , mesh , time , externalStateProcedure)
!
!        *****************************************************
!           BaseClass computes Local Gradients by default
!              Do not change.. Used by ComputeTimeDerivativeIsolated
!        *****************************************************
!           
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(ViscousMethod_t), intent(in) :: self
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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------------------------------
!     Subroutine to get the Jacobian of the contravariant fluxes (with respect to ∇q) on the inner Gauss points of an element
!     -> The contravariant Jacobian is a 4th order tensor in every Gauss point
!     -> dFdQ (:,:,:,i,j,k,dim)
!              |_| | |_|_| |
!               |  |  |    |_flux direction: f(1), g(2), h(3)
!               |  |  |______Coordinate indexes in element 
!               |  |_________∇q component: 1, 2, 3
!               |____________Jacobian for this component
!     --------------------------------------------------------------------------------------------------
      subroutine BaseClass_ComputeInnerFluxJacobian( self, e, dFdQ) 
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(ViscousMethod_t), intent(in)  :: self
         type(Element)         , intent(in)  :: e                                                                             !<  This element
         real(kind=RP)         , intent(out) :: dFdQ( NCONS, NCONS, NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3), NDIM ) !>  Contravariant Jacobians for each Gaus point in the element
         !--------------------------------------------
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM,NDIM) :: df_dgradq   ! Cartesian Jacobian tensor
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM,NDIM) :: df_dgradq_  ! Intermediate Jacobian tensor (between the cartesian and the contravariant)
         integer :: i,j,k     ! Coordinate counters
         integer :: i1, i2    ! Index of G_xx
         !--------------------------------------------
#if defined(NAVIERSTOKES)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            
            call ViscousJacobian(e % storage % Q(:,i,j,k), df_dgradq)
            
!           Fill intermedate Jacobian tensor
!           --------------------------------
            df_dgradq_ = 0._RP
            do i1 = 1, NDIM ; do i2 = 1, NDIM
               df_dgradq_(:,:,i1,1) = df_dgradq_(:,:,i1,1) + df_dgradq(:,:,i2,i1) * e % geom % jGradXi  (i2,i,j,k)
            end do          ; end do
            do i1 = 1, NDIM ; do i2 = 1, NDIM
               df_dgradq_(:,:,i1,2) = df_dgradq_(:,:,i1,2) + df_dgradq(:,:,i2,i1) * e % geom % jGradEta (i2,i,j,k)
            end do          ; end do
            do i1 = 1, NDIM ; do i2 = 1, NDIM
               df_dgradq_(:,:,i1,3) = df_dgradq_(:,:,i1,3) + df_dgradq(:,:,i2,i1) * e % geom % jGradZeta(i2,i,j,k)
            end do          ; end do
            df_dgradq_ = df_dgradq_ * e % geom % invJacobian(i,j,k)
            
!           Fill Finished Jacobian contravariant tensor
!           *******************************************
            
!           Xi-flux
!           -------
            dFdQ(:,:,:,i,j,k,1) = 0._RP
            do i1 = 1, NDIM ; do i2 = 1, NDIM
               dFdQ(:,:,i1,i,j,k,1) = dFdQ(:,:,i1,i,j,k,1) + df_dgradq_(:,:,i2,i1) * e % geom % jGradXi  (i2,i,j,k)
            end do          ; end do
            
!           Eta-flux
!           --------
            dFdQ(:,:,:,i,j,k,2) = 0._RP
            do i1 = 1, NDIM ; do i2 = 1, NDIM
               dFdQ(:,:,i1,i,j,k,2) = dFdQ(:,:,i1,i,j,k,2) + df_dgradq_(:,:,i2,i1) * e % geom % jGradEta (i2,i,j,k)
            end do          ; end do
            
!           Zeta-flux
!           ---------
            dFdQ(:,:,:,i,j,k,3) = 0._RP
            do i1 = 1, NDIM ; do i2 = 1, NDIM
               dFdQ(:,:,i1,i,j,k,3) = dFdQ(:,:,i1,i,j,k,3) + df_dgradq_(:,:,i2,i1) * e % geom % jGradZeta(i2,i,j,k)
            end do          ; end do
            
         end do                ; end do                ; end do
#endif
      end subroutine BaseClass_ComputeInnerFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeInnerFluxesWithSGS( self , e , contravariantFlux )
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

      end subroutine BaseClass_ComputeInnerFluxesWithSGS

      subroutine BaseClass_RiemannSolver ( self, f, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                           nHat, dWall, flux )
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
         real(kind=RP)                        :: dWall
         real(kind=RP), dimension(N_EQN)      :: flux
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
         real(kind=RP)                        :: dWall
         real(kind=RP), dimension(N_EQN)      :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolverWithSGS
end module ViscousMethodClass
