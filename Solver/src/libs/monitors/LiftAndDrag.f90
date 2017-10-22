!module LiftAndDrag
!   use SMConstants
!   use ZoneClass
!   use ElementClass
!   use FaceClass
!   use Physics
!   use HexMeshClass
!   use NodalStorageClass
!   use ProlongToFacesProcedures
!   implicit none
!   
!contains
!!
!!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!!
!!  ------------------------------------
!!  Computes lift and drag for all zones
!!  ------------------------------------
!   subroutine calc_LiftDrag(mesh,spA)
!      implicit none
!      !------------------------------------------------
!      type(HexMesh)  , intent(inout)   :: mesh
!      type(NodalStorage), intent(in)   :: spA(0:,0:,0:)
!      !------------------------------------------------ 
!      integer       :: iZone                    ! Zone counter
!      real(kind=RP) :: Fx, Fy, Fz               ! Nondimensional forces in cartesian directions
!      real(kind=RP) :: Fx_abs, Fy_abs, Fz_abs   ! Absolute forces in cartesian directions
!      real(kind=RP) :: Fd, Fl                   ! Drag and lift (nondimensional) forces
!      real(kind=RP) :: Fd_abs, Fl_abs           ! Absolute forces in cartesian directions
!      real(kind=RP) :: Cd, Cl                   ! Drag and lift coefficients
!      !------------------------------------------------
!      
!      WRITE(STD_OUT,'(A45)') "*********************************************"
!      WRITE(STD_OUT,'(A1,43X,A1)') "*","*"
!      WRITE(STD_OUT,'(A1,15X,A13,15X,A1)') "*","Zone Analysis","*"
!      WRITE(STD_OUT,'(A1,43X,A1)') "*","*"
!
!      do iZone = 1, size(mesh % zones)
!         
!         call calc_LiftDrag_zone(mesh, spA, mesh % zones(iZone),Fx, Fy, Fz)
!!
!!        ---------------------------------------------------------
!!        Get the actual forces from the nondimensional computation 
!!        ---------------------------------------------------------
!!
!         Fx_abs = Fx * refValues % rho * refValues % V**2 * refValues % L**2
!         Fy_abs = Fy * refValues % rho * refValues % V**2 * refValues % L**2
!         Fz_abs = Fz * refValues % rho * refValues % V**2 * refValues % L**2
!
!!
!!        ---------------------------------------------------------
!!        Get drag and lift forces
!!          TODO: Add routine to project forces onto flow direction and some perpendicular (or other user specified direction), for drag and lift
!!                Meanwhile, we simply assign Fd = Fx and Fl = Fy
!!        ---------------------------------------------------------
!!
!         Fd = Fx
!         Fl = Fy
!         
!         Fd_abs = Fd * refValues % rho * refValues % V**2 * refValues % L**2
!         Fl_abs = Fl * refValues % rho * refValues % V**2 * refValues % L**2
!         
!!
!!        ------------------------------------------------------
!!        Get drag and lift coefficients. Remember:
!!        
!!                        2 * Fd
!!           Cd = -------------------------
!!                 rho_inf * V_inf^2 * A
!!
!!         Here, we will suppose A = L * 1. (drag per unit span)
!!        ------------------------------------------------------
!!
!         
!         Cd = Fd * 2._RP * refValues % L
!         Cl = Fl * 2._RP * refValues % L
!         
!         
!         WRITE(STD_OUT,'(A45)') "*********************************************"
!         WRITE(STD_OUT,'(A8,I3,A12,A20,A2)') "* ZONE: ", iZone, ". Boundary: ", trim(mesh % zones(iZone) % Name) , " *"
!         WRITE(STD_OUT,'(A7,ES12.5,A2,22X,A2)') "* Fx = ", Fx_abs," N"," *"
!         WRITE(STD_OUT,'(A7,ES12.5,A2,22X,A2)') "* Fy = ", Fy_abs," N"," *"
!         WRITE(STD_OUT,'(A7,ES12.5,A2,22X,A2)') "* Fz = ", Fz_abs," N"," *"
!         WRITE(STD_OUT,'(A7,ES12.5,A2,A10,ES12.5,A2)') "* Fd = ", Fd_abs," N","  |  Cd = ",Cd," *"
!         WRITE(STD_OUT,'(A7,ES12.5,A2,A10,ES12.5,A2)') "* Fl = ", Fl_abs," N","  |  Cl = ",Cl," *"
!      end do
!      
!      WRITE(STD_OUT,'(A45)') "*********************************************"
!      
!   end subroutine
!!
!!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!!
!!  -----------------------------------
!!  Computes lift and drag for one zone
!!  -----------------------------------
!   subroutine calc_LiftDrag_zone(mesh,spA,zone,Fx, Fy, Fz)
!      implicit none
!      !------------------------------------------------
!      type(HexMesh), intent(inout)      :: mesh
!      type(NodalStorage), intent(in)    :: spA(0:,0:,0:)
!      type(Zone_t) , intent(in)         :: zone
!      real(kind=RP), intent(out)        :: Fx, Fy, Fz
!      !------------------------------------------------
!      integer       :: iface       ! Face counter
!      integer       :: fID, eID    ! Face and element indexes
!      real(kind=RP) :: FxFace,FyFace,FzFace
!      !------------------------------------------------
!      
!      Fx = 0._RP
!      Fy = 0._RP
!      Fz = 0._RP
!      
!!$omp parallel do reduction(+:Fx,Fy,Fz) private(fID,eID,FxFace,FyFace,FzFace) schedule(guided)
!      do iface = 1, zone % no_of_faces
!         fID = zone % faces(iface)
!         eID = mesh % faces(fID) % elementIDs(1)
!         
!         call calc_LiftDrag_face (mesh % elements(eID), &
!                                  spA, &
!                                  mesh % faces(fID) % elementSide(1),&
!                                  FxFace,FyFace,FzFace)
!         Fx = Fx + FxFace
!         Fy = Fy + FyFace
!         Fz = Fz + FzFace
!      end do
!!$omp end parallel do
!      
!   end subroutine calc_LiftDrag_zone
!!
!!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!!
!!  -----------------------------------------------
!!  Computes lift and drag contribution of one face
!!  -----------------------------------------------
!   subroutine calc_LiftDrag_face(eL,spA,elSide,Fx,Fy,Fz)
!      implicit none
!      !------------------------------------------------
!      type(Element)     , intent(inout) :: eL        !<> Element adjacent to face
!      type(NodalStorage), TARGET, intent(in) :: spA(0:,0:,0:)
!      integer           , intent(in)    :: elSide    !<  element side: Index of face
!      real(kind=RP)     , intent(out)   :: Fx,Fy,Fz  !>  Cartesian components of the force
!      !------------------------------------------------
!      integer                :: N(2)               ! polynomial orders on face
!      integer                :: Nx, Ny, Nz         ! polynomial orders of adjacent element
!      integer                :: i,j                ! coordinate counters
!      real(kind=RP)          :: ds                 ! mapping Jacobian term on each node of the surface
!      real(kind=RP)          :: p                  ! pressure
!      real(kind=RP)          :: norm(3)            ! unit normal vector
!      real(kind=RP), pointer :: wx(:), wy(:)       ! Weights on the face (relative coordinates)
!      real(kind=RP)          :: Fx_p, Fy_p, Fz_p   ! Pressure forces
!      real(kind=RP)          :: Fx_v, Fy_v, Fz_v   ! Viscous forces
!      real(kind=RP)          :: F(N_EQN,3)         ! Viscous fluxes
!      real(kind=RP)          :: Qb(N_EQN)          ! Conserved variables at the boundary
!      !------------------------------------------------
!      
!!
!!     -----------------
!!     Basic definitions
!!     -----------------
!!
!      Nx = eL % Nxyz(1)
!      Ny = eL % Nxyz(2)
!      Nz = eL % Nxyz(3)
!      N = eL % Nxyz (axisMap(:,elSide))
!      
!      Fx_p = 0._RP
!      Fy_p = 0._RP
!      Fz_p = 0._RP
!      Fx_v = 0._RP
!      Fy_v = 0._RP
!      Fz_v = 0._RP
!      
!      ! Get the appropriate weights for this face
!      select case (elSide)
!         case(1,2)
!            wx => spA(Nx,Ny,Nz) % wx
!            wy => spA(Nx,Ny,Nz) % wz
!         case(3,5)
!            wx => spA(Nx,Ny,Nz) % wx
!            wy => spA(Nx,Ny,Nz) % wy
!         case(4,6)
!            wx => spA(Nx,Ny,Nz) % wy
!            wy => spA(Nx,Ny,Nz) % wz
!      end select
!      
!!
!!     --------------------------
!!     Prolong variables to faces
!!     --------------------------
!!      
!      call ProlongToFaces( eL, spA(Nx,Ny,Nz) )
!      if (flowIsNavierStokes) call ProlongGradientToFaces( eL, spA(Nx,Ny,Nz) )
!      
!      ! Numerical integration over surface
!      DO j = 0, N(2)
!         DO i = 0, N(1)
!            ds = eL % geom % scal (i,j,elSide)
!            Qb = eL % Qb (:,i,j,elSide)
!            norm = eL % geom % normal (:,i,j,elSide)
!            
!!           ---------------------------------
!!           Add pressure contribution
!!           (p at infinity is 1.0_RP/gammaM2)
!!           ---------------------------------
!            p = (Pressure(Qb) - 1.0_RP/dimensionless % gammaM2)  ! *2._RP !????
!            
!            Fx_p = Fx_p + p * norm(1) * wx(i)*wy(j)*ds
!            Fy_p = Fy_p + p * norm(2) * wx(i)*wy(j)*ds
!            Fz_p = Fz_p + p * norm(3) * wx(i)*wy(j)*ds
!            
!            
!!           ---------------------------------
!!           Add viscous contribution
!!           ---------------------------------
!            
!            if (flowIsNavierStokes) then
!               F = ViscousFlux0D (Qb , &
!                                  eL % U_xb (:,i,j,elSide), &
!                                  eL % U_yb (:,i,j,elSide), &
!                                  eL % U_zb (:,i,j,elSide)) 
!               
!               ! Negative sign because the boundary "feels" a force equal and opposite to the one that the fluid "feels" (3rd law)
!               Fx_v = Fx_v - (F(IRHOU,1) * norm(1) + F(IRHOU,2) * norm(2) + F(IRHOU,3) * norm(3) ) * wx(i)*wy(j)*ds
!               Fy_v = Fy_v - (F(IRHOV,1) * norm(1) + F(IRHOV,2) * norm(2) + F(IRHOV,3) * norm(3) ) * wx(i)*wy(j)*ds
!               Fz_v = Fz_v - (F(IRHOW,1) * norm(1) + F(IRHOW,2) * norm(2) + F(IRHOW,3) * norm(3) ) * wx(i)*wy(j)*ds
!               
!            end if
!            
!         end do
!      end do
!      
!      nullify(wx,wy)
!      
!      Fx = Fx_p + Fx_v
!      Fy = Fy_p + Fy_v
!      Fz = Fz_p + Fz_v
!      
!   end subroutine calc_LiftDrag_face
!   
!end module LiftAndDrag
