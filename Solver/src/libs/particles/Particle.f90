!
!//////////////////////////////////////////////////////
!
!   @File:    Particle.f90
!   @Author:  Gonzalo (g.rubio@upm.es)
!   @Created: Tue Apr 10 17:31:21 2018
!   @Last revision date: Mon Jul  2 14:17:27 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7af1f42fb2bc9ea3a0103412145f2a925b4fac5e
!
!//////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES) || defined(INCNS)
module ParticleClass
use SMConstants
! use NodalStorageClass
 use HexMeshClass
 use ElementClass
 use PhysicsStorage
! use MonitorDefinitions
! use ResidualsMonitorClass
! use StatisticsMonitor
! use ProbeClass
! use SurfaceMonitorClass
! use VolumeMonitorClass
implicit none
!
#include "Includes.h"

private
public  Particle_t
!
!  ******************************
!  Main particle class definition
!  ******************************
!  
type Particle_t
    type(HexMesh), pointer :: mesh 
    real(KIND=RP)   :: pos(3)
    real(KIND=RP)   :: pos_old(3)    
    real(KIND=RP)   :: vel(3)
    real(KIND=RP)   :: vel_old(3)
    real(KIND=RP)   :: temp
    real(KIND=RP)   :: temp_old
    real(KIND=RP)   :: fluidVel(3)
    real(KIND=RP)   :: fluidTemp
    ! Physical properties of the particles (independent)
    real(KIND=RP)  :: D            !Particle diameter
    real(KIND=RP)  :: rho          !Particle density   
    real(KIND=RP)  :: h            !Convective heat transfer coefficient (2 x k / D_p) (Assuming Nu = 2) 
    real(KIND=RP)  :: cv           !Particle specific heat        

    ! Impact parameters 
    real(KIND=RP)  :: rcoeff    ! Restitution coefficient
    real(KIND=RP)  :: delta     ! Roughness angle (sommerfeld model)

    ! Physical properties of the particles (dependent)
    real(KIND=RP)  :: m            !Particle mass ( rho_p x PI x D_p^3 / 6 )
    real(KIND=RP)  :: tau          !Aerodynamic response time ( rho_p x D_p^2 / 18 mu )
    logical                         :: active
    integer                         :: rank
    integer                         :: ID
    integer                         :: eID
    integer                         :: eID_old
    real(kind=RP)                   :: x(3)
    real(kind=RP)                   :: xi(3)
    real(kind=RP)                   :: xi_old(3)
    real(kind=RP), allocatable      :: lxi(:) , leta(:), lzeta(:)
    logical                         :: lost 
!    real(kind=RP)                   :: values(BUFFER_SIZE)    
!    character(len=STR_LEN_MONITORS) :: fileName
!    character(len=STR_LEN_MONITORS) :: monitorName
!    character(len=STR_LEN_MONITORS) :: variable
    contains
        procedure   :: init                => particle_init
        procedure   :: set_pos             => particle_set_pos
        procedure   :: set_vel             => particle_set_vel
        procedure   :: set_temp            => particle_set_temp    
        procedure   :: setGlobalPos        => particle_setGlobalPos
        procedure   :: getFluidVelandTemp  => particle_getFluidVelandTemp
        procedure   :: show                => particle_show
        procedure   :: integrate           => particle_integrate
        procedure   :: source              => particle_source
        procedure   :: updateVelRK3        => particle_updateVelRK3
        procedure   :: updateTempRK3       => particle_updateTempRK3        
end type Particle_t

real(kind=RP)  :: DEG2RAD                 = PI/180._RP
!
!  ========
contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_init(self, mesh)
   class(particle_t)       :: self
   type(HexMesh), target   :: mesh

    self % mesh        => mesh
    self % pos(:)      = 0.0_RP
    self % vel(:)      = 0.0_RP
    self % temp        = 0.0_RP
    self % fluidVel(:) = 0.0_RP
    self % fluidTemp   = 0.0_RP    
    self % active      = .true.
    self % eID         = -1
    self % D           = 1.0_RP
    self % rho         = 1.0_RP
    self % h           = 1.0_RP
    self % cv          = 1.0_RP
    self % m           = self % rho * PI * POW3(self % D) / 6
    self % lost        = .false. 
    self % rcoeff      = 1.0_RP
    self % delta       = 0.0_RP
    !self % tau         = self % rho * POW2(self % D) / (18 * mu) 
    ! DEPENDS ON THE FLUID VISCOSITY. IT SHOULD BE UPDATED ACCORDINGLY
end subroutine particle_init
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_set_pos(self, pos)
    implicit none
    class(particle_t)  :: self
    real(KIND=RP)      :: pos(3)

    self % pos(:) = pos

end subroutine particle_set_pos
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_set_vel(self, vel)
    implicit none
    class(particle_t)  :: self
    real(KIND=RP)      :: vel(3)

    self % vel(:) = vel

end subroutine particle_set_vel
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_set_temp(self, temp)
    implicit none
    class(particle_t)  :: self
    real(KIND=RP)      :: temp

    self % temp = temp

end subroutine particle_set_temp
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_show ( self ) 
    implicit none
    class(Particle_t), intent(inout)  :: self   

    write(*,*) "eID         = ", self % eID          
    write(*,*) "active      = ", self % active       
    write(*,*) "vel(:)      = ", self % vel(:)       
    write(*,*) "temp        = ", self % temp         
    write(*,*) "pos(:)      = ", self % pos(:)       
    write(*,*) "fluidVel(:) = ", self % fluidVel(:)  
    write(*,*) "fluidTemp   = ", self % fluidTemp        

end subroutine     

!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_setGlobalPos ( self, mesh )

    implicit none
    class(Particle_t)       , intent(inout)  :: self    
    class(HexMesh)          , intent(in)     :: mesh    
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    integer                 :: i ,j, k 
    integer, dimension(1)   :: elementwas 
    
    elementwas(1)  = self % eID
    self % eID_old = self % eID 
    self % xi_old  = self % xi 
    !
    ! Find the requested point in the mesh
    ! ------------------------------------
    self % active = &
    mesh % FindPointWithCoords(self % pos, & ! physical position of the particle
                               self % eID, & ! element in which the particle is
                               self % xi,  & ! computational position of the particle (rel to eID)
                               elementwas  & ! element in which the particle was
                               )

    if ( .not. self % active) return 
    !
    ! Check whether the probe is located in other partition
    ! -----------------------------------------------------
    !  call self % LookInOtherPartitions
    !
    ! Disable the probe if the point is not found
    ! -------------------------------------------
    !  if ( .not. self % active ) then
    !     if ( MPI_Process % isRoot ) then
    !        write(STD_OUT,'(A,I0,A)') "Probe ", ID, " was not successfully initialized."
    !        print*, "Probe is set to inactive."
    !     end if

    !     return
    !  end if    
    !
    ! Get the Lagrange interpolants
    ! -----------------------------

    associate(e => mesh % elements(self % eID))
        if ( elementwas(1) /= self % eID ) then 
            if (allocated(self % lxi)) deallocate(self % lxi)
            allocate( self % lxi(0 : e % Nxyz(1)) )
            if (allocated(self % leta)) deallocate(self % leta)
            allocate( self % leta(0 : e % Nxyz(2)) )
            if (allocated(self % lzeta)) deallocate(self % lzeta)        
            allocate( self % lzeta(0 : e % Nxyz(3)) )
            self % lxi = e % spAxi % lj(self % xi(1))
            self % leta = e % spAeta % lj(self % xi(2))
            self % lzeta = e % spAzeta % lj(self % xi(3))
        endif 
        ! 
        ! Recover the coordinates from direct projection. These will be the real coordinates
        ! ----------------------------------------------
        self % x = 0.0_RP
        do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            self % x = self % x + e % geom % x(:,i,j,k) * self % lxi(i) * self % leta(j) * self % lzeta(k)
        end do                  ; end do                ; end do      
    end associate   
#endif                      
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_getFluidVelandTemp ( self, mesh ) 
    use Physics
    use MPI_Process_Info
    use VariableConversion
    implicit none
    class(Particle_t), intent(inout)  :: self   
    class(HexMesh)                    :: mesh
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    integer        :: i, j, k, ierr
    integer        :: dir
    real(kind=RP)  :: value
    real(kind=RP)  :: vel(3,&
                          0:mesh % elements(self % eID) % Nxyz(1),&
                          0:mesh % elements(self % eID) % Nxyz(2),&
                          0:mesh % elements(self % eID) % Nxyz(3)  )
    real(kind=RP)  :: temp(0:mesh % elements(self % eID) % Nxyz(1),&
                           0:mesh % elements(self % eID) % Nxyz(2),&
                           0:mesh % elements(self % eID) % Nxyz(3)  )

    ! if ( MPI_Process % rank .eq. self % rank ) then
!
!           Update the probe
!           ----------------
       associate( e => mesh % elements(self % eID) )
       associate( Q => e % storage % Q )

        do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1) 
            ! get u velocity
            vel(1,i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k) 
            ! get v velocity
            vel(2,i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
            ! get w velocity            
            vel(3,i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
            ! get temp
            temp(i,j,k) = Temperature( Q(:,i,j,k) )
        end do            ; end do             ; end do

        self % fluidVel(:) = 0.0_RP
        do k = 0, e % Nxyz(3)    ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
            self % fluidVel(:) = self % fluidVel(:) + vel(:,i,j,k) * self % lxi(i) * self % leta(j) * self % lzeta(k)
        end do               ; end do             ; end do

        self % fluidTemp = 0.0_RP
        do k = 0, e % Nxyz(3)    ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
            self % fluidTemp = self % fluidTemp + temp(i,j,k) * self % lxi(i) * self % leta(j) * self % lzeta(k)
        end do               ; end do             ; end do            
       end associate
       end associate
! #ifdef _HAS_MPI_            
!        if ( MPI_Process % doMPIAction ) then
! !
! !              Share the result with the rest of the processes
! !              -----------------------------------------------         
!           call mpi_bcast(value, 1, MPI_DOUBLE, self % rank, MPI_COMM_WORLD, ierr)

!        end if
! #endif
!     else
! !
! !           Receive the result from the rank that contains the probe
! !           --------------------------------------------------------
! #ifdef _HAS_MPI_
!        if ( MPI_Process % doMPIAction ) then
!           call mpi_bcast(self % values(bufferPosition), 1, MPI_DOUBLE, self % rank, MPI_COMM_WORLD, ierr)
!        end if
! #endif
!     end if                    
#endif
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_integrate ( self, dt, St, Nu, phim, cvpdivcv, I0 ) 
#if defined(NAVIERSTOKES)
    use Physics,   only : sutherlandsLaw
    use FluidData, only : dimensionless, thermodynamics
#endif
    implicit none
    class(Particle_t)               , intent(inout)  :: self   
    real(KIND=RP)                   , intent(in)     :: dt 
    real(KIND=RP), intent(in) :: St         ! Particle non dimensional number
    real(KIND=RP), intent(in) :: Nu         ! Particle non dimensional number 
    real(KIND=RP), intent(in) :: phim       ! Particle non dimensional number 
    real(KIND=RP), intent(in) :: cvpdivcv   ! Particle non dimensional number 
    real(KIND=RP), intent(in) :: I0         ! Particle non dimensional number (radiation intensity)
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    real(KIND=RP) :: mu
    real(KIND=RP) :: invFr2
    real(KIND=RP) :: gamma
    real(KIND=RP) :: Pr
    real(KIND=RP) :: gravity(3) 
    logical       :: bounce 

    self % lost = .false. 

    if ( .not. self % active ) then 
      call compute_bounce_parameters(self, bounce)
    endif 

    mu              = SutherlandsLaw(self % fluidTemp)    ! Non dimensional viscosity mu(T)
    invFr2          = dimensionless  % invFr2             ! Fluid non dimensional number
    gravity         = dimensionless  % gravity_dir        ! Direction of gravity vector (normalised)    
    gamma           = thermodynamics % gamma              ! Fluid non dimensional number
    Pr              = dimensionless  % Pr                 ! Fluid non dimensional number

    ! St              = dimensionlessParticles % St        ! Particle non dimensional number
    ! Nu              = dimensionlessParticles % Nu        ! Particle non dimensional number 
    ! phim            = dimensionlessParticles % phim      ! Particle non dimensional number     
    ! cvpdivcv        = dimensionlessParticles % cvpdivcv  ! Particle non dimensional number

!
! RK3 method for now
! --------------------------- 
!TDG: esto no es realmente correcto. Hay dos aproximaciones intrínsecas en este modelo:
!   1.- Si avanza el tiempo, habría que actualizar el flujo (self % fluidVel). 
!   2.- Si avanza el tiempo, cambia la posición de la partícula, que habría que actualizar, y el flujo
!en ese punto (self % fluidVel).
!   Al hacerlo de esta forma, el código es más eficiente. Habría que cuantificar el error cometido.

    self % vel_old  = self % vel 
    self % pos_old  = self % pos 
    self % temp_old = self % temp 

    ! VELOCITY
    self % vel = self % updateVelRK3 ( dt, mu, St, invFr2, gravity  ) 
    ! POSITION
    self % pos = self % pos + dt * self % vel
    ! TEMPERATURE
    self % temp = self % updateTempRK3 ( dt, gamma, cvpdivcv, St, Pr, I0, Nu ) 

    !print*,  self % pos, "//", self % vel, "//", self % temp 
#endif
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine compute_bounce_parameters(p, bounce)
    class(particle_t)   , intent(inout)      :: p
    
    real(kind=RP)        :: GEOM_TOL = 1e-5_RP
    integer              :: i, j, k, fdir, face, eID, faceID 
    integer, parameter   :: ploc(3)=[4,2,5] ! Face elemente mapping
    integer, parameter   :: nloc(3)=[6,1,3]
    integer, parameter   :: fdirmap(6) = [2,2,3,1,3,1]
    real(kind=RP)        :: pos(3), pos_out(3), pos_in(3), xi_init(3), xi(3)
    real(kind=RP)        :: normal(3),v_in(3), d_in(3), d_out(3), f_xi(2)
    real(kind=RP)        :: dt_after_collision, normal_rough(3)
    real(kind=RP)        :: CollisionPlaneNormal(3), ImpactAngle,RoughAngle
    logical              :: inside, bounce, fine_search
    integer, parameter   :: max_iter = 50
    integer              :: Nxyz(3), neighbours(6)
    real(kind=RP),allocatable  :: lxi(:), leta(:), lzeta(:)


    eID         = p%eID_old
    pos_out     = p%pos
    xi_init     = p%xi_old
    bounce      = .false.
    fine_search = .false.
    
    ! Search the face intersected by the trace
    ! -----

    !inside = p%mesh%elements(eID)%FindPointWithCoords(pos_out, xi, xi_init)
    fdir = maxloc (abs(xi_init), dim=1)

    if (xi_init(fdir) > 0._RP) then
       face = ploc(fdir)
    else
       face = nloc(fdir)
    end if

    ! Check face boundary type
    if ( (trim(p%mesh%elements(eID)%boundaryName(face)) == "wall") .or. &
        & (trim(p%mesh%elements(eID)%boundaryName(face)) == "pipe") ) then
       bounce = .true.
       fine_search = .false.
    else if (trim(p%mesh%elements(eID)%boundaryName(face)) == "---") then
       ! Particle is closer to an interior face than a boundary, 
       ! Collision may occur in other element
       bounce = .true.
       fine_search = .true.

    else
       ! The particle has reached a permeable boundary
       bounce = .false.
    end if

    if (.not. bounce) return


    ! Compute intersection between particle trace and face
    ! ----
    
    ! The fast way
    if (.not. fine_search) then

       pos_in   = p%pos_old
       pos_out  = p%pos
       xi_init  = 0._RP
       
       ! Try locationg the collision pont using a bisection method
       do i = 1, max_iter
          pos = 0.5_RP * (pos_in + pos_out)
          !inside = p%mesh%elements(eID)%FindPointWithCoords(pos, xi, xi_init)
          inside = p%mesh%elements(eID)%FindPointWithCoords(pos, xi)
       
          if (inside) then
             ! Check tolerance in the face direction
             if (abs(abs(xi(fdir))-1._RP) < GEOM_TOL) then
                exit
             else
                pos_in = pos
                xi_init = xi
             end if
          else
             pos_out = pos
          end if
       end do
       ! If collision point has not been found, activate the fine search 
       if (i >= max_iter)  fine_search = .true.
    end if


    
    ! Try to find the collision point using a more expensive method in case
    ! the fast method fails

    if (fine_search) then
       
       pos_in   = p%pos_old
       pos_out  = p%pos
       xi_init  = p%xi_old
       
       ! Find the neighbours of the element in which the particle was

       neighbours = -1
       do j = 1, 6
          if (p%mesh%elements(eID)%NumberOfConnections(j) > 0) then
                neighbours(j) = p%mesh%elements(eID)%Connection(j)%ElementIDs(1) 
          else 
                neighbours(j) = -1
          end if
       end do

       iter_loop:do i = 1, max_iter
          
          pos = 0.5_RP * (pos_in + pos_out)
          
          ! Include neigbours when searching the point coordinates
          call FindPointWithCoords_Neighbours(p, pos,eID,neighbours, xi,inside)
          
          if (inside) then

             if  (any( abs(xi) > 1._RP - GEOM_TOL) ) then
                exit iter_loop
             end if
             

             pos_in = pos
          else
             pos_out = pos
          end if
          
          if ((i >= max_iter))  then
             bounce = .false.
             !Write (*,*) "Warning, particle lost"
             p%lost = .true.
             return
          end if
       end do iter_loop

    
       ! Re check the collision face
       fdir = maxloc (abs(xi), dim=1)
       if (xi(fdir) > 0._RP) then
          face = ploc(fdir)
       else
          face = nloc(fdir)
       end if
       
       if ( (trim(p%mesh%elements(eID)%boundaryName(face)) == "wall") .or. &
        & (trim(p%mesh%elements(eID)%boundaryName(face)) == "pipe") ) then
          bounce = .true.
       
       else if (trim(p%mesh%elements(eID)%boundaryName(face)) == "---") then
          ! Particle has left the domain through an intertal face, the
          !  particle is lost
          bounce = .false.
          !Write (*,*) "Warning, particle lost after crossing internal face"
          p%lost = .true.
          return
       
       else
          ! The particle has reached a permeable boundary
          bounce = .false.
          return
       end if
    
    end if
    
    ! Get face normal at collision point
    ! -----
    
    faceID = p%mesh%elements(eID)%faceIDs(face)
    
    Nxyz        = p % mesh % elements(eID) % Nxyz
    allocate(lxi(0:Nxyz(1)))
    allocate(leta(0:Nxyz(2)))
    allocate(lzeta(0:Nxyz(3)))
    lxi(0:)     = p % mesh % elements(eID) % spAxi % lj   (xi(1))
    leta(0:)     = p % mesh % elements(eID) % spAeta % lj  (xi(2))
    lzeta(0:)   = p % mesh % elements(eID) % spAzeta % lj (xi(3))
    
    normal = 0.0_RP
       select case (face)
          case (1,2)
             do k = 0,Nxyz(3)  ; do i = 0, Nxyz(1)
                normal = normal + p%mesh%faces(faceID)%geom%normal(:,i,k)  * lxi(i) * lzeta(k)
             end do            ; end do
             f_xi = [xi(1),xi(3)]
          case (3,5)
             do j = 0,Nxyz(2)  ; do i = 0, Nxyz(1)
                normal = normal + p%mesh%faces(faceID)%geom%normal(:,i,j) * lxi(i) * leta(j)
             end do            ; end do
             f_xi = [xi(1),xi(2)]
          case (4,6)
             do j = 0,Nxyz(2)  ; do k = 0, Nxyz(3)
                normal = normal + p%mesh%faces(faceID)%geom%normal(:,j,k) * leta(j) * lzeta(k)
             end do            ; end do
             f_xi = [xi(2),xi(3)]
       end select
    
    normal = - normal
    

    ! Compute velocity and position after collsion
    ! -----

!     ! Out direction
!     print*, 595
! !    print*, ImpactAngle
!     print*, p%delta 
!     read(*,*)
    d_in  = -p%vel/norm2(p%vel)
    normal = normal/norm2(normal)

    CollisionPlaneNormal = cross(normal,d_in)
    CollisionPlaneNormal = CollisionPlaneNormal / norm2(CollisionPlaneNormal)
    

   !  ! Modify Colision plane using Gaussian Noise
   !  RoughAngle = GaussianNoise() * p%delta
   !  CollisionPlaneNormal = rotate_vector(CollisionPlaneNormal, d_in, RoughAngle)
    
   !  ! Modify Normal using Surface Rough Angle
   !  print*, ImpactAngle
   !  print*, p%delta 
   !  read(*,*)
   !  RoughAngle = WallRoughnessAngleFast(ImpactAngle, p%delta)
   !  normal_rough = rotate_vector(normal, CollisionPlaneNormal, p%delta)
    normal_rough = normal    
    ! Specular rebound with modified normal
    d_out = 2._RP * dot_product(normal_rough,d_in)*normal_rough - d_in
   
    ImpactAngle = PI/2._RP - acos(dot_product(normal_rough,d_in))

    ! Velocity vector
    ! -----
    v_in = p%vel
    p%vel = norm2(v_in) * d_out/norm2(d_out) * p%rcoeff
    
    ! New position
    ! -----
    dt_after_collision = norm2(p%pos-p%pos_old)/norm2(p%vel)

    p%pos = pos + p%vel * dt_after_collision
    
    ! Set element and particle as active  
    ! -----

    p%active = .true.
    p%eID = eID

    ! store collision parameters
    ! -----
    ! if (size(p%collisions,dim=1) == p%ncollisions ) then
    !    call p%increase_collisions_size
    ! end if
    
    ! p%ncollisions = p%ncollisions + 1
    ! p%collisions(p%ncollisions)%pos           = pos
    ! p%collisions(p%ncollisions)%normal        = normal
    ! p%collisions(p%ncollisions)%vel_in        = v_in * refValues%V
    ! p%collisions(p%ncollisions)%vel_out       = p%vel * refValues%V
    ! p%collisions(p%ncollisions)%eID           = eID
    ! p%collisions(p%ncollisions)%faceID        = faceID
    ! p%collisions(p%ncollisions)%ImpactAngle   = ImpactAngle * 180._RP / PI
    ! p%collisions(p%ncollisions)%xi            = f_xi

 end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////
!
 subroutine  FindPointWithCoords_Neighbours(p, pos,eID,neighbours, xi,inside)
   class(Particle_t), intent(in)      :: p
   real(kind=RP)  , intent(in)      :: pos(3)
   integer        , intent(inout)   :: eID
   integer        , intent(in)      :: neighbours(6)
   real(kind=RP)  , intent(out)     :: xi(3)
   logical        , intent(out)     :: inside

   integer  :: i
   
   inside = p%mesh%elements(eID)%FindPointWithCoords(pos, xi)
   if (inside) return

   ! Else check if the point resides iniside a neigbour

   do i = 1, 6
      if (neighbours(i) <= 0) cycle
      inside = p%mesh%elements(neighbours(i))%FindPointWithCoords(pos, xi)
      if (inside) then
         eID = neighbours(i)
         return
      end if
   end do

   inside = .false.


end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////
!
function particle_updateVelRK3 (self, dt, mu, St, invFr2, gravity) result(Q)
    implicit none 

    class(Particle_t), intent(in)     :: self 
    real(KIND=RP),     intent(in)     :: dt
    real(KIND=RP),     intent(in)     :: mu
    real(KIND=RP),     intent(in)     :: St
    real(KIND=RP),     intent(in)     :: invFr2
    real(KIND=RP),     intent(in)     :: gravity(3)
    REAL(KIND=RP)                     :: Q(3)
#if defined(NAVIERSTOKES)
!
!   ---------------
!   Local variables
!   ---------------
!
    INTEGER        :: i, j, k    
    REAL(KIND=RP), DIMENSION(3) :: G, Qdot
    REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
    REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
    REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)   

    G = 0.0_RP
    Q = self % vel
    DO k = 1,3
        Qdot = mu / St * ( self % fluidVel - self % vel) + invFr2 * gravity
        G = a(k) * G  + Qdot
        Q = Q         + c(k) * dt * G
    END DO
#endif
end function 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
function particle_updateTempRK3 (self, dt, gamma, cvpdivcv, St, Pr, I0, Nu ) result(Q)
        implicit none 
    
        class(Particle_t), intent(in)     :: self 
        real(KIND=RP),     intent(in)     :: dt
        real(KIND=RP),     intent(in)     :: gamma
        real(KIND=RP),     intent(in)     :: cvpdivcv
        real(KIND=RP),     intent(in)     :: St
        real(KIND=RP),     intent(in)     :: Pr
        real(KIND=RP),     intent(in)     :: I0
        real(KIND=RP),     intent(in)     :: Nu
        real(KIND=RP)                     :: Q
#if defined(NAVIERSTOKES)    
    !
    !   ---------------
    !   Local variables
    !   ---------------
    !
        INTEGER        :: i, j, k    
        REAL(KIND=RP)               :: G, Qdot
        REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
        REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
        REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)   
    
        !TDG: convert gamma / (3 * cvpdivcv * St * Pr) into a variable
        G    = 0.0_RP
        Q = self % temp
        DO k = 1,3
            Qdot = gamma / (3 * cvpdivcv * St * Pr) * &
                ( I0 - Nu * ( self % temp - self % fluidTemp ) )
            G = a(k) * G  + Qdot
            Q = Q         + c(k) * dt * G
        END DO
#endif    
    end function 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_source ( self, e, source )
    implicit none
    class(Particle_t)       , intent(in)    :: self    
    class(element)          , intent(in)    :: e    
    real(kind=RP)           , intent(inout) :: source(:,0:,0:,0:)    
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    integer         :: i ,j, k                     
    real(kind=RP)   :: delta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3) ) 
    real(kind=RP)   :: x(3) 
    real(kind=RP)   :: val 
    real(kind=RP)   :: dx 

    if ( .not. self % active ) return 

    !**************************************
    ! COMPUTE SCALING OF THE ELEMENT (dx) 
    !**************************************
    val = 0.0_RP
    do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
      val = val + e % spAxi % w(i) * e % spAeta % w(j) * e % spAzeta % w(k) * e % geom % jacobian(i,j,k)
    end do            ; end do           ; end do
    dx = val ** (1.0_RP/3.0_RP)
    !*********************************
    ! CONSTRUCT DISCRETE DIRAC DELTA
    !*********************************
    do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
      x = e % geom % x(:,i,j,k)
      !Compute value of approximation of delta Diract (exp)
      delta(i,j,k) = deltaDirac( x(1) - self % x(1), &
                          x(2) - self % x(2), &
                          x(3) - self % x(3), &
                          dx )
    enddo ; enddo ; enddo    
    !*********************************************************************************
    ! DISCRETE NORMALIZATION OF DIRAC DELTA. DISCRETE INTEGRAL IN ELEMENT EQUAL TO 1 
    !*********************************************************************************
    val = 0.0_RP
    do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
      val = val + e % spAxi % w(i) * e % spAeta % w(j) * e % spAzeta % w(k) * e % geom % jacobian(i,j,k) * delta(i,j,k)
    end do            ; end do           ; end do
    delta = delta / val 
    !*****************************
    ! COMPUTATION OF SOURCE TERM 
    !*****************************
    do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)              
      source(2,i,j,k) = source(2,i,j,k) - ( self % fluidVel(1) - self % Vel(1) ) * delta(i,j,k)   
      source(3,i,j,k) = source(3,i,j,k) - ( self % fluidVel(2) - self % Vel(2) ) * delta(i,j,k)
      source(4,i,j,k) = source(4,i,j,k) - ( self % fluidVel(3) - self % Vel(3) ) * delta(i,j,k)
      source(5,i,j,k) = source(5,i,j,k) - ( self % fluidTemp   - self % temp   ) * delta(i,j,k) 
    enddo ; enddo ; enddo     
#endif     
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
function deltaDirac(x,y,z,dx)
    !This function creates an approximation of a delta Dirac
    !   centered at x,y,z. 
    implicit none 
    real(kind=RP)            :: deltaDirac
    real(kind=RP)            :: x,y,z
    real(kind=RP)            :: dx

    real(kind=RP)            :: sigma 

    !sigma = 1.5 * D
    sigma = dx / 10 
    ! This value has been taken from Maxey, Patel, Chang and Wang 1997
    ! According to them, sigma > 1.5 grid spacing AND sigma > 1.3 D for stability
    ! They use a value of sigma ~ D 

    ! This approximation has the same triple integral as the exact 
    ! Dirac delta

    deltaDirac = ( 2 * pi * POW2( sigma ) ) ** ( -3.0_RP / 2.0_RP ) * &
        exp( - ( POW2(x) + POW2(y) + POW2(z) ) / ( 2 * POW2( sigma ) ) )

end function 
!
!///////////////////////////////////////////////////////////////////////////////////////
!    
!
!///////////////////////////////////////////////////////////////////////////////////////
!
function WallRoughnessAngleFast(alpha1, delta)
    real(kind=RP), intent(in)  :: alpha1
    real(kind=RP), intent(in)  :: delta

    real(kind=RP)  :: WallRoughnessAngleFast

    real(kind=RP)  :: g, rnd(2), rnd_normal
    integer        :: i

    

    do i = 1, 100 ! 
       call random_number(rnd)
       rnd_normal = sqrt(-2._RP * log(rnd(1))) * cos(2._RP * PI * rnd(2)) !Miller-Box formula
       g = rnd_normal* delta
       if ((alpha1 + g) > 1e-3_RP )  then
          WallRoughnessAngleFast = g
          return
       end if
    end do

    ! Just in case...
    WallRoughnessAngleFast = 0._RP

 end function
!--
!--
 function GaussianNoise()
    real(kind=RP)  :: GaussianNoise 
    real(kind=RP)  :: rnd(2)
    call random_number(rnd)
    GaussianNoise = sqrt(-2._RP * log(rnd(1))) * cos(2._RP * PI * rnd(2)) !Miller-Box formula

 end function 
!--
!--
 function cross(a,b) result(c)
    real(kind=RP), dimension(3), intent(in)  :: a,b
    real(kind=RP), dimension(3)               :: c
       c(1) = a(2) * b(3) - a(3) * b(2)
       c(2) = a(3) * b(1) - a(1) * b(3)
       c(3) = a(1) * b(2) - a(2) * b(1)
 end function cross
!--
!--
 function rotate_vector(vector, rotation_vector, angle)
    real(kind=RP), dimension(3), intent(in)   :: vector,rotation_vector
    real(kind=RP)              , intent(in)   :: angle       
    real(kind=RP), dimension(3)               :: rotate_vector

    real(kind=RP)  :: c, s, onemc,  u, v, w, rmat(3,3), rot_norm
    
    c = cos(angle*DEG2RAD)
    s = sin(angle*DEG2RAD)
    onemc = 1._RP - c

    ! Normalize rotation vector
    rot_norm = norm2(rotation_vector)
    u = rotation_vector(1)/rot_norm
    v = rotation_vector(2)/rot_norm
    w = rotation_vector(3)/rot_norm

    ! Rotation matrix
    rmat(1,:) = [c+u*u*onemc  , u*v*onemc-w*s, u*w*onemc+v*s]
    rmat(2,:) = [u*v*onemc+w*s, c+v*v*onemc  , v*w*onemc-u*s]
    rmat(3,:) = [u*w*onemc-v*s, v*w*onemc+u*s, c+w*w*onemc  ]

    rotate_vector = matmul(rmat,vector)


 end function
end module 
#endif
