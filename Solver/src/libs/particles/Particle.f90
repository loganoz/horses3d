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
    real(KIND=RP)   :: pos(3)
    real(KIND=RP)   :: vel(3)
    real(KIND=RP)   :: temp
    real(KIND=RP)   :: fluidVel(3)
    real(KIND=RP)   :: fluidTemp
    ! Physical properties of the particles (independent)
    real(KIND=RP)  :: D            !Particle diameter
    real(KIND=RP)  :: rho          !Particle density   
    real(KIND=RP)  :: h            !Convective heat transfer coefficient (2 x k / D_p) (Assuming Nu = 2) 
    real(KIND=RP)  :: cv           !Particle specific heat        

    ! Physical properties of the particles (dependent)
    real(KIND=RP)  :: m            !Particle mass ( rho_p x PI x D_p^3 / 6 )
    real(KIND=RP)  :: tau          !Aerodynamic response time ( rho_p x D_p^2 / 18 mu )
    logical                         :: active
    integer                         :: rank
    integer                         :: ID
    integer                         :: eID
    real(kind=RP)                   :: x(3)
    real(kind=RP)                   :: xi(3)
    real(kind=RP), allocatable      :: lxi(:) , leta(:), lzeta(:)
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
!
!  ========
contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_init(self)
    implicit none
    class(particle_t)  :: self

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
    
    elementwas(1) = self % eID
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

    if ( .not. self % active ) return 

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
subroutine particle_integrate ( self, dt, St, Nu, phim, cvpdivcv, I0, gravity ) 
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
    real(KIND=RP), intent(in) :: gravity(3) ! Particle non dimensional vector (gravity)
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


    if ( .not. self % active ) return 

    mu              = SutherlandsLaw(self % fluidTemp)    ! Non dimensional viscosity mu(T)
    invFr2 = dimensionless  % invFr2    ! Fluid non dimensional number
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

    ! VELOCITY
    self % vel = self % updateVelRK3 ( dt, mu, St, invFr2, gravity  ) 
    ! POSITION
    self % pos = self % pos + dt * self % vel
    ! TEMPERATURE
    self % temp = self % updateTempRK3 ( dt, gamma, cvpdivcv, St, Pr, I0, Nu ) 

    print*,  self % pos, "//", self % vel, "//", self % temp 
#endif
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
        Qdot = mu / St * ( self % fluidVel - self % vel) - invFr2 * gravity
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
    real(kind=RP)   :: delta 

    if ( .not. self % active ) return 

    do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
        !Compute value of approximation of delta Diract (exp)
        delta = deltaDirac( e % spAxi % x(i)      - self % xi(1), &
                            e % spAeta % x(j)     - self % xi(2), &
                            e % spAzeta % x(k)    - self % xi(3) )

        !Compute source term                    
        source(2,i,j,k) = source(2,i,j,k) - ( self % fluidVel(1) - self % Vel(1) ) * delta   
        source(3,i,j,k) = source(3,i,j,k) - ( self % fluidVel(2) - self % Vel(2) ) * delta
        source(4,i,j,k) = source(4,i,j,k) - ( self % fluidVel(3) - self % Vel(3) ) * delta
        source(5,i,j,k) = source(5,i,j,k) - ( self % fluidTemp   - self % temp   ) * delta 
    enddo ; enddo ; enddo     
#endif     
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
function deltaDirac(x,y,z)
    !This function creates an approximation of a delta Dirac
    !   centered at x,y,z. 
    !   c -> 0 == delta dirac
    !   The value of the parameter c should be adapted depending 
    !   on the polynomial order used in the approximation. 
    !TDG: this constant c should depend on the pol order or the diameter of the particle.
    implicit none 
    real(kind=RP)            :: deltaDirac
    real(kind=RP)            :: x,y,z
    real(kind=RP), parameter :: c = 0.02_RP

    deltaDirac = exp( - (POW2(x) + POW2(y) + POW2(z)) / c )

end function 
!
!///////////////////////////////////////////////////////////////////////////////////////
!    
end module 
#endif
