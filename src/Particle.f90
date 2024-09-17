#include "Includes.h"
#ifdef FLOW
module ParticleClass
 use SMConstants
 use HexMeshClass
 use ElementClass
 use PhysicsStorage
 use NodalStorageClass, only: NodalStorage
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
    real(KIND=RP)   :: pos_old(3)    
    real(KIND=RP)   :: vel(3)
    real(KIND=RP)   :: vel_old(3)
    real(KIND=RP)   :: temp
    real(KIND=RP)   :: temp_old
    real(KIND=RP)   :: fluidVel(3)
    real(KIND=RP)   :: fluidTemp

    logical         :: active
    integer                         :: eID
    integer                         :: eID_old
    real(kind=RP)                   :: x(3)
    real(kind=RP)                   :: xi(3)
    real(kind=RP)                   :: xi_old(3)
    real(kind=RP), allocatable      :: lxi(:), leta(:), lzeta(:)

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
subroutine particle_init(self, mesh)
   class(particle_t)       :: self
   type(HexMesh), target   :: mesh

    self % pos(:)      = 0.0_RP
    self % vel(:)      = 0.0_RP
    self % temp        = 0.0_RP
    self % fluidVel(:) = 0.0_RP
    self % fluidTemp   = 0.0_RP    
    self % active      = .true.
    self % eID         = -1

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
    associate(spAxi   => NodalStorage(e % Nxyz(1)), &
              spAeta  => NodalStorage(e % Nxyz(2)), &
              spAzeta => NodalStorage(e % Nxyz(3)) )
        if ( elementwas(1) /= self % eID ) then 
            if (allocated(self % lxi)) deallocate(self % lxi)
            allocate( self % lxi(0 : e % Nxyz(1)) )
            if (allocated(self % leta)) deallocate(self % leta)
            allocate( self % leta(0 : e % Nxyz(2)) )
            if (allocated(self % lzeta)) deallocate(self % lzeta)        
            allocate( self % lzeta(0 : e % Nxyz(3)) )
            self % lxi = spAxi % lj(self % xi(1))
            self % leta = spAeta % lj(self % xi(2))
            self % lzeta = spAzeta % lj(self % xi(3))
        endif 
        ! 
        ! Recover the coordinates from direct projection. These will be the real coordinates
        ! ----------------------------------------------
        self % x = 0.0_RP
        do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            self % x = self % x + e % geom % x(:,i,j,k) * self % lxi(i) * self % leta(j) * self % lzeta(k)
        end do                  ; end do                ; end do      
    end associate   
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
subroutine particle_integrate ( self, mesh, dt, St, Nu, phim, I0, gammaDiv3cvpdivcvStPr, minbox, maxbox , bcbox ) 
#if defined(NAVIERSTOKES)
    use VariableConversion,   only : sutherlandsLaw
    use FluidData, only : dimensionless
#endif
    implicit none
    
    class(Particle_t)               , intent(inout)  :: self   
    type(HexMesh)                   , intent(in)     :: mesh
    real(KIND=RP)                   , intent(in)     :: dt 
    real(KIND=RP), intent(in) :: St         ! Particle non dimensional number
    real(KIND=RP), intent(in) :: Nu         ! Particle non dimensional number 
    real(KIND=RP), intent(in) :: phim       ! Particle non dimensional number 
    real(KIND=RP), intent(in) :: I0         ! Particle non dimensional number (radiation intensity)
    real(KIND=RP), intent(in) :: gammaDiv3cvpdivcvStPr ! Particle and fluid comb of properties to increase performance
    real(KIND=RP), intent(in)       :: minbox(3)  ! Minimum value of box for particles    
    real(KIND=RP), intent(in)       :: maxbox(3)  ! Maximum value of box for particles
    integer, intent(in)       :: bcbox(3)   ! Boundary conditions of box for particles [0,1,2] [inflow/outflow, wall, periodic]      

#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    real(KIND=RP) :: mu
    real(KIND=RP) :: invFr2
    real(KIND=RP) :: gravity(3) 

    if ( .not. self % active ) then 
      call compute_bounce_parameters(self, mesh, minbox, maxbox, bcbox)
    endif 

    mu              = SutherlandsLaw(self % fluidTemp)    ! Non dimensional viscosity mu(T)
    invFr2          = dimensionless  % invFr2             ! Fluid non dimensional number
    gravity         = dimensionless  % gravity_dir        ! Direction of gravity vector (normalised)    

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
    self % temp = self % updateTempRK3 ( dt, I0, Nu, gammaDiv3cvpdivcvStPr ) 

#endif
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine compute_bounce_parameters(p, mesh, minbox, maxbox, bcbox)
   class(particle_t)   , intent(inout)      :: p
   type(HexMesh), intent(in)       :: mesh
   real(KIND=RP), intent(in)       :: minbox(3)  ! Minimum value of box for particles    
   real(KIND=RP), intent(in)       :: maxbox(3)  ! Maximum value of box for particles
   integer, intent(in)       :: bcbox(3)   ! Boundary conditions of box for particles [0,1,2] [inflow/outflow, wall, periodic] 


   integer              :: eID
   real(kind=RP)        :: pos_out(3), pos_in(3)

   eID         = p%eID_old
   pos_in      = p%pos_old
   pos_out     = p%pos
   
      ! Check exit face of the box [i+-,j+-,k+-]

         ! bcbox[yz,xz,xy] 0 is inflow/outflow, 1 is wall, 2 is periodic  
         if ( p % pos(1) < minbox(1) ) then 
            if     ( bcbox(1) == 0 ) then 
               ! Particle abandoned the domain through inflow or outflow
               ! Reinject at outflow or inflow
               ! call p % set_vel  ( v )
               ! call p % set_temp ( T )
               ! p % pos(1) = p % pos(1) - ( minbox(1) - maxbox(1) )
            elseif ( bcbox(1) == 1 ) then 
               !print*, "Warning. Particle lost through wall."
               !pos_in(1) 
               p % pos(1) = minbox(1) - ( pos_out(1) - minbox(1) )  
               p % vel(1)   = - p % vel(1)
            elseif ( bcbox(1) == 2 ) then  
               p % pos(1) = p % pos(1) - ( minbox(1) - maxbox(1) )
            endif
         endif 

         if ( p % pos(2) < minbox(2) ) then 
            if     ( bcbox(2) == 0 ) then 
               ! Particle abandoned the domain through inflow or outflow
               ! Reinject at outflow or inflow
               ! call p % set_vel  ( v )
               ! call p % set_temp ( T )
               ! p % pos(2) = p % pos(2) - ( minbox(2) - maxbox(2) )
            elseif ( bcbox(2) == 1 ) then 
               !print*, "Warning. Particle lost through wall."   
               p % pos(2) = minbox(2) - ( pos_out(2) - minbox(2) )  
               p % vel(2)   = - p % vel(2)
            elseif ( bcbox(2) == 2 ) then  
               p % pos(2) = p % pos(2) - ( minbox(2) - maxbox(2) )
            endif
         endif 

         if ( p % pos(3) < minbox(3) ) then 
            if     ( bcbox(3) == 0 ) then 
               ! Particle abandoned the domain through inflow or outflow
               ! Reinject at outflow or inflow
               ! call p % set_vel  ( v )
               ! call p % set_temp ( T )
               ! p % pos(3) = p % pos(3) - ( minbox(3) - maxbox(3) )
            elseif ( bcbox(3) == 1 ) then 
               !print*, "Warning. Particle lost through wall."
               p % pos(3) = minbox(3) - ( pos_out(3) - minbox(3) )  
               p % vel(3)   = - p % vel(3)
            elseif ( bcbox(3) == 2 ) then  
               p % pos(3) = p % pos(3) - ( minbox(3) - maxbox(3) )
            endif
         endif 

         if ( p % pos(1) > maxbox(1) ) then 
            if     ( bcbox(1) == 0 ) then 
               ! Particle abandoned the domain through inflow or outflow
               ! Reinject at outflow or inflow
               ! call p % set_vel  ( v )
               ! call p % set_temp ( T )
               ! p % pos(1) = p % pos(1) + ( minbox(1) - maxbox(1) ) 
            elseif ( bcbox(1) == 1 ) then 
               !print*, "Warning. Particle lost through wall."
               p % pos(1) = maxbox(1) - ( pos_out(1) - maxbox(1) )  
               p % vel(1)   = - p % vel(1)
            elseif ( bcbox(1) == 2 ) then  
               p % pos(1) = p % pos(1) + ( minbox(1) - maxbox(1) ) 
            endif
         endif 

         if ( p % pos(2) > maxbox(2) ) then 
            if     ( bcbox(2) == 0 ) then 
               ! Particle abandoned the domain through inflow or outflow
               ! Reinject at outflow or inflow
               ! call p % set_vel  ( v )
               ! call p % set_temp ( T )
               ! p % pos(2) = p % pos(2) + ( minbox(2) - maxbox(2) )
            elseif ( bcbox(2) == 1 ) then 
               !print*, "Warning. Particle lost through wall."
               p % pos(2) = maxbox(2) - ( pos_out(2) - maxbox(2) )  
               p % vel(2)   = - p % vel(2)
            elseif ( bcbox(2) == 2 ) then  
               p % pos(2) = p % pos(2) + ( minbox(2) - maxbox(2) ) 
            endif
         endif 

         if ( p % pos(3) > maxbox(3) ) then 
            if     ( bcbox(3) == 0 ) then 
               ! Particle abandoned the domain through inflow or outflow
               ! Reinject at outflow or inflow
               ! call p % set_vel  ( v )
               ! call p % set_temp ( T )
               ! p % pos(3) = p % pos(3) + ( minbox(3) - maxbox(3) )
            elseif ( bcbox(3) == 1 ) then 
               !print*, "Warning. Particle lost through wall."
               p % pos(3) = maxbox(3) - ( pos_out(3) - maxbox(3) )  
               p % vel(3)   = - p % vel(3)        
            elseif ( bcbox(3) == 2 ) then  
               p % pos(3) = p % pos(3) + ( minbox(3) - maxbox(3) )
            endif
         endif 

         p % eID = p%eID_old 
         call p % setGlobalPos ( mesh )

 end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////
!
 subroutine  FindPointWithCoords_Neighbours(p, mesh, pos,eID,neighbours, xi,inside)
   class(Particle_t), intent(in)      :: p
   class(HexMesh)   , intent(in)      :: mesh    
   real(kind=RP)  , intent(in)      :: pos(3)
   integer        , intent(inout)   :: eID
   integer        , intent(in)      :: neighbours(6)
   real(kind=RP)  , intent(out)     :: xi(3)
   logical        , intent(out)     :: inside

   integer  :: i
   
   inside = mesh%elements(eID)%FindPointWithCoords(pos, mesh % dir2D, xi)
   if (inside) return

   ! Else check if the point resides inside a neighbour

   do i = 1, 6
      if (neighbours(i) <= 0) cycle
      inside = mesh%elements(neighbours(i))%FindPointWithCoords(pos, mesh % dir2D, xi)
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
function particle_updateTempRK3 (self, dt, I0, Nu, gammaDiv3cvpdivcvStPr ) result(Q)
        implicit none
    
        class(Particle_t), intent(in)     :: self 
        real(KIND=RP),     intent(in)     :: dt
        real(KIND=RP),     intent(in)     :: I0
        real(KIND=RP),     intent(in)     :: Nu
        real(KIND=RP),     intent(in)     :: gammaDiv3cvpdivcvStPr
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

        G    = 0.0_RP
        Q = self % temp
        DO k = 1,3
            Qdot = gammaDiv3cvpdivcvStPr * &
                ( I0 - Nu * ( self % temp - self % fluidTemp ) )
            G = a(k) * G  + Qdot
            Q = Q         + c(k) * dt * G
        END DO
#endif    
    end function 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_source ( self, e, source, highordersource )
    implicit none
    class(Particle_t)       , intent(in)    :: self    
    class(element)          , intent(in)    :: e    
    real(kind=RP)           , intent(inout) :: source(:,0:,0:,0:)    
    logical                                 :: highordersource
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    integer         :: i ,j, k                     
    real(kind=RP)   :: vol  

    associate(spAxi   => NodalStorage(e % Nxyz(1)), &
      spAeta  => NodalStorage(e % Nxyz(2)), &
      spAzeta => NodalStorage(e % Nxyz(3)) )

    if ( .not. self % active ) return 

    !**************************************
    ! COMPUTE SCALING OF THE ELEMENT (vol) 
    !**************************************

    ! This information could be stored in the element so it does not have to be recomputed
    !vol = 0.0_RP
    !do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
    !  vol = vol + e % spAxi % w(i) * e % spAeta % w(j) * e % spAzeta % w(k) * e % geom % jacobian(i,j,k)
    !end do            ; end do           ; end do
    
    ! This is the same as the lines commented at top. Update this for particle_source_gaussian if I recover it.
    vol = e % geom % volume

    !*****************************
    ! COMPUTATION OF SOURCE TERM 
    !*****************************
    if (highordersource) then 
      do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)  
         ! New implementation - Kopriva's idea of dealing with dirac function in weak form 
         ! This can be optimized by precomputing:
         ! self % lxi(i) * self % leta(j) * self % lzeta(k) / ( e % geom % jacobian(i,j,k) * spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) )          
         source(2,i,j,k) = source(2,i,j,k) - ( self % fluidVel(1) - self % Vel(1) ) * self % lxi(i) * self % leta(j) * self % lzeta(k) &
         / ( e % geom % jacobian(i,j,k) * spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) )
         source(3,i,j,k) = source(3,i,j,k) - ( self % fluidVel(2) - self % Vel(2) ) * self % lxi(i) * self % leta(j) * self % lzeta(k) &
         / ( e % geom % jacobian(i,j,k) * spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) ) 
         source(4,i,j,k) = source(4,i,j,k) - ( self % fluidVel(3) - self % Vel(3) ) * self % lxi(i) * self % leta(j) * self % lzeta(k) &
         / ( e % geom % jacobian(i,j,k) * spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) ) 
         source(5,i,j,k) = source(5,i,j,k) - ( self % fluidTemp   - self % temp   ) * self % lxi(i) * self % leta(j) * self % lzeta(k) &
         / ( e % geom % jacobian(i,j,k) * spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) ) 
      enddo ; enddo ; enddo     
    else
      do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)  
         ! Old implementation - Low order approximation of the source term inside the element (constant)
         source(2,i,j,k) = source(2,i,j,k) - ( self % fluidVel(1) - self % Vel(1) ) * 1.0_RP / vol 
         source(3,i,j,k) = source(3,i,j,k) - ( self % fluidVel(2) - self % Vel(2) ) * 1.0_RP / vol 
         source(4,i,j,k) = source(4,i,j,k) - ( self % fluidVel(3) - self % Vel(3) ) * 1.0_RP / vol 
         source(5,i,j,k) = source(5,i,j,k) - ( self % fluidTemp   - self % temp   ) * 1.0_RP / vol 
      enddo ; enddo ; enddo        
    endif 
      
    end associate
#endif     
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine particle_source_gaussian ( self, e, source )
   implicit none
   class(Particle_t)       , intent(in)    :: self    
   class(element)          , intent(in)    :: e    
   real(kind=RP)           , intent(inout) :: source(:,0:,0:,0:)    
   ! This subroutine is not used. It could replace particle_source if I want.
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
   integer         :: i ,j, k                     
   real(kind=RP)   :: delta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3) ) 
   real(kind=RP)   :: x(3) 
   real(kind=RP)   :: val, vol  
   real(kind=RP)   :: dx 

   if ( .not. self % active ) return 

   !**************************************
   ! COMPUTE SCALING OF THE ELEMENT (dx) 
   !**************************************
   associate(spAxi   => NodalStorage(e % Nxyz(1)), &
              spAeta  => NodalStorage(e % Nxyz(2)), &
              spAzeta => NodalStorage(e % Nxyz(3)) )
              
   ! only once! 
   vol = 0.0_RP
   do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
     vol = vol + spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) * e % geom % jacobian(i,j,k)
   end do            ; end do           ; end do
   dx = vol ** (1.0_RP/3.0_RP)

   delta = 1.0_RP / vol 


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
     val = val + spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) * e % geom % jacobian(i,j,k) * delta(i,j,k)
   end do            ; end do           ; end do

   ! CON ESTO DESACTIVADO PIERDO ENERGÍA. TENGO QUE EXTENDER ESTA IDEA A LOS VECINOS. 
   ! HACER SOLO UNA VEZ POR PARTICULA
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
   
   end associate
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
    sigma = 1.5 * dx
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

end module 
#endif