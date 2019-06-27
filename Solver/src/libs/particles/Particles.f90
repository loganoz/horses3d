!
!//////////////////////////////////////////////////////
!
!   @File:    Particles.f90
!   @Author:  Gonzalo (g.rubio@upm.es)
!   @Created: Tue Apr 10 17:31:22 2018
!   @Last revision date: Thu Sep 27 16:42:15 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 5ab4fc5764dead65069a92d809d881f964ea4900
!
!//////////////////////////////////////////////////////
!
module ParticlesClass
#if defined(NAVIERSTOKES) || defined(INCNS)
   use SMConstants
   use ParticleClass
   use FluidData
   use PhysicsStorage
   use HexMeshClass
   use FileReadingUtilities      , only: getRealArrayFromString, getIntArrayFromString
   implicit none
!
#include "Includes.h"

private
public  Particles_t
!
!  *******************************
!  Main particles class definition
!  *******************************
!  
type DimensionlessParticles_t
    real(kind=RP) :: St
    real(kind=RP) :: Nu
    real(kind=RP) :: phim
    real(kind=RP) :: cvpdivcv
    real(kind=RP) :: I0
end type DimensionlessParticles_t

! The implementation of particles is limited to boxes right now. It should be easy
! to couple this with the boundary conditions of the solver. I do not have the time
! to do it properly right now. 
type pMesh_t 
    real(kind=RP)  :: min(3) ! minimum dimension of box
    real(kind=RP)  :: max(3) ! maximum dimension of box
    integer  :: bc(3)  ! boundary condition (inflow/outflow [0], wall [1], periodic [2])
end type pMesh_t

type injection_t 
    logical  :: active
    integer  :: axis(3)  ! [ , , ] direction of injection
    integer  :: number   ! particles injected per step
    integer  :: period   ! period of iterations between injections 
    integer  :: injected ! number of particles injected
    real(KIND=RP) :: v(3) ! Injection velocity
    real(KIND=RP) :: T    ! Injection temperature
end type injection_t

type Particles_t
    integer                             :: no_of_particles
    type(Particle_t),       allocatable :: particle(:)
    type(dimensionlessParticles_t)      :: dimensionless
    logical                             :: active 
    type(pMesh_t)                       :: pMesh
    type(injection_t)                   :: injection 
    contains
        procedure   :: Construct      => ConstructParticles    
        procedure   :: Integrate      => IntegrateParticles
        procedure   :: ExportToVTK    => ExportToVTKParticles
        procedure   :: AddSource      => AddSourceParticles
        procedure   :: ComputeSourceTerm => ComputeSourceTermParticles
        procedure   :: Inject   => InjectParticles
end type Particles_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine ConstructParticles( self, mesh, controlVariables )
    use FTValueDictionaryClass
#if defined(NAVIERSTOKES)
    use Physics_NSKeywordsModule
#endif
    use headers
    implicit none
    class(Particles_t)      , intent(inout) :: self
    class(HexMesh)          , intent(in)    :: mesh
    class(FTValueDictionary), intent(in)    :: controlVariables        
#if defined(NAVIERSTOKES)
    !TDG: understand difference between class and type    
    ! http://www.pgroup.com/lit/articles/insider/v3n1a3.htm
!
!        ---------------
!        Local variables
!        ---------------
!
    integer            :: i 
    real(KIND=RP)      :: pos(3)
    real(KIND=RP)      :: vel(3)
    real(KIND=RP)      :: temp
    real(kind=RP)      :: dy, dz, y, z, Ly, Lz
    character(LEN=132) :: partFile 
    
    self % no_of_particles = controlVariables % integerValueForKey(numberOfParticlesKey)
    
    if (self % no_of_particles == huge(1)) then 
        self % no_of_particles = 0
    else
        write(STD_OUT,'(/)')
        call Section_Header('Loading particles')
        write(STD_OUT,'(/)')
        write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of particles: " , self % no_of_particles
    endif 

    allocate( self % particle( self % no_of_particles) ) 

    do i = 1, self % no_of_particles 
        call self%particle(i)%init(mesh)
    enddo 

    self % dimensionless % St       = controlVariables % doublePrecisionValueForKey(STOKES_NUMBER_PART_KEY)    
    self % dimensionless % phim     = controlVariables % doublePrecisionValueForKey(PHI_M_PART_KEY) 
    self % dimensionless % cvpdivcv = controlVariables % doublePrecisionValueForKey(GAMMA_PART_KEY) 
    self % dimensionless % Nu       = 2.0_RP !Stokeian flow
    self % dimensionless % I0       = controlVariables % doublePrecisionValueForKey(I0_PART_KEY) 

    partFile = controlVariables % StringValueForKey(key = PART_FILE_KEY, requestedLength = 132)

    self % pMesh % min = getRealArrayFromString( controlVariables % StringValueForKey(key = MIN_BOX_KEY,&
    requestedLength = 132))
    self % pMesh % max = getRealArrayFromString( controlVariables % StringValueForKey(key = MAX_BOX_KEY,&
    requestedLength = 132))
    self % pMesh % bc  = getIntArrayFromString( controlVariables % StringValueForKey(key = BC_BOX_KEY,&
    requestedLength = 132))

    self % injection % active = controlVariables % logicalValueForKey("injection")

    if (self % injection % active) then 
        self % injection % axis   = getIntArrayFromString( controlVariables % StringValueForKey(key = PART_INJ_KEY,&
    requestedLength = 132))
        self % injection % number = controlVariables % integerValueForKey(key = PART_NUMB_PER_STEP_KEY)
        self % injection % period = controlVariables % integerValueForKey(key = PART_PERIOD_KEY)
        self % injection % v      = getRealArrayFromString( controlVariables % StringValueForKey(key = INJ_VEL_KEY,&
        requestedLength = 132))
        self % injection % T      = controlVariables % doublePrecisionValueForKey(INJ_TEMP_KEY)        
    else 
        self % injection % axis   = [0,0,0]
        self % injection % number = 0
        self % injection % period = 0
        self % injection % v      = 0.0_RP
        self % injection % T      = 0.0_RP
    endif 

    write(STD_OUT,'(30X,A,A28,A30)') "->" , "Initialization file: " , partFile

! print*, " self % dimensionless % St ",  self % dimensionless % St !, STOKES_NUMBER_PART_KEY
! print*, " self % dimensionless % phim ",  self % dimensionless % phim !, STOKES_NUMBER_PART_KEY
! print*, " self % dimensionless % cvpdivcv ",  self % dimensionless % cvpdivcv !, STOKES_NUMBER_PART_KEY
! print*, " self % dimensionless % I0 ", self % dimensionless % I0 !, I0_PART_KEY
! print*, " self % dimensionless % g ",  self % dimensionless % g !, STOKES_NUMBER_PART_KEY
! print*, "self % no_of_particles", self % no_of_particles
    !TDG: if particles go out of the domain, they should stop being tracked down and 
    !       new particles should be created. The implementation performed right now
    !       makes that the particles that go out of the domain are not tracked down.


    ! vel  = (/0.d0, 0.d0, 0.d0/)
    ! temp = 1.d0  

    ! dy = 3e-2_RP
    ! dz = 3e-2_RP   
    ! y = 1e-2_RP
    ! z = 1e-2_RP  
    ! Ly = pi
    ! Lz = 2.0_RP
    ! pos(1) = 1.e-2_RP !1.e-2_RP
    ! do i = 1, self % no_of_particles
    !     y = y + dy
    !     if ( y < Ly ) then 
    !         pos(2) = y
    !         pos(3) = z
    !     else 
    !         y = dy
    !         z = z + dz
    !     endif
    !     if ( z > Lz) then 
    !         write(*,*) "There is no space for so many particles"
    !         write(*,*) i
    !         stop
    !     endif 
    
    !     call self % particle(i) % set_pos ( pos )
    !     call self % particle(i) % set_vel ( vel )
    !     call self % particle(i) % set_temp( temp )        
    ! enddo 
    
    ! vel  = (/0.d0, 0.d0, 0.0d0/)
    ! temp = 1.d0

    if (self % injection % active) then 
        !print*, "Particle initialization file missing. No particles in the domain at t=0."
        do i = 1, self % no_of_particles 
            self % particle (i) % active = .false. 
        enddo 
        self % injection % injected = 0
    else 
        self % injection % injected = self % no_of_particles - 1
        open(UNIT=10, FILE=partFile)
        read(10,*)
        do i = 1, self % no_of_particles 
            ! Read position of the particles from RandomParticles.txt
            read(10,*) pos(1), pos(2), pos(3)
            call self % particle(i) % set_pos ( pos )

            ! Position the particle in the computational mesh (get element)
            call self % particle(i) % setGlobalPos( mesh )

            ! Get Fluid velocity and temperature for the random positions of each particle
            call self % particle(i) % getFluidVelandTemp( mesh )
            
            ! Initialise particle velocity and temperature with the fluid values
            call self % particle(i) % set_vel  ( self % particle(i) % fluidVel  )
            call self % particle(i) % set_temp ( self % particle(i) % fluidTemp )
        enddo 
        close(10)
    endif 
    ! pos (1) = 1.d0
    ! pos (2) = 1.d0
    ! pos (3) = 1.d0 
    ! vel  = (/0.d0, 0.d0, 0.d0/)
    ! temp = 1.d0  
    ! call self % particle(1) % set_pos ( pos )
    ! call self % particle(1) % set_vel ( vel )
    ! call self % particle(1) % set_temp( temp )        
#endif
end subroutine ConstructParticles
!
!///////////////////////////////////////////////////////////////////////////////////
!
subroutine IntegrateParticles( self, mesh, dt )
    implicit none
    class(HexMesh)          , intent(in)     :: mesh     
    class(Particles_t)      , intent(inout)  :: self
    real(KIND=RP)           , intent(in)     :: dt
#if defined(NAVIERSTOKES)
    !Debugging variables
    logical :: debug  = .false.
    logical :: debug2 = .false. 
    real    :: t1, t2, t3, t4, t5, t6
    real    :: gt1, gt2, gt3
!
!        ---------------
!        Local variables
!        ---------------
!
    integer     :: i
    !integer     :: no_of_particles

    !no_of_particles = size( self % particle )


    !GTD: Change to for all particles
    !GTD: Adapt for MPI compatibility
    !GTD: Performance can be improved (a lot)
    !GTD: OpenMP parallelization
    !GTD: Stop following particle if it gets out of the domain 
    !GTD: Fix makefile dependencies (now I have to sometimes do make clean)
gt1 = 0.d0
gt2 = 0.d0
gt3 = 0.d0
    !$omp parallel do schedule(runtime)       
    do i = 1, self % injection % injected + 1
        if (self % particle(i) % active) then 
        !    
        ! Get particle global position and set up interpolation
        !------------------------------------------------------
        if (debug) call cpu_time(t1)
            call self % particle(i) % setGlobalPos( mesh )
            ! call self % particle(i) % show() 
            ! read(*,*)
        if (debug) call cpu_time(t2)
        if (debug2) gt1 = gt1 + (t2 - t1)
        !    
        ! Get fluid velocity and temperature at that position
        !------------------------------------------------------   
        if (debug) call cpu_time(t3)
        if ( self % particle(i) % active ) then 
            call self % particle(i) % getFluidVelandTemp( mesh )
        endif 
        if (debug) call cpu_time(t4)
        if (debug2) gt2 = gt2 + (t4 - t3)        
        !        
        ! Integrate in time to get new particle velocity and temperature
        !------------------------------------------------------   
        if (debug) call cpu_time(t5)
            call self % particle(i) % integrate( dt, &
                                                    self % dimensionless % St, &
                                                    self % dimensionless % Nu, &
                                                    self % dimensionless % phim, &
                                                    self % dimensionless % cvpdivcv, &
                                                    self % dimensionless % I0, &
                                                    self % pMesh % min, & 
                                                    self % pMesh % max, &
                                                    self % pMesh % bc ) 
        if (debug) call cpu_time(t6)
        if (debug2) gt3 = gt3 + (t6 - t5)
        !    
        ! Print particle position (only for debugging)
        !------------------------------------------------------           
!                 call self % particle(i) % show()  
        endif 
!        call self % particle(i) % show()  
        ! if (debug) then 
        !     print*, "Particle number", i,"/",self % no_of_particles
        !     print*, "Set Global Position", t2 - t1, "seconds"
        !     print*, "Get fluid v and T  ", t4 - t3, "seconds"
        !     print*, "Integrate particle ", t6 - t5, "seconds"       
        !     print*, "Press space bar to continue"     
        !     read(*,*)       
        ! endif 

    enddo 
    !$omp end parallel do
        !    
        ! Print particle position (only for debugging)
        !------------------------------------------------------           
                ! call self % particle(1000) % show()  
if (debug2) then 
    print*, "overall time", gt1+gt2+gt3, "seconds"
    print*, "set pos", gt1 / (gt1+gt2+gt3) * 100, "%"
    print*, "get vel", gt2 / (gt1+gt2+gt3) * 100, "%"
    print*, "integra", gt3 / (gt1+gt2+gt3) * 100, "%"
endif 
#endif
end subroutine IntegrateParticles
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine AddSourceParticles( self, iP, e, time, thermodynamics_, dimensionless_, refValues_ )
    USE ElementClass
#if defined(NAVIERSTOKES)
    use Physics,            only : sutherlandsLaw
    use VariableConversion, only : temperature 
#endif
    IMPLICIT NONE
    class(Particles_t)      , intent(in)    :: self
    integer                 , intent(in)    :: iP 
    CLASS(element)          , intent(inout) :: e
    REAL(KIND=RP)           , intent(in)    :: time
    type(Thermodynamics_t)  , intent(in)    :: thermodynamics_
    type(Dimensionless_t)   , intent(in)    :: dimensionless_
    type(RefValues_t)       , intent(in)    :: refValues_
#if defined(NAVIERSTOKES)    
!
!        ---------------
!        Local variables
!        ---------------
!
    integer                     :: i, j, k, eID
    real(KIND=RP)               :: Source( NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3) )

    real :: t1,t2

!    call cpu_time(t1)
!    if ( self % no_of_particles > 0 ) then
        associate ( d => self % dimensionless )
        !allocate( Source( NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3) ) )  
        call self % computeSourceTerm( e, iP, Source )                          
        !do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

            !   Compute source term in coordinate i, j, k
            !   Contributions of all the particles are taken into account
            !   TDG: these numbers that are going to be used a lot of times. Store them somewhere

            Source(1, :, :, : ) = 0.0_RP

            Source(2:4, :, :, :) = Source(2:4, :, :, :) * d % phim / ( self % no_of_particles * d % St) * &
                SutherlandsLaw( Temperature( e % storage % Q(:,i,j,k) ) ) 

            Source(5, :, :, :)   = Source(5, :, :, :)   * d % phim / (3 * self % no_of_particles) * d % Nu &
                                        / ( thermodynamics_ % gammaminus1 * dimensionless_ % Pr * &
                                        dimensionless_ % Mach ** 2 * d % St )
        !end do                  ; end do                ; end do

            ! Add to NS source term
            !e % storage % S_NS(:,i,j,k) = e % storage % S_NS(:,i,j,k) + Source(:,i,j,k)
!$omp critical
        e % storage % S_NSP = e % storage % S_NSP + Source
!$omp end critical
            !   Add source term to Qdot                           
            !e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source(:,i,j,k)

        !deallocate(Source)
        end associate 
!    else 
        !Nothing
!    endif 
!    call cpu_time(t2)

!    print*, "Compute and add source", t2-t1, "seconds"
#endif
end subroutine AddSourceParticles    
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine ComputeSourceTermParticles(self, e, iP, Source)
    USE ElementClass    
    IMPLICIT NONE
    class(Particles_t)      , intent(in)    :: self
    class(element)          , intent(in)    :: e     
    integer                 , intent(in)    :: iP 
    real(KIND=RP)           , intent(out)   :: Source(:,0:,0:,0:)
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
! 
    integer       :: neighbours(6),j 

    ! Find the neighbours of the element in which the particle is
    ! neighbours = -1
    ! do j = 1, 6
    !    if (self % particle(iP) % mesh%elements(e % eID)%NumberOfConnections(j) > 0) then
    !          neighbours(j) = self % particle(iP) % mesh%elements(e % eID)%Connection(j)%ElementIDs(1) 
    !    else 
    !          neighbours(j) = -1
    !    end if
    ! end do     

    Source = 0.0_RP
!    do iP = 1, self % injection % injected + 1
!        if (self % particle(iP) % active ) then
            !if (self % particle(iP) % eID == e % eID) then 
                call self % particle(iP) % Source( e, Source )
            !else
                ! do j = 1,6
                !     if ( neighbours(j) == e % eID ) then 
                !         call self % particle(iP) % Source( e, Source )    
                !     endif 
                ! enddo 
                ! if it is not in the element it does not affect
                ! the element 
           

            !endif 
            ! print*, "e % eID", e % eID
            ! print*, "neighbours", neighbours
            ! read(*,*)
            !Source = Source + self % particle(iP) % Source()
 !       endif 
!    enddo 
#endif
end subroutine 
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine ExportToVTKParticles( self, iter, solution_file )
    implicit none
    class(Particles_t)      , intent(in)  :: self
    character(len=LINE_LENGTH)            :: solution_file
    integer                               :: iter 
#if defined(NAVIERSTOKES)

!
!        ---------------
!        Local variables
!        ---------------
!
    character(len=LINE_LENGTH)   :: filename 
    integer     :: i
    integer     :: no_of_particles

    
    no_of_particles = size( self % particle )

    write(fileName,'(A,A,I10.10,A)') trim(solution_file),'.parts.',iter,'.csv'

    open(FILE = fileName, UNIT=10)

    write(10,*) "i", ",", "x coord", ",", "y coord", ",", "z coord", ",", "u", ",", "v", ",", "w", ",", "T"
    do i = 1, no_of_particles
        if ( self % particle (i) % active ) then 
            write( 10, '(i10,7(A,E12.6))' ) i, ",", &
            self % particle (i) % pos(1), ",", self % particle (i) % pos(2), ",", self % particle (i) % pos(3), ",",&
            self % particle (i) % vel(1), ",", self % particle (i) % vel(2), ",", self % particle (i) % vel(3),",", &
            self % particle (i) % temp
        endif 
    enddo 

    close(10)
#endif
end subroutine ExportToVTKParticles
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine InjectParticles( self, mesh  )
    implicit none
    class(Particles_t)      , intent(inout)  :: self
    class(HexMesh)          , intent(in)     :: mesh
#if defined(NAVIERSTOKES)

!
!        ---------------
!        Local variables
!        ---------------
!
    integer           :: i, k 
    real(KIND=RP)     :: pos(3)
    real(KIND=RP)     :: v(3), T
!    real(KIND=RP)     :: eps

    ! ! This is a harcoded value that makes sure that the particle is inside the domain.
    ! eps = 1.d-1 

    if ( self % injection % injected + self % injection % number + 1  > self % no_of_particles ) then 
        !print*, "Maximum number of particles reached."
        !print*, "No more particles will be added to this simulation."
        return 
    endif 

    ! Injection velocity
    v = self % injection % v 
    ! Injection temperature
    T = self % injection % T
    
!$omp do schedule(runtime) private(k, pos)
    do i = self % injection % injected, self % injection % injected + self % injection % number

        do k = 1,3
            call random_number( pos(k) )
            pos(k) = ( pos(k) - 0.5_RP ) * ( self % pMesh % max(k) - self % pMesh % min(k) ) + &
                ( self % pMesh % max(k) + self % pMesh % min(k)  ) / 2 
        enddo
    
        do k = 1, 3
            if ( self % injection % axis (k) /= 0 ) then 
                if ( self % injection % axis (k) == 1 ) then 
                    pos(k) = self % pMesh % min(k) + ( self % pMesh % max(k) - self % pMesh % min(k) ) / 100
                elseif ( self % injection % axis (i) == - 1 ) then 
                    pos(k) = self % pMesh % max(k) - ( self % pMesh % max(k) - self % pMesh % min(k) ) / 100
                endif 
            endif 
        enddo      

        call self % particle(i+1) % set_pos ( pos )

        ! Position the particle in the computational mesh (get element)
        call self % particle(i+1) % setGlobalPos( mesh )

        ! Get Fluid velocity and temperature for the random positions of each particle
        ! It will be required for the computation of the source term.
        call self % particle(i+1) % getFluidVelandTemp( mesh )
        ! Initialise particle velocity and temperature with the fluid values
        !call self % particle(i+1) % set_vel  ( self % particle(i+1) % fluidVel  )
        !call self % particle(i+1) % set_temp ( self % particle(i+1) % fluidTemp )        

        call self % particle(i+1) % set_vel  ( v )
        call self % particle(i+1) % set_temp ( T )      

    enddo 
!$omp end do 
    self % injection % injected = self % injection % injected + self % injection % number

#endif
end subroutine InjectParticles
!
!///////////////////////////////////////////////////////////////////////////////////////
!
#endif
end module ParticlesClass
