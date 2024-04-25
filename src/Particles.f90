#include "Includes.h"
module ParticlesClass
#ifdef FLOW
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
    ! Mixed fluid particles to improve performance
    real(kind=RP) :: gammaDiv3cvpdivcvStPr
    real(kind=RP) :: phimDivNo_of_particlesSt
    real(kind=RP) :: phimDiv3No_of_particlesNuDivgammaminus1PrM2St
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
    integer                             :: part_per_parcel
    type(Particle_t),       allocatable :: particle(:)
    type(dimensionlessParticles_t)      :: dimensionless
    logical                             :: active 
    logical                             :: highordersource
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
subroutine ConstructParticles( self, mesh, controlVariables, solution_file )
    use FTValueDictionaryClass
    use FluidData, only : dimensionless, thermodynamics 
#if defined(SPALARTALMARAS)
    use Physics_NSSAKeywordsModule
#elif defined(NAVIERSTOKES)
    use Physics_NSKeywordsModule
#endif
    use headers
    implicit none
    class(Particles_t)        , intent(inout) :: self
    class(HexMesh)            , intent(in)    :: mesh
    class(FTValueDictionary)  , intent(in)    :: controlVariables     
    character(len=LINE_LENGTH), intent(in)    :: solution_file   
#if defined(NAVIERSTOKES)
!
!        ---------------
!        Local variables
!        ---------------
!
    integer            :: i 
    real(KIND=RP)      :: pos(3)
    real(KIND=RP)      :: vel(3)
    real(KIND=RP)      :: temp
!    real(kind=RP)      :: dy, dz, y, z, Ly, Lz
    character(LEN=LINE_LENGTH) :: partFile  
    character(LEN=1)   :: trash
    integer            :: itrash
    logical            :: velAndTempFromFile
    
    !    
    ! Read information related to particles from control file
    !--------------------------------------------------------
    
    self % no_of_particles = controlVariables % integerValueForKey(numberOfParticlesKey)
    
    self % part_per_parcel = controlVariables % integerValueForKey(particlesPerParcelKey)

    if ( controlVariables % ContainsKey(sourceTermKey) ) then
        self % highordersource = controlVariables % logicalValueForKey(sourceTermKey)
    else
        self % highordersource = .false.
    end if

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

    !    
    ! Collapse variables to improve performance
    !--------------------------------------------------------

        ! gamma / (3 * cvpdivcv * St * Pr)
    self % dimensionless % gammaDiv3cvpdivcvStPr = thermodynamics % gamma / &
        ( 3 * self % dimensionless % cvpdivcv * & 
        self % dimensionless % St * dimensionless  % Pr)

        ! phim / ( no_of_particles * St)
    self % dimensionless % phimDivNo_of_particlesSt = self % part_per_parcel * self % dimensionless % phim &
        / (self % no_of_particles * self % dimensionless % St)

        ! phim / (3 * no_of_particles) * Nu / ( gammaminus1 * Pr * Mach ** 2 * St )
    self % dimensionless % phimDiv3No_of_particlesNuDivgammaminus1PrM2St = &
        self % part_per_parcel * self % dimensionless % phim / (3 * self % no_of_particles) * self % dimensionless % Nu / &
        ( thermodynamics % gammaminus1 * dimensionless % Pr * dimensionless % Mach ** 2 * self % dimensionless % St )


    partFile = controlVariables % StringValueForKey(key = PART_FILE_KEY, requestedLength = LINE_LENGTH)
    velAndTempFromFile = controlVariables % logicalValueForKey(PART_LOG_FILE_KEY)

    self % pMesh % min = getRealArrayFromString( controlVariables % StringValueForKey(key = MIN_BOX_KEY,&
    requestedLength = 132))
    self % pMesh % max = getRealArrayFromString( controlVariables % StringValueForKey(key = MAX_BOX_KEY,&
    requestedLength = 132))
    self % pMesh % bc  = getIntArrayFromString( controlVariables % StringValueForKey(key = BC_BOX_KEY,&
    requestedLength = 132))

    self % injection % active = controlVariables % logicalValueForKey(PART_LOG_INJ_KEY)

    !    
    ! Set up injection if required
    !--------------------------------------------------------

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

    !    
    ! Initialize particles
    !--------------------------------------------------------

    if (self % injection % active) then 
        do i = 1, self % no_of_particles 
            self % particle (i) % active = .false. 
        enddo 
        self % injection % injected = 0
    else 
        self % injection % injected = self % no_of_particles - 1
        open(UNIT=10, FILE=partFile)
        read(10,*)
        do i = 1, self % no_of_particles 
            ! Read position of the particles from file
            read(10,*) itrash, &
                       pos(1), pos(2), pos(3), &
                       vel(1), vel(2), vel(3), temp

            if (itrash .ne. i ) then 
                write(*,*) "Particles missing in the initialization from file, reducing number of particles."
                write(*,*) "This functionality is in beta mode. Check code for more details."
                itrash = itrash + 1
                self % injection % injected = self % injection % injected - 1
                ! This has not been tested. 
                ! Includes this if + the update of line 238 of self % no_of_particles
            endif 

            call self % particle(i) % set_pos ( pos )

            ! Position the particle in the computational mesh (get element)
            call self % particle(i) % setGlobalPos( mesh )

            if (velAndTempFromFile) then 
                ! Initialise particle velocity and temperature from file
                call self % particle(i) % set_vel  ( vel  )
                call self % particle(i) % set_temp ( temp )
            else 

                ! Get Fluid velocity and temperature for the random positions of each particle
                call self % particle(i) % getFluidVelandTemp( mesh )
                
                ! Initialise particle velocity and temperature with the fluid values
                call self % particle(i) % set_vel  ( self % particle(i) % fluidVel  )
                call self % particle(i) % set_temp ( self % particle(i) % fluidTemp )
            endif 
        enddo 
        self % no_of_particles      = self % injection % injected + 1
        close(10)
    endif 

    call ExportToVTKParticles( self, 0, solution_file )

    !    
    ! Show particles at screen
    !--------------------------------------------------------

    write(STD_OUT,'(30X,A,A28,L)')   "->" , "Injection active: " , self % injection % active
    write(STD_OUT,'(30X,A,A28,L)')   "->" , "High order source: " , self % highordersource
    write(STD_OUT,'(30X,A,A28,A132)')   "->" , "Initialization file: " , partFile
    write(STD_OUT,'(30X,A,A30,E10.3)') "->", "Stokes number: ", self % dimensionless % St
    write(STD_OUT,'(30X,A,A30,E10.3)') "->", "phim: ", self % dimensionless % phim
    write(STD_OUT,'(30X,A,A30,E10.3)') "->", "Gamma (cvpdivcv): ", self % dimensionless % cvpdivcv
    write(STD_OUT,'(30X,A,A30,E10.3)') "->", "Nusselt: ", self % dimensionless % Nu
    write(STD_OUT,'(30X,A,A30,E10.3)') "->", "I0: ", self % dimensionless % I0

    write(STD_OUT,'(30X,A,A20,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "minimum box: ","[", &
    self % pMesh % min(1), ", ", &
    self % pMesh % min(2), ", ", &
    self % pMesh % min(3), "]"

    write(STD_OUT,'(30X,A,A20,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "maximum box: ","[", &
    self % pMesh % max(1), ", ", &
    self % pMesh % max(2), ", ", &
    self % pMesh % max(3), "]"

    write(STD_OUT,'(30X,A,A20,A,i4.1,A,i4.1,A,i4.1,A, A61)') "->" , "bc box: ","[", &
    self % pMesh % bc(1), ", ", &
    self % pMesh % bc(2), ", ", &
    self % pMesh % bc(3), "]",  "     // [i,j,k] 0 is inflow/outflow, 1 is wall, 2 is periodic"

    if ( self % injection % active ) then 
        write(STD_OUT,'(30X,A,A40,A,i4.1,A,i4.1,A,i4.1,A)') "->" , "Injection axis: ","[", &
        self % injection % axis(1), ", ", &
        self % injection % axis(2), ", ", &
        self % injection % axis(3), "]"
        write(STD_OUT,'(30X,A,A40,I7)') "->", "Injection particles per step: ", self % injection % number 
        write(STD_OUT,'(30X,A,A40,I7)') "->", "Injection iter period: ", self % injection % period     
        write(STD_OUT,'(30X,A,A40,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "Dimensionless Injection velocity: ","[", &
        self % injection % v (1), ", ", &
        self % injection % v (2), ", ", &
        self % injection % v (3), "]"        
        write(STD_OUT,'(30X,A,A40,E10.3)') "->", "Dimensionless Injection temp: ", self % injection % T   
    endif    

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
!
!        ---------------
!        Local variables
!        ---------------
!
    integer     :: i


    !GTD: Change to for all particles
    !GTD: Adapt for MPI compatibility
    !GTD: Performance can be improved (a lot)
    !GTD: Fix makefile dependencies (now I have to sometimes do make clean)

    !$omp parallel do schedule(runtime)       
    do i = 1, self % injection % injected + 1
        if (self % particle(i) % active) then 
        !    
        ! Get particle global position and set up interpolation
        !------------------------------------------------------
            !call self % particle(i) % show
            call self % particle(i) % setGlobalPos( mesh )
            !call self % particle(i) % show
        !    
        ! Get fluid velocity and temperature at that position
        !------------------------------------------------------ 
            ! PREVENTS TRYING TO GET VEL AND TEMP OUTSIDE DOMAIN  
            if ( self % particle(i) % active ) then 
                call self % particle(i) % getFluidVelandTemp( mesh )
            endif    
        !        
        ! Integrate in time to get new particle velocity and temperature
        !------------------------------------------------------   
            call self % particle(i) % integrate( mesh, dt, &
                                                    self % dimensionless % St, &
                                                    self % dimensionless % Nu, &
                                                    self % dimensionless % phim, &
                                                    self % dimensionless % I0, &
                                                    self % dimensionless % gammaDiv3cvpdivcvStPr, &
                                                    self % pMesh % min, & 
                                                    self % pMesh % max, &
                                                    self % pMesh % bc ) 
        !    
        ! Print particle position (only for debugging)
        !------------------------------------------------------           
!                 call self % particle(i) % show()  
        endif 

    enddo 
    !$omp end parallel do
        !    
        ! Print particle position (only for debugging)
        !------------------------------------------------------           
  
#endif
end subroutine IntegrateParticles
!
!///////////////////////////////////////////////////////////////////////////////////////
!
subroutine AddSourceParticles( self, iP, e, time, thermodynamics_, dimensionless_, refValues_ )
    USE ElementClass
#if defined(NAVIERSTOKES)
    use VariableConversion, only : temperature, sutherlandsLaw
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


        associate ( d => self % dimensionless )

        call self % computeSourceTerm( e, iP, Source )                          

        !   Compute source term in coordinate i, j, k

        do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            Source(1, i, j, k ) = 0.0_RP

            Source(2:4, i, j, k ) = Source(2:4, i, j, k ) * d % phimDivNo_of_particlesSt * &
                    SutherlandsLaw( Temperature( e % storage % Q(:,i,j,k) ) ) 

            Source(5, i, j, k )   = Source(5, i, j, k )   * d % phimDiv3No_of_particlesNuDivgammaminus1PrM2St
        enddo ; enddo ; enddo 
!$omp critical
        do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            e % storage % S_NSP(1:5,i,j,k) = e % storage % S_NSP(1:5,i,j,k) + Source(1:5,i,j,k)
        enddo ; enddo ; enddo 
!$omp end critical

        end associate 
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

    Source = 0.0_RP
    call self % particle(iP) % Source( e, Source, self % highordersource )

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
    real(KIND=RP)     :: eps

    ! This is a hardcoded value that makes sure that the particle is inside the domain.
    ! Injection position is Min + (Max - Min) * eps or Max - (Max - Min) * eps depending if it is in positive or negative directions
    eps = 1.d-2 

    if ( self % injection % injected + self % injection % number + 1  > self % no_of_particles ) then 
        return 
    endif 

    ! Injection velocity
    v = self % injection % v 
    ! Injection temperature
    T = self % injection % T
    
!!!!$omp do schedule(runtime) private(k, pos)
!$omp single
    do i = self % injection % injected, self % injection % injected + self % injection % number

        do k = 1,3
            call random_number( pos(k) )
            pos(k) = ( pos(k) - 0.5_RP ) * ( self % pMesh % max(k) - self % pMesh % min(k) ) + &
                ( self % pMesh % max(k) + self % pMesh % min(k)  ) / 2 
        enddo
    
        do k = 1, 3
            if ( self % injection % axis (k) /= 0 ) then 
                if ( self % injection % axis (k) == 1 ) then 
                    pos(k) = self % pMesh % min(k) + ( self % pMesh % max(k) - self % pMesh % min(k) ) * eps
                elseif ( self % injection % axis (i) == - 1 ) then 
                    pos(k) = self % pMesh % max(k) - ( self % pMesh % max(k) - self % pMesh % min(k) ) * eps
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
!$omp end single
!!!!!$omp end do 
    self % injection % injected = self % injection % injected + self % injection % number

#endif
end subroutine InjectParticles
!
!///////////////////////////////////////////////////////////////////////////////////////
!
#endif
end module ParticlesClass