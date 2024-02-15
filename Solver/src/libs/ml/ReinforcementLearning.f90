!
!//////////////////////////////////////////////////////
!
!      Module containing classes and subroutines for Reinforcement Learning
!      Currently implemented:
!        * p-Adaptation agent based on Value Iteration.
!
!////////////////////////////////////////////////////////////////////////
!

module ReinforcementLearning

    implicit none

    private 

    public pAdaptationAgent_t

    type, abstract :: MatrixND_t
        contains
            procedure(MatrixND_ConstructInterface), deferred :: MatrixND_Construct
            generic :: construct => MatrixND_Construct
            procedure(MatrixND_DestructInterface), deferred :: MatrixND_Destruct
            generic :: destruct  => MatrixND_Destruct
            procedure(MatrixND_LoadInterface), deferred :: MatrixND_Load
            generic :: load  => MatrixND_Load
            procedure(MatrixND_GetDataInterface), deferred :: MatrixND_GetData
            generic :: getData  => MatrixND_GetData
    end type MatrixND_t

    type, extends(MatrixND_t) :: Matrix2D_t
        integer(kind=1), allocatable :: data(:,:)
        contains
            procedure :: MatrixND_Construct => Matrix2D_Construct
            procedure :: MatrixND_Destruct  => Matrix2D_Destruct
            procedure :: MatrixND_Load      => Matrix2D_Load
            procedure :: MatrixND_GetData   => Matrix2D_GetData
    end type Matrix2D_t

    type, extends(MatrixND_t) :: Matrix3D_t
        integer(kind=1), allocatable :: data(:,:,:)
        contains
            procedure :: MatrixND_Construct => Matrix3D_Construct
            procedure :: MatrixND_Destruct  => Matrix3D_Destruct
            procedure :: MatrixND_Load      => Matrix3D_Load
            procedure :: MatrixND_GetData   => Matrix3D_GetData
    end type Matrix3D_t

    type, extends(MatrixND_t) :: Matrix4D_t
        integer(kind=1), allocatable :: data(:,:,:,:)
        contains
            procedure :: MatrixND_Construct => Matrix4D_Construct
            procedure :: MatrixND_Destruct  => Matrix4D_Destruct
            procedure :: MatrixND_Load      => Matrix4D_Load
            procedure :: MatrixND_GetData   => Matrix4D_GetData
    end type Matrix4D_t

    type, extends(MatrixND_t) :: Matrix5D_t
        integer(kind=1), allocatable :: data(:,:,:,:,:)
        contains
            procedure :: MatrixND_Construct => Matrix5D_Construct
            procedure :: MatrixND_Destruct  => Matrix5D_Destruct
            procedure :: MatrixND_Load      => Matrix5D_Load
            procedure :: MatrixND_GetData   => Matrix5D_GetData
    end type Matrix5D_t

    type, extends(MatrixND_t) :: Matrix6D_t
        integer(kind=1), allocatable :: data(:,:,:,:,:,:)
        contains
            procedure :: MatrixND_Construct => Matrix6D_Construct
            procedure :: MatrixND_Destruct  => Matrix6D_Destruct
            procedure :: MatrixND_Load      => Matrix6D_Load
            procedure :: MatrixND_GetData   => Matrix6D_GetData
    end type Matrix6D_t

    type, extends(MatrixND_t) :: Matrix7D_t
        integer(kind=1), allocatable :: data(:,:,:,:,:,:,:)
        contains
            procedure :: MatrixND_Construct => Matrix7D_Construct
            procedure :: MatrixND_Destruct  => Matrix7D_Destruct
            procedure :: MatrixND_Load      => Matrix7D_Load
            procedure :: MatrixND_GetData   => Matrix7D_GetData
    end type Matrix7D_t

    type :: adaptiveArray_t
        class(MatrixND_t), allocatable :: matrix
    end type adaptiveArray_t

    type :: pAdaptationAgent_t
        integer                            :: smax      ! Maximum value of the state
        integer                            :: pmin      ! Minimum polynomial order
        integer                            :: pmax      ! Maximum polynomial order
        type(adaptiveArray_t), allocatable :: policy(:) ! Policy
        contains
            procedure :: construct => pAdaptationAgent_Construct
            procedure :: destruct => pAdaptationAgent_Destruct
    end type pAdaptationAgent_t

    !  ------------------------------------------
!  Interface for constructing the MatrixND classs
!  ----------------------------------------------
   interface
   subroutine MatrixND_ConstructInterface(self, Npoints)
      import MatrixND_t
      !-------------------------------------------------
      implicit none
      !-------------------------------------------------
      class(MatrixND_t) , intent(inout) :: self   
      integer           , intent(in)    :: Npoints      
   end subroutine MatrixND_ConstructInterface
   end interface

!  --------------------------------------------
!  Interface for destructing the MatrixND class
!  --------------------------------------------
   interface
   subroutine MatrixND_DestructInterface(self)
      import MatrixND_t
      !--------------------------------------
      implicit none
      !--------------------------------------
      class(MatrixND_t) , intent(inout) :: self
   end subroutine MatrixND_DestructInterface
   end interface

!  -----------------------------------------------------------
!  Interface for loading the MatrixND class from a binary file
!  -----------------------------------------------------------
   interface
   subroutine MatrixND_LoadInterface(self, agentFile)
      import MatrixND_t
      !--------------------------------------
      implicit none
      !--------------------------------------
      class(MatrixND_t) , intent(inout) :: self
      character(len=*)  , intent(in)    :: agentFile
   end subroutine MatrixND_LoadInterface
   end interface

!  -------------------------------------------------------
!  Interface for retrieving the data from a MatrixND class
!  -------------------------------------------------------
   interface
   function MatrixND_GetDataInterface(self, indices) result(val)
      import MatrixND_t
      !--------------------------------------
      implicit none
      !--------------------------------------
      class(MatrixND_t)     , intent(in)    :: self
      integer               , intent(in)    :: indices(:)
      integer(kind=1)                       :: val
   end function MatrixND_GetDataInterface
   end interface

    !========
    contains
    !========
    !
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  -------------------------------------------
    !  Routine for constructing the MatrixND class
    !  -------------------------------------------
    subroutine Matrix2D_Construct(self, Npoints)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix2D_t), intent(inout) :: self
        integer          , intent(in)    :: Npoints
        !----------------------------------------------------
        allocate(self % data(Npoints,Npoints))
    end subroutine Matrix2D_Construct

    subroutine Matrix3D_Construct(self, Npoints)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix3D_t), intent(inout) :: self
        integer          , intent(in)    :: Npoints
        !----------------------------------------------------
        allocate(self % data(Npoints,Npoints,Npoints))
    end subroutine Matrix3D_Construct

    subroutine Matrix4D_Construct(self, Npoints)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix4D_t), intent(inout) :: self
        integer          , intent(in)    :: Npoints
        !----------------------------------------------------
        allocate(self % data(Npoints,Npoints,Npoints,Npoints))
    end subroutine Matrix4D_Construct

    subroutine Matrix5D_Construct(self, Npoints)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix5D_t), intent(inout) :: self
        integer          , intent(in)    :: Npoints
        !----------------------------------------------------
        allocate(self % data(Npoints,Npoints,Npoints,Npoints,Npoints))
    end subroutine Matrix5D_Construct

    subroutine Matrix6D_Construct(self, Npoints)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix6D_t), intent(inout) :: self
        integer          , intent(in)    :: Npoints
        !----------------------------------------------------
        allocate(self % data(Npoints,Npoints,Npoints,Npoints,Npoints,Npoints))
    end subroutine Matrix6D_Construct

    subroutine Matrix7D_Construct(self, Npoints)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix7D_t), intent(inout) :: self
        integer          , intent(in)    :: Npoints
        !----------------------------------------------------
        allocate(self % data(Npoints,Npoints,Npoints,Npoints,Npoints,Npoints,Npoints))
    end subroutine Matrix7D_Construct
    !
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  -------------------------------------------
    !  Routine for destructing the MatrixNd class
    !  -------------------------------------------
    subroutine Matrix2D_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix2D_t), intent(inout) :: self
        !----------------------------------------------------
        deallocate(self % data)
    end subroutine Matrix2D_Destruct

    subroutine Matrix3D_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix3D_t), intent(inout) :: self
        !----------------------------------------------------
        deallocate(self % data)
    end subroutine Matrix3D_Destruct

    subroutine Matrix4D_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix4D_t), intent(inout) :: self
        !----------------------------------------------------
        deallocate(self % data)
    end subroutine Matrix4D_Destruct

    subroutine Matrix5D_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix5D_t), intent(inout) :: self
        !----------------------------------------------------
        deallocate(self % data)
    end subroutine Matrix5D_Destruct

    subroutine Matrix6D_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix6D_t), intent(inout) :: self
        !----------------------------------------------------
        deallocate(self % data)
    end subroutine Matrix6D_Destruct

    subroutine Matrix7D_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix7D_t), intent(inout) :: self
        !----------------------------------------------------
        deallocate(self % data)
    end subroutine Matrix7D_Destruct
    !
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  ---------------------------------------------------------
    !  Routine for loading the MatrixNd class from a binary file
    !  ---------------------------------------------------------
    subroutine Matrix2D_Load(self, agentFile)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix2D_t), intent(inout)  :: self
        character(len=*)  , intent(in)    :: agentFile
        !----------------------------------------------------
        integer :: fd
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//"_p1.bin", form="unformatted", action="read", status="old", access='stream') 
        read(fd) self % data
        close(UNIT=fd)
    end subroutine Matrix2D_Load

    subroutine Matrix3D_Load(self, agentFile)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix3D_t), intent(inout)  :: self
        character(len=*)  , intent(in)    :: agentFile
        !----------------------------------------------------
        integer :: fd
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//"_p2.bin", form="unformatted", action="read", status="old", access='stream') 
        read(fd) self % data
        close(UNIT=fd)
    end subroutine Matrix3D_Load

    subroutine Matrix4D_Load(self, agentFile)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix4D_t), intent(inout)  :: self
        character(len=*)  , intent(in)    :: agentFile
        !----------------------------------------------------
        integer :: fd
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//"_p3.bin", form="unformatted", action="read", status="old", access='stream') 
        read(fd) self % data
        close(UNIT=fd)
    end subroutine Matrix4D_Load

    subroutine Matrix5D_Load(self, agentFile)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix5D_t), intent(inout)  :: self
        character(len=*)  , intent(in)    :: agentFile
        !----------------------------------------------------
        integer :: fd
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//"_p4.bin", form="unformatted", action="read", status="old", access='stream') 
        read(fd) self % data
        close(UNIT=fd)
    end subroutine Matrix5D_Load

    subroutine Matrix6D_Load(self, agentFile)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix6D_t), intent(inout)  :: self
        character(len=*)  , intent(in)    :: agentFile
        !----------------------------------------------------
        integer :: fd
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//"_p5.bin", form="unformatted", action="read", status="old", access='stream') 
        read(fd) self % data
        close(UNIT=fd)
    end subroutine Matrix6D_Load

    subroutine Matrix7D_Load(self, agentFile)
        implicit none
        !-arguments-------------------------------------------
        class(Matrix7D_t), intent(inout)  :: self
        character(len=*)  , intent(in)    :: agentFile
        !----------------------------------------------------
        integer :: fd
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//"_p6.bin", form="unformatted", action="read", status="old", access='stream') 
        read(fd) self % data
        close(UNIT=fd)
    end subroutine Matrix7D_Load
    !
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  -------------------------------------------------------
    !  Routine for retrieving the data from the MatrixNd class
    !  -------------------------------------------------------
    function Matrix2D_GetData(self, indices) result(val)
        !--------------------------------------
        implicit none
        !--------------------------------------
        class(Matrix2D_t)     , intent(in)    :: self
        integer               , intent(in)    :: indices(:)
        integer(kind=1)                       :: val
        !--------------------------------------
        val = self % data(indices(1), indices(2))
    end function Matrix2D_GetData

    function Matrix3D_GetData(self, indices) result(val)
        !--------------------------------------
        implicit none
        !--------------------------------------
        class(Matrix3D_t)     , intent(in)    :: self
        integer               , intent(in)    :: indices(:)
        integer(kind=1)                       :: val
        !--------------------------------------
        val = self % data(indices(1), indices(2), indices(3))
    end function Matrix3D_GetData

    function Matrix4D_GetData(self, indices) result(val)
        !--------------------------------------
        implicit none
        !--------------------------------------
        class(Matrix4D_t)     , intent(in)    :: self
        integer               , intent(in)    :: indices(:)
        integer(kind=1)                       :: val
        !--------------------------------------
        val = self % data(indices(1), indices(2), indices(3), indices(4))
    end function Matrix4D_GetData

    function Matrix5D_GetData(self, indices) result(val)
        !--------------------------------------
        implicit none
        !--------------------------------------
        class(Matrix5D_t)     , intent(in)    :: self
        integer               , intent(in)    :: indices(:)
        integer(kind=1)                       :: val
        !--------------------------------------
        val = self % data(indices(1), indices(2), indices(3), indices(4), indices(5))
    end function Matrix5D_GetData

    function Matrix6D_GetData(self, indices) result(val)
        !--------------------------------------
        implicit none
        !--------------------------------------
        class(Matrix6D_t)     , intent(in)    :: self
        integer               , intent(in)    :: indices(:)
        integer(kind=1)                       :: val
        !--------------------------------------
        val = self % data(indices(1), indices(2), indices(3), indices(4), indices(5), indices(6))
    end function Matrix6D_GetData

    function Matrix7D_GetData(self, indices) result(val)
        !--------------------------------------
        implicit none
        !--------------------------------------
        class(Matrix7D_t)     , intent(in)    :: self
        integer               , intent(in)    :: indices(:)
        integer(kind=1)                       :: val
        !--------------------------------------
        val = self % data(indices(1), indices(2), indices(3), indices(4), indices(5), indices(6), indices(7))
    end function Matrix7D_GetData

    !
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  ----------------------------------------
    !  Routine for constructing the p-adaptator
    !  ----------------------------------------
    subroutine pAdaptationAgent_Construct(self, agentFile)
        use, intrinsic :: iso_fortran_env
        implicit none
        !-arguments-------------------------------------------
        class(pAdaptationAgent_t), intent(inout) :: self
        character(len=*)         , intent(in)    :: agentFile
        !-local variables-------------------------------------
        integer          :: fd, NPoints, policySize, i, j, k, p, header, indices3(3), indices2(2)
        !----------------------------------------------------
        open(newunit = fd, FILE = agentFile//".bin", form="unformatted", action="read", status="old", access='stream')   
            READ(fd) NPoints, header
            self % smax = NPoints / 2
            READ(fd) self % pmin, header
            READ(fd) self % pmax, header
        close(UNIT=fd)
        
        policySize = self % pmax - self % pmin + 1
        allocate(self % policy(policySize))
        do i = 1, policySize
            p = self % pmin + i - 1
            if (p == 1) then
                allocate(Matrix2D_t::self % policy(i) % matrix)
            else if (p == 2) then
                allocate(Matrix3D_t::self % policy(i) % matrix)
            else if (p == 3) then
                allocate(Matrix4D_t::self % policy(i) % matrix)
            else if (p == 4) then
                allocate(Matrix5D_t::self % policy(i) % matrix)
            else if (p == 5) then
                allocate(Matrix6D_t::self % policy(i) % matrix)
            else if (p == 6) then
                allocate(Matrix7D_t::self % policy(i) % matrix)
            else 
                error stop 'The maximum polynomial order is 6 for p-adaptation with Value Iteration RL'
            end if
            call self % policy(i) % matrix % construct(NPoints)

            ! Read the policy from the file
            call self % policy(i) % matrix % load(agentFile)
        enddo

    end subroutine pAdaptationAgent_Construct
    !
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !  ----------------------------------------
    !  Routine for destructing the p-adaptator
    !  ----------------------------------------
    subroutine pAdaptationAgent_Destruct(self)
        implicit none
        !-arguments-------------------------------------------
        class(pAdaptationAgent_t), intent(inout) :: self
        !-local variables-------------------------------------
        integer :: policySize, i
        !----------------------------------------------------     
        policySize = self % pmax - self % pmin + 1
        do i = 1, policySize
            call self % policy(i) % matrix % destruct()
        enddo

        deallocate(self % policy)

    end subroutine pAdaptationAgent_Destruct

end module ReinforcementLearning