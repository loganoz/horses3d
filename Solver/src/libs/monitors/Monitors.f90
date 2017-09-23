module MonitorsClass
   use SMConstants
   use Physics
   use QuadMeshClass
   implicit none
!
   private
   public      Monitor_t , ConstructMonitors
!
!                                ***********************
   integer, parameter         :: BUFFER_SIZE      = 1000
   integer, parameter         :: STR_LEN_MONITORS = 128
   integer, parameter         :: MONITOR_LENGTH   = 10
!                                ***********************
!
!  **********************
!  Probe class definition
!  **********************
!
   type Probe_t
      logical                         :: active
      integer                         :: eID
      integer                         :: ID
      real(kind=RP)                   :: x(NDIM)
      real(kind=RP)                   :: xi, eta
      real(kind=RP)                   :: values(BUFFER_SIZE)
      real(kind=RP), allocatable      :: l_xi(:) , l_eta(:)
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS) :: variable
      contains
         procedure   :: Initialization => Probe_Initialization
         procedure   :: Update         => Probe_Update
         procedure   :: WriteLabel     => Probe_WriteLabel
         procedure   :: WriteValues    => Probe_WriteValue
         procedure   :: WriteToFile    => Probe_WriteToFile
   end type Probe_t
!
!  ********************************
!  Surface monitor class definition
!  ********************************
!
   type SurfaceMonitor_t

   end type SurfaceMonitor_t
!
!  *******************************
!  Volume monitor class definition
!  *******************************
!
   type VolumeMonitor_t

   end type VolumeMonitor_t

!
!
!  ************************
!  Monitor class definition
!  ************************
!  
   type Monitor_t
      integer                                  :: no_of_probes
      integer                                  :: no_of_surfaceMonitors
      integer                                  :: no_of_volumeMonitors
      integer                                  :: bufferLine
      real(kind=RP)                            :: t(BUFFER_SIZE)
      class ( Probe_t          ) , allocatable :: probes          ( : )
      class ( SurfaceMonitor_t ) , allocatable :: surfaceMonitors ( : )
      class ( VolumeMonitor_t  ) , allocatable :: volumeMonitors  ( : )
      contains
         procedure   :: WriteLabel => Monitor_WriteLabel
         procedure   :: WriteUnderlines => Monitor_WriteUnderlines
         procedure   :: WriteValues    => Monitor_WriteValues
         procedure   :: UpdateValues => Monitor_UpdateValues
         procedure   :: WriteToFile => Monitor_WriteToFile
   end type Monitor_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      function ConstructMonitors( mesh ) result(Monitors)
         use Setup_Class
         use ParamfileIO
         implicit none
         type(Monitor_t)         :: Monitors
         class(QuadMesh_t)       :: mesh
!        -------------------------------------------------
         integer                         :: fID , io
         integer                         :: i
         character(len=STR_LEN_MONITORS) :: line
!
!        Search in case file for probes, surface monitors, and volume monitors
!        ---------------------------------------------------------------------
         Monitors % no_of_probes          = 0
         Monitors % no_of_surfaceMonitors = 0
         Monitors % no_of_volumeMonitors  = 0

         open ( newunit = fID , file = trim(Setup % case_file) , status = "old" , action = "read" )

readloop:do 
            read ( fID , '(A)' , iostat = io ) line

            if ( io .lt. 0 ) then
!
!              End of file
!              -----------
               line = ""
               exit readloop

            elseif ( io .gt. 0 ) then
!
!              Error
!              -----
               print*, "Error reading case file, in line 96 of Monitors.f90"
               stop "Stopped."

            else
!
!              Succeeded
!              ---------
               line = getSquashedLine( line )

               if ( index ( line , '#defineprobe' ) .gt. 0 ) then
                  Monitors % no_of_probes = Monitors % no_of_probes + 1

               elseif ( index ( line , '#definesurfacemonitor' ) .gt. 0 ) then
                  Monitors % no_of_surfaceMonitors = Monitors % no_of_surfaceMonitors + 1

               elseif ( index ( line , '#definevolumemonitor' ) .gt. 0 ) then
                  Monitors % no_of_volumeMonitors = Monitors % no_of_volumeMonitors + 1

               end if
               
            end if

         end do readloop
!
!        Close case file
!        ---------------
         close ( fID )
!
!        Allocate monitors
!        -----------------
         allocate ( Monitors % probes          ( Monitors % no_of_probes          )  ) 
         allocate ( Monitors % surfaceMonitors ( Monitors % no_of_surfaceMonitors )  ) 
         allocate ( Monitors % volumeMonitors  ( Monitors % no_of_volumeMonitors  )  ) 
!
!        Initialize
!        ----------
         do i = 1 , Monitors % no_of_probes
            call Monitors % probes(i) % Initialization(mesh , i)
         end do

         Monitors % bufferLine = 0

      end function ConstructMonitors

      subroutine Monitor_WriteLabel ( self )
         implicit none
         class(Monitor_t)              :: self
         integer                       :: i 

         do i = 1 , self % no_of_probes
            call self % probes(i) % WriteLabel
         end do

         do i = 1 , self % no_of_surfaceMonitors
!            call self % probe(i) % WriteLab
         end do

         do i = 1 , self % no_of_volumeMonitors
!            call self % probe(i) % WriteLabel
         end do
         

      end subroutine Monitor_WriteLabel

      subroutine Monitor_WriteUnderlines( self ) 
         implicit none
         class(Monitor_t)             :: self
         integer                      :: i
         character(len=MONITOR_LENGTH), parameter :: dashes = "----------"
   
         do i = 1 , self % no_of_probes
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % probes(i) % monitorName ) + 2 ) )
         end do

      end subroutine Monitor_WriteUnderlines

      subroutine Monitor_WriteValues ( self )
         implicit none
         class(Monitor_t)           :: self
         integer                    :: i

         do i = 1 , self % no_of_probes
            call self % probes(i) % WriteValues ( self % bufferLine )
         end do

      end subroutine Monitor_WriteValues

      subroutine Monitor_UpdateValues ( self , mesh , t )
         implicit none
         class(Monitor_t)              :: self
         class(QuadMesh_t)             :: mesh
         real(kind=RP)                 :: t
         integer                       :: i 

         self % bufferLine = self % bufferLine + 1

         self % t ( self % bufferLine ) = t

         do i = 1 , self % no_of_probes
            call self % probes(i) % Update( mesh , self % bufferLine )

         end do

      end subroutine Monitor_UpdateValues

      subroutine Monitor_WriteToFile ( self , force) 
         implicit none
         class(Monitor_t)        :: self
         logical, optional       :: force
         integer                 :: i 
         logical                 :: forceVal

         if ( present ( force ) ) then
            forceVal = force

         else
            forceVal = .false.

         end if

         if ( forceVal ) then 
            do i = 1 , self % no_of_probes
               call self % probes(i) % WriteToFile ( self % t , self % bufferLine )

            end do

            self % bufferLine = 0

         else

            if ( self % bufferLine .eq. BUFFER_SIZE ) then
               do i = 1 , self % no_of_probes
                  call self % probes(i) % WriteToFile ( self % t , BUFFER_SIZE ) 

               end do

               self % bufferLine = 0
   
            end if


         end if

      end subroutine Monitor_WriteToFile
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           PROBE ROUTINES
!           --------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Probe_Initialization( self , mesh , ID ) 
         use ParamfileIO
         use Setup_Class
         implicit none
         class(Probe_t)          :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: ID
!        ----------------------------------------------
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         real(kind=RP), allocatable       :: x,y
         integer                          :: pos
         integer                          :: fID

         self % ID = ID

         write(in_label , '(A,I0)') "#define probe " , self % ID
         
         call readValueInRegion ( trim ( Setup % case_file )  , "Name"       , self % monitorName , in_label , "# end" ) 
         call readValueInRegion ( trim ( Setup % case_file )  , "x position" , x                  , in_label , "# end" ) 
         call readValueInRegion ( trim ( Setup % case_file )  , "Variable"   , self % variable    , in_label , "# end" ) 
         call readValueInRegion ( trim ( Setup % case_file )  , "y position" , y                  , in_label , "# end" ) 

         call mesh % findElementWithCoords([x,y] , self % eID , self % xi , self % eta )

         if ( self % eID .eq. -1 ) then
            self % active = .false.
            return

         else
            self % active = .true.

         end if

         self % l_xi  = mesh % elements(self % eID) % spA % lj( self % xi  )
         self % l_eta = mesh % elements(self % eID) % spA % lj( self % eta )

         fileName = trim(Setup % solution_file)
   
         pos = index(trim(fileName) , '.HiORst' )
            
         write( self % fileName , '(A,A,A,A)') fileName(1:pos-1) , "." , trim(self % monitorName) , ".probe"  
!
!        Create file
!        -----------
         open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
         write( fID , '(A24,2X,A24)' ) "time" , trim(self % variable)
         close ( fID ) 
              
      end subroutine Probe_Initialization

      subroutine Probe_Update ( self , mesh , bufferPosition )
         use MatrixOperations
         implicit none
         class(Probe_t)          :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: bufferPosition
!        ------------------------------------------------------
         real(kind=RP)           :: rho , rhou , rhov , rhoe
         real(kind=RP)           :: rhot , rhout , rhovt , rhoet
      
         associate ( N => mesh % elements( self % ID ) % spA % N )

         call BilinearForm ( mesh % elements ( self % eID )  % Q    ( 0:N , 0:N , IRHO  )  , self % l_xi , self % l_eta , rho   ) 
         call BilinearForm ( mesh % elements ( self % eID )  % Q    ( 0:N , 0:N , IRHOU )  , self % l_xi , self % l_eta , rhou  ) 
         call BilinearForm ( mesh % elements ( self % eID )  % Q    ( 0:N , 0:N , IRHOV )  , self % l_xi , self % l_eta , rhov  ) 
         call BilinearForm ( mesh % elements ( self % eID )  % Q    ( 0:N , 0:N , IRHOE )  , self % l_xi , self % l_eta , rhoe  ) 

         call BilinearForm ( mesh % elements ( self % eID )  % QDot ( 0:N , 0:N , IRHO  )  , self % l_xi , self % l_eta , rhot  ) 
         call BilinearForm ( mesh % elements ( self % eID )  % QDot ( 0:N , 0:N , IRHOU )  , self % l_xi , self % l_eta , rhout ) 
         call BilinearForm ( mesh % elements ( self % eID )  % QDot ( 0:N , 0:N , IRHOV )  , self % l_xi , self % l_eta , rhovt ) 
         call BilinearForm ( mesh % elements ( self % eID )  % QDot ( 0:N , 0:N , IRHOE )  , self % l_xi , self % l_eta , rhoet ) 

         end associate 

         select case ( trim( self % variable ) )

            case ("rho")
               self % values(bufferPosition) = rho * refValues % rho
 
            case ("rhou")
               self % values(bufferPosition) = rhou * refValues % rho * refValues % a
 
            case ("rhov")
               self % values(bufferPosition) = rhov * refValues % rho * refValues % a
 
            case ("rhoe")
               self % values(bufferPosition) = rhoe * refValues % rho * refValues % p
 
            case ("rhot")
               self % values(bufferPosition) = rhot * refValues % rho / refValues % tc
 
            case ("rhout")
               self % values(bufferPosition) = rhout * refValues % rho * refValues % a / refValues % tc
 
            case ("rhovt")
               self % values(bufferPosition) = rhovt * refValues % rho * refValues % a / refValues % tc
 
            case ("rhoet")
               self % values(bufferPosition) = rhoet * refValues % rho * refValues % p / refValues % tc
 
            case ("u")
               self % values(bufferPosition) = rhou / rho * refValues % a
 
            case ("v")
               self % values(bufferPosition) = rhov / rho * refValues % a
    
            case ("p")
               self % values(bufferPosition) = Thermodynamics % gm1 * ( rhoe - 0.5_RP * ( rhou * rhou + rhov * rhov ) / rho ) * refValues % p
       
            case ("Mach")
               self % values(bufferPosition) = sqrt(rhou * rhou + rhov * rhov) / rho / sqrt(Thermodynamics % gamma)
 
         end select                        

      end subroutine Probe_Update

      subroutine Probe_WriteLabel ( self )
         implicit none
         class(Probe_t)             :: self

         write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))

      end subroutine Probe_WriteLabel
   
      subroutine Probe_WriteValue ( self , bufferLine ) 
         implicit none
         class(Probe_t)             :: self
         integer                    :: bufferLine

         write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values ( bufferLine ) 

      end subroutine Probe_WriteValue 

      subroutine Probe_WriteToFile ( self , t , no_of_lines)
         implicit none  
         class(Probe_t)             :: self
         real(kind=RP)              :: t(:)
         integer                    :: no_of_lines
!        -------------------------------------------
         integer                    :: i
         integer                    :: fID

         open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
         
         do i = 1 , no_of_lines
            write( fID , '(ES24.16,2X,ES24.16)' ) t(i) , self % values(i)

         end do
        
         close ( fID )
      
      end subroutine Probe_WriteToFile

end module MonitorsClass
