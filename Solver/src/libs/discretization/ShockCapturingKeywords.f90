!
!//////////////////////////////////////////////////////
!
!   @File:    ShockCapturingKeywords.f90
!   @Author:  Andrés Mateo (andres.mgabin@upm.es)
!   @Created: Thu Jun 17 2021
!   @Last revision date: Thu Jun 17 2021
!   @Last revision author: Andrés Mateo (andres.mgabin@upm.es)
!
!//////////////////////////////////////////////////////
!
module ShockCapturingKeywords
!
!  Keywords
!  --------
   character(len=*), parameter :: SC_KEY          = "enable shock-capturing"
   character(len=*), parameter :: SC_SENSOR_KEY   = "shock-capturing sensor"
   character(len=*), parameter :: SC_VARIABLE_KEY = "shock-capturing variable"
   character(len=*), parameter :: SC_METHOD       = "shock-capturing method"

   character(len=*), parameter :: SC_MODAL_KEY = "modal"
   character(len=*), parameter :: SC_RHO_KEY   = "rho"
   character(len=*), parameter :: SC_RHOU_KEY  = "rhou"
   character(len=*), parameter :: SC_RHOV_KEY  = "rhov"
   character(len=*), parameter :: SC_RHOW_KEY  = "rhow"
   character(len=*), parameter :: SC_RHOE_KEY  = "rhoe"
   character(len=*), parameter :: SC_U_KEY     = "u"
   character(len=*), parameter :: SC_V_KEY     = "v"
   character(len=*), parameter :: SC_W_KEY     = "w"
   character(len=*), parameter :: SC_P_KEY     = "p"
   character(len=*), parameter :: SC_RHOP_KEY  = "rhop"
!
!  Sensor types
!  ------------
   enum, bind(C)
      enumerator :: SC_MODAL
   end enum
!
!  Sensed variables
!  ----------------
   enum, bind(C)
      enumerator :: SC_RHO, SC_RHOU, SC_RHOV, SC_RHOW, SC_RHOE
      enumerator :: SC_U, SC_V, SC_W, SC_P
      enumerator :: SC_RHOP
   end enum 

end module ShockCapturingKeywords
