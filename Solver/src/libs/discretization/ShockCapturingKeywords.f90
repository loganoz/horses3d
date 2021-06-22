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
   character(len=*), parameter :: SC_KEY            = "enable shock-capturing"
   character(len=*), parameter :: SC_SENSOR_KEY     = "shock-capturing sensor"
   character(len=*), parameter :: SC_METHOD_KEY     = "shock-capturing method"
   character(len=*), parameter :: SC_VARIABLE_KEY   = "sensor variable"
   character(len=*), parameter :: SC_LOW_THRES_KEY  = "sensor lower limit"
   character(len=*), parameter :: SC_HIGH_THRES_KEY = "sensor higher limit"
   character(len=*), parameter :: SC_THRES_1_KEY    = "sensor threshold 1"
   character(len=*), parameter :: SC_THRES_2_KEY    = "sensor threshold 2"
   character(len=*), parameter :: SC_AVIS_KEY       = "shock-capturing art. viscosity"
!
!  Sensor types
!  ------------
   character(len=*), parameter :: SC_MODAL_KEY = "modal"
!
!  Shock-capturing methods
!  -----------------------
   character(len=*), parameter :: SC_SVV_KEY   = "svv"
   character(len=*), parameter :: SC_SSFV_KEY  = "ssfv-svv"
   character(len=*), parameter :: SC_NOSVV_KEY = "no svv"
!
!  Sensed variables
!  ----------------
   character(len=*), parameter :: SC_RHO_KEY  = "rho"
   character(len=*), parameter :: SC_RHOU_KEY = "rhou"
   character(len=*), parameter :: SC_RHOV_KEY = "rhov"
   character(len=*), parameter :: SC_RHOW_KEY = "rhow"
   character(len=*), parameter :: SC_RHOE_KEY = "rhoe"
   character(len=*), parameter :: SC_U_KEY    = "u"
   character(len=*), parameter :: SC_V_KEY    = "v"
   character(len=*), parameter :: SC_W_KEY    = "w"
   character(len=*), parameter :: SC_P_KEY    = "p"
   character(len=*), parameter :: SC_RHOP_KEY = "rhop"

end module ShockCapturingKeywords
