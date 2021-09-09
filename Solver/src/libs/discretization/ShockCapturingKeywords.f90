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
   character(len=*), parameter :: SC_SENSOR_KEY     = "shock sensor"
   character(len=*), parameter :: SC_METHOD_KEY     = "shock method"
   character(len=*), parameter :: SC_VISC_FLUX_KEY  = "shock viscous flux"
   character(len=*), parameter :: SC_UPDATE_KEY     = "shock update viscosity"
   character(len=*), parameter :: SC_MU1_KEY        = "shock mu 1"
   character(len=*), parameter :: SC_ALPHA1_KEY     = "shock alpha 1"
   character(len=*), parameter :: SC_MU2_KEY        = "shock mu 2"
   character(len=*), parameter :: SC_ALPHA2_KEY     = "shock alpha 2"
   character(len=*), parameter :: SC_ALPHA_MU_KEY   = "shock alpha/mu"
   character(len=*), parameter :: SC_VARIABLE_KEY   = "sensor variable"
   character(len=*), parameter :: SC_LOW_THRES_KEY  = "sensor lower limit"
   character(len=*), parameter :: SC_HIGH_THRES_KEY = "sensor higher limit"
   character(len=*), parameter :: SC_THRES_1_KEY    = "sensor threshold 1"
   character(len=*), parameter :: SC_THRES_2_KEY    = "sensor threshold 2"
!
!  Sensor types
!  ------------
   character(len=*), parameter :: SC_ZERO_VAL    = "zeros"
   character(len=*), parameter :: SC_ONE_VAL     = "ones"
   character(len=*), parameter :: SC_GRADRHO_VAL = "grad rho"
   character(len=*), parameter :: SC_MODAL_VAL   = "modal"

   integer, parameter :: SC_ZERO_ID    = 1
   integer, parameter :: SC_ONE_ID     = 2
   integer, parameter :: SC_GRADRHO_ID = 3
   integer, parameter :: SC_MODAL_ID   = 4
!
!  Shock-capturing methods
!  -----------------------
   character(len=*), parameter :: SC_NOSVV_VAL = "non-filtered"
   character(len=*), parameter :: SC_SVV_VAL   = "svv"
   character(len=*), parameter :: SC_SSFV_VAL  = "ssfv-svv"
!
!  Artificial viscosity fluxes
!  ---------------------------
   character(len=*), parameter :: SC_PHYS_VAL = "physical"
   character(len=*), parameter :: SC_GP_VAL   = "guermond-popov"

   integer, parameter :: SC_PHYS_ID = 1
   integer, parameter :: SC_GP_ID   = 2
!
!  Viscosity update method
!  -----------------------
   character(len=*), parameter :: SC_CONST_VAL  = "constant"
   character(len=*), parameter :: SC_SENSOR_VAL = "sensor"
   character(len=*), parameter :: SC_SMAG_VAL   = "smagorinsky"

   integer, parameter :: SC_CONST_ID  = 1
   integer, parameter :: SC_SENSOR_ID = 2
   integer, parameter :: SC_SMAG_ID   = 3
!
!  Sensed variables
!  ----------------
   character(len=*), parameter :: SC_RHO_VAL  = "rho"
   character(len=*), parameter :: SC_RHOU_VAL = "rhou"
   character(len=*), parameter :: SC_RHOV_VAL = "rhov"
   character(len=*), parameter :: SC_RHOW_VAL = "rhow"
   character(len=*), parameter :: SC_RHOE_VAL = "rhoe"
   character(len=*), parameter :: SC_U_VAL    = "u"
   character(len=*), parameter :: SC_V_VAL    = "v"
   character(len=*), parameter :: SC_W_VAL    = "w"
   character(len=*), parameter :: SC_P_VAL    = "p"
   character(len=*), parameter :: SC_RHOP_VAL = "rhop"

   integer, parameter :: SC_RHO_ID  = 1
   integer, parameter :: SC_RHOU_ID = 2
   integer, parameter :: SC_RHOV_ID = 3
   integer, parameter :: SC_RHOW_ID = 4
   integer, parameter :: SC_RHOE_ID = 5
   integer, parameter :: SC_U_ID    = 6
   integer, parameter :: SC_V_ID    = 7
   integer, parameter :: SC_W_ID    = 8
   integer, parameter :: SC_P_ID    = 9
   integer, parameter :: SC_RHOP_ID = 10

end module ShockCapturingKeywords
