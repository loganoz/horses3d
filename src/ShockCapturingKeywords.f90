module ShockCapturingKeywords
!
!  Keywords
!  --------
   character(len=*), parameter :: SC_KEY                = "enable shock-capturing"
   character(len=*), parameter :: SC_SENSOR_KEY         = "shock sensor"
   character(len=*), parameter :: SC_METHOD1_KEY        = "shock first method"
   character(len=*), parameter :: SC_METHOD2_KEY        = "shock second method"
   character(len=*), parameter :: SC_VISC_FLUX1_KEY     = "shock viscous flux 1"
   character(len=*), parameter :: SC_VISC_FLUX2_KEY     = "shock viscous flux 2"
   character(len=*), parameter :: SC_UPDATE_KEY         = "shock update strategy"
   character(len=*), parameter :: SC_MU1_KEY            = "shock mu 1"
   character(len=*), parameter :: SC_ALPHA1_KEY         = "shock alpha 1"
   character(len=*), parameter :: SC_MU2_KEY            = "shock mu 2"
   character(len=*), parameter :: SC_ALPHA2_KEY         = "shock alpha 2"
   character(len=*), parameter :: SC_ALPHA_MU_KEY       = "shock alpha/mu"
   character(len=*), parameter :: SC_VARIABLE_KEY       = "sensor variable"
   character(len=*), parameter :: SC_LOW_THRES_KEY      = "sensor lower limit"
   character(len=*), parameter :: SC_HIGH_THRES_KEY     = "sensor higher limit"
   character(len=*), parameter :: SC_TE_NMIN_KEY        = "sensor te min n"
   character(len=*), parameter :: SC_TE_DELTA_KEY       = "sensor te delta n"
   character(len=*), parameter :: SC_TE_DTYPE_KEY       = "sensor te derivative"
   character(len=*), parameter :: SC_NUM_CLUSTERS_KEY   = "sensor number of clusters"
   character(len=*), parameter :: SC_SENSOR_INERTIA_KEY = "sensor min. timesteps"
!
!  Sensor types
!  ------------
   character(len=*), parameter :: SC_ZERO_VAL          = "zeros"
   character(len=*), parameter :: SC_ONE_VAL           = "ones"
   character(len=*), parameter :: SC_INTEGRAL_VAL      = "integral"
   character(len=*), parameter :: SC_INTEGRAL_SQRT_VAL = "integral with sqrt"
   character(len=*), parameter :: SC_MODAL_VAL         = "modal"
   character(len=*), parameter :: SC_TE_VAL            = "truncation error"
   character(len=*), parameter :: SC_ALIAS_VAL         = "aliasing error"
   character(len=*), parameter :: SC_GMM_VAL           = "gmm"

   integer, parameter :: SC_ZERO_ID          = 1
   integer, parameter :: SC_ONE_ID           = 2
   integer, parameter :: SC_INTEGRAL_ID      = 3
   integer, parameter :: SC_INTEGRAL_SQRT_ID = 4
   integer, parameter :: SC_MODAL_ID         = 5
   integer, parameter :: SC_TE_ID            = 6
   integer, parameter :: SC_ALIAS_ID         = 7
   integer, parameter :: SC_GMM_ID           = 8
!
!  Shock-capturing methods
!  -----------------------
   character(len=*), parameter :: SC_NO_VAL    = "none"
   character(len=*), parameter :: SC_NOSVV_VAL = "non-filtered"
   character(len=*), parameter :: SC_SVV_VAL   = "svv"
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
   character(len=*), parameter :: SC_RHO_VAL         = "rho"
   character(len=*), parameter :: SC_RHOU_VAL        = "rhou"
   character(len=*), parameter :: SC_RHOV_VAL        = "rhov"
   character(len=*), parameter :: SC_RHOW_VAL        = "rhow"
   character(len=*), parameter :: SC_RHOE_VAL        = "rhoe"
   character(len=*), parameter :: SC_U_VAL           = "u"
   character(len=*), parameter :: SC_V_VAL           = "v"
   character(len=*), parameter :: SC_W_VAL           = "w"
   character(len=*), parameter :: SC_P_VAL           = "p"
   character(len=*), parameter :: SC_RHOP_VAL        = "rhop"
   character(len=*), parameter :: SC_RHO_GRAD_VAL    = "grad rho"
   character(len=*), parameter :: SC_DIVV_VAL        = "div v"
   character(len=*), parameter :: SC_DIVV_P_GRAD_VAL = "div v + grad p"

   integer, parameter :: SC_RHO_ID         = 1
   integer, parameter :: SC_RHOU_ID        = 2
   integer, parameter :: SC_RHOV_ID        = 3
   integer, parameter :: SC_RHOW_ID        = 4
   integer, parameter :: SC_RHOE_ID        = 5
   integer, parameter :: SC_U_ID           = 6
   integer, parameter :: SC_V_ID           = 7
   integer, parameter :: SC_W_ID           = 8
   integer, parameter :: SC_P_ID           = 9
   integer, parameter :: SC_RHOP_ID        = 10
   integer, parameter :: SC_RHO_GRAD_ID    = 11
   integer, parameter :: SC_DIVV_ID        = 12
   integer, parameter :: SC_DIVV_P_GRAD_ID = 13
!
!  Derivative types
!  ----------------
   character(len=*), parameter :: SC_ISOLATED_KEY     = "isolated"
   character(len=*), parameter :: SC_NON_ISOLATED_KEY = "non-isolated"

end module ShockCapturingKeywords