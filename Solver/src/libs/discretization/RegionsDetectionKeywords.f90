module RegionsDetectionKeywords

    character(len=*), parameter :: VIS_KEY                = "enable viscous regions detection" 
    character(len=*), parameter :: VIS_SENSOR_KEY         = "viscous regions sensor" 
    character(len=*), parameter :: VIS_GMM_VAL            = "gmm based viscous regions detection" 
    character(len=*), parameter :: VIS_NUM_CLUSTERS_KEY   = "sensor number of clusters" 
    character(len=*), parameter :: VIS_ITER_JUMP          = "interval of iterations to activate the viscous sensor"
    character(len=*), parameter :: VIS_SENSOR_INERTIA_KEY = "sensor min. timesteps"
    character(len=*), parameter :: VIS_SENSOR_ELEM_CLUST  = "method to cluster mesh elements" 
    character(len=*), parameter :: VIS_TO_HYBRID          = "set off viscous fluxes in the inviscid region" 
    character(len=*), parameter :: VIS_TO_ADAPT           = "perform p-adaptation with viscous sensor"
    character(len=*), parameter :: VIS_ITER_MIN           = "minimum number of iterations to activate the sensor" 
    
    integer, parameter :: VIS_GMM_ID          = 1
    character(len=*), parameter :: VIS_INV_VAL        = "Q_s +R_s +Q_omega"
    integer, parameter :: VIS_INV_ID        = 1
end module RegionsDetectionKeywords