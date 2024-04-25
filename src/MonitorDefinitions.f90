module MonitorDefinitions
   implicit none
   private
   public BUFFER_SIZE_DEFAULT, STR_LEN_MONITORS
   public MONITOR_LENGTH, VOLUME_UNDEFINED, VOLUME_INTEGRAL
   public BUFFER_SIZE

   integer, parameter         :: BUFFER_SIZE_DEFAULT = 100
   integer                    :: BUFFER_SIZE = BUFFER_SIZE_DEFAULT
   integer, parameter         :: STR_LEN_MONITORS  = 128
   integer, parameter         :: MONITOR_LENGTH    = 10
   integer, parameter         :: VOLUME_UNDEFINED  = 0
   integer, parameter         :: VOLUME_INTEGRAL   = 1

end module MonitorDefinitions