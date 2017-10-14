!
!//////////////////////////////////////////////////////
!
!   @File:    OutputVariables.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sat Oct 14 20:44:38 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module OutputVariables
   use SMConstants

   private
!   public   

   integer, parameter   :: STR_VAR_LEN = 8
!
!  ***************************
!  Variables without gradients
!  ***************************
!
   integer, parameter :: Q_V    = 0
   integer, parameter :: RHO_V  = 1
   integer, parameter :: U_V    = 2
   integer, parameter :: V_V    = 3
   integer, parameter :: W_V    = 4
   integer, parameter :: P_V    = 5
   integer, parameter :: T_V    = 6
   integer, parameter :: Mach_V = 7
   integer, parameter :: S_V    = 8
   integer, parameter :: Vabs_V = 9
   integer, parameter :: Vvec_V = 10
   integer, parameter :: Ht_V   = 11

   character(len = STR_VAR_LEN), parameter  :: QKey    = "Q"
   character(len = STR_VAR_LEN), parameter  :: RHOKey  = "rho"
   character(len = STR_VAR_LEN), parameter  :: UKey    = "u"
   character(len = STR_VAR_LEN), parameter  :: VKey    = "v"
   character(len = STR_VAR_LEN), parameter  :: WKey    = "w"
   character(len = STR_VAR_LEN), parameter  :: PKey    = "p"
   character(len = STR_VAR_LEN), parameter  :: TKey    = "T"
   character(len = STR_VAR_LEN), parameter  :: MachKey = "Mach"
   character(len = STR_VAR_LEN), parameter  :: SKey    = "s"
   character(len = STR_VAR_LEN), parameter  :: VabsKey = "Vabs"
   character(len = STR_VAR_LEN), parameter  :: VvecKey = "V"
   character(len = STR_VAR_LEN), parameter  :: HtKey   = "Ht"

   character(len=STR_VAR_LEN), dimension(12), parameter  :: variableNames = (/ QKey, RHOKey, UKey, VKey, WKey, &
                                                                            PKey, TKey, MachKey, SKey, VabsKey, &
                                                                            VvecKey, HtKey /)
                                                               

   integer                :: no_of_outputVariables
   integer, allocatable   :: outputVariableNames(:)

   contains
   
      integer function outputVariablesForVariable(iVar)
!
!        ************************************************
!           This subroutine specifies if a variable
!           implies more than one variable, e.g.,
!           Q = [rho,rhou,rhov,rhow,rhoe], or
!           V = [u,v,w]
!        ************************************************
!
         implicit none
         integer,    intent(in)     :: iVar
   
         select case(iVar)
   
         case(Q_V)
            outputVariablesForVariable = 5

         case(Vvec_V)
            outputVariablesForVariable = 3

         case default
            outputVariablesForVariable = 1
      
         end select

      end function outputVariablesForVariable

end module OutputVariables
