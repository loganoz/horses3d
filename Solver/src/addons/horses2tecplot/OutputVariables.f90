!
!//////////////////////////////////////////////////////
!
!   @File:    OutputVariables.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sat Oct 14 20:44:38 2017
!   @Last revision date: Sun Oct 15 23:09:45 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fffbd669a72f9b212fc48c02a022b43cfc365ef2
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module OutputVariables
!
!//////////////////////////////////////////////////////////////////////////////////
!
!        This module selects user-defined output variables for the 
!     tecplot file. The user may add extra variables, this can be done by
!     following the steps:
!
!           * Increase NO_OF_VARIABLES in 1.
!           * Adding a new ID for the variable NEWVAR_V (following the given order).
!           * Adding a new key for the variable, NEWVARKey.
!           * Adding the new key to the keys array.
!           * Adding the "select case" procedure that computes the output
!                 variable.
!
!//////////////////////////////////////////////////////////////////////////////////
!
   use SMConstants
   use PhysicsStorage

   private
   public   no_of_outputVariables
   public   getOutputVariables, ComputeOutputVariables, getOutputVariablesLabel

   integer, parameter   :: STR_VAR_LEN = 8
   integer, parameter   :: NO_OF_VARIABLES = 16
!
!  ***************************
!  Variables without gradients
!  ***************************
!
   integer, parameter :: Q_V    = 1
   integer, parameter :: RHO_V  = 2
   integer, parameter :: U_V    = 3
   integer, parameter :: V_V    = 4
   integer, parameter :: W_V    = 5
   integer, parameter :: P_V    = 6
   integer, parameter :: T_V    = 7
   integer, parameter :: Mach_V = 8
   integer, parameter :: S_V    = 9
   integer, parameter :: Vabs_V = 10
   integer, parameter :: Vvec_V = 11
   integer, parameter :: Ht_V   = 12
   integer, parameter :: RHOU_V = 13
   integer, parameter :: RHOV_V = 14
   integer, parameter :: RHOW_V = 15
   integer, parameter :: RHOE_V = 16

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
   character(len = STR_VAR_LEN), parameter  :: RHOUKey = "rhou"
   character(len = STR_VAR_LEN), parameter  :: RHOVKey = "rhov"
   character(len = STR_VAR_LEN), parameter  :: RHOWKey = "rhow"
   character(len = STR_VAR_LEN), parameter  :: RHOEKey = "rhoe"

   character(len=STR_VAR_LEN), dimension(16), parameter  :: variableNames = (/ QKey, RHOKey, UKey, VKey, WKey, &
                                                                            PKey, TKey, MachKey, SKey, VabsKey, &
                                                                            VvecKey, HtKey, RHOUKey, RHOVKey, RHOWKey, &
                                                                            RHOEKey /)
                                                               

   integer                :: no_of_outputVariables
   integer, allocatable   :: outputVariableNames(:)

   contains
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine getOutputVariables(flag_name)
         implicit none
         character(len=*), intent(in)  :: flag_name
!
!        ---------------
!        Local variables
!        ---------------
!
         logical                       :: flagPresent
         integer                       :: pos, pos2, i, preliminarNoOfVariables
         integer, allocatable          :: preliminarVariables(:)
         character(len=LINE_LENGTH)    :: flag   
         character(len=STR_VAR_LEN)   :: inputVar

         flagPresent = .false.

         do i = 1, command_argument_count()
            call get_command_argument(i, flag)
            pos = index(trim(flag),trim(flag_name))

            if ( pos .ne. 0 ) then
               flagPresent = .true.
               exit
            end if 
         end do
!
!        ***********************************************************
!              Read the variables. They are first loaded into
!           a "preliminar" variables, prior to introduce them in 
!           the real outputVariables. This is because some of
!           the output variables lead to multiple variables (e.g.
!           Q or V).
!        ***********************************************************
!
         if ( .not. flagPresent ) then
!
!           Default: export Q
!           -------
            preliminarNoOfVariables = 1
            allocate( preliminarVariables(preliminarNoOfVariables) )
            preliminarVariables(1) = Q_V
   
         else
!
!           Get the variables from the command argument
!           -------------------------------------------
            pos = index(trim(flag),"=")

            if ( pos .eq. 0 ) then
               print*, 'Missing "=" operator in --output-variables flag'
               errorMessage(STD_OUT)
               stop
            end if
!
!           Prepare to read the variable names
!           ----------------------------------
            preliminarNoOfVariables = getNoOfCommas(trim(flag)) + 1
            allocate( preliminarVariables(preliminarNoOfVariables) )
            
            if ( preliminarNoOfVariables .eq. 1 ) then
               read(flag(pos+1:len_trim(flag)),*) inputVar
               preliminarVariables(1) = outputVariableForName(trim(inputVar))
            else
               do i = 1, preliminarNoOfVariables-1
                  pos2 = index(trim(flag(pos+1:)),",") + pos
                  read(flag(pos+1:pos2),*) inputVar
                  preliminarVariables(i) = outputVariableForName(trim(inputVar))
                  pos = pos2
               end do
            
               pos = index(trim(flag),",",BACK=.true.)
               preliminarVariables(preliminarNoOfVariables) = outputVariableForName(flag(pos+1:))
               
            end if
         end if
!
!        *******************************************************
!        Convert the preliminary variables into output variables
!        *******************************************************
!
         no_of_outputVariables = 0    
         
         do i = 1, preliminarNoOfVariables
            no_of_outputVariables = no_of_outputVariables + outputVariablesForVariable(preliminarVariables(i))
         end do

         allocate( outputVariableNames(no_of_outputVariables) )

         pos = 1
         do i = 1, preliminarNoOfVariables
            pos2 = pos + outputVariablesForVariable(preliminarVariables(i))
            call outputVariablesForPreliminarVariable(preliminarVariables(i), outputVariableNames(pos:pos2) )
      
            pos = pos + outputVariablesForVariable(preliminarVariables(i))
         end do

      end subroutine getOutputVariables

      subroutine ComputeOutputVariables(N, Q, output, refs)
         use SolutionFile
         implicit none
         integer, intent(in)     :: N(3)
         real(kind=RP), intent(in)  :: Q(0:N(1),0:N(2),0:N(3),1:5)
         real(kind=RP), intent(out) :: output(0:N(1),0:N(2),0:N(3),1:no_of_outputVariables)
         real(kind=RP), intent(in)  :: refs(NO_OF_SAVED_REFS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: var, i, j, k
         
         do var = 1, no_of_outputVariables

            select case (outputVariableNames(var))

            case(RHO_V)
               output(:,:,:,var) = Q(:,:,:,IRHO)
            case(U_V)
               output(:,:,:,var) = Q(:,:,:,IRHOU) / Q(:,:,:,IRHO)

            case(V_V)
               output(:,:,:,var) = Q(:,:,:,IRHOV) / Q(:,:,:,IRHO)

            case(W_V)
               output(:,:,:,var) = Q(:,:,:,IRHOW) / Q(:,:,:,IRHO)

            case(P_V)
               output(:,:,:,var) = (refs(GAMMA_REF) - 1.0_RP)*(Q(:,:,:,IRHOE) - 0.5_RP*&
                                 ( POW2(Q(:,:,:,IRHOU)) + POW2(Q(:,:,:,IRHOV)) + POW2(Q(:,:,:,IRHOW))) /Q(:,:,:,IRHO)) 
            case(T_V)
               output(:,:,:,var) = (refs(GAMMA_REF) - 1.0_RP) * refs(GAMMA_REF) * POW2(refs(MACH_REF)) * (Q(:,:,:,IRHOE) - 0.5_RP*&
                                 ( POW2(Q(:,:,:,IRHOU)) + POW2(Q(:,:,:,IRHOV)) + POW2(Q(:,:,:,IRHOW))) /POW2(Q(:,:,:,IRHO))) 

            case(MACH_V)

            case(S_V)
   
            case(Vabs_V)
               output(:,:,:,var) = sqrt(POW2(Q(:,:,:,IRHOU)) + POW2(Q(:,:,:,IRHOV)) + POW2(Q(:,:,:,IRHOW)))/Q(:,:,:,IRHO)

            case(Ht_V)

            case(RHOU_V)
               output(:,:,:,var) = Q(:,:,:,IRHOU) 

            case(RHOV_V)
               output(:,:,:,var) = Q(:,:,:,IRHOV) 
      
            case(RHOW_V)
               output(:,:,:,var) = Q(:,:,:,IRHOW) 

            case(RHOE_V)
               output(:,:,:,var) = Q(:,:,:,IRHOE) 

            end select

         end do

   
      end subroutine ComputeOutputVariables

      character(len=1024) function getOutputVariablesLabel()
         implicit none
!
!        ---------------
!        Local variables         
!        ---------------
!
         integer  :: iVar

         getOutputVariablesLabel = ""
         do iVar = 1, no_of_outputVariables
            write(getOutputVariablesLabel,'(A,A,A,A)') trim(getOutputVariablesLabel), ',"',trim(variableNames(outputVariableNames(iVar))),'"'
         end do

      end function getOutputVariablesLabel

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

      subroutine OutputVariablesForPreliminarVariable(iVar, output)
         implicit none
         integer, intent(in)     :: iVar
         integer, intent(out)    :: output(:)

         select case(iVar)

         case(Q_V)
            output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V/)

         case(Vvec_V)
            output = (/U_V, V_V, W_V/)

         case default
            output = iVar

         end select

      end subroutine OutputVariablesForPreliminarVariable

      integer function outputVariableForName(variableName)
         implicit none
         character(len=*),    intent(in)     :: variableName

         do outputVariableForName = 1, NO_OF_VARIABLES
            if ( trim(variableNames(outputVariableForName)) .eq. trim(variableName) ) return
         end do
!
!        Return "-1" if not found
!        ------------------------
         outputVariableForName = -1
      
      end function outputVariableForName

      integer function getNoOfCommas(input_line)
         implicit none
         character(len=*), intent(in)     :: input_line
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i

         getNoOfCommas = 0
         do i = 1, len_trim(input_line)
            if ( input_line(i:i) .eq. "," ) getNoOfCommas = getNoOfCommas + 1 
         end do

      end function getNoOfCommas

end module OutputVariables
