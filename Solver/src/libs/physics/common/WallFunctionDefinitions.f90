!//////////////////////////////////////////////////////
!
!This module stores the definitions of the wall function and update the default values based on the controlVariables

#include "Includes.h"
#if defined(NAVIERSTOKES)
Module WallFunctionDefinitions  !

    use SMConstants
    Implicit None
!   
!  *****************************
!  Default everything to private
!  *****************************
!
    private
!
!  ****************
!  Public variables
!  ****************
!
    public STD_WALL, ABL_WALL
    public wallFuncIndex
    public kappa, WallC
    public y0, d
    public u_tau0, newtonAlpha, newtonTol, newtonMaxIter
    public useAverageV
!
!  ******************
!  Public definitions
!  ******************
!
    public Initialize_Wall_Function
!
    integer,       parameter                                  :: STD_WALL = 1
    integer,       parameter                                  :: ABL_WALL = 2

    integer                                                   :: wallFuncIndex

    real(kind=RP), parameter                                  :: DEFAULT_PLANE_DISPLACEMENT = 0.0_RP
    real(kind=RP), parameter                                  :: DEFAULT_VON_KARMAN = 0.38_RP, DEFAULT_WALL_C = 4.1_RP ! Constants are taken from Frere 2007
    real(kind=RP), parameter                                  :: DEFAULT_NEWTON_DAMP = 1.0_RP, DEFAULT_NEWTON_SEED = 1.0_RP, DEFAULT_NEWTON_TOL = 1.0E-10_RP
    integer,       parameter                                  :: DEFAULT_NEWTON_INTER = 100

    real(kind=RP)                                             :: y0, d                            ! for ABL wall func
    real(kind=RP)                                             :: u_tau0, newtonAlpha, newtonTol   ! for standrd wall func
    integer                                                   :: newtonMaxIter                    ! for standrd wall func
    real(kind=RP)                                             :: kappa, WallC                     ! for either

    logical                                                   :: useAverageV = .false.

    contains 
!   
!------------------------------------------------------------------------------------------------------------------------
!
    Subroutine Initialize_Wall_Function(controlVariables, correct)
        use FTValueDictionaryClass
        implicit none
        class(FTValueDictionary),  intent(in)  :: controlVariables
        logical,  intent(out)                  :: correct

!       ---------------
!       Local variables
!       ---------------
!
        character(len=STRING_CONSTANT_LENGTH)  :: wallFuncType
        real(kind=RP)                          :: Lref

        ! is false by default so it returns false if stopped
        correct = .false.

        if (.not. controlVariables % containsKey("wall function")) then
            return
        end if         

        wallFuncType = controlVariables % stringValueForKey("wall function", STRING_CONSTANT_LENGTH)
      
        select case (trim(wallFuncType))
        case ("standard")
            wallFuncIndex = STD_WALL
        case ("abl")
            wallFuncIndex = ABL_WALL
        case default
            write(STD_OUT,'(A,A,A)') 'Requested wall function "',trim(wallFuncType),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented functions are:"
            write(STD_OUT,'(A)') "  * Standard"
            write(STD_OUT,'(A)') "  * ABL"
            write(STD_OUT,'(A)') "Wall function will not be activated"
            return
        end select

        y0 = controlVariables % getValueOrDefault("wall roughness", 0.0_RP)
        if (wallFuncIndex .eq. ABL_WALL .and. y0 .eq. 0.0_RP) then
            write(STD_OUT,'(A)') "Wall roughness is mandatory for ABL wall function"
            write(STD_OUT,'(A)') "Wall function will not be activated"
            return
        end if 
        Lref = controlVariables % getValueOrDefault("reference length (m)", 1.0_RP)
        y0 = y0 / Lref

        d               = controlVariables % getValueOrDefault("wall plane displacement", DEFAULT_PLANE_DISPLACEMENT)
        u_tau0          = controlVariables % getValueOrDefault("wall function seed", DEFAULT_NEWTON_SEED)
        newtonAlpha     = controlVariables % getValueOrDefault("wall function damp", DEFAULT_NEWTON_DAMP)
        newtonTol       = controlVariables % getValueOrDefault("wall function tolerance", DEFAULT_NEWTON_TOL)
        newtonMaxIter   = controlVariables % getValueOrDefault("wall function max iter", DEFAULT_NEWTON_INTER)
        kappa           = controlVariables % getValueOrDefault("wall function kappa", DEFAULT_VON_KARMAN)
        WallC           = controlVariables % getValueOrDefault("wall function c", DEFAULT_WALL_C)

        d = d / Lref

        useAverageV = controlVariables%logicalValueForKey("wall function use average")

        !todo: see if there are negative values and return if that's the case
            ! write(STD_OUT,'(A)') "Wall function will not be activated"

        ! if it arrives here everything has gone well
        correct = .true.

    End Subroutine Initialize_Wall_Function
!   
!------------------------------------------------------------------------------------------------------------------------
!
End Module WallFunctionDefinitions
#endif