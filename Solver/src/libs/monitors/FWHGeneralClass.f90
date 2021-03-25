!
!   @File:    ObserverClass.f90
!   @Author:  Oscar Marino (oscar.marino@upm.es)
!   @Created: Mar 25 2020
!   @Last revision date: 
!   @Last revision author: 
!   @Last revision commit: 
!
!//////////////////////////////////////////////////////
!
!This class represents the general behaviour of the Fwoc Williams and Hawckings aero accoustic analogy

#include "Includes.h"
Module FWHGeneralClass  !

    use SMConstants
    Implicit None

!
!   *****************************
!   Main FWH class definition
!   *****************************
    type FWHClass

    end type FWHClass
           ! se debe construir desde la clase general de FW, esta debe hacer algo similar a la de monitores: crear update, escribir,
           ! crear archivo de escritura, allocar, leer de control file, etc...
           ! debe crear la clase zona que tenga todas las caras de source (aglomerar todas las zone de BC dadas)
           ! debe tener como atributos t, iter, zones, buffer, dt_update, y similar a monitor


End Module FWHGeneralClass
