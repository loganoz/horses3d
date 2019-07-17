!
!//////////////////////////////////////////////////////
!
!   @File:    JacobianDefinitions.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Wed Jul 17 11:53:06 2019
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
!  Module containing types and definitions for the Jacobian matrices
!
!//////////////////////////////////////////////////////
!
module JacobianDefinitions
   use SMConstants
   implicit none
!
!  **********
!  Parameters
!  **********
!
   real(kind=RP), parameter :: JACEPS = 1.e-9_RP ! Minimum value of a Jacobian entry (smaller values are considered as 0._RP)
   
end module JacobianDefinitions
   
