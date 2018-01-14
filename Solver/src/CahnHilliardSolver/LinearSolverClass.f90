!
!//////////////////////////////////////////////////////
!
!   @File:    LinearSolverClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:40 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!////////////////////////////////////////////////////////////////////////
!
!      LinearSolverClass.f90
!      Created: 2017-04-12 00:20:00 +0100 
!      By: Andrés Rueda
!
!      Module to load all linear solvers
!
!////////////////////////////////////////////////////////////////////////
MODULE LinearSolverClass
   USE PetscSolverClass
   USE MKLPardisoSolverClass
   USE IterativeSolverClass
   USE MultigridSolverClass
   IMPLICIT NONE
END MODULE LinearSolverClass
