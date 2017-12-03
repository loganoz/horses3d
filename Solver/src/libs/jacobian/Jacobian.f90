!
!////////////////////////////////////////////////////////////////////////
!
!      Jacobian.f90
!      Created: 2011-09-27 16:21:15 +0200 
!      By: Gonzalo Rubio Calzado
!          Carlos Redondo
!          Andrés Rueda  
!
!////////////////////////////////////////////////////////////////////////
!
Module Jacobian 

USE SMConstants
USE HexMeshClass
USE PhysicsStorage 
IMPLICIT NONE 

   private
   public Look_for_neighbour, ijk2local, local2ijk
!
!-------------------------------------------------------------------
! This module contains subroutines to calculate and print the 
! Jacobian matrix using sparse notation
!-------------------------------------------------------------------
!

!========
 CONTAINS
!========
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Look_for_neighbour(this, mesh)  !arueda: ready for 3D
      IMPLICIT NONE 
!
!     -----------------------------------------------------------------------------
!     This subroutine finds the neighbours of all elements (conforming meshes only)
!        current implementation only finds neighbour 3D elements, but not faces 
!        (modify commented part to perform this)
!     -----------------------------------------------------------------------------
!
!
!     ----------------
!     Input parameters
!     ----------------
!
      TYPE(Neighbour)                   :: this(:)
      TYPE(HexMesh)                     :: mesh
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                         :: i,j,iEl  ! k, iEdge, iSlave
      
      DO iEl = 1, SIZE(mesh%elements)
!~       k = 1      
!         PRINT*, "Element = ", iEl
         this(iEl)%elmnt(7) = iEl  ! The last one is itself
         DO j = 1, 6
            IF (mesh%elements(iEl)%NumberOfConnections(j) == 0) THEN
               this(iEl)%elmnt(j) = 0
            ELSE
               this(iEl)%elmnt(j) = mesh%elements(iEl)%Connection(j)%ElementIDs(1)
            ENDIF
         ENDDO
!~         DO iEdge = 1, SIZE(sem%mortar)
!~            IF (mesh%masters(iEdge)%elementID == iEl) THEN
!~               this(iEl)%edge(k) = iEdge
!~               k = k + 1
!~!            PRINT*, "Element = ", iEl, "conection with edge", iEdge             
!~            ELSE
!~               DO iSlave = 1, SIZE(sem%mesh%masters(iEdge)%sElementIDs)
!~                  IF (sem%mesh%masters(iEdge)%sElementIDs(iSlave) == iEl) THEN 
!~!            PRINT*, "Element = ", iEl, "conection with edge", iEdge                   
!~                     this(iEl)%edge(k) = iEdge
!~                     k = k + 1
!~                  ENDIF
!~               ENDDO
!~            ENDIF
!~         ENDDO 
         !print*, "iEl", iEl
         !print*, "this(iEl)%elmnt", this(iEl)%elmnt
         !print*, "this(iEl)%edge", this(iEl)%edge
         !read(*,*)
      ENDDO
    
    END SUBROUTINE 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
!  Returns the local index relative to an element from the local coordinates: i(lagrange node x), j(lagrange node y), 
!  k(lagrange node z), l(equation number)
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
   
   FUNCTION ijk2local(i,j,k,l,N_EQN,Nx,Ny,Nz) RESULT(idx)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)   :: i, j, k, l, Nx, Ny, Nz, N_EQN
      INTEGER               :: idx
      
      IF (l < 1 .OR. l > N_EQN)  STOP 'error in ijk2local, l has a wrong value'
      IF (i < 0 .OR. i > Nx)     STOP 'error in ijk2local, i has a wrong value'
      IF (j < 0 .OR. j > Ny)     STOP 'error in ijk2local, j has a wrong value'
      IF (k < 0 .OR. k > Nz)     STOP 'error in ijk2local, k has a wrong value'
      
      idx = k*(Nx+1)*(Ny+1)*N_EQN + j*(Nx+1)*N_EQN + i*N_EQN + l
   END FUNCTION
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Returns the coordinates relative to an element: l(equation number), i(lagrange node x), j(lagrange node y), k(lagrange node z)
!  from the local index  
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
   
   FUNCTION local2ijk(idx,N_EQN,Nx,Ny,Nz) RESULT (indices)
   
      INTEGER, INTENT(IN)   :: idx, Nx, Ny, Nz, N_EQN
      INTEGER               :: indices(4)
      INTEGER               :: tmp1, tmp2
      
      IF (idx < 1 .OR. idx > (Nx+1)*(Ny+1)*(Nz+1)*N_EQN) STOP 'error in local2ijk, idx has wrong value'
      
      indices(4) = (idx-1) / ((Nx+1)*(Ny+1) * N_EQN)
      tmp1       = MOD((idx-1),((Nx+1)*(Ny+1) * N_EQN) )
      indices(3) = tmp1 / ((Nx+1)*N_EQN)
      tmp2       = MOD(tmp1,((Nx+1)*N_EQN))
      indices(2) = tmp2 / (N_EQN)
      indices(1) = MOD(tmp2, N_EQN) + 1
   END FUNCTION




!
!////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE JacobianCalculation(sem, t, ExternalState, ExternalGradient, jacobianFileName) 
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Calculation of the Jacobian matrix of the equation from the steady
!~ !     state in SD. The precission is specified using the variable eps
!~ !     -------------------------------------------------------------------
!~ !
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: t
!~       EXTERNAL                          :: ExternalState, ExternalGradient
!~       CHARACTER(LEN=LINE_LENGTH)        :: jacobianFileName    
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER                         :: i,j,k, elmnt, N, M, counter, nvpn
!~       REAL(KIND=RP),DIMENSION(N_eqn)   :: eps                 
      
!~       INTEGER                         :: outputUnit = 22
!~       CHARACTER(LEN=LINE_LENGTH)      :: outputFile     
!~       INTEGER                         :: nUnit = 23
!~       CHARACTER(LEN=LINE_LENGTH)      :: nFile         
!~       INTEGER                         :: iUnit = 24
!~       CHARACTER(LEN=LINE_LENGTH)      :: iFile 
!~       INTEGER                         :: tUnit = 25
!~       CHARACTER(LEN=LINE_LENGTH)      :: tFile 

!~        outputFile = TRIM(jacobianFileName)//TRIM("Adense.dat")
!~        nFile = TRIM(jacobianFileName)//TRIM("NVPN.dat")
!~        iFile = TRIM(jacobianFileName)//TRIM("NNVI.dat")
!~        tFile = TRIM(jacobianFileName)//TRIM("CONT.dat")
!~ !
!~ !     ---------
!~ !     OPEN UNIT
!~ !     ---------
!~ !            
!~       OPEN(UNIT=outputUnit,FILE=outputFile)
!~       OPEN(UNIT=nUnit,FILE=nFile)
!~       OPEN(UNIT=iUnit,FILE=iFile)
!~       OPEN(UNIT=tUnit,FILE=tFile)
      
!~       counter = 0
!~       nvpn    = 1
!~       WRITE(nUnit, *) nvpn

!~       CALL ComputeTimeDerivative( sem )
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%QDot = sem%dgS(elmnt)%QDot
!~       END DO                                 
!~       DO k = 1,N_EQN
!~          eps(k) = JacEpsilon(sem,k)
!~       END DO 

!~       PRINT*, "PlotJacDense has been activated"
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          N       = sem%mesh%elements(elmnt)%Nx
!~          M       = sem%mesh%elements(elmnt)%Ny
!~          PRINT*, "Element", elmnt, "of", SIZE(sem%mesh%elements)
!~          DO j = 0,M
!~             DO i = 0,N
!~                 DO k = 1,N_EQN
!~                    !counter = k + j*4 + i*4*(N+1) + (elmnt-1)*4*(N+1)**2
!~                    counter = counter + 1
!~                    CALL Jac_ith(sem,eps(k),elmnt,i,j,k &
!~                                ,ExternalState, ExternalGradient)
!~                    CALL PlotJacDense ( sem, outputUnit )
!~                    !CALL PlotJacSparse( sem, nvpn, counter, nUnit, iUnit, tUnit )
!~                 END DO 
!~             END DO
!~          END DO 
!~       END DO
!~ !
!~ !     ----------
!~ !     CLOSE UNIT
!~ !     ----------
!~ !
!~       CLOSE(outputUnit)
!~       CLOSE(nUnit)
!~       CLOSE(iUnit)
!~       CLOSE(tUnit)
        
    
!~     END SUBROUTINE 
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE JacobianCalculationAllocate(sem, t, ExternalState, ExternalGradient, Jac) 
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Calculation of the Jacobian matrix of the equation from the steady
!~ !     state in SD. The precission is specified using the variable eps
!~ !     -------------------------------------------------------------------
!~ !
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: t
!~       EXTERNAL                          :: ExternalState, ExternalGradient
      
!~       REAL(KIND=RP)     :: Jac(:,:)    
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER                          :: i,j,k, elmnt, N, M, counter
!~       REAL(KIND=RP),DIMENSION(N_eqn)   :: eps                 
      
      
!~       counter = 0

!~       CALL ComputeTimeDerivative( sem )
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%QDot = sem%dgS(elmnt)%QDot
!~       END DO                                 
!~       DO k = 1,N_EQN
!~          eps(k) = JacEpsilon(sem,k)
!~       END DO 
!~       !!eps = 1.0_RP
      
!~       !!!I have to check this for the Jacobian!!!
!~       PRINT*, "Calculation of the Jacobian. Storage in Jac(:,:)"
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          N       = sem%mesh%elements(elmnt)%Nx
!~          M       = sem%mesh%elements(elmnt)%Ny
!~          PRINT*, "Element", elmnt, "of", SIZE(sem%mesh%elements)
!~          DO j = 0,M
!~             DO i = 0,N
!~                 DO k = 1,N_EQN
!~                    counter = counter + 1                                       
!~                    CALL Jac_ith(sem,eps(k),elmnt,i,j,k &
!~                                ,ExternalState, ExternalGradient)
!~                    CALL WriteAithInJac(sem, Jac(:,counter)) 
!~                 END DO 
!~             END DO
!~          END DO 
!~       END DO    
!~       PRINT*, "Cest finit"
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%QDot = sem%JAC(elmnt)%QDot
!~       END DO   
!~     END SUBROUTINE JacobianCalculationAllocate     
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE LocalJacobianCalculationAllocate(sem, t, ExternalState, ExternalGradient, Jac) 
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Calculation of the Jacobian matrix of the equation from the steady
!~ !     state in SD. The precission is specified using the variable eps
!~ !     -------------------------------------------------------------------
!~ !
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: t
!~       EXTERNAL                          :: ExternalState, ExternalGradient
      
!~       REAL(KIND=RP)     :: Jac(:,:)    
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER                          :: i,j,k, elmnt, N, M, counter
!~       REAL(KIND=RP),DIMENSION(N_eqn)   :: eps                 
      
      
!~       counter = 0

!~       CALL ComputeTimeDerivativeLTE( sem, 0.0_RP, &
!~                                       ExternalState, &
!~                                       ExternalGradient )
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%QDot = sem%dgS(elmnt)%QDot
!~       END DO                                 
!~       DO k = 1,N_EQN
!~          eps(k) = JacEpsilon(sem,k)
!~       END DO 
!~       !!eps = 1.0_RP
      
!~       !!!I have to check this for the Jacobian!!!
!~       PRINT*, "Calculation of the Jacobian. Storage in Jac(:,:)"
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          N       = sem%mesh%elements(elmnt)%Nx
!~          M       = sem%mesh%elements(elmnt)%Ny
!~          PRINT*, "Element", elmnt, "of", SIZE(sem%mesh%elements)
!~          DO j = 0,M
!~             DO i = 0,N
!~                 DO k = 1,N_EQN
!~                    counter = counter + 1                                       
!~                    CALL Local_Jac_ith(sem,eps(k),elmnt,i,j,k &
!~                                ,ExternalState, ExternalGradient)
!~                    CALL WriteAithInJac(sem, Jac(:,counter)) 
!~                 END DO 
!~             END DO
!~          END DO 
!~       END DO    
!~       PRINT*, "Cest finit"
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%QDot = sem%JAC(elmnt)%QDot
!~       END DO   
!~     END SUBROUTINE LocalJacobianCalculationAllocate         
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE ResidualCalculationAllocate(sem, t, ExternalState, ExternalGradient, res) 
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Residual calculation. Use in the TE estimation with NC solution.
!~ !     -------------------------------------------------------------------
!~ !
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: t
!~       EXTERNAL                          :: ExternalState, ExternalGradient
      
!~       REAL(KIND=RP)                     :: res(:)    
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER                         :: i,j,k, elmnt, N, M, counter
!~       REAL(KIND=RP),DIMENSION(N_eqn)   :: eps                 
!~ !
!~ !     ---------
!~ !     OPEN UNIT
!~ !     ---------
!~ !      
      
!~       counter = 0

!~       CALL ComputeTimeDerivative( sem )
!~       DO elmnt = 1, SIZE(sem%mesh%elements)  
!~          sem%JAC(elmnt)%Aith = sem%dgS(elmnt)%QDot
!~       END DO                                 

!~       CALL WriteAithInJac(sem, res)          

!~     END SUBROUTINE ResidualCalculationAllocate    
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE TEvectorAllocate(sem, t, ExternalState, ExternalGradient, TEvector) 
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Residual calculation. Use in the TE estimation with NC solution.
!~ !     -------------------------------------------------------------------
!~ !
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: t
!~       EXTERNAL                          :: ExternalState, ExternalGradient
      
!~       REAL(KIND=RP)                     :: TEvector(:)    
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER                          :: i,j,k, elmnt, N, M, counter
!~       REAL(KIND=RP),DIMENSION(N_eqn)   :: eps                 
!~ !
!~ !     ---------
!~ !     OPEN UNIT
!~ !     ---------
!~ !      
      
!~       counter = 0

!~       DO elmnt = 1, SIZE(sem%mesh%elements)  
!~          sem%JAC(elmnt)%Aith = sem%dgS(elmnt)%QDot
!~       END DO                                 

!~       CALL WriteAithInJac(sem, TEvector)          

!~     END SUBROUTINE TEvectorAllocate            
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE Jac_ith( sem, eps, Pelmnt, iP, jP, kP, &
!~                         ExternalState, ExternalGradient )
!~ !
!~ !     ------------------------------------------
!~ !     Calculation of the ith row of the Jacobian
!~ !     ------------------------------------------
!~ !
!~       IMPLICIT NONE
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: eps
!~       EXTERNAL                          :: ExternalState, ExternalGradient   
!~       INTEGER                           :: Pelmnt, iP, jP, kP
!~ !
!~ !     ---------------
!~ !     Local variables
!~ !     ---------------
!~ !

!~       INTEGER                           :: elmnt, i, j, k

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%eps    = 0.0_RP
!~       END DO 
!~       !IF (JAC%u(i).ge.0.0_RP) JAC%eps(i) =  eps
!~       !IF (JAC%u(i) <  0.0_RP) JAC%eps(i) = -eps        
!~       sem%JAC(Pelmnt)%eps(kP,iP,jP) = eps

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%Q    = sem%dgS(elmnt)%Q
!~       END DO 
      
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%Q = sem%dgS(elmnt)%Q + sem%JAC(elmnt)%eps
!~       END DO 

!~       CALL ComputeTimeDerivative( sem )
                                      
!~      ! DO elmnt = 1, SIZE(sem%mesh%elements)
!~      !    sem%JAC(elmnt)%UPQDot = sem%dgS(elmnt)%QDot
!~      ! END DO                                       

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%Aith = (sem%dgS(elmnt)%QDot-sem%JAC(elmnt)%QDot)/eps!-sem%JAC(elmnt)%Aith
!~       END DO  

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%Q = sem%JAC(elmnt)%Q
!~       END DO  

      
!~    END SUBROUTINE Jac_ith
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE Local_Jac_ith( sem, eps, Pelmnt, iP, jP, kP, &
!~                         ExternalState, ExternalGradient )
!~ !     ------------------------------------------
!~ !     Calculation of the ith row of the Jacobian
!~ !     ------------------------------------------
!~       IMPLICIT NONE
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~       TYPE(DGSem)                       :: sem
!~       REAL(KIND=RP)                     :: eps
!~       EXTERNAL                          :: ExternalState, ExternalGradient   
!~       INTEGER                           :: Pelmnt, iP, jP, kP
!~ !     ---------------
!~ !     Local variables
!~ !     ---------------
!~       INTEGER                           :: elmnt, i, j, k

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%eps    = 0.0_RP
!~       END DO 
!~       !IF (JAC%u(i).ge.0.0_RP) JAC%eps(i) =  eps
!~       !IF (JAC%u(i) <  0.0_RP) JAC%eps(i) = -eps        
!~       sem%JAC(Pelmnt)%eps(kP,iP,jP) = eps

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%Q    = sem%dgS(elmnt)%Q
!~       END DO 
      
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%Q = sem%dgS(elmnt)%Q + sem%JAC(elmnt)%eps
!~       END DO 

!~       CALL ComputeTimeDerivativeLTE( sem, 0.0_RP, &
!~                                       ExternalState, &
!~                                       ExternalGradient )
                                      
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%UPQDot = sem%dgS(elmnt)%QDot
!~       END DO                                       

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%Aith = (sem%JAC(elmnt)%UPQDot-sem%JAC(elmnt)%QDot)/eps
!~       END DO  

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%Q = sem%JAC(elmnt)%Q
!~       END DO  

      
!~    END SUBROUTINE Local_Jac_ith
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !

!~    FUNCTION JacEpsilon(sem,k) RESULT(eps)
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Calculation of the perturbation to calculate the Jacobian matrix 
!~ !     -------------------------------------------------------------------
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)    ,INTENT(IN)  :: sem
!~       REAL(KIND=RP)               :: eps
!~       INTEGER                     :: k
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER         :: i,j,N, M, elmnt, NTot
!~       REAL(KIND=RP)   :: b
      
!~       b = 1e-6_RP       
!~       eps = 0.0_RP
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~       N       = sem%mesh%elements(elmnt)%Nx
!~       M       = sem%mesh%elements(elmnt)%Ny
!~          DO j = 0,M
!~             DO i = 0,N
!~                eps = eps + b*sem%dgS(elmnt)%Q(k,i,j)
!~             END DO
!~          END DO 
!~       END DO
!~       NTot = N_eqn*SIZE(sem%mesh%elements)*(N+1)*(M+1)
!~       eps = eps + b
!~       eps = eps/(NTot)       
    
!~    END FUNCTION  
   
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE PlotJacDense( sem, outputUnit )
!~ !
!~ !     ------------------------------------------
!~ !     Calculation of the ith row of the Jacobian
!~ !     ------------------------------------------
!~ !
!~       IMPLICIT NONE
!~ !
!~ !     -----------------
!~ !     Input parameters:
!~ !     -----------------
!~ !
!~       TYPE(DGSem)     :: sem
!~       INTEGER         :: outputUnit
!~ !
!~ !     ---------------
!~ !     Local variables
!~ !     ---------------

!~       INTEGER              :: i, j, k, elmnt
!~       INTEGER              :: N, M, NTot
     
       
!~       DO elmnt = 1, SIZE(sem%mesh%elements) 
!~           N = sem%mesh%elements(elmnt)%Nx
!~           M = sem%mesh%elements(elmnt)%Ny
!~           DO j = 0,M
!~              DO i = 0,N
!~                 DO k = 1,N_EQN                
!~                   WRITE(outputUnit, *) sem%JAC(elmnt)%Aith(k,i,j)
!~ !                  IF (ABS(sem%JAC(elmnt)%Aith(k,i,j)) > 1.d-4) THEN 
!~ !                     PRINT*, elmnt, k, j, i, sem%JAC(elmnt)%Aith(k,i,j)
!~ !                  ENDIF
!~                 END DO
!~              END DO
!~           END DO
!~       END DO 
      
      
!~    END SUBROUTINE PlotJacDense   
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE WriteAithInJac( sem, Jac )
!~ !
!~ !     ------------------------------------------
!~ !     Calculation of the ith row of the Jacobian
!~ !     ------------------------------------------
!~ !
!~       IMPLICIT NONE
!~ !
!~ !     -----------------
!~ !     Input parameters:
!~ !     -----------------
!~ !
!~       TYPE(DGSem)     :: sem
!~       REAL(KIND = RP) :: Jac(:)
!~ !
!~ !     ---------------
!~ !     Local variables
!~ !     ---------------

!~       INTEGER              :: i, j, k, elmnt
!~       INTEGER              :: N, M, NTot
!~       INTEGER              :: counter
!~       !PRINT*, "i   j   Aith(k,i,j)"
!~       counter = 0 
!~       DO elmnt = 1, SIZE(sem%mesh%elements) 
!~           N = sem%mesh%elements(elmnt)%Nx
!~           M = sem%mesh%elements(elmnt)%Ny
!~           DO j = 0,M
!~              DO i = 0,N
!~                 DO k = 1,N_EQN
!~                   counter = counter + 1                 
!~                   Jac(counter) = sem%JAC(elmnt)%Aith(k,i,j)   
!~                   !PRINT*, i, j,  sem%JAC(elmnt)%Aith(k,i,j)             
!~                 END DO
!~              END DO
!~           END DO
!~       END DO 
!~       !PRINT*, "____________________"
      
      
!~    END SUBROUTINE WriteAithInJac     
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE PlotJacSparse( sem, nvpn, last, nUnit, iUnit, tUnit )
!~ !
!~ !     ------------------------------------------
!~ !     Calculation of the ith row of the Jacobian
!~ !     ------------------------------------------
!~ !
!~       IMPLICIT NONE
!~ !
!~ !     -----------------
!~ !     Input parameters:
!~ !     -----------------
!~ !
!~       TYPE(DGSem)     :: sem
!~       INTEGER         :: last, nUnit, iUnit, tUnit 
!~       INTEGER         :: nvpn
!~ !
!~ !     ---------------
!~ !     Local variables
!~ !     ---------------
!~ !

      

!~       INTEGER              :: i, j, k, elmnt
!~       INTEGER              :: N, M, NTot, counter

      
!~       !NTot = 4*SIZE(sem%mesh%elements)*(N+1)*(N+1) 

!~       counter = 0
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~           N = sem%mesh%elements(elmnt)%Nx
!~           M = sem%mesh%elements(elmnt)%Ny
!~           DO j = 0,M 
!~              DO i = 0,N 
!~                 DO k = 1,4 
!~                    counter = counter + 1
!~                    IF ( ABS(sem%JAC(elmnt)%Aith(k,i,j)) > 1.0e-8_RP  ) THEN
!~                       nvpn = nvpn + 1
!~                       WRITE(iUnit, *) counter !k + j*4 + i*(4*(N+1)) + (elmnt-1)*4*(N+1)**2
!~                       WRITE(tUnit, *) sem%JAC(elmnt)%Aith(k,i,j)
!~                    END IF
!~                 END DO
!~              END DO
!~           END DO
      
!~       END DO 
      
!~       WRITE(nUnit, *) nvpn
      
!~    END SUBROUTINE PlotJacSparse   
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~       SUBROUTINE ExportMesh( sem ) 
!~          IMPLICIT NONE
!~ !
!~ !        ---------
!~ !        Arguments
!~ !        ---------
!~ !
!~          TYPE(DGSem)   :: sem
!~ !
!~ !        ---------------
!~ !        Local variables
!~ !        ---------------
!~ !
!~          INTEGER                                  :: i, j, k, N, M
         
!~          OPEN(UNIT=11,FILE="Mesh.dat")
!~          DO k = 1, SIZE(sem%mesh%elements) 
!~             N = sem%mesh%elements(k)%Nx
!~             M = sem%mesh%elements(k)%Ny
!~             !WRITE(11,*) "ZONE I=", N+1, ",J=",N+1,", F=POINT"
!~             DO j= 0, M 
!~                DO i = 0, N                  
!~                   WRITE(11,*) sem%mesh%elements(k)%geom%x(i,j), &
!~                                               sem%mesh%elements(k)%geom%y(i,j)
!~                END DO
!~             END DO
!~          END DO
!~          CLOSE(11)
         
!~       END SUBROUTINE ExportMesh   
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!
   
!~ !
!~     SUBROUTINE JacobianCalculationNeighbour(sem, semOptJac, nbr, t, ExternalState, ExternalGradient, jacobianFileName) 
!~       IMPLICIT NONE 
!~ !
!~ !     -------------------------------------------------------------------
!~ !     Calculation of the Jacobian matrix of the equation from the steady
!~ !     state in SD. The precission is specified using the variable eps.
!~ !     This version is optimized for the Euler case so only the neighbours
!~ !     are calculated
!~ !     -------------------------------------------------------------------
!~ !
!~ !
!~ !     ----------------
!~ !     Input parameters
!~ !     ----------------
!~ !
!~       TYPE(DGSem)                       :: sem, semOptJac
!~       TYPE(Neighbour)                   :: nbr(:)
!~       REAL(KIND=RP)                     :: t
!~       EXTERNAL                          :: ExternalState, ExternalGradient
!~       CHARACTER(LEN=LINE_LENGTH)        :: jacobianFileName    
!~ !
!~ !     ---------------
!~ !     Local Variables
!~ !     ---------------
!~ !
!~       INTEGER                         :: i,j,k, elmnt, N, M, counter, nvpn
!~       REAL(KIND=RP),DIMENSION(N_eqn)   :: eps                 
      
!~       INTEGER                         :: outputUnit = 22
!~       CHARACTER(LEN=LINE_LENGTH)      :: outputFile     
!~       INTEGER                         :: nUnit = 23
!~       CHARACTER(LEN=LINE_LENGTH)      :: nFile         
!~       INTEGER                         :: iUnit = 24
!~       CHARACTER(LEN=LINE_LENGTH)      :: iFile 
!~       INTEGER                         :: tUnit = 25
!~       CHARACTER(LEN=LINE_LENGTH)      :: tFile 

!~        outputFile = TRIM(jacobianFileName)//TRIM("Adense.dat")
!~        nFile = TRIM(jacobianFileName)//TRIM("NVPN.dat")
!~        iFile = TRIM(jacobianFileName)//TRIM("NNVI.dat")
!~        tFile = TRIM(jacobianFileName)//TRIM("CONT.dat")


!~ !
!~ !     ---------
!~ !     OPEN UNIT
!~ !     ---------
!~ !            
!~       OPEN(UNIT=outputUnit,FILE=outputFile)
!~       OPEN(UNIT=nUnit,FILE=nFile)
!~       OPEN(UNIT=iUnit,FILE=iFile)
!~       OPEN(UNIT=tUnit,FILE=tFile)
      
!~       counter = 0
!~       nvpn    = 1
!~       WRITE(nUnit, *) nvpn

!~       CALL ComputeTimeDerivative( sem )
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%QDot = sem%dgS(elmnt)%QDot
!~       END DO                                 
      
!~          semOptJac = sem
      
!~       DO k = 1,N_EQN
!~          eps(k) = JacEpsilon(sem,k)
!~       END DO 
!~       PRINT*, "PlotJacDense has been activated"
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          N       = sem%mesh%elements(elmnt)%Nx
!~          M       = sem%mesh%elements(elmnt)%Ny
!~          PRINT*, "Element", elmnt, "of", SIZE(sem%mesh%elements)
!~          DO j = 0,M
!~             DO i = 0,N
!~                 DO k = 1,N_EQN
!~                    !counter = k + j*4 + i*4*(N+1) + (elmnt-1)*4*(N+1)**2
!~                    counter = counter + 1              
!~                    CALL Jac_ith_Neigh(sem, semOptJac, nbr, eps(k),elmnt,i,j,k &
!~                                ,ExternalState, ExternalGradient)
!~                    !CALL Jac_ith(sem, eps(k),elmnt,i,j,k &
!~                    !            ,ExternalState, ExternalGradient)                                                       
!~                    CALL PlotJacDense ( sem, outputUnit )
!~                    !CALL PlotJacSparse( sem, nvpn, counter, nUnit, iUnit, tUnit )
!~                 END DO 
!~             END DO
!~          END DO 
!~       END DO
!~ !
!~ !     ----------
!~ !     CLOSE UNIT
!~ !     ----------
!~ !
!~       CLOSE(outputUnit)
!~       CLOSE(nUnit)
!~       CLOSE(iUnit)
!~       CLOSE(tUnit)
        
    
!~     END SUBROUTINE 
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
!~     SUBROUTINE Jac_ith_Neigh( sem, semOptJac, nbr, eps, Pelmnt, iP, jP, kP, &
!~                         ExternalState, ExternalGradient )
!~ !     ------------------------------------------
!~ !     Calculation of the ith row of the Jacobian
!~ !     Optimized version for the Euler case.
!~ !     ------------------------------------------
!~       IMPLICIT NONE
!~       TYPE(DGSem)                       :: sem, semOptJac
!~       TYPE(Neighbour)                   :: nbr(:)
!~       REAL(KIND=RP)                     :: eps
!~       EXTERNAL                          :: ExternalState, ExternalGradient   
!~       INTEGER                           :: Pelmnt, iP, jP, kP

!~       INTEGER                           :: elmnt, i, j, k
!~       LOGICAL, SAVE                     :: isfirst = .TRUE.

!~       IF (ISFIRST) THEN
!~          isfirst = .FALSE.
!~          DO elmnt = 1, SIZE(sem%mesh%elements)
!~             sem%JAC(elmnt)%eps    = 0.0_RP
!~             semOptJac%JAC(elmnt)%eps    = 0.0_RP
!~          END DO 
!~       ENDIF
      
!~       sem%JAC(Pelmnt)%eps(kP,iP,jP) = eps

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%Q    = sem%dgS(elmnt)%Q
!~       END DO 
      
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%QDot =  sem%JAC(elmnt)%QDot
!~          sem%dgS(elmnt)%bound =  semOptJac%dgS(elmnt)%bound 
!~       ENDDO
!~       !sem%dgS = semOptJac%dgS
      
!~       sem%dgS(Pelmnt)%Q = sem%dgS(Pelmnt)%Q + sem%JAC(Pelmnt)%eps

         

!~       CALL ComputeTimeDerivativeJac1( sem, nbr, pElmnt, 0.0_RP, ExternalState,ExternalGradient )

!~       !DO elmnt = 1, SIZE(sem%mesh%elements)
!~       !   sem%JAC(elmnt)%UPQDot = sem%dgS(elmnt)%QDot
!~       !END DO  
                                     
!~       !CALL ComputeTimeDerivative( sem, 0.0_RP,ExternalState,ExternalGradient )

!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%JAC(elmnt)%Aith = (sem%dgS(elmnt)%QDot-sem%JAC(elmnt)%QDot)/eps !- sem%JAC(elmnt)%Aith
!~       END DO  

!~ !     Restaurar todo para el cálculo de la siguiente columna
      
!~       DO elmnt = 1, SIZE(sem%mesh%elements)
!~          sem%dgS(elmnt)%Q = sem%JAC(elmnt)%Q
!~          sem%dgs(elmnt)%QDot = sem%JAC(elmnt)%QDot
!~       END DO  
!~       !CALL ComputeTimeDerivative( sem, 0.0_RP,ExternalState,ExternalGradient )
!~       sem%JAC(Pelmnt)%eps(kP,iP,jP) = 0.0_RP
      
!~    END SUBROUTINE Jac_ith_Neigh
!~ !
!~ !////////////////////////////////////////////////////////////////////////
!~ !
END MODULE  


