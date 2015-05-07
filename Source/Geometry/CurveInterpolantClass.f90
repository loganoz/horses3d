!
!////////////////////////////////////////////////////////////////////////
!
!      CurveInterpolantClass.f
!      Created: 2007-06-12 15:37:20 -0400 
!      By: David Kopriva
!
!      Contains: 
!         ALGORITHM 96: CurveInterpolant (Class)
!         ALGORITHM 97: CurveInterpolantProcedures
!            SUBROUTINE ConstructCurveInterpolant( this, N, nodes, values)
!            SUBROUTINE EvaluateAt( this, atLocation, givingResult )
!            SUBROUTINE DerivativeAt( this, atLocation, givingResult )
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE CurveInterpolantClass
      USE PolynomialInterpAndDerivsModule
      IMPLICIT NONE 
!
!---------------------------------------------------------------------
! This module defines a class that creates a curve from an interpolant
!---------------------------------------------------------------------
!
      TYPE CurveInterpolant
         INTEGER                    :: numberOfNodes
         REAL(KIND=RP), ALLOCATABLE :: nodes(:)
         REAL(KIND=RP), ALLOCATABLE :: values(:,:)
         REAL(KIND=RP), ALLOCATABLE :: bWeights(:)
      END TYPE CurveInterpolant
!
!     --------
!     Generics
!     --------
!
      INTERFACE Construct
         MODULE PROCEDURE ConstructCurveInterpolant
      END INTERFACE Construct
      INTERFACE Destruct
         MODULE PROCEDURE DestructCurveInterpolant
      END INTERFACE Destruct
!
!     ========
      CONTAINS 
!     ========
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructCurveInterpolant( this, N, nodes, values)
!
!-------------------------------------------------------------------
! Constructor pfor a curve interpolant
!-------------------------------------------------------------------
!
         TYPE(CurveInterpolant)         , INTENT(OUT) :: this
         INTEGER                        , INTENT(IN)  :: N
         REAL(KIND=RP), DIMENSION(0:N)  , INTENT(IN)  :: nodes
         REAL(KIND=RP), DIMENSION(0:N,2), INTENT(IN)  :: values
         
         ALLOCATE( this%nodes(0:N) )
         ALLOCATE( this%values(0:N,2) )
         ALLOCATE( this%bWeights(0:N) )
         
         this%numberOfNodes = N
         this%nodes         = nodes
         this%values        = values
         CALL BarycentricWeights( N, nodes, this%bWeights )
      
      END SUBROUTINE ConstructCurveInterpolant
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetValues( this, values )
!
!-------------------------------------------------------------------
! Constructor pfor a curve interpolant
!-------------------------------------------------------------------
!
         TYPE(CurveInterpolant)        , INTENT(INOUT) :: this
         REAL(KIND=RP), DIMENSION(0:,:), INTENT(IN)  :: values
         
         this%values        = values
      
      END SUBROUTINE SetValues
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructCurveInterpolant( this ) 
         TYPE(CurveInterpolant)       :: this
         IF(ALLOCATED(this%nodes))     DEALLOCATE( this%nodes )
         IF(ALLOCATED(this%values))    DEALLOCATE( this%values )
         IF(ALLOCATED(this%bWeights))  DEALLOCATE( this%bWeights )
      END SUBROUTINE DestructCurveInterpolant
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE EvaluateAt( this, atLocation, givingResult ) 
!
!     -------------------------------------------------------------------
!     Evaluate the interpolant of this at s returning v
!     -------------------------------------------------------------------
!
         TYPE(CurveInterpolant), INTENT(IN)  :: this
         REAL(KIND=RP)         , INTENT(IN)  :: atLocation
         REAL(KIND=RP)         , INTENT(OUT) :: givingResult(2)
         
         givingResult(1) = LagrangeInterpolation( atLocation, this%numberOfNodes, &
                           this%nodes, this%values(:,1), this%bWeights)
         givingResult(2) = LagrangeInterpolation( atLocation, this%numberOfNodes, &
                           this%nodes, this%values(:,2), this%bWeights)
   
      END SUBROUTINE EvaluateAt
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DerivativeAt( this, atLocation, givingResult ) 
!
!     -------------------------------------------------------------------
!     Evaluate the derivative of the interpolant of this at atLocation
!     returning givingResult.
!     -------------------------------------------------------------------
!
         TYPE(CurveInterpolant), INTENT(IN)  :: this
         REAL(KIND=RP)         , INTENT(IN)  :: atLocation
         REAL(KIND=RP)         , INTENT(OUT) :: givingResult(2)
         
         givingResult(1) = LagrangeInterpolantDerivative( atLocation, this%numberOfNodes, &
                           this%nodes, this%values(:,1), this%bWeights)
         givingResult(2) = LagrangeInterpolantDerivative( atLocation, this%numberOfNodes, &
                           this%nodes, this%values(:,2), this%bWeights)
   
      END SUBROUTINE DerivativeAt

      END MODULE CurveInterpolantClass
