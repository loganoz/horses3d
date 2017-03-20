!
!////////////////////////////////////////////////////////////////////////
!
!      Implicit_NJ.f90
!      Created: 2017-03-17 15:21:00 +0100 
!      By:  Carlos Redondo (module for 2D) 
!           Andrés Rueda   (3D implementation and changes) 
!      Module for computing colorings
!           In the moment it is optimized for computing colorings of Legendre-Gauss 
!           discretizations for the Euler or Navier-Stokes equations.
!           Implement Legendre-Gauss-Lobatto colorings!  
!
!////////////////////////////////////////////////////////////////////////
MODULE ColorsClass
   
   USE DGSEMClass

   IMPLICIT NONE

   TYPE Colors
         INTEGER                             :: ncolors ! number of colors
         INTEGER,DIMENSION(:), ALLOCATABLE   :: elmnts  ! color ordered elements
         INTEGER,DIMENSION(:), ALLOCATABLE   :: bounds  ! idx of the first element on a color
         INTEGER                             :: ntotal  ! total  numer of elements
      CONTAINS
         PROCEDURE                           :: construct   => constructcolors !construct(nbr)
         PROCEDURE                           :: info        => getcolorsinfo   ! prints coloring info  
   END TYPE Colors
   
   PRIVATE
   PUBLIC Colors
   PUBLIC constructcolors, getcolorsinfo  

   CONTAINS
      SUBROUTINE constructcolors(this, nbr,flowIsNavierStokes)
         CLASS(Colors), INTENT(OUT)          :: this
         TYPE(neighbour),DIMENSION(:)        :: nbr
         LOGICAL                             :: flowIsNavierStokes
         
         INTEGER                             :: ncolored = 0
         LOGICAL, DIMENSION(:), ALLOCATABLE  :: colored, used
         LOGICAL                             :: allcolored = .FALSE.
         INTEGER                             :: i, j, counter, idx
         INTEGER                             :: ntotal, ncolors, maxcolor
         INTEGER, DIMENSION(:), ALLOCATABLE  :: colors
         
         ntotal = SIZE(nbr)
         this%ntotal = ntotal
         ALLOCATE(used(0:ntotal)) !0 correspond to boundary "neighbour"
         ALLOCATE(colored(ntotal))
         ALLOCATE(colors(ntotal))
         colored(:) = .FALSE.
         maxcolor = 0
         
         IF (flowIsNavierStokes) THEN
            ! Color elements using 2 neighbours (works for BR1 using LG)
            DO WHILE (.NOT. allcolored)
               used(:) = .FALSE.
               maxcolor = maxcolor + 1
               elemloop: DO i = 1, ntotal
                  IF (.NOT. colored(i)) THEN
                     used(0) = .FALSE. !boundary neigbour is empty
                     
                     DO j=1, 7   !arueda: hardcoded 7 for conforming hexahedral meshes
                        IF (nbr(i)%elmnt(j) /= 0) THEN
                           IF (ANY(used(nbr(nbr(i)%elmnt(j))%elmnt))) CYCLE elemloop
                        END IF
                     END DO
                     
                     !used(nbr(i)%elmnt)= .TRUE.
                     
                     DO j=1, 7   !arueda: hardcoded 7 for conforming hexahedral meshes
                        IF (nbr(i)%elmnt(j) /= 0) used(nbr(nbr(i)%elmnt(j))%elmnt) = .TRUE.
                     END DO
                     
                     colored(i) = .TRUE.
                     colors(i) = maxcolor
                     ncolored = ncolored + 1
                     IF (ncolored == ntotal) allcolored = .TRUE.
                     
                  ENDIF
               ENDDO elemloop
            ENDDO
         
         ELSE
         
            !Color elements using  1 neighbour (works for Euler and NS -IP using LG)
            DO WHILE (.NOT. allcolored)
               used(:) = .FALSE.
               maxcolor = maxcolor + 1
               DO i = 1, ntotal
                  IF (.NOT. colored(i)) THEN
                     used(0) = .FALSE. !boundary neigbour is empty
                     IF (.NOT. ANY(used(nbr(i)%elmnt))) THEN
                        used(nbr(i)%elmnt)= .TRUE.
                        colored(i) = .TRUE.
                        colors(i) = maxcolor
                        ncolored = ncolored + 1
                        IF (ncolored == ntotal) allcolored = .TRUE.
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         
         END IF
         
         this%ncolors = maxcolor
         ALLOCATE(this%elmnts(ntotal))
         ALLOCATE(this%bounds(this%ncolors + 1))
         
         ! order elements acording to colors
         idx = 1
         DO i = 1, this%ncolors
            this%bounds(i)= idx
            counter = 0
            DO j = 1, ntotal
               IF (colors(j) == i) THEN
                  this%elmnts(idx) = j
                  counter = counter + 1
                  idx = idx + 1
               ENDIF      
            ENDDO
         ENDDO
         this%bounds(i)= idx
      END SUBROUTINE constructcolors
      
      SUBROUTINE getcolorsinfo(this)
         CLASS(Colors)        :: this
         INTEGER              :: i
       
         WRITE(*,'(A13,I2)') "# of colors: ", this%ncolors
         WRITE(*,*) "Element list:"         
         WRITE(*,'(*(I5,1X))') (this%elmnts(i), i = 1,this%ntotal)
         WRITE(*,*) "New color indexes"
         WRITE(*,'(*(I5,1X))') (this%bounds(i), i= 1,this%ncolors + 1)
      END SUBROUTINE getcolorsinfo
      
  
END MODULE ColorsClass
