!
! ////////////////////////////////////////////////////////////////////
!   HORSES3D to Foam Result - readHorses2Foam Module
!
!      This module read horses2foam.convert input file
!
!////////////////////////////////////////////////////////////////////////////////////////

MODULE readHorses2Foam
      USE SMConstants
      IMPLICIT NONE
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
!  Read File Input "horses2foam.convert"
!  ------------------------
    SUBROUTINE readHorses2FoamConvert(filename, meshFile, boundaryFile, resultFile, outputVariables, polyOrder, generateMesh)
        IMPLICIT NONE
        CHARACTER(LEN=LINE_LENGTH), INTENT(IN)      :: filename
        CHARACTER(LEN=LINE_LENGTH), INTENT(OUT)     :: meshFile, boundaryFile, resultFile, outputVariables
        INTEGER                   , INTENT(OUT)     :: polyOrder(3)
		LOGICAL                   , INTENT(OUT)     :: generateMesh

        !///////////////////////////
        ! LOCAL VARIABLES
        !///////////////////////////
        INTEGER                     :: fID, io, ist, cStart, cEnd
        CHARACTER(LEN=LINE_LENGTH)  :: inputLine
		logical :: fileExist
		
		generateMesh =.true.

        ! OPEN FILE

		inquire(file=filename, exist=fileExist)
		if (fileExist) then
			OPEN(newunit=fID,file=filename,status="old",action="read", iostat=io)
			DO
				READ(fid,'(A132)',IOSTAT = ist) inputLine
				! Exit at the end of File
				IF(ist /= 0 ) EXIT

				cStart = INDEX(inputLine,'=')
				cEnd   = INDEX(inputLine, ' ', .true. )

				IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Mesh Filename") THEN
					meshFile=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Boundary Filename") THEN
					boundaryFile=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Result") THEN
					resultFile=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Polynomial Order") THEN
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
					cStart = INDEX(inputLine,'(')
					cEnd   = INDEX(inputLine, ')', .true. )
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
					cStart = INDEX(inputLine,',')
					cEnd   = LEN(inputLine)
					read (inputLine(1: cStart-1),'(I10)') polyOrder(1) ! Convert to Integer
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd )))
					cStart = INDEX(inputLine,',')
					cEnd   = LEN(inputLine)
					read (inputLine(1: cStart-1),'(I10)') polyOrder(2)
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd )))
					read (inputLine,'(I10)') polyOrder(3)
				 ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Output Variables ") THEN
					outputVariables=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				 ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Create Foam Mesh ") THEN
					READ( inputLine( cStart+1:cEnd-1 ), * ) generateMesh
				END IF
			END DO
			CLOSE(fID) ! Close File
		else
			write(STD_OUT,'(10X,A,A)') "ERROR::horses2foam.convert does not exist. File is generated"
			OPEN(newunit=fID,file=filename,status="new",action="write")
			WRITE(fid,'(A)')'!//////////////////////////////////////////////////////////////////////////////////////////////////////////!'
			WRITE(fid,'(A)')'!////////////////////////////////// HORSES3D to FOAM format for PARAVIEW //////////////////////////////////!'
			WRITE(fid,'(A)')'!////////////////////////////////// ---------Control Input File--------- //////////////////////////////////!'
			WRITE(fid,'(A)')'!//////////////////////////////////////////////////////////////////////////////////////////////////////////!'
			WRITE(fid,'(A)')'Mesh Filename     = ../MESH/'
			WRITE(fid,'(A)')'Boundary Filename = ../MESH/'
			WRITE(fid,'(A)')'Result            = Result.hsol'
			WRITE(fid,'(A)')'Polynomial Order  = ( 1 , 1 , 1 )'
			WRITE(fid,'(A)')'Output Variables  = V, Qcrit, p, pt, T, rho, Mach, gradV, omega'
			WRITE(fid,'(A)')'Create Foam Mesh  = .true.'
			CLOSE(fID)
		    call EXIT(0)
		end if 

    END SUBROUTINE readHorses2FoamConvert

END MODULE readHorses2Foam
