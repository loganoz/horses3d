!
! ////////////////////////////////////////////////////////////////////
! HORSES3D Converter
!     Main program horsesConverter
!	This module read control file "horses2Converter.convert"
!
! ////////////////////////////////////////////////////////////////////
!
MODULE getTaskConverter
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
!  Read Input Task "horses2Converter.convert"
!  ------------------------------------------
    SUBROUTINE readGetTaskInput(filename, meshFile1, boundaryFile1, resultFile1, meshFile2, boundaryFile2, polyOrder2, taskType, VTKfile, Ref)
        IMPLICIT NONE
        CHARACTER(LEN=LINE_LENGTH), INTENT(IN)      :: filename
        CHARACTER(LEN=LINE_LENGTH), INTENT(OUT)     :: meshFile1, boundaryFile1, resultFile1, VTKfile
		CHARACTER(LEN=LINE_LENGTH), INTENT(OUT)     :: meshFile2, boundaryFile2
		INTEGER                   , INTENT(OUT)     :: polyOrder2(3)
		INTEGER					  , INTENT(OUT)     :: taskType
		REAL(KIND=RP)             , INTENT(OUT)     :: Ref(4)

!
!       ---------------
!       Local variables
!       ---------------
!
        INTEGER                     :: fID, io, ist, cStart, cEnd, cDot
        CHARACTER(LEN=LINE_LENGTH)  :: inputLine
		logical :: fileExist
!
!  		Default Value
!  		-------------	
		Ref=0
		meshFile1=""
		boundaryFile1=""
		resultFile1=""
		meshFile2=""
		boundaryFile2=""
		VTKfile=""
		taskType=0
		polyOrder2=0
!
!  		Check and open control file
!  		---------------------------
		inquire(file=filename, exist=fileExist)
		if (fileExist) then
			OPEN(newunit=fID,file=filename,status="old",action="read", iostat=io)
			DO
				READ(fid,'(A132)',IOSTAT = ist) inputLine
				! Exit at the end of File
				IF(ist /= 0 ) EXIT

				cStart = INDEX(inputLine,'=')
				cEnd   = INDEX(inputLine, ' ', .true. )

				IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Mesh Filename 1") THEN
					meshFile1=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Boundary Filename 1") THEN
					boundaryFile1=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Result 1") THEN
					resultFile1=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Mesh Filename 2") THEN
					meshFile2=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Boundary Filename 2") THEN
					boundaryFile2=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Polynomial Order") THEN
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
					cStart = INDEX(inputLine,'(')
					cEnd   = INDEX(inputLine, ')', .true. )
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
					cStart = INDEX(inputLine,',')
					cEnd   = LEN(inputLine)
					read (inputLine(1: cStart-1),'(I10)') polyOrder2(1) ! Convert to Integer
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd )))
					cStart = INDEX(inputLine,',')
					cEnd   = LEN(inputLine)
					read (inputLine(1: cStart-1),'(I10)') polyOrder2(2)
					inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd )))
					read (inputLine,'(I10)') polyOrder2(3)
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Task") THEN
					IF (TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))=="meshInterpolation") THEN
						taskType=1
					ELSE IF (TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))=="horsesMesh2OF") THEN
						taskType=2
					ELSE IF (TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))=="OF2Horses") THEN
						taskType=3
					END IF
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="VTK file") THEN
					VTKfile=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Reynolds Number") THEN
					cDot=0
					cDot = INDEX(inputLine(cStart+1: cEnd-1),'.', .true. )
					if (cDot .eq. 0) then
						read (inputLine(cStart+1: cEnd-1),'(ES8.0)')  Ref(1)
					else 
						read (inputLine(cStart+1: cEnd-1),'(ES10.2)')  Ref(1) 
					end if 
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Mach Number") THEN
					cDot=0
					cDot = INDEX(inputLine(cStart+1: cEnd-1),'.', .true. )
					if (cDot .eq. 0) then
						read (inputLine(cStart+1: cEnd-1),'(ES8.0)')  Ref(2)
					else 
						read (inputLine(cStart+1: cEnd-1),'(ES10.2)')  Ref(2) 
					end if  
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Reference pressure (Pa)") THEN
					cDot=0
					cDot = INDEX(inputLine(cStart+1: cEnd-1),'.', .true. )
					if (cDot .eq. 0) then
						read (inputLine(cStart+1: cEnd-1),'(ES8.0)')  Ref(3)
					else 
						read (inputLine(cStart+1: cEnd-1),'(ES10.2)')  Ref(3) 
					end if 
				ELSE IF (TRIM(ADJUSTL(inputLine(1:cStart-1)))=="Reference temperature (K)") THEN
					cDot=0
					cDot = INDEX(inputLine(cStart+1: cEnd-1),'.', .true. )
					if (cDot .eq. 0) then
						read (inputLine(cStart+1: cEnd-1),'(ES8.0)')  Ref(4)
					else 
						read (inputLine(cStart+1: cEnd-1),'(ES10.2)')  Ref(4) 
					end if 
				END IF
				
			END DO
			CLOSE(fID) ! Close File
		else
!
!  			Create template for control input file if does not exist
!  			--------------------------------------------------------
			write(STD_OUT,'(10X,A,A)') "ERROR::horses2foam.convert does not exist. File is generated"
			OPEN(newunit=fID,file=filename,status="new",action="write")
			WRITE(fid,'(A)')'!//////////////////////////////////////////////////////////////////////////////////////////////////////////!'
			WRITE(fid,'(A)')'!////////////////////////////////// HORSES3D to FOAM format for PARAVIEW //////////////////////////////////!'
			WRITE(fid,'(A)')'!////////////////////////////////// ---------Control Input File--------- //////////////////////////////////!'
			WRITE(fid,'(A)')'!//////////////////////////////////////////////////////////////////////////////////////////////////////////!'
			WRITE(fid,'(A)')'Task                = meshInterpolation/horsesMesh2OF/OF2Horses'
			WRITE(fid,'(A)')'Mesh Filename 1     = ../MESH/'
			WRITE(fid,'(A)')'Boundary Filename 1 = ../MESH/'
			WRITE(fid,'(A)')'Result 1            = Result.hsol'
			WRITE(fid,'(A)')'Mesh Filename 2     = ../MESH_2/'
			WRITE(fid,'(A)')'Boundary Filename 2 = ../MESH_2/'
			WRITE(fid,'(A)')'Polynomial Order    = ( 1 , 1 , 1 )'
			WRITE(fid,'(A)')'--------------- For OF2Horses ---------------' 
			WRITE(fid,'(A)')'VTK file                   = foamToVTK.vtk'
			WRITE(fid,'(A)')'Reynolds Number            = '
			WRITE(fid,'(A)')'Mach Number                = ' 
			WRITE(fid,'(A)')'Reference pressure (Pa)    = 101325'
			WRITE(fid,'(A)')'Reference temperature (K)  = 288.889'
			CLOSE(fID)
		    call EXIT(0)
		end if 

    END SUBROUTINE readGetTaskInput
!
!////////////////////////////////////////////////////////////////////////
!    

END MODULE getTaskConverter
