!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!   HORSES3D to Foam Result - foamMeshStorage Module
!
!      This module convert Horses3D mesh storage to foamMesh Storage
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE foamResultVTKStorageConverter
    USE SMConstants
    IMPLICIT NONE
	
	
	TYPE VTK1D_t
		REAL(KIND=RP),ALLOCATABLE   :: data(:)
	END TYPE VTK1D_t
		
	TYPE VTK2D_t
		REAL(KIND=RP),ALLOCATABLE   :: data(:,:)
	END TYPE VTK2D_t
	
    !
    TYPE VTKResult_t
		CHARACTER(len=LINE_LENGTH)  :: resultName
        INTEGER                     :: nPoints
		TYPE(VTK2D_t)           	:: x
		TYPE(VTK2D_t)				:: U
		TYPE(VTK1D_t)           	:: p
        TYPE(VTK1D_t)           	:: rho
		REAL(KIND=RP),ALLOCATABLE   :: Q(:,:)			!Q(1:5, nPoints) 
		REAL(KIND=RP)               :: Re				!Reynolds Number
		REAL(KIND=RP)               :: Mach				!Mach Number
        REAL(KIND=RP)               :: rhoRef
		REAL(KIND=RP)               :: VRef
        REAL(KIND=RP)               :: pRef
        REAL(KIND=RP)               :: TRef
		contains
		PROCEDURE   :: Construct     => constructVTKResult
    END TYPE VTKResult_t
    !

!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!     This subroutine construct VTKResult from VTKfile 
!
      SUBROUTINE constructVTKResult(self, VTKfile, Ref)
		IMPLICIT NONE
		CLASS(VTKResult_t)               :: self
		CHARACTER(len=*), INTENT(in)     :: VTKfile
		REAL(KIND=RP)   , INTENT(in)     :: Ref(4)
		
!
!        	---------------
!        	Local variables
!        	---------------
!
		integer :: fid, flag, ist, cStart, cEnd, cDot
		integer	:: i, j, k, ii, jj, nData
		CHARACTER(LEN=LINE_LENGTH)  :: inputLine, input
		REAL(KIND=RP)               :: R, gamma, pRefInv, a, vrefInv, rhoInv
		
		R=287.15_RP
		gamma = 1.4_RP 
		
		
		
		self % resultName = trim(VTKfile)
		self % Re		  = Ref(1)
		self % Mach       = Ref(2)
		self % pRef       = Ref(3)
		self % TRef       = Ref(4)
		
		write(STD_OUT,'(10X,A,A)') "Loading VTK Result File:"
		write(STD_OUT,'(10X,A,A)') "---------------------"
		write(STD_OUT,'(30X,A,A30,A)') "->","VTK File: ", trim(VTKfile)	
		
!
!       Open file
!       ---------
		fid = putVTKFileInReadDataMode(VTKfile)

!
!       Get number of points
!       --------------------
		j=0
		 DO WHILE (j.eq.0)
			READ(fid,'(A132)',IOSTAT = ist) inputLine
			IF (TRIM(ADJUSTL(inputLine(1:6)))=="POINTS") THEN
				cStart = INDEX(inputLine,'S')
				cEnd   = INDEX(inputLine, 'f', .true. )
				inputLine=TRIM(ADJUSTL(inputLine( cStart+1: cEnd-1 )))
				read (inputLine,'(I10)') self % nPoints ! Convert to Integer
				j=1
			END IF
		 END DO 
		
		 write(STD_OUT,'(30X,A,A30,I10)') "->","Number of Points: ", self % nPoints	

!
!        Read Points Location
!        --------------------
		 CALL read2DVTKData(self % x, fid, self % nPoints, 3)		
		 write(STD_OUT,'(30X,A,A30)') "->","Points data loaded"	
!
!        Read Data 
!        ---------			 
		 j=0
		 DO WHILE (j.eq.0)
			READ(fid,'(A132)',IOSTAT = ist) inputLine
			cEnd   = INDEX(inputLine,' ', .false. )
			inputLine=inputLine(1:cEnd-1)
			IF ((TRIM(ADJUSTL(inputLine))=="POINT_DATA")) THEN
				j=1
			end if 
		 END DO 
!
!        Obtain p, rho, and U from VTK Files Points data
!        -----------------------------------------------					 		 
		 DO i=1,3       
			j=0
			 DO WHILE (j.eq.0)
				READ(fid,'(A132)',IOSTAT = ist) inputLine
				cEnd   = INDEX(inputLine,' ', .false. )
				inputLine=inputLine(1:cEnd-1)
				IF (TRIM(ADJUSTL(inputLine))=="rho") THEN
					CALL read1DVTKData(self % rho, fid, self % nPoints)	
					write(STD_OUT,'(30X,A,A30)') "->","rho data loaded"	
					j=1
				ELSE IF (TRIM(ADJUSTL(inputLine))=="p") THEN
					CALL read1DVTKData(self % p, fid, self % nPoints)	
					write(STD_OUT,'(30X,A,A30)') "->","p data loaded"	
					j=1
				ELSE IF (TRIM(ADJUSTL(inputLine))=="U") THEN
					CALL read2DVTKData(self % U, fid, self % nPoints, 3)	
					write(STD_OUT,'(30X,A,A30)') "->","U data loaded"	
					j=1
				END IF 
				if (ist.gt.0) then
					write(STD_OUT,'(30X,A,A30)') "->","ERROR - Check Input File !!!"
				    call EXIT(0)
				else if (ist.lt.0) then 
					write(STD_OUT,'(30X,A,A30)') "->","ERROR - end of file reached. rho, p, or U might not in the files !!!"	
					call EXIT(0)
				end if 
			 END DO 
		  END DO 
!
!        Close file
!        ----------
		 close(fid)
		 
		 ! Normalized with reference values
		   self % rhoRef = self % pRef / (R * self % TRef)
		   rhoInv = 1_RP/self % rhoRef
		   a=sqrt(gamma* self%pRef /self % rhoRef)
		   self % VRef = self % Mach *a
		   vrefInv = 1_RP/self % VRef
		   pRefInv = 1_RP/(((self % Mach)**2) * gamma * self % rhoRef * R * self % TRef)
		   
		   self % rho % data = self % rho % data * rhoInv
		   self % U % data = self % U % data * vrefInv
		   self % p % data = self % p % data * pRefInv
		   
		   allocate(self % Q(1:5, self % nPoints))
		   
		   do i=1, self % nPoints
				self % Q(1,i) = self % rho % data(i)
				self % Q(2,i) = self % rho % data(i) * self % U % data(1,i) 
				self % Q(3,i) = self % rho % data(i) * self % U % data(2,i) 
				self % Q(4,i) = self % rho % data(i) * self % U % data(3,i) 
				self % Q(5,i) = self % p % data (i) / (gamma-1.0_RP) + 0.5_RP*(self % Q(2,i) * self % U % data(1,i) + self % Q(3,i) * self % U % data(2,i) & 
									+ self % Q(4,i) * self % U % data(3,i) )
		   end do 
		
		write(STD_OUT,'(10X,A,A)') "Finish Loading VTK Results"

      END SUBROUTINE constructVTKResult
!
!////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!     This routine open VTK file and put it in the read mode
!     ------------------------------------------------------
      FUNCTION putVTKFileInReadDataMode(VTKfile) RESULT (fid)
		IMPLICIT NONE
		CHARACTER(len=*), INTENT(in)     :: VTKfile
		INTEGER                          :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(VTKfile), status="old", action="read", &  
                           form="formatted" , access="stream")
        END FUNCTION putVTKFileInReadDataMode
!
!     --------------------------------------
!     This routine read 1D data from VTK
!     --------------------------------------
      SUBROUTINE read1DVTKData(self, fid, nPoints)
		IMPLICIT NONE
		CLASS(VTK1D_t) 	   , INTENT(INOUT)		:: self
		INTEGER            , INTENT(IN)         :: fid
		INTEGER		       , INTENT(IN)         :: nPoints
!
!        	---------------
!        	Local variables
!        	---------------
!
		integer :: ist, cDot, j, k, ii, jj, nData, cEnd
		CHARACTER(LEN=LINE_LENGTH)  :: inputLine, input
		
		 ALLOCATE(self % data (nPoints))
			
		 j=1
		 ii=0
		 jj=1
		DO WHILE (j.le. nPoints)
		READ(fid,'(A132)',IOSTAT = ist) inputLine
		nData=10
		if (j.eq.1) nData=11
		DO k=1, nData
			inputLine=trim(inputLine)
			ii=ii+1
			if (ii.gt.nPoints) exit
			cEnd   = INDEX(inputLine,' ', .false. )
			input  = inputLine(1:cEnd-1)
			inputLine = inputLine(cEnd+1:len(inputLine))
			cDot = INDEX(input,'.', .true. )
			if (cDot .eq. 0) then
				read (input,'(ES7.0)') self % data(ii) 
			else 
				read (input,'(ES14.7)') self % data(ii) 
			end if 
		END DO 
		j=j+nData
		END DO
		END SUBROUTINE read1DVTKData
!
!     --------------------------------------
!     This routine read 2D data from VTK
!     --------------------------------------
      SUBROUTINE read2DVTKData(self, fid, nPoints, nsize)
		IMPLICIT NONE
		CLASS(VTK2D_t) 	   , INTENT(INOUT)		:: self
		INTEGER            , INTENT(IN)         :: fid
		INTEGER		       , INTENT(IN)         :: nPoints
		INTEGER		       , INTENT(IN)         :: nsize
!
!        	---------------
!        	Local variables
!        	---------------
!
		 integer :: ist, cDot, j, k, ii, jj, nData, cEnd
		 CHARACTER(LEN=LINE_LENGTH)  :: inputLine, input
		
		 ALLOCATE(self % data (nsize, nPoints))
			
		 j=1
		 ii=0
		 jj=1
		 DO WHILE (j.le.nPoints*nsize)
			READ(fid,'(A132)',IOSTAT = ist) inputLine
			nData=10
			if (j.eq.1) nData=11
			DO k=1, nData
				inputLine=trim(inputLine)
				ii=ii+1
				if (ii.eq.nsize+1) then
					ii=1
					jj=jj+1
				end if
				if (jj.gt.nPoints) exit
				cEnd   = INDEX(inputLine,' ', .false. )
				input  = inputLine(1:cEnd-1)
				inputLine = inputLine(cEnd+1:len(inputLine))
				cDot = INDEX(input,'.', .true. )
				if (cDot .eq. 0) then
					read (input,'(ES7.0)') self % data(ii,jj) ! Convert to Integer
				else 
					read (input,'(ES14.7)') self % data(ii,jj) ! Convert to Integer
				end if 
			END DO 
			j=j+nData
		 END DO 
		END SUBROUTINE read2DVTKData
!
END MODULE foamResultVTKStorageConverter
!
!////////////////////////////////////////////// END OF FILE //////////////////////////////////////////////
!
