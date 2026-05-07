#include "Includes.h"
module RigidMotionFrameClass
   use SMConstants
   use Utilities, only: toLower
   use ParamfileRegions, only: getSquashedLine
   use MPI_Process_Info, only: MPI_Process
   implicit none

   private
   public :: GetRigidMotionFramePose
   public :: GetRigidMotionFrameVelocity

   integer, parameter :: MAX_RIGID_MOTION_FRAMES = 16

   type RigidMotionFrame_t
      logical                         :: configured = .false.
      logical                         :: loaded     = .false.
      character(len=LINE_LENGTH)      :: name       = ""
      character(len=LINE_LENGTH)      :: tableFile  = ""
      logical                         :: anglesInDegrees = .true.
      real(kind=RP)                   :: referenceCenter(NDIM) = 0.0_RP
      real(kind=RP)                   :: referenceRotation(NDIM,NDIM) = 0.0_RP
      integer                         :: n = 0
      real(kind=RP), allocatable      :: time(:)
      real(kind=RP), allocatable      :: center(:,:)
      real(kind=RP), allocatable      :: angles(:,:)
   end type RigidMotionFrame_t

   type(RigidMotionFrame_t), save :: frames(MAX_RIGID_MOTION_FRAMES)

contains

   subroutine GetRigidMotionFramePose(frameName, t, center, R)
      implicit none

      character(len=*), intent(in)  :: frameName
      real(kind=RP),    intent(in)  :: t
      real(kind=RP),    intent(out) :: center(NDIM)
      real(kind=RP),    intent(out) :: R(NDIM,NDIM)

      integer       :: frameID
      real(kind=RP) :: angles(NDIM)
      real(kind=RP) :: tableCenter(NDIM), tableRotation(NDIM,NDIM)

      call EnsureFrameReady(frameName, frameID)
      call PoseAtTime(frames(frameID), t, tableCenter, angles)
      call EulerToRotationMatrix(angles, tableRotation)

      center = frames(frameID)%referenceCenter + matmul(frames(frameID)%referenceRotation, tableCenter)
      R      = matmul(frames(frameID)%referenceRotation, tableRotation)
   end subroutine GetRigidMotionFramePose

   subroutine GetRigidMotionFrameVelocity(frameName, t, centerDot, omega)
      implicit none

      character(len=*), intent(in)  :: frameName
      real(kind=RP),    intent(in)  :: t
      real(kind=RP),    intent(out) :: centerDot(NDIM), omega(NDIM)

      integer       :: frameID
      real(kind=RP) :: tableCenterDot(NDIM), tableOmega(NDIM)

      call EnsureFrameReady(frameName, frameID)
      call VelocityAtTime(frames(frameID), t, tableCenterDot, tableOmega)

      centerDot = matmul(frames(frameID)%referenceRotation, tableCenterDot)
      omega     = matmul(frames(frameID)%referenceRotation, tableOmega)
   end subroutine GetRigidMotionFrameVelocity

   subroutine EnsureFrameReady(frameName, frameID)
      implicit none

      character(len=*), intent(in)  :: frameName
      integer,          intent(out) :: frameID

      frameID = GetFrameID(frameName)

      if (.not. frames(frameID)%configured) call ConfigureFrame(frames(frameID), frameName)
      if (.not. frames(frameID)%loaded)     call LoadFrameTable(frames(frameID))
   end subroutine EnsureFrameReady

   integer function GetFrameID(frameName) result(frameID)
      implicit none

      character(len=*), intent(in) :: frameName

      integer                    :: i, emptySlot
      character(len=LINE_LENGTH) :: requestedName, storedName

      emptySlot = 0
      requestedName = frameName
      call toLower(requestedName)

      do i = 1, MAX_RIGID_MOTION_FRAMES
         if (frames(i)%configured) then
            storedName = frames(i)%name
            call toLower(storedName)
            if (trim(storedName) == trim(requestedName)) then
               frameID = i
               return
            end if
         else if (emptySlot == 0) then
            emptySlot = i
         end if
      end do

      if (emptySlot == 0) then
         call MotionFrameError("ERROR: Too many rigid motion frames configured.")
      end if

      frameID = emptySlot
   end function GetFrameID

   subroutine ConfigureFrame(frame, frameName)
      implicit none

      type(RigidMotionFrame_t), intent(inout) :: frame
      character(len=*),         intent(in)    :: frameName

      character(len=LINE_LENGTH) :: caseFile
      character(len=1024)        :: line, lineNoComment, lowerLine
      character(len=1024)        :: key, value, lowerKey, lowerValue
      character(len=1024)        :: inLabel, lowerInLabel
      real(kind=RP)              :: xAxis(NDIM), yAxis(NDIM), zAxis(NDIM)
      integer                    :: fid, io, pos
      logical                    :: inside, found, hasXAxis, hasYAxis, hasZAxis

      call get_command_argument(1, caseFile)
      if (len_trim(caseFile) == 0) then
         call MotionFrameError("ERROR: Cannot configure rigid motion frame without a control file argument.")
      end if

      write(inLabel,'(A,A)') "#define rigid motion frame ", trim(frameName)
      lowerInLabel = inLabel
      call toLower(lowerInLabel)

      inside = .false.
      found  = .false.
      hasXAxis = .false.
      hasYAxis = .false.
      hasZAxis = .false.
      call SetIdentity(frame%referenceRotation)

      open(newunit=fid, file=trim(caseFile), status="old", action="read", iostat=io)
      if (io /= 0) call MotionFrameError("ERROR: Could not open control file for rigid motion frame setup.")

      do
         read(fid,'(A)',iostat=io) line
         if (io /= 0) exit

         lineNoComment = StripInlineComment(line)
         lowerLine = lineNoComment
         call toLower(lowerLine)

         if (getSquashedLine(lowerLine) == getSquashedLine(lowerInLabel)) then
            inside = .true.
            found  = .true.
            cycle
         end if

         if (inside .and. getSquashedLine(lowerLine) == getSquashedLine("#end")) exit
         if (.not. inside) cycle

         pos = index(lineNoComment, "=")
         if (pos == 0) pos = index(lineNoComment, ":")
         if (pos == 0) cycle

         key = adjustl(lineNoComment(1:pos-1))
         value = adjustl(lineNoComment(pos+1:))

         lowerKey = key
         lowerValue = value
         call toLower(lowerKey)
         call toLower(lowerValue)

         select case (trim(getSquashedLine(lowerKey)))
         case ("type")
            if (index(trim(lowerValue), "prescribed") == 0 .or. index(trim(lowerValue), "table") == 0) then
               call MotionFrameError("ERROR: Rigid motion frame type must be prescribed table.")
            end if

         case ("tablefile", "motiontable", "file")
            frame%tableFile = TrimQuotes(value)

         case ("referencecenter", "referenceorigin", "origin", "center")
            call ReadVector3(value, "reference center", frame%referenceCenter)

         case ("referencexaxis", "xaxis", "axisx")
            call ReadVector3(value, "reference x axis", xAxis)
            hasXAxis = .true.

         case ("referenceyaxis", "yaxis", "axisy")
            call ReadVector3(value, "reference y axis", yAxis)
            hasYAxis = .true.

         case ("referencezaxis", "zaxis", "axisz")
            call ReadVector3(value, "reference z axis", zAxis)
            hasZAxis = .true.

         case ("angleunits", "angleunit")
            select case (trim(getSquashedLine(lowerValue)))
            case ("degrees", "degree", "deg")
               frame%anglesInDegrees = .true.
            case ("radians", "radian", "rad")
               frame%anglesInDegrees = .false.
            case default
               call MotionFrameError("ERROR: Rigid motion frame angle units must be degrees or radians.")
            end select

         case ("rotationconvention")
            if (trim(getSquashedLine(lowerValue)) /= "zyx") then
               call MotionFrameError("ERROR: Only ZYX rigid motion frame rotation convention is implemented.")
            end if
         end select
      end do

      close(fid)

      if (.not. found) then
         call MotionFrameError("ERROR: Rigid motion frame '"//trim(frameName)//"' was not found in the control file.")
      end if

      if (len_trim(frame%tableFile) == 0) then
         call MotionFrameError("ERROR: Rigid motion frame '"//trim(frameName)//"' needs a Table file entry.")
      end if

      call BuildReferenceRotation(frame%referenceRotation, hasXAxis, xAxis, hasYAxis, yAxis, hasZAxis, zAxis)

      frame%name = trim(frameName)
      frame%configured = .true.
   end subroutine ConfigureFrame

   subroutine LoadFrameTable(frame)
      implicit none

      type(RigidMotionFrame_t), intent(inout) :: frame

      character(len=1024) :: line
      integer             :: fid, io, nRows, row, i, j
      real(kind=RP)       :: values(7)
      logical             :: ok
      real(kind=RP)       :: angleOffset, currentAngle, deltaAngle

      open(newunit=fid, file=trim(frame%tableFile), status="old", action="read", iostat=io)
      if (io /= 0) then
         call MotionFrameError("ERROR: Could not open rigid motion table: "//trim(frame%tableFile))
      end if

      nRows = 0
      do
         read(fid,'(A)',iostat=io) line
         if (io /= 0) exit
         call ParseMotionLine(line, values, ok)
         if (ok) nRows = nRows + 1
      end do
      close(fid)

      if (nRows < 2) call MotionFrameError("ERROR: A rigid motion table needs at least two numeric rows.")

      allocate(frame%time(nRows))
      allocate(frame%center(NDIM,nRows))
      allocate(frame%angles(NDIM,nRows))

      open(newunit=fid, file=trim(frame%tableFile), status="old", action="read", iostat=io)
      if (io /= 0) then
         call MotionFrameError("ERROR: Could not reopen rigid motion table: "//trim(frame%tableFile))
      end if

      row = 0
      do
         read(fid,'(A)',iostat=io) line
         if (io /= 0) exit

         call ParseMotionLine(line, values, ok)
         if (.not. ok) cycle

         row = row + 1
         frame%time(row)     = values(1)
         frame%center(:,row) = values(2:4)
         frame%angles(:,row) = values(5:7)
         if (frame%anglesInDegrees) frame%angles(:,row) = frame%angles(:,row) * DEG2RAD
      end do
      close(fid)

      frame%n = nRows
      if (row /= frame%n) call MotionFrameError("ERROR: Rigid motion table changed while reading.")

      do i = 2, frame%n
         if (frame%time(i) <= frame%time(i-1)) then
            call MotionFrameError("ERROR: Rigid motion table times must be strictly increasing.")
         end if
      end do

      do j = 1, NDIM
         angleOffset = 0.0_RP
         currentAngle = frame%angles(j,1)
         do i = 2, frame%n
            currentAngle = frame%angles(j,i) + angleOffset
            deltaAngle = currentAngle - frame%angles(j,i-1)

            do while (deltaAngle > PI)
               angleOffset = angleOffset - 2.0_RP * PI
               currentAngle = frame%angles(j,i) + angleOffset
               deltaAngle = currentAngle - frame%angles(j,i-1)
            end do

            do while (deltaAngle < -PI)
               angleOffset = angleOffset + 2.0_RP * PI
               currentAngle = frame%angles(j,i) + angleOffset
               deltaAngle = currentAngle - frame%angles(j,i-1)
            end do

            frame%angles(j,i) = currentAngle
         end do
      end do

      frame%loaded = .true.
      if (MPI_Process%isRoot) then
         write(STD_OUT,'(30X,A,A,A,A,I0,A)') "->", "Rigid motion frame ", trim(frame%name), &
                                             ": loaded ", frame%n, " poses."
      end if
   end subroutine LoadFrameTable

   subroutine ParseMotionLine(line, values, ok)
      implicit none

      character(len=*), intent(in)  :: line
      real(kind=RP),    intent(out) :: values(7)
      logical,          intent(out) :: ok

      character(len=1024) :: cleanLine
      integer             :: io

      ok = .false.
      cleanLine = adjustl(StripInlineComment(line))
      if (len_trim(cleanLine) == 0) return
      if (cleanLine(1:1) == "#") return

      read(cleanLine, *, iostat=io) values
      ok = (io == 0)
   end subroutine ParseMotionLine

   subroutine PoseAtTime(frame, tq, center, angles)
      implicit none

      type(RigidMotionFrame_t), intent(in)  :: frame
      real(kind=RP),            intent(in)  :: tq
      real(kind=RP),            intent(out) :: center(NDIM), angles(NDIM)

      integer       :: interval
      real(kind=RP) :: alpha

      call FindInterval(frame, tq, interval, alpha)

      center = (1.0_RP - alpha) * frame%center(:,interval) + alpha * frame%center(:,interval+1)
      angles = (1.0_RP - alpha) * frame%angles(:,interval) + alpha * frame%angles(:,interval+1)
   end subroutine PoseAtTime

   subroutine VelocityAtTime(frame, tq, centerDot, omega)
      implicit none

      type(RigidMotionFrame_t), intent(in)  :: frame
      real(kind=RP),            intent(in)  :: tq
      real(kind=RP),            intent(out) :: centerDot(NDIM), omega(NDIM)

      integer       :: interval
      real(kind=RP) :: alpha, dtSegment
      real(kind=RP) :: angles(NDIM), angleRates(NDIM)

      centerDot = 0.0_RP
      omega = 0.0_RP

      if (tq < frame%time(1) .or. tq > frame%time(frame%n)) return

      call FindInterval(frame, tq, interval, alpha)
      dtSegment = frame%time(interval+1) - frame%time(interval)

      centerDot = (frame%center(:,interval+1) - frame%center(:,interval)) / dtSegment
      angles = (1.0_RP - alpha) * frame%angles(:,interval) + alpha * frame%angles(:,interval+1)
      angleRates = (frame%angles(:,interval+1) - frame%angles(:,interval)) / dtSegment

      call EulerRatesToOmega(angles, angleRates, omega)
   end subroutine VelocityAtTime

   subroutine FindInterval(frame, tq, interval, alpha)
      implicit none

      type(RigidMotionFrame_t), intent(in)  :: frame
      real(kind=RP),            intent(in)  :: tq
      integer,                  intent(out) :: interval
      real(kind=RP),            intent(out) :: alpha

      integer :: low, high, mid

      if (tq <= frame%time(1)) then
         interval = 1
         alpha = 0.0_RP
         return
      else if (tq >= frame%time(frame%n)) then
         interval = frame%n - 1
         alpha = 1.0_RP
         return
      end if

      low = 1
      high = frame%n
      do while (high - low > 1)
         mid = (low + high) / 2
         if (frame%time(mid) <= tq) then
            low = mid
         else
            high = mid
         end if
      end do

      interval = low
      alpha = (tq - frame%time(interval)) / (frame%time(interval+1) - frame%time(interval))
   end subroutine FindInterval

   subroutine EulerToRotationMatrix(angles, R)
      implicit none

      real(kind=RP), intent(in)  :: angles(NDIM)
      real(kind=RP), intent(out) :: R(NDIM,NDIM)

      real(kind=RP) :: phi, theta, psi
      real(kind=RP) :: cphi, sphi, ctheta, stheta, cpsi, spsi

      phi   = angles(IX)
      theta = angles(IY)
      psi   = angles(IZ)

      cphi = cos(phi)
      sphi = sin(phi)
      ctheta = cos(theta)
      stheta = sin(theta)
      cpsi = cos(psi)
      spsi = sin(psi)

      R(IX,IX) =  cpsi * ctheta
      R(IX,IY) =  cpsi * stheta * sphi - spsi * cphi
      R(IX,IZ) =  cpsi * stheta * cphi + spsi * sphi

      R(IY,IX) =  spsi * ctheta
      R(IY,IY) =  spsi * stheta * sphi + cpsi * cphi
      R(IY,IZ) =  spsi * stheta * cphi - cpsi * sphi

      R(IZ,IX) = -stheta
      R(IZ,IY) =  ctheta * sphi
      R(IZ,IZ) =  ctheta * cphi
   end subroutine EulerToRotationMatrix

   subroutine EulerRatesToOmega(angles, angleRates, omega)
      implicit none

      real(kind=RP), intent(in)  :: angles(NDIM), angleRates(NDIM)
      real(kind=RP), intent(out) :: omega(NDIM)

      real(kind=RP) :: theta, psi, phiDot, thetaDot, psiDot

      theta = angles(IY)
      psi = angles(IZ)
      phiDot = angleRates(IX)
      thetaDot = angleRates(IY)
      psiDot = angleRates(IZ)

      omega(IX) = -thetaDot * sin(psi) + phiDot * cos(psi) * cos(theta)
      omega(IY) =  thetaDot * cos(psi) + phiDot * sin(psi) * cos(theta)
      omega(IZ) =  psiDot - phiDot * sin(theta)
   end subroutine EulerRatesToOmega

   subroutine SetIdentity(A)
      implicit none

      real(kind=RP), intent(out) :: A(NDIM,NDIM)

      integer :: i

      A = 0.0_RP
      do i = 1, NDIM
         A(i,i) = 1.0_RP
      end do
   end subroutine SetIdentity

   subroutine ReadVector3(value, label, vector)
      implicit none

      character(len=*), intent(in)  :: value, label
      real(kind=RP),    intent(out) :: vector(NDIM)

      character(len=1024) :: cleanValue
      integer             :: i, io, pos1, pos2, commaCount

      cleanValue = TrimQuotes(value)
      pos1 = index(cleanValue, "[")
      pos2 = index(cleanValue, "]")

      if (pos1 == 0 .or. pos2 == 0 .or. pos2 <= pos1) then
         call MotionFrameError("ERROR: "//trim(label)//" must be a vector like [x,y,z].")
      end if

      cleanValue = cleanValue(pos1+1:pos2-1)
      commaCount = 0
      do i = 1, len_trim(cleanValue)
         if (cleanValue(i:i) == ",") then
            commaCount = commaCount + 1
            cleanValue(i:i) = " "
         end if
      end do

      if (commaCount /= NDIM - 1) then
         call MotionFrameError("ERROR: "//trim(label)//" must have exactly three comma-separated values.")
      end if

      read(cleanValue, *, iostat=io) vector
      if (io /= 0) then
         call MotionFrameError("ERROR: "//trim(label)//" must have exactly three values.")
      end if
   end subroutine ReadVector3

   subroutine BuildReferenceRotation(referenceRotation, hasXAxis, xAxis, hasYAxis, yAxis, hasZAxis, zAxis)
      implicit none

      real(kind=RP), intent(out) :: referenceRotation(NDIM,NDIM)
      logical,       intent(in)  :: hasXAxis, hasYAxis, hasZAxis
      real(kind=RP), intent(in)  :: xAxis(NDIM), yAxis(NDIM), zAxis(NDIM)

      integer       :: nAxes
      real(kind=RP) :: ex(NDIM), ey(NDIM), ez(NDIM)

      call SetIdentity(referenceRotation)

      nAxes = 0
      if (hasXAxis) nAxes = nAxes + 1
      if (hasYAxis) nAxes = nAxes + 1
      if (hasZAxis) nAxes = nAxes + 1

      if (nAxes == 0) return
      if (nAxes < 2) then
         call MotionFrameError("ERROR: A rigid motion reference frame needs at least two axis directions.")
      end if

      if (hasXAxis .and. hasYAxis) then
         ex = NormalizeAxis(xAxis, "reference x axis")
         ey = yAxis - dot_product(yAxis, ex) * ex
         ey = NormalizeAxis(ey, "reference y axis")
         ez = NormalizeAxis(CrossProduct(ex, ey), "reference z axis")

         if (hasZAxis) call CheckAxisSense(ez, zAxis, "reference z axis")

      else if (hasXAxis .and. hasZAxis) then
         ex = NormalizeAxis(xAxis, "reference x axis")
         ez = zAxis - dot_product(zAxis, ex) * ex
         ez = NormalizeAxis(ez, "reference z axis")
         ey = NormalizeAxis(CrossProduct(ez, ex), "reference y axis")

      else
         ey = NormalizeAxis(yAxis, "reference y axis")
         ez = zAxis - dot_product(zAxis, ey) * ey
         ez = NormalizeAxis(ez, "reference z axis")
         ex = NormalizeAxis(CrossProduct(ey, ez), "reference x axis")
      end if

      referenceRotation(:,IX) = ex
      referenceRotation(:,IY) = ey
      referenceRotation(:,IZ) = ez
   end subroutine BuildReferenceRotation

   function NormalizeAxis(axis, label) result(unitAxis)
      implicit none

      real(kind=RP),    intent(in) :: axis(NDIM)
      character(len=*), intent(in) :: label
      real(kind=RP)                :: unitAxis(NDIM)

      real(kind=RP) :: axisNorm

      axisNorm = sqrt(sum(axis * axis))
      if (axisNorm <= 100.0_RP * epsilon(1.0_RP)) then
         call MotionFrameError("ERROR: "//trim(label)//" cannot be zero or parallel to another axis.")
      end if

      unitAxis = axis / axisNorm
   end function NormalizeAxis

   subroutine CheckAxisSense(referenceAxis, candidateAxis, label)
      implicit none

      real(kind=RP),    intent(in) :: referenceAxis(NDIM), candidateAxis(NDIM)
      character(len=*), intent(in) :: label

      real(kind=RP) :: candidateUnit(NDIM)

      candidateUnit = NormalizeAxis(candidateAxis, label)
      if (dot_product(referenceAxis, candidateUnit) <= 0.0_RP) then
         call MotionFrameError("ERROR: "//trim(label)//" gives a left-handed rigid motion reference frame.")
      end if
   end subroutine CheckAxisSense

   function CrossProduct(a, b) result(c)
      implicit none

      real(kind=RP), intent(in) :: a(NDIM), b(NDIM)
      real(kind=RP)             :: c(NDIM)

      c(IX) = a(IY) * b(IZ) - a(IZ) * b(IY)
      c(IY) = a(IZ) * b(IX) - a(IX) * b(IZ)
      c(IZ) = a(IX) * b(IY) - a(IY) * b(IX)
   end function CrossProduct

   function StripInlineComment(line) result(cleanLine)
      implicit none

      character(len=*), intent(in) :: line
      character(len=len(line))     :: cleanLine

      integer :: pos

      cleanLine = line
      pos = index(cleanLine, "!")
      if (pos > 0) cleanLine = cleanLine(1:pos-1)
   end function StripInlineComment

   function TrimQuotes(value) result(cleanValue)
      implicit none

      character(len=*), intent(in) :: value
      character(len=LINE_LENGTH)   :: cleanValue

      integer :: n

      cleanValue = adjustl(value)
      cleanValue = trim(cleanValue)
      n = len_trim(cleanValue)

      if (n >= 2) then
         if ((cleanValue(1:1) == '"' .and. cleanValue(n:n) == '"') .or. &
             (cleanValue(1:1) == "'" .and. cleanValue(n:n) == "'")) then
            if (n == 2) then
               cleanValue = ""
            else
               cleanValue = cleanValue(2:n-1)
            end if
         end if
      end if
   end function TrimQuotes

   subroutine MotionFrameError(message)
      implicit none

      character(len=*), intent(in) :: message

      write(STD_OUT,'(A)') trim(message)
      stop 1
   end subroutine MotionFrameError

end module RigidMotionFrameClass
