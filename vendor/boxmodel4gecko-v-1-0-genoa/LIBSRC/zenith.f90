MODULE zenith
  IMPLICIT NONE

CONTAINS
SUBROUTINE CZEN(LAT,LONG,IIYEAR,IMTH,IDAY,GMT,ZENITH)
! ----------------------------------------------------------------------
! Input/output of the original subroutine was adapted to the boxmodel
! and the code switched to fortran 90 by B. Aumont, summer 2016
! ----------------------------------------------------------------------
! Written by Sasha Madronich at NCAR 23 February 1989.
! Based on equations given by Paltridge and Platt [1976] "Radiative 
! Processes in Meteorology and Climatology", Elsevier, pp. 62,63.
! Originally from Spencer, J.W., 1972, Fourier series representation 
! of the  position of the sun, Search, 2:172.
! This subroutine calculates solar zenith and azimuth angles for a 
! particular time and location.  Must specify:
!
! INPUT:
!  LAT - latitude in decimal degrees 
!  LONG - longitude in decimal degrees
!  IDATE (not longer used) - Date at Greenwich 
!                          - specify year (19yy), month (mm), day (dd)
!                            format is six-digit integer:  yymmdd)
!  GMT  - Greenwich mean time - decimal military eg.
!              22.75 = 45 min after ten pm gmt
! OUTPUT  
!  Zenith
!  Azimuth
! NOTE:  this approximate program has no changes from year to year.
! ----------------------------------------------------------------------
  IMPLICIT NONE
  
  REAL,INTENT(IN)     :: LAT,LONG
  REAL,INTENT(IN)     :: GMT
  INTEGER, INTENT(IN) :: IIYEAR
  INTEGER, INTENT(IN) :: IMTH
  INTEGER, INTENT(IN) :: IDAY
  REAL, INTENT(OUT)   :: ZENITH
!  REAL, INTENT(OUT)   :: AZIMUTH  ! uncomment if azimuth is wanted

  INTEGER :: IMN(12)
  REAL    :: LBGMT,LZGMT

  REAL    :: RLT, D, TZ, RDECL, EQR, EQH, ZPT
  REAL    :: CSZ, ZR
  INTEGER :: IIY, IJD
  INTEGER :: I

  REAL, PARAMETER :: PI = 3.1415926535590
  REAL, PARAMETER :: DR = PI/180.

  IMN = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! day / month
  RLT = LAT*DR                                  ! convert to radians

! parse date
!  IIYEAR = IDATE/10000
!!!!  IYEAR = 19*100 + IIYEAR  ! IYEAR is not used
!  IMTH = (IDATE - IIYEAR*10000)/100
!  IDAY = IDATE - IIYEAR*10000 - IMTH*100

! identify and correct leap years
  IIY = (IIYEAR/4)*4
  IF (IIY == IIYEAR) IMN(2) = 29

! compute current day of year IJD = 1 to 365
  IJD = 0
  DO I=1,IMTH - 1
    IJD = IJD + IMN(I)
  ENDDO
  IJD = IJD + IDAY

! calculate fractional day from start of year:
  D = FLOAT(IJD-1) + GMT/24.

! Equation 3.8 for "day-angle"
  TZ = 2.*PI*D/365.

! Equation 3.7 for declination in radians
  RDECL = 0.006918 - 0.399912*COS(   TZ) + 0.070257*SIN(   TZ) &
                   - 0.006758*COS(2.*TZ) + 0.000907*SIN(2.*TZ) &  
                   - 0.002697*COS(3.*TZ) + 0.001480*SIN(3.*TZ)

! Equation 3.11 for Equation of time  in radians
  EQR = 0.000075 + 0.001868*COS(   TZ) - 0.032077*SIN(   TZ) &
                 - 0.014615*COS(2.*TZ) - 0.040849*SIN(2.*TZ)

! convert equation of time to hours:
  EQH = EQR*24./(2.*PI) 

! calculate local hour angle (hours):
  LBGMT = 12. - EQH - LONG*24./360 

! convert to angle from GMT
  LZGMT = 15.*(GMT - LBGMT)
  ZPT = LZGMT*DR

! Equation 2.4 for cosine of zenith angle 
  CSZ = SIN(RLT)*SIN(RDECL) + COS(RLT)*COS(RDECL)*COS(ZPT)
  ZR = ACOS(CSZ)
  ZENITH = ZR/DR

! calc local solar azimuth (ba: uncomment if wanted)
!  CAZ = (SIN(RDECL) - SIN(RLT)*COS(ZR))/(COS(RLT)*SIN(ZR))
!  RAZ = ACOS(CAZ)
!  AZIMUTH = RAZ/DR

  RETURN
END SUBROUTINE CZEN

END MODULE zenith
