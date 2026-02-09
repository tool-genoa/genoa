MODULE phototool
!$ USE omp_lib
IMPLICIT NONE
INTEGER, PARAMETER :: mxsza=20        ! max # of tabulated solar zenith angle
INTEGER, PARAMETER :: mxcsc=mxsza*3   ! max # of cubic spline coef. (interpol. J vs sza)
INTEGER            :: nsza            ! # of tabulated solar zenith angle (sza)
REAL               :: xsza(mxsza)     ! value of the tabulated sza

INTEGER            :: njlab           ! # of distinct jlabel
REAL,ALLOCATABLE   :: ratj(:,:)       ! values of tabulated J values
REAL,ALLOCATABLE   :: cscj(:,:)       ! cubic spline coef. for J values

! structure for storing J labels and reaction numbers
TYPE :: jrxn
  INTEGER             :: label        ! J label 
  INTEGER             :: nrxn         ! # of rxn having the corresponding label
  INTEGER,ALLOCATABLE :: idrxn(:)     ! ID of the rxn having corresponding J label
END TYPE
TYPE(jrxn),ALLOCATABLE :: jmap(:)     ! map of all J labels

CONTAINS
!=======================================================================
! PURPOSE: Allocate memory to the jmap(:)%idrxn field and set the 
! overall fields of jmap. 
!=======================================================================
SUBROUTINE make_jmap(numhv,rxhvtag)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: numhv      ! number of reaction using HV keyword
  INTEGER,INTENT(IN) :: rxhvtag(:) ! tag (label) used for the ith reaction  

  INTEGER :: labelcount(njlab)     ! number of records of each label
  INTEGER :: jlabel(njlab)         ! value of the label recorded
  INTEGER :: lab,i,j,ndat
  
  labelcount(:)=0 ; jlabel(:)=0 ; jmap%label = 0
  
! record and count all labels 
! --------------------------- 
  ndat = 0
  hvloop: DO i=1,numhv
    lab=rxhvtag(i)

    ! check if label already exist. If yes, count and go to next
    DO j=1,ndat
      IF (lab == jlabel(j)) THEN
        labelcount(j)=labelcount(j)+1
        CYCLE hvloop 
      ENDIF
    ENDDO

    ! if that point is reached, label not recorded yet (=> add label)
    ndat=ndat+1 ; jlabel(ndat)=lab ; labelcount(ndat)=1

  ENDDO hvloop              
  
! allocate the required memory to jmap
! ------------------------------------ 
  DO i=1,njlab
    ALLOCATE(jmap(i)%idrxn(labelcount(i)))
  ENDDO

! fill jmap with the mechanism data
! --------------------------------- 
  jmap%label = 0 ; jmap%nrxn = 0
  DO i=1,njlab
    jmap(i)%label = jlabel(i)
  ENDDO

  DO i=1,numhv
    lab=rxhvtag(i)
    DO j=1,njlab
      IF (lab==jmap(j)%label) THEN
        jmap(j)%nrxn = jmap(j)%nrxn + 1
        jmap(j)%idrxn(jmap(j)%nrxn) = i
        EXIT
      ENDIF
    ENDDO
  ENDDO
  
END SUBROUTINE make_jmap

!=======================================================================  
! PURPOSE: Read tabulated J values as a function of the solar zenith 
! angle (sza) for each "label" associated to the photolytic reactions 
! in the mechanism. A cubic spline interpolation is performed, to be  
! used during time integration.
! NOTE: The tabulated sza must be identical for each label. Value are
! stored in xsza(:), declared on top on the module
!=======================================================================  
SUBROUTINE readjvalue(lout)
  USE boxtool, ONLY: stoperr
  USE stringdata, ONLY: stringarg
  USE cubicspline, ONLY: cspline
  USE parameter_mod, ONLY: phot_file
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: lout          ! file unit (report)

  CHARACTER(LEN=76)    :: line,cpline
  CHARACTER(LEN=4)     :: keyword
  REAL    :: val1(SIZE(xsza,1))
  REAL    :: val2(SIZE(xsza,1))
  REAL    :: yy(SIZE(xsza,1)),cpho(SIZE(cscj,2))
  LOGICAL :: lostop
  INTEGER :: ii,filu,ilabel,i,ij,ios
  INTEGER :: isza
  LOGICAL :: lofirst    
  INTEGER :: chklabel(njlab)

! parameter for subroutine stringarg (see stringdata module)   
  INTEGER, PARAMETER :: mxarg=10      ! maximum # of arguments in a line
  INTEGER :: posi(mxarg,2)
  INTEGER :: narg

  CHARACTER(LEN=10),PARAMETER :: progname='readjvalue'
  CHARACTER(LEN=80) :: mesg1, mesg2  

  lostop=.FALSE. ; lofirst=.TRUE.
  chklabel(:)=0  ; ratj(:,:)=0.
  val1(:)=0.     ; val2(:)=0.
  nsza=0         ; xsza(:)=0.

! open file of input data
! ------------------------
  filu=12 
  OPEN(filu,FILE=phot_file,FORM='formatted',STATUS='old',IOSTAT=ios)
  IF (ios /= 0) THEN
    mesg1="File "//TRIM(phot_file)//" cannot be open"
    CALL stoperr(progname,mesg1)
  ENDIF

! -----------------------------
! READ NEXT LINE
! -----------------------------
  rdline: & 
  DO                 ! escape the loop if keyword 'END ' is found 
    READ(filu,'(a4,(a))',IOSTAT=ios) keyword, line
    IF (ios/=0) THEN 
      mesg1="End of file reached, keyword 'END' missing? "
      CALL stoperr(progname,mesg1)
    ENDIF

! comment line - read next
    IF (keyword(1:1) == '/' .OR. keyword(1:1) == '!') CYCLE rdline

! CASE SELECTOR FOR KEYWORDS
! -----------------------------
    SELECT CASE (keyword)   

    ! End of file - escape
    CASE('END ')
      EXIT rdline 

    ! Read the photolytic data
    CASE('PHOT')

      CALL stringarg(line,narg,posi)
      IF (narg /= 2 .AND. narg /= 3) THEN    ! species label nteta
        lostop=.TRUE.
        WRITE(lout,*) "--error--, in readjvalue while reading phot"
        WRITE(lout,*) "           # of arg not equal 3 or 4:"
        WRITE(lout,*) keyword,' ',TRIM(line) 
        CYCLE rdline
      ENDIF
      cpline=line(posi(narg-1,1):)   ! species is not read

      READ(cpline,*,IOSTAT=ios) ilabel, isza
      IF (ios/=0) THEN
        WRITE(lout,*) "--error--, in readjvalue while reading phot:"
        WRITE(lout,*) keyword,' ',TRIM(line) 
        lostop=.TRUE.
      ENDIF

      ! read the data for the set
      DO i=1,isza
        READ(filu,*) val1(i),val2(i)
      ENDDO

      ! check values of provided (must be identical in each set)
      IF (lofirst) THEN
        nsza=isza ; lofirst=.FALSE. ; xsza(:)=val1(:)
        IF (nsza>mxsza) THEN
          mesg1="number of tabulated sza are greater than mxsza"
          CALL stoperr(progname,mesg1)
        ENDIF
      ELSE
        IF (isza /= nsza) THEN
          mesg1="number of tabulated sza are not identical in "//TRIM(phot_file)
          CALL stoperr(progname,mesg1)
        ENDIF
        DO i=1,nsza
          IF (val1(i)/=xsza(i)) THEN
            WRITE(lout,*) "--error--, in readjvalue while reading "//TRIM(phot_file)
            WRITE(lout,*) "sza does not match previous set at: ",TRIM(line)
            lostop=.TRUE.
          ENDIF
        ENDDO
      ENDIF

      ! add J & check that each label is given only once 
      labloop: DO ij=1,njlab
        IF (ilabel == jmap(ij)%label) THEN
          IF (chklabel(ij) == 0) THEN
            ratj(ij,:)=val2(:) ; chklabel(ij)=1 ; EXIT labloop
          ELSE
            WRITE(lout,'(a19)') "-- error-- in readjvalue"
            WRITE(lout,'(i5,a25)') ilabel," ident. present + 1 time"
            lostop=.TRUE.
          ENDIF
        ENDIF
      ENDDO labloop

    ! Default: keyword unknown (error)
    CASE DEFAULT
      WRITE(lout,*) "--error--,  while reading in "//TRIM(phot_file), keyword
      WRITE(lout,*) "            keyword unkown"
      lostop=.TRUE.

    END SELECT  
  ENDDO rdline
  CLOSE(filu)

! -----------------------------
! CHECK LABELS AVAILABILITY
! -----------------------------

! check that each label in the mechanism has J values
  DO ii=1,njlab
    IF (chklabel(ii) == 0) THEN
      WRITE(lout,*) "--error--, in readjvalue, label not found!" 
      WRITE(lout,*) "  => missing label: ",jmap(ii)%label
      lostop=.TRUE.
    ENDIF
  ENDDO

! stop if error
  IF (lostop) THEN
    mesg1="Miscellaneous errors found while reading "//TRIM(phot_file)
    mesg2="see the *.out file for details"
    CALL stoperr(progname,mesg1,mesg2)
  ENDIF

! ------------------------------------------------ 
! GET COEFFICIENTS FOR CUBIC SPLINE INTERPOLATION 
! ------------------------------------------------ 
  DO i=1,njlab
    yy(:)=ratj(i,:)
    CALL cspline(nsza,xsza,yy,cpho)
    cscj(i,:)=cpho(:)
  ENDDO

END SUBROUTINE readjvalue

!=======================================================================  
! PURPOSE: Set the J values of the photolytic reactions. J values are  
! tabulated for each labels. A cubic spline interpolation is performed  
! to get J value at the corresponding solar zenith angle (sza). 
! Rate coefficients (i.e. J values) are set as the 1st arrhenius 
! coefficient of the reactions.  
!=======================================================================  
SUBROUTINE getjvalue(numhv,idhv,hvfact,time,arrhcf)
  USE initnbound, ONLY: iyear,imonth,iday,sla,slo,tz
  USE cubicspline, ONLY: csplint
  USE zenith, ONLY: czen
  USE parameter_mod, ONLY:flag_ZA,fixedZA
  USE boxtool, ONLY: stoperr
  USE frozenstuff,   ONLY: envfix    ! genoa
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: numhv             ! # of 'HV' rx
  INTEGER,INTENT(IN) :: idhv(:)           ! mechanism rx index of the ith 'HV' reaction
  REAL,INTENT(IN)    :: hvfact(:)         ! scaling factor applied to reference J
  REAL,INTENT(IN)    :: time              ! time (modulo 1 day)
  REAL,INTENT(INOUT) :: arrhcf(:,:)       ! [i,3] arrhenius coef. for the ith rx

  REAL    :: ti                           ! greenwich mean time in hour (e.g. 13h30 = 13.5) 
  REAL    :: sza                          ! solar zenith angle
  REAL    :: cck(SIZE(cscj,2))
  REAL    :: yy(SIZE(xsza,1))
  REAL    :: ratecoef(SIZE(idhv,1))
  INTEGER :: i,j,ire,ierr
  REAL    :: valj
  REAL    :: cor                          ! scaling factor applied to all J values

  CHARACTER(LEN=10),PARAMETER :: progname='getjvalue'
  CHARACTER(LEN=80) :: mesg1 

! ---------------------------------
! compute the solar zenith angle (sza)
! ---------------------------------

! ti is the greenwich mean time (in hour) - use tz to change to local time.    
  ti=time/3600. + tz
  !CALL czen(sla,slo,iyear,imonth,iday,ti,sza)
  
  IF (envfix(1,4) > 0.) THEN
    i = MAX(1,MIN(24, CEILING(time/3600.))) ! time in hour
    sza = envfix(4,i) ! Read fixed value
  ELSE IF (flag_ZA) THEN
    sza=fixedZA   ! fixed value
  ELSE ! Calculate the solar zenith angle
    CALL czen(sla,slo,iyear,imonth,iday,ti,sza)
  ENDIF

! check whether night or day (no photolysis at night!)
  IF (sza > 90.0) THEN 
    DO i=1,numhv
      ire=idhv(i)
      arrhcf(ire,1)=0. ; arrhcf(ire,2)=0. ; arrhcf(ire,3)=0.
    ENDDO
    RETURN
  ENDIF

! -----------------------------------------
! GET J VALUES (cubic spline interpolation) 
! -----------------------------------------

! set J value for each label
  DO i=1,njlab
    yy(:) = ratj(i,:)     ! set J values @ tabulated sza
    cck(:)=cscj(i,:)      ! set the polynomial coef. in J value for label i
    CALL csplint(nsza,xsza,yy,cck,sza,valj,ierr) ! compute J value at the given sza
    IF (ierr /= 0) THEN
      mesg1="interpolation issue while computing J"
      CALL stoperr(progname,mesg1)
    ENDIF
    valj = ABS(valj)

    ! set the rate coef (J) for each rxn having the corresponding label
    DO j=1,jmap(i)%nrxn
      ratecoef(jmap(i)%idrxn(j)) = valj
    ENDDO

  ENDDO

! --------------------------------------
! SET J VALUES (1st arrhenius parameter)
! --------------------------------------
!  cor=0.252  ! scaling factor for blacklight in Camergie mellon chamber
  cor=1.0  ! scaling factor for no corrections  
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire) SHARED(arrhcf,idhv,ratecoef,hvfact,cor,numhv)
  DO i=1,numhv
    ire=idhv(i)
    arrhcf(ire,1) = MAX(ratecoef(i)*hvfact(i)*cor,0.)
    arrhcf(ire,2) = 0.  ; arrhcf(ire,3) = 0.
  ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE getjvalue

END MODULE phototool
