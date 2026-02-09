MODULE auxtool
IMPLICIT NONE
CONTAINS

! SUBROUTINE readaux(line,lout,loaux,xauxcf,auxflag,lostop)
! SUBROUTINE chkkeypar (lout, numre, auxflag, xauxcf, loaux, etc ...

!=======================================================================
! PURPOSE: Read the line providing the information (numbers) related  
! to a keyword. The format of the line is : 
!   KEYWORD / value1 value2 ... /
!
! NOTE : Only one keyword per line is allowed
!=======================================================================
SUBROUTINE readaux(line,lout,loaux,xauxcf,auxflag,lostop)
  USE keywordlist, ONLY: extraidx,foidx,hvidx,isoidx,tabcfidx,winidx, &
                         aouidx,wouidx,nauxarg                
  USE intptool, ONLY: stringarg  
  IMPLICIT NONE

  CHARACTER(lEN=*),INTENT(IN) :: line  ! line to be read
  INTEGER,INTENT(IN)    :: lout        ! unit output file
  LOGICAL,INTENT(INOUT) :: loaux       ! flag raised if aux. info provided
  REAL,INTENT(OUT)      :: xauxcf(:,:) ! auxiliary coefficients 
  INTEGER,INTENT(OUT)   :: auxflag(:)  ! flag related to each keyword
  LOGICAL,INTENT(OUT)   :: lostop      ! error flag

  INTEGER :: i,ii,ios,ipos
  INTEGER :: lenkey,idaux,narg
  REAL    :: xcoeff
  LOGICAL :: loslash,loflag 
  CHARACTER(lEN=LEN(line)) :: wline, chkline, emptyline
  INTEGER :: posi(SIZE(xauxcf,2),2)
  INTEGER :: mxaux                     ! max # of data in auxiliary info 

! Initialize
  mxaux=SIZE(xauxcf,2) 
  emptyline=' '
  IF (.NOT.loaux) THEN
    loaux=.TRUE.  ;  xauxcf(:,:)=0.  ;  auxflag(:)=0
  ENDIF

! cp the working line (left adjust) 
  wline=ADJUSTL(line)

! -----------------------------
! Find the keyword involved
! -----------------------------

! Set the "keyword ID number" (idaux) - see keywordlist
  lenkey=0
  IF (wline(1:2)=='HV')           THEN ; lenkey=2 ; idaux=hvidx
  ELSE IF (wline(1:7)=='FALLOFF') THEN ; lenkey=7 ; idaux=foidx
  ELSE IF (wline(1:5)=='EXTRA')   THEN ; lenkey=5 ; idaux=extraidx
  ELSE IF (wline(1:5)=='TABCF')   THEN ; lenkey=5 ; idaux=tabcfidx
  ELSE IF (wline(1:3)=='AOU')     THEN ; lenkey=3 ; idaux=aouidx
  ELSE IF (wline(1:3)=='WIN')     THEN ; lenkey=3 ; idaux=winidx
  ELSE IF (wline(1:3)=='WOU')     THEN ; lenkey=3 ; idaux=wouidx
  ELSE IF (wline(1:4)=='ISOM')    THEN ; lenkey=4 ; idaux=isoidx
  ENDIF

! Keyword not found => error.
  IF (lenkey == 0) THEN
    WRITE (lout,*) 'wline: ', TRIM(wline)
    WRITE (lout,*)' --error--, unidentified (or unexpected) keyword'
    lostop=.TRUE. ; RETURN
  ENDIF

! check that the keyword is not already used for the same reaction.
! Note : use of multiple keywords is not checked (done in "check key" routine) 
  IF (auxflag(idaux) /= 0) THEN
    WRITE(lout,*) 'wline: ', TRIM(wline)
    WRITE(lout,*)  '--error--, keyword already used for the reaction'
    lostop=.TRUE. ; RETURN
  ELSE
    auxflag(idaux)=1
  ENDIF

! ------------------------
! Read the numerical data 
! ------------------------

! find the data block (between 2 slahes) given after any keyword. 
  ii=lenkey+1
  wline=wline(ii:) ! clean line 
  ipos=INDEX(wline,'/') 
  loslash =.FALSE.
  IF (ipos == 0) THEN
    loslash=.TRUE.
  ELSE
    wline(ipos:ipos)=' ' ! rm 1st '/'
  ENDIF

  ipos=INDEX(wline,'/') 
  IF (ipos == 0) THEN
    loslash=.TRUE.
  ELSE
    wline(ipos:ipos)=' ' ! rm 2nd '/'
  ENDIF
 
! slashes not found => error
  IF (loslash) THEN
    WRITE(lout,*) 'line: ', TRIM(line)
    WRITE(lout,*)'--error--, expected "/" not found'
    lostop=.TRUE. ; RETURN
  ENDIF

! unexpected data remaining after the second slash (e.g. a 2nd data block)
  chkline(1:)=wline(ipos:)
  IF (chkline /= emptyline) THEN
    WRITE(lout,*) 'line: ', TRIM(line)
    WRITE(lout,*)'--error--, unexpected argument after the 2nd slash'
    lostop=.TRUE. ; RETURN
  ENDIF

! find the number of parameters between the slashes (narg). 
  CALL stringarg(wline,narg,posi)

! check parameters matches the expected numbers.
  IF (nauxarg(idaux) /= narg) THEN
    loflag=.TRUE.
! number of parameter for the keyword "EXTRA" is free (at least 1)
    IF (idaux==extraidx .AND. narg>=1 .AND. narg<=mxaux) loflag=.FALSE.
    IF (loflag) THEN
      WRITE(lout,*) 'line: ', TRIM(line)
      WRITE(lout,*)'--error--, inappropriate number of parameter'
      lostop=.TRUE. ; RETURN
    ENDIF
  ENDIF  

! read the numerical values between the slashes
  DO i=1,narg
    READ(wline(posi(i,1):posi(i,2)),*,IOSTAT=ios)  xcoeff
    IF (ios/=0) THEN
      WRITE(lout,*) 'line: ', TRIM(line)
      WRITE(lout,'(a,i3)')'--error--, while reading coef number:', i
      lostop=.TRUE. ; RETURN
    ENDIF
    xauxcf(idaux,i)=xcoeff
  ENDDO

END SUBROUTINE readaux

!=======================================================================
! PURPOSE: Check the consistency between the keyword used in the  
! reaction and auxiliary info then provided. Raise counters for  
! reaction type, store the labels, id numbers and auxiliary data in the 
! corresponding tables. 
! If an error is found => lostop is turned to true.     
!=======================================================================
SUBROUTINE chkkeypar (lout, numre, auxflag, xauxcf, loaux, &
                      reactkey, itype, num_n, numretype, id_n, &
                      idhv, ido2, id_m, idfo, idextra, &
                      idmeo2, idpero, idtabcf, &
                      idain, idaou, idwin, idwou, idiso, &
                      hvcf, focf, extracf, labtabcf, &
                      wincf,woucf,isocf,&
                      lostop)
  USE keywordlist, ONLY: nauxarg,hvidx,o2idx,tbidx,foidx,extraidx,meidx,&
                         p1idx,p9idx,tabcfidx,ainidx,aouidx,winidx, &
                         wouidx,isoidx
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: lout        ! unit output file
  INTEGER,INTENT(IN) :: numre       ! number of reactions
  INTEGER,INTENT(IN) :: auxflag(:)  ! keywords used in the auxiliary lines of "current" rxn
  REAL,INTENT(IN)    :: xauxcf(:,:) ! auxiliary coefficients set for the current reaction
  LOGICAL,INTENT(IN) :: loaux       ! flag raised if aux. info provided
  INTEGER,INTENT(IN) :: reactkey(:) ! flag to check if keyword is used

  INTEGER,INTENT(OUT)   :: itype(:)   ! "type" of the reaction (linked to keywords)
  INTEGER,INTENT(INOUT) :: num_n      ! # of "normal" reaction (i.e. not using a keyword) 
  INTEGER,INTENT(INOUT) :: numretype(:) ! # of reaction per keyword type
  INTEGER,INTENT(INOUT) :: id_n(:)    ! ID# of "normal" thermal reactions
  INTEGER,INTENT(OUT) :: idhv(:)      ! ID# of reactions using the keyword HV
  INTEGER,INTENT(OUT) :: ido2(:)      ! ID# of reactions using the keyword OXYGEN
  INTEGER,INTENT(OUT) :: id_m(:)      ! ID# of reactions using the keyword TBODY
  INTEGER,INTENT(OUT) :: idfo(:)      ! ID# of reactions using the keyword FALLOFF
  INTEGER,INTENT(OUT) :: idextra(:)   ! ID# of reactions using the keyword EXTRA
  INTEGER,INTENT(OUT) :: idmeo2(:)    ! ID# of reactions using the keyword MEPERO
  INTEGER,INTENT(OUT) :: idpero(:,:)  ! ID# of reactions using the keyword PEROx
  INTEGER,INTENT(OUT) :: idtabcf(:)   ! ID# of reactions using the keyword TABCF
  INTEGER,INTENT(OUT) :: idain(:)     ! ID# of reactions using the keyword AIN
  INTEGER,INTENT(OUT) :: idaou(:)     ! ID# of reactions using the keyword AOU
  INTEGER,INTENT(OUT) :: idwin(:)     ! ID# of reactions using the keyword WIN
  INTEGER,INTENT(OUT) :: idwou(:)     ! ID# of reactions using the keyword WOU
  INTEGER,INTENT(OUT) :: idiso(:)     ! isomerization reactions (vereckeen SAR) 
  REAL,INTENT(OUT) :: hvcf(:,:)       ! auxiliary info for "HV" reactions    
  REAL,INTENT(OUT) :: focf(:,:)       ! auxiliary info for "FALLOFF" reactions
  REAL,INTENT(OUT) :: extracf(:,:)    ! auxiliary info for "EXTRA" reactions
  REAL,INTENT(OUT) :: labtabcf(:,:)   ! auxiliary info for "TABCF" reactions
  REAL,INTENT(OUT) :: wincf(:,:)      ! auxiliary info for "WIN" reactions
  REAL,INTENT(OUT) :: woucf(:,:)      ! auxiliary info for "WOU" reactions
  REAL,INTENT(OUT) :: isocf(:,:)      ! auxiliary info for isomerization reaction
  LOGICAL,INTENT(OUT) :: lostop       ! error flag


! local
  INTEGER ikey, rerow
  INTEGER, DIMENSION(SIZE(nauxarg,1)) :: chkaux, chktype
  INTEGER :: sumchkaux, sumkey
  INTEGER :: rtype, idummy(1)

! --------------------------------------------------
! PRELIMINARY CHECKS
! --------------------------------------------------

! stop if no reaction provided while auxiliary information is used
  IF (numre == 0) THEN
    IF (loaux) THEN
      WRITE(lout,*) '--error-- auxiliary keyword provided before 1st reaction' 
      STOP '--error-- keyword found before 1st reaction'
    ENDIF
    RETURN ! 1st entry in the routine (no current reaction) -> escape
  ENDIF

! no keyword used - raise counter and return
  sumkey=SUM(reactkey)     ! does not count keyword NOTHING
  sumchkaux=SUM(auxflag)
  IF (sumkey==0) THEN
    IF (sumchkaux==0) THEN
      num_n=num_n+1
      id_n(num_n)=numre
      RETURN
    ELSE
      WRITE(lout,*) '--error-- auxiliary info provided without keyword in the reaction'
      lostop=.TRUE. ; RETURN
    ENDIF
  ENDIF

! current version allows only 1 keyword per reactions
  IF (sumkey>1) THEN
    WRITE(lout,*)'--error-- only 1 keyword allowed per reaction'
    lostop=.TRUE. ; RETURN
  ENDIF

! check consistency between reaction type and auxiliary information
  chkaux(:)=0
  WHERE (reactkey>0 .AND. nauxarg>0) 
    chkaux(:)=1
  END WHERE
  sumchkaux=SUM(chkaux)

  IF (sumchkaux == 0) THEN
    IF (loaux) THEN  ! auxiliary info provided but not expected 
      WRITE(lout,*)'--error--  unexpected auxiliary info provided'
      lostop=.TRUE. ; RETURN
    ENDIF
  ELSE
    IF (.NOT.loaux) THEN  ! auxiliary info expected but not found
      WRITE(lout,*)'--error--  reaction must have auxiliary info'
      lostop=.TRUE. ; RETURN
    ENDIF
  ENDIF

  chktype(:)=0
  WHERE (chkaux /= auxflag) ! aux. info is not for the reaction type
    chktype(:)=1
  END WHERE 
  sumchkaux=SUM(chktype)
  IF (sumchkaux/=0) THEN
    WRITE(lout,*)'--error-- mismatch in auxiliary and reaction keywords'
    lostop=.TRUE. ; RETURN
  ENDIF


! --------------------------------------------------
! SET REACTION TYPE AND STORE AUXILIARY INFORMATION
! --------------------------------------------------


! find which index equal 1 (only one left)
  idummy=MAXLOC(reactkey) ; rtype=idummy(1) ! idummy is of rank 1

  SELECT CASE (rtype)   ! rtype number provide reaction type 

! HV reaction
! ---------------
    CASE(hvidx)   
      IF (numretype(rtype) >= SIZE(idhv,1)) THEN
        WRITE(lout,*)'--error--   numhv > maxhv'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idhv(rerow)=numre
      hvcf(rerow,1)=xauxcf(rtype,1) ; hvcf(rerow,2)=xauxcf(rtype,2)
      RETURN

! OXYGEN reaction
! ---------------
    CASE(o2idx)
      IF (numretype(rtype) >= SIZE(ido2,1)) THEN
        WRITE(lout,*)'   --error-- numo2 > maxo2'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      ido2(rerow)=numre
      RETURN

! TBODY (3rd body) reaction
! --------------------------
    CASE(tbidx)
      IF (numretype(rtype) >= SIZE(id_m,1)) THEN
        WRITE(lout,*)'--error--   num_m > max_m'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      id_m(rerow)=numre
      RETURN

! FALL OFF reaction
! --------------------------
    CASE(foidx)
      IF (numretype(rtype) >= SIZE(idfo,1)) THEN
        WRITE(lout,*)'--error-- numfo > maxfo'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idfo(rerow)=numre
      focf(rerow,:)=xauxcf(rtype,:)
      RETURN

! EXTRA reaction
! --------------------------
    CASE(extraidx)
      IF (numretype(rtype) >= SIZE(idextra,1)) THEN
        WRITE(lout,*)'--error-- numextra > maxextra'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idextra(rerow)=numre
      extracf(rerow,:)=xauxcf(rtype,:)  
      RETURN

! MEPERO (RO2+CH3O2) reactions
! -----------------------------
    CASE(meidx)
      IF (numretype(rtype) >= SIZE(idmeo2,1)) THEN
        WRITE(lout,*)'--error-- nummeo2 > mxrpero'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idmeo2(rerow)=numre
      RETURN

! PERO-1-9 (RO2+PEROX) reactions
! -----------------------------
    CASE(p1idx:p9idx)
      ikey=rtype-p1idx+1 ! warning: pero type start @ index p1idx
      IF (numretype(rtype) >= SIZE(idpero,1)) THEN
        WRITE(lout,*)'--error-- numpero > mxrpero, pero=',ikey
        lostop=.TRUE. ; RETURN
      ENDIF

      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idpero(rerow,ikey)=numre
      itype(numre)=rtype
      RETURN

! TABCF (tabulated stoi coef)
! -----------------------------
    CASE(tabcfidx)
      IF (numretype(rtype) >= SIZE(idtabcf,1)) THEN
        WRITE(lout,*)'--error-- numtabcf > maxtabcf'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idtabcf(rerow)=numre
      labtabcf(rerow,1)=xauxcf(rtype,1)
      labtabcf(rerow,2)=xauxcf(rtype,2)
      RETURN

! AIN (gas->aerosol mass transfer)
! -----------------------------
    CASE(ainidx)
      IF (numretype(rtype) >= SIZE(idain,1)) THEN
        WRITE(lout,*)'--error-- numain > maxeq'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idain(rerow)=numre
      RETURN

! AOU (aerosol->gas mass transfer)
! -----------------------------
    CASE(aouidx)
      IF (numretype(rtype) >= SIZE(idaou,1)) THEN
        WRITE(lout,*)'--error-- numaou > maxeq'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idaou(rerow)=numre
!      aoucf(rerow,:)=xauxcf(rtype,:)
      RETURN

! WIN (gas->wall mass transfer)
! -----------------------------
    CASE(winidx)
      IF (numretype(rtype) >= SIZE(idwin,1)) THEN
        WRITE(lout,*)'--error-- numwin > maxeq'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idwin(rerow)=numre
      wincf(rerow,:)=xauxcf(rtype,:)
      RETURN

! WOU (wall->gas mass transfer)
! -----------------------------
    CASE(wouidx)
      IF (numretype(rtype) >= SIZE(idwou,1)) THEN
        WRITE(lout,*)'--error-- numwou > maxeq'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idwou(rerow)=numre
      woucf(rerow,:)=xauxcf(rtype,:)
      RETURN

! ISOM (isomerization reactions, Vereecken based)
! -----------------------------
    CASE(isoidx)
      IF (numretype(rtype) >= SIZE(idiso,1)) THEN
        WRITE(lout,*)'--error-- # isomerization > maxiso'
        lostop=.TRUE. ; RETURN
      ENDIF

      itype(numre)=rtype
      numretype(rtype)=numretype(rtype)+1
      rerow=numretype(rtype)
      idiso(rerow)=numre
      isocf(rerow,:)=xauxcf(rtype,:)
      RETURN

! Unidentified reaction
! -----------------------------
    CASE DEFAULT
      WRITE(lout,*)'--error-- reaction type not identified'
      WRITE(lout,*)'unexpected rtype= ',rtype
      lostop=.TRUE. ; RETURN
   
  END SELECT  

END SUBROUTINE chkkeypar

END MODULE auxtool
