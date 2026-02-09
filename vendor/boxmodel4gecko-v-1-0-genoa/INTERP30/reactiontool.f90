MODULE reactiontool
IMPLICIT NONE

! structure for storing reagents of each reactions
TYPE :: rgntmap
  INTEGER             :: nrg       ! number of reagents (in current rx)
  INTEGER,ALLOCATABLE :: idrg(:)   ! ID of the reagents (in current rx)
  REAL,ALLOCATABLE    :: stoirg(:) ! stoichiometric coef. of each product
END TYPE
TYPE(rgntmap),ALLOCATABLE :: rxrgnt(:) ! reagent list for rx [i]


! structure for storing products of each reactions
TYPE :: pdctmap
  INTEGER             :: npd       ! number of products (in current rx)
  INTEGER,ALLOCATABLE :: idpd(:)   ! ID of the products (in current rx)
  REAL,ALLOCATABLE    :: stoipd(:) ! stoichiometric coef. of each product
END TYPE
TYPE(pdctmap),ALLOCATABLE :: rxpdct(:) ! product list for rx [i]

CONTAINS

!SUBROUTINE readreac(lout,line,numsp,chrsp,isortlk,sortsp, etc ...
!=======================================================================
! PURPOSE: Read a reaction string to find reactants, products,  
! stoichiometric coefficients and Arrhenius parameters. Store the data 
! in the corresponding tables. 
!=======================================================================
SUBROUTINE readreac(lout,line,numsp,chrsp,isortlk,sortsp, &
                    numre,maxrg,maxpd,itype,arrhcf,reactkey,lostop)
  USE keywordlist, ONLY:mxlkey,mxkeywd,keywdlist,tabcfidx
  USE intptool, ONLY: stringarg, search
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: lout            ! unit output file
  CHARACTER(LEN=*),INTENT(IN) :: line   ! reaction string to be read
  INTEGER,INTENT(IN) :: numsp           ! number of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)  ! list of species
  INTEGER,INTENT(IN) :: isortlk(:)         ! index link "sorted" -> "mechanism" species
  CHARACTER(LEN=*),INTENT(IN) :: sortsp(:) ! sorted list of species
  INTEGER,INTENT(INOUT) :: numre       ! number of reactions
  INTEGER,INTENT(INOUT) :: maxrg       ! max # of reagent identified in a reaction
  INTEGER,INTENT(INOUT) :: maxpd       ! max # of product identified in a reaction
  INTEGER,INTENT(OUT) :: itype(:)      ! "type" of the reaction (i.e. keyword used)
  REAL,INTENT(OUT)    :: arrhcf(:,:)   ! arrhenius coef. for the reactions
  INTEGER,INTENT(OUT) :: reactkey(:)   ! flag to check if keyword is used
  LOGICAL,INTENT(OUT) :: lostop        ! error flag

  CHARACTER(LEN=mxlkey) :: tempkey      ! mxlkey: see keywordlist
  CHARACTER(LEN=LEN(chrsp(1))) :: tempsp
  CHARACTER(LEN=1) :: tchar
  INTEGER :: i,j,iloc,ios,iarg,iterm
  INTEGER :: ik, isep, isize, idspe
  INTEGER :: chknum, nleft, nright
  INTEGER :: spemax(1), speind
  REAL    :: xvalue, xstoi
  LOGICAL :: lonumber

  INTEGER, PARAMETER :: maxarg=100   ! maximum number of argument in a line
  INTEGER :: posi(maxarg,2)          ! position of the arguments on string
  INTEGER :: narg                    ! # of arguments

  CHARACTER(LEN=LEN(line)) :: lline, wlline ! left line
  CHARACTER(LEN=LEN(line)) :: rline, wrline ! right line
 
  ! maximum number of rgnt in rx to be read 
  ! (keep mxrgnt=2, unless you know what you are doing)
  INTEGER, PARAMETER :: mxrgnt=2 
  REAL    :: stoileft(mxrgnt), cpstoileft(mxrgnt)
  INTEGER :: speleft(mxrgnt), cpspeleft(mxrgnt)

  ! maximum number of pdct in rx to be read (arbitrarily large)
  INTEGER, PARAMETER :: mxpdct=50 
  REAL    :: stoiright(mxpdct), cpstoiright(mxpdct)
  INTEGER :: speright(mxpdct), cpsperight(mxpdct)

  INTEGER :: tabcfid(mxpdct), cktabcf(mxpdct)

  INTEGER :: mxre      ! max # of reactions
  INTEGER :: nidleft, nidright ! number of distinct ID
  INTEGER :: mlspe     ! max length of a species

! Initialize
!------------------
  mlspe=LEN(chrsp(1)) ! maximum length of a species
  mxre=SIZE(arrhcf,1) ! maximum # of reaction
  stoileft(:)=0.  ;  speleft(:)=0   
  stoiright(:)=0. ;  speright(:)=0  
  tabcfid(:)=0    ;  cktabcf(:)=0

! New reaction: increase reaction counter 
  numre=numre+1
  itype(numre)=0 ! default is no type
  IF (numre > mxre) THEN
    WRITE(lout,*) TRIM(line)
    WRITE(lout,*)'--error-- number of reactions > mxre, i.e. ',mxre
    lostop=.TRUE. ; RETURN
  ENDIF

! --------------------------
! READ ARRHENIUS PARAMETER
! --------------------------
  CALL stringarg(line,narg,posi)
  IF (narg < 6) THEN              ! 3 arrh. coefs + arrow + 1 reactant + 1 product
    WRITE(lout,*)'  --error--  missing arguments in the reactions'
    lostop=.TRUE. ; RETURN
  ENDIF

 ! Arrhenius coefficients are the last 3 arguments of the string
  DO i=1,3
    iarg=narg-3+i
    READ(line(posi(iarg,1):posi(iarg,2)),*,IOSTAT=ios)  xvalue
    IF (ios/=0) THEN
      WRITE(lout,*) 'line: ', TRIM(line)
      WRITE(lout,*)'--error--, while reading arrhenius coef number: ', i
      WRITE(lout,*)'          or a "/" might be missing for aux. info. '
      lostop=.TRUE. ; RETURN
    ENDIF
    arrhcf(numre,i)=xvalue
  ENDDO
  iterm=posi(narg-2,1)-1 ! set "pointer" 1 character before 1st Arrh. Coef.
 
! ----------------------------------------------
! READ THE REACTION (SPECIES AND STOICHIOMETRIE)
! ----------------------------------------------

! find the reaction arrow "=>" 
  iloc=INDEX(line,'=>')
  IF (iloc==0) THEN
    WRITE(lout,*)'--error-- line is not a reaction - missing "=>" '
    WRITE(lout,*)'          or a "/" might be missing for aux. info. '
    WRITE(lout,*) TRIM(line)
    lostop=.TRUE. ; RETURN
  ENDIF

! check that duplicate sign does not exist
  ik=INDEX(line(iloc+2:),'=>')
  IF (ik/=0) THEN
    WRITE(lout,*)'--error--  duplicate "=>" sign in line: '
    WRITE(lout,*) TRIM(line)
    lostop=.TRUE. ; RETURN
  ENDIF

! copy left and right reaction line (without kinetic parameters)
  lline=ADJUSTL(line(1:iloc-1))      ! reactive (left) side 
  rline=ADJUSTL(line(iloc+2:iterm))  ! product (right) side

! ----------------------------------
! READ THE LEFT SIDE OF THE REACTION  
! ----------------------------------

! make a working copy
  wlline=lline

! switch '+' into ' ' and find arguments location
  isep=LEN_TRIM(lline)
  DO i=1,isep
    IF (lline(i:i)=='+') wlline(i:i)=' '
  ENDDO

  CALL stringarg(wlline,narg,posi)
  IF (narg==0) THEN
    WRITE(lout,*)'  --error--  no species on the left side'
    WRITE(lout,*) TRIM(line)
    lostop=.TRUE. ; RETURN
  ENDIF

! initialize counters 
  nleft=0          ! number of species on the left side
  stoileft(:)=1.   ! stoi. coef. (default=1)
  speleft(:)=0     ! ID of the species on the left side
  reactkey(:)=0    ! keyword flag
  lonumber=.FALSE. ! logical for a number in the "stack" (a stoi. coef.)


! --- ENTRY LOOP - decrypt the arguments
  leftread:  DO iarg=1,narg

! check if argument is a number 
    tchar=wlline(posi(iarg,1):posi(iarg,1)) ! get 1st character
    chknum=SCAN(tchar,'0123456789.')        ! no negative number allowed
    IF (chknum.EQ.1) THEN  ! if argument is a number ...
      IF (lonumber) THEN   ! check that numbers are allowed
        WRITE(lout,*)'--error--, unexpected successive numbers: ', iarg
        WRITE(lout,*) 'left side: ', TRIM(lline)
        lostop=.TRUE. ; RETURN
      ELSE 
        lonumber=.TRUE.    ! raise flag for stoi. coef.
      ENDIF

      READ(wlline(posi(iarg,1):posi(iarg,2)),*,IOSTAT=ios)  xstoi
      IF (ios/=0) THEN
        WRITE(lout,*) 'left side: ', TRIM(lline)
        WRITE(lout,*)'--error--, while reading coef. number: ', iarg
        lostop=.TRUE. ; RETURN
      ENDIF
      CYCLE leftread
    ENDIF

! check if argument is a keyword 
    tempkey=wlline(posi(iarg,1):posi(iarg,2))
    DO ik=1,mxkeywd                   ! mxkeywd defined in "keywordlist"
      IF (tempkey==keywdlist(ik)) THEN
        IF (lonumber) THEN
          WRITE(lout,*) '--error--  stoi. coef. forbidden before keyword'
          WRITE(lout,*) 'left side: ', TRIM(lline)
          lostop=.TRUE. ; RETURN
        ELSE
          reactkey(ik)=reactkey(ik)+1 ! raise flag for the keyword
          CYCLE leftread
        ENDIF
      ENDIF
    ENDDO      

    IF (tempkey(1:8)=='NOTHING ') THEN  ! keyword 'NOTHING' not allowed on the left side
      WRITE(lout,*) '--error-- "NOTHING" not allowed on the left side:'
      WRITE(lout,*) 'left side: ', TRIM(lline)
      lostop=.TRUE. ; RETURN
    ENDIF

! check length of the species
    isize=posi(iarg,2)-posi(iarg,1)+1    
    IF (isize > mlspe) THEN
      WRITE(lout,*) '--error--  species length > mlspe, i.e. ', mlspe
      WRITE(lout,*) 'check left side species : ', TRIM(lline)
      lostop=.TRUE. ; RETURN
    ENDIF

! get species ID number
    tempsp=wlline(posi(iarg,1):posi(iarg,2))
    iloc = search(tempsp,sortsp,numsp)
    IF (iloc <= 0) THEN
      WRITE(lout,*) '--error--, species unknown: ',tempsp
      WRITE(lout,*) 'left side: ', TRIM(lline)
      lostop=.TRUE. ; RETURN
    ENDIF
    idspe=isortlk(iloc) 

! add newspecies
    nleft=nleft+1
    IF (nleft>SIZE(stoileft)) THEN
      WRITE(lout,*) '--error-- #species exceed max number: ',SIZE(stoileft)
      WRITE(lout,*) 'on the left side: ', TRIM(lline)
      lostop=.TRUE. ; RETURN
    ENDIF
    IF (lonumber) stoileft(nleft)=xstoi  ! store the coefficient
    speleft(nleft)= idspe                ! store the species number
    lonumber=.FALSE.

  ENDDO leftread  ! ----- end reading loop

! check that the stack for number is empty
  IF (lonumber) THEN
    WRITE(lout,*) '--error--, the last stoi. coef. is not linked to a species '
    lostop=.TRUE. ; RETURN
  ENDIF


! ------------------------------------
! READ THE RIGHT SIDE OF THE REACTION  
! ------------------------------------

! make a working copy
  wrline=rline

! switch '+' into ' ' and find arguments location
  isep=LEN_TRIM(rline)
  DO i=1,isep
    IF (rline(i:i)=='+') wrline(i:i)=' '
  ENDDO

  CALL stringarg(wrline,narg,posi)
  IF (narg==0) THEN
    WRITE(lout,*)'  --error--  no species on the right side'
    WRITE(lout,*)'             use keyword NOTHING if apropriate'
    WRITE(lout,*) TRIM(line)
    lostop=.TRUE. ; RETURN
  ENDIF

! reset counters 
  nright=0         ! number of species on the right side
  stoiright(:)=1.  ! stoi. coef. (default=1) 
  speright(:)=0    ! ID of the species on the right side
  lonumber=.FALSE. ! logical for a number in the "stack" (a stoi. coef.)

! --- ENTRY LOOP - decrypt the arguments
  rightread:  DO iarg=1,narg

! check if argument is a number 
    tchar=wrline(posi(iarg,1):posi(iarg,1))
    chknum=SCAN(tchar,'0123456789-.') ! negative stoi. coef. allowed
    IF (chknum==1) THEN    ! argument is a number
      IF (lonumber) THEN     ! check that a number is allowed
        WRITE(lout,*)'--error--, unexpected successive numbers: ', iarg
        WRITE(lout,*) 'right side: ', TRIM(rline)
        lostop=.TRUE. ; RETURN
      ELSE 
        lonumber=.TRUE. ! raise flag for stoi. coef.
      ENDIF

      READ(wrline(posi(iarg,1):posi(iarg,2)),*,IOSTAT=ios)  xstoi
      IF (ios/=0) THEN
        WRITE(lout,*)'--error--, while reading coef number: ', iarg
        WRITE(lout,*) 'right side: ', TRIM(rline)
        lostop=.TRUE. ; RETURN
      ENDIF
      CYCLE rightread
    ENDIF

! check if argument is NOTHING. Other keywords forbidden (interpreted as species)
    tempkey=wrline(posi(iarg,1):posi(iarg,2))
    IF (tempkey=='NOTHING ') THEN
      IF (lonumber) THEN
        WRITE(lout,*) '--error--  stoi. coef. not allowed before keyword'
        WRITE(lout,*) 'right side: ', TRIM(rline)
        lostop=.TRUE. ; RETURN
      ELSE
        CYCLE rightread  ! do nothing
      ENDIF
    ENDIF

! check length of the species
    isize=posi(iarg,2)-posi(iarg,1)+1
    IF (isize > mlspe) THEN
      WRITE(lout,*) '--error--  species length > mlspe, i.e. ', mlspe
      WRITE(lout,*) 'right side: ', TRIM(rline)
      lostop=.TRUE. ; RETURN
    ENDIF

! get species ID number
    tempsp=wrline(posi(iarg,1):posi(iarg,2))
    iloc = search(tempsp,sortsp,numsp)
    IF (iloc <= 0) THEN
      WRITE(lout,*) '--error--, unknow species: ',tempsp
      WRITE(lout,*) 'right side: ', TRIM(rline)
      lostop=.TRUE. ; RETURN
    ENDIF
    idspe=isortlk(iloc)

! add newspecies
    nright=nright+1
    IF (nright>SIZE(stoiright)) THEN
      WRITE(lout,*) '--error-- #species exceed max number: ',SIZE(stoiright)
      WRITE(lout,*) 'on the right side: ', TRIM(rline)
      lostop=.TRUE. ; RETURN
    ENDIF
    IF (lonumber) stoiright(nright)=xstoi  ! store the coefficient
    speright(nright)= idspe                ! store the species number
    lonumber=.FALSE.

  ENDDO rightread ! end of the reading loop

! check that the stack for number is empty
  IF (lonumber) THEN  
    WRITE(lout,*) '--error--, the last stoi. coef. is not linked to a species'
    lostop=.TRUE. ; RETURN
  ENDIF

! -----------------------------------------------
! CLEAN AND SORT THE DATA 
! -----------------------------------------------

! for tabulated stoi. coef. (TABCF keyword), the order of the species
! must not be changed. Make a copy of the right ID order before sorting
  IF (reactkey(tabcfidx)>0) tabcfid(:)=speright(:)  ! tabcfidx: see keywordlist

! search duplicate species (e.g. A+A => X+X)
  nidleft=nleft 
  DO i=1,nleft-1
    DO j=i+1,nleft
      IF (speleft(i) == speleft(j)) THEN
        speleft(j)=0
        stoileft(i)=stoileft(i)+stoileft(j)
        nidleft=nidleft-1
      ENDIF
    ENDDO
  ENDDO 

  nidright=nright
  DO i=1,nright-1
    DO j=i+1,nright
      IF (speright(i) == speright(j)) THEN
        IF (speright(i) /= 0) THEN
          speright(j)=0
          stoiright(i)=stoiright(i)+stoiright(j)
          nidright=nidright-1
        ENDIF
      ENDIF
    ENDDO
  ENDDO 

! sort the species according to order in "SPECIES" declaration
! lowest ID number set to the lowest index
  cpspeleft(:)=0 ; cpstoileft(:)=0.
  DO i=nidleft,1,-1
    spemax=MAXLOC(speleft) ; speind=spemax(1) ! spemax is of rank 1
    IF (speind == 0) THEN
      WRITE(lout,*) '--error--, speind=0 while sorting left side ID numbers'
      lostop=.TRUE. ; RETURN
    ENDIF 
    cpspeleft(i)=speleft(speind)
    cpstoileft(i)=stoileft(speind)
    speleft(speind)=0 ! reset species
  ENDDO
  IF (MAXVAL(speleft) > 0) THEN
    WRITE(lout,*) '--error--, species remain after sorting left side ID numbers'
    WRITE(lout,*) '# ID left side=',nidleft
    WRITE(lout,*) 'ID remaining=', speleft
    lostop=.TRUE. ; RETURN
  ENDIF

  cpsperight(:)=0 ; cpstoiright(:)=0.
  DO i=nidright,1,-1
    spemax=MAXLOC(speright) ; speind=spemax(1) ! spemax is of rank 1
    IF (speind == 0) THEN
      WRITE(lout,*) '--error--, speind=0 after sorting right side ID numbers'
      lostop=.TRUE. ; RETURN
    ENDIF 
    cpsperight(i)=speright(speind)
    cpstoiright(i)=stoiright(speind)
    speright(speind)=0 ! reset species
  ENDDO
  IF (MAXVAL(speright) > 0) THEN
    WRITE(lout,*) '--error--, species remain while sorting right side ID numbers'
    WRITE(lout,*) '# ID right side=',nidright
    WRITE(lout,*) 'ID remaining=', speright
    lostop=.TRUE. ; RETURN
  ENDIF

! check that ID order was not changed for reactions with tabulated stoi. coef.
  IF (reactkey(tabcfidx)>0) THEN
    cktabcf(:)=0
    WHERE (tabcfid /= cpsperight) cktabcf=1
    IF (SUM(cktabcf) /= 0) THEN
      WRITE(lout,*) '--error--, order of the species in a reaction using'
      WRITE(lout,*) 'TABCF keyword must match order given in SPECIES'
      lostop=.TRUE. ; RETURN
    ENDIF
  ENDIF

! -----------------------------------------------
! SET THE VALUES TO THE OUTPUT VARIABLES 
! -----------------------------------------------

  ! reagent stuff
  ALLOCATE(rxrgnt(numre)%idrg(nidleft), rxrgnt(numre)%stoirg(nidleft) )
  rxrgnt(numre)%nrg = nidleft
  rxrgnt(numre)%idrg(:) = cpspeleft(1:nidleft)
  rxrgnt(numre)%stoirg(:) = cpstoileft(1:nidleft)
  IF (nidleft>maxrg)  maxrg=nidleft    ! save mx reagent encountered

  ! product stuff
  ALLOCATE(rxpdct(numre)%idpd(nidright), rxpdct(numre)%stoipd(nidright) )
  rxpdct(numre)%npd = nidright
  rxpdct(numre)%idpd(:) = cpsperight(1:nidright)
  rxpdct(numre)%stoipd(:) = cpstoiright(1:nidright)
  IF (nidright>maxpd)  maxpd=nidright  ! save mx product encountered

END SUBROUTINE readreac

END MODULE reactiontool
