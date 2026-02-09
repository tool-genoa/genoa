! ----------------------------------------------------------------------
! GECKO2BOX  -  main program 
! ----------------------------------------------------------------------
! Purpose: 
! provide the interface between geckoa mechanisms (or mechanisms 
! formatted in an identical way) and the box model. 
!                               -----
! The mechanism must have the following structure:
! SPECIES
!    list of species ... (1 name + 1 molec. weight per line)
! END
! REACTIONS
!    list of reactions ... 
! END
!                               -----
! the format of a reaction is:
!    A + B + KEYWORD (optional) => X + Y      a  b  c 
!            KEYWORD / data data .../
!
! Available keywords are in the module "keywordlist"
! a, b and c are the arrhenius coefficients
! ----------------------------------------------------------------------
PROGRAM gecko2box
  USE keywordlist, ONLY: mxkeywd,hvidx,nauxarg
  USE intptool, ONLY: quicksort
  USE speciestool, ONLY: ptrspecies,readspecies,chkspecies
  USE auxtool, ONLY: chkkeypar,readaux
  USE reactiontool, ONLY: readreac,rxpdct,rxrgnt
  USE mecatool, ONLY: chkmeca,hvinfo,countself
  IMPLICIT NONE

! size of tables (as read from the input file)
  INTEGER :: maxsp     ! max # of species 
  INTEGER :: maxre     ! max # of reactions
  INTEGER :: maxhv     ! max # of react. using 'HV' keyword
  INTEGER :: max_m     ! max # of react. using 'TBODY' keyword
  INTEGER :: maxextra  ! max # of react. using 'EXTRA' keyword
  INTEGER :: maxo2     ! max # of react. using 'OXYGEN' keyword
  INTEGER :: mxrpero   ! max # of react. using 'PEROX' or 'MEPERO' keyword
  INTEGER :: maxfo     ! max # of react. using 'FALLOFF' keyword
  INTEGER :: maxaeq    ! max # of gas <-> particles equilibrium
  INTEGER :: maxweq    ! max # of gas <-> wall equilibrium
  INTEGER :: maxiso    ! max # of react. using 'ISOM' keyword
  INTEGER :: maxtabcf  ! max # of react. using 'TABCF' keyword

! size of tables (set as parameters)
  INTEGER, PARAMETER :: maxlsp=8        ! max length of species names
  INTEGER, PARAMETER :: maxaux=7        ! max # of arguments in the auxiliary info. lines
  INTEGER, PARAMETER :: mxcpero=9       ! max # of peroxy classes
  INTEGER, PARAMETER :: mchromo=140     ! max # of labels in HV reactions
  
! length of string to be read (parameters)
  INTEGER, PARAMETER :: mxline=8        ! max # of lines to be used for a reaction
  INTEGER, PARAMETER :: lenread=160     ! string length (max) in the mechanism file
  INTEGER, PARAMETER :: lenline=lenread + (mxline-1)*100 ! string length after concatenation

! file units (parameters)
  INTEGER, PARAMETER :: lin=15      ! unit input file
  INTEGER, PARAMETER :: lout=16     ! unit output file
  INTEGER, PARAMETER :: bof=11      ! unit output (boxmod link) binary file

! species tables
  CHARACTER(LEN=maxlsp),ALLOCATABLE :: chrsp(:)   ! [size=maxsp] list of species (order = mechanism list)
  CHARACTER(LEN=maxlsp),ALLOCATABLE :: sortsp(:)  ! [size=maxsp] list of species (sorted)
  INTEGER,ALLOCATABLE               :: isortlk(:) ! [size=maxsp] corresponding link "sorted" -> "mechanism" index
  REAL,ALLOCATABLE                  :: wmol(:)    ! [size=maxsp] molar mass of the species

! index per reaction type (reaction number in the full mechanism)
  REAL,ALLOCATABLE    :: arrhcf(:,:) ! [size=maxre,3] arrhenius coef. for the reactions
  INTEGER,ALLOCATABLE :: itype(:)    ! [size=maxre] "type" of the reaction in the mechanism (0 is used for "normal")
  INTEGER,ALLOCATABLE :: id_n(:)     ! [size=maxre] "normal" thermal reactions
  INTEGER,ALLOCATABLE :: idhv(:)     ! [size=maxhv] photolytic reactions (keyword HV)
  INTEGER,ALLOCATABLE :: ido2(:)     ! [size=maxo2] reactions using O2 as reactant (keyword OXYGEN)
  INTEGER,ALLOCATABLE :: id_m(:)     ! [size=max_m] 3rd body reactions (keyword TBODY)
  INTEGER,ALLOCATABLE :: idfo(:)     ! [size=maxfo] reactions in the FALLOFF regime (keyword FALLOFF)
  INTEGER,ALLOCATABLE :: idextra(:)  ! [size=maxextra] special reaction (particular reaction rate - keyword EXTRA)
  INTEGER,ALLOCATABLE :: idmeo2(:)   ! [size=mxrpero] reactions using CH3O2 as reactant (keyword MEPERO)
  INTEGER,ALLOCATABLE :: idpero(:,:) ! [size=mxrpero,mxcpero] reactions using PEROx as reactant (keyword PEROx)
  INTEGER,ALLOCATABLE :: idtabcf(:)  ! [size=maxtabcf] reactions with tabuled stoi. coef.
  INTEGER,ALLOCATABLE :: idain(:)    ! [size=maxaeq] gas->aerosol mass transfer "reaction"
  INTEGER,ALLOCATABLE :: idaou(:)    ! [size=maxaeq] aerosol->gas mass transfer "reaction"
  INTEGER,ALLOCATABLE :: idwin(:)    ! [size=maxweq] gas->wall mass transfer "reaction"
  INTEGER,ALLOCATABLE :: idwou(:)    ! [size=maxweq] wall->gas mass transfer "reaction"
  INTEGER,ALLOCATABLE :: idiso(:)    ! [size=maxiso] isomerization reactions (vereckeen SAR) 

! auxiliary information (provided between / / and linked to a keyword)
  LOGICAL :: loaux                  ! flag (=true) if auxiliary info is available
  REAL,ALLOCATABLE :: hvcf(:,:)     ! [size=(maxhv,maxaux)] auxiliary info for "hv" reactions
  REAL,ALLOCATABLE :: focf(:,:)     ! [size=(maxfo,maxaux)] auxiliary info for "fall off" reactions 
  REAL,ALLOCATABLE :: extracf(:,:)  ! [size=(maxextra,maxaux)]  auxiliary info for "extra" reactions
  REAL,ALLOCATABLE :: labtabcf(:,:) ! [size=(maxtabcf,maxaux)]  auxiliary info for "tabulated stoi coef" reactions
  REAL,ALLOCATABLE :: wincf(:,:)    ! [size=(maxweq,maxaux)]  auxiliary info for gas->wall mass transfer
  REAL,ALLOCATABLE :: woucf(:,:)    ! [size=(maxweq,maxaux)]  auxiliary info for wall->gas mass transfer
  REAL,ALLOCATABLE :: isocf(:,:)    ! [size=(maxiso,maxaux)]  auxiliary info for isomerization reaction

! set names and properties of the keywords available to design the mechanism
! See also the subroutine keywordlist providing parameters
  INTEGER :: reactkey(mxkeywd)         ! identify which keywords are used in a given reaction
  REAL    :: xauxcf(mxkeywd,maxaux)=0. ! auxiliary coefficient of the "current" reaction
  INTEGER :: auxflag(mxkeywd)=0        ! keywords in auxiliary line of the "current" reaction

! Number of species and reactions 
  INTEGER :: numsp=0    ! number of species (incremented while reading the mechanism)
  INTEGER :: numre=0    ! number of reactions (incremented while reading the mechanism)
  INTEGER :: num_n=0    
  INTEGER :: numretype(mxkeywd)=0  ! number of reaction per keyword type

! index for starting/ending of phase contributors in the list of species
  INTEGER :: ngas                ! # of gas phase species
  INTEGER :: npart               ! # of part. phase species
  INTEGER :: nwall               ! # of wall phase species
  INTEGER :: ptrgas(2)           ! pointer (1st and last) gas phase species
  INTEGER :: ptrpart(2)          ! pointer (1st and last) part. phase species
  INTEGER :: ptrwall(2)          ! pointer (1st and last) wall phase species

! interface with the box model (not used within the gecko2box program)
  INTEGER :: numo2,nummeo2,num_m,numfo,numhv,numextra
  INTEGER :: numtabcf,numpero(mxcpero)
  INTEGER :: numain,numaou,numwin,numwou,numiso
  INTEGER :: numself
  INTEGER :: maxrg=0                 ! track max # of reagents found in a reaction
  INTEGER :: maxpd=0                 ! track max # of reagents found in a reaction


  INTEGER :: fg_fname=0              ! flag for reading filename
  CHARACTER(LEN=40) :: filename='-'

! other operating variables
  LOGICAL :: lostop         ! true if error found in the mechanism

! "local" (working) variable
  INTEGER :: nlines, ire, ios, last
  INTEGER :: i, k
  LOGICAL :: lolin, lopen
  CHARACTER(LEN=lenline) :: line, newline, oldline 

! ----------------
! USER INTERFACE
! ----------------
! Provide output information (may output large file for big mechanism)
! Turn to "false" if no information (or no "debugging" needed) 
  LOGICAL, PARAMETER :: lospeech=.FALSE. ! print info on screen 
  LOGICAL, PARAMETER :: loreport=.FALSE. ! print additional info in the output file 

! Checking extra large mechanism is time consuming (in chkmeca subroutine) and usually
! not necessary. Default is no deep check. Turn to false if check is wanted
  LOGICAL, PARAMETER :: lonodeep=.TRUE.

! ------------
! INITIALIZE
! ------------

! genoa read input argument
  i = IARGC()
  IF (i > 0) THEN
    CALL GETARG(1,line)
    ! Read filename - 1st argument
    READ(line,*,IOSTAT=ios) filename
    IF (ios /= 0) THEN
      PRINT*, 'error reading input argument from: ', TRIM(line)
      STOP
    ENDIF
    IF (TRIM(filename) /= '-') THEN ! check if filename is provided
      PRINT*, 'filename = ', filename
      fg_fname = 1
    ELSE
      PRINT*, 'filename not provided. Using default.'
    ENDIF

  ENDIF

! check that maxaux is greater than declaration in nauxarg ...
  IF (maxaux < MAXVAL(nauxarg)) THEN   
    PRINT*, '--error--, maxaux < # in nauxarg'
    STOP
  ENDIF
 
  loaux = .FALSE. 

! Open input and output (report) files
  IF (fg_fname > 0) THEN
    OPEN(lin , FILE=TRIM(filename)//'.mech', STATUS='OLD')
    OPEN(lout, FILE=TRIM(filename)//'.out')
  ELSE ! default
    OPEN(lin , FILE='indat.mech', STATUS='OLD')
    OPEN(lout, FILE='outdat.out') 
  ENDIF

! ------------
! READ VARIOUS SIZE
! ------------

! search keyword 'SIZE' 
  line=' '
  DO WHILE (line(1:5) /= 'SIZE ')
    READ(lin,'(a)', IOSTAT=ios) line
    IF (ios/=0) THEN 
      WRITE(*,*) '--error--, keyword SIZE not found'
      STOP 'error while reading input file - keyword SIZE missing ?'
    ENDIF
  ENDDO

  IF (lospeech) WRITE(*,*) 'reading sizes of tables ...'
  READ(lin,*) maxsp     ! ; WRITE(*,*) 'maxsp=',maxsp
  READ(lin,*) maxre     ! ; WRITE(*,*) 'maxre=',maxre   
  READ(lin,*) maxhv     ! ; WRITE(*,*) 'maxhv=',maxhv   
  READ(lin,*) max_m     ! ; WRITE(*,*) 'max_m=',max_m   

  READ(lin,*) maxo2     ! ; WRITE(*,*) 'maxo2=',maxo2   
  READ(lin,*) maxextra  ! ; WRITE(*,*) 'maxextra=',maxextra
  READ(lin,*) mxrpero   ! ; WRITE(*,*) 'mxrpero=',mxrpero 
  READ(lin,*) maxfo     ! ; WRITE(*,*) 'maxfo=',maxfo   

  READ(lin,*) maxiso    ! ; WRITE(*,*) 'maxiso=',maxiso  
  READ(lin,*) maxaeq    ! ; WRITE(*,*) 'maxaeq=',maxaeq
  READ(lin,*) maxweq    ! ; WRITE(*,*) 'maxweq=',maxweq
  READ(lin,*) maxtabcf  ! ; WRITE(*,*) 'maxtabcf=',maxtabcf

! overwrite 0 sizes (avoid to allocate table without size)
  IF (max_m==0)   max_m=1 
  IF (maxo2==0)   maxo2=1
  IF (maxextra==0) maxextra=1
  IF (maxtabcf==0) maxtabcf=1
  IF (mxrpero==0) mxrpero=1
  IF (maxfo==0)   maxfo=1
  IF (maxiso==0)  maxiso=1
  IF (maxweq==0)   maxweq=1
  IF (maxaeq==0)   maxaeq=1
    
! ---------------------------------------------------------------------
! ------      "SPECIES" DATA BLOCK - READ THE LIST OF SPECIES   -------  
! ---------------------------------------------------------------------

! allocate space needed for the "species" tables 
  ALLOCATE(chrsp(maxsp))   ; chrsp(:)=' '
  ALLOCATE(sortsp(maxsp))  ; sortsp(:)=' '
  ALLOCATE(isortlk(maxsp)) ; isortlk(:)=0
  ALLOCATE(wmol(maxsp))    ; wmol(:)=0.

  lostop=.FALSE.  ! if true: program stop  after reading full data block 

! search keyword 'SPECIES' 
  line=' '
  DO WHILE (line(1:7) /= 'SPECIES')
    READ(lin,'(a)', IOSTAT=ios) line
    IF (ios/=0) THEN 
      WRITE(lout,*) '--error--, keyword SPECIES not found'
      STOP 'error while reading input file - keyword SPECIES missing ?'
    ENDIF
  ENDDO

  lopen=.FALSE.
  ptrgas(:)=0 ; ptrpart(:)=0 ; ptrwall(:)=0

  IF (lospeech) PRINT*, 'reading the species ...'

! Start loop to read list of species 
  getspecies : DO                 ! escape loop if keyword END found 
    READ(lin,'(a)',IOSTAT=ios) line
    IF (ios/=0) THEN 
      WRITE(*,*) '--error--, keyword END not found'
      STOP 'error while reading SPECIES block - keyword END missing ?'
    ENDIF

    IF (line(1:11)=='END SPECIES') EXIT getspecies 

    IF (line(1:1)=='!') CYCLE       ! comment in the list of species 

    IF (line(1:6)=='PHASE:') THEN   ! start/stop phase membership
      CALL ptrspecies(lout,line,numsp,lopen,ptrgas,ptrpart,ptrwall)
      CYCLE
    ENDIF

    IF (lopen) THEN
      CALL readspecies(lout,line,numsp,chrsp,wmol,lostop)
    ELSE
      WRITE(*,*) 'expecting next species but no phase open'
      STOP ''
    ENDIF
  ENDDO getspecies 

  IF (lospeech) PRINT*, 'numsp= ',numsp

! sort the list of species (efficient search needed for large mechanism)
  IF (lospeech) PRINT*, 'sort the list of species ...'

  CALL quicksort(numsp,chrsp,isortlk,sortsp)
! other options for sorting (old version) 
!  CALL linuxsort(numsp,chrsp,isortlk,sortsp) ! efficient - but require file writing
!  CALL bubblesort(numsp,chrsp,isortlk,sortsp) ! too slow if numsp>1000
 
! check species names
  IF (lospeech) WRITE(*,*) 'checking the species ...'
  CALL chkspecies(lout,lospeech,numsp,chrsp,sortsp,isortlk,wmol, &
                  ptrgas,ptrpart,ptrwall, &
                  ngas,npart,nwall,lostop)
  IF (lospeech) THEN
    WRITE(*,*) 'ptrgas: ' ,ptrgas(1),  ptrgas(2),  ngas
    WRITE(*,*) 'ptrpart: ',ptrpart(1), ptrpart(2), npart
    WRITE(*,*) 'ptrwall: ',ptrwall(1), ptrwall(2), nwall
  ENDIF

! stop if error found while reading the list of species
  IF (lostop) THEN
    WRITE(lout,*)'--error-- while reading the list of species'
    CLOSE(lin) ; CLOSE(lout)
    STOP ' ERROR LEVEL 1 while reading species list => see *.out file'
  ENDIF

! ----------------------------------------------------------------------
! --------   "REACTIONS" DATA RECORD - read list of reactions   -------- 
! ----------------------------------------------------------------------

! allocate space needed for the "reaction" tables 
  IF (lospeech) WRITE(*,*) 'allocate reaction tables ...'
  ALLOCATE(arrhcf(maxre,3))        ; arrhcf(:,:)=0.    ! REAL, arrhenius coef. for the reactions

  ALLOCATE(rxrgnt(maxre))                              ! Structure storing the reaction reagent data 
  ALLOCATE(rxpdct(maxre))                              ! Structure storing the reaction product data 

  ALLOCATE(itype(maxre))             ; itype(:)=-1     ! INTEGER, rxn "type" in the mechanism (0 is used for "normal")
  ALLOCATE(id_n(maxre))              ; id_n(:)=0       ! INTEGER, "normal" thermal reactions
  ALLOCATE(idhv(maxhv))              ; idhv(:)=0       ! INTEGER, photolytic reactions (keyword HV)
  ALLOCATE(hvcf(maxhv,maxaux))       ; hvcf(:,:)=0.    ! REAL, auxiliary info for "hv" reactions   
  ALLOCATE(id_m(max_m))              ; id_m(:)=0       ! INTEGER, 3rd body reactions (keyword TBODY)
  ALLOCATE(ido2(maxo2))              ; ido2(:)=0       ! INTEGER, reactions using O2 as reactant (keyword OXYGEN)
  ALLOCATE(idextra(maxextra))        ; idextra(:)=0    ! INTEGER, special reaction (particular reaction rate - keyword EXTRA)
  ALLOCATE(extracf(maxextra,maxaux)) ; extracf(:,:)=0. ! REAL, auxiliary info for "tabulated stoi coef" reactions       
  ALLOCATE(idtabcf(maxtabcf))        ; idtabcf(:)=0    ! INTEGER, reactions with tabuled stoi. coef.
  ALLOCATE(idmeo2(mxrpero))          ; idmeo2(:)=0     ! INTEGER, reactions using CH3O2 as reactant (keyword MEPERO)
  ALLOCATE(idpero(mxrpero,mxcpero))  ; idpero(:,:)=0   ! INTEGER, reactions using PEROx as reactant (keyword PEROx)
  ALLOCATE(idfo(maxfo))              ; idfo(:)=0       ! INTEGER, reactions in the FALLOFF regime (keyword FALLOFF)
  ALLOCATE(focf(maxfo,maxaux))       ; focf(:,:)=0.    ! REAL, auxiliary info for "fall off" reactions       
  ALLOCATE(idiso(maxiso))            ; idiso(:)=0      ! INTEGER, isomerization reactions (vereckeen SAR) 
  ALLOCATE(isocf(maxiso,maxaux))     ; isocf(:,:)=0.   ! REAL, auxiliary info for isomerization reaction        
  ALLOCATE(idain(maxaeq))            ; idain(:)=0      ! INTEGER, gas->aerosol mass transfer "reaction"
  ALLOCATE(idaou(maxaeq))            ; idaou(:)=0      ! INTEGER, aerosol->gas mass transfer "reaction"
  ALLOCATE(idwin(maxweq))            ; idwin(:)=0      ! INTEGER, gas->wall mass transfer "reaction"
  ALLOCATE(idwou(maxweq))            ; idwou(:)=0      ! INTEGER, wall->gas mass transfer "reaction"
  ALLOCATE(wincf(maxweq,maxaux))     ; wincf(:,:)=0.   ! REAL, auxiliary info for gas->wall mass transfer      
  ALLOCATE(woucf(maxweq,maxaux))     ; woucf(:,:)=0.   ! REAL, auxiliary info for wall->gas mass transfer     

  lostop=.FALSE. ! stop if turned to true after reading full data block 

! search for keyword 'REACTIONS'
  line=' '
  DO WHILE (line(1:9) /= 'REACTIONS')
    READ(lin,'(a)', IOSTAT=ios) line
    IF (ios/=0) THEN 
      WRITE(*,*) '--error--, keyword REACTIONS not found'
      STOP '--error-- while reading input file - keyword REACTIONS missing?'
    ENDIF
  ENDDO

  IF (lospeech) PRINT*, 'reading the reactions ...'

! START LOOP to read reactions
  loopreac : DO                   ! escape the loop if keyword END found 
    READ(lin,'(a)',IOSTAT=ios) line
    IF (ios/=0) THEN 
      WRITE(*,*) '--error--, keyword END not found in REACTIONS data'
      STOP 'error while reading REACTIONS data - keyword END missing ?'
    ENDIF

! comment line in the mechanism
    IF (line(1:1)=='!') CYCLE loopreac
      
! keyword end => exit the loop
    IF (line(1:4)=='END ') THEN 
      CALL chkkeypar(lout, numre, auxflag, xauxcf, loaux, &
                      reactkey, itype, num_n, numretype, id_n, &
                      idhv, ido2, id_m, idfo, idextra, &
                      idmeo2, idpero, idtabcf, &
                      idain, idaou, idwin, idwou, idiso, &
                      hvcf, focf, extracf, labtabcf, &
                      wincf, woucf, isocf,&
                      lostop)

      EXIT loopreac 
    ENDIF

! check for auxiliary information => read data and cycle
    last=LEN_TRIM(line)
    IF ( line(last:last) == '/' ) THEN 
      IF (loreport) WRITE(lout,'(a)') TRIM(line)
      CALL readaux(line,lout,loaux,xauxcf,auxflag,lostop)
      CYCLE loopreac

! else read a reaction
    ELSE

! concatenation of a reaction written on 2 lines (or more if allowed - see mxline)
      lolin=.FALSE. ;  nlines=1
      IF ( line(last:last) == '+' ) lolin=.TRUE. 
      DO WHILE (lolin)
        nlines=nlines+1
        IF (nlines > mxline) THEN
          WRITE(lout,*)'--error-- too many input lines for a reaction'
          WRITE(lout,*)'nlines > mxline'
          lostop=.TRUE.
          CYCLE loopreac
        ENDIF

        READ(lin,'(a)', IOSTAT=ios) newline
        IF (ios/=0) THEN 
          WRITE(lout,*)'--error-- while reading reactions on multiple lines '
          STOP 'error while reading reactions on multiple lines '
        ENDIF
        oldline=line
        line=TRIM(oldline)//' '//TRIM(ADJUSTL(newline))
        last=LEN_TRIM(line) 
        IF ( line(last:last) /= '+' ) lolin=.FALSE. 
      ENDDO

! raise counters and give the attributes of the current reaction 
! (i.e. NOT the reaction currently "in line"). 
      CALL chkkeypar (lout, numre, auxflag, xauxcf, loaux, &
                      reactkey, itype, num_n, numretype, id_n, &
                      idhv, ido2, id_m, idfo, idextra, &
                      idmeo2, idpero, idtabcf, &
                      idain, idaou, idwin, idwou, idiso, &
                      hvcf, focf, extracf, labtabcf, &
                      wincf, woucf, isocf, &
                      lostop)

! "close" the reaction before processing the "new" (in "line") reaction 
! reset flag for auxiliary info before interpreting the next reaction
      auxflag(:)=0 ; reactkey(:)=0 ;  xauxcf(:,:)=0. ; loaux=.FALSE. 
      
! "start" the new reaction, i.e interpret the reaction in "line"
      IF (loreport) WRITE(lout,'(a16,i7)') 'REACTION NUMBER ',numre+1
      IF (loreport) WRITE(lout,'(a)') TRIM(line)

      CALL readreac(lout,line,numsp,chrsp,isortlk,sortsp, &
                    numre,maxrg,maxpd,itype,arrhcf,reactkey,lostop)
      CYCLE loopreac
    ENDIF
  ENDDO loopreac

! close mechanism file
  CLOSE(lin) 

! stop if errors are found in the mechanism
  IF (lostop) THEN
    WRITE(*,*)'-- error exit -- in the mechanism'
    CLOSE(lout)
    STOP 'ERROR LEVEL 2: in the list of reaction => check *.out file '
  ENDIF

! ----------------------------------------------------------------------
! ------           FINAL CHECK OF THE CHEMICAL SCHEMES           -------  
! ----------------------------------------------------------------------

! scroll the mechanism (if lonodeep is true => simple check)
  CALL chkmeca(lout,num_n,numre,itype,numretype,maxrg,maxpd,lonodeep,lostop)

! count self reactions
  CALL countself(lout,numre,numself,lostop)

! stop if error found
  IF (lostop) THEN
    WRITE(*,*)'   -- error exit -- after reading reaction data'
    CLOSE(lin) ; CLOSE(lout)
    STOP 'ERROR LEVEL 3: in the list of reactions => check file '
  ENDIF

! get info for hv reaction (nothing returned) - for boxmodeling only
  IF (lospeech) PRINT*, 'checking hv reaction ....'
  numhv=numretype(hvidx)
  CALL hvinfo(lout,mchromo,numhv,hvcf)


! ----------------------------------------------------------------------
! --  WRITE THE MECHANISM IN A BINARY FILE TO LINK WITH THE BOXMODEL  --
! ----------------------------------------------------------------------
  numhv = numretype(1)
  numo2 = numretype(2)
  num_m = numretype(3)
  numfo = numretype(4)
  numextra = numretype(5)
  nummeo2 = numretype(6)
  numpero(1:9) = numretype(7:15)
  numtabcf = numretype(16)
  numain = numretype(17)
  numaou = numretype(18)
  numwin = numretype(19)
  numwou = numretype(20)
  numiso = numretype(21)

  ! open the binary output file
  IF (lospeech) PRINT*, 'writing the output ...'
  IF (fg_fname > 0) THEN
    !OPEN(bof,FILE=TRIM(filename)//'.bin', FORM='UNFORMATTED')
    OPEN(bof,FILE=TRIM(filename)//'.li', FORM='UNFORMATTED')
  ELSE
    !OPEN(bof,FILE='outdat.bin', FORM='UNFORMATTED')
    OPEN(bof,FILE='indat.li', FORM='UNFORMATTED')
  ENDIF

  ! write info about the size of the mechanism
  WRITE(bof) maxlsp,maxrg,maxpd,maxaux
  WRITE(bof) numsp,numre,num_n, numhv,numo2,num_m,numfo,numextra,&
               nummeo2, numtabcf,numain,numaou,numwin,numwou,numiso, &
               numself
  WRITE(bof) mxcpero

  WRITE(bof) (numpero(i),i=1,mxcpero)

  ! write species names and info
  WRITE(bof) (chrsp(i),i=1,numsp)
  WRITE(bof) (wmol(i),i=1,numsp)
  WRITE(bof) ngas, (ptrgas(i),i=1,2)
  WRITE(bof) npart, (ptrpart(i),i=1,2)
  WRITE(bof) nwall, (ptrwall(i),i=1,2)

  ! write reaction numbers (id) for each reaction type
  WRITE(bof) (id_n(i),i=1,num_n)
  WRITE(bof) (idhv(i),i=1,numhv)
  WRITE(bof) (ido2(i),i=1,numo2)
  WRITE(bof) (id_m(i),i=1,num_m)
  WRITE(bof) (idfo(i),i=1,numfo)
  WRITE(bof) (idextra(i),i=1,numextra)
  WRITE(bof) (idmeo2(i),i=1,nummeo2)
  WRITE(bof) ((idpero(i,k),i=1,numpero(k)),k=1,mxcpero)
  WRITE(bof) (idain(i),i=1,numain)
  WRITE(bof) (idaou(i),i=1,numaou)
  WRITE(bof) (idwin(i),i=1,numwin)
  WRITE(bof) (idwou(i),i=1,numwou)
  WRITE(bof) (idiso(i),i=1,numiso)
  
  ! write reaction detail - arrhenius coef.
  WRITE(bof)((arrhcf(ire,k),k=1,3),ire=1,numre)

  ! write # of reagents and products in each reaction
  WRITE(bof) (rxrgnt(ire)%nrg,rxpdct(ire)%npd, ire=1,numre)

  ! write for each rxn reagents, products and corresponding stoi. coef.
  ! NOTE: stoi. coef are written as integer type for reagents (checked 
  ! as appropriate - i.e. 1 or 2 - in chkmeca.
  WRITE(bof)(                                                               &
    (rxrgnt(ire)%idrg(i),NINT(rxrgnt(ire)%stoirg(i)), i=1,rxrgnt(ire)%nrg), &
    (rxpdct(ire)%idpd(i),rxpdct(ire)%stoipd(i), i=1,rxpdct(ire)%npd),       &
    ire=1,numre )

  ! write auxiliary information
  WRITE(bof) (NINT(hvcf(k,1)),k=1,numhv)
  WRITE(bof) (hvcf(k,2),k=1,numhv)
  WRITE(bof) ((focf(k,i),i=1,maxaux),k=1,numfo)  
  WRITE(bof) ((extracf(k,i),i=1,maxaux),k=1,numextra)
  WRITE(bof) ((woucf(k,i),i=1,1),k=1,numwou) 
  WRITE(bof) ((wincf(k,i),i=1,1),k=1,numwin)
  WRITE(bof) ((isocf(k,i),i=1,maxaux),k=1,numiso)

  ! not used but might be usefull for futur release (keep as commented)
  !  WRITE(bof) (idtabcf(i),i=1,numtabcf)    
  !  WRITE(bof) ((aoucf(k,i),i=1,2),k=1,numaou) 
  !  WRITE(bof) ((woucf(k,i),i=1,3),k=1,numwou) 
  !  WRITE(bof) ((wincf(k,i),i=1,3),k=1,numwin)

  ! done writing - close files
  CLOSE(bof)

! ----------------------------------------------------------------------
!--                        FINAL WORDS ...                            --
! ----------------------------------------------------------------------
  WRITE(lout,*)
  WRITE(lout,*)'   -- NORMAL EXIT --  LINK-FILE HAS BEEN CREATED'
  WRITE(lout,*)
  CLOSE(lout)

  PRINT*, 'number of species : ', numsp
  PRINT*, 'number of reaction : ', numre
  PRINT*, '-- no pb found in the chemical scheme --'

END PROGRAM gecko2box
 
