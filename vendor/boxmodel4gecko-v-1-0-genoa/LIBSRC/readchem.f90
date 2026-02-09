MODULE readchem
IMPLICIT NONE
CONTAINS

!=======================================================================  
! Purpose : 
! 1. Read the binary file produced by the interpreter. 
! 2. Put away some reactions (to decrease memory storage and to increase
!    efficiency in solving the system). A special storage is performed
!    for "self reactions (i.e. reaction X+X => products)
! 3. Count the number of distinct HV label (i.e. "njlab" to be stored 
!    in phototool). Needed here as output to allocate next the memory 
!    required to compute the J values.
!=======================================================================  
SUBROUTINE readmech(lout,numreacro2, &
                    nspg,nspp,nspw,ptrgas,ptrpart,ptrwall, &
                    idhv,ido2,id_m,idfo,idextra,idmeo2,idreacro2,&
                    idain,idaou,idwin,idwou,idiso,&
                    arrhcf,rxhvtag,hvfact,focf,extracf,&
                    woucf,wincf,isocf,&
                    idselfreac,wmol,chrsp,njlab)
  USE boxtool, ONLY: stoperr
  USE mechdata, ONLY: rxpdct,rxrgnt
  USE parameter_mod, ONLY: input_dir_chem, fg_output
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: lout              ! file unit for information outputs   
  INTEGER, INTENT(OUT) :: numreacro2(:)     ! # of rx with "PEROj" in class [j]
  INTEGER, INTENT(OUT) :: nspg              ! # of gas phase species
  INTEGER, INTENT(OUT) :: nspp              ! # of part. phase species
  INTEGER, INTENT(OUT) :: nspw              ! # of wall phase species
  INTEGER, INTENT(OUT) :: ptrgas(:)         ! pointer (1st and last) gas phase species
  INTEGER, INTENT(OUT) :: ptrpart(:)        ! pointer (1st and last) part. phase species
  INTEGER, INTENT(OUT) :: ptrwall(:)        ! pointer (1st and last) wall phase species
  INTEGER, INTENT(OUT) :: idhv(:)           ! mechanism rx index of the ith 'HV' reaction 
  INTEGER, INTENT(OUT) :: ido2(:)           ! mechanism rx index of the ith 'OXYGEN' reaction
  INTEGER, INTENT(OUT) :: id_m(:)           ! mechanism rx index of the ith 'TBODY' reaction
  INTEGER, INTENT(OUT) :: idfo(:)           ! mechanism rx index of the ith 'FALLOFF' reaction
  INTEGER, INTENT(OUT) :: idextra(:)        ! mechanism rx index of the ith 'EXTRA' reaction
  INTEGER, INTENT(OUT) :: idiso(:)          ! mechanism rx index of the ith 'ISOM' reaction
  INTEGER, INTENT(OUT) :: idmeo2(:)         ! mechanism rx index of the ith 'MEPERO' reaction
  INTEGER, INTENT(OUT) :: idreacro2(:,:)    ! [i,j] mechanism rx index of the ith 'PEROj' reaction
  INTEGER, INTENT(OUT) :: idain(:)          ! mechanism rx index for ith gas -> part. "reaction"
  INTEGER, INTENT(OUT) :: idaou(:)          ! mechanism rx index for ith part. -> gas "reaction"
  INTEGER, INTENT(OUT) :: idwin(:)          ! mechanism rx index for ith gas -> wall "reaction"
  INTEGER, INTENT(OUT) :: idwou(:)          ! mechanism rx index for ith wall -> gas "reaction"
  REAL, INTENT(OUT)    :: arrhcf(:,:)       ! [i,3] arrhenius parameter for the ith mech. rx index
  INTEGER, INTENT(OUT) :: rxhvtag(:)        ! reference (tag) number for J(xs,qy) of ith hv reaction 
  REAL, INTENT(OUT)    :: hvfact(:)         ! scaling factor applied to reference J 
  REAL, INTENT(OUT)    :: focf(:,:)         ! [i,j] j coef. (k0, Fc...) for the ith falloff rx  
  REAL, INTENT(OUT)    :: extracf(:,:)      ! [i,j] j coef. (k...) for the ith "extra" rx
  REAL, INTENT(OUT)    :: woucf(:,:)        ! [i,j] j extra coef. for the ith w->g rx (WOU)
  REAL, INTENT(OUT)    :: wincf(:,:)        ! [i,j] j extra coef. for the ith g->w rx (WIN)
  REAL, INTENT(OUT)    :: isocf(:,:)        ! [i,j] j extra coef. for the ith 'ISOM' rx
  INTEGER, INTENT(OUT) :: idselfreac(:,:)   ! [i,1] mech. rx index and species ID [i,2] for self 
  CHARACTER(LEN=*),INTENT(OUT) :: chrsp(:)  ! list (names) of the species
  REAL, INTENT(OUT)    :: wmol(:)           ! molar mass of the ith species 
  INTEGER, INTENT(OUT) :: njlab             ! # of distinct jlabel  

! internal
  LOGICAL :: lostop
  INTEGER :: i,j,k,ire,filin,loc,label,nreagent,nproduct
  INTEGER,ALLOCATABLE :: stoirg(:,:)        ! [i,j] stoi. coef. for species j in rxn i

  INTEGER :: numsp,numre,num_n,numhv,numo2,num_m,numfo,numextra,nummeo2
  INTEGER :: numain,numaou,numwin,numwou,numiso,ncpero,numself
  INTEGER :: dummyid,dummyid2,mx1spe,mx2spe,mxaux,mxlsp,nself

! J table memory allocated later, need a temporary max # of J labels
  INTEGER, PARAMETER  :: mxjlabel=250 
  INTEGER :: hvlabel(mxjlabel)

  CHARACTER(LEN=10),PARAMETER :: progname='readmech'
  CHARACTER(LEN=80) :: mesg1, mesg2

! ----------------------------- 
! initialization of all values
! ----------------------------- 
  numsp=0 ; chrsp(:)=' '
  numre=0 ; num_n=0 ; num_m=0 ; numfo=0 ; numhv=0 
  numextra=0 ; numo2=0 ; numiso=0 ; nummeo2=0 ; mx1spe=0 ; mx2spe=0
      
  id_m(:)=0 ; idfo(:)=0 ; idhv(:)=0 ; idextra(:)=0 ; ido2(:)=0 
  idmeo2(:)=0 ; ptrgas(:)=0 ; ptrpart(:)=0 ; ptrwall(:)=0 

  arrhcf(:,:)=0. ; focf(:,:)=0. ; rxhvtag(:)=0 ; hvfact(:)=0. 
  extracf(:,:)=0. ; isocf(:,:)=0. ; wmol(:)=0.

  lostop=.FALSE.

! ------------------------------------------
! Start reading the chemical scheme (binary file)
! ------------------------------------------

! open the chemical file
  filin=12
  i=LEN_TRIM(input_dir_chem)
  IF (input_dir_chem(i:i)=='/') THEN
    OPEN (filin, FILE=TRIM(input_dir_chem)//'indat.li',FORM='unformatted',STATUS='OLD')
  ELSE
    OPEN (filin, FILE=TRIM(input_dir_chem)//'.li',FORM='unformatted',STATUS='OLD')
  ENDIF
  
  ! following numbers were already read (see readmechsize) to allocate memory 
  READ(filin) mxlsp,mx1spe,mx2spe,mxaux
  READ(filin) numsp,numre,num_n,numhv,numo2,num_m,numfo,numextra,&
              nummeo2,dummyid,numain,numaou,numwin,numwou,numiso,&
              numself
  READ(filin) ncpero

  ! not read before
  READ(filin) (numreacro2(i),i=1,ncpero)

  ALLOCATE(stoirg(numre,mx1spe)) ; stoirg(:,:) = 0

! -----------------------------------------------
! Check size allocate for the various table 
! -----------------------------------------------

! double check for the length of the tables
  IF (mx1spe > 2) THEN
    mesg1="Unexpected number of reagents in a reaction (larger than 2)."
    CALL stoperr(progname,mesg1)
  ENDIF
  IF (numre > SIZE(rxrgnt) ) lostop=.TRUE.                       ! # of reactants per rxn 
  IF (mx2spe > SIZE(rxpdct) ) lostop=.TRUE.                      ! # of products per rxn
  IF (mxaux > SIZE(extracf,2) ) lostop=.TRUE.                    ! # of auxiliary info
  IF (mxlsp > LEN(chrsp(1))) lostop=.TRUE.                       ! length of species names
  IF (numsp > SIZE(chrsp,1)) lostop=.TRUE.                       ! # of species
  IF (numre > SIZE(arrhcf,1)) lostop=.TRUE.                      ! # of rxn
  IF (numhv > SIZE(idhv,1)) lostop=.TRUE.                        ! # of hv rxn
  IF (numo2 > SIZE(ido2,1)) lostop=.TRUE.                        ! # of O2 rxn
  IF (num_m > SIZE(id_m,1)) lostop=.TRUE.                        ! # third body rxn
  IF (numfo > SIZE(idfo,1)) lostop=.TRUE.                        ! # of falloff rxn
  IF (numextra > SIZE(idextra,1)) lostop=.TRUE.                  ! # of extra rxn
  IF (nummeo2 > SIZE(idmeo2,1)) lostop=.TRUE.                    ! # of ch3o2 rxn
  IF (numain > SIZE(idain,1)) lostop=.TRUE.                      ! # of gas -> aero rxn
  IF (numaou > SIZE(idaou,1)) lostop=.TRUE.                      ! # of aero -> gas rxn
  IF (numwin > SIZE(idwin,1)) lostop=.TRUE.                      ! # of gas -> wall rxn
  IF (numwou > SIZE(idwou,1)) lostop=.TRUE.                      ! # of wall -> gas rxn
  IF (numiso > SIZE(idiso,1) ) lostop=.TRUE.                     ! # of isomerisation rxn
  IF (ncpero > SIZE(numreacro2,1)) lostop=.TRUE.                 ! # of RO2 classes
  IF (MAXVAL(numreacro2) > SIZE(idreacro2,1) ) lostop=.TRUE.     ! # of RO2+PERO rxn

! stop on error
  IF (lostop) THEN
    mesg1="Sizes allocated to the tables does not match the mechanism data"
    CALL stoperr(progname,mesg1)
  ENDIF

! ---------------------------------------------------
! Read the mechanism 
! ---------------------------------------------------

! read species names and info
  READ(filin) (chrsp(i),i=1,numsp)       ! all species (all phases)
  READ(filin) (wmol(i),i=1,numsp)        ! molar masses 
  READ(filin) nspg, (ptrgas(i),i=1,2)    ! gas phase species
  READ(filin) nspp, (ptrpart(i),i=1,2)   ! particle phase species
  READ(filin) nspw, (ptrwall(i),i=1,2)   ! wall phase species

! read reaction numbers (id) for each reaction type
  READ(filin) (dummyid,i=1,num_n)
  READ(filin) (idhv(i),i=1,numhv)
  READ(filin) (ido2(i),i=1,numo2)
  READ(filin) (id_m(i),i=1,num_m)
  READ(filin) (idfo(i),i=1,numfo)
  READ(filin) (idextra(i),i=1,numextra)
  READ(filin) (idmeo2(i),i=1,nummeo2)
  READ(filin) ((idreacro2(i,k),i=1,numreacro2(k)),k=1,ncpero)
  READ(filin) (idain(i),i=1,numain)
  READ(filin) (idaou(i),i=1,numaou)
  READ(filin) (idwin(i),i=1,numwin)
  READ(filin) (idwou(i),i=1,numwou)
  READ(filin) (idiso(i),i=1,numiso)

! read reaction detail (species, stoi. coef. and arrhenius coef.)
  READ(filin) ((arrhcf(ire,k),k=1,3),ire=1,numre)
  
  ! read # of reagents and products in each rxn and allocate memory
  READ(filin) (rxrgnt(ire)%nrg, rxpdct(ire)%npd, ire=1,numre)
  DO ire=1,numre
    nreagent = rxrgnt(ire)%nrg 
    ALLOCATE(rxrgnt(ire)%idrg(nreagent))
    nproduct = rxpdct(ire)%npd 
    ALLOCATE(rxpdct(ire)%idpd(nproduct), rxpdct(ire)%stoipd(nproduct) ) 
  ENDDO
  
  ! read for each rxn reagents, products and corresponding stoi. coef.
  ! note: stoi. coef are integer type for reagent, real type for products
  READ(filin) (                                                        &
    (rxrgnt(ire)%idrg(i), stoirg(ire,i), i=1,rxrgnt(ire)%nrg),         & 
    (rxpdct(ire)%idpd(i), rxpdct(ire)%stoipd(i), i=1,rxpdct(ire)%npd), &
    ire=1,numre)

  ! read the auxiliary information (hv, falloff ...)
  READ(filin) (rxhvtag(k),k=1,numhv)
  READ(filin) (hvfact(k),k=1,numhv)
  READ(filin) ((focf(k,i),i=1,mxaux),k=1,numfo) 
  READ(filin) ((extracf(k,i),i=1,mxaux),k=1,numextra)
  READ(filin) ((woucf(k,i),i=1,1),k=1,numwou) 
  READ(filin) ((wincf(k,i),i=1,1),k=1,numwin) 
  READ(filin) ((isocf(k,i),i=1,mxaux),k=1,numiso)

! nothing left to read. Close the file.      
  CLOSE(filin)

! -----------------------------------------------------
! Check stoichiometric coefficient on the reagent side
! -----------------------------------------------------

! The default stoi. coef. is 1. If equal 2 (self reaction), then store 
! the info in a table.  Stoi. coef. with values not equal to 1 or 2 are 
! not allowed in this version of the program.
  nself=0
  reactloop: DO ire=1,numre
    loc=0
    DO j=1,rxrgnt(ire)%nrg
      IF (stoirg(ire,j) == 1) THEN 
        CYCLE
      ELSE IF (stoirg(ire,j) == 2) THEN
        nself=nself+1
        IF (nself > SIZE(idselfreac,1)) THEN
          mesg1="The number of self reactions exceed the size allowed."
          CALL stoperr(progname,mesg1)
        ENDIF
        loc=loc+1
        IF (loc > 1) THEN   ! 2A + 2B => product rxn (not allowed)
          mesg1="Two species has a stoi. coef. of 2 in the same rxn."
          CALL stoperr(progname,mesg1)
        ENDIF
        idselfreac(nself,1)=ire              ! store rx number
        idselfreac(nself,2)=rxrgnt(ire)%idrg(j)  ! store species ID
      ELSE                    ! error, stoi. coef. not allowed
        mesg1="Stoi. coef. at the reactant side is not 1 or 2 " 
        WRITE(mesg2,'(a,i7)') "in reaction number: ",ire 
        CALL stoperr(progname,mesg1,mesg2)
      ENDIF
    ENDDO
  ENDDO reactloop
  IF (nself/=numself) THEN
    mesg1="expect nself/=numself"
    CALL stoperr(progname,mesg1)
  ENDIF
  DEALLOCATE(stoirg)

! ----------------------------------------
! get the # of distinct 'HV' label (njlab) 
! ----------------------------------------

  njlab=0                  ! # of distinct hv label (chromophore) used 
  hvlabel(:)=0             ! "temporary" table (needed to compute njlab)

  ! loop over HV rxn, count each new label 
  hvreact: DO i=1,numhv
    label=rxhvtag(i)

    ! check if label already recorded. If yes, goto next
    DO j=1,njlab
      IF (label == hvlabel(j))  CYCLE hvreact 
    ENDDO

    ! if that point is reached, the label is not recorded yet 
    njlab=njlab+1
    IF (njlab > mxjlabel) THEN
      mesg1="number of J labels > mxjlabel"
       CALL stoperr(progname,mesg1)
    ENDIF
    hvlabel(njlab)=label    ! record the new label
  ENDDO hvreact              

! -------------------------
! Write some data 
! -------------------------
  IF (fg_output>0) THEN
    WRITE(lout,*) '------ INFO ABOUT THE MECHANISM USED -------'
    WRITE(lout,*) '  numsp       =',numsp
    WRITE(lout,*) '  numre       =',numre
    WRITE(lout,*) '  num_n       =',num_n
    WRITE(lout,*) '  numhv       =',numhv
    WRITE(lout,*) '  numo2       =',numo2
    WRITE(lout,*) '  num_m       =',num_m
    WRITE(lout,*) '  numfo       =',numfo
    WRITE(lout,*) '  numextra    =',numextra
    WRITE(lout,*) '  nummeo2     =',nummeo2
    WRITE(lout,*) '  numain      =',numain
    WRITE(lout,*) '  numaou      =',numaou
    WRITE(lout,*) '  numwin      =',numwin
    WRITE(lout,*) '  numwou      =',numwou
    WRITE(lout,*) '  numiso      =',numiso
    WRITE(lout,*) '  numreacro2  =',numreacro2
    WRITE(lout,*) '  njlab       =',njlab
    WRITE(lout,*) ''

  ENDIF
END SUBROUTINE readmech

END MODULE readchem
