MODULE mecatool
IMPLICIT NONE
CONTAINS
!***********************************************************************
! PURPOSE: Scroll the mechanism and check for unexpected reactions, in 
! particular duplicate reactions. Check the stoichiometric on the 
! reagent side ("integer" 1 or 2). 
! Warning : the mechanism check can be time consuming and is usually 
! not needed. If "lonodeep" is turned to false, only a simple 
! reaction counting is performed (useful to optimize the size of the 
! tables in the boxmodel).
! Nothing is returned by this routine, but the "lostop" flag is raised
! if an error is found 
!***********************************************************************
SUBROUTINE chkmeca(lout,num_n,numre,itype,numretype,maxrg,maxpd,lonodeep,lostop)
  USE keywordlist, ONLY: mxkeywd,keywdlist
  USE reactiontool, ONLY: rxrgnt,rxpdct
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: lout         ! unit output file
  INTEGER,INTENT(IN)  :: num_n        ! # of "normal" reaction (i.e. not using a keyword)
  INTEGER,INTENT(IN)  :: numre        ! # of reactions
  INTEGER,INTENT(IN)  :: itype(:)     ! "type" of the reactions  
  INTEGER,INTENT(IN)  :: numretype(:) ! # of reaction per keyword type
  INTEGER,INTENT(IN)  :: maxrg        ! max # of reagent identified in a reaction
  INTEGER,INTENT(IN)  :: maxpd        ! max # of product identified in a reaction
  LOGICAL,INTENT(IN)  :: lonodeep     ! flag (deeper search if true)
  LOGICAL,INTENT(OUT) :: lostop       ! stop error flag

  INTEGER i, j, l, nchk

! -------------------------------------------
! CHECK THE NUMBER OF REAGENTS IN A REACTION 
! -------------------------------------------
  IF (maxrg>2) THEN
    WRITE(lout,*) '--error--, reaction with reagents>2 is unexpected'
    lostop=.TRUE.
  ENDIF
  WRITE(lout,*) 
  WRITE(lout,*) 'max # of pdct found in a rxn:',maxpd ; WRITE(lout,*)

! -------------------------------------------
! CHECK THAT REAGENT STOI COEF ARE 1 OR 2  
! -------------------------------------------
  DO i=1,numre
    rgloop: DO j=1,rxrgnt(i)%nrg
      IF ( (abs(rxrgnt(i)%stoirg(j)) - 1.0) < 1e-4) CYCLE rgloop
      IF ( (abs(rxrgnt(i)%stoirg(j)) - 2.0) < 1e-4) CYCLE rgloop
      WRITE(lout,*) '--error--, reaction with unexpected stoi coef, reagent side'
      WRITE(lout,*) '           stoi.coef is:',rxrgnt(i)%stoirg(j)
      lostop=.TRUE.
    ENDDO rgloop
  ENDDO

! ---------------------------
! CHECK THE SUM OF REACTIONS 
! ---------------------------

! check that no reaction was forgotten ...
  WRITE(lout,*)
  WRITE(lout,*) 'num_n=',num_n
  DO i=1,mxkeywd 
    WRITE(lout,*) keywdlist(i),'=',numretype(i)
  ENDDO
  WRITE(lout,*) '------------------------------------------------'
  WRITE(lout,*) 'numre=',numre ; WRITE(lout,*)

  nchk=SUM(numretype)+num_n
  IF (nchk /= numre) THEN
    WRITE(lout,*) '--error--  in the equation count'
    lostop=.TRUE.
  ENDIF

! exit if no deep check is wanted (usually not needed)
  IF (lonodeep) THEN
    WRITE(lout,*) 'THE SET OF REACTIONS WAS NOT DEEPLY CHECKED'
    RETURN
  ENDIF


! -----------------------------------------------
! IDENTIFY NO BUDGET REACTIONS (REACTANTS = PRODUCTS)
! -----------------------------------------------
  nobud: DO i=1,numre
    IF (rxrgnt(i)%nrg == rxpdct(i)%npd) THEN
      DO j=1,rxrgnt(i)%nrg
        IF (rxrgnt(i)%idrg(j) /= rxpdct(i)%idpd(j)) CYCLE nobud
      ENDDO 
    ELSE
      CYCLE
    ENDIF

    WRITE(lout,*)'--note-- equation number ',i,' is trivial'
    WRITE(lout,*)
    lostop=.TRUE.
  ENDDO nobud 

! ---------------------------------------------------------------
! IDENTIFY DUPLICATE REACTIONS (DOES NOT PRODUCE A STOP ERROR)
! --------------------------------------------------------------
  nodup : DO i=1,numre-1
   chkloop : DO j=i+1,numre

      IF ( rxrgnt(i)%nrg == rxrgnt(j)%nrg .AND.   &
           rxpdct(i)%npd == rxpdct(j)%npd ) THEN
        DO l=1,rxrgnt(i)%nrg
          IF (rxrgnt(i)%idrg(l) /= rxrgnt(j)%idrg(l)) CYCLE chkloop
        ENDDO
        DO l=1,rxpdct(i)%npd
          IF (rxpdct(i)%idpd(l) /= rxpdct(j)%idpd(l)) CYCLE chkloop
        ENDDO

        IF (itype(i) == itype(j)) THEN
          WRITE(lout,*)'--note-- equations number ',i, ' and ',j, 'are identical'
          PRINT*,      '--note-- equations number ',i, ' and ',j, 'are identical'
          lostop=.TRUE.
        ENDIF
      ENDIF

    ENDDO chkloop
  ENDDO nodup

END SUBROUTINE chkmeca


!***********************************************************************
! PURPOSE: provide information about hv reactions. 
! No output for this subroutine - just numbers that may be useful to
! optimize the size of the tables once used in the boxmodel.
! Indeed, to decrease CPU time and limit the size of the tables,
! HV reactions are categorized (in the boxmodel) according to the 
! occurence of each label. This subroutine provides this 
! information (in the *.out file) to then optimize the size of 
! the tables in the boxmodel. 
!***********************************************************************
SUBROUTINE hvinfo(lout,mchromo,numhv,hvcf)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lout        ! unit for the output file
  INTEGER, INTENT(IN) :: mchromo     ! maximum number of chromophore (label)
  INTEGER, INTENT(IN) :: numhv       ! number of "hv" reactions
  REAL, INTENT(IN)    :: hvcf(:,:)   ! hv auxiliary info

  INTEGER :: ntchromo, mxlab, chromocf(mchromo), numchromo(mchromo)
  INTEGER :: i,j,idcf

! Initialize
! ------------
  ntchromo=0 ; mxlab=0 ; chromocf(:)=0 ; numchromo(:)=0


! scroll the hv reactions - identify # of distinct hv label (i.e. chromo)
! -----------------------------------------------------------------
  looprehv : DO i=1,numhv
    idcf=NINT(hvcf(i,1))

! check if hv label already exist. If yes, store new # and goto next reaction
    DO j=1,ntchromo
      IF (idcf == chromocf(j)) THEN
        numchromo(j)=numchromo(j)+1
        CYCLE looprehv
      ENDIF
    ENDDO

! if that point is reached, hv label does not exist => add to the list
    ntchromo=ntchromo+1
    IF (ntchromo > mchromo) THEN
      WRITE(lout,*) '--error--, number of chromophore (i.e. labels)'
      WRITE(lout,*) '           exceed mchromo. Change akparameter.h'
      STOP '--error-- in hvinfo subroutine. ntchromo > mchromo' 
    ENDIF
    chromocf(ntchromo)=idcf
    numchromo(ntchromo)=1        
  ENDDO looprehv

! write info about hv reactions 
! ------------------------------
  mxlab=MAXVAL(numchromo)
  WRITE(lout,*) ' '
  WRITE(lout,*) 'total number of HV reaction :',numhv
  WRITE(lout,*) 'number of chromophore (label):',ntchromo
  WRITE(lout,*) 'max number of species in a given label:',mxlab
  WRITE(lout,*) ' '
  WRITE(lout,*) 'occurence (column 3) of each label (column 2):'
  DO i=1,ntchromo
    WRITE (lout,'(i3,2x,i6,2x,i6)') i,chromocf(i), numchromo(i)
  ENDDO

END SUBROUTINE hvinfo

! ======================================================================
! PURPOSE : count the number of self reaction and check the stoi. coef.
! on the reactant side. 
! ======================================================================
SUBROUTINE countself(lout,numre,numself,lostop)
  USE reactiontool, ONLY: rxrgnt
  IMPLICIT NONE
  INTEGER,INTENT(IN)  :: lout          ! unit output file
  INTEGER,INTENT(IN)  :: numre         ! # of reactions
  INTEGER,INTENT(OUT) :: numself       ! # of self reactions
  LOGICAL,INTENT(OUT) :: lostop        ! stop error flag

  INTEGER :: ire,loc, j
  
  numself=0

! For the reactant side: the default stoi. coef. is 1. If equal 2 (self
! reaction), then store the info in a table.  Stoi. coef. with values 
! not equal to 1. or 2. are not allowed in this version of the program. 
! Purpose is to win time in the routine involved in time integration.
  numself=0
  reactloop: DO ire=1,numre
    loc=0
    DO j=1,rxrgnt(ire)%nrg

      IF (rxrgnt(ire)%stoirg(j) == 1.) THEN 
        CYCLE

      ELSE IF (rxrgnt(ire)%stoirg(j) == 2.) THEN
        numself=numself+1
        loc=loc+1
        IF (loc > 1) THEN   ! only 1 "self reaction" per reaction
          WRITE(lout,*) '--error--, in countself. Two species have a value '
          WRITE(lout,*) '           of 2 for the stoi. coef. in the same'
          WRITE(lout,*) '           reaction (i.e. 4 order reaction).' 
          WRITE(lout,*) '           This case is not allowed.'
          lostop=.TRUE.
        ENDIF
!        idselfreac(numself,1)=ire
!        idselfreac(numself,2)=j

      ELSE                    ! error, stoi. coef. not allowed
        WRITE(lout,*) '--error--, in readmech. Stoi. coef. at the reactant'
        WRITE(lout,*) 'side is not 1 or 2 in the reaction #:',ire 
        lostop=.TRUE.
      ENDIF
    ENDDO
  ENDDO reactloop
END SUBROUTINE countself


END MODULE mecatool
