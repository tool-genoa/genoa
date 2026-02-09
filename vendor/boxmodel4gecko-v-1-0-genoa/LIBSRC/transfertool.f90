MODULE transfertool
IMPLICIT NONE
CONTAINS 

! ======================================================================
! Purpose: Find  the ID number of the corresponding species in the 
! gas (idgsat - not in this routine), aerosol (idasat) and wall (idwsat) 
! phases. 
! ======================================================================
SUBROUTINE readtransfer(chrsp,numsp,numain,numwin, &
                        nsat,namsat,idasat,idwsat,mech2sat)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE
      
  INTEGER,INTENT(IN)   :: numsp             ! # of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)   ! name of the species
  INTEGER,INTENT(IN)   :: numain            ! # of gas -> part. rxn
  INTEGER,INTENT(IN)   :: numwin            ! # of gas -> wall rxn
  CHARACTER(LEN=*),INTENT(IN):: namsat(:)   ! name of the partitioning species (with 'G')
  INTEGER,INTENT(IN)   :: nsat              ! # of partitioning species
  INTEGER,INTENT(OUT)  :: idasat(:)         ! ID of partitioning species - part. phase
  INTEGER,INTENT(OUT)  :: idwsat(:)         ! ID of partitioning species - wall phase
  INTEGER,INTENT(INOUT):: mech2sat(:)       ! [j] index of mech species j in the pvap table

  INTEGER              :: i,j,k
  CHARACTER(LEN=LEN(chrsp(1))) :: xspe
  LOGICAL              :: loerr
  CHARACTER(LEN=14),PARAMETER :: progname='readtransfer'
  CHARACTER(LEN=80) :: mesg1

! link species in namsat with particle phase species (idasat)
! -----------------------------------------------------------
  IF (numain > 0) THEN
    k=1
    loerr=.FALSE.

! fast (sorted) seek
    sataloop1:&
    DO i=1,nsat
      xspe=namsat(i)
      xspe(1:1)="A"      ! overwrite 'G' to 'A'
      DO j=k,numsp
        IF (xspe == chrsp(j)) THEN
          idasat(i)=j
          mech2sat(j)=i
          k=j+1          ! species order in "sat" and "mech" is the same
          CYCLE sataloop1
        ENDIF
      ENDDO

      ! if that point is reached then species not found ... 
      loerr=.TRUE. 
      EXIT sataloop1
    ENDDO sataloop1

! try slow (unsorted) seek
    IF (loerr) THEN  
      PRINT*, "   => try slow seek in readtransfer"
      sataloop2:&
      DO i=1,nsat
        xspe=namsat(i)
        xspe(1:1)="A"
        DO j=1,numsp
          IF (xspe == chrsp(j)) THEN
            idasat(i)=j
            mech2sat(j)=i
            CYCLE sataloop2
          ENDIF
        ENDDO
        
        ! species not found ... 
        mesg1="The (part. phase) species is not identified: "//TRIM(xspe)
        CALL stoperr(progname,mesg1)
      ENDDO sataloop2
    ENDIF

  ENDIF

! link species in namsat with wall phase species (idwsat)
! -----------------------------------------------------------
  IF (numwin > 0) THEN
    k=1
    satwloop:&
    DO i=1,nsat
      xspe=namsat(i)
      xspe(1:1)="W"   ! overwrite 'G' to 'A'
      DO j=k,numsp
        IF (xspe == chrsp(j)) THEN
          idwsat(i)=j
          mech2sat(j)=i
          k=j+1
          CYCLE satwloop
        ENDIF
      ENDDO

      ! if that point is reached then species not found ...
      mesg1="The (wall phase) species is not identified: "//TRIM(xspe)
      CALL stoperr(progname,mesg1)
    ENDDO satwloop
  ENDIF

END SUBROUTINE readtransfer

! ======================================================================
! Purpose: return (1) ID # of the "forward" reaction for a "reverse" 
! condensed phase -> gas phase reaction and (2) the label (row in the
! pvap file) to compute the vapor pressure of the species involved. 
! ======================================================================
SUBROUTINE rev2forward(lout,chrsp,numaou,numwin,numwou,nsat,namsat, &
                       idaou,idwin,idwou,gairefor,gwirefor)
  USE boxtool, ONLY: stoperr
  USE mechdata, ONLY: rxpdct,rxrgnt
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: lout               ! file unit for information outputs 
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)   ! name of the species
  INTEGER,INTENT(IN)  :: numaou             ! # of part. -> gas rxn
  INTEGER,INTENT(IN)  :: numwin             ! # of gas -> wall rxn
  INTEGER,INTENT(IN)  :: numwou             ! # of wall -> gas rxn
  CHARACTER(LEN=*),INTENT(IN):: namsat(:)   ! name of the partitioning species ('G')
  INTEGER,INTENT(IN)  :: nsat               ! # of partitioning species
  INTEGER, INTENT(IN) :: idaou(:)           ! mechanism rx index for ith p->g (AOU)
  INTEGER, INTENT(IN) :: idwin(:)           ! mechanism rx index for ith g->w (WIN)
  INTEGER, INTENT(IN) :: idwou(:)           ! mechanism rx index for ith w->g (WOU)
  INTEGER,INTENT(OUT) :: gairefor(:)        ! rxn ID of forward g->p matching the reverse p->g
  INTEGER,INTENT(OUT) :: gwirefor(:)        ! rxn ID of forward g->w matching the reverse w->g

  INTEGER :: g2apvaplab(SIZE(namsat))       ! species ID in the pvap list involved in the p->g 
  INTEGER :: w2gpvaplab(SIZE(namsat))       ! species ID in the pvap list involved in the w->g
  INTEGER :: g2wpvaplab(SIZE(namsat))       ! species ID in the pvap list involved in the g->w
  INTEGER :: ire, i, j, k
  INTEGER :: idaero, idgas, idwall 
  INTEGER :: ichk(SIZE(namsat))
  LOGICAL :: loerr
  CHARACTER(LEN=14),PARAMETER :: progname='rev2forward'
  CHARACTER(LEN=80) :: mesg1, mesg2
  

  gairefor(:)=0  ;  g2apvaplab(:)=0
  gwirefor(:)=0  ;  g2wpvaplab(:)=0  ;  w2gpvaplab(:)=0
  loerr=.FALSE.
  
! -----------------   
! GAS <-> AEROSOL
! -----------------   

! scroll the aerosol -> gas reactions
  k=1
  DO i=1,numaou
    ire=idaou(i)
    idaero=rxrgnt(ire)%idrg(1)
    idgas=rxpdct(ire)%idpd(1)

    ! find the forward g->p rxn matching the current reverse p->g rxn
    DO j=ire-1,ire+1,2     ! forward just before or after reverse 
      IF ( (rxrgnt(j)%idrg(1)==idgas) .AND. (rxpdct(j)%idpd(1)==idaero) ) THEN
        gairefor(i)=j      ! index of the forward p->g found
        EXIT
      ENDIF
    ENDDO

    ! find the 'Gwxyz' species in pvap list matching 'wxyz' species
    DO j=k,nsat
      IF (chrsp(idgas) == namsat(j)) THEN
        g2apvaplab(i)=j     
        k=j+1         ! species order is the same
        EXIT
      ENDIF
    ENDDO
  ENDDO 

! check that all reactions were found
  ichk=0  
  WHERE (gairefor(1:numaou) == 0) ichk(1:numaou)=1
  IF (SUM(ichk(1:numaou)) /= 0) THEN
    loerr=.TRUE.
    WRITE(lout,*) "--error--, in rev2forward, while trying to identify" 
    WRITE(lout,*) "2 successive rxn in gas/aero equilibrium. Species are:"
    DO i=1,numaou
      IF (ichk(i) /= 0 ) THEN
        WRITE(lout,*) i, ' ', chrsp(rxrgnt(idaou(i))%idrg(1))
        WRITE(lout,*) i, ' ', chrsp(rxpdct(idaou(i))%idpd(1))
      ENDIF
    ENDDO
  ENDIF

! check that all species were found
  ichk=0  
  WHERE (g2apvaplab(1:numaou) == 0) ichk(1:numaou)=1
  IF (SUM(ichk(1:numaou)) /= 0) THEN
    loerr=.TRUE.
    WRITE(lout,*) "--error--, vapor pressure not found for species:"
    DO i=1,numaou
      IF (ichk(i) /= 0 ) THEN
        WRITE(lout,*) i, ' ', chrsp(rxpdct(idaou(i))%idpd(1))
      ENDIF
    ENDDO
  ENDIF

! -----------------   
! GAS -> WALL
! -----------------   

! scroll the gas -> wall reactions
  k=1
  DO i=1,numwin
    ire=idwin(i)
    idgas=rxrgnt(ire)%idrg(1)
    idwall=rxpdct(ire)%idpd(1)
    ! find index of the 'Gwxyz' species in the pvap list 
    DO j=k,nsat
      IF (chrsp(idgas) == namsat(j)) THEN
        g2wpvaplab(i)=j
        k=j+1         ! species order is the same
        EXIT
      ENDIF
    ENDDO
  ENDDO 

! check that all species were found
  ichk(:)=0  
  WHERE (g2wpvaplab(1:numwin) == 0) ichk(1:numwin)=1
  IF (SUM(ichk(1:numwin)) /= 0) THEN
    loerr=.TRUE.
    WRITE(lout,*) "--error--, vapor pressure not found for species:"
    DO i=1,numwin
      IF (ichk(i) /= 0 ) THEN
        WRITE(lout,*) i, ' ', chrsp(rxpdct(idwin(i))%idpd(1))
      ENDIF
    ENDDO
  ENDIF

! -----------------   
! WALL -> GAS
! -----------------   

! scroll the wall -> gas reactions
  k=1
  DO i=1,numwou
    ire=idwou(i)
    idwall=rxrgnt(ire)%idrg(1)
    idgas=rxpdct(ire)%idpd(1)

    ! find the forward g->w rxn matching the current reverse w->g rxn
    DO j=ire-1,ire+1,2   ! forward just before or after reverse
      IF ( (rxrgnt(j)%idrg(1)==idgas) .AND. (rxpdct(j)%idpd(1)==idwall) ) THEN
        gwirefor(i)=j    ! index of the forward p->g found
        EXIT
      ENDIF
    ENDDO

    ! find the 'Gwxyz' species in pvap list matching 'wxyz' species
    DO j=k,nsat
      IF (chrsp(idgas) == namsat(j)) THEN
        w2gpvaplab(i)=j
        k=j+1         ! species order is the same
        EXIT
      ENDIF
    ENDDO
  ENDDO 

! check that all reactions were found
  ichk=0  
  WHERE (gwirefor(1:numwou) == 0) ichk(1:numwou)=1
  IF (SUM(ichk(1:numwou)) /= 0) THEN
    loerr=.TRUE.
    WRITE(lout,*) "--error--, in rev2forwardnot while trying to identify" 
    WRITE(lout,*) "  2 successive reactions in gas/wall equilibrium. Species are:"
    DO i=1,numwou
      IF (ichk(i) /= 0 ) THEN
        WRITE(lout,*) i, ' ', chrsp(rxrgnt(idwou(i))%idrg(1))
        WRITE(lout,*) i, ' ', chrsp(rxpdct(idwou(i))%idpd(1))
      ENDIF
    ENDDO
  ENDIF

! check that all species were found
  ichk=0  
  WHERE (w2gpvaplab(1:numwou) == 0) ichk(1:numwou)=1
  IF (SUM(ichk(1:numwou)) /= 0) THEN
    loerr=.TRUE.
    WRITE(lout,*) "--error--, vapor pressure not found for species:"
    DO i=1,numwou
      IF (ichk(i) /= 0 ) THEN
        WRITE(lout,*) i, ' ', chrsp(rxpdct(idwou(i))%idpd(1))
      ENDIF
    ENDDO
  ENDIF

! -----------------   
! STOP if error 
! -----------------   
  IF (loerr) THEN 
    mesg1="reverse reactions and/or pvap label not found"
    mesg2="See the *.out file for the details" 
    CALL stoperr(progname,mesg1,mesg2)
  ENDIF

END SUBROUTINE rev2forward


! ======================================================================
! Purpose: return the particle number concentration (pnc) for the 
! non-volatile organic species (i.e. the seed species here) provided 
! as input. The properties of the seed (M, density, radius of the 
! particles) are provided by the parameter module.
! ======================================================================
SUBROUTINE get1stpnc(lout,pnc)
  USE parameter_mod, ONLY: Mseed,avogadro,pi,densseed,Rp0,nvoc,nvic,fg_output
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: lout      ! file unit for information outputs
  REAL,INTENT(OUT)   :: pnc       ! particle number concentration

  CHARACTER(LEN=12),PARAMETER :: progname='get1stpnc'
  CHARACTER(LEN=80) :: mesg1

  ! compute particle number concentration given seed (nvoc) 
  ! concentration(molec/cm3) and initial radius   
  IF (nvoc > 100.) THEN     ! arbitrarily low number here 
    pnc=(Mseed*(nvoc+nvic))/(avogadro*(4./3.)*pi*densseed*(Rp0**3)) 
  ELSE
    mesg1="Initial particle concentration requested without seed (nvoc) prescribed"
    CALL stoperr(progname,mesg1)
  ENDIF
 
  IF (fg_output>0) THEN
    WRITE (lout,'(a24,1P,E10.2,0P)') 'Particle number (/cm3):', pnc
    WRITE (lout,*) '---------------'
  ENDIF
END SUBROUTINE get1stpnc


! ======================================================================
! Purpose: compute particle loss due to dilution or any other loss 
! process (e.g. to wall). The subroutine should only be used if those
! processes are not included the set of ODE (see iteration). This is the
! case when equilibrium it forced at each time step (i.e. if soa_fg=1).
! Concentration at t+deltat is computed assuming a simple 1st order
! loss. 
! ======================================================================
SUBROUTINE particleloss(soa_fg,kdilu,kaerloss,nsat,deltat,caer)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: soa_fg      ! SOA flag (method used for gas/particle partioning)
  REAL,INTENT(IN)    :: kdilu       ! first order dilution rate constant 
  REAL,INTENT(IN)    :: kaerloss    ! first order loss rate constant for particles
  INTEGER,INTENT(IN) :: nsat        ! number of species in the condensed phase
  REAL,INTENT(IN)    :: deltat      ! time step
  REAL,INTENT(INOUT) :: caer(:)     ! species concentration in the condensed phase
  
  REAL :: kloss, scalefac
  CHARACTER(LEN=12),PARAMETER :: progname='particleloss'
  CHARACTER(LEN=80) :: mesg1
  
  IF (soa_fg/=1) THEN
    mesg1="Unexpected flag number used (must be 1 in the routine)."
    CALL stoperr(progname,mesg1)
  ENDIF
  
  kloss=kdilu+kaerloss
  IF (kloss > 1E-30) THEN
    scalefac=EXP(-kloss*deltat)
    caer(1:nsat)=caer(1:nsat)*scalefac
  ENDIF
END SUBROUTINE particleloss

END MODULE transfertool
