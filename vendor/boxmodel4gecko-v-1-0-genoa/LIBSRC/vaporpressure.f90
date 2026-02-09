MODULE vaporpressure
!$ USE omp_lib
IMPLICIT NONE
LOGICAL :: lo1pvap
CONTAINS 

!=======================================================================
! PURPOSE: Read the saturation vapor pressure and the latent heat of 
! vaporization of the species undergoing gas/particle partitioning in
! the mechanism. Properties are provided for a reference temperature,
! which is 1st value to be read. The format of the file to be read is:
!   comment lines (as many as wanted)
!   Tref   data
!   name(1) pvap(1) eheat(1)
!   name(i) pvap(i) eheat(i)
!   ...
!   END
!=======================================================================
SUBROUTINE readpvap(chrsp,numsp,nsat,namsat,idgsat,rt4psat,rpsat, &
                    eheat,mech2sat,lodicnam)
  USE boxtool, ONLY: stoperr
  USE parameter_mod, ONLY: Rgas,input_dir_chem
  IMPLICIT NONE
  INTEGER,INTENT(IN)   :: numsp               ! # of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)     ! name of the species
  INTEGER,INTENT(OUT)  :: nsat                ! # of partitioning species 
  CHARACTER(LEN=*),INTENT(OUT)::  namsat(:)   ! name of the partitioning species
  INTEGER,INTENT(OUT)  :: idgsat(:)           ! ID of the partitioning species
  REAL,INTENT(OUT)     :: rt4psat             ! reference T for vapor pressure
  REAL,INTENT(OUT)     :: rpsat(:)            ! vapor pressure @ ref T
  REAL,INTENT(OUT)     :: eheat(:)            ! latent heat of evaporation (@ ref T in K-1)
  INTEGER,INTENT(INOUT):: mech2sat(:)         ! [j] index of mech species j in the pvap table
  LOGICAL,INTENT(IN),OPTIONAL :: lodicnam     ! if true and present then chrsp is dictionary names

  CHARACTER(LEN=160) :: line
  INTEGER            :: i, j, iloc, ios, filu
  LOGICAL            :: loerr
  INTEGER            :: ibeg 

  CHARACTER(LEN=10),PARAMETER :: progname='readpvap'
  CHARACTER(LEN=80) :: mesg1, mesg2
  
! chrsp is either a "gasphase name" with 1st letter G (boxmodel) or the
! dictionary name (no "G" for postprocessing). Default: gasphase names
  ibeg=1 
  IF (PRESENT(lodicnam)) THEN
    IF (lodicnam) THEN ; ibeg=2 ; ELSE; ibeg=1 ; ENDIF
  ENDIF

! open the file
! -------------
  nsat=0 ; idgsat(:)=0
  filu=12
  i=LEN_TRIM(input_dir_chem)
  IF (input_dir_chem(i:i)=='/') THEN
    OPEN (filu, FILE=TRIM(input_dir_chem)//'pvap.sat', STATUS='OLD')
  ELSE
    OPEN (filu, FILE=TRIM(input_dir_chem)//'.pvap', STATUS='OLD')
  ENDIF

! read the file
! -------------

! read comment and exit once Tref found
  DO                             
    READ (filu, '(a)', IOSTAT=ios) line
    IF (ios/=0) THEN 
      mesg1="ios issue while reading file pvap.sat"
      mesg2="keyword 'TREF' missing ?" 
      CALL stoperr(progname,mesg1,mesg2)
    ENDIF
    IF (line(1:4) /= 'TREF') CYCLE   ! get next line 
    READ (line(5:),*) rt4psat        ! keywd TREF found
    EXIT 
  ENDDO

! read pvap in the file
  nsat=0
  DO                             
    READ (filu, '(a)', IOSTAT=ios) line
 
    ! exit if keyword END found
    IF (ios/=0) THEN 
      mesg1="ios issue while reading file pvap.sat"
      mesg2="keyword 'END' missing ?" 
      CALL stoperr(progname,mesg1,mesg2)
    ENDIF  
    IF (line(1:3) == 'END') EXIT        
    
   ! get next input
    nsat=nsat+1 
    READ(line,*,IOSTAT=ios) namsat(nsat),rpsat(nsat),eheat(nsat)
    IF (ios/=0) THEN
      IF (LEN_TRIM(line) /= 0) THEN 
        mesg1="while reading file pvap.sat (unexpected format?)"
        mesg2="at line: "//TRIM(line) 
        CALL stoperr(progname,mesg1,mesg2)
      ELSE
        mesg1="empty line - keyword END misssing?"
        CALL stoperr(progname,mesg1)
      ENDIF
    ENDIF  
  ENDDO
  CLOSE(filu)

  IF (nsat > SIZE(namsat)) THEN
    mesg1="nsat > SIZE(namsat) - size of tables too small"
    CALL stoperr(progname,mesg1)
  ENDIF

! get the ID for the species in namsat (store ID in idgsat)
! ---------------------------------------------------------
  loerr=.FALSE. ;  iloc=1

! fast seek - require same order in tables (OK in gecko)    
  seekloop1:&    
  DO i=1,nsat
    DO j=iloc,numsp
      IF (namsat(i)(ibeg:) == chrsp(j)) THEN
        iloc=j
        idgsat(i)=j
        mech2sat(j)=i
        CYCLE seekloop1
      ENDIF
    ENDDO
    
    ! if that point is reached then species not found (try unsorted loop)
    loerr=.TRUE. 
    EXIT seekloop1
  ENDDO seekloop1

! if fast seek did not succeed, then try slow seek (no order required - OK for MCM)
  IF (loerr) THEN  
    PRINT*, "       => try slow seek in readpvap"
    seekloop2: &   
    DO i=1,nsat
      DO j=1,numsp
        IF (namsat(i)(ibeg:) == chrsp(j)) THEN
          idgsat(i)=j
          mech2sat(j)=i
          CYCLE seekloop2
        ENDIF
      ENDDO
      mesg1="Species in the input file not found: "//TRIM(namsat(i))
      mesg2="The species is unknown in the mechanism."
      CALL stoperr(progname,mesg1,mesg2)
    ENDDO seekloop2
  ENDIF

! change evaporation heat from KJ/mol to K-1
! ------------------------------------------------------
  eheat(:)=eheat(:)*1000./Rgas
  lo1pvap = .TRUE. ! raise flag for 1st check of vapor pressure 

END SUBROUTINE readpvap

! ======================================================================
! Purpose : Compute the vapor pressure of every species in the 
! tables @ temp provided is input. Vapor pressure are computed based 
! on the Clausius Clapeyron equation.                              
! ======================================================================
SUBROUTINE getpvap(nsat,temp,rt4psat,rpsat,eheat,psat)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE
  INTEGER,INTENT(IN):: nsat             ! # of partitioning species 
  REAL,INTENT(IN)   :: temp             ! T for which Psat is computed
  REAL,INTENT(IN)   :: rt4psat          ! reference T for vapor pressure
  REAL,INTENT(IN)   :: rpsat(:)         ! vapor pressure @ ref T
  REAL,INTENT(IN)   :: eheat(:)         ! latent heat of evaporation (@ ref T in K-1)
  REAL,INTENT(OUT)  :: psat(:)          ! species saturation vapor pressure (atm)

  REAL    :: multifac
  INTEGER :: i 
  CHARACTER(LEN=10),PARAMETER :: progname='getpvap'
  CHARACTER(LEN=80) :: mesg1

! compute the saturation vapor pressure (atm) at temp. Use Clausius 
! Clapeyron to compute pvap @ temp if temp /= Tref
  IF (temp==rt4psat) THEN
      psat(1:nsat)=rpsat(1:nsat)

  ELSE                 
    multifac = ( (1./rt4psat)-(1./temp) )

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(psat,rpsat,eheat,multifac,nsat)
    DO i=1,nsat
      psat(i) = rpsat(i) * EXP(eheat(i)*multifac)
    ENDDO
!$OMP END PARALLEL DO 

  ENDIF

! check value if getpvap was called for the 1st time
  IF (lo1pvap) THEN
    DO i=1,nsat
      IF (psat(i) <= 0.) THEN
        mesg1="unexpected values found is computed pvap"
        CALL stoperr(progname,mesg1)
      ENDIF
    ENDDO
    lo1pvap=.FALSE.    
  ENDIF

END SUBROUTINE getpvap

! ======================================================================
! Purpose : Compute gas/aerosol partitioning assuming absorption into 
! an ideal  organic phase. The routine also returns the saturaration  
! vapor pressure.                                                     
! ======================================================================
SUBROUTINE soapartition(lout,nsat,idgsat,temp,nvoc,rt4psat,rpsat, &
                        eheat,cgas,caer)
  USE parameter_mod, ONLY: avogadro
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE
  INTEGER,INTENT(IN):: lout             ! file unit (report)
  INTEGER,INTENT(IN):: nsat             ! # of partitioning species 
  INTEGER,INTENT(IN):: idgsat(:)        ! ID of the partitioning species
  REAL,INTENT(IN)   :: temp             ! T for which Psat is computed
  REAL,INTENT(IN)   :: nvoc             ! Non Volatile Organic Compound (molec/cm3)
  REAL,INTENT(IN)   :: rt4psat          ! reference T for vapor pressure
  REAL,INTENT(IN)   :: rpsat(:)         ! vapor pressure @ ref T
  REAL,INTENT(IN)   :: eheat(:)         ! latent heat of evaporation (@ ref T in K-1)
  REAL,INTENT(INOUT):: cgas(:)          ! species gas phase concentration
  REAL,INTENT(INOUT):: caer(:)          ! species particle phase concentration 

  REAL    :: tratio(SIZE(caer))
  REAL    :: psat(SIZE(rpsat,1))        ! species saturation vapor pressure (atm)
  REAL    :: cmin, cmax, cinf, csup, ctotit
  REAL    :: conv, ctot, ratio, div
  REAL    :: multifac
  INTEGER :: it
  LOGICAL :: loconverge
  REAL, PARAMETER :: Rgas=0.0820578   ! Gas const. in L.atm./(K.mol)
  INTEGER, PARAMETER :: niter=100     ! maximum number of iteration allowed
  REAL, PARAMETER    :: cseed=1E3     ! min (seed) conc. - to start iteration
  LOGICAL, PARAMETER :: debug=.FALSE. ! turn to true if report wanted          

  CHARACTER(LEN=14),PARAMETER :: progname='soapartition'
  CHARACTER(LEN=80) :: mesg1, mesg2

! --------------------
! COMPUTE Pvap
! --------------------

! compute the saturation vapor pressure (atm) at temp. 
  CALL getpvap(nsat,temp,rt4psat,rpsat,eheat,psat)

! compute total concentration for each species (gas+aerosol) and provide
! a first estimate for the total number of molecule in the aerosol phase 
! based on the previous time step 
  ctotit=SUM(caer(1:nsat)) + nvoc
  cgas(idgsat(1:nsat)) = cgas(idgsat(1:nsat)) + caer(1:nsat)

! --------------------
! CONVERGENCE LOOP
! --------------------
! intialize parameters for the iteration loops
  loconverge=.FALSE.   ! to be turned to true once convergence achieved
  cmin = cseed + nvoc  ! min conc. for the sum of condensed species
  cmax = 1E20          ! max conc. for the sum of condensed species
  cinf = cmin          ! min. conc obtained during the iteration process
  csup = cmax          ! min. conc obtained during the iteration process
  IF (ctotit <= cmin) ctotit = cmin
  multifac = avogadro / (1000.*temp*Rgas)      ! Rgas here in L.atm./(K.mol)

! start iterations
  IF (debug) WRITE(lout,*) '--- convergence loop in soapartitioning --'    
  DO it = 1,niter
    ctot = 0. ; ratio = 0.
    div=multifac/ctotit
    caer(1:nsat) = cgas(idgsat(1:nsat)) / (1.+(psat(1:nsat)*div) )
    tratio(1:nsat) = cgas(idgsat(1:nsat)) / (ctotit+(psat(1:nsat)*multifac))
    ctot = SUM(caer(1:nsat)) + nvoc
    ratio = SUM(tratio(1:nsat)) + (nvoc/ctotit)

    IF (debug) WRITE(39,'(i3,3(1pe16.6))') it,ctotit,ctot, ratio

    ! Check convergence: reached if total aerosol conc. is modified by 
    ! by less than 1/1000 between 2 successive iterations 
    conv = (ctotit-ctot)/ctotit
    IF (ABS(conv) < 1E-3) THEN
      loconverge=.TRUE.
      EXIT
    ENDIF 

    ! ctot overestimated, concentration must be between cinf and ctot
    IF (ratio < 1.) THEN
      IF (ctotit <= (cseed+nvoc) ) THEN
        ctot = 0.
        loconverge = .TRUE.
        EXIT
      ENDIF 
      csup = ctotit       ! cinf is not changed 
      ctotit = (cinf*csup)**0.5

    ! ctot underestimated, concentration must be between ctot and csup
    ELSE
      cinf = ctotit       ! csup is not changed  
      ctotit = (cinf*csup)**0.5
    ENDIF
  ENDDO

! -----------------------
! CONVERGENCE  SUCCEEDED
! -----------------------
! If particle phase concentration is below a threshold of 1E3 molec/cm3 
! then reset aerosols to 0 else compute the gas/particle partioning.
  IF (loconverge) THEN
    IF (ctot == 0.) THEN
      caer(:) = 0.
    ELSE
      div=multifac/ctot
      caer(1:nsat) = cgas(idgsat(1:nsat)) / ( 1. +(psat(1:nsat)*div) )
      cgas(idgsat(1:nsat)) = cgas(idgsat(1:nsat)) - caer(1:nsat)
    ENDIF

! --------------------
! CONVERGENCE  FAILED
! --------------------
! Two cases can be considered for non convergence (Stop in any cases!):
! 1- If csup==cmax after iteration, then Cmax might be too low. 
! 2- The number of iteration (niter) is too low. 
  ELSE
    IF (csup == cmax) THEN
      mesg1="Convergence failed with csup=cmax"
      mesg2="Change the value for cmax?"
      CALL stoperr(progname,mesg1,mesg2)
    ELSE
      mesg1="Convergence failed it=niter"
      mesg2="Change the value for niter?"
      CALL stoperr(progname,mesg1,mesg2)
    ENDIF
  ENDIF

END SUBROUTINE soapartition

! ======================================================================
! PURPOSE: compute rate coefficient for mass transfer "reaction".
! ======================================================================
SUBROUTINE mtrat(numain,numaou,numwin,numwou,idain,idaou,idwin,idwou, &
                 woucf,temp,wmol,pnc,nsat,rt4psat,rpsat,eheat, &
                 mech2sat,gairefor,gwirefor,nvoc,cbox,idasat, &
                 rh,qfor)
  USE parameter_mod, ONLY: winfac,weqfac,ratfac,lifetimegw,fg_alpha,&
                           pi,avogadro,Rgas,fg_iterKgp,denssoa,Rp0
  USE viscositytool, ONLY: get_alpha,alpha               
  USE diffusivitytool, ONLY: dvol
  USE mechdata, ONLY: rxrgnt
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: numain        ! # of gas -> part. rxn
  INTEGER,INTENT(IN) :: numaou        ! # of part. -> gas rxn
  INTEGER,INTENT(IN) :: numwin        ! # of gas -> wall rxn
  INTEGER,INTENT(IN) :: numwou        ! # of wall -> gas rxn
  INTEGER,INTENT(IN) :: idain(:)      ! mechanism rx index for ith g->p (AIN)
  INTEGER,INTENT(IN) :: idaou(:)      ! mechanism rx index for ith p->g (AOU)
  INTEGER,INTENT(IN) :: idwin(:)      ! mechanism rx index for ith g->w (WIN)
  INTEGER,INTENT(IN) :: idwou(:)      ! mechanism rx index for ith w->g (WOU)
  REAL,INTENT(IN)    :: woucf(:,:)    ! [i,j] j extra coef. for the ith w->g rx (WOU)
  REAL,INTENT(IN)    :: temp          ! temperature
  REAL,INTENT(IN)    :: wmol(:)       ! molar mass of the ith species
  REAL,INTENT(IN)    :: pnc           ! particle number concentration
  INTEGER,INTENT(IN) :: nsat          ! # of partitioning species 
  REAL,INTENT(IN)    :: rt4psat       ! reference T for vapor pressure
  REAL,INTENT(IN)    :: rpsat(:)      ! vapor pressure @ ref T
  REAL,INTENT(IN)    :: eheat(:)      ! latent heat of evaporation (@ ref T in K-1)
  INTEGER,INTENT(IN) :: mech2sat(:)   ! [j] index of mech species j in the pvap table
  INTEGER,INTENT(IN) :: gairefor(:)   ! rxn ID of forward g->p matching the reverse p->g
  INTEGER,INTENT(IN) :: gwirefor(:)   ! rxn ID of forward g->w matching the reverse w->g
  REAL,INTENT(IN)    :: nvoc          ! Non Volatile OC (molec/cm3) (seed)
  REAL,INTENT(IN)    :: cbox(:)       ! concentration of species
  INTEGER,INTENT(IN) :: idasat(:)     ! ID of partitioning species - part. phase
  REAL,INTENT(IN)    :: rh            ! relative humidity
  REAL,INTENT(INOUT) :: qfor(:)       ! rate coefficient for all reactions

  REAL, PARAMETER :: RatmL=0.0820578  ! (atm L K-1 mol-1)
  REAL    :: psat(SIZE(rpsat,1))      ! species saturation vapor pressure (atm)
  REAL    :: Rp                       ! radius of the particles (cm)
  REAL    :: cmean(nsat)
  REAL    :: Dg(nsat)
  REAL    :: pressure
  REAL    :: soamass,mfac
  INTEGER :: ire,i,irefor
  REAL    :: multifac,keq,pvap,rate,tpfact
  REAL    :: ctotaer, caer
  INTEGER :: idsat,idspe

! return if no mass transfer equations
  IF (numain+numwin == 0) RETURN

! Compute vapor pressure
  CALL getpvap(nsat,temp,rt4psat,rpsat,eheat,psat)

! Compute mean molecular velocity (m/s)  
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(cmean,wmol,idasat,temp,nsat)
  DO i=1,nsat
    cmean(i) = SQRT((8.*Rgas*temp)/(pi*wmol(idasat(i))*0.001)) 
  ENDDO
!$OMP END PARALLEL DO

! Compute diffusivity (cm2/s)  
  pressure = 1.            ! pressure unit in Dg formula is bar. Assume 1 bar
  tpfact=temp**(7./4.)/pressure
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(tpfact,nsat,Dg,dvol)
  DO i=1,nsat
    Dg(i)=dvol(i)*tpfact
  ENDDO
!$OMP END PARALLEL DO

! --------------------------
! gas <-> aerosol parameters
! --------------------------
  IF (fg_alpha==1) THEN
    CALL get_alpha(cbox,nvoc,idasat,nsat,rh,psat,wmol,temp,pnc,cmean,alpha)
  ELSE
    alpha(:)=1. 
  ENDIF
  
! --------------
! gas -> aerosol
! --------------

! compute the radius (cm) of the particles (assume pnc is unchanged)
  mfac=1E12/avogadro
  soamass= SUM(cbox(idasat(1:nsat))*wmol(idasat(1:nsat)))*mfac ! ug/m3
  !seedmass=seed*Mseed*mfac                                     ! ug/m3

  !! option 1: mix seed and SOA to compute radius (
  !Rp = ( ((soamass+seedmass)*1E-12/(denssoa*pnc))*(3./(4.*pi)) )**(1./3.)
  !! option 2: make SOA shell around seed (of radius Rp0)
  Rp = (Rp0**3 + (soamass*1E-12/(denssoa*pnc))*(3./(4.*pi)))**(1./3.)
  !Rp = 1E-5   ! overwrite (assume 100 nm, in cm)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire,idspe,idsat,rate) &
!$OMP SHARED(numain,idain,rxrgnt,Dg,mech2sat,qfor,alpha,pnc,Rp,cmean)
  DO i=1,numain
    ire=idain(i)
    idspe=rxrgnt(ire)%idrg(1) 
    idsat=mech2sat(idspe)
    rate = g2prate(Dg(idsat),pnc,Rp,cmean(idsat),alpha(idsat))
    qfor(ire) = rate*ratfac
  ENDDO
!$OMP END PARALLEL DO 

! --------------
! aerosol -> gaz
! --------------
! Reverse reaction based on the equilibrium constant Keq=k_in/k_out where
! k_in is the forward reaction (see gas -> aerosol). The equilibrium 
! constant Keq=Caer/Cgas (concentration in molecule/cm3 of air) is: 
! Keq=(RT*Caer)/(Nav*Pvap).

  ! compute species concentration in particles (per cm3 of air)
  caer=SUM(cbox(idasat(1:nsat)))
  ctotaer=nvoc+caer         ! seed = nvoc (non volatile organic compound)
  IF (ctotaer<1E8) ctotaer=1E8  ! overwrite "tiny seed" if ctotaer is too small 

  multifac = RatmL*1000.*temp/avogadro     ! R in atm.L.K-1.mol-1 and 1000 cm3/L

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire,idspe,idsat,pvap,keq,irefor) &
!$OMP SHARED(numaou,idaou,rxrgnt,mech2sat,psat,multifac,gairefor,qfor,ctotaer)
  DO i=1,numaou
    ire=idaou(i)
    idspe=rxrgnt(ire)%idrg(1) 
    idsat=mech2sat(idspe)
    pvap=psat(idsat)
    
    ! Compute K equilibrium (see also iteration.f90)
    IF (fg_iterKgp) THEN ; keq=multifac/pvap          ! scaling to ctotaer performed in iteration
    ELSE                 ; keq=multifac*ctotaer/pvap  ! constant ctotaer between time steps
    ENDIF

    ! compute kout as kin/keq
    irefor=gairefor(i)            ! irefor: ID # for the forward rxn
    qfor(ire)=qfor(irefor)/keq      
  ENDDO
!$OMP END PARALLEL DO 

! --------------
! gas -> wall
! --------------
! Lifetime related to the gas -> wall loss is provided as "lifetimegw".
! All species are assumed here to have the same rate coefficient. 
  IF (numwin>0) rate=(1./lifetimegw)*winfac  
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire) SHARED(qfor,numwin,idwin,rate)
  DO i=1,numwin
    ire=idwin(i)
    qfor(ire) = rate
  ENDDO
!$OMP END PARALLEL DO 
      
! --------------
! wall -> gaz
! --------------
! Reverse reaction based on the equilibrium constant Keq=kgw/kwg where
! kwg is the forward reaction (see gas -> wall). The equilibrium 
! constant Keq=Cwall/Cgas (concentration in molecule/cm3 of air) is: 
!   Keq=(RT/Pvap)*(cmw)
! where cmw=(Cw/Mw*gamma) [mole/m3] is provided in the woucf.
! See Ziemann et al., AST, 2010, p887).
  multifac = RatmL*1E-3*temp ! RT in atm.m3.mol-1 

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire,idspe,idsat,pvap,keq,irefor) &
!$OMP SHARED(numwou,idwou,rxrgnt,mech2sat,psat,multifac,woucf,gwirefor,qfor)
  DO i=1,numwou
    ire=idwou(i)
    idspe=rxrgnt(ire)%idrg(1) 
    idsat=mech2sat(idspe)
    pvap=psat(idsat)

    ! compute the rate constant for the transformation
    keq=multifac*woucf(i,1)*weqfac/pvap

    ! compute kout as kin/keq
    irefor=gwirefor(i)            ! irefor is the ID # for the forward rxn
    qfor(ire)=qfor(irefor)/keq      
  ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE mtrat

! ======================================================================
! PURPOSE: compute the gas -> particle rate coefficient.
! The diffusion constant Dg is provided as input.
! ======================================================================
REAL FUNCTION g2prate(Dg,pnc,Rp,mspeed,alpha)
  USE parameter_mod, ONLY: Cair,imtr,Mair,pi,Rgas,Sair
  IMPLICIT NONE

  REAL,INTENT(IN) :: Dg        ! diffusivity (cm2/s) the species in air
  REAL,INTENT(IN) :: pnc       ! particle number concentration
  REAL,INTENT(IN) :: Rp        ! radius of the particles (cm)
  REAL,INTENT(IN) :: mspeed    ! mean molecular speed (m/s) 
  REAL,INTENT(IN) :: alpha     ! mass accommodation coefficient 

  REAL :: ktd,kti,ktot,Kn,beta
 
! Kn Knudsen number (mean free path: Lambda=3Dg/mspeed)
  Kn = (3.*Dg)/(100.*mspeed*Rp)
  beta = (0.75*alpha*(1+Kn))/(Kn**2. + Kn + 0.283*Kn*alpha +0.75*alpha)

! 1st order rate coef for diffusion limited mass transfer (continuous regime)
  ktd = 4.*pi*Dg*pnc*Rp

! 1st order rate coef for collision limited mass transfer (kinetic regime)
  kti = alpha*mspeed*100.*pi*(Rp**2)*pnc

! return the rate coef
  IF      (imtr == 1) THEN ; g2prate=ktd       ! limited by gas-phase diffusion
  ELSE IF (imtr == 2) THEN ; g2prate=kti       ! limited by gas/particles collision
  ELSE IF (imtr == 3) THEN ; g2prate=ktd*beta  ! Fuchs-Sutugin approach
  ELSE                                         ! resistance approach
    ktot=(1./ktd)+(1./kti)  
    g2prate=1./ktot
  ENDIF
      
END FUNCTION g2prate

END MODULE vaporpressure
