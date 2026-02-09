MODULE viscositytool
!$ USE omp_lib
IMPLICIT NONE

INTEGER             :: ntg         ! # of partitioning species
REAL,ALLOCATABLE    :: Tg(:)       ! Tg
REAL,ALLOCATABLE    :: alpha(:)    ! accommodation coefficient  

CONTAINS 
! SUBROUTINE get_alpha(cbox,seed,idasat,nsat,rh,psat,wmol,temp,pnc,alpha)
! SUBROUTINE readTg(chrsp,numsp,lodicnam)

!=======================================================================
! PURPOSE: compute the mass accomodation coefficient alpha for each 
! species.
!=======================================================================
SUBROUTINE get_alpha(cbox,seed,idasat,nsat,rh,psat,wmol,temp,pnc,cmean,alpha)
  USE parameter_mod, ONLY: alphas,avogadro,D,denssoa,densw,fg_alpha_bulk, &
                           fg_tgseed,kappa,kb,kGT,pi,rad,Rgas,Rp0,Tgseed, &
                           Mseed,tgw
  IMPLICIT NONE

  REAL,INTENT(IN)    :: cbox(:)    ! species concentration
  REAL,INTENT(IN)    :: seed       ! Non Volatile OC (molec/cm3) (NVOC)
  INTEGER,INTENT(IN) :: idasat(:)  ! ID of partitioning species - part. phase
  INTEGER,INTENT(IN) :: nsat       ! # of partitioning species
  REAL,INTENT(IN)    :: rh         ! relative humidity
  REAL,INTENT(IN)    :: psat(:)    ! species saturation vapor pressure (atm)
  REAL,INTENT(IN)    :: wmol(:)    ! molar weight of species
  REAL,INTENT(IN)    :: temp       ! temperature
  REAL,INTENT(IN)    :: pnc        ! particle number concentration
  REAL,INTENT(IN)    :: cmean(:)   ! mean molecular velocity (m/s)
  REAL,INTENT(OUT)   :: alpha(:)   ! accommodation coefficient

  INTEGER :: i
  REAL    :: Zgeq(SIZE(psat))      ! equilibrium gas phase concentration (molec.cm-3)
  REAL    :: visc                  ! viscosity (Pa.s)
  REAL    :: Db                    ! bulk diffusivity (cm2.s-1)
  REAL    :: Rp                    ! particle radius (cm)
  REAL    :: soamass               ! SOA mass concentration (µg/m3)
  REAL    :: wmass                 ! water mass concentration (condensed)
  REAL    :: seedmass              ! seed mass concentration (µg/m3)
  REAL    :: fsoa,mspeed
  REAL    :: Tgdry,Tgwet,T0
  REAL,PARAMETER :: mfac=1E12/avogadro

  alpha(:)=1.
  
! compute SOA mass concentration
  soamass= SUM(cbox(idasat(1:nsat))*wmol(idasat(1:nsat)))*mfac
  seedmass=seed*Mseed*mfac
  IF (soamass==0) THEN ; visc=0. ; RETURN ; ENDIF ! exit if no SOA
  
! compute Tg for dry SOA mixture
  Tgdry= SUM(cbox(idasat(1:nsat))*wmol(idasat(1:nsat))*Tg(1:nsat))*mfac/soamass
  IF (fg_tgseed==1) THEN
    Tgdry= (soamass/(soamass+seedmass))*Tgdry + (seedmass/(soamass+seedmass))*Tgseed
  ENDIF

! compute water mass concentration at given RH 
  IF (rh==0.) THEN ; wmass=0.
  ELSE             ; wmass=(kappa*densw*soamass)/(denssoa*((100./rh)-1))
  ENDIF

! calculate Tg for SOA+water mixture
  IF (fg_tgseed==0) THEN ; fsoa = soamass / (soamass+wmass)
  ELSE                   ; fsoa = (soamass+seedmass) / (soamass+wmass+seedmass)
  ENDIF
  Tgwet = ((1-fsoa)*Tgw  + (1/kGT)*fsoa*Tgdry) / ((1-fsoa) + (1/kGT) *fsoa)

! calculate viscosity
  T0= (39.17 * Tgwet) / (D+39.17)
  visc = 10**(-5 + 0.434 * T0*D / (temp-T0))
  
! convert viscosity to bulk diffusivity (cm2.s-1)
  Db = (kb*temp*1E4) / (6*pi*rad*visc)

! compute equilibrium gas phase concentration for each species
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(Zgeq,psat,temp,nsat)
  DO i=1,nsat
    Zgeq(i)= (psat(i)*1.001325E5*avogadro) / (Rgas*temp*1E6)
  ENDDO
!$OMP END PARALLEL DO
  
! compute particle radius (cm)
  Rp = (Rp0**3 + (soamass*1E-12/(denssoa*pnc))*(3./(4.*pi)))**(1./3.)

! compute alpha
  IF (fg_alpha_bulk==0) THEN

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,mspeed) SHARED(alpha,Zgeq,cmean,nsat,Db,Rp)
    DO i=1,nsat
      mspeed=cmean(i)*100.   ! change unit from m/s to cm/s
      alpha(i) = alphas / (1+ ( ((rad*100)**4 *alphas*mspeed*Zgeq(i)) / (2*Db) ) )
    ENDDO
!$OMP END PARALLEL DO

  ELSE

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,mspeed) SHARED(alpha,Zgeq,cmean,nsat,Db,Rp)
    DO i=1,nsat
      mspeed=cmean(i)*100.  
      alpha(i) = alphas / (1+ ( ((rad*100)**4 *alphas*mspeed*Zgeq(i))/(4.*Db) ) * &
                              (1.+(Rp/(2.*rad*100))) )
    ENDDO
!$OMP END PARALLEL DO

  ENDIF
  
END SUBROUTINE get_alpha

!=======================================================================
! Purpose: read Tg to calculate viscosity, accomodation coefficient and 
! rate for mass tranfer.
!=======================================================================
SUBROUTINE readTg(chrsp,numsp,namsat,idgsat,lodicnam)
  USE boxtool, ONLY: stoperr
  USE parameter_mod, ONLY: mxlsp, input_dir_chem
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: numsp             ! # of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:) ! name of the species
  CHARACTER(LEN=*),INTENT(IN):: namsat(:) ! name of the partitioning species (with 'G')
  INTEGER,INTENT(IN) :: idgsat(:)         ! ID of the partitioning species
  LOGICAL,INTENT(IN),OPTIONAL :: lodicnam ! if true and present then chrsp is dictionary names

  CHARACTER(LEN=mxlsp) :: namtg           ! name of the partitioning species
  CHARACTER(LEN=160) :: line
  INTEGER            :: i,j,iloc,ios,filu,ibeg,nrec,itg
  LOGICAL            :: loerr
  CHARACTER(LEN=15),PARAMETER :: progname='readTg'
  CHARACTER(LEN=80) :: mesg1, mesg2
  
! chrsp is either a "gas phase name" with 1st letter G (boxmodel) or the
! dictionary name (no "G" for postprocessing). Default is gasphase names.
  ibeg=1 
  IF (PRESENT(lodicnam)) THEN
    IF (lodicnam) THEN ; ibeg=2 ; ELSE; ibeg=1 ; ENDIF
  ENDIF
  
! open the file
! -------------
  filu=12
  i=LEN_TRIM(input_dir_chem)
  IF (input_dir_chem(i:i)=='/') THEN
    OPEN (filu, FILE=TRIM(input_dir_chem)//'Tg.dat', STATUS='OLD')
  ELSE
    OPEN (filu, FILE=TRIM(input_dir_chem)//'.Tg', STATUS='OLD')
  ENDIF
  
! read the number of record 
  READ (filu,*, IOSTAT=ios) nrec
  IF (ios/=0) THEN 
    mesg1="# record cannot be read - keyword END missing?"
    CALL stoperr(progname,mesg1)
  ENDIF  
  IF (nrec/=ntg) THEN 
    mesg1="# record does not match size used to allocate table memory"
    CALL stoperr(progname,mesg1)
  ENDIF  

! read the file
! -------------
  itg=0
  DO     
    READ (filu, '(a)', IOSTAT=ios) line

    ! exit if keyword END found
    IF (ios/=0) THEN 
      mesg1="keyword END missing in the input file?"
      CALL stoperr(progname,mesg1)
    ENDIF  
    IF (line(1:3) == 'END') EXIT        
    
    ! get next input 
    itg=itg+1                  
    READ(line,*,IOSTAT=ios) namtg,Tg(itg)
    IF (ios/=0) THEN
      IF (LEN_TRIM(line) /= 0) THEN 
        mesg1="unexpected format? - see out file"
        mesg2="at line: "//TRIM(line) 
        CALL stoperr(progname,mesg1,mesg2)
      ELSE
        mesg1="empty line - keyword END misssing?"
        CALL stoperr(progname,mesg1)
      ENDIF
    ENDIF
    
    ! check that sorting in Tg and Pvap tables are identical
    IF (namtg/=namsat(itg)) THEN
      mesg1="species in namsat (pvap) and namtg (Tg) are not identical"
      mesg2="same order is required. Mismatch at species: "//TRIM(namtg)
      CALL stoperr(progname,mesg1,mesg2)
    ENDIF  
  ENDDO
  CLOSE(filu)

  IF (itg/=nrec) THEN
    mesg1="number of record /= number of data"
    CALL stoperr(progname,mesg1)
  ENDIF
  IF (itg > SIZE(Tg)) THEN
    mesg1="itg > SIZE(Tg) - The size of the tables is too small"
    CALL stoperr(progname,mesg1)
  ENDIF

! check the ID for the species in namsat 
! --------------------------------------
  loerr=.FALSE.
  iloc=1

! fast seek - require same order in tables (OK in gecko)    
  seekloop1:&    
  DO i=1,ntg
    DO j=iloc,numsp
      IF (namsat(i)(ibeg:) == chrsp(j)) THEN
        iloc=j
        ! check the consistency with the Pvap data (must be identical)
        IF (j/=idgsat(i)) THEN
          mesg1="mismatch in the ID between Tg and Pvap data"
          mesg2="Mismatch at species: "//namsat(i)
          CALL stoperr(progname,mesg1,mesg2)
        ENDIF
        CYCLE seekloop1
      ENDIF
    ENDDO
    
    ! if that point is reached then species not found (try unsorted loop)
    loerr=.TRUE. 
    EXIT seekloop1
  ENDDO seekloop1

! if fast seek did not succeed, then try slow seek (no order required - OK for MCM)
  IF (loerr) THEN  
    WRITE(6,*) "       => try slow seek in readTg"
    seekloop2: &   
    DO i=1,ntg
      DO j=1,numsp
        IF (namsat(i)(ibeg:) == chrsp(j)) THEN
          ! check the consistency with the Pvap data (must be identical)
          IF (j/=idgsat(i)) THEN
            mesg1="mismatch in the ID between Tg and Pvap data"
            mesg2="Mismatch at species: "//namsat(i)
            CALL stoperr(progname,mesg1,mesg2)
          ENDIF
          CYCLE seekloop2
        ENDIF
      ENDDO
      mesg1="Species in the input file not found: "//TRIM(namsat(i))
      mesg2="The species is unknown in the mechanism."
      CALL stoperr(progname,mesg1,mesg2)
    ENDDO seekloop2
  ENDIF

END SUBROUTINE readTg

END MODULE viscositytool
