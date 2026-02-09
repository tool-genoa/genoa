MODULE initnbound
IMPLICIT NONE
REAL    :: temp0
REAL    :: rh0
INTEGER :: iyear     ! date - year
INTEGER :: imonth    ! date - month
INTEGER :: iday      ! date - day
REAL    :: sla       ! latitude,
REAL    :: slo       ! longitude
REAL    :: tz        ! time zone

CONTAINS
!=======================================================================  
! Purpose: return T, RH, water for the time provided as input. User 
! can set prescribed function, as needed.
!=======================================================================  
SUBROUTINE getconstrain(time,temp,sumc,rh,water,envfix)
  USE parameter_mod, ONLY: avogadro,Rgas
  IMPLICIT NONE

  REAL, INTENT(IN)    :: time        ! time (warning: modulo 1 day)
  REAL, INTENT(IN)    :: envfix(:,:) ! Fixed environmental data read from key file
  REAL, INTENT(OUT)   :: temp        ! temperature
  REAL, INTENT(OUT)   :: sumc        ! M (in molec/cm3)
  REAL, INTENT(OUT)   :: rh          ! relative humidity
  REAL, INTENT(OUT)   :: water       ! water concentration (molec/cm3)
  
! "Magnus formula": psat(H2O) = c1 * EXP(c2*T/(c3+T)) (units Pa)
!  where: T is expressed in degrees ! : T(K)-273.16
!         c1 = 610.94, c2 = 17.625, c3 = 243.04
! Ref: Alduchov & Eskridge 1996, J.Appl.Meteorol. 35 601-609.
! quoted by: Lawrence, 2005, BAMS DOI:10.1175/BAMS-86-2-225
  REAL,PARAMETER :: c1 = 610.94, c2 = 17.625,  c3 = 243.04
  REAL           :: TC, psat_H2O, p_H2O
  INTEGER        :: i

  ! Ger index for fixed environmental parameters - genoa
  i = MAX(1, MIN(24, CEILING(time/3600.))) ! time in hours

! set temperature
! --------------------
  ! As provided in the input file:
  IF (envfix(1,2) > 0.) THEN ! Use the fixed value
    temp=envfix(2,i)
  ELSE ! Read from the input file
    temp=temp0 
  ENDIF

! set pressure (sumc in molec/cm3) 
! --------------------------------
  ! Use perfect gas law :sumc=(pres*6.022E+23)/(8.32*temp)
  ! Constant value (forced by user):
  IF (envfix(1,5) > 0.) THEN ! Use the fixed value
    sumc=envfix(5,i)
  ELSE ! Read from the input file
    sumc=2.5E19
  ENDIF

! set relative humidity (in %) and the water content 
! -------------------------------------
  ! As provided in the input file:
  IF (envfix(1,3) > 0.) THEN ! Use the fixed value
    rh=envfix(3,i)
  ELSE ! Read from the input file
    rh=rh0
  ENDIF

! Water concentration (molec/cm3)
! -------------------------------
  TC = temp-273.16                   ! temperature in Celsius
  psat_H2O = c1*EXP(c2*TC/(c3+TC))   ! saturation water pressure
  p_H2O = psat_H2O*rh/100.           ! water pressure
  water = p_H2O / (Rgas*temp*1.e6/avogadro)

END SUBROUTINE getconstrain

!=======================================================================  
! PURPOSE: Read the initial concentrations of some species
!=======================================================================  
SUBROUTINE readinitial(lout,numsp,chrsp,cbox,nfix,idfix,cfix,noxfix,cnox,envfix)
  USE stringdata, ONLY: stringarg
  USE boxtool, ONLY: spe2id,stoperr
  USE parameter_mod, ONLY: smallc,input_cond
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: lout              ! file unit (report)
  INTEGER,INTENT(IN)  :: numsp             ! # of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)  ! name of the species
  REAL,INTENT(OUT)    :: cbox(:)           ! species concentration
  INTEGER,INTENT(OUT) :: nfix              ! number of species having fixed concentration
  INTEGER,INTENT(OUT) :: idfix(:)          ! id of the species having fixed concentration
  REAL,INTENT(OUT)    :: cfix(:,:)         ! value of the fixed concentration
  REAL,INTENT(OUT)    :: envfix(:,:)       ! value of the fixed environmental parameters - genoa
                                           ! (temperature, relative humidity, solar zenith angle, pressure)
  INTEGER,INTENT(OUT) :: noxfix            ! flag to fix NOx (NO+O2) concentration 
  REAL,INTENT(OUT)    :: cnox              ! value of the fixed NOx concentration

  CHARACTER(LEN=LEN(chrsp(1))) :: species
  LOGICAL            :: lostop
  CHARACTER(LEN=76)  :: line
  CHARACTER(LEN=400) :: longline
  CHARACTER(LEN=4)   :: keyword
  INTEGER            :: i,isp,filu,ios
  REAL               :: cin                ! genoa - input concs

! parameter for subroutine stringarg (see stringdata module)   
  INTEGER, PARAMETER :: mxarg=10      ! maximum # of arguments in a line
  INTEGER :: posi(mxarg,2)
  INTEGER :: narg

  CHARACTER(LEN=15),PARAMETER :: progname='readinitial'
  CHARACTER(LEN=80) :: mesg1, mesg2

! -----------------------------
! INITIALIZE               
! -----------------------------
  lostop=.FALSE.

! set default value (likely overwritten next)
  cbox(:) = smallc                    ! lower limit for a concentration (avoid 0.) 
  nfix=0  ; idfix(:)=0 ; cfix(:,:)=0. ; envfix(:,:)=-1. ! genoa
  noxfix=0 ; cnox=0.
  
! -----------------------------
! OPEN THE FILE 
! -----------------------------
  filu=12 ! 'indat.key'
  OPEN(filu,FILE=input_cond,FORM='formatted',STATUS='old',IOSTAT=ios)
  IF (ios /= 0) THEN
    mesg1="file "//TRIM(input_cond)//" not found"
    CALL stoperr(progname,mesg1)
  ENDIF

  !WRITE(lout,*)
  !WRITE(lout,*)'--------- Keyword input -----------'

! -----------------------------
! READ NEXT LINE
! -----------------------------
  rdline: & 
  DO                 ! escape the loop if keyword 'END ' is found 
    READ(filu,'(a4,(a))',IOSTAT=ios) keyword, line
    IF (ios/=0) THEN 
      mesg1="End of file reached - keyword 'END' missing ?"
      CALL stoperr(progname,mesg1)
    ENDIF

! comment line - read next
    IF (keyword(1:1)=='/' .OR. keyword(1:1)=='!') CYCLE rdline

! CASE SELECTOR FOR KEYWORDS
! -----------------------------
    SELECT CASE (keyword)   

     ! End of file - escape
     CASE ('END ') 
       EXIT rdline

     ! Initial concentrations
     CASE('REAC')
       CALL stringarg(line,narg,posi)
       IF (narg /= 2) THEN            ! species + c*inbox
         WRITE(lout,*) " --error--, in readinitial while reading reac."
         WRITE(lout,*) "            # of arg is not appropriate:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       species=line(posi(1,1):posi(1,2))
       isp=spe2id(lout,chrsp,numsp,species,csub='readinitial.f90',loerr=.FALSE.)
       IF (isp == 0) THEN
         WRITE(lout,*) " --error--, in readinitial while reading reac."
         WRITE(lout,*) "            species is unkown in line :"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       READ(line(posi(2,1):),*,IOSTAT=ios) cbox(isp)
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading reac:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       CYCLE rdline

     ! fixed concentrations
     CASE('CFIX')
       nfix = nfix+1
       IF (nfix > SIZE(idfix)) THEN
         WRITE(lout,*) " --error--, in readinitial while reading cfix."
         WRITE(lout,*) "            # of species exceed max allowed: ",SIZE(idfix)
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF

       CALL stringarg(line,narg,posi)
       IF (narg /= 2) THEN            
         WRITE(lout,*) " --error--, in readinitial while reading cfix."
         WRITE(lout,*) "            # of arg is not appropriate:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       species=line(posi(1,1):posi(1,2))
       isp=spe2id(lout,chrsp,numsp,species,csub='readinitial.f90',loerr=.FALSE.)
       IF (isp == 0) THEN
         WRITE(lout,*) " --error--, in readinitial while reading reac."
         WRITE(lout,*) "            species is unkown in line :"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       idfix(nfix)=isp
       ! Read fixed concentration and copy ! genoa
       READ(line(posi(2,1):),*,IOSTAT=ios) cin
       DO i=1, 24
         cfix(nfix,i) = cin
       ENDDO
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading reac:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       CYCLE rdline

     ! fixed hourly concentrations - genoa
     CASE('CF24')
       nfix = nfix+1
       IF (nfix > SIZE(idfix)) THEN
         WRITE(lout,*) " --error--, in readinitial while reading cf24."
         WRITE(lout,*) "            # of species exceed max allowed: ",SIZE(idfix)
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF

       CALL stringarg(line,narg,posi)
       IF (narg /= 1) THEN            
         WRITE(lout,*) " --error--, in readinitial while reading cf24."
         WRITE(lout,*) "            # of arg is not appropriate:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       species=line(posi(1,1):posi(1,2))
       isp=spe2id(lout,chrsp,numsp,species,csub='readinitial.f90',loerr=.FALSE.)
       IF (isp == 0) THEN
         WRITE(lout,*) " --error--, in readinitial while reading cf24."
         WRITE(lout,*) "            species is unkown in line :"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       idfix(nfix)=isp
       ! Read hourly fixed concentration from the next line ! genoa
       READ(filu,'(a)',IOSTAT=ios) longline
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading concs from cf24:"
         WRITE(lout,*) keyword,' ',TRIM(line), ' ',TRIM(longline)
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       longline=TRIM(ADJUSTL(longline))
       ! Read the 24 hourly fixed concentrations 
       READ(longline,*,IOSTAT=ios) (cfix(nfix,i),i=1,24)
       !PRINT*, "Read ", TRIM(line), " ", nfix, cfix(nfix,:)
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading reac from cf24 at hour ",i
         WRITE(lout,*) keyword,' ',TRIM(line), ' ',TRIM(longline)
         lostop=.TRUE.
         CYCLE rdline
       ELSE
         PRINT*, "Read diurnal fixed concs for ", TRIM(line), " ", nfix, cfix(nfix,:)
       ENDIF
       CYCLE rdline

     ! fixed environmental parameters - genoa
     CASE('TEMP') ! temperature (K)
       READ(filu,'(a)',IOSTAT=ios) longline
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, while reading the 2nd line for temperature TEMP (K):"
         WRITE(lout,*) keyword,' ',TRIM(longline)
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       longline=TRIM(ADJUSTL(longline)) ! remove leading and trailing spaces
       READ(longline,*,IOSTAT=ios) (envfix(2,i),i=1,24)
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, while reading values for temperature TEMP (K):"
         WRITE(lout,*) keyword, ' ',TRIM(longline), envfix(2,:)
         lostop=.TRUE.
       ELSE
         envfix(1,2) = 1. ! Set Flag
         PRINT*, "Read diurnal temperature ", envfix(2,:)
       ENDIF 
       CYCLE rdline

     CASE('RELH') ! relative humidity (%)
        READ(filu,'(a)',IOSTAT=ios) longline
        IF (ios/=0) THEN
          WRITE(lout,*) "--error--, while reading the 2nd line for relative humidity RELH (%):"
          WRITE(lout,*) keyword, ' ',TRIM(longline)
          lostop=.TRUE.
          CYCLE rdline
        ENDIF
        longline=TRIM(ADJUSTL(longline))
        READ(longline,*,IOSTAT=ios) (envfix(3,i),i=1,24)
        IF (ios/=0) THEN
          WRITE(lout,*) "--error--, while reading values for relative humidity RELH (%):"
          WRITE(lout,*) keyword,' ',TRIM(longline), envfix(3,:)
          lostop=.TRUE.
       ELSE
         envfix(1,3) = 1. ! Set Flag
         PRINT*, "Read diurnal relative humidity ", envfix(3,:)
       ENDIF 
        CYCLE rdline

     CASE('SZAS') ! solar zenith angle (degree)
        READ(filu,'(a)',IOSTAT=ios) longline
        IF (ios/=0) THEN
          WRITE(lout,*) "--error--, while reading the 2nd line for solar zenith angle SZAS (degree):"
          WRITE(lout,*) keyword,' ',TRIM(longline)
          lostop=.TRUE.
          CYCLE rdline
        ENDIF
        longline=TRIM(ADJUSTL(longline))
        READ(longline,*,IOSTAT=ios) (envfix(4,i),i=1,24)
        IF (ios/=0) THEN
          WRITE(lout,*) "--error--, while reading values for solar zenith angle SZAS (degree):"
          WRITE(lout,*) keyword,' ',TRIM(longline), envfix(4,:)
          lostop=.TRUE.
       ELSE
         envfix(1,4) = 1. ! Set Flag
         PRINT*, "Read diurnal solar zenith angle ", envfix(4,:)
       ENDIF 
        CYCLE rdline

     CASE('PRES') ! pressure (molec/cm3)
        READ(filu,'(a)',IOSTAT=ios) longline
        IF (ios/=0) THEN
          WRITE(lout,*) "--error--, while reading the 2nd line for pressure PRES (molec/cm3):"
          WRITE(lout,*) keyword,' ',TRIM(longline)
          lostop=.TRUE.
          CYCLE rdline
        ENDIF
        longline=TRIM(ADJUSTL(longline))
        READ(longline,*,IOSTAT=ios) (envfix(5,i),i=1,24)
        IF (ios/=0) THEN
          WRITE(lout,*) "--error--, while reading values for pressure PRES (molec/cm3):"
          WRITE(lout,*) keyword,TRIM(longline), envfix(5,:)
          lostop=.TRUE.
       ELSE
         envfix(1,5) = 1. ! Set Flag
         PRINT*, "Read diurnal pressure ", envfix(5,:)
       ENDIF 
       CYCLE rdline

     ! Fixed NOx concentration
     CASE('CNOX')
       READ(line,*,IOSTAT=ios) cnox
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading CNOX:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
       ENDIF
       noxfix=1                       ! raise the "NOx fixed" flag
       CYCLE rdline

    ! Unidentified keyword
    CASE DEFAULT
      WRITE(lout,*)'--error--, in readinitial. Keyword unknown:'
      WRITE(lout,*) keyword,' ',TRIM(line) 
      lostop=.TRUE.
      CYCLE rdline

    END SELECT  
  ENDDO rdline
  CLOSE(filu)

! Check inconsistencies in inputs
! -------------------------------

  ! if NOx is fixed then NO and NO2 must be "free" 
  IF (noxfix/=0) THEN
    IF (nfix/=0) THEN
      DO i=1,nfix
        IF (chrsp(idfix(i))=='GNO ') THEN
          WRITE(lout,*) "--error--, in readinitial. Inconsistency in"
          WRITE(lout,*) " input data: both NO and NOx are fixed."
          lostop=.TRUE.
        ENDIF
        IF (chrsp(idfix(i))=='GNO2 ') THEN
          WRITE(lout,*) "--error--, in readinitial. Inconsistency in"
          WRITE(lout,*) " input data: both NO2 and NOx are fixed."
          lostop=.TRUE.
        ENDIF
      ENDDO
    ENDIF
  ENDIF

! -------------------
! STOP if error found
! -------------------
  IF (lostop) THEN
    mesg1="Miscellaneous errors found while reading "//TRIM(input_cond) &
    //" (see *.out file)"
    mesg2="See details in the *.out file"
    CALL stoperr(progname,mesg1,mesg2)
  ENDIF

END SUBROUTINE readinitial

!=======================================================================  
! PURPOSE: Read the initial conditions (temperature, initial C ...) 
! for the simulation and parameters for the time integration (start and 
! stop time, # of time steps, solver tolerance ...) in the namelist
!=======================================================================  
SUBROUTINE rd_nml(filein)
  USE parameter_mod, ONLY: flag_ZA,fixedZA,tstart,tstop,ntstep,nskip,&
                           rtol,atol,numit,dtmin,&
                           nvoc,nvic,Rp0,denssoa,densseed,Mseed, &
                           ! genoa
                           fg_out_col,err_sps_list_in, &
                           fg_output,input_dir_chem,input_cond, &
                           output_dir,phot_file,ref_conc_list_in, &
                           init_conc_list_in

  USE mechdata,      ONLY: kdilu, kaerloss
  IMPLICIT NONE

  CHARACTER(LEN=200), INTENT(IN) :: filein ! 'simu.nml'
  INTEGER,PARAMETER           :: nmlu=45
  INTEGER                     :: ierr
  
  NAMELIST /simu_conditions/ tstart,tstop,ntstep,nskip,rtol,atol,numit,dtmin
  NAMELIST /env_conditions/  temp0,rh0,flag_ZA,fixedZA, kdilu, kaerloss,&
                             iday,imonth,iyear,sla,slo,tz
  NAMELIST /mass_transfert_parameters/ Rp0,nvoc,nvic,denssoa,densseed,Mseed
  
  ! genoa
  NAMELIST /genoa_options/ fg_out_col,err_sps_list_in,fg_output, &
                           input_dir_chem,input_cond,output_dir,phot_file, &
                           ref_conc_list_in,init_conc_list_in

  OPEN(nmlu, file=filein, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading namelist file'
      STOP
    ENDIF

    READ(nmlu,nml=simu_conditions, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "simu_conditions" arguments in namelist'
      STOP
    ENDIF

    READ(nmlu,nml=env_conditions, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "env_conditions" arguments in namelist'
      STOP
    ENDIF

    READ(nmlu,nml=mass_transfert_parameters, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "env_conditions" arguments in namelist'
      STOP
    ENDIF
    
    READ(nmlu,nml=genoa_options, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "genoa_options" arguments in namelist'
      STOP
    ENDIF
  CLOSE(nmlu)

  IF ( (imonth <1 .OR. imonth > 12) .OR. &
       (iday  < 1 .OR. iday   > 31)       ) THEN 
    PRINT*, 'problem with the date: ',iday,'/',imonth,'/',iyear
    STOP
  ENDIF
  
  ! genoa
  IF (fg_output>0) THEN
    IF (fg_out_col > 0) print*, "Read genoa settings from nml: ", fg_out_col
    IF (LEN_TRIM(err_sps_list_in) > 0) print*, "err_sps_list_in: ", trim(err_sps_list_in)
    IF (LEN_TRIM(init_conc_list_in) > 0) print*, "init_conc_list_in: ", trim(init_conc_list_in)
    IF (LEN_TRIM(ref_conc_list_in) > 0) print*, "ref_conc_list_in: ", trim(ref_conc_list_in)
  ENDIF

END SUBROUTINE rd_nml

!-----------------------------------
! default  parameters for bomxmodel simulations
! => values to use if not supplied by user
!-----------------------------------
SUBROUTINE define_defaults
  USE parameter_mod, ONLY: tstart,tstop,ntstep,nskip,rtol,atol,numit,dtmin,&
                           Rp0,nvoc,nvic,denssoa,densseed,Mseed
  USE mechdata,      ONLY: kdilu, kaerloss
  IMPLICIT NONE

! simulation parameters
    tstart = 0. ; tstop  = 3600.
    ntstep = 12 ; nskip  = 1         
    rtol=0.01   ; atol=100. 
    numit=10    ; dtmin=0.1

    iyear=2000 ; imonth=3 ; iday=21       
    sla=45.    ; slo=0.   ; tz=0.     

    rh0      = 3.
    temp0    = 298.
    kdilu    = 0.      ! dilution rate constant (apply to all species in all phase)
    kaerloss = 0.      ! particle loss rate constant (apply to species in particles only) 

    nvoc     = 0.      ! Non Volatile Organic Compound (concentration)
    nvic     = 0.      ! Non Volatile Inorganic Compound (concentration)
    Rp0      = 1E-5    ! Initial radius of the seed particles 
    denssoa  = 1.06    ! soa density (g.cm-3)    
    densseed = 1.06    ! seed density (g.cm-3)    
    Mseed    = 427.    ! seed molar mass (DOS=427 g/mol)
    
END SUBROUTINE define_defaults


!=======================================================================  
! PURPOSE: Read species name for error computations - genoa
!=======================================================================  
SUBROUTINE read4genoa(chrsp,err_sps_inds,refconc,cinit,fg_init,refname)
  USE parameter_mod, ONLY: fg_ref_conc,fg_out_col,err_sps_list_in, &
                           ref_conc_list_in, init_conc_list_in, fg_output
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)
  INTEGER,INTENT(IN) :: fg_init ! Flag to read the initial concentration
  CHARACTER(LEN=30),INTENT(IN) :: refname ! Name of the reference file
  INTEGER,ALLOCATABLE,INTENT(OUT) :: err_sps_inds(:,:) ! no.group, no.species
  REAL,ALLOCATABLE,INTENT(OUT) :: cinit(:,:)  ! Save the index and initial concentration
  REAL,ALLOCATABLE,INTENT(OUT) :: refconc(:,:,:) ! Real concentration from the file

  INTEGER :: i,j,k,n,m,l,ios,ln
  CHARACTER(LEN=LEN(chrsp(1))) :: sname
  CHARACTER(LEN=80) :: mesg1, mesg2
  CHARACTER(LEN=100) :: filename

!=== Read the species for error computation ===
  IF (fg_out_col.GT.1) THEN
    l=LEN_TRIM(err_sps_list_in) ! Length of the string
    IF (l==0) THEN
      PRINT*, "Error: fg_out_col > 1 but err_sps_list_in is empty"
      STOP
    ENDIF
    IF (fg_output>0) PRINT*, "Read the species for error computation: ",TRIM(err_sps_list_in)
    n=2; m=1 ! Number of group (1 preserved for total soas), number of totoal species
    DO i=1,l
      IF (err_sps_list_in(i:i)==';') THEN
        n=n+1 ! Number of group
        m=m+1 ! Number of total species
      ELSE IF (err_sps_list_in(i:i)==',') THEN
        m=m+1 ! Number of total species
      ENDIF
    ENDDO

    ! Check if the number of group is the same as fg_out_col
    IF (n/=fg_out_col) THEN
      WRITE(*,*) 'Number of group in err_sps_list_in is different from fg_out_col',n,fg_out_col
      STOP
    ENDIF
  
    ! Allocate the array
    ALLOCATE(err_sps_inds(m,2)); err_sps_inds=0

    ! Reread err_sps_list_in
    j=1; n=2; m=0; sname=''
    loop_err_sps: &
    DO i=1,l
      IF (err_sps_list_in(i:i)==';'.OR.err_sps_list_in(i:i)==',') THEN
        sname=TRIM(ADJUSTL(err_sps_list_in(j:i-1)))
      ELSE IF (i==l) THEN
        sname=TRIM(ADJUSTL(err_sps_list_in(j:i)))
      ENDIF
      IF (LEN_TRIM(sname)>0) THEN ! Read the species
        m=m+1 ! Number of species
        j=i+1 ! Start of the next species
        DO k=1,SIZE(chrsp) ! Find the index of the species in chrsp
          !print*, trim(sname), trim(chrsp(k))
          IF (sname==chrsp(k)) THEN
            err_sps_inds(m,1)=n
            err_sps_inds(m,2)=k
            ! Print to check
            IF (fg_output>0) PRINT*, "Read species ",TRIM(sname)," in group ",n &
                                    ," at index ",k," in chrsp ",k
            IF (err_sps_list_in(i:i)==';') n=n+1 ! Update number of group
            sname='' ! Reset species name
            CYCLE loop_err_sps
          ENDIF
        ENDDO
        ! Stop if the species is not found in chrsp
        PRINT*, 'Error Species ' // TRIM(sname) // ' not found in the species list'
        STOP
      ENDIF
    ENDDO loop_err_sps
  ENDIF

!=== Read the reference concentration ===
  IF (LEN_TRIM(ref_conc_list_in) .GT. 0) THEN
    ! Update fg_ref_conc by "," in ref_conc_list_in
    fg_ref_conc=1; ln=0 ! Number of lines(timesteps) in the file
    DO i=1, LEN_TRIM(ref_conc_list_in)
      IF (ref_conc_list_in(i:i)==',') fg_ref_conc=fg_ref_conc+1
    ENDDO
    IF (fg_output>0) PRINT*, "Read reference concs from ", fg_ref_conc, &
            " files in ",TRIM(ADJUSTL(ref_conc_list_in))
    
    ! Read concs from files
    k=0; j=1; m=1 ! End of the file name, Start of the file name, No. of files
    DO i=1, LEN_TRIM(ref_conc_list_in)
      ! Find the end of the file name
      IF (ref_conc_list_in(i:i)==',') THEN
        k=i-1
      ELSE IF (i==LEN_TRIM(ref_conc_list_in)) THEN
        k=i
      ENDIF
      ! Read the file
      IF (k/=0) THEN
        IF (j>=k) THEN
          PRINT*, 'Error: Empty file name in ref_conc_list_in'
          STOP
        ENDIF
        ! Open the file
        filename=TRIM(ADJUSTL(ref_conc_list_in(j:k)))
        IF (TRIM(refname)/='outdat') THEN
          l=LEN_TRIM(refname)  ! Add '/' to the directory if needed
          IF (filename(l:l)/='/') filename=TRIM(filename)//'/'
          filename=TRIM(filename)//TRIM(refname)//'.genoa'
        ENDIF
        IF (fg_output>0) PRINT*,"Reading file: ",TRIM(filename)," ..."
        OPEN(UNIT=30+m,FILE=TRIM(filename),STATUS='OLD',IOSTAT=ios)
        IF (ios/=0) THEN
          PRINT*, 'Error: Ref conc File ',TRIM(filename),' not found'
          STOP
        ENDIF
        ! Read no. lines(timesteps) in the file
        IF (m==1) THEN
          DO
            READ(30+m,*,IOSTAT=ios)
            IF (ios/=0) EXIT
            ln=ln+1
          ENDDO
          ! Allocate the array
          ALLOCATE(refconc(fg_ref_conc,fg_out_col,ln)); refconc=0.
          REWIND(30+m)
          ! Print to check
          IF (fg_output>0) PRINT*, "Read ",ln," lines from reference files"
        ENDIF

        ! Read the concentration
        DO n=1, ln
          READ(30+m,*,IOSTAT=ios) (refconc(m,j,n),j=1,fg_out_col)
          IF (ios/=0) THEN
            PRINT*, 'Error: Can not read the concentration from line',n &
                    ,' in file ',TRIM(filename), &
                    'Assuming file contains ',ln,' lines.'
            STOP
          ENDIF
        ENDDO

        ! Close the file
        CLOSE(30+m)
        m=m+1 ! Update the counter
        j=i+1 ! Start of the next file
        k=0   ! Reset the end of the file name
        IF (fg_output>0) PRINT*, "Finished reading ",TRIM(filename)
      ENDIF
    ENDDO
    ! check m
    IF (m-1/=fg_ref_conc) THEN
      PRINT*, 'Error: No. of files read ',m-1,' not equal to fg_ref_conc ',fg_ref_conc
      STOP
    ENDIF
  ENDIF

!=== Read the initial concentrations ===
  IF (fg_init > 0) THEN
    l=LEN_TRIM(init_conc_list_in) ! Length of the string
    IF (l==0) THEN
      PRINT*, "Error: fg_init > 0 but init_conc_list_in is empty"
      STOP
    ENDIF
    ! Get the part contains info
    n=1; m=1; j=1; k=1 ! No.species in the group, No.group, group start/end index
    rd_init: &
    DO i=1,l
      IF (init_conc_list_in(i:i)==',') THEN
        n=n+1 ! Read No.species in the group
      ELSE IF (init_conc_list_in(i:i)==';') THEN ! Find new group
        IF (m==fg_init) THEN
          k=i-1 ! End of the group
          EXIT rd_init
        ENDIF
        m=m+1
        n=1 ! Reset No.species
        j=i+1 ! Start of the next group
      ENDIF
    ENDDO rd_init
    ! Check if found fg_init
    IF (m < fg_init) THEN
      PRINT*, "Error: fg_init > no.groups in init_conc_list_in: ",fg_init,m," from: ",trim(init_conc_list_in)
      STOP
    ENDIF
    IF (k<=j) k=l   ! No ";"
    ! Allocate the conc array
    ALLOCATE(cinit(n,2)); cinit=0.
    ! Reread and find species with conc. in init_conc_list_in(j:k)
    sname=''; mesg1=''; l=1
    rd_init2: &
    DO i=j, k
      IF (init_conc_list_in(i:i)==',') THEN
        mesg1=TRIM(ADJUSTL(init_conc_list_in(j:i-1)))
      ELSE IF (i==k) THEN
        mesg1=TRIM(ADJUSTL(init_conc_list_in(j:i)))
      ENDIF
      IF (LEN_TRIM(mesg1)==0) CYCLE rd_init2
      ! Read species and conc. from mesg1
      m=INDEX(mesg1,' ') ! Find 1st space
      IF (m==0) THEN
        PRINT*, 'Error: No space in the string: ',mesg1
        STOP
      ENDIF
      READ(mesg1(1:m-1),*,IOSTAT=ios) sname
      IF (ios/=0) THEN
        PRINT*, 'Error: Can not read species name: ',TRIM(mesg1(1:m-1))," from ",trim(init_conc_list_in(j:k))
        STOP
      ENDIF
      sname=TRIM(ADJUSTL(sname))
      DO n=1,SIZE(chrsp) ! Find the index of the species in chrsp
        IF (sname==chrsp(n)) THEN
          cinit(l,1)=n ! Save the index in chrsp
          READ(mesg1(m+1:),*,IOSTAT=ios) cinit(l,2) ! Read the concentration
          IF (ios/=0) THEN
            PRINT*, 'species: ',TRIM(sname)
            PRINT*, 'Error: Can not read conc.: ',TRIM(mesg1(m+1:))," from ",trim(init_conc_list_in(j:k))
            STOP
          ENDIF
          l=l+1 ! Update counter
          IF (fg_output>0) PRINT*, "Read species ",TRIM(sname)," at index ",n," with conc. ",cinit(l-1,2)
          CYCLE rd_init2
        ENDIF
      ENDDO
      PRINT*, 'Error: Species ',TRIM(sname),' not found in the species list'
      STOP
    ENDDO rd_init2
  ENDIF
END SUBROUTINE read4genoa

END MODULE initnbound
