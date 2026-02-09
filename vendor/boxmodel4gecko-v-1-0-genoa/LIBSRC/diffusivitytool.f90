MODULE diffusivitytool
IMPLICIT NONE

INTEGER             :: ndv         ! # of partitioning species
REAL,ALLOCATABLE    :: dvol(:)     ! [i] diffusion volume of species i

CONTAINS 
! SUBROUTINE readdv(chrsp,numsp,lodicnam)

!=======================================================================
! Purpose: change diffusion volume into the required numbers to compute 
! diffusivity. After the call, dvol(i) is transformed into:
!     0.00143/[ [(M_AB)^0.5 x [dvol(A)^(1/3)+dvol_B^(1/3)]^2]
! where A stand for the species i and B is air. This numbers are 
! next used in the Fuller SAR, as described in: "The Properties of 
! Gases and Liquids", 5th Edition, B.E. Poling, J.M. Prausnitz, J.P. 
! O’Connell, The McGraw-Hill, 2004. See equation 11-4.4 in that book 
!=======================================================================
SUBROUTINE chge_dvol(idasat,wmol)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: idasat(:)         ! ID of partitioning species - part. phase
  REAL,INTENT(IN)    :: wmol(:)           ! molar mass of the ith species

  INTEGER            :: i
  REAL               :: d_air,wfact

  CHARACTER(LEN=15),PARAMETER :: progname='chge_dvol'

  ! diffusion volume for air ("the properties of gases and liquids", 
  ! Poling et al., 2005, Table 11.1) 
  REAL,PARAMETER :: air_dvol=19.7 

! store diffusion volume as a pre-computed parameter, i.e. :
! 0.00143/[ [(M_AB)**0.5 * [dvol**(1/3)+dvol_air(1/3)]**2]
! See equation (11-4.4) in Poling, 2005, book (see above))
! Note : M_AB=2*(1/M_A+1/M_B), as described in Poling
  d_air=air_dvol**(1./3.)
  dvol(:)=dvol(:)**(1./3.) 
  dvol(:)=(dvol(:)+d_air)**2.
  DO i=1,ndv
    wfact=1./(wmol(idasat(i)))+1./29.     ! M in g/mol 
    wfact=SQRT(2./wfact)        
    dvol(i)=dvol(i)*wfact
    dvol(i)=0.00143/dvol(i)
  ENDDO
END SUBROUTINE chge_dvol

!=======================================================================
! Purpose: read diffusion volumes.
! BA : another similar subroutine is readTg. Consider merging both
! routines to improve code efficiency... 
!=======================================================================
SUBROUTINE readdv(chrsp,numsp,namsat,idgsat,idasat,wmol,lodicnam)
  USE boxtool, ONLY: stoperr
  USE parameter_mod, ONLY: mxlsp,input_dir_chem
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: numsp             ! # of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:) ! name of the species
  CHARACTER(LEN=*),INTENT(IN):: namsat(:) ! name of the partitioning species (with 'G')
  INTEGER,INTENT(IN) :: idgsat(:)         ! ID of the partitioning species
  INTEGER,INTENT(IN) :: idasat(:)         ! ID of partitioning species - part. phase
  REAL,INTENT(IN)    :: wmol(:)           ! molar mass of the ith species
  LOGICAL,INTENT(IN),OPTIONAL :: lodicnam ! if true and present then chrsp is dictionary names

  CHARACTER(LEN=mxlsp) :: namdv           ! name of the partitioning species
  CHARACTER(LEN=160) :: line
  INTEGER            :: i,j,iloc,ios,filu,ibeg,nrec,idv
  LOGICAL            :: loerr

  CHARACTER(LEN=15),PARAMETER :: progname='readdv'
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
    OPEN (filu, FILE=TRIM(input_dir_chem)//'mdv.dat', STATUS='OLD')
  ELSE
    OPEN (filu, FILE=TRIM(input_dir_chem)//'.mdv', STATUS='OLD')
  ENDIF

! read the number of record 
  READ (filu,*, IOSTAT=ios) nrec
  IF (ios/=0) THEN 
    mesg1="# record cannot be read - keyword END missing?"
    CALL stoperr(progname,mesg1)
  ENDIF  
  IF (nrec/=ndv) THEN 
    mesg1="# record does not match size used to allocate table memory"
    CALL stoperr(progname,mesg1)
  ENDIF  

! read the file
! -------------
  idv=0
  DO     
    READ (filu, '(a)', IOSTAT=ios) line
    IF (line(1:1)=='!') CYCLE

    ! exit if keyword END found
    IF (ios/=0) THEN 
      mesg1="keyword END missing in the input file?"
      CALL stoperr(progname,mesg1)
    ENDIF  
    IF (line(1:3) == 'END') EXIT        
    
    ! get next input 
    idv=idv+1                  
    READ(line,*,IOSTAT=ios) namdv,dvol(idv)
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
    
    ! check that sorting in mdv and Pvap tables are identical
    IF (namdv/=namsat(idv)) THEN
      mesg1="species in namsat (pvap) and namdv (mdv) are not identical"
      mesg2="same order is required. Mismatch at species: "//TRIM(namdv)
      CALL stoperr(progname,mesg1,mesg2)
    ENDIF  
  ENDDO
  CLOSE(filu)

  IF (idv/=nrec) THEN
    mesg1="number of record /= number of data"
    CALL stoperr(progname,mesg1)
  ENDIF
  IF (idv > SIZE(dvol)) THEN
    mesg1="idv > SIZE(dvol) - The size of the tables is too small"
    CALL stoperr(progname,mesg1)
  ENDIF

! check the ID for the species in namsat 
! --------------------------------------
  loerr=.FALSE.
  iloc=1

! fast seek - require same order in tables (OK in gecko)    
  seekloop1:&    
  DO i=1,ndv
    DO j=iloc,numsp
      IF (namsat(i)(ibeg:) == chrsp(j)) THEN
        iloc=j
        ! check the consistency with the Pvap data (must be identical)
        IF (j/=idgsat(i)) THEN
          mesg1="mismatch in the ID between dvol and Pvap data"
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
    WRITE(6,*) "       => try slow seek in readdv"
    seekloop2: &   
    DO i=1,ndv
      DO j=1,numsp
        IF (namsat(i)(ibeg:) == chrsp(j)) THEN
          ! check the consistency with the Pvap data (must be identical)
          IF (j/=idgsat(i)) THEN
            mesg1="mismatch in the ID between dvol and Pvap data"
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

END SUBROUTINE readdv

END MODULE diffusivitytool
