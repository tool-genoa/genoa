MODULE speciestool
IMPLICIT NONE
CONTAINS
! SUBROUTINE ptrspecies(lout,line,numsp,lopen,ptrgas,ptrpart,ptrwall)          
!***********************************************************************
! PURPOSE: set pointer in the list of species to identify contributors 
! to the composition of each phase (or block of species).
! The program stop whenever a problem is found. 
!***********************************************************************
SUBROUTINE ptrspecies(lout,line,numsp,lopen,ptrgas,ptrpart,ptrwall)          
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: lout           ! unit output file
  CHARACTER(LEN=*),INTENT(IN) :: line     ! line to be checked
  INTEGER,INTENT(IN)    :: numsp          ! number of species
  INTEGER,INTENT(INOUT) :: ptrgas(:)      ! pointer (1st and last) gas phase species
  INTEGER,INTENT(INOUT) :: ptrpart(:)     ! pointer (1st and last) part. phase species
  INTEGER,INTENT(INOUT) :: ptrwall(:)     ! pointer (1st and last) wall phase species
  LOGICAL,INTENT(INOUT) :: lopen          ! true is phase "open"

  IF (line(1:16)=='PHASE: START GAS') THEN
    IF (lopen) STOP 'start gas phase while another phase still not closed'
    IF (ptrgas(1)/=0) STOP 'keyword "start gas" used more than once'
    ptrgas(1)=numsp+1 ; lopen=.TRUE. ; RETURN
  ENDIF
  IF (line(1:14)=='PHASE: END GAS') THEN
    IF (.not.lopen) STOP 'close gas phase while no phase open'
    IF (ptrgas(2)/=0) STOP 'keyword "end gas" used more than once'
    ptrgas(2)=numsp ; lopen=.FALSE. ; RETURN
  ENDIF

  IF (line(1:18)=='PHASE: START PART.') THEN
    IF (lopen) STOP 'start part. phase while another phase still not closed'
    IF (ptrpart(1)/=0) STOP 'keyword "start part" used more than once'
    ptrpart(1)=numsp+1 ; lopen=.TRUE. ; RETURN
  ENDIF
  IF (line(1:16)=='PHASE: END PART.') THEN
    IF (.not.lopen) STOP 'close part. phase while no phase open'
    IF (ptrpart(2)/=0) STOP 'keyword "end part" used more than once'
    ptrpart(2)=numsp ; lopen=.FALSE. ; RETURN
  ENDIF

  IF (line(1:17)=='PHASE: START WALL') THEN
    IF (lopen) STOP 'start wall phase while another phase still not closed'
    IF (ptrwall(1)/=0) STOP 'keyword "start wall" used more than once'
    ptrwall(1)=numsp+1 ; lopen=.TRUE. ; RETURN
  ENDIF
  IF (line(1:15)=='PHASE: END WALL') THEN
    IF (.not.lopen) STOP 'close wall phase while no phase open'
    IF (ptrwall(2)/=0) STOP 'keyword "end wall" used more than once'
    ptrwall(2)=numsp ; lopen=.FALSE. ; RETURN
  ENDIF

END SUBROUTINE ptrspecies

!***********************************************************************
! PURPOSE: Given an input line of the form : 
!     species_name  /weight/  
! add the species in the chrsp table and the corresponding molecular
! weight in the wmol table. The line must contain one species only.
! Molecular weight is optional here (set to 0 if not given). 
!***********************************************************************
SUBROUTINE readspecies(lout,line,numsp,chrsp,wmol,lostop)          
  USE intptool, ONLY: stringarg
  IMPLICIT NONE

  INTEGER,INTENT(IN)    :: lout            ! unit output file
  CHARACTER(LEN=*),INTENT(INOUT) :: line   ! line to be read
  INTEGER,INTENT(INOUT) :: numsp           ! number of species
  CHARACTER(LEN=*),INTENT(OUT) :: chrsp(:) ! list of species
  REAL,INTENT(OUT)      :: wmol(:)         ! molar mass of the species
  LOGICAL,INTENT(OUT)   :: lostop          ! error flag

  INTEGER,PARAMETER :: mxarg=10            ! max # of arguments in a line
  INTEGER :: posi(mxarg,2)
  INTEGER :: narg
  INTEGER :: ios, lspe, mlspe, ipos

! max length of the species
  mlspe=LEN(chrsp(1))

! ------------------
! READ THE LINE 
! ------------------
! Molecular weight is optional in this routine. Slashes is optional 
! Warning : only 1 species per line!
  
! remove the slashes (if any)
  ipos=-1
  DO WHILE (ipos /= 0)  
    ipos=INDEX(line,'/') 
    IF (ipos /= 0) line(ipos:ipos)=' '           ! rm '/'
  ENDDO

! get positions of the arguments
  CALL stringarg(line,narg,posi)
  IF (narg == 0) THEN
    WRITE(lout,*)'--error--,  unexpected empty line in SPECIES data '
    lostop=.TRUE. ; RETURN
  ENDIF
  IF (narg > 2) THEN
    WRITE(lout,*) TRIM(line) 
    WRITE(lout,*)'--error--,  only 2 arguments expected'
    lostop=.TRUE. ; RETURN
  ENDIF

! raise number of species - Stop if exceed the size of the table
  numsp=numsp+1
  IF (numsp > SIZE(chrsp,1)) THEN
    WRITE(lout,*) '--error--, numsp > max species, i.e.', SIZE(chrsp,1)
    STOP 'in readspecies - too many species (see out file)' 
  ENDIF

! add new species name
  lspe=posi(1,2)-posi(1,1) + 1
  IF (lspe > mlspe) THEN
    WRITE(lout,*) TRIM(line) 
    WRITE(lout,*)'--error--, length of species name > mlspe'
    WRITE(lout,*)'           mslpe =', mlspe
    lostop=.TRUE. ; RETURN
  ENDIF
  chrsp(numsp)=line(posi(1,1):posi(1,2))

! add molecular weight
  wmol(numsp)=0.  ! default value
  IF (narg==2) THEN
    READ(line(posi(2,1):posi(2,2)),*,IOSTAT=ios)  wmol(numsp)
    IF (ios/=0) THEN
      WRITE(lout,*) TRIM(line) 
      WRITE(lout,*)'--error--, while reading MW'
      lostop=.TRUE.
    ENDIF
  ENDIF

END SUBROUTINE readspecies

!***********************************************************************
! PURPOSE: Check that the species given as input are correctly set. 
! Outputs the number of species in each phase. 
! "lostop" flag is raised if an error is found. 
!***********************************************************************
SUBROUTINE chkspecies(lout,lospeech,numsp,chrsp,sortsp,isortlk,wmol, &
                  ptrgas,ptrpart,ptrwall, &
                  ngas,npart,nwall,lostop)
  USE keywordlist, ONLY: mxkeywd,keywdlist
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: lout        ! unit output file
  LOGICAL,INTENT(IN)  :: lospeech    ! print info on screen
  INTEGER,INTENT(IN)  :: numsp       ! number of species
  CHARACTER(lEN=*),INTENT(IN) ::  chrsp(:) ! list of species
  CHARACTER(lEN=*),INTENT(IN) :: sortsp(:) ! sorted list of species
  INTEGER,INTENT(IN)  :: isortlk(:)  ! index link "sorted" -> "mechanism" species
  REAL,INTENT(IN)     :: wmol(:)     ! molar mass of the species
  INTEGER,INTENT(IN)  :: ptrgas(:)   ! pointer (1st and last) gas phase species
  INTEGER,INTENT(INOUT) :: ptrpart(:)  ! pointer (1st and last) part. phase species
  INTEGER,INTENT(INOUT) :: ptrwall(:)  ! pointer (1st and last) wall phase species
  INTEGER,INTENT(OUT) :: ngas        ! # of gas phase species
  INTEGER,INTENT(OUT) :: npart       ! # of part. phase species
  INTEGER,INTENT(OUT) :: nwall       ! # of wall phase species
  LOGICAL,INTENT(OUT) :: lostop      ! error flag

  INTEGER :: i, j, mlspe, nall
  LOGICAL :: lochar

  mlspe=LEN(chrsp(1)) ! maximum length of the species
  
! ------------------------------------------------------------
! Check if the name of the species contains illegal character
! ------------------------------------------------------------
  IF (lospeech) PRINT*, '.... check characters in species'
  
  speloop: DO  i=1,numsp

! check the first character: must start with Aa-Zz or '(<[{'        
    lochar=.TRUE.      
    IF      (chrsp(i)(1:1)>='A' .AND. chrsp(i)(1:1)<='Z') THEN ; lochar=.FALSE.
    ELSE IF (chrsp(i)(1:1)>='a' .AND. chrsp(i)(1:1)<='z') THEN ; lochar=.FALSE.
    ELSE IF (chrsp(i)(1:1)=='(') THEN ; lochar=.FALSE.
    ELSE IF (chrsp(i)(1:1)=='<') THEN ; lochar=.FALSE.
    ELSE IF (chrsp(i)(1:1)=='[') THEN ; lochar=.FALSE.
    ELSE IF (chrsp(i)(1:1)=='{') THEN ; lochar=.FALSE.
    ENDIF
    IF (lochar) THEN 
      WRITE(lout,*) '--error--, illegal 1st character in: ', chrsp(i)
       lostop=.TRUE.
    ENDIF
  
! check the other characters: must not contain '+-='         
    j=SCAN(chrsp(i),'=+-')
    IF (j /= 0) THEN 
      WRITE(lout,*) '--error--, character not allowed in:', chrsp(i)
      lostop=.TRUE.
    ENDIF
  ENDDO speloop

! ------------------------------------------------------------
! check for forbidden names (keywords)
! ------------------------------------------------------------
  IF (lospeech) PRINT*, '.... search forbidden names'
  DO i=1,numsp
    lochar=.FALSE.
    DO j=1,mxkeywd  
      IF (chrsp(i)==keywdlist(j)) lochar=.TRUE.
    ENDDO
    IF (chrsp(i)=='NOTHING') lochar=.TRUE. ! keywd not in keywdlist 
    IF (lochar) THEN
      WRITE(lout,*)'--error--, species is also a keyword: ',chrsp(i)
      lostop=.TRUE.
    ENDIF
  ENDDO

! ------------------------------------------------------------
! Check for duplicate
! ------------------------------------------------------------
! Sorted list so duplicates must be neighbours
  IF (lospeech) PRINT*, '.... search for duplicate species'
  DO i=1,numsp-1
    IF (sortsp(i)==sortsp(i+1)) THEN
      WRITE(lout,*)'--error--, duplicate species: ', sortsp(i)
      lostop=.TRUE.
    ENDIF
  ENDDO

! -------------------------------------------------------
! Check molecular weight  
! -------------------------------------------------------
  IF (lospeech) PRINT*, '.... check molecular weight'
  DO i=1,numsp
    IF (wmol(i) <= 0.0) THEN
      WRITE(lout,*)'--error--, no molecular weight set for: ', chrsp(i)
      lostop=.TRUE.
    ENDIF 
  ENDDO 

! -------------------------------------------------------
! Check molecular weight  
! -------------------------------------------------------
  IF (lospeech) PRINT*, '.... check number of species in each phase ...'

! check that closed phase were first open
  IF (ptrgas(1)==0) THEN;
    IF (ptrgas(2)/=0) STOP 'in chkspecies, gas phase closed but not open'
  ENDIF
  IF (ptrpart(1)==0) THEN;
    IF (ptrpart(2)/=0) STOP 'in chkspecies, part. phase closed but not open'
  ENDIF
  IF (ptrwall(1)==0) THEN;
    IF (ptrwall(2)/=0) STOP 'in chkspecies, wall phase closed but not open'
  ENDIF
   
! check consistency of the pointer
  ngas=0 ; npart=0 ; nwall=0

  IF (ptrgas(1)>0) THEN
    IF (ptrgas(2)<ptrgas(1)) THEN
      WRITE(lout,*)'--error--, # of "gas phase" species <= 0'
      lostop=.TRUE.
    ELSE  
      ngas=ptrgas(2)-ptrgas(1)+1
    ENDIF
  ENDIF

  IF (ptrpart(1)>0) THEN
    IF (ptrpart(2)<ptrpart(1)) THEN
      IF (ptrpart(1)==ptrpart(2)+1) THEN ! pointer found but empty list of species
        npart=0
        ptrpart(1)=0
        ptrpart(2)=0
      ELSE
        WRITE(lout,*)'--error--, # of "particulate phase" species < 0'
        lostop=.TRUE.
      ENDIF
    ELSE  
      npart=ptrpart(2)-ptrpart(1)+1
    ENDIF
  ENDIF

  IF (ptrwall(1)>0) THEN
    IF (ptrwall(2)<ptrwall(1)) THEN
      IF (ptrwall(1)==ptrwall(2)+1) THEN ! pointer found but empty list of species
        nwall=0
        ptrwall(1)=0
        ptrwall(2)=0
      ELSE
        WRITE(lout,*)'--error--, # of "wall phase" species < 0'
        lostop=.TRUE.
      ENDIF
    ELSE  
      nwall=ptrwall(2)-ptrwall(1)+1
    ENDIF
  ENDIF

  nall=ngas+npart+nwall
  IF (nall /= numsp) THEN
    WRITE(lout,*)'--error--, num of species is each phase /= numsp'
    lostop=.TRUE.
  ENDIF 
                
! -------------------------------------------------------
! EXIT
! -------------------------------------------------------
  IF (lospeech) PRINT*,'EXIT chkspecies ...'

END SUBROUTINE chkspecies

END MODULE speciestool
