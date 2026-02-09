MODULE ro2maptool

! structure for storing RO2 species in the corresponding classes 
  TYPE :: clasro2
    INTEGER             :: nro2          ! number of RO2 in the class
    INTEGER,ALLOCATABLE :: idro2(:)      ! ID of the RO2 species
  END TYPE
  TYPE(clasro2),ALLOCATABLE :: ro2map(:) ! RO2 species map (size: # of pero class)

CONTAINS  
!=======================================================================
! PURPOSE: RO2+RO2 reactions use "counters", i.e. the sum of the 
! concentrations of RO2 species belonging to a given "class". This 
! routine read the files providing the list of species in each class,
! allocate the required memory for each class and store the species 
! that belong to each class.
!=======================================================================
SUBROUTINE readro2(lout,chrsp,ncpero)
  USE boxtool, ONLY: stoperr
  USE parameter_mod, ONLY: input_dir_chem, fg_output
  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: lout      ! file unit for information outputs
  CHARACTER(LEN=*), INTENT(IN) :: chrsp(:)  ! list (names) of the species
  INTEGER, INTENT(IN)          :: ncpero    ! # of peroxy (RO2) classes

  CHARACTER(LEN=LEN(chrsp(1))),ALLOCATABLE :: ro2sp(:)
  CHARACTER(LEN=LEN(chrsp(1))) :: sname
  CHARACTER(LEN=50) :: liner
  CHARACTER(LEN=200) :: filnam
  INTEGER :: i,j,nclas,nclas0,ispe,ios,iloc,nspeall
  INTEGER :: nspe(ncpero)
  LOGICAL :: loerr
  INTEGER, ALLOCATABLE :: ro2id(:), flg(:) ! RO2 class
  
  INTEGER,PARAMETER :: filu=12                          ! unit of the files to be read 
  CHARACTER(LEN=10),PARAMETER :: progname='readro2'
  CHARACTER(LEN=80) :: mesg1, mesg2

! Minor change required in this routine if ncpero > 9 (see open below)
  IF (ncpero  > 9) THEN
    mesg1="# of RO2 classes > 9."
    mesg2="Change the program to read more classes." 
    CALL stoperr(progname,mesg1,mesg2)
  ENDIF

  !WRITE(lout,*) ' '
  !WRITE(lout,*) '---- RO2+RO2 reactions -----'

! LOOP OVER THE RO2 CLASS
! -----------------------

  IF (fg_output>0) PRINT*, "   read RO2 counter file: indat.ro2"

! open file to be read
! --------------------
  i=LEN_TRIM(input_dir_chem)
  IF (input_dir_chem(i:i) == '/') THEN
    filnam=TRIM(input_dir_chem)//'indat.ro2'
  ELSE
    filnam=TRIM(input_dir_chem)//'.ro2'
  ENDIF
  OPEN(filu, FILE=filnam, STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
    mesg1="file not found/open: "//TRIM(filnam)
    CALL stoperr(progname,mesg1)
  ENDIF

! read RO2 species that belong to the class
! -----------------------------------------
! Read first line: numbers for each class

  READ(filu,*,IOSTAT=ios) (nspe(i), i=1, ncpero)
  IF (ios/=0) THEN 
    mesg1="cannot read the list of # species in file: "//TRIM(filnam)
    CALL stoperr(progname,mesg1)
  ENDIF
   ! allocate required memory to the ro2map
  DO nclas=1,ncpero
      ro2map(nclas)%nro2 = nspe(nclas)
      ALLOCATE(ro2map(nclas)%idro2(nspe(nclas)+1))  ! add 1 slot (just in case)
      ro2map(nclas)%idro2(:)=0
  ENDDO
  nspeall = sum(nspe)      ! total number of ro2 species
  ALLOCATE(ro2sp(nspeall))     ! not used out of this subroutine
  ALLOCATE(ro2id(nspeall))     ! not used out of this subroutine
  ALLOCATE(flg(nspeall))       ! not used out of this subroutine
   ! loop over the species in the file
  ispe=0
  speloop: & 
  DO                  
    READ(filu, '(a)', IOSTAT=ios) liner ! Read species name and ro2 class
    IF (ios/=0) THEN 
      mesg1="nothing left to read in file: "//TRIM(filnam)
      mesg2="keyword 'END' missing?" 
      CALL stoperr(progname,mesg1,mesg2)
    ENDIF
     ! escape the loop if keyword 'END' is found
    IF (liner(1:3)=='END') THEN   
      IF (ispe/=nspeall) THEN
        mesg1="number of species mismatches expected record in: "//TRIM(filnam)
        CALL stoperr(progname,mesg1)
      ENDIF
      EXIT speloop
    ENDIF 
     ! check memory (just in case)
    ispe=ispe+1
    IF (ispe > nspeall) THEN
      mesg1="# of RO2 species > expected record in: "//TRIM(filnam)
      CALL stoperr(progname,mesg1)
    ENDIF
     ! check empty line
    IF (LEN_TRIM(liner) == 0) THEN
      mesg1="Unexpected empty line read in: "//TRIM(filnam)
      CALL stoperr(progname,mesg1)
    ENDIF
    ! store (temporarily) the species
    READ(liner, *, IOSTAT=ios) sname, nclas
    ! check before storing
    if (ios/=0) THEN
      mesg1="cannot read the species from line: "//TRIM(liner)
      CALL stoperr(progname,mesg1)
      EXIT speloop
    ENDIF
    ro2sp(ispe)=adjustl(sname)
    ro2id(ispe)=nclas
    !print*, "READ: ", trim(ro2sp(ispe)), "  ",ro2id(ispe), ispe
  ENDDO speloop
  CLOSE(filu)

! get the ID for the RO2 species & store ID in idchemro2
! ------------------------------------------------------
    loerr=.FALSE. ;  iloc=1 ; nspe=0; flg=0; nclas0=1

! fast seek - require same order in tables (OK in gecko)    
    seekloop1:&    
    DO i=1,nspeall
      sname=ro2sp(i) ! RO2 species name 
      nclas=ro2id(i) ! RO2 class index
      IF (nclas0 .ne. nclas) THEN
        iloc = 1       ! Change of group - start from beginning
        nclas0 = nclas ! Reset
      ENDIF
      DO j=iloc,SIZE(chrsp)
        IF (sname==chrsp(j)) THEN

          ! check for duplicate
          IF (flg(i) .gt. 0) THEN          
            mesg1="The following species was provided twice: "//TRIM(sname)
            mesg2="last found in file: "//TRIM(filnam)
            CALL stoperr(progname,mesg1,mesg2)
          ENDIF

          ! store ID
          nspe(nclas) = nspe(nclas) + 1 ! Count current species number in RO2 class
          ro2map(nclas)%idro2(nspe(nclas))=j
          flg(i)=j ; iloc=j+1
          CYCLE seekloop1
        ENDIF
      ENDDO

      ! if that point is reached then species not found => try unsorted loop
      loerr=.TRUE. 
      IF (fg_output>0) PRINT*, TRIM(sname), " not found in fast seek."
    ENDDO seekloop1

! if fast seek did not succeed, then try slow seek (no order required - OK for MCM)
    IF (loerr) THEN  
      seekloop2:&  
      DO i=1,nspeall
        if (flg(i) > 0) CYCLE seekloop2 ! Already found
        sname=ro2sp(i) ! RO2 species name 
        nclas=ro2id(i) ! RO2 class index
        IF (fg_output>0) PRINT*, "   => try slow seek for: ro2 species ", TRIM(sname)
        DO j=1,SIZE(chrsp)
          IF (sname==chrsp(j)) THEN
              ! store ID
              nspe(nclas) = nspe(nclas) + 1 ! Count current species number in RO2 class
              ro2map(nclas)%idro2(nspe(nclas))=j
              flg(i) = j
              CYCLE seekloop2
          ENDIF
        ENDDO
        ! if that point is reached then species not found
        mesg1="The following species not identified: "//TRIM(sname)
        mesg2="but provided in file: "//TRIM(filnam)
        CALL stoperr(progname,mesg1,mesg2)
      ENDDO seekloop2
    ENDIF
    DEALLOCATE(ro2sp)
    DEALLOCATE(ro2id)
    DEALLOCATE(flg)
      
END SUBROUTINE readro2

END MODULE ro2maptool
