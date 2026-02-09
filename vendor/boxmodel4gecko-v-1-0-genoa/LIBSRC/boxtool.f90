MODULE boxtool
  IMPLICIT NONE
  CONTAINS
!-----------------------------------------------------------------------
! Purpose : return the ID number of the species (chem) provided as 
! input. Two arguments are optional to handle error :
! - csub: the name of the calling subroutine (to be printed if provided
! - loerr: stop the program if set to true (also the default value) if
!          not provided as argument
! ----------------------------------------------------------------------    
FUNCTION spe2id(lout,chrsp,numsp,chem,csub,loerr)
  IMPLICIT NONE

! Input/output
  INTEGER, INTENT(IN)         :: lout
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)
  INTEGER, INTENT(IN)         :: numsp
  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: csub    ! name of the calling subroutine 
  LOGICAL,INTENT(IN),OPTIONAL :: loerr            ! if true then stop if the function return 0 
  INTEGER                     :: spe2id

! Local
  CHARACTER(LEN=LEN(chrsp(1))):: ichrsp
  INTEGER :: isp,lensp
  LOGICAL :: lostop

! default : stop if the species is not found  
  IF (PRESENT(loerr)) THEN; lostop=loerr ; ELSE ; lostop=.TRUE. ; ENDIF 

  ichrsp=' ' ; spe2id=0 ; lensp=LEN(chrsp(1))
  
  IF (chem(1:1)==' ') THEN
    WRITE(lout,*) '--error-- in spe2id. 1st character is " "'
    WRITE(lout,*) '          for species:',TRIM(chem)
    IF (PRESENT(csub)) WRITE(lout,*) ' => calling subroutine was: ',csub
    IF (lostop) STOP 'in spe2id => species start with " " character' 
  ENDIF
  IF (LEN_TRIM(chem) > lensp) THEN
    WRITE(lout,*) '--error-- in spe2id. Species length is too long'
    WRITE(lout,*) '          for species:',TRIM(chem)
    WRITE(lout,*) '          max length is:',lensp
    IF (PRESENT(csub)) WRITE(lout,*) ' => calling subroutine was: ',csub
    IF (lostop) STOP 'in spe2id => species too long' 
  ENDIF

  ichrsp=chem
  DO isp=1,numsp
    IF (chrsp(isp) == ichrsp) THEN
      spe2id=isp
      RETURN
    ENDIF
  ENDDO
  
  IF (spe2id==0) THEN     
    WRITE(lout,*) '--error-- in spe2id. '
    WRITE(lout,*) ' => species unindentified:',ichrsp
    IF (PRESENT(csub)) WRITE(lout,*) ' => calling subroutine was: ',csub
    IF (lostop) STOP 'in spe2id => species not found' 
  ENDIF

END FUNCTION spe2id

!=======================================================================
! PURPOSE: return the number of record in a file
!=======================================================================
FUNCTION getnrec(lout,filename)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: lout 
  INTEGER :: getnrec
  CHARACTER(LEN=*) :: filename

  INTEGER :: filu,ios

  getnrec=0
  filu=12
  OPEN (filu, FILE=filename, STATUS='OLD')
  READ (filu,*,IOSTAT=ios) getnrec
  IF (ios/=0) THEN 
    WRITE(lout,*) '--error--, while reading number of record'
    WRITE(lout,*) '         , in file:',TRIM(filename)
    WRITE(lout,*) '           in function getnrec'
    STOP 'Stop error in getnrec'
  ENDIF
  CLOSE(12)
  RETURN
END FUNCTION getnrec

!=======================================================================
! PURPOSE: print error messages and stop
!=======================================================================
SUBROUTINE stoperr(prog,mesg1,mesg2)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: prog   ! name of calling prog routine
  CHARACTER(LEN=*),INTENT(IN) :: mesg1  ! nature of problem
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: mesg2  ! 2nd message (optional)
  
  CHARACTER(LEN=200) :: addmesg 
!  LOGICAL,PARAMETER :: kill=.FALSE.
  LOGICAL,PARAMETER :: kill=.TRUE.

  IF (PRESENT(mesg2)) THEN ; addmesg=mesg2
  ELSE                     ; addmesg=' '
  ENDIF 

  WRITE(6,'(a)') ' '
  WRITE(6,'(a)') '--error-- in: '//TRIM(prog)
  WRITE(6,'(a)') TRIM(mesg1)
  IF (addmesg/=' ') WRITE(6,'(a)') TRIM(addmesg)

  IF (kill) STOP "in stoperr"

END SUBROUTINE stoperr


END MODULE boxtool
