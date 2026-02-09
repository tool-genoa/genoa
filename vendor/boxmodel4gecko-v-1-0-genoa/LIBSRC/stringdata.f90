! ********************************************************************** 
! This subroutine returns the number of arguments (narg) and 
! position (posi) in the string (line) provided as input.
! ********************************************************************** 
MODULE stringdata
  IMPLICIT NONE
  CONTAINS
   
SUBROUTINE stringarg(line,narg,posi)
  IMPLICIT NONE

! input/output
  CHARACTER(LEN=*), INTENT(IN) :: line ! line to be analyzed

! output
  INTEGER, INTENT(OUT) :: narg       ! number of arguments in the line
  INTEGER, INTENT(OUT) :: posi(:,:)  ! position of the arguments
                                     ! posi(i,1): start the ith argument
                                     ! posi(i,2): end the ith argument

! local
  INTEGER :: istart,istop, ipos, lenstr
  INTEGER :: maxarg  ! max # of arguments in the string (size of "posi")
  LOGICAL :: loarg

! initialize
  istart=0
  narg=0 ; posi(:,:)=0 ; lenstr=LEN(line)
  maxarg=SIZE(posi,1)

! find last character - Return if empty line 
  istop=LEN_TRIM(line)
  IF (istop == 0) RETURN ! (no argument found in line)
   
! find first character 
  DO ipos=1,lenstr
    IF (line(ipos:ipos) /= ' ') THEN
      istart=ipos
      EXIT
    ENDIF
  ENDDO
 
! initialize the loop
  narg=1
  posi(narg,1)=istart
  posi(narg,2)=istart ! likely overwritten next
  loarg=.TRUE.
  DO ipos=istart,istop ! scroll the line
    IF (loarg) THEN
      IF (line(ipos:ipos)==' ') THEN
        loarg=.FALSE.
      ELSE
        posi(narg,2)=ipos
        CYCLE
      ENDIF
    ELSE
      IF (line(ipos:ipos)==' ') THEN
        CYCLE
      ELSE
        loarg=.TRUE.
        narg=narg+1
        posi(narg,1)=ipos
        posi(narg,2)=ipos
        CYCLE
      ENDIF
    ENDIF
  ENDDO
  posi(narg,2)=istop

! check that the number of arguments is below threshold (size of posi)
  IF (narg > maxarg) THEN
    PRINT*, '--error-- in stringarg. To many arguments in :'
    PRINT*, TRIM(line)
    PRINT*, 'Maximum number of arguments set to:', maxarg
    STOP 'in stringarg routine - too many arguments'
  ENDIF

END SUBROUTINE stringarg

END MODULE stringdata
