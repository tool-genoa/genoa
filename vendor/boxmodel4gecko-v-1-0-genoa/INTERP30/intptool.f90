MODULE intptool
IMPLICIT NONE
CONTAINS

!SUBROUTINE stringarg(line,narg,posi)
!INTEGER FUNCTION search(aseek,alist,nlist)
!SUBROUTINE bubblesort(numsp,chrsp,isortlk,sortsp)
!SUBROUTINE linuxsort(numsp,chrsp,isortlk,sortsp)
!SUBROUTINE quicksort(numsp,chrsp,isortlk,sortsp)

!======================================================================= 
! PURPOSE: This subroutine returns the number of arguments (narg) and 
! position (posi) in the string (line) provided as input.
!======================================================================= 
SUBROUTINE stringarg(line,narg,posi)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: line ! line to be analyzed
  INTEGER,INTENT(OUT) :: narg       ! number of arguments in the line
  INTEGER,INTENT(OUT) :: posi(:,:)  ! position of the arguments
                                    ! posi(i,1): start the ith argument
                                    ! posi(i,2): end the ith argument

! local
  INTEGER :: istart,istop, ipos, lenstr
  INTEGER :: maxarg  ! max # of arguments in the string (size of "posi")
  LOGICAL :: loarg

! initialize
  istart=0 ; narg=0 ; posi(:,:)=0 ; lenstr=LEN(line)
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
        posi(narg,1)=ipos ; posi(narg,2)=ipos
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

!======================================================================= 
! PURPOSE: search if a name (aseek) exist in a sorted list (alist). 
! If found : return ! the index in the list - otherwise return a 
! negative index (where aseek should have been found) 
!======================================================================= 
INTEGER FUNCTION search(aseek,alist,nlist)
  IMPLICIT NONE

  INTEGER,INTENT(IN)          :: nlist        ! # of element in the list                    
  CHARACTER(LEN=*),INTENT(IN) :: aseek        ! string to be find in the list
  CHARACTER(LEN=*),INTENT(IN) :: alist(nlist) ! sorted list of strings

  INTEGER jhi, jlo, jold, j

! initialize:
  search = 0 ; jold = 0 ; jlo = 1 ; jhi = nlist + 1

! search loop
  searchloop : &
  DO                              ! exit the loop if same index is found
    j = (jhi+jlo)/2               ! new index
       
    IF (j == jold) THEN           ! new=old index => not found
      search = -j
      RETURN
    ENDIF

    jold = j                      ! store last position
    IF (aseek > alist(j)) THEN    ! index too small
      jlo = j
      CYCLE searchloop
    ELSE IF (aseek == alist(j)) THEN  ! Bingo
      search = j
      RETURN
    ELSE                          ! index too large
      jhi  = j
      CYCLE searchloop
    ENDIF

  ENDDO  searchloop
END FUNCTION search

!======================================================================= 
! PURPOSE: sort the list of species based on simple bubble sort. 
! Sorting is not efficient => to be used only for small mechanism 
!======================================================================= 
SUBROUTINE bubblesort(numsp,chrsp,isortlk,sortsp)
  IMPLICIT NONE

  INTEGER,INTENT(IN)           :: numsp      ! # of species
  CHARACTER(LEN=*),INTENT(IN)  :: chrsp(:)   ! list of species
  INTEGER,INTENT(OUT)          :: isortlk(:) ! sorted list of species
  CHARACTER(LEN=*),INTENT(OUT) :: sortsp(:)  ! index link "sorted" -> "mechanism" species

  INTEGER :: i,j
  CHARACTER(LEN=LEN(sortsp(1)))  :: temp
  
  isortlk(:)=0  ; sortsp(:)=chrsp(:)
  
  DO i=1,numsp
    DO j=numsp,i+1,-1
      IF (sortsp(j-1) > sortsp(j) )THEN
        temp=sortsp(j-1)
        sortsp(j-1)=sortsp(j)
        sortsp(j)=temp
      ENDIF
    ENDDO
  ENDDO

! get the index of species
  DO i=1,numsp
    j = search(chrsp(i),sortsp,numsp)
    IF (j <= 0) STOP 'in bubblesort, species not found'
    isortlk(j)=i
  ENDDO
END SUBROUTINE bubblesort

!======================================================================= 
! PURPOSE: sort the list of species based on the linux "sort" command.
! Linux sort appears to be very efficient but requires to temporary
! write a file with the list of species and perform a system call
!======================================================================= 
SUBROUTINE linuxsort(numsp,chrsp,isortlk,sortsp)
  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: numsp      ! # of species
  CHARACTER(LEN=*),INTENT(IN)   :: chrsp(:)   ! list of species
  INTEGER, INTENT(OUT)          :: isortlk(:) ! sorted list of species
  CHARACTER(LEN=*), INTENT(OUT) :: sortsp(:)  ! index link "sorted" -> "mechanism" species

  INTEGER :: i
    
  OPEN(55, FILE='dummy') 
  DO i=1,numsp
    WRITE(55,'(a,i9)') chrsp(i),i
  ENDDO 
  CLOSE(55)
  CALL SYSTEM('LC_ALL=C sort -k 1,1 dummy -o sortedlist')
  CALL SYSTEM('rm -f dummy')

! Read the sorted list and store in sortsp. "isortlk" provide the  
! link (index) between the unsorted and sorted lists
  OPEN(55,FILE='sortedlist',STATUS='old')
  DO i=1,numsp
    READ(55,'(a,i9)') sortsp(i),isortlk(i)
  ENDDO
  CLOSE(55)
END SUBROUTINE linuxsort 

!======================================================================= 
! PURPOSE: sort the list of species based on an "efficient" sorting 
! routine.
!======================================================================= 
SUBROUTINE quicksort(numsp,chrsp,isortlk,sortsp)
  USE sortstring
  IMPLICIT NONE

  INTEGER,INTENT(IN)           :: numsp      ! # of species
  CHARACTER(LEN=*),INTENT(IN)  :: chrsp(:)   ! list of species
  INTEGER,INTENT(OUT)          :: isortlk(:) ! sorted list of species
  CHARACTER(LEN=*),INTENT(OUT) :: sortsp(:)  ! index link "sorted" -> "mechanism" species

! local
  INTEGER :: i,j

  isortlk(:)=0 ;  sortsp(:)=chrsp(:)
  
  CALL sort_string(sortsp(1:numsp))
    
! get the index of species
  DO i=1,numsp
    j = search(chrsp(i),sortsp,numsp)
    IF (j <= 0) STOP 'in quicksort, species not found'
    isortlk(j)=i
  ENDDO
END SUBROUTINE quicksort 

END MODULE intptool
 
