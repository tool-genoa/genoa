MODULE frozenstuff
IMPLICIT NONE

! concentration externally fixed
INTEGER, PARAMETER   :: mxfix=10        ! max number of species for which concentration can be fixed
INTEGER              :: nfix            ! number of species having fixed concentration
INTEGER              :: idfix(mxfix)    ! id of the species having fixed concentration
REAL                 :: cfix(mxfix,24)  ! value of the fixed concentration
REAL                 :: envfix(5,24)    ! value of the fixed environmental variables
                                        ! (flags, temperature, relative humidity, solar zenith angle, pressure)
INTEGER              :: noxfix          ! flag to fix NOx (NO+NO2) concentration (0: not fixed, 1:fixed) 
REAL                 :: cnox            ! value of the fixed NOx concentration

! id of some commonly used species
INTEGER :: idch3o2     ! ID for CH3O2
INTEGER :: idno        ! ID for NO
INTEGER :: idno2       ! ID for NO2
INTEGER :: idno3       ! ID for NO3
INTEGER :: idoh        ! ID for OH
INTEGER :: idho2       ! ID for HO2
INTEGER :: ido3        ! ID for O3
INTEGER :: idhono      ! ID for HONO (HNO2)

CONTAINS
!=======================================================================
! Purpose: fix the concentration for the species having a fixed 
! (externally prescribed) value.
!=======================================================================
SUBROUTINE resetc(c,tnow)
  IMPLICIT NONE
  REAL, INTENT(INOUT) :: c(:)
  REAL, INTENT(IN)    :: tnow
  INTEGER :: i, t ! loop index

  IF (nfix/=0) THEN
    t = MAX(1, MIN(24, CEILING(MODULO(tnow, 86400.)/3600.)))  ! time in hours
    DO i=1,nfix
      c(idfix(i))=cfix(i,t)
    ENDDO
  ENDIF
END SUBROUTINE  resetc

!=======================================================================
! Purpose: set the production and removal rate to 0 for the species 
! having fixed concentration (externally prescribed value).
!=======================================================================
SUBROUTINE resetrate(xfr,xfl)
  IMPLICIT NONE
  REAL, INTENT(INOUT) :: xfr(:)
  REAL, INTENT(INOUT) :: xfl(:)
  
  INTEGER :: i

  IF (nfix/=0) THEN
    DO i=1,nfix
      xfr(idfix(i))=0.
      xfl(idfix(i))=0.
    ENDDO
  ENDIF
END SUBROUTINE  resetrate

!=======================================================================
! Purpose: fix the NOx (NO+NO2) concentration to the prescribed value. 
! The  NO/NO2 ratio is unchanged. 
!=======================================================================
SUBROUTINE resetnox(c)
  IMPLICIT NONE
  REAL,INTENT(INOUT) :: c(:)
  
  REAL :: rap

  IF (noxfix==1) THEN
    IF (c(idno2)==0.) THEN  
      c(idno)=cnox
    ELSE
      rap=c(idno)/c(idno2)
      c(idno2)=cnox*(1./(rap+1.)) 
      c(idno)=cnox*(rap/(rap+1.))
    ENDIF
  ENDIF
END SUBROUTINE resetnox

!=======================================================================
! Purpose: store ID for some commonly used species 
!=======================================================================
SUBROUTINE storeid(lout,numsp,chrsp)
  USE boxtool, ONLY: spe2id
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: lout
  INTEGER,INTENT(IN) :: numsp
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)

  idch3o2=spe2id(lout,chrsp,numsp,chem='GCH3O2',csub='storeid')
  idno=spe2id(lout,chrsp,numsp,chem='GNO',csub='storeid')  
  idno2=spe2id(lout,chrsp,numsp,chem='GNO2',csub='storeid')
  idno3=spe2id(lout,chrsp,numsp,chem='GNO3',csub='storeid')
  idho2=spe2id(lout,chrsp,numsp,chem='GHO2',csub='storeid')
  ido3=spe2id(lout,chrsp,numsp,chem='GO3',csub='storeid')
  idoh=spe2id(lout,chrsp,numsp,chem='GHO',csub='storeid') 
  idhono=spe2id(lout,chrsp,numsp,chem='GHNO2',csub='storeid')
END SUBROUTINE storeid

END MODULE frozenstuff
