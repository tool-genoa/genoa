MODULE iteration
!$ USE omp_lib
IMPLICIT NONE
CONTAINS 

!=======================================================================
! PURPOSE: compute the production and removal rates of the species. 
!  => Gauss-Seidel iterations for twostep.
!=======================================================================
SUBROUTINE iterfirst(rk,prevc,gdt,soa_fg,nsat,idasat,numaou,idaou,nvoc,&
                     c,xfr,xfl,tnow)
  USE parameter_mod, ONLY: smallc,fg_iterKgp
  USE mechdata, ONLY: numsp,numre,numself,idselfreac,kdilu,kaerloss
  USE frozenstuff, ONLY: resetc,resetrate,nfix
  USE mechdata, ONLY: rxrgnt
  USE lpmaptool, ONLY: lpmap
  IMPLICIT NONE

  REAL, INTENT(IN)    :: gdt
  REAL, INTENT(IN)    :: tnow
  REAL, INTENT(IN)    :: rk(:)
  REAL, INTENT(IN)    :: prevc(:)

  INTEGER, INTENT(IN) :: soa_fg
  INTEGER, INTENT(IN) :: nsat
  INTEGER, INTENT(IN) :: idasat(:)
  INTEGER, INTENT(IN) :: numaou
  INTEGER, INTENT(IN) :: idaou(:)
  REAL, INTENT(IN)    :: nvoc

  REAL,INTENT(INOUT)  :: c(:)
  REAL,INTENT(OUT)    :: xfr(:)
  REAL,INTENT(OUT)    :: xfl(:)

! local variable
  REAL    :: rate(SIZE(rk))    ! reaction rate
  INTEGER :: i, j, ire
  REAL    :: ctotaer

! Initialize
  xfr(:)=0.  ;  xfl(:)=0.  ;  rate(:)=rk(:)

! reset concentration if requested
  IF (nfix/=0) CALL resetc(c, tnow)

! compute the reaction rate. The stoi. coef. of the reactants is 1,
! except for the species given in idselfreac where it is 2.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ire,i) SHARED(rate,c,rxrgnt,numre)
  DO ire=1,numre
    DO i=1,rxrgnt(ire)%nrg
      rate(ire)=rate(ire)*c(rxrgnt(ire)%idrg(i))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  DO i=1,numself
    ire=idselfreac(i,1) ; j=idselfreac(i,2)
    rate(ire)=rate(ire)*c(j)
  ENDDO

! compute rate for part. -> gas at each time step (if requested)
  IF (fg_iterKgp) THEN
    ctotaer = nvoc + SUM(c(idasat(1:nsat)))
    IF (ctotaer<1E8) ctotaer=1E8  ! overwrite "tiny seed" if ctotaer is too small 

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ire,i) SHARED(rate,idaou,ctotaer,numaou)
    DO i=1,numaou
      ire=idaou(i)
      rate(ire)=rate(ire)/ctotaer 
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

! Compute the production and loss expression for each species. For the loss term, 
! stoi. coef. is 1, expect for the species given in idselfreac where it is 2.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j) SHARED(lpmap,rate,xfl,xfr,numsp) SCHEDULE(DYNAMIC,numsp/50)
  DO i=1,numsp
    DO j=1,SIZE(lpmap(i)%idl)
      xfl(i) = xfl(i) + rate(lpmap(i)%idl(j))
    ENDDO

    DO j=1,SIZE(lpmap(i)%idp)
      xfr(i) = xfr(i) + rate(lpmap(i)%idp(j))*lpmap(i)%stpd(j)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

! correct the loss term for self reaction
  DO i=1,numself
    ire=idselfreac(i,1) ; j=idselfreac(i,2)
    xfl(j) = xfl(j) + rate(ire)
  ENDDO

! correct the loss rate and add dilution to loss
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(xfl,c,kdilu,numsp) 
  DO i=1,numsp
    xfl(i) = xfl(i)/c(i) + kdilu 
  ENDDO
!$OMP END PARALLEL DO

! add aerosol loss term (e.g. to the wall) 
  IF (soa_fg == 2) THEN
    IF (kaerloss > 1E-30) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(xfl,idasat,kaerloss,nsat) 
      DO i=1,nsat
        xfl(idasat(i)) = xfl(idasat(i)) + kaerloss
      ENDDO 
!$OMP END PARALLEL DO
    ENDIF
  ENDIF

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(c,xfr,xfl,gdt,prevc,numsp) 
  DO i=1,numsp
    c(i)=MAX(smallc, (prevc(i)+gdt*xfr(i)) / (1.+gdt*xfl(i)) )
  ENDDO
!$OMP END PARALLEL DO

! reset concentrations and rates if requested 
  IF (nfix/=0) THEN
    CALL resetc(c, tnow)
    CALL resetrate(xfr,xfl)
  ENDIF 

END SUBROUTINE iterfirst
END MODULE iteration
