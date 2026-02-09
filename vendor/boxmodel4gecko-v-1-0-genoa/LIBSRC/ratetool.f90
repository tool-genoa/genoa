MODULE ratetool
!$ USE omp_lib
IMPLICIT NONE
CONTAINS
!=======================================================================  
! PURPOSE: Compute the rate coefficient of the reactions in the 
! mechanism as a function T, as well as other parameters for "special"  
! reaction type (e.g. falloff, extra, RO2+RO2 ...)
!=======================================================================  
SUBROUTINE ratecoef(numre,numo2,num_m,numfo,numextra,nummeo2,&
                    numiso,numreacro2,ncpero,&
                    ido2,id_m,idfo,idextra,idmeo2,idreacro2,idiso,&
                    arrhcf,focf,extracf,isocf,&
                    cmeo2,cro2,temp,sumc,water,qfor)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: numre             ! number of reaction in the mechanism
  INTEGER,INTENT(IN) :: numo2             ! # of 'OXYGEN' rx
  INTEGER,INTENT(IN) :: num_m             ! # of 'TBODY' rx
  INTEGER,INTENT(IN) :: numfo             ! # of 'FALLOFF' rx
  INTEGER,INTENT(IN) :: numextra          ! # of 'EXTRA' rx
  INTEGER,INTENT(IN) :: nummeo2           ! # of 'MEPERO' (CH3O2) rx
  INTEGER,INTENT(IN) :: numiso            ! # of 'ISOM' rx
  INTEGER,INTENT(IN) :: numreacro2(:)     ! # of rx with "PEROj" in class [j]
  INTEGER,INTENT(IN) :: ncpero            ! # of peroxy classes (for RO2+RO2 reactions)
  INTEGER,INTENT(IN) :: ido2(:)           ! mechanism rx index of the ith 'OXYGEN' reaction
  INTEGER,INTENT(IN) :: id_m(:)           ! mechanism rx index of the ith 'TBODY' reaction
  INTEGER,INTENT(IN) :: idfo(:)           ! mechanism rx index of the ith 'FALLOFF' reaction
  INTEGER,INTENT(IN) :: idextra(:)        ! mechanism rx index of the ith 'EXTRA' reaction
  INTEGER,INTENT(IN) :: idmeo2(:)         ! mechanism rx index of the ith 'MEPERO' reaction 
  INTEGER,INTENT(IN) :: idreacro2(:,:)    ! [i,j] mechanism rx index of the ith 'PEROj' reaction
  INTEGER,INTENT(IN) :: idiso(:)          ! mechanism rx index of the ith 'ISOM' reaction
  REAL,INTENT(IN)  :: arrhcf(:,:)         ! [i,3] arrhenius parameter for the ith mech. rx index
  REAL,INTENT(IN)  :: focf(:,:)           ! [i,j] j coef. (k0, Fc...) for the ith falloff rx
  REAL,INTENT(IN)  :: extracf(:,:)        ! [i,j] j coef. (k...) for the ith "extra" rx)
  REAL,INTENT(IN)  :: isocf(:,:)          ! [i,j] j extra coef. for the ith 'ISOM' rx
  REAL,INTENT(IN)  :: cmeo2               ! CH3O2 concentration 
  REAL,INTENT(IN)  :: cro2(:)             ! RO2 concentration for 'PERO' in class [j]
  REAL,INTENT(IN)  :: temp                ! temperature
  REAL,INTENT(IN)  :: sumc                ! M (in molec/cm3)
  REAL,INTENT(IN)  :: water               ! water concentration (molec/cm3)
  REAL,INTENT(OUT) :: qfor(:)             ! rate coefficient for the ith reaction

  INTEGER :: ire, i,j
  REAL    :: c_o2,t4,t3,t2

  CHARACTER(LEN=10),PARAMETER :: progname='ratecoef'
  CHARACTER(LEN=80) :: mesg1 

  c_o2=0.2*sumc ; qfor(:)=0.0

  DO i=1,ncpero
    IF (cro2(i) < 0.) THEN 
      WRITE(mesg1,'(a,i2)') "unexpected negative c for ro2 class: ", i
      CALL stoperr(progname,mesg1)      
    ENDIF
  ENDDO

! compute rate constant
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ire) SHARED(qfor,arrhcf,temp,numre)
  DO ire=1,numre
    qfor(ire)=arrhcf(ire,1)*(temp**(arrhcf(ire,2)))*exp(-arrhcf(ire,3)/temp)
  ENDDO
!$OMP END PARALLEL DO 

! reaction with third body M
  DO i=1,num_m
    ire=id_m(i)
    qfor(ire)=qfor(ire)*sumc
  ENDDO

! reaction with O2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire) SHARED(qfor,c_o2,numo2,ido2)
  DO i=1,numo2
    ire=ido2(i)
    qfor(ire)=qfor(ire)*c_o2
  ENDDO
!$OMP END PARALLEL DO 

! reaction with CH3O2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire) SHARED(qfor,cmeo2,nummeo2,idmeo2)
  DO i=1,nummeo2
    ire=idmeo2(i)
    qfor(ire)=qfor(ire)*cmeo2
  ENDDO
!$OMP END PARALLEL DO 

! reaction with RO2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,ire) SHARED(qfor,cro2,idreacro2,numreacro2,ncpero)
  DO i=1,ncpero
    DO j=1,numreacro2(i)
      ire=idreacro2(j,i)
      qfor(ire)=qfor(ire)*cro2(i)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO 

! isomerisation reaction
  t4=temp**4 ; t3=temp**3 ; t2=temp**2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ire) SHARED(qfor,temp,t4,t3,t2,isocf,numiso,idiso)
  DO i=1,numiso
    ire=idiso(i)
    qfor(ire) = qfor(ire)*(isocf(i,1)*t4 + isocf(i,2)*t3 + isocf(i,3)*t2 + &
                           isocf(i,4)*temp + isocf(i,5) )
  ENDDO
!$OMP END PARALLEL DO 

! fall off reaction
  DO i=1,numfo
    ire=idfo(i)
    qfor(ire) = forate(focf(i,:),temp,sumc,arrhcf(ire,:))
  ENDDO

! extra reaction
  DO i=1,numextra
    ire=idextra(i)
    qfor(ire) = extrarate(extracf(i,:),temp,sumc,water,qfor(ire))
  ENDDO

END SUBROUTINE ratecoef

!=======================================================================
! Purpose: Compute the rate coefficient for the reactions using the
! keyword EXTRA in the mechanism, i.e. 
!      A + B + EXTRA => X +Y
!              EXTRA /label data1 data2 .../
! where data_x are used to compute these "special" reaction rate coef.
!=======================================================================
REAL FUNCTION extrarate(coef,temp,sumc,water,qrxn)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE

  REAL,INTENT(IN) :: coef(:)     ! 'extra' coef. for the current rxn
  REAL,INTENT(IN) :: temp        ! temperature
  REAL,INTENT(IN) :: sumc        ! M (in molec/cm3)
  REAL,INTENT(IN) :: water       ! water concentration (molec/cm3)
  REAL,INTENT(IN) :: qrxn         ! rate coefficient for the current rxn

  INTEGER :: label
  REAL    :: xk,xk0,xk1,xk2,xk3,kd,water_dimer
  CHARACTER(LEN=10),PARAMETER :: progname='extrarate'
  CHARACTER(LEN=80) :: mesg1

! CASE SELECTOR FOR LABEL
! -----------------------------
  label=NINT(coef(1))
  SELECT CASE(label)   
  
! label 100: O+O2+M=>O3 // 3rd order reaction 
    CASE (100)
      extrarate=qrxn*sumc*0.2*sumc  

! label 500: reaction with water vapor  
    CASE (500)
      extrarate=qrxn*water

! label 501: reaction with H20 (xk1) and H2O+M (xk2)  // specific to HO2+HO2 => H2O2
    CASE (501)
      xk1=qrxn*water                                      
      xk2=coef(2)*EXP(-coef(4)/temp)*water*sumc
      extrarate = xk1+xk2                                      

! label 502: reaction with water dimer (Kd according to Scribano et al., 2006)
    CASE (502)
      kd = 4.7856e-4*exp(1851/temp-5.10485e-3*temp)    ! atm^(-1)
      kd = kd*8.314*temp*1e6/(1.01325e5*6.02E23)       ! in molec-1 cm3
      water_dimer=kd*(water**2.)
      extrarate = qrxn * water_dimer

! label 550: OH+HNO3 reaction: k=k0 + k3*M/(1+K3*M/K2)
    CASE (550) 
      xk0=qrxn
      xk2=coef(2)*(temp**coef(3))*exp(-coef(4)/temp)
      xk3=coef(5)*(temp**coef(6))*exp(-coef(7)/temp)
      xk=xk0 + xk3*sumc / (1.+ ((xk3*sumc)/xk2) ) 
      extrarate=xk

! unidentified label
    CASE DEFAULT
      WRITE(mesg1,*) "Label unknown:",label
      CALL stoperr(progname,mesg1)

  END SELECT
END FUNCTION extrarate

!=======================================================================
! Purpose: Compute reaction rate coefficient for "falloff" type  
! reaction. A "Troe" expression is used :
! k=[koM/(1+koM/ki))]*Fc^[(1+log10(koM/ki))^2]^-1
! where ko is the low pressure rate constant and ki is the high 
! pressure rate constant. 
! ki: provided as the "reaction" arrhenius parameter (arrh)
! ko: provided as the 1st, 2nd and 3rd data in table coef
! Fc: provided as the 4th date in table coef
!=======================================================================
REAL FUNCTION forate(coef,temp,sumc,arrh)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE

  REAL,INTENT(IN) :: temp       ! temperature
  REAL,INTENT(IN) :: sumc       ! M (in molec/cm3)
  REAL,INTENT(IN) :: coef(:)    ! coef. (k0, Fc...) for the current rxn
  REAL,INTENT(IN) :: arrh(:)    ! arrhenius coef for the current rxn

  REAL :: xkom,xki,kratio,factor
  CHARACTER(LEN=10),PARAMETER :: progname='forate'
  CHARACTER(LEN=80) :: mesg1
  
! patch MCM - see http://mcm.york.ac.uk/parameters/complex.htt for details
  REAL :: fcent,factor2  

! xkom: ko*M  (effective low pressure rate coef.) /// xki: k_infiny   
  xki = arrh(1)*((temp/300.)**arrh(2))*exp(-arrh(3)/temp)
  xkom = coef(1)*((temp/300.)**coef(2))*exp(-coef(3)/temp)*sumc
  kratio = xkom/xki
  factor = 1./(1.+(log10(kratio))**2.)

! First IF uses the usual Troe equation to compute rate constant
  IF ((coef(4)/=0.) .AND. (coef(5)==0.)) THEN
    forate = (xkom/(1.+kratio)) * (coef(4)**factor)

! The following IF use a different equation for MCM rates constant
  ELSE IF ((coef(4)/=0.) .AND. (coef(5)==1.)) THEN
    factor2 = 1./(1.+(log10(kratio)/(0.75-1.27*log10(coef(4))))**2.)
    forate = (xkom/(1.+kratio)) * (coef(4)**factor2)

  ELSE IF ((coef(4)==0.00) .AND. (coef(5)==204.) .AND. &
           (coef(6)==0.17) .AND. (coef(7)==51.)) THEN
    fcent = exp(-temp/coef(5))+coef(6)+exp(-coef(7)/temp)
    factor2 = 1./(1.+(log10(kratio)/(0.75-1.27*log10(fcent)))**2.)
    forate = (xkom/(1.+kratio)) * (fcent**factor2)
  ELSE
    mesg1="unexpected setup of the auxiliary coef in falloff rxn" 
    CALL stoperr(progname,mesg1)
  ENDIF

END FUNCTION forate

END MODULE ratetool
