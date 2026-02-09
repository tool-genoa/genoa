MODULE wrtboftool
  IMPLICIT NONE
  INTEGER :: numsav       ! # of outputs (i.e. time step) written in bof
  CONTAINS
!=======================================================================
! PURPOSE: write the species in the binary output files (bof).
!=======================================================================
SUBROUTINE wrtspecies(lbof,lgbof,lpbof,lwbof,&
                      numsp,mxlsp,chrsp,wmol,soa_fg,nsat,idgsat, &
                      nspg,ptrgas,nspp,ptrpart,nspw,ptrwall,out_fg)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: lbof                 ! binary output file - all species
  INTEGER,INTENT(IN) :: lgbof                ! binary output file - gas phase only
  INTEGER,INTENT(IN) :: lpbof                ! binary output file - part. phase only
  INTEGER,INTENT(IN) :: lwbof                ! binary output file - wall phase only
  INTEGER,INTENT(IN) :: numsp                ! number of species in the mechanism
  INTEGER,INTENT(IN) :: mxlsp                ! max length of species names
  CHARACTER(LEN=*), INTENT(IN) :: chrsp(:)   ! name of the species
  REAL,INTENT(IN)    :: wmol(:)              ! molar mass of the ith species
  INTEGER,INTENT(IN) :: soa_fg               ! SOA flag (method used for gas/particle partioning)
  INTEGER,INTENT(IN) :: nsat                 ! # of partitioning species
  INTEGER,INTENT(IN) :: idgsat(:)            ! ID of the partitioning species
  INTEGER,INTENT(IN) :: nspg                 ! # of gas phase species
  INTEGER,INTENT(IN) :: ptrgas(:)            ! pointer (1st and last) gas phase species
  INTEGER,INTENT(IN) :: nspp                 ! # of part. phase species
  INTEGER,INTENT(IN) :: ptrpart(:)           ! pointer (1st and last) part. phase species
  INTEGER,INTENT(IN) :: nspw                 ! # of wall phase species
  INTEGER,INTENT(IN) :: ptrwall(:)           ! pointer (1st and last) wall phase species
  INTEGER,INTENT(IN) :: out_fg               ! flag for output: if 1: only out lbof, 2: out all
  
  INTEGER :: i

  IF (out_fg < 1) RETURN ! output nothing

  numsav = 0  ! initialize numsav 
 
  ! write all species included in the mechanism
  WRITE(lbof) numsp,mxlsp,numsav
  WRITE(lbof) (chrsp(i),i=1,numsp)
  WRITE(lbof) (wmol(i),i=1,numsp)

  IF (out_fg == 1) RETURN ! only output lbof
  
  ! soa_fg=1 ==> equilibrium at each time step. 
  IF (soa_fg == 1) THEN
    WRITE(lpbof) nsat,mxlsp,numsav
    WRITE(lpbof) (chrsp(idgsat(i)),i=1,nsat)

  ! soa_fg=2 ==> solve mass transport equations.
  ELSE IF (soa_fg == 2) THEN
    
    ! gas phase only
    WRITE(lgbof) nspg,mxlsp,numsav
    WRITE(lgbof) (chrsp(i),i=ptrgas(1),ptrgas(2))

    ! particle phase only
    IF (nspp>0) THEN
      WRITE(lpbof) nspp,mxlsp,numsav
      WRITE(lpbof) (chrsp(i),i=ptrpart(1),ptrpart(2))
    ENDIF
    
    ! "wall" phase only
    IF (nspw>0) THEN
      WRITE(lwbof) nspw,mxlsp,numsav
      WRITE(lwbof) (chrsp(i),i=ptrwall(1),ptrwall(2))
    ENDIF
  ENDIF
  
END SUBROUTINE wrtspecies

!=======================================================================
! PURPOSE: write concentration in the binary output files (bof).
!=======================================================================
SUBROUTINE wrtresu(time,temp,lbof,lgbof,lpbof,lwbof,&
                   numsp,cbox,caer,soa_fg,nsat,ptrgas, &
                   nspp,ptrpart,nspw,ptrwall,out_fg)
  IMPLICIT NONE
  REAL,INTENT(IN)    :: time        ! current time
  REAL,INTENT(IN)    :: temp        ! Ttemperature
  INTEGER,INTENT(IN) :: lbof        ! binary output file - all species
  INTEGER,INTENT(IN) :: lgbof       ! binary output file - gas phase only
  INTEGER,INTENT(IN) :: lpbof       ! binary output file - part. phase only
  INTEGER,INTENT(IN) :: lwbof       ! binary output file - wall phase only
  INTEGER,INTENT(IN) :: numsp       ! number of species in the mechanism
  REAL,INTENT(IN)    :: cbox(:)     ! concentration of species 
  REAL,INTENT(IN)    :: caer(:)     ! particle phase concentration (molec/cm3)
  INTEGER,INTENT(IN) :: soa_fg      ! SOA flag (method used for gas/particle partioning)
  INTEGER,INTENT(IN) :: nsat        ! # of partitioning species
  INTEGER,INTENT(IN) :: ptrgas(:)   ! pointer (1st and last) gas phase species
  INTEGER,INTENT(IN) :: nspp        ! # of part. phase species
  INTEGER,INTENT(IN) :: ptrpart(:)  ! pointer (1st and last) part. phase species
  INTEGER,INTENT(IN) :: nspw        ! # of wall phase species
  INTEGER,INTENT(IN) :: ptrwall(:)  ! pointer (1st and last) wall phase species
  INTEGER,INTENT(IN) :: out_fg      ! flag for output: if 1: only out lbof, 2: out all

  INTEGER :: isp

  IF (out_fg < 1) RETURN ! output nothing

  numsav = numsav + 1     !  raise the # of outputs counters
  
  ! all species included in the mechanism
  WRITE(lbof) time,temp,(cbox(isp),isp=1,numsp)

  IF (out_fg == 1) RETURN ! only output lbof

  ! soa_fg=1 ==> equilibrium at each time step. 
  IF (soa_fg == 1) THEN
    WRITE(lpbof) time,temp,(caer(isp),isp=1,nsat)

  ! soa_fg=2 ==> equilibrium at each time step. 
  ELSE IF (soa_fg == 2) THEN
    WRITE(lgbof) time,temp,(cbox(isp),isp=ptrgas(1),ptrgas(2))
    IF (nspp>0) WRITE(lpbof) time,temp,(cbox(isp),isp=ptrpart(1),ptrpart(2))
    IF (nspw>0) WRITE(lwbof) time,temp,(cbox(isp),isp=ptrwall(1),ptrwall(2))
  ENDIF
END SUBROUTINE wrtresu

!=======================================================================
! PURPOSE: overwrite the header of the binary output files (bof) using 
! actual number of time datasets (numsav) and close all bof. 
!=======================================================================
SUBROUTINE bofheader(lout,lbof,lgbof,lpbof,lwbof,numsp,mxlsp, &
                     nsat,nspg,nspp,nspw,soa_fg,out_fg)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: lout     ! file unit for information outputs
  INTEGER,INTENT(IN) :: lbof     ! binary output file - all species
  INTEGER,INTENT(IN) :: lgbof    ! binary output file - gas phase only
  INTEGER,INTENT(IN) :: lpbof    ! binary output file - part. phase only
  INTEGER,INTENT(IN) :: lwbof    ! binary output file - wall phase only
  INTEGER,INTENT(IN) :: numsp    ! number of species in the mechanism
  INTEGER,INTENT(IN) :: mxlsp    ! max length of species names
  INTEGER,INTENT(IN) :: nsat     ! # of partitioning species
  INTEGER,INTENT(IN) :: nspg     ! # of gas phase species
  INTEGER,INTENT(IN) :: nspp     ! # of part. phase species
  INTEGER,INTENT(IN) :: nspw     ! # of wall phase species
  INTEGER,INTENT(IN) :: soa_fg   ! SOA flag (method used for gas/particle partioning)
  INTEGER,INTENT(IN) :: out_fg   ! flag for output: if 1: only out lbof, 2: out all

  IF (out_fg < 1) RETURN ! output nothing

  ! write info to log file
  WRITE(lout,*)
  WRITE(lout,*) ' Binariy file has ',numsav,' time datasets'

  OPEN(lbof,POSITION='REWIND') 
  WRITE(lbof) numsp,mxlsp,numsav
  CLOSE(lbof)

  IF (out_fg == 1) RETURN ! only output lbof

  ! soa_fg=1 ==> equilibrium at each time step. 
  IF (soa_fg == 1) THEN
    OPEN(lpbof,POSITION='REWIND')  
    WRITE(lpbof) nsat,mxlsp,numsav
    CLOSE(lpbof)

  ! soa_fg=2 ==> solve mass transport equations.
  ELSE IF (soa_fg == 2) THEN
    
    ! gas phase only
    OPEN(lgbof,POSITION='REWIND')  
    WRITE(lgbof) nspg,mxlsp,numsav
    CLOSE(lgbof)

    ! particle phase only
    IF (nspp>0) THEN
      OPEN(lpbof,POSITION='REWIND')  
      WRITE(lpbof) nspp,mxlsp,numsav
      CLOSE(lpbof)
    ENDIF
    
    ! "wall" phase only
    IF (nspw>0) THEN
      OPEN(lwbof,POSITION='REWIND')  
      WRITE(lwbof) nspw,mxlsp,numsav
      CLOSE(lwbof)
    ENDIF
  ENDIF

END SUBROUTINE bofheader

END MODULE wrtboftool
