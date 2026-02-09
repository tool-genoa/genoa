!-----------------------------------------------------------------
!------------- MODULE FOR NETCDF OUTPUT FILES --------------------
!-----------------------------------------------------------------
!SUBROUTINE opennc(filenam,numsp,mxlsp) 
!SUBROUTINE wrtnc_const(chrsp) 
!SUBROUTINE wrtnc_tchge(time,cbox) 
!SUBROUTINE ncdf_close
!-----------------------------------------------------------------
MODULE ncdf_saving
  USE netcdf               ! use netcdf library
  IMPLICIT NONE
  
  ! ID for the nc file
  INTEGER :: ncid          
  
  ! ID for the various dimension
  INTEGER :: dimid_time    ! number of record (i.e. time steps)
  INTEGER :: dimid_numsp   ! number of species
  INTEGER :: dimid_mxlsp   ! max length of species names
  INTEGER :: dimid_nsat    ! number of species in aerosol phase (only with equilibrium mass transfer)
  
  ! ID for the variables
  INTEGER :: varid_chrsp   ! lst (name) of species
  INTEGER :: varid_namsat  ! lst (name) of species in aerosol (only with equilibrium mass transfer)
  INTEGER :: varid_time    ! time (at record)
  INTEGER :: varid_conc    ! concentration
  INTEGER :: varid_caer    ! aerosol concentration (only with equilibrium mass transfer)

  ! other stuff
  CHARACTER(LEN=64)   :: title 
  CHARACTER(LEN=132)  :: institute 
  CHARACTER(LEN=64)   :: source 
  CHARACTER(LEN=32)   :: systime
  CHARACTER(LEN=32)   :: usrname  
  CHARACTER(LEN=1024) :: history
 
  INTEGER :: nrec    ! number of record (time step) 

CONTAINS

!-----------------------------------------------------------------
! open the nc file
!-----------------------------------------------------------------
SUBROUTINE opennc(filenam,numsp,mxlsp,nsat) 
  IMPLICIT NONE
  INTEGER :: retval
  CHARACTER(LEN=*),INTENT(IN)   :: filenam   ! name of the output file
  INTEGER,INTENT(IN) :: numsp                ! number of species in the mechanism
  INTEGER,INTENT(IN) :: mxlsp                ! max length of species names
  INTEGER,INTENT(IN),OPTIONAL :: nsat        ! number of species in aerosol phase (only with equilibrium mass transfer)

  nrec = 0   ! initialize the number of record stored 
  CALL fdate(systime)
  CALL getenv('USER',usrname)
  history = 'created by '//usrname(1:len_trim(usrname))// &
            ' on '// systime(1:len_trim(systime))
  title = 'box model simulation results using gecko-a mechanism'
  institute = 'LISA, UMR CNRS 7583, Univ Paris Est Creteil, Universite Paris Cite'
  
!-----------------------------
!----- open the file
!-----------------------------
  retval = nf90_create(filenam, NF90_NETCDF4, ncid)
  IF (retval /= NF90_NOERR) CALL handle_error(retval)
  
  retval = nf90_put_att(ncid, NF90_GLOBAL, 'title',trim(title))
  if (retval /= 0) call handle_error(retval)
  retval = nf90_put_att(ncid, NF90_GLOBAL, 'history',trim(history))
  if (retval /= 0) call handle_error(retval)
  retval = nf90_put_att(ncid, NF90_GLOBAL, 'institute', trim(institute))
  if (retval /= 0) call handle_error(retval)
  retval = nf90_put_att(ncid, NF90_GLOBAL, 'source', trim(source))
  if (retval /= 0) call handle_error(retval)

!-----------------------------
!----- declare dimension -----
!-----------------------------

  ! max length species
  retval = nf90_def_dim(ncid, 'lenspecies', mxlsp, dimid_mxlsp)
  if (retval /= NF90_NOERR) CALL handle_error(retval)
  print*,'mxlsp=',mxlsp

  ! species number
  retval = nf90_def_dim(ncid, 'numspecies', numsp, dimid_numsp)
  if (retval /= NF90_NOERR) CALL handle_error(retval)
  print*,'numsp=',numsp

  IF (PRESENT(nsat)) THEN
  ! species number in aerosol phase
    retval = nf90_def_dim(ncid, 'nsat', nsat, dimid_nsat)
    if (retval /= NF90_NOERR) CALL handle_error(retval)
    print*,'nsat=',nsat
  ENDIF
  
  ! number of time steps
  retval = nf90_def_dim(ncid,'ntstep',NF90_UNLIMITED,dimid_time)
  if (retval /= NF90_NOERR) CALL handle_error(retval)

!-----------------------------
!----- declare variables -----
!-----------------------------
  ! species ID  
  IF (PRESENT(nsat)) THEN
    retval = nf90_def_var(ncid,'Gaseous_Species', NF90_CHAR,[dimid_mxlsp, dimid_numsp], varid_chrsp)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
    retval=nf90_put_att(ncid,varid_chrsp,'comment','all gaseous species name in the mechanism')
    IF (retval /= NF90_NOERR) CALL handle_error(retval)

    retval = nf90_def_var(ncid,'Aerosol_Species', NF90_CHAR,[dimid_mxlsp, dimid_nsat], varid_namsat)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
    retval=nf90_put_att(ncid,varid_namsat,'comment','all aerosol species name in the mechanism')
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
  ELSE
    retval = nf90_def_var(ncid,'Species', NF90_CHAR,[dimid_mxlsp, dimid_numsp], varid_chrsp)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
    retval=nf90_put_att(ncid,varid_chrsp,'comment','all species name in the mechanism')
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
  ENDIF


  ! time at time steps
  retval = nf90_def_var(ncid,'time',NF90_REAL,[dimid_time], varid_time)
  IF (retval /= NF90_NOERR) CALL handle_error(retval)
  retval=nf90_put_att(ncid,varid_time,'units','second')
  IF (retval /= NF90_NOERR) CALL handle_error(retval)
  
  ! species concentration at time steps
  IF (PRESENT(nsat)) THEN
    retval = nf90_def_var(ncid,'gaseous_concentration',NF90_REAL,[dimid_numsp,dimid_time], varid_conc)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
    retval=nf90_put_att(ncid,varid_conc,'units','molec/cm3')

    retval = nf90_def_var(ncid,'aerosol_concentration',NF90_REAL,[dimid_nsat,dimid_time], varid_caer)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
    retval=nf90_put_att(ncid,varid_caer,'units','molec/cm3')
  ELSE
    retval = nf90_def_var(ncid,'concentration',NF90_REAL,[dimid_numsp,dimid_time], varid_conc)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
    retval=nf90_put_att(ncid,varid_conc,'units','molec/cm3')
  ENDIF
 
!-----------------------------
!----- close declaration mode
!-----------------------------
  retval = nf90_enddef(ncid)
  IF (retval /= NF90_NOERR) CALL handle_error(retval)

END SUBROUTINE opennc 

!=======================================================================
! write time independent variables (e.g. species name)
!=======================================================================
SUBROUTINE wrtnc_const(chrsp,namsat) 
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)  ! list (names) of the species
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL::  namsat(:)   ! name of the partitioning species

  INTEGER :: retval

  ! write name of species 
  retval = nf90_put_var(ncid, varid_chrsp, chrsp)
  IF (retval /= NF90_NOERR) CALL handle_error(retval)

  IF (PRESENT(namsat)) THEN 
    retval = nf90_put_var(ncid, varid_namsat, namsat)
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
  ENDIF

END SUBROUTINE wrtnc_const

!=======================================================================
! write time dependent variables (e.g. concentrations)
!=======================================================================
SUBROUTINE wrtnc_tchge(time,cbox,caer) 
  IMPLICIT NONE
  REAL,INTENT(IN) :: time
  REAL,INTENT(IN) :: cbox(:)
  REAL,INTENT(IN),OPTIONAL :: caer(:)

  INTEGER :: nspe,nsat
  INTEGER :: retval

  nrec=nrec+1   ! add a record
  nspe = SIZE(cbox)
  
  ! write time at time step
  retval = nf90_put_var(ncid,varid_time,time, start = [nrec])
  IF (retval /= NF90_NOERR) CALL handle_error(retval)
  
  ! write concentration at time step
  retval = nf90_put_var(ncid,varid_conc,cbox,start=(/1,nrec/),count=(/nspe,1/))
  IF (retval /= NF90_NOERR) CALL handle_error(retval)

  IF (PRESENT(caer)) THEN
    nsat=SIZE(caer) 
    retval = nf90_put_var(ncid,varid_caer,caer,start=(/1,nrec/),count=(/nsat,1/))
    IF (retval /= NF90_NOERR) CALL handle_error(retval)
  ENDIF

END SUBROUTINE wrtnc_tchge


! ----------------------------------------------------------------------
! Subroutine to close ncdf files
! ----------------------------------------------------------------------
SUBROUTINE ncdf_close
  IMPLICIT NONE
  INTEGER :: retval

  retval = nf90_close(ncid)
  IF (retval /= 0) CALL handle_error(retval)

END SUBROUTINE ncdf_close

!==================================
SUBROUTINE handle_error(retval)
  INTEGER, INTENT(IN) :: retval
  CHARACTER(LEN=256) :: ncstring

  ! print the error msg provided by netcdf
  ncstring = nf90_strerror(retval)
  PRINT *, trim(ncstring)
  STOP 'Error in NetCDF operation'
END SUBROUTINE handle_error

END MODULE ncdf_saving
