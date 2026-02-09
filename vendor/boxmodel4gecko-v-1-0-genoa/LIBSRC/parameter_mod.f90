MODULE parameter_mod

! parameters for boxmodel
  REAL, PARAMETER     :: smallc=1E-40        ! lower limit for a concentration (avoid 0.)
  INTEGER, PARAMETER  :: mxlsp=8             ! max length of species names 
  REAL                :: tstart,tstop        ! starting time and stopping time for simulation
  INTEGER             :: ntstep              ! # of printouts (binary file) & "refresh" time (T, rate, ...)
  INTEGER             :: nskip               ! skip printout in the binary file
  REAL                :: rtol, atol          ! relative & absolute tolerance - solver parameters
  INTEGER             :: numit               ! # number of iteration - solver parameters
  REAL                :: dtmin               ! minimum delta time - solver parameters

! Fixed parameters for properties calculation
  REAL, PARAMETER     :: avogadro=6.02214E23 ! avogadro number
  REAL, PARAMETER     :: Rgas=8.3144621      ! (J.K-1.mol-1)
  REAL, PARAMETER     :: kb=1.38E-23         ! Boltzman constant (J K-1)
  REAL, PARAMETER     :: pi=3.14159
  REAL, PARAMETER     :: multimass=1.66E-12  ! 1E12/avogadro: factor involved for molec/cm3 to microg/m3 conversion

  REAL, PARAMETER     :: Sair=4.E-19         ! sigma: collision cross section (m2)
  REAL, PARAMETER     :: Mair=28.6           ! mean molar mass of air
  REAL, PARAMETER     :: Cair=2.5E25         ! air concentration (molec/m3) @ room T and P

  REAL, PARAMETER     :: kappa=0.1           ! Hygroscopicity
  REAL, PARAMETER     :: densw=1.0           ! water density (g.cm-3)
  REAL, PARAMETER     :: Tgw=136.            ! glass transition temperature of water in K
  REAL, PARAMETER     :: kGT=2.5             ! Gordon Taylor constant
  REAL, PARAMETER     :: D=10.               ! fragility
  REAL, PARAMETER     :: alphas=1            ! surface accomodation coefficient
  REAL, PARAMETER     :: rad=5E-10           ! effective molecular radius (m)
 
  REAL, PARAMETER     :: Tgseed=163.         ! Tg of seed
  
  REAL, PARAMETER     :: lifetimegw=900.     ! 15 min

! Non Fixed parameters for properties calculation
  REAL                :: Rp0                 ! Initial radius of the seed particles 
  REAL                :: nvoc                ! Non Volatile Organic Compound (concentration)
  REAL                :: nvic                ! Non Volatile Inorganic Compound (concentration)
  REAL                :: denssoa             ! soa density (g.cm-3)
  REAL                :: densseed            ! seed density (g.cm-3)
  REAL                :: Mseed               ! seed molar mass
  
! Options for the simulation
  INTEGER, PARAMETER  :: fg_tgseed = 0       ! equal 1 to take into account the seed in Tg calculation
  INTEGER, PARAMETER  :: fg_alpha  = 1       ! equal 1 to calculate alpha with viscosity, fg=0 means alpha=1 for all species
  INTEGER, PARAMETER  :: fg_alpha_bulk  = 0  ! equal 1 to calculate alpha bulk, 0 for alpha to near surface
  LOGICAL, PARAMETER  :: fg_iterKgp =.FALSE. ! compute Keq for gas <-> particle at each time step (within the solver)

  INTEGER, PARAMETER  :: imtr=3              ! 0 mt approach, 1 diff, 2 coll, 3 Fuchs-Sutugin
  REAL, PARAMETER     :: winfac=1.0          ! simple multiplicative factor for gas->wall rate (by default, lifetime=15min, for 30min, put 0.5)
  REAL, PARAMETER     :: weqfac=1.0          ! simple multiplicative factor for gas->wall rate
  REAL, PARAMETER     :: ratfac=1.0          ! simple multiplicative factor gas-> aero rate

  LOGICAL             :: flag_ZA             ! flag to fix the zenithal angle
  REAL                :: fixedZA             ! fixed zenithal angle to set between 0 and 90 (put a value over 90 to perform simulations with no light)

! Options add for genoa - default not used
  INTEGER  :: fg_ref_conc = 0    ! if > 0: read concentrations from ref_conc_list_in
  INTEGER  :: fg_out_col = 0     ! if > 0: output only list species
  INTEGER  :: fg_output = 2      ! 2/1/0: with all/with only lbof/without(no print) regular outputs 
  CHARACTER (len=200) :: err_sps_list_in=""      ! Read by namelist in the format "A,B;C;D" to compute group concentrations and output in a file named outconc.dat
  CHARACTER (len=200) :: input_dir_chem="./"     ! folder contains chemistry related files
  CHARACTER (len=200) :: input_cond="indat.key"  ! File contains initial/constant condition
  CHARACTER (len=200) :: ref_conc_list_in=""       ! folder contains reference concentrations
  CHARACTER (len=200) :: phot_file="jfile.phot"  ! folder contains condition related files
  CHARACTER (len=200) :: output_dir="./"         ! folder to output results
  CHARACTER (len=200) :: init_conc_list_in=""    ! Read initial concentrations for species A 100,B 100; C for initID 1;2
  
END MODULE parameter_mod
