MODULE mechdata
IMPLICIT NONE
INTEGER              :: lenspe          ! (max) length of species names
INTEGER              :: numsp           ! number of species in the mechanism
INTEGER              :: numre           ! number of reaction in the mechanism
INTEGER              :: numself         ! number of "self reaction"
INTEGER,ALLOCATABLE  :: idselfreac(:,:) ! ID of reaction 2X => pdct, SIZE=[numself,2]

! structure for storing products of each reactions
TYPE :: pdctmap
  INTEGER             :: npd       ! number of products (in current rx)
  INTEGER,ALLOCATABLE :: idpd(:)   ! ID of the products (in current rx)
  REAL,ALLOCATABLE    :: stoipd(:) ! stoichiometric coef. of each product
END TYPE
TYPE(pdctmap),ALLOCATABLE :: rxpdct(:) ! product list for rx [i]

! structure for storing reagents of each reactions
TYPE :: rgntmap
  INTEGER             :: nrg       ! number of reagents (in current rx)
  INTEGER,ALLOCATABLE :: idrg(:)   ! ID of the reagents (in current rx)
END TYPE
TYPE(rgntmap),ALLOCATABLE :: rxrgnt(:) ! reagent list for rx [i]


REAL                 :: kdilu     ! first order dilution (loss) rate constant (all species in all phases)
REAL                 :: kaerloss  ! first order loss rate constant for particles (e.g. aerosol wall loss)

END MODULE mechdata
