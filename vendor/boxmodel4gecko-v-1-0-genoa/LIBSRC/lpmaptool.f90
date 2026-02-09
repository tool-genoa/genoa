MODULE lpmaptool

! structure for storing loss/production reactions of each species 
  TYPE :: spec_reac_map
    INTEGER             :: nloss   ! number of loss term
    INTEGER             :: nprod   ! number of production term
    INTEGER,ALLOCATABLE :: idl(:)  ! ID of loss reactions
    INTEGER,ALLOCATABLE :: idp(:)  ! ID of production reactions
    REAL,ALLOCATABLE    :: stpd(:) ! prod related stoichiometric coef.
  END TYPE
  TYPE(spec_reac_map),ALLOCATABLE :: lpmap(:) ! loss_production map (size=numsp)

CONTAINS
!=======================================================================
! PURPOSE: map the reactions contributing to the loss and production of 
! each species. Allocate memory to the lpmap(:)%idl, lpmap(:)%idp and 
! lpmap(:)%stpd fields. 
!=======================================================================
SUBROUTINE make_lpmap(numsp,numre)
  USE mechdata, ONLY: rxpdct,rxrgnt
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: numsp          ! number of species in the mechanism
  INTEGER,INTENT(IN) :: numre          ! number of reaction in the mechanism
  
  INTEGER :: lossterm(numsp),prodterm(numsp)
  INTEGER :: i,ire,ispec,nloss,nprod

  lossterm(:)=0 ; prodterm(:)=0

! count the number of loss & production rxn (i.e. term) for each species 
  DO ire=1,numre      
    DO i=1,rxrgnt(ire)%nrg
      ispec = rxrgnt(ire)%idrg(i)
      lossterm(ispec) = lossterm(ispec) + 1
    ENDDO
    DO i=1,rxpdct(ire)%npd
      ispec = rxpdct(ire)%idpd(i)
      prodterm(ispec) = prodterm(ispec) + 1
    ENDDO
  ENDDO

! allocate the space required for loss/production map
  DO ispec = 1, numsp
    ALLOCATE(lpmap(ispec)%idl(lossterm(ispec)),      &
             lpmap(ispec)%idp(prodterm(ispec)),      &
             lpmap(ispec)%stpd(prodterm(ispec)) )
  ENDDO

! fill lpmap with the mechanism data
  lpmap%nloss = 0 ; lpmap%nprod = 0  
  DO ire=1,numre
    DO i=1,rxrgnt(ire)%nrg             ! loss side of the reaction
      ispec = rxrgnt(ire)%idrg(i)
      nloss = lpmap(ispec)%nloss + 1
      lpmap(ispec)%nloss = nloss
      lpmap(ispec)%idl(nloss) = ire
    ENDDO
    DO i=1,rxpdct(ire)%npd             ! production side of the reaction
      ispec = rxpdct(ire)%idpd(i)
      nprod = lpmap(ispec)%nprod + 1
      lpmap(ispec)%nprod = nprod
      lpmap(ispec)%idp(nprod) = ire
      lpmap(ispec)%stpd(nprod) = rxpdct(ire)%stoipd(i)
    ENDDO
  ENDDO
    
END SUBROUTINE make_lpmap

END MODULE lpmaptool

