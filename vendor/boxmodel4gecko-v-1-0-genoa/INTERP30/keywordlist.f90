! Set the list of keyword available to design the mechanism.
MODULE keywordlist
  IMPLICIT NONE

! number of keyword available
  INTEGER, PARAMETER :: mxkeywd=21

! maximum length of a keyword string  
  INTEGER, PARAMETER :: mxlkey=9

! list of keywords     
  CHARACTER(LEN=mxlkey), DIMENSION(mxkeywd), PARAMETER :: keywdlist=  & 
    (/ 'HV       ', &  ! index  1 
       'OXYGEN   ', &  ! index  2  
       'TBODY    ', &  ! index  3
       'FALLOFF  ', &  ! index  4
       'EXTRA    ', &  ! index  5
       'MEPERO   ', &  ! index  6
       'PERO1    ', &  ! index  7
       'PERO2    ', &  ! index  8
       'PERO3    ', &  ! index  9
       'PERO4    ', &  ! index 10
       'PERO5    ', &  ! index 11 
       'PERO6    ', &  ! index 12 
       'PERO7    ', &  ! index 13 
       'PERO8    ', &  ! index 14 
       'PERO9    ', &  ! index 15
       'TABCF    ', &  ! index 16
       'AIN      ', &  ! index 17
       'AOU      ', &  ! index 18
       'WIN      ', &  ! index 19
       'WOU      ', &  ! index 20
       'ISOM     '  &  ! index 21
        /)                                                     

! set numbers to some indexes (check consistency with "keywdlist")
  INTEGER, PARAMETER :: hvidx=1     ! index for HV reaction
  INTEGER, PARAMETER :: o2idx=2     ! index for OXYGEN reaction
  INTEGER, PARAMETER :: tbidx=3     ! index for TBODY reaction
  INTEGER, PARAMETER :: foidx=4     ! index for FALLOFF reaction
  INTEGER, PARAMETER :: extraidx=5  ! index for EXTRA reaction
  INTEGER, PARAMETER :: meidx=6     ! index for MEPERO reaction
  INTEGER, PARAMETER :: p1idx=7     ! index for PERO1 reaction
  INTEGER, PARAMETER :: p2idx=8     ! index for PERO2 reaction
  INTEGER, PARAMETER :: p3idx=9     ! index for PERO3 reaction
  INTEGER, PARAMETER :: p4idx=10    ! index for PERO4 reaction
  INTEGER, PARAMETER :: p5idx=11    ! index for PERO5 reaction
  INTEGER, PARAMETER :: p6idx=12    ! index for PERO6 reaction
  INTEGER, PARAMETER :: p7idx=13    ! index for PERO7 reaction
  INTEGER, PARAMETER :: p8idx=14    ! index for PERO8 reaction
  INTEGER, PARAMETER :: p9idx=15    ! index for PERO9 reaction
  INTEGER, PARAMETER :: tabcfidx=16 ! index for TABCF reaction
  INTEGER, PARAMETER :: ainidx=17   ! index for AIN reaction
  INTEGER, PARAMETER :: aouidx=18   ! index for AOU reaction
  INTEGER, PARAMETER :: winidx=19   ! index for WIN reaction
  INTEGER, PARAMETER :: wouidx=20   ! index for WOU reaction
  INTEGER, PARAMETER :: isoidx=21   ! index for ISOM reaction
  
! Number of arguments (i.e. auxiliary info) expected for each keywords
  INTEGER, DIMENSION(mxkeywd), PARAMETER :: nauxarg= &
    (/ 2,   & ! index  1 - HV reaction (label, coef.)
       0,   & ! index  2 - OXYGEN reaction 
       0,   & ! index  3 - TBODY reaction 
       7,   & ! index  4 - FALLOFF reaction (3 stoi coef, 1 Fc, 3 others)
       7,   & ! index  5 - EXTRA reaction (depends on reaction)
       0,   & ! index  6 - MEPERO reaction
       0,   & ! index  7 - PERO1 reaction  
       0,   & ! index  8 - PERO2 reaction  
       0,   & ! index  9 - PERO3 reaction  
       0,   & ! index 10 - PERO4 reaction 
       0,   & ! index 11 - PERO5 reaction 
       0,   & ! index 12 - PERO6 reaction  
       0,   & ! index 13 - PERO7 reaction  
       0,   & ! index 14 - PERO8 reaction  
       0,   & ! index 15 - PERO9 reaction 
       2,   & ! index 16 - TABCF reaction (label, tab size)
       0,   & ! index 17 - AIN reaction
       0,   & ! index 18 - AOU reaction
       1,   & ! index 19 - WIN reaction
       1,   & ! index 20 - WOU reaction
       5    & ! index 21 - ISOM reaction
        /)    
  
END MODULE keywordlist
