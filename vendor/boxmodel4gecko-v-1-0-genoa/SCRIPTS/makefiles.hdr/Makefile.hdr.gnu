#-*-makefile-*-
SHELL := /bin/bash

# select mode (MODE=DEVEL for development, anything else for regular compilation (e.g. PROD) 
ifeq ($(MODE),DEVEL)
     FFLAGS = -I${my_netcdf_fortran_inc} -g -O0 -fbounds-check -fbacktrace -finit-integer=-100000 -ffpe-trap=zero,overflow,invalid -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wunused-variable -Wuninitialized -fdefault-real-8 -mcmodel=large
else
     FFLAGS = -I${my_netcdf_fortran_inc} -fopenmp -fdefault-real-8 -mcmodel=large
endif

# Add for docker compatibility: ensure no PIE
FFLAGS += -fno-pie

#########################################################################
### In principle, you should not have to modify too many things below ...
# NetCDF config for various possible cases (nc-config/nf-config, hdf5, etc.)
INFCONFIG               =       $(shell [ -e $(my_netcdf_fortran_bin)/nf-config ]&& echo yes)
INCCONFIG               =       $(shell [ -e $(my_netcdf_c_bin)/nc-config ]&& echo yes)
ifeq ($(INFCONFIG),yes)
        NCCONFIG          =       $(my_netcdf_fortran_bin)/nf-config
else ifeq ($(INCCONFIG),yes)
        NCCONFIG          =       $(my_netcdf_c_bin)/nc-config
else
        NCCONFIG          =       none
endif
#
ifeq ($(NCCONFIG),none)
        NCFLIB          =       $(shell [ -e $(my_netcdf_fortran_lib)/libnetcdff.a ]&& echo twolibs)
        CULIB           =       $(shell nm $(my_netcdf_fortran_lib)/libnetcdf.a | grep -q curl && echo need_curl)
        HDLIB           =       $(shell [ -e $(my_netcdf_c_bin)/nc-config ]&& $(my_netcdf_c_bin)/nc-config --has-hdf5)
        CLIB            =       $(shell [ -e $(my_netcdf_c_bin)/nc-config ]&& $(my_netcdf_c_bin)/nc-config --libs)
        ifeq ($(NCFLIB),twolibs)
        CDFLIB          =       -lnetcdff -lnetcdf
        else
        CDFLIB          =       -lnetcdf
        endif
        ifeq ($(HDLIB),yes)
        CDFLIB1          =       $(CDFLIB) -lhdf5 -lhdf5_hl
        else
        CDFLIB1          =       $(CDFLIB)
        endif
        ifeq ($(CULIB),need_curl)
        CDFLIBS         =       $(CDFLIB1) -lcurl
        else
        CDFLIBS         =       $(CDFLIB1)
        endif
        NETCDFLIBS      =       $(CDFLIBS) -L${my_netcdf_fortran_lib} -L${my_hdf5_lib} ${CLIB}
else
        NETCDFLIBS_F    =       $(shell $(NCCONFIG) --flibs | gawk '{for (i=1;i<=NF;i++) if(substr($$i,1,2)=="-L" || substr($$i,1,2)=="-l") st=st" "$$i}END{print st}')
        NETCDFLIBS_C    =       $(shell $(my_netcdf_c_bin)/nc-config --libs | gawk '{for (i=1;i<=NF;i++) if(substr($$i,1,2)=="-L" || substr($$i,1,2)=="-l") st=st" "$$i}END{print st}')
        NETCDFLIBS      =       $(NETCDFLIBS_F) $(NETCDFLIBS_C)
endif




