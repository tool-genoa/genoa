#!/bin/bash
#---------------------------------------------------------------------------------
#	Architecture file for compiling and running BOXMOD	
#	Specify path to libraries 
#---------------------------------------------------------------------------------
module load XXXXX

export my_netcdf_fortran_lib='XXXX'
export my_netcdf_fortran_inc='XXXX'
export my_netcdf_fortran_bin='XXXX'
export my_hdf5_lib='XXXX'
export my_netcdf_c_lib='XXXX'
export my_netcdf_c_inc='XXXX'
export my_netcdf_c_bin='XXXX'

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

export my_hdr=Makefile.hdr.gnu
export FC=gfortran
