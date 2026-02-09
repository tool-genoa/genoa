#!/bin/bash
#---------------------------------------------------------------------------------
#	Architecture file for compiling and running BOXMOD	
#	Specify path to libraries 
#---------------------------------------------------------------------------------
# NetCDF-Fortran

export my_netcdf_fortran_lib="/usr/lib64"
# export my_netcdf_fortran_inc="/usr/include"
export my_netcdf_fortran_inc="/usr/lib64/gfortran/modules"
export my_netcdf_fortran_bin="/usr/bin"

# NetCDF-C
export my_netcdf_c_lib="/usr/lib64"
export my_netcdf_c_inc="/usr/include"
export my_netcdf_c_bin="/usr/bin"

# HDF5
export my_hdf5_lib="/usr/lib64"

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

export my_hdr=Makefile.hdr.gnu
export FC=gfortran
