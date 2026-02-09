#!/bin/bash
#---------------------------------------------------------------------------------
#	Architecture file for compiling and running BOXMOD	
#	Specify path to libraries 
#---------------------------------------------------------------------------------

# NetCDF Fortran library paths
export my_netcdf_fortran_lib='/usr/lib/aarch64-linux-gnu'
export my_netcdf_fortran_inc='/usr/include'
export my_netcdf_fortran_bin='/usr/bin'

# NetCDF C library paths
export my_netcdf_c_lib='/usr/local/lib'
export my_netcdf_c_inc='/usr/local/include'
export my_netcdf_c_bin='/usr/local/bin'

# HDF5 library paths
export my_hdf5_lib='/usr/lib/aarch64-linux-gnu/hdf5/serial'
export my_hdf5_inc='/usr/include/hdf5/serial'
export my_hdf5_bin='/usr/bin'

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

export my_hdr=Makefile.hdr.gnu
export FC=/usr/bin/gfortran
