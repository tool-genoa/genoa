#!/bin/bash
#---------------------------------------------------------------------------------
#	Architecture file for compiling and running BOXMOD	
#	Specify path to libraries 
#---------------------------------------------------------------------------------
# NetCDF-Fortran
export my_netcdf_fortran_lib="$CONDA_PREFIX/lib"
export my_netcdf_fortran_inc="$CONDA_PREFIX/include"
export my_netcdf_fortran_bin="$CONDA_PREFIX/bin"

# NetCDF-C
export my_netcdf_c_lib="$CONDA_PREFIX/lib"
export my_netcdf_c_inc="$CONDA_PREFIX/include"
export my_netcdf_c_bin="$CONDA_PREFIX/bin"

# HDF5
export my_hdf5_lib="$CONDA_PREFIX/lib"

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

export my_hdr=Makefile.hdr.gnu
export FC=gfortran
