#!/bin/bash
module purge
module load gnu/11.3.0
module load openmpi/4.1.4
module load hdf5/1.12.1
module load netcdf-c/4.9.0  
module load netcdf-f/4.6.0

export my_netcdf_fortran_lib=/opt/spack_soft/netcdf-fortran/4.6.0/gcc-11.3.0-sgxjb6/lib
export my_netcdf_fortran_inc=/opt/spack_soft/netcdf-fortran/4.6.0/gcc-11.3.0-sgxjb6/include
export my_netcdf_fortran_bin=/opt/spack_soft/netcdf-fortran/4.6.0/gcc-11.3.0-sgxjb6/bin
export my_hdf5_lib=/opt/spack_soft/hdf5/1.12.1/gcc-11.3.0-pjyyu4/lib
export my_netcdf_c_lib=/opt/spack_soft/netcdf-c/4.9.0/gcc-11.3.0-pzomjb/lib
export my_netcdf_c_inc=/opt/spack_soft/netcdf-c/4.9.0/gcc-11.3.0-pzomjb/include
export my_netcdf_c_bin=/opt/spack_soft/netcdf-c/4.9.0/gcc-11.3.0-pzomjb/bin

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

export my_hdr=Makefile.hdr.gnu
export FC=gfortran

