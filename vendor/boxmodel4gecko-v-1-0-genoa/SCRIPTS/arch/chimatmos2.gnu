#!/bin/bash
#---------------------------------------------------------------------------------
#	Architecture file for compiling and running BOXMOD	
#	Specify path to libraries 
#---------------------------------------------------------------------------------
module purge
module load slurm
module load openmpi/gnu/64/4.0.5
module load hdf5/parallel/gnu/1.12.0
module load netcdf/netcdf-c/gnu/4.7.4
module load netcdf/netcdf-f/gnu/4.5.3

export my_netcdf_fortran_lib=${NETCDFFORTRAN_LIBDIR}
export my_netcdf_fortran_inc=${NETCDFFORTRAN_INCDIR}
export my_netcdf_fortran_bin=${NETCDFFORTRAN_EXEDIR}
export my_hdf5_lib=${HDF5_LIBDIR}
export my_netcdf_c_lib=${NETCDFC_LIBDIR}
export my_netcdf_c_inc=${NETCDFC_INCDIR}
export my_netcdf_c_bin=${NETCDFC_LEXEDIR}

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${my_netcdf_c_lib}:${my_netcdf_fortran_lib}

export my_hdr=Makefile.hdr.gnu
export FC=gfortran
