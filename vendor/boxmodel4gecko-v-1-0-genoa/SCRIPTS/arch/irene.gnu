#!/bin/bash
module purge
module load flavor/buildcompiler/gcc/8
module load flavor/buildmpi/openmpi/4.0
module load gnu/8.3.0
module load mpi/openmpi/4.0.5
module switch feature/openmpi/net/ib/ucx-rc
module load flavor/hdf5/parallel
module load hdf5/1.8.20
module load netcdf-fortran/4.4.4

export my_netcdf_fortran_lib=${NETCDFFORTRAN_LIBDIR}
export my_netcdf_fortran_inc=${NETCDFFORTRAN_INCDIR}
export my_netcdf_fortran_bin=${NETCDFFORTRAN_EXEDIR}
export my_hdf5_lib=${HDF5_LIBDIR}
export my_netcdf_c_lib=${NETCDFC_LIBDIR}
export my_netcdf_c_inc=${NETCDFC_INCDIR}
export my_netcdf_c_bin=${NETCDFC_EXEDIR}

export LD_LIBRARY_PATH=${my_netcdf_c_lib}:${my_netcdf_fortran_lib}:$LD_LIBRARY_PATH

export my_hdr=Makefile.hdr.gnu
export FC=${C_GNU_ROOT}/bin/gfortran

