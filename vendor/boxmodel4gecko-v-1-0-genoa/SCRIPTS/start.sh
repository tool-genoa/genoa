#!/bin/bash
source my_param.sh

ARCH=$(cat "${PATH_BOX}/SIMU/compilation_arch")
source ${PATH_BOX}/SCRIPTS/arch/${ARCH}

ulimit -s unlimited
export OMP_NUM_THREADS=${nb_proc}

################################################
# Copy needed files for simulation             #
################################################
out='xxxxx'
workdir='WORK'

[ -e $workdir ] || mkdir $workdir
cd $workdir
rm *	

cp ${PATH_BOX}/PROG/boxmod ./ || { echo 'program has not started' ; exit 1; }

cp ../concentrations.init ./indat.key || { echo 'program has not started' ; exit 1; }
cp ../simu.nml ./simu.nml || { echo 'program has not started' ; exit 1; }

cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero1.dat ./indat1.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero2.dat ./indat2.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero3.dat ./indat3.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero4.dat ./indat4.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero5.dat ./indat5.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero6.dat ./indat6.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero7.dat ./indat7.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero8.dat ./indat8.ro2 || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/COUNTERS/XXXXX/pero9.dat ./indat9.ro2 || { echo 'program has not started' ; exit 1; }

cp ${PATH_BOX}/CHEMDAT/xxxxx.bin ./indat.li || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/xxxxx.sat ./pvap.sat || { echo 'program has not started' ; exit 1; }
#cp ${PATH_BOX}/INPUT/xxxxx.dep ./vfile.dep || { echo 'program has not started' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/xxxxx.Tg ./Tg.dat || { echo 'program has not started, Tg.dat missing' ; exit 1; }
cp ${PATH_BOX}/CHEMDAT/xxxxx.mdv ./mdv.dat || { echo 'program has not started, mdv.dat missing' ; exit 1; }

cp ${PATH_BOX}/INPUT/solarlight.phot ./jfile.phot || { echo 'program has not started' ; exit 1; }
#cp ${PATH_BOX}/INPUT/blacklight.phot ./jfile.phot || { echo 'program has not started' ; exit 1; }

time ./boxmod

for file in `ls outdat.*` ; do
mv $file ../RESU/${out}.${file#*.}
done

for file in `ls fort.*` ; do
mv $file ../RESU/${out}.${file#*.}
done


