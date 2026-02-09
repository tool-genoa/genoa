#!/bin/bash

if [ $# -eq 0 ]; then
  echo "--ERROR-- this unix procedure must be run with some arguments"
  exit 1
fi

##################################################
# copy the files in the working directory
# stop if required files are not found
##################################################
cd ./WORK || exit 1

savestatus=0
cp ../../INTERP30/gecko2box ./intp.exe
if [ $? -ne 0 ]; then
  savestatus=$?
fi

cp ../"$1".mech ./indat.mech
if [ $? -ne 0 ]; then
  savestatus=$?
fi

if [ $savestatus -ne 0 ]; then
  echo
  echo "program has not been started"
  cd ../ || exit 1
  exit 1
fi

##################################################
#  run the interpretor
##################################################
echo "running program"
./intp.exe

##################################################
#  copy the output files
##################################################
for i in outdat.*; do
  mv "$i" ../"$1"."${i##*.}"
done

##################################################
#  clean and exit
##################################################
rm -f intp.exe indat.mech
cd ../
