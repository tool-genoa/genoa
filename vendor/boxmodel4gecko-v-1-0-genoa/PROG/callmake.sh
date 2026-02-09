#/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[1;34m'
NC='\033[0m' # No Color

if [ "$#" == 0 ] ; then 
  echo -e "${RED}--ERROR-- this unix procedure must be run with one argument (name of fortran file without extension) ${NC}"
  exit 1
fi

echo -e "${BLUE} ============================ ${NC}"
echo -e "${BLUE}   Making $1 executable  ${NC}"
echo -e "${BLUE} ============================ ${NC}"

cd ./WORK
rm -f *

ln -sf $hdr_file Makefile.hdr

cp ../prog.makefile ./
make -f ./prog.makefile prog=$1 
if [ $? != 0 ] ;  then
  echo -e "${RED}-- COMPILATION ERROR --${NC}"
  exit 1
else
  echo -e "${GREEN}-- NO ERROR --${NC}"
fi

#rm *.mod *.o *.a *.makefile 
cd ..

exit

