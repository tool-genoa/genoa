#!/bin/bash
ulimit -s unlimited

compil_mode=PROD
arch_file=''
export BOXDIR=`pwd`

if [ $# -eq 0 ]
then
  echo 'argument missing : type --help'
  exit
fi

# Deal with prompt
while (($# > 0))
do
    case $1 in
        "-h"|"--h"|"--help"|"-help")
            echo "build.sh - Compile Boxmodel for GECKO"
            echo "build.sh [options]"
            echo "      [ -h | --h | -help | --help ] : this help message"
            echo "      [--arch] : specify file with your architecture paths"
            echo "      [--dev | --devel] : compilation in development/debug mode (default : production mode)"
            exit ;;
        "--devel"|"--dev") compil_mode=DEVEL ; shift ;;
        "--prod")          compil_mode=PROD  ; shift ;;
        "--arch") arch_file=$2 ; shift ; shift ;;
        *) echo "Wrong restart option $1"; exit 1 ;;
    esac
done

echo '================================='
echo '  Loading environment variables '
echo '================================='
pwdir=`pwd`
cd SCRIPTS/arch/
if [ -e "${arch_file}" ]; then
  source ${arch_file}
  export MODE=${compil_mode}
  export hdr_file=${pwdir}/SCRIPTS/makefiles.hdr/$my_hdr
  [ -s ${hdr_file} ]|| { echo "non existent $hdr_file makefile header ; exiting." ; exit 1; }
  cd ../..
  echo ${arch_file} > SIMU/compilation_arch
else
    echo "arch file doesn't exist in SCRIPTS/arch/ - see ./build.sh --help"
    ls *
    exit
fi

echo '============================'
echo '  Compiling files  '
echo '============================'
echo "INTERP30"
cd INTERP30/
ln -sf $hdr_file Makefile.hdr
make clean
make || exit 1
echo ""
echo "OBJ"
echo ""
cd ../OBJ/
ln -sf $hdr_file Makefile.hdr
make clean
make || exit 1
echo ""
echo "PROG"
echo ""
cd ../PROG/
./callmake.sh boxmod || exit 1
