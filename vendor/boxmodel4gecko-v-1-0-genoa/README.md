## GECKO-A box model wiki
A detailed wiki is avalaible on GECKO-a box model gitlab at the following URL : https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/boxmodel4gecko/-/wikis/home



## Structure of the repository

##### CHEMDAT

Contains the chemical mechanims, saturation vapore pressure files, ...

##### INPUT

Contains input files such as photolysis tables.

##### INTERP

Contains the interpretor used to translaste chemical mechanisms into binary file used by the boxmodel.

##### LIBSRC

Contains all the fortran subroutines.

##### OBJ

This is the compilation directory.

##### SIMU

This the folder where the simulations are.

##### PROG

Contains main boxmodel program.


## Compile the code
Requirements for compilation: fortran compiler (e.g. gfortran), netcdf library for output files.
To compile the GECKO-A box model, you need to use *build.sh* script. this script must be launched with options


[ -h | --h | -help | --help ] : this help message

[--arch] : specify file with your architecture paths

[--dev | --devel] : compilation in development/debug mode (default : production mode)


At least --arch need to be specified with the name of the arch file containing the path to your netcdf library. When the repo is cloned, you need to create your own file corresponding to your machine:

go into scripts/arch folder.
copy template_for_other_machine.gnu (or .intel)  to name_of_your_machine.gnu (or .intel).
edit the file and add the paths to your netcdf library.

Then you can go back at the root of the repo and run:
./build.sh --arch name_of_your_machine.gnu


## prepare a simulation
To import GECKO files, run the following script:

./prepare_simu.sh --import_mech nameforyourmechanism --geckooutdir PATH/TO/GECKO/RUN/OUT

Replace *nameforyourmechanism* with the name you want for files in the box model.

If the script runs successfully, it will copy files into the *BOXMODEL/CHEMDAT* folder and create the simulation folder in *BOXMODEL/SIMU*.


## launch a simulation
To run a simulation, go to the simulation folder and run the *start.sh* script.
