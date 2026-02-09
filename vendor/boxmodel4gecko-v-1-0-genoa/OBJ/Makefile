#-*-makefile-*-

include Makefile.hdr

# set the path to search for *f90 files
vpath %.f90 ../LIBSRC 

MODULE=boxtool.mod mechdata.mod ncdf_saving.mod lpmaptool.mod parameter_mod.mod ro2maptool.mod mechdata.mod readchem.mod stringdata.mod \
       cubicspline.mod  initnbound.mod frozenstuff.mod ratetool.mod iteration.mod solver.mod zenith.mod phototool.mod\
       transfertool.mod viscositytool.mod diffusivitytool.mod vaporpressure.mod wrtboftool.mod

OBJECT=boxtool.o ncdf_saving.o lpmaptool.o parameter_mod.o ro2maptool.o mechdata.mod readchem.o stringdata.o \
       cubicspline.o  initnbound.o frozenstuff.o  ratetool.o iteration.o solver.o zenith.o phototool.o\
       viscositytool.o diffusivitytool.o vaporpressure.o transfertool.o wrtboftool.o

# object files needed for the boxmodel
all: $(MODULE) $(OBJECT)

# make the modules 
%.mod : %.f90
	$(FC) -c $(FFLAGS) $(NETCDFLIBS) -s $<
	ar -r boxlib.a $*.o
	touch $*.mod $*.o

# make the objects 
%.o : %.f90 
	$(FC) -c $(FFLAGS) $(NETCDFLIBS) -s $<
	ar -r boxlib.a $@
	
clean :
	rm *.o *.mod *.a

