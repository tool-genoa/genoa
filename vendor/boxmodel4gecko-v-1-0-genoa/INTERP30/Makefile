#-*-makefile-*-

include Makefile.hdr

# define modules (to be produced 1st)
MODULE=sortstring.mod intptool.mod keywordlist.mod speciestool.mod \
       auxtool.mod reactiontool.mod mecatool.mod
OBJECT=intptool.o keywordlist.o speciestool.o auxtool.o reactiontool.o\
       mecatool.o gecko2box.o

# 
gecko2box : $(MODULE) $(OBJECT)
	$(FC) $(FFLAGS) *.o -o gecko2box
	rm *.o *.mod
	
# make the modules 
%.mod : %.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.mod $*.o

# make the objects 
%.o : %.f90
	$(FC) $(FFLAGS) -c $<

# make the modules first (use update of the main program)
#gecko2box.o :: gecko2box.f90 $(MODULE)
#	$(FC) -c $(MODULE)
#	$(FC) -c gecko2box.f90

clean :
	rm -f *.o *.mod *.a

