#-*-makefile-*-

include Makefile.hdr

$(BOXDIR)/PROG/$(prog): $(BOXDIR)/PROG/$(prog).f90 $(BOXDIR)/OBJ/boxlib.a
	cp $(BOXDIR)/PROG/$(prog).f90 ./
	cp $(BOXDIR)/OBJ/boxlib.a ./
	cp $(BOXDIR)/OBJ/*.mod ./
	$(FC) $(prog).f90 $(FFLAGS) boxlib.a $(NETCDFLIBS) -o ../../PROG/$(prog)

