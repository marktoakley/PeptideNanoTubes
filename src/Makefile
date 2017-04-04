SOURCE = ../source
BINDIR = ../bin
F90 = gfortran 
FFLAGS = -w -O3 -Wall

exe :  dir 
	$(F90) $(FFLAGS) tube_module.f90 metric_pdb.f90 -o $(BINDIR)/PDBHBond
	$(F90) $(FFLAGS) tube_module.f90 metric_pathsample.f90  -o $(BINDIR)/PSHBond

dir :
	if [ ! -d $(BINDIR) ] ; then mkdir $(BINDIR) ;fi

clean :
	rm -rv   $(BINDIR)/* *.mod

