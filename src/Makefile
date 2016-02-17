SOURCE = ../source
BINDIR = ../bin
EXEC = $(BINDIR)/NanoTubeMaker
F90 = gfortran 
FFLAGS = -w -O3 -Wall


exe :  dir 
	$(F90) tubehbond.f90 pdb_hbond.f90 -o $(BINDIR)/TubeHBond
	$(F90) tubehbond.f90 metric_pathsample.f90  -o $(BINDIR)/PSHBond

dir :
	if [ ! -d $(BINDIR) ] ; then mkdir $(BINDIR) ;fi

clean :
	rm -rv   $(BINDIR)/* *.mod

