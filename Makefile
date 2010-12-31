#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
## Description:   makefile for SAMRAI FAC Stokes
##
#########################################################################

SAMRAI	      = ../AMR/SAMRAI-v3.1.0-beta
SRCDIR        = .
VPATH         = 
TESTTOOLS     = ../testtools
OBJECT        = ../AMR/SAMRAI-v3.1.0-beta-build

default:      main

include $(OBJECT)/config/Makefile.config

CPPFLAGS_EXTRA= -DTESTING=0

NUM_TESTS = 2

TEST_NPROCS = 0,2

CXX_OBJS      = main.o FACStokes/FACStokes.o \
	FACStokes/initializeLevelData.o \
	FACStokes/packDerivedDataIntoDoubleBuffer.o \
	FACStokes/resetHierarchyConfiguration.o \
	FACStokes/setupPlotter.o \
	FACStokes/solveStokes.o StokesFACOps.o \
	StokesHypreSolver.o StokesSpecifications.o StokesFACSolver.o
F_OBJS      = facpoisson2d.o facpoisson3d.o

main:	$(CXX_OBJS) $(F_OBJS) $(LIBSAMRAIDEPEND)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXX_OBJS) $(F_OBJS) \
	$(LIBSAMRAI3D) $(LIBSAMRAI3D) $(LIBSAMRAI) $(LDLIBS) -o main

check:
	$(MAKE) check2d
	$(MAKE) check3d

check2d:	main
	@for i in test_inputs/*2d.input ; do	\
	  $(OBJECT)/config/serpa-run $(TEST_NPROCS) \
		./main $${i};	\
	done

check3d:	main
	@for i in test_inputs/*3d.input ; do	\
	  $(OBJECT)/config/serpa-run $(TEST_NPROCS) \
		./main $${i};	\
	done

checkcompile: main

checktest:
	rm -f makecheck.logfile
	$(MAKE) check 2>&1 | $(TEE) makecheck.logfile
	$(TESTTOOLS)/testcount.sh $(TEST_NPROCS) $(NUM_TESTS) makecheck.logfile
	rm -f makecheck.logfile

examples2d: main
	@for i in $(SRCDIR)/example_inputs/*.2d.input ; do	\
	  $(OBJECT)/config/serpa-run $(TEST_NPROCS) \
		./main $${i};	\
	done

examples3d: main
	@for i in $(SRCDIR)/example_inputs/*.3d.input ; do	\
	  $(OBJECT)/config/serpa-run $(TEST_NPROCS) \
		./main $${i};	\
	done

examples:
	$(MAKE) examples2d
	$(MAKE) examples3d

clean-check:
	$(SAMCLEAN)

clean:		clean-check
	$(RM) main *.f *.o */*.o

redo:
	$(RM) core main

include $(SRCDIR)/Makefile.depend

FORTRAN       = $(SRCDIR)/fortran
M4DIRS          = -DFORTDIR=$(FORTRAN) $(SAMRAI_M4_FLAGS)

facpoisson2d.o:	$(FORTRAN)/facpoisson2d.m4
	$(M4) $(M4DIRS) $(FORTRAN)/facpoisson2d.m4 > facpoisson2d.f
	$(F77) $(FFLAGS) -c facpoisson2d.f -o $@

facpoisson3d.o:	$(FORTRAN)/facpoisson3d.m4
	$(M4) $(M4DIRS) $(FORTRAN)/facpoisson3d.m4 > facpoisson3d.f
	$(F77) $(FFLAGS) -c facpoisson3d.f -o $@
