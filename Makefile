#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
## Description:   makefile for SAMRAI FAC Stokes
##
#########################################################################

SAMRAI	      = ../SAMRAI-hg
SRCDIR        = .
VPATH         = 
TESTTOOLS     = ../testtools
OBJECT        = ../SAMRAI-hg-build

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
	FACStokes/solveStokes.o \
	StokesFACOps/StokesFACOps.o \
	StokesFACOps/checkInputPatchDataIndices.o \
	StokesFACOps/computeCompositeResidualOnLevel.o \
	StokesFACOps/computeFluxOnPatch.o \
	StokesFACOps/computeResidualNorm.o \
	StokesFACOps/computeResidualOnPatch.o \
	StokesFACOps/computeVectorWeights.o \
	StokesFACOps/deallocateOperatorState.o \
	StokesFACOps/ewingFixFlux.o \
	StokesFACOps/finalizeCallback.o \
	StokesFACOps/initializeOperatorState.o \
	StokesFACOps/postprocessOneCycle.o \
	StokesFACOps/prolongErrorAndCorrect.o \
	StokesFACOps/redOrBlackSmoothingOnPatch.o \
	StokesFACOps/relax.o \
	StokesFACOps/restrictResidual.o \
	StokesFACOps/restrictSolution.o \
	StokesFACOps/smoothError.o \
	StokesFACOps/smoothErrorByRedBlack.o \
	StokesFACOps/solveCoarsestLevel.o \
	StokesFACOps/solveCoarsestLevel_HYPRE.o \
	StokesFACOps/xeqScheduleFluxCoarsen.o \
	StokesFACOps/xeqScheduleGhostFill.o \
	StokesFACOps/xeqScheduleGhostFillNoCoarse.o \
	StokesFACOps/xeqScheduleProlongation.o \
	StokesFACOps/xeqScheduleRRestriction.o \
	StokesFACOps/xeqScheduleURestriction.o \
	StokesFACSolver/StokesFACSolver.o \
	StokesFACSolver/StokesFACSolver_Destructor.o \
	StokesFACSolver/createVectorWrappers.o \
	StokesFACSolver/deallocateSolverState.o \
	StokesFACSolver/destroyVectorWrappers.o \
	StokesFACSolver/enableLogging.o \
	StokesFACSolver/getFromInput.o \
	StokesFACSolver/initializeSolverState.o \
	StokesFACSolver/initializeStatics.o \
	StokesFACSolver/setBcObject.o \
	StokesFACSolver/setBoundaries.o \
	StokesFACSolver/solveSystem.o \
	StokesHypreSolver.o StokesSpecifications.o
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
