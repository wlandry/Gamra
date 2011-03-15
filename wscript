# top = '.'
# out = 'build'

def options(opt):
    opt.load('compiler_cxx')

def configure(conf):
    conf.load('compiler_cxx')
def build(bld):
    # bld.recurse('src/libtcod-1.5.0')
    bld.program(
        features     = ['cxx','cprogram'],
        source       = ['main.C',
                        'FACStokes/FACStokes.C',
                        'FACStokes/fix_viscosity.C',
                        'FACStokes/initializeLevelData.C',
                        'FACStokes/packDerivedDataIntoDoubleBuffer.C',
                        'FACStokes/resetHierarchyConfiguration.C',
                        'FACStokes/setupPlotter.C',
                        'FACStokes/solveStokes.C',
                        'P_Refine.C',
                        'P_MDPI_Refine.C',
                        'V_Refine.C',
                        'V_Coarsen.C',
                        'P_Boundary_Refine.C',
                        'V_Boundary_Refine/refine.C',
                        'V_Boundary_Refine/Update_V.C',
                        'P_Refine_Patch_Strategy.C',
                        'V_Refine_Patch_Strategy.C',
                        'V_Coarsen_Patch_Strategy.C',
                        'Cell_Viscosity_Coarsen.C',
                        'Edge_Viscosity_Coarsen.C',
                        'set_V_boundary.C',
                        'StokesFACOps/StokesFACOps.C',
                        'StokesFACOps/checkInputPatchDataIndices.C',
                        'StokesFACOps/computeCompositeResidualOnLevel.C',
                        'StokesFACOps/residual_2D.C',
                        'StokesFACOps/residual_3D.C',
                        'StokesFACOps/computeFluxOnPatch.C',
                        'StokesFACOps/computeResidualNorm.C',
                        'StokesFACOps/computeVectorWeights.C',
                        'StokesFACOps/deallocateOperatorState.C',
                        'StokesFACOps/ewingFixFlux.C',
                        'StokesFACOps/finalizeCallback.C',
                        'StokesFACOps/initializeOperatorState.C',
                        'StokesFACOps/postprocessOneCycle.C',
                        'StokesFACOps/prolongErrorAndCorrect.C',
                        'StokesFACOps/restrictResidual.C',
                        'StokesFACOps/restrictSolution.C',
                        'StokesFACOps/smoothError.C',
                        'StokesFACOps/smooth_Tackley.C',
                        'StokesFACOps/smooth_Gerya.C',
                        'StokesFACOps/set_boundaries.C',
                        'StokesFACOps/solveCoarsestLevel.C',
                        'StokesFACOps/solveCoarsestLevel_HYPRE.C',
                        'StokesFACOps/smooth_V.C',
                        'StokesFACOps/xeqScheduleFluxCoarsen.C',
                        'StokesFACOps/xeqScheduleGhostFill.C',
                        'StokesFACOps/xeqScheduleGhostFillNoCoarse.C',
                        'StokesFACOps/xeqScheduleProlongation.C',
                        'StokesFACOps/xeqScheduleRRestriction.C',
                        'StokesFACOps/xeqScheduleURestriction.C',
                        'StokesFACSolver/StokesFACSolver.C',
                        'StokesFACSolver/StokesFACSolver_Destructor.C',
                        'StokesFACSolver/createVectorWrappers.C',
                        'StokesFACSolver/deallocateSolverState.C',
                        'StokesFACSolver/destroyVectorWrappers.C',
                        'StokesFACSolver/enableLogging.C',
                        'StokesFACSolver/getFromInput.C',
                        'StokesFACSolver/initializeSolverState.C',
                        'StokesFACSolver/initializeStatics.C',
                        'StokesFACSolver/setBcObject.C',
                        'StokesFACSolver/setBoundaries.C',
                        'StokesFACSolver/solveSystem.C',
                        'StokesHypreSolver.C',
                        'StokesSpecifications.C'],

        target       = 'gamr',
        # cxxflags      = ['-std=c++0x','-g','-D_GLIBCXX_DEBUG'],
        cxxflags      = ['-g', '-Wall', '-Wextra', '-Wconversion',
                         '-DTESTING=0'],
        lib          = ['SAMRAI_appu', 'SAMRAI_algs', 'SAMRAI_solv',
                        'SAMRAI_geom', 'SAMRAI_mesh', 'SAMRAI_math',
                        'SAMRAI_pdat', 'SAMRAI_xfer', 'SAMRAI_hier',
                        'SAMRAI_tbox',   'petsc', 'HYPRE',
                        'HYPRE_sstruct_ls', 'HYPRE_sstruct_mv',
                        'HYPRE_struct_ls', 'HYPRE_struct_mv',
                        'HYPRE_parcsr_ls',
                        'HYPRE_DistributedMatrixPilutSolver',
                        'HYPRE_ParaSails', 'HYPRE_Euclid',
                        'HYPRE_MatrixMatrix', 'HYPRE_DistributedMatrix',
                        'HYPRE_IJ_mv', 'HYPRE_parcsr_mv', 'HYPRE_seq_mv',
                        'HYPRE_krylov', 'HYPRE_utilities', 'hdf5',
                        'dl', 'm', 'gfortranbegin', 'gfortran', 'm',
                        'gfortranbegin', 'gfortran', 'm', 'gfortranbegin',
                        'gfortran', 'm', 'sundials_cvode', 'sundials_kinsol'],
        libpath      = ['../../SAMRAI-hg-build/lib',
                        '/usr/lib/petsc/linux-gnu-c-opt/lib'],
        includes = ['.','../SAMRAI-hg-build/include',
                    '../SAMRAI-hg/source',
                    '/usr/lib/petsc/include',
                    '/usr/lib/petsc/bmake/linux-gnu-c-opt',
                    '/usr/lib/petsc/src/vec',
                    '/usr/include']
        )
