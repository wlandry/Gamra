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
        source       = ['src/main.C',
                        'src/FACStokes/FACStokes.C',
                        'src/FACStokes/fix_viscosity.C',
                        'src/FACStokes/applyGradientDetector.C',
                        'src/FACStokes/initializeLevelData.C',
                        'src/FACStokes/packDerivedDataIntoDoubleBuffer.C',
                        'src/FACStokes/resetHierarchyConfiguration.C',
                        'src/FACStokes/setupPlotter.C',
                        'src/FACStokes/solveStokes.C',
                        'src/P_Refine.C',
                        'src/P_MDPI_Refine.C',
                        'src/V_Refine.C',
                        'src/V_Coarsen/coarsen_2D.C',
                        'src/V_Coarsen/coarsen_3D.C',
                        'src/P_Boundary_Refine/refine.C',
                        'src/P_Boundary_Refine/Update_P_2D.C',
                        'src/P_Boundary_Refine/Update_P_3D.C',
                        'src/V_Boundary_Refine/refine.C',
                        'src/V_Boundary_Refine/Update_V_2D.C',
                        'src/V_Coarsen_Patch_Strategy/postprocessCoarsen_2D.C',
                        'src/V_Coarsen_Patch_Strategy/postprocessCoarsen_3D.C',
                        'src/Cell_Viscosity_Coarsen.C',
                        'src/Edge_Viscosity_Coarsen.C',
                        'src/set_boundary.C',
                        'src/StokesFACOps/StokesFACOps.C',
                        'src/StokesFACOps/computeCompositeResidualOnLevel.C',
                        'src/StokesFACOps/residual_2D.C',
                        'src/StokesFACOps/residual_3D.C',
                        'src/StokesFACOps/computeResidualNorm.C',
                        'src/StokesFACOps/computeVectorWeights.C',
                        'src/StokesFACOps/deallocateOperatorState.C',
                        'src/StokesFACOps/finalizeCallback.C',
                        'src/StokesFACOps/initializeOperatorState.C',
                        'src/StokesFACOps/postprocessOneCycle.C',
                        'src/StokesFACOps/prolongErrorAndCorrect.C',
                        'src/StokesFACOps/restrictResidual.C',
                        'src/StokesFACOps/restrictSolution.C',
                        'src/StokesFACOps/smoothError.C',
                        'src/StokesFACOps/smooth_Tackley_2D.C',
                        'src/StokesFACOps/smooth_Tackley_3D.C',
                        'src/StokesFACOps/smooth_Gerya.C',
                        'src/StokesFACOps/set_boundaries.C',
                        'src/StokesFACOps/solveCoarsestLevel.C',
                        'src/StokesFACOps/solveCoarsestLevel_HYPRE.C',
                        'src/StokesFACOps/smooth_V_2D.C',
                        'src/StokesFACOps/smooth_V_3D.C',
                        'src/StokesFACOps/xeqScheduleGhostFill.C',
                        'src/StokesFACOps/xeqScheduleGhostFillNoCoarse.C',
                        'src/StokesFACOps/xeqScheduleProlongation.C',
                        'src/StokesFACOps/xeqScheduleRRestriction.C',
                        'src/StokesFACOps/xeqScheduleURestriction.C',
                        'src/StokesFACSolver/StokesFACSolver.C',
                        'src/StokesFACSolver/StokesFACSolver_Destructor.C',
                        'src/StokesFACSolver/createVectorWrappers.C',
                        'src/StokesFACSolver/deallocateSolverState.C',
                        'src/StokesFACSolver/destroyVectorWrappers.C',
                        'src/StokesFACSolver/enableLogging.C',
                        'src/StokesFACSolver/getFromInput.C',
                        'src/StokesFACSolver/initializeSolverState.C',
                        'src/StokesFACSolver/initializeStatics.C',
                        'src/StokesFACSolver/setBcObject.C',
                        'src/StokesFACSolver/setBoundaries.C',
                        'src/StokesFACSolver/solveSystem.C',
                        'src/StokesHypreSolver.C'],

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
        includes = ['src','../SAMRAI-hg-build/include',
                    '../SAMRAI-hg/source',
                    '/usr/lib/petsc/include',
                    '/usr/lib/petsc/bmake/linux-gnu-c-opt',
                    '/usr/lib/petsc/src/vec',
                    '/usr/include']
        )
