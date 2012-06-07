# top = '.'
# out = 'build'

def options(opt):
    opt.load('compiler_cxx')

def configure(conf):
    conf.load('compiler_cxx')
def build(bld):
    bld.program(
        features     = ['cxx','cprogram'],
        source       = ['src/main.C',
                        'src/Stokes/FAC/FAC.C',
                        'src/Stokes/FAC/fix_viscosity.C',
                        'src/Stokes/FAC/applyGradientDetector.C',
                        'src/Stokes/FAC/initializeLevelData.C',
                        'src/Stokes/FAC/packDerivedDataIntoDoubleBuffer.C',
                        'src/Stokes/FAC/resetHierarchyConfiguration.C',
                        'src/Stokes/FAC/setupPlotter.C',
                        'src/Stokes/FAC/solve.C',
                        'src/Stokes/P_Refine.C',
                        'src/Stokes/V_Refine.C',
                        'src/Stokes/Resid_Coarsen.C',
                        'src/Stokes/V_Coarsen/coarsen_2D.C',
                        'src/Stokes/V_Coarsen/coarsen_3D.C',
                        'src/Stokes/P_Boundary_Refine/refine.C',
                        'src/Stokes/P_Boundary_Refine/Update_P_2D.C',
                        'src/Stokes/P_Boundary_Refine/Update_P_3D.C',
                        'src/Stokes/V_Boundary_Refine/refine.C',
                        'src/Stokes/V_Boundary_Refine/Update_V_2D.C',
                        'src/Stokes/V_Boundary_Refine/Update_V_3D.C',
                        'src/Stokes/V_Coarsen_Patch_Strategy/postprocessCoarsen_2D.C',
                        'src/Stokes/V_Coarsen_Patch_Strategy/postprocessCoarsen_3D.C',
                        'src/Stokes/set_boundary.C',
                        'src/Stokes/FACOps/FACOps.C',
                        'src/Stokes/FACOps/computeCompositeResidualOnLevel.C',
                        'src/Stokes/FACOps/residual_2D.C',
                        'src/Stokes/FACOps/residual_3D.C',
                        'src/Stokes/FACOps/computeResidualNorm.C',
                        'src/Stokes/FACOps/computeVectorWeights.C',
                        'src/Stokes/FACOps/deallocateOperatorState.C',
                        'src/Stokes/FACOps/finalizeCallback.C',
                        'src/Stokes/FACOps/initializeOperatorState.C',
                        'src/Stokes/FACOps/postprocessOneCycle.C',
                        'src/Stokes/FACOps/prolongErrorAndCorrect.C',
                        'src/Stokes/FACOps/restrictResidual.C',
                        'src/Stokes/FACOps/restrictSolution.C',
                        'src/Stokes/FACOps/smoothError.C',
                        'src/Stokes/FACOps/smooth_Tackley_2D.C',
                        'src/Stokes/FACOps/smooth_Tackley_3D.C',
                        'src/Stokes/FACOps/smooth_Gerya.C',
                        'src/Stokes/FACOps/set_boundaries.C',
                        'src/Stokes/FACOps/solveCoarsestLevel.C',
                        'src/Stokes/FACOps/solveCoarsestLevel_HYPRE.C',
                        'src/Stokes/FACOps/smooth_V_2D.C',
                        'src/Stokes/FACOps/smooth_V_3D.C',
                        'src/Stokes/FACOps/xeqScheduleGhostFill.C',
                        'src/Stokes/FACOps/xeqScheduleGhostFillNoCoarse.C',
                        'src/Stokes/FACOps/xeqScheduleProlongation.C',
                        'src/Stokes/FACOps/xeqScheduleRRestriction.C',
                        'src/Stokes/FACOps/xeqScheduleURestriction.C',
                        'src/Stokes/FACSolver/FACSolver.C',
                        'src/Stokes/FACSolver/FACSolver_Destructor.C',
                        'src/Stokes/FACSolver/createVectorWrappers.C',
                        'src/Stokes/FACSolver/deallocateSolverState.C',
                        'src/Stokes/FACSolver/destroyVectorWrappers.C',
                        'src/Stokes/FACSolver/enableLogging.C',
                        'src/Stokes/FACSolver/getFromInput.C',
                        'src/Stokes/FACSolver/initializeSolverState.C',
                        'src/Stokes/FACSolver/initializeStatics.C',
                        'src/Stokes/FACSolver/setBcObject.C',
                        'src/Stokes/FACSolver/setBoundaries.C',
                        'src/Stokes/FACSolver/solveSystem.C',
                        'src/Stokes/HypreSolver.C',
                        'src/Elastic/FAC/FAC.C',
                        'src/Elastic/FAC/fix_moduli.C',
                        'src/Elastic/FAC/applyGradientDetector.C',
                        'src/Elastic/FAC/initializeLevelData.C',
                        'src/Elastic/FAC/packDerivedDataIntoDoubleBuffer.C',
                        'src/Elastic/FAC/resetHierarchyConfiguration.C',
                        'src/Elastic/FAC/setupPlotter.C',
                        'src/Elastic/FAC/solve.C',
                        'src/Elastic/P_Refine.C',
                        'src/Elastic/V_Refine/refine.C',
                        'src/Elastic/V_Refine/refine_along_line.C',
                        'src/Elastic/Resid_Coarsen.C',
                        'src/Elastic/V_Coarsen/coarsen_2D.C',
                        'src/Elastic/V_Coarsen/coarsen_3D.C',
                        'src/Elastic/P_Boundary_Refine/refine.C',
                        'src/Elastic/P_Boundary_Refine/Update_P_2D.C',
                        'src/Elastic/P_Boundary_Refine/Update_P_3D.C',
                        'src/Elastic/V_Boundary_Refine/refine.C',
                        'src/Elastic/V_Boundary_Refine/Update_V_2D.C',
                        'src/Elastic/V_Boundary_Refine/Update_V_3D.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/postprocessCoarsen_2D.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/postprocessCoarsen_3D.C',
                        'src/Elastic/set_boundary.C',
                        'src/Elastic/FACOps/FACOps.C',
                        'src/Elastic/FACOps/computeCompositeResidualOnLevel.C',
                        'src/Elastic/FACOps/residual_2D.C',
                        'src/Elastic/FACOps/residual_3D.C',
                        'src/Elastic/FACOps/computeResidualNorm.C',
                        'src/Elastic/FACOps/computeVectorWeights.C',
                        'src/Elastic/FACOps/deallocateOperatorState.C',
                        'src/Elastic/FACOps/finalizeCallback.C',
                        'src/Elastic/FACOps/initializeOperatorState.C',
                        'src/Elastic/FACOps/postprocessOneCycle.C',
                        'src/Elastic/FACOps/prolongErrorAndCorrect.C',
                        'src/Elastic/FACOps/restrictResidual.C',
                        'src/Elastic/FACOps/restrictSolution.C',
                        'src/Elastic/FACOps/smoothError.C',
                        'src/Elastic/FACOps/smooth_Tackley_2D.C',
                        'src/Elastic/FACOps/smooth_Tackley_3D.C',
                        'src/Elastic/FACOps/set_boundaries.C',
                        'src/Elastic/FACOps/solveCoarsestLevel.C',
                        'src/Elastic/FACOps/solveCoarsestLevel_HYPRE.C',
                        'src/Elastic/FACOps/smooth_V_2D.C',
                        'src/Elastic/FACOps/smooth_V_3D.C',
                        'src/Elastic/FACOps/xeqScheduleGhostFill.C',
                        'src/Elastic/FACOps/xeqScheduleGhostFillNoCoarse.C',
                        'src/Elastic/FACOps/xeqScheduleProlongation.C',
                        'src/Elastic/FACOps/xeqScheduleRRestriction.C',
                        'src/Elastic/FACOps/xeqScheduleURestriction.C',
                        'src/Elastic/FACSolver/FACSolver.C',
                        'src/Elastic/FACSolver/FACSolver_Destructor.C',
                        'src/Elastic/FACSolver/createVectorWrappers.C',
                        'src/Elastic/FACSolver/deallocateSolverState.C',
                        'src/Elastic/FACSolver/destroyVectorWrappers.C',
                        'src/Elastic/FACSolver/enableLogging.C',
                        'src/Elastic/FACSolver/getFromInput.C',
                        'src/Elastic/FACSolver/initializeSolverState.C',
                        'src/Elastic/FACSolver/initializeStatics.C',
                        'src/Elastic/FACSolver/setBcObject.C',
                        'src/Elastic/FACSolver/setBoundaries.C',
                        'src/Elastic/FACSolver/solveSystem.C',
                        'src/Elastic/HypreSolver.C'],

        target       = 'gamra',
        cxxflags      = ['-O3', '-Wall', '-Wextra', '-Wconversion',
                        '-DTESTING=0', '-Drestrict='],
        # cxxflags      = ['-g', '-Wall', '-Wextra', '-Wconversion',
        #                  '-Drestrict='],
        lib          = ['SAMRAI_appu', 'SAMRAI_algs', 'SAMRAI_solv',
                        'SAMRAI_geom', 'SAMRAI_mesh', 'SAMRAI_math',
                        'SAMRAI_pdat', 'SAMRAI_xfer', 'SAMRAI_hier',
                        'SAMRAI_tbox', 'petsc', 'hdf5',
                        'HYPRE',
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
                        'gfortran', 'm'],
        libpath      = ['/home/boo/cig/SAMRAI-hg-build/lib',
                        '/usr/lib/petsc/linux-gnu-c-opt/lib'],
        includes = ['src','../FTensor','../../SAMRAI-hg-build/include',
                    '../../SAMRAI-hg/source',
                    '/usr/lib/petsc/include',
                    '/usr/lib/petsc/bmake/linux-gnu-c-opt',
                    '/usr/lib/petsc/src/vec',
                    '/usr/include']
        )
