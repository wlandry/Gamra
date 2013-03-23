# top = '.'
# out = 'build'

def options(opt):
    opt.load('compiler_cxx')

    samrai=opt.add_option_group('SAMRAI Options')
    samrai.add_option('--samrai-dir',
                   help='Base directory where samrai is installed')
    samrai.add_option('--samrai-incdir',
                   help='Directory where samrai include files are installed')
    samrai.add_option('--samrai-libdir',
                   help='Directory where samrai library files are installed')
    samrai.add_option('--samrai-libs',
                   help='Names of the samrai libraries without prefix or suffix\n'
                   '(e.g. "SAMRAI_appu SAMRAI_algs SAMRAI_solv"')

    muparser=opt.add_option_group('muParser Options')
    muparser.add_option('--muparser-dir',
                   help='Base directory where muparser is installed')
    muparser.add_option('--muparser-incdir',
                   help='Directory where muparser include files are installed')
    muparser.add_option('--muparser-libdir',
                   help='Directory where muparser library files are installed')
    muparser.add_option('--muparser-libs',
                   help='Names of the muparser libraries without prefix or suffix\n'
                   '(e.g. "muparser"')

    hdf5=opt.add_option_group('HDF5 Options')
    hdf5.add_option('--hdf5-dir',
                   help='Base directory where HDF5 is installed')
    hdf5.add_option('--hdf5-incdir',
                   help='Directory where HDF5 include files are installed')
    hdf5.add_option('--hdf5-libdir',
                   help='Directory where HDF5 library files are installed')
    hdf5.add_option('--hdf5-libs',
                   help='Names of the HDF5 libraries without prefix or suffix\n'
                   '(e.g. "hdf5"')

def configure(conf):
    conf.setenv('debug')
    configure_variant(conf);
    conf.setenv('release')
    configure_variant(conf);

def configure_variant(conf):
    conf.load('compiler_cxx')

    # Find SAMRAI
    if conf.options.samrai_dir:
        if not conf.options.samrai_incdir:
            conf.options.samrai_incdir=conf.options.samrai_dir + "/include"
        if not conf.options.samrai_libdir:
            conf.options.samrai_libdir=conf.options.samrai_dir + "/lib"
    frag="#include \"SAMRAI/SAMRAI_config.h\"\n" + 'int main()"\n' \
        + "{}\n"
    if not conf.options.samrai_incdir:
        conf.options.samrai_incdir='/usr/include'

    conf.check_cxx(msg="Checking for SAMRAI",
                  header_name='SAMRAI/SAMRAI_config.h',
                  includes=[conf.options.samrai_incdir], uselib_store='samrai',
                  libpath=[conf.options.samrai_libdir],
                  rpath=[conf.options.samrai_libdir],
                  lib=['SAMRAI_appu', 'SAMRAI_algs', 'SAMRAI_solv',
                       'SAMRAI_geom', 'SAMRAI_mesh', 'SAMRAI_math',
                       'SAMRAI_pdat', 'SAMRAI_xfer', 'SAMRAI_hier',
                       'SAMRAI_tbox'])

    # Find muParser
    if conf.options.muparser_dir:
        if not conf.options.muparser_incdir:
            conf.options.muparser_incdir=conf.options.muparser_dir + "/include"
        if not conf.options.muparser_libdir:
            conf.options.muparser_libdir=conf.options.muparser_dir + "/lib"
    frag="#include \"muParser.h\"\n" + 'int main()"\n' \
        + "{}\n"

    conf.check_cxx(msg="Checking for muParser",
                  header_name='muParser.h',
                  includes=[conf.options.muparser_incdir], uselib_store='muparser',
                  libpath=[conf.options.muparser_libdir],
                  rpath=[conf.options.muparser_libdir],
                  lib=['muparser'])

    # Find Parallel HDF5
    if conf.options.hdf5_dir:
        if not conf.options.hdf5_incdir:
            conf.options.hdf5_incdir=conf.options.hdf5_dir + "/include"
        if not conf.options.hdf5_libdir:
            conf.options.hdf5_libdir=conf.options.hdf5_dir + "/lib"
    frag="#include \"mpi.h\"\n#include \"hdf5.h\"\n" + 'int main()"\n' \
        + "{}\n"

    conf.check_cxx(msg="Checking for hdf5",
                  header_name='hdf5.h',
                  includes=[conf.options.hdf5_incdir], uselib_store='hdf5',
                  libpath=[conf.options.hdf5_libdir],
                  rpath=[conf.options.hdf5_libdir],
                  lib=['hdf5'])

def build(bld):
    if not bld.variant:
        bld.fatal('call "waf build_debug" or "waf build_release", and try "waf --help"')
    default_flags=['-Wall', '-Wextra', '-Wconversion', '-Drestrict=']
    variant_flags= {'release' : ['-O3', '-DTESTING=0'],
                    'debug' : ['-g']}
    bld.program(
        features     = ['cxx','cprogram'],
        source       = ['src/main.C',
                        'src/Elastic/FAC/FAC.C',
                        'src/Elastic/FAC/fix_moduli.C',
                        'src/Elastic/FAC/intersection.C',
                        'src/Elastic/FAC/applyGradientDetector.C',
                        'src/Elastic/FAC/initializeLevelData.C',
                        'src/Elastic/FAC/pack_strain.C',
                        'src/Elastic/FAC/pack_level_set.C',
                        'src/Elastic/FAC/pack_v_v_rhs.C',
                        'src/Elastic/FAC/resetHierarchyConfiguration.C',
                        'src/Elastic/FAC/setupPlotter.C',
                        'src/Elastic/FAC/solve.C',
                        'src/Elastic/V_Refine/refine.C',
                        'src/Elastic/V_Refine/refine_along_line.C',
                        'src/Elastic/V_Boundary_Refine/V_Boundary_Refine.C',
                        'src/Elastic/V_Boundary_Refine/refine.C',
                        'src/Elastic/V_Boundary_Refine/Update_V_2D.C',
                        'src/Elastic/V_Boundary_Refine/Correction_2D.C',
                        'src/Elastic/V_Boundary_Refine/Update_V_3D.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/postprocessCoarsen.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/coarsen_2D.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/coarsen_3D.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/fix_boundary_elements_2D.C',
                        'src/Elastic/V_Coarsen_Patch_Strategy/fix_boundary_elements_3D.C',
                        'src/Elastic/Boundary_Conditions/Boundary_Conditions.C',
                        'src/Elastic/Boundary_Conditions/set_boundary.C',
                        'src/Elastic/Boundary_Conditions/set_embedded_boundary.C',
                        'src/Elastic/Boundary_Conditions/set_regular_boundary/set_regular_boundary.C',
                        'src/Elastic/Boundary_Conditions/set_regular_boundary/set_dirichlet.C',
                        'src/Elastic/Boundary_Conditions/set_regular_boundary/set_shear_derivs.C',
                        'src/Elastic/FACOps/FACOps.C',
                        'src/Elastic/FACOps/v_level_set_operator_2D.C',
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
                        'src/Elastic/FACOps/smooth_2D.C',
                        'src/Elastic/FACOps/smooth_3D.C',
                        'src/Elastic/FACOps/set_boundaries.C',
                        'src/Elastic/FACOps/solveCoarsestLevel.C',
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
                        'src/Elastic/FACSolver/solveSystem.C',
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
                        'src/Stokes/HypreSolver.C'],
        target       = 'gamra',
        cxxflags = variant_flags[bld.variant] + default_flags,
        lib          = ['dl','gfortranbegin', 'gfortran', 'm'],
        libpath      = ['/sw/lib/gcc4.7/lib','/sw/lib/gcc4.7/lib/gcc/x86_64-apple-darwin11.4.0/4.7.2'],
        includes = ['src','../FTensor'],
        use=['samrai','muparser','hdf5']
        )

from waflib.Build import BuildContext, CleanContext, \
    InstallContext, UninstallContext

for x in 'debug release'.split():
    for y in (BuildContext, CleanContext, InstallContext, UninstallContext):
        name = y.__name__.replace('Context','').lower()
        class tmp(y):
            cmd = name + '_' + x
            variant = x
