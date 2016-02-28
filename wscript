# top = '.'
# out = 'build'

def options(opt):
    opt.load('compiler_cxx FTensor okada boost')

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
    conf.setenv('prof')
    configure_variant(conf);

def configure_variant(conf):
    conf.load('compiler_cxx FTensor okada boost')
    # conf.env.SHLIB_MARKER = '-Wl,-Bstatic'
    conf.check_boost()

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
                        'SAMRAI_tbox'],
                   use=['BOOST'])

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

    # Optimization flags
    optimize_msg="Checking for optimization flag "
    optimize_fragment="int main() {}\n"
    optimize_flags=['-Ofast','-O3','-O2','-O']
    for flag in optimize_flags:
        try:
            conf.check_cxx(msg=optimize_msg+flag, fragment=optimize_fragment,
                           cxxflags=flag, uselib_store='optimize')
        except conf.errors.ConfigurationError:
            continue
        else:
            found_optimize=True
            break

def build(bld):
    default_flags=['-Wall', '-Wextra', '-Wconversion', '-Wvla']
    cxxflags_variant= {'release' : ['-Ofast', '-DTESTING=0'],
                    'prof' : ['-pg','-Ofast', '-DTESTING=0'],
                    'debug' : ['-g']}
    linkflags_variant={'release' : [],
                       'prof' : ['-pg'],
                       'debug' : []}

    use_array=['samrai','muparser','hdf5','FTensor','okada','BOOST']
    if bld.variant=='release':
        use_array.append('optimize')

    bld.program(
        features     = ['cxx','cprogram'],
        source       = ['src/main.cxx',
                        'src/Elastic/FAC/FAC.cxx',
                        'src/Elastic/FAC/fix_moduli.cxx',
                        'src/Elastic/FAC/intersection.cxx',
                        'src/Elastic/FAC/applyGradientDetector.cxx',
                        'src/Elastic/FAC/initializeLevelData.cxx',
                        'src/Elastic/FAC/pack_strain.cxx',
                        'src/Elastic/FAC/pack_level_set.cxx',
                        'src/Elastic/FAC/pack_v_v_rhs.cxx',
                        'src/Elastic/FAC/pack_v_initial.cxx',
                        'src/Elastic/FAC/resetHierarchyConfiguration.cxx',
                        'src/Elastic/FAC/setupPlotter.cxx',
                        'src/Elastic/FAC/solve.cxx',
                        'src/Elastic/V_Refine/V_Refine.cxx',
                        'src/Elastic/V_Refine/refine.cxx',
                        'src/Elastic/V_Refine/refine_box.cxx',
                        'src/Elastic/V_Refine/refine_along_line.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/Coarse_Fine_Boundary_Refine.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/refine.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/refine_box.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/Update_V_2D.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/Update_V_embedded_2D.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/Correction_2D.cxx',
                        'src/Elastic/Coarse_Fine_Boundary_Refine/Update_V_3D.cxx',
                        'src/Elastic/V_Coarsen_Patch_Strategy/postprocessCoarsen.cxx',
                        'src/Elastic/V_Coarsen_Patch_Strategy/coarsen_2D.cxx',
                        'src/Elastic/V_Coarsen_Patch_Strategy/coarsen_3D.cxx',
                        'src/Elastic/V_Coarsen_Patch_Strategy/fix_boundary_elements_2D.cxx',
                        'src/Elastic/V_Coarsen_Patch_Strategy/fix_boundary_elements_3D.cxx',
                        'src/Elastic/Boundary_Conditions/Boundary_Conditions.cxx',
                        'src/Elastic/Boundary_Conditions/set_physical_boundary.cxx',
                        'src/Elastic/Boundary_Conditions/set_embedded_boundary.cxx',
                        'src/Elastic/Boundary_Conditions/set_regular_boundary/set_regular_boundary.cxx',
                        'src/Elastic/Boundary_Conditions/set_regular_boundary/set_dirichlet.cxx',
                        'src/Elastic/Boundary_Conditions/set_regular_boundary/set_shear_stress.cxx',
                        'src/Elastic/FACOps/FACOps.cxx',
                        'src/Elastic/FACOps/v_level_set_operator_2D.cxx',
                        'src/Elastic/FACOps/computeCompositeResidualOnLevel.cxx',
                        'src/Elastic/FACOps/residual_2D.cxx',
                        'src/Elastic/FACOps/residual_3D.cxx',
                        'src/Elastic/FACOps/computeResidualNorm.cxx',
                        'src/Elastic/FACOps/deallocateOperatorState.cxx',
                        'src/Elastic/FACOps/finalizeCallback.cxx',
                        'src/Elastic/FACOps/initializeOperatorState.cxx',
                        'src/Elastic/FACOps/postprocessOneCycle.cxx',
                        'src/Elastic/FACOps/prolongErrorAndCorrect.cxx',
                        'src/Elastic/FACOps/restrictResidual.cxx',
                        'src/Elastic/FACOps/restrictSolution.cxx',
                        'src/Elastic/FACOps/smoothError.cxx',
                        'src/Elastic/FACOps/Gauss_Seidel_red_black_2D.cxx',
                        'src/Elastic/FACOps/Gauss_Seidel_red_black_3D.cxx',
                        'src/Elastic/FACOps/set_physical_boundaries.cxx',
                        'src/Elastic/FACOps/solveCoarsestLevel.cxx',
                        'src/Elastic/FACOps/update_V_2D.cxx',
                        'src/Elastic/FACOps/update_V_3D.cxx',
                        'src/Elastic/FACOps/ghostfill.cxx',
                        'src/Elastic/FACOps/ghostfill_nocoarse.cxx',
                        'src/Elastic/FACOps/refine.cxx',
                        'src/Elastic/FACOps/coarsen_resid.cxx',
                        'src/Elastic/FACOps/coarsen_solution.cxx',
                        'src/Elastic/FACSolver/FACSolver.cxx',
                        'src/Elastic/FACSolver/FACSolver_Destructor.cxx',
                        'src/Elastic/FACSolver/createVectorWrappers.cxx',
                        'src/Elastic/FACSolver/deallocateSolverState.cxx',
                        'src/Elastic/FACSolver/getFromInput.cxx',
                        'src/Elastic/FACSolver/initializeSolverState.cxx',
                        'src/Stokes/FAC/FAC.cxx',
                        'src/Stokes/FAC/fix_viscosity.cxx',
                        'src/Stokes/FAC/applyGradientDetector.cxx',
                        'src/Stokes/FAC/initializeLevelData.cxx',
                        'src/Stokes/FAC/packDerivedDataIntoDoubleBuffer.cxx',
                        'src/Stokes/FAC/resetHierarchyConfiguration.cxx',
                        'src/Stokes/FAC/setupPlotter.cxx',
                        'src/Stokes/FAC/solve.cxx',
                        'src/Stokes/P_Refine.cxx',
                        'src/Stokes/V_Refine.cxx',
                        'src/Stokes/Resid_Coarsen.cxx',
                        'src/Stokes/V_Coarsen/coarsen_2D.cxx',
                        'src/Stokes/V_Coarsen/coarsen_3D.cxx',
                        'src/Stokes/P_Boundary_Refine/refine.cxx',
                        'src/Stokes/P_Boundary_Refine/Update_P_2D.cxx',
                        'src/Stokes/P_Boundary_Refine/Update_P_3D.cxx',
                        'src/Stokes/V_Boundary_Refine/refine.cxx',
                        'src/Stokes/V_Boundary_Refine/Update_V_2D.cxx',
                        'src/Stokes/V_Boundary_Refine/Update_V_3D.cxx',
                        'src/Stokes/V_Coarsen_Patch_Strategy/postprocessCoarsen_2D.cxx',
                        'src/Stokes/V_Coarsen_Patch_Strategy/postprocessCoarsen_3D.cxx',
                        'src/Stokes/set_boundary.cxx',
                        'src/Stokes/FACOps/FACOps.cxx',
                        'src/Stokes/FACOps/computeCompositeResidualOnLevel.cxx',
                        'src/Stokes/FACOps/residual_2D.cxx',
                        'src/Stokes/FACOps/residual_3D.cxx',
                        'src/Stokes/FACOps/computeResidualNorm.cxx',
                        'src/Stokes/FACOps/deallocateOperatorState.cxx',
                        'src/Stokes/FACOps/finalizeCallback.cxx',
                        'src/Stokes/FACOps/initializeOperatorState.cxx',
                        'src/Stokes/FACOps/postprocessOneCycle.cxx',
                        'src/Stokes/FACOps/prolongErrorAndCorrect.cxx',
                        'src/Stokes/FACOps/restrictResidual.cxx',
                        'src/Stokes/FACOps/restrictSolution.cxx',
                        'src/Stokes/FACOps/smoothError.cxx',
                        'src/Stokes/FACOps/smooth_Tackley_2D.cxx',
                        'src/Stokes/FACOps/smooth_Tackley_3D.cxx',
                        'src/Stokes/FACOps/smooth_Gerya.cxx',
                        'src/Stokes/FACOps/set_physical_boundaries.cxx',
                        'src/Stokes/FACOps/solveCoarsestLevel.cxx',
                        'src/Stokes/FACOps/solveCoarsestLevel_HYPRE.cxx',
                        'src/Stokes/FACOps/smooth_V_2D.cxx',
                        'src/Stokes/FACOps/smooth_V_3D.cxx',
                        'src/Stokes/FACOps/xeqScheduleGhostFill.cxx',
                        'src/Stokes/FACOps/xeqScheduleGhostFillNoCoarse.cxx',
                        'src/Stokes/FACOps/xeqScheduleProlongation.cxx',
                        'src/Stokes/FACOps/xeqScheduleRRestriction.cxx',
                        'src/Stokes/FACOps/xeqScheduleURestriction.cxx',
                        'src/Stokes/FACSolver/FACSolver.cxx',
                        'src/Stokes/FACSolver/FACSolver_Destructor.cxx',
                        'src/Stokes/FACSolver/createVectorWrappers.cxx',
                        'src/Stokes/FACSolver/deallocateSolverState.cxx',
                        'src/Stokes/FACSolver/destroyVectorWrappers.cxx',
                        'src/Stokes/FACSolver/enableLogging.cxx',
                        'src/Stokes/FACSolver/getFromInput.cxx',
                        'src/Stokes/FACSolver/initializeSolverState.cxx',
                        'src/Stokes/FACSolver/setBcObject.cxx',
                        'src/Stokes/FACSolver/setBoundaries.cxx',
                        'src/Stokes/FACSolver/solveSystem.cxx',
                        'src/Stokes/HypreSolver.cxx'],
        target       = 'gamra',
        cxxflags     = cxxflags_variant[bld.variant] + default_flags,
        lib          = ['dl','gfortranbegin', 'gfortran', 'm'],
        libpath      = ['/sw/lib/gcc4.8/lib','/sw/lib/gcc4.8/lib/gcc/x86_64-apple-darwin13.0.0/4.8.2'],
        linkflags    = linkflags_variant[bld.variant],
        includes = ['src'],
        use=use_array
        )

from waflib.Build import BuildContext, CleanContext, \
    InstallContext, UninstallContext

for x in 'debug prof release'.split():
    for y in (BuildContext, CleanContext, InstallContext, UninstallContext):
        name = y.__name__.replace('Context','').lower()
        class tmp(y):
            if x=='release':
                cmd = name
            else:
                cmd = name + '_' + x
            variant = x
