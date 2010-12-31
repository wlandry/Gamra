/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "FACStokes.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "StokesSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

extern "C" {
void F77_FUNC(setexactandrhs2d, SETEXACTANDRHS2D) (const int & ifirst0,
   const int & ilast0,
   const int & ifirst1,
   const int & ilast1,
   double* exact,
   double* rhs,
   const double* dx,
   const double* xlower);

void F77_FUNC(setexactandrhs3d, SETEXACTANDRHS3D) (const int & ifirst0,
   const int & ilast0,
   const int & ifirst1,
   const int & ilast1,
   const int & ifirst2,
   const int & ilast2,
   double* exact,
   double* rhs,
   const double* dx,
   const double* xlower);
}

namespace SAMRAI {

/*
 *************************************************************************
 * Constructor creates a unique context for the object and register      *
 * all its internal variables with the variable database.                *
 *************************************************************************
 */
FACStokes::FACStokes(
   const std::string& object_name,
   const tbox::Dimension& dim,
   tbox::Pointer<tbox::Database> database):
   d_object_name(object_name),
   d_dim(dim),
   d_hierarchy(NULL),
   d_stokes_fac_solver((d_dim),
                        object_name + "::stokes_hypre",
                        (!database.isNull() &&
                         database->isDatabase("fac_solver")) ?
                        database->getDatabase("fac_solver"):
                           tbox::Pointer<tbox::Database>(NULL)),
   d_bc_coefs(d_dim,
              object_name + "::bc_coefs",
              (!database.isNull() &&
               database->isDatabase("bc_coefs")) ?
              database->getDatabase("bc_coefs"):
                 tbox::Pointer<tbox::Database>(NULL)),
   d_context()
{

   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();

   /*
    * Get a unique context for variables owned by this object.
    */
   d_context = vdb->getContext(d_object_name + ":Context");

   /*
    * Register variables with hier::VariableDatabase
    * and get the descriptor indices for those variables.
    */

   tbox::Pointer<pdat::CellVariable<double> >
     p(new pdat::CellVariable<double>(dim, object_name + ":p", 1));
   p_id = vdb->registerVariableAndContext(p, d_context, hier::IntVector(dim, 1)
                                          /* ghost cell width is 1 for
                                             stencil widths */);

   tbox::Pointer<pdat::CellVariable<double> >
     p_exact(new pdat::CellVariable<double>(dim, object_name + ":p exact"));
   p_exact_id = vdb->registerVariableAndContext(p_exact,d_context,
                                                hier::IntVector(dim, 1)
                                                /* ghost cell width is
                                                   1 in case needed */);

   tbox::Pointer<pdat::CellVariable<double> >
     p_rhs(new pdat::CellVariable<double>(dim,object_name
                                          + ":p right hand side"));
                                          
   p_rhs_id = vdb->registerVariableAndContext(p_rhs,d_context,
                                              hier::IntVector(dim, 0)
                                              /* ghost cell width is 0 */);

   tbox::Pointer<pdat::FaceVariable<double> >
     v(new pdat::FaceVariable<double>(dim, object_name + ":v", 1));
   v_id = vdb->registerVariableAndContext(v, d_context, hier::IntVector(dim, 1)
                                          /* ghost cell width is 1 for
                                             stencil widths */);

   tbox::Pointer<pdat::FaceVariable<double> >
     v_rhs(new pdat::FaceVariable<double>(dim,object_name
                                          + ":v right hand side"));
   v_rhs_id = vdb->registerVariableAndContext(v_rhs,d_context,
                                              hier::IntVector(dim, 0)
                                              /* ghost cell width is 0 */);

   /*
    * Specify an implementation of solv::RobinBcCoefStrategy for the
    * solver to use.  We use the implementation
    * solv::LocationIndexRobinBcCoefs, but other implementations are
    * possible, including user-implemented.
    */
   d_stokes_fac_solver.setBcObject(&d_bc_coefs);
}

/*
 *************************************************************************
 * Destructor does nothing interesting                                   *
 *************************************************************************
 */
FACStokes::~FACStokes()
{
}

/*
 *************************************************************************
 * Initialize data on a level.                                           *
 *                                                                       *
 * Allocate the solution, exact solution and rhs memory.                 *
 * Fill the rhs and exact solution.                                      *
 *************************************************************************
 */
void FACStokes::initializeLevelData(
   const tbox::Pointer<hier::BasePatchHierarchy> patch_hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer<hier::BasePatchLevel> old_level,
   const bool allocate_data)
{

   (void)init_data_time;
   (void)can_be_refined;
   (void)initial_time;
   (void)old_level;

   tbox::Pointer<hier::PatchHierarchy> hierarchy = patch_hierarchy;
   tbox::Pointer<geom::CartesianGridGeometry> grid_geom =
      hierarchy->getGridGeometry();

   tbox::Pointer<hier::PatchLevel> level =
      hierarchy->getPatchLevel(level_number);

   if (allocate_data) {
      level->allocatePatchData(p_id);
      level->allocatePatchData(p_rhs_id);
      level->allocatePatchData(p_exact_id);
      level->allocatePatchData(v_id);
      level->allocatePatchData(v_rhs_id);
   }

   /*
    * Initialize data in all patches in the level.
    */
   hier::PatchLevel::Iterator pi(*level);
   for (pi.initialize(*level); pi; pi++) {

      tbox::Pointer<hier::Patch> patch = *pi;
      if (patch.isNull()) {
         TBOX_ERROR(d_object_name
            << ": Cannot find patch.  Null patch pointer.");
      }
      hier::Box pbox = patch->getBox();
      tbox::Pointer<geom::CartesianPatchGeometry> patch_geom =
         patch->getPatchGeometry();

      tbox::Pointer<pdat::CellData<double> > exact_data =
         patch->getPatchData(p_exact_id);
      tbox::Pointer<pdat::CellData<double> > rhs_data =
         patch->getPatchData(p_rhs_id);

      /*
       * Set source function and exact solution.
       */
      if (d_dim == tbox::Dimension(2)) {
         F77_FUNC(setexactandrhs2d, SETEXACTANDRHS2D) (
            pbox.lower()[0],
            pbox.upper()[0],
            pbox.lower()[1],
            pbox.upper()[1],
            exact_data->getPointer(),
            rhs_data->getPointer(),
            patch_geom->getDx(),
            patch_geom->getXLower());
      } else if (d_dim == tbox::Dimension(3)) {
         F77_FUNC(setexactandrhs3d, SETEXACTANDRHS3D) (
            pbox.lower()[0],
            pbox.upper()[0],
            pbox.lower()[1],
            pbox.upper()[1],
            pbox.lower()[2],
            pbox.upper()[2],
            exact_data->getPointer(),
            rhs_data->getPointer(),
            patch_geom->getDx(),
            patch_geom->getXLower());
      }

   }    // End patch loop.
}

/*
 *************************************************************************
 * Reset the hierarchy-dependent internal information.                   *
 *************************************************************************
 */
void FACStokes::resetHierarchyConfiguration(
   tbox::Pointer<hier::BasePatchHierarchy> new_hierarchy,
   int coarsest_level,
   int finest_level)
{
   (void)coarsest_level;
   (void)finest_level;

   d_hierarchy = new_hierarchy;
}

/*
 *************************************************************************
 * Set up the initial guess and problem parameters                       *
 * and solve the Stokes problem.  We explicitly initialize and          *
 * deallocate the solver state in this example.                          *
 *************************************************************************
 */
int FACStokes::solveStokes()
{

   if (d_hierarchy.isNull()) {
      TBOX_ERROR(d_object_name
         << "Cannot solve using an uninitialized object.\n");
   }

   int ln;
   /*
    * Fill in the initial guess.
    */
   for (ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::Iterator ip(*level);
      for ( ; ip; ip++) {
         tbox::Pointer<hier::Patch> patch = *ip;
         tbox::Pointer<pdat::CellData<double> > data = patch->getPatchData(
               p_id);
         data->fill(0.0);
         tbox::Pointer<pdat::FaceData<double> > vdata = patch->getPatchData(
               v_id);
         vdata->fill(0.0);
      }
   }

   /*
    * Set the parameters for the Stokes equation.
    * See classes solv::CellStokesFACSolver or
    * solv::StokesSpecifications.
    * (D is the diffusion coefficient.
    * C is the source term which is not used in this example.)
    */
   d_stokes_fac_solver.setDConstant(1.0);
   d_stokes_fac_solver.setCConstant(0.0);

   d_stokes_fac_solver.initializeSolverState(
      p_id,
      p_rhs_id,
      d_hierarchy,
      0,
      d_hierarchy->getFinestLevelNumber());

   tbox::plog << "solving..." << std::endl;
   int solver_ret;
   solver_ret = d_stokes_fac_solver.solveSystem(p_id,
         p_rhs_id);
   /*
    * Present data on the solve.
    */
   double avg_factor, final_factor;
   d_stokes_fac_solver.getConvergenceFactors(avg_factor, final_factor);
   tbox::plog << "\t" << (solver_ret ? "" : "NOT ") << "converged " << "\n"
              << "	iterations: "
              << d_stokes_fac_solver.getNumberOfIterations() << "\n"
              << "	residual: "<< d_stokes_fac_solver.getResidualNorm()
              << "\n"
              << "	average convergence: "<< avg_factor << "\n"
              << "	final convergence: "<< final_factor << "\n"
              << std::flush;

   d_stokes_fac_solver.deallocateSolverState();

   return 0;
}

#ifdef HAVE_HDF5
/*
 *************************************************************************
 * Set up external plotter to plot internal data from this class.        *
 * Register variables appropriate for plotting.                          *
 *************************************************************************
 */
int FACStokes::setupPlotter(
   appu::VisItDataWriter& plotter) const {
   if (d_hierarchy.isNull()) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                               << " FACStokes::setupPlotter\n"
                               << "The hierarchy must be set before calling\n"
                               << "this function.\n");
   }
   plotter.registerPlotQuantity("Computed solution",
      "SCALAR",
      p_id);
   plotter.registerDerivedPlotQuantity("Error",
      "SCALAR",
      (appu::VisDerivedDataStrategy *)this);
   plotter.registerPlotQuantity("Exact solution",
      "SCALAR",
      p_exact_id);
   plotter.registerPlotQuantity("Stokes source",
      "SCALAR",
      p_rhs_id);

   return 0;
}
#endif

/*
 *************************************************************************
 * Write derived data to the given stream.                               *
 *************************************************************************
 */
bool FACStokes::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id) const
{
   (void)depth_id;

   pdat::CellData<double>::Iterator icell(region);

   if (variable_name == "Error") {
      tbox::Pointer<pdat::CellData<double> > current_solution_ =
         patch.getPatchData(p_id);
      tbox::Pointer<pdat::CellData<double> > exact_solution_ =
         patch.getPatchData(p_exact_id);
      pdat::CellData<double>& current_solution = *current_solution_;
      pdat::CellData<double>& exact_solution = *exact_solution_;
      for ( ; icell; icell++) {
         double diff = (current_solution(*icell) - exact_solution(*icell));
         *buffer = diff;
         buffer = buffer + 1;
      }
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
         <<
         "FACStokesX::packDerivedDataIntoDoubleBuffer");

   }
   // Return true if this patch has derived data on it.
   // False otherwise.
   return true;
}

}
