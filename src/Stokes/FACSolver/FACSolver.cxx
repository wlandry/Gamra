/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#include <SAMRAI/pdat/CellVariable.h>
#include "Stokes/FACSolver.hxx"
#include <SAMRAI/tbox/PIO.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/tbox/StartupShutdownManager.h>

#include IOMANIP_HEADER_FILE

/*
*************************************************************************
*                                                                       *
* Constructor sets uninitialized solver state.                          *
* Set default iteration and convergence parameters.                     *
*                                                                       *
* By default settings:                                                  *
*   - Stokes equation specified has D=1, C=0.                          *
*   - State is uninitialized                                            *
*   - Logging is disabled                                               *
*   - Context for internal data is set based on object name.            *
*                                                                       *
*************************************************************************
*/

Stokes::FACSolver::FACSolver(const SAMRAI::tbox::Dimension& dim,
                             const std::string& object_name,
                             boost::shared_ptr<SAMRAI::tbox::Database> database):
  d_dim(dim),
  d_object_name(object_name),
  d_fac_ops(boost::make_shared<FACOps>(d_dim, object_name + "::fac_ops",
                                       database)),
  d_fac_precond(object_name + "::fac_precond",d_fac_ops,database),
  d_bc_object(),
  d_simple_bc(d_dim, object_name + "::bc"),
  d_hierarchy(),
  d_level_min(-1),
  d_level_max(-1),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::CONTEXT")),
  d_uv(),
  d_fv(),
  d_solver_is_initialized(false),
  d_enable_logging(false)
{
  // #ifdef HAVE_HYPRE
  //       setCoarsestLevelSolverChoice("hypre");
  //       setCoarsestLevelSolverTolerance(1e-10);
  //       setCoarsestLevelSolverMaxIterations(20);
  //       setUseSMG(true);
  // #else
  setCoarsestLevelSolverChoice("Tackley");
  setCoarsestLevelSolverTolerance(1e-8);
  setCoarsestLevelSolverMaxIterations(10);
  // #endif

  // /*
  //  * The default RobinBcCoefStrategy used,
  //  * SimpleCellRobinBcCoefs only works with constant refine
  //  * for prolongation.  So we use constant refinement
  //  * for prolongation by default.
  //  */
  // setProlongationMethod("CONSTANT_REFINE");

  /*
   * The FAC operator optionally uses the preconditioner
   * to get data for logging.
   */
  d_fac_ops->setPreconditioner((const SAMRAI::solv::FACPreconditioner *)
                               (&d_fac_precond));

  if (database) {
    getFromInput(database);
  }
}
