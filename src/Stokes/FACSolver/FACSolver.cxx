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
* Initialize the static data members.                                   *
*                                                                       *
*************************************************************************
*/

bool Stokes::FACSolver::s_initialized = 0;
int Stokes::FACSolver::s_weight_id[SAMRAI::MAX_DIM_VAL];
int Stokes::FACSolver::s_instance_counter[SAMRAI::MAX_DIM_VAL];

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
  d_ln_min(-1),
  d_ln_max(-1),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::CONTEXT")),
  d_uv(),
  d_fv(),
  d_solver_is_initialized(false),
  d_enable_logging(false)
{

  if (!s_initialized) {
    initializeStatics();
  }

  setCoarseFineDiscretization("Ewing");
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

  /*
   * Construct integer tag variables and add to variable database.  Note that
   * variables and patch data indices are shared among all instances.
   * The VariableDatabase holds the variables, once contructed and
   * registered via the VariableDatabase::registerInternalSAMRAIVariable()
   * function call.  Note that variables are registered and patch data indices
   * are made only for the first time through the constructor.
   */
  SAMRAI::hier::VariableDatabase* var_db = SAMRAI::hier::VariableDatabase::getDatabase();

  {
    static std::string cell_weight_name("Stokes::FACSolver_cell_weight");

    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > weight =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
      (var_db->getVariable(cell_weight_name));
    if (!weight) {
      weight = boost::make_shared<SAMRAI::pdat::CellVariable<double> >
        (d_dim, cell_weight_name, 1);
    }

    if (s_weight_id[d_dim.getValue() - 1] < 0) {
      s_weight_id[d_dim.getValue() - 1] =
        var_db->registerInternalSAMRAIVariable
        (weight,SAMRAI::hier::IntVector::getZero(d_dim));
    }
  }

  {
    static std::string side_weight_name("Stokes::FACSolver_side_weight");

    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > weight =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
      (var_db->getVariable(side_weight_name));
    if (!weight) {
      weight = boost::make_shared<SAMRAI::pdat::SideVariable<double> >
        (d_dim, side_weight_name,SAMRAI::hier::IntVector::getOne(d_dim),1);
    }

    if (s_weight_id[d_dim.getValue() - 2] < 0) {
      s_weight_id[d_dim.getValue() - 2] =
        var_db->registerInternalSAMRAIVariable
        (weight,SAMRAI::hier::IntVector::getZero(d_dim));
    }
  }

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

  s_instance_counter[d_dim.getValue() - 1]++;
}
