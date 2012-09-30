/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar elastic equation. 
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellVariable.h"
#include "Elastic/FACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

#ifndef SAMRAI_INLINE
#include "Elastic/FACSolver.I"
#endif

/*
*************************************************************************
*                                                                       *
* Initialize the static data members.                                   *
*                                                                       *
*************************************************************************
*/

bool Elastic::FACSolver::s_initialized = 0;
int Elastic::FACSolver::s_weight_id[SAMRAI::tbox::Dimension::
                                    MAXIMUM_DIMENSION_VALUE];
int Elastic::FACSolver::s_instance_counter[SAMRAI::tbox::Dimension::
                                           MAXIMUM_DIMENSION_VALUE];

/*
*************************************************************************
*                                                                       *
* Constructor sets uninitialized solver state.                          *
* Set default iteration and convergence parameters.                     *
*                                                                       *
* By default settings:                                                  *
*   - Elastic equation specified has D=1, C=0.                          *
*   - State is uninitialized                                            *
*   - Logging is disabled                                               *
*   - Context for internal data is set based on object name.            *
*                                                                       *
*************************************************************************
*/

Elastic::FACSolver::FACSolver
(const SAMRAI::tbox::Dimension& dim,
 const std::string& object_name,
 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database,
 Boundary_Conditions &bc):
  d_dim(dim),
  d_object_name(object_name),
  d_boundary_conditions(bc),
  d_fac_ops(d_dim, object_name + "::fac_ops",database,bc),
  d_fac_precond(object_name + "::fac_precond", d_fac_ops),
  d_hierarchy(NULL),
  d_ln_min(-1),
  d_ln_max(-1),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::CONTEXT")),
  d_uv(NULL),
  d_fv(NULL),
  d_solver_is_initialized(false),
  d_enable_logging(false)
{

  if (!s_initialized) {
    initializeStatics();
  }

  setMaxCycles(10);
  setResidualTolerance(1e-6);
  setPresmoothingSweeps(1);
  setPostsmoothingSweeps(1);
  setCoarseFineDiscretization("Ewing");
  setCoarsestLevelSolverChoice("Tackley");
  setCoarsestLevelSolverTolerance(1e-8);
  setCoarsestLevelSolverMaxIterations(10);

  /*
   * Construct integer tag variables and add to variable database.  Note that
   * variables and patch data indices are shared among all instances.
   * The VariableDatabase holds the variables, once contructed and
   * registered via the VariableDatabase::registerInternalSAMRAIVariable()
   * function call.  Note that variables are registered and patch data indices
   * are made only for the first time through the constructor.
   */
  SAMRAI::hier::VariableDatabase*
    var_db = SAMRAI::hier::VariableDatabase::getDatabase();

  {
    static std::string cell_weight_name("Elastic::FACSolver_cell_weight");

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<double> >
      weight = var_db->getVariable(cell_weight_name);
    if (weight.isNull()) {
      weight = new SAMRAI::pdat::CellVariable<double>(d_dim, cell_weight_name, 1);
    }

    if (s_weight_id[d_dim.getValue() - 1] < 0) {
      s_weight_id[d_dim.getValue() - 1] =
        var_db->registerInternalSAMRAIVariable
        (weight,SAMRAI::hier::IntVector::getZero(d_dim));
    }
  }

  {
    static std::string side_weight_name("Elastic::FACSolver_side_weight");

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<double> >
      weight = var_db->getVariable(side_weight_name);
    if (weight.isNull()) {
      weight = new SAMRAI::pdat::SideVariable<double>(d_dim,side_weight_name,1);
    }

    if (s_weight_id[d_dim.getValue() - 2] < 0) {
      s_weight_id[d_dim.getValue() - 2] =
        var_db->registerInternalSAMRAIVariable
        (weight,SAMRAI::hier::IntVector::getZero(d_dim));
    }
  }

  d_fac_precond.setAlgorithmChoice("fas");

  /*
   * The FAC operator optionally uses the preconditioner
   * to get data for logging.
   */
  d_fac_ops.setPreconditioner(&d_fac_precond);

  if (database) {
    getFromInput(database);
  }

  s_instance_counter[d_dim.getValue() - 1]++;
}
