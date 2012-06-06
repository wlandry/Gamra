/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#ifndef included_solv_StokesFACSolver_C
#define included_solv_StokesFACSolver_C

#include "SAMRAI/pdat/CellVariable.h"
#include "Stokes/FACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

#ifndef SAMRAI_INLINE
#include "Stokes/FACSolver.I"
#endif

namespace SAMRAI {
  namespace solv {

    /*
*************************************************************************
*                                                                       *
* Initialize the static data members.                                   *
*                                                                       *
*************************************************************************
*/

    bool Stokes::FACSolver::s_initialized = 0;
    int Stokes::FACSolver::s_weight_id[SAMRAI::tbox::Dimension::
                                     MAXIMUM_DIMENSION_VALUE];
    int Stokes::FACSolver::s_instance_counter[SAMRAI::tbox::Dimension::
                                            MAXIMUM_DIMENSION_VALUE];

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

    Stokes::FACSolver::FACSolver(const tbox::Dimension& dim,
                                 const std::string& object_name,
                                 tbox::Pointer<tbox::Database> database):
      d_dim(dim),
      d_object_name(object_name),
      d_fac_ops(d_dim, object_name + "::fac_ops",database),
      d_fac_precond(object_name + "::fac_precond", d_fac_ops),
      d_bc_object(NULL),
      d_simple_bc(d_dim, object_name + "::bc"),
      d_hierarchy(NULL),
      d_ln_min(-1),
      d_ln_max(-1),
      d_context(hier::VariableDatabase::getDatabase()
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
      hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

      {
        static std::string cell_weight_name("Stokes::FACSolver_cell_weight");

        tbox::Pointer<pdat::CellVariable<double> >
          weight = var_db->getVariable(cell_weight_name);
        if (weight.isNull()) {
          weight = new pdat::CellVariable<double>(d_dim, cell_weight_name, 1);
        }

        if (s_weight_id[d_dim.getValue() - 1] < 0) {
          s_weight_id[d_dim.getValue() - 1] =
            var_db->registerInternalSAMRAIVariable
            (weight,hier::IntVector::getZero(d_dim));
        }
      }

      {
        static std::string side_weight_name("Stokes::FACSolver_side_weight");

        tbox::Pointer<pdat::SideVariable<double> >
          weight = var_db->getVariable(side_weight_name);
        if (weight.isNull()) {
          weight = new pdat::SideVariable<double>(d_dim, side_weight_name, 1);
        }

        if (s_weight_id[d_dim.getValue() - 2] < 0) {
          s_weight_id[d_dim.getValue() - 2] =
            var_db->registerInternalSAMRAIVariable
            (weight,hier::IntVector::getZero(d_dim));
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
      d_fac_ops.setPreconditioner((const FACPreconditioner *)(&d_fac_precond));

      if (database) {
        getFromInput(database);
      }

      s_instance_counter[d_dim.getValue() - 1]++;
    }

  }
}
#endif
