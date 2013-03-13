/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar Elastic equation. 
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellVariable.h"
#include "Elastic/FACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

/*
*************************************************************************
*                                                                       *
* Prepare internal data for solve.                                      *
* Allocate scratch data.  Create vectors for u and f                    *
* required by the FACPreconditioner interface.                    *
* Set up internal boundary condition object.                            *
* Share data to coordinate with FAC preconditioner and                  *
* Elastic FAC operator.                                                 *
*                                                                       *
*************************************************************************
*/

void Elastic::FACSolver::initializeSolverState
(const int cell_moduli_id,
 const int edge_moduli_id,
 const int dv_diagonal_id,
 const int dv_mixed_id, 
 const int level_set_id,
 const int v_id,
 const int v_rhs_id,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
 const int coarse_level,
 const int fine_level)
{
  TBOX_ASSERT(hierarchy);
  TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

#ifdef DEBUG_CHECK_ASSERTIONS
  if (v_id < 0 || v_rhs_id < 0) {
    TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
  }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
  if (!hierarchy) {
    TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed\n"
               << "in inititialization.");
  }
#endif
  d_hierarchy = hierarchy;

  d_ln_min = coarse_level;
  d_ln_max = fine_level;
  if (d_ln_min == -1) {
    d_ln_min = 0;
  }
  if (d_ln_max == -1) {
    d_ln_max = d_hierarchy->getFinestLevelNumber();
  }

#ifdef DEBUG_CHECK_ASSERTIONS
  if (d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max) {
    TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
               << "inititialization.\n");
  }
#endif

  int ln;
  for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
    d_hierarchy->getPatchLevel(ln)->allocatePatchData(s_weight_id[d_dim.getValue() - 1]);
  }

  d_fac_ops.computeVectorWeights(d_hierarchy,
                                 s_weight_id[d_dim.getValue() - 1],
                                 d_ln_min,
                                 d_ln_max);

  d_fac_ops.set_extra_ids(cell_moduli_id,edge_moduli_id,dv_diagonal_id,
                          dv_mixed_id,level_set_id);

  createVectorWrappers(v_id, v_rhs_id);

  d_fac_precond.initializeSolverState(*d_uv, *d_fv);

  d_solver_is_initialized = true;
}
