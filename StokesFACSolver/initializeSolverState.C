/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellVariable.h"
#include "StokesFACSolver.h"
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
* Stokes FAC operator.                                                 *
*                                                                       *
*************************************************************************
*/

void SAMRAI::solv::StokesFACSolver::initializeSolverState
(const int p,
 const int p_rhs,
 const int v,
 const int v_rhs,
 tbox::Pointer<hier::PatchHierarchy> hierarchy,
 const int coarse_level,
 const int fine_level)
{
  TBOX_ASSERT(!hierarchy.isNull());
  TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

  if (d_bc_object == NULL) {
    TBOX_ERROR(
               d_object_name << ": No BC coefficient strategy object!\n"
               <<
               "Use either setBoundaries or setPhysicalBcCoefObject\n"
               << "to specify the boundary conidition.\n");
  }

#ifdef DEBUG_CHECK_ASSERTIONS
  if (p < 0 || p_rhs < 0) {
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

  if (d_bc_object == &d_simple_bc) {
    d_simple_bc.setHierarchy(d_hierarchy,
                             d_ln_min,
                             d_ln_max);
    if (d_stokes_spec.dIsConstant()) {
      d_simple_bc.setDiffusionCoefConstant(d_stokes_spec.getDConstant());
    } else {
      d_simple_bc.setDiffusionCoefId(d_stokes_spec.getDPatchDataId());
    }
  }

  d_fac_ops.setStokesSpecifications(d_stokes_spec);

  createVectorWrappers(p, p_rhs, v, v_rhs);

  d_fac_precond.initializeSolverState(*d_uv, *d_fv);

  d_solver_is_initialized = true;
}
