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
* Prepare internal data for solve.                                      *
* Allocate scratch data.  Create vectors for u and f                    *
* required by the FACPreconditioner interface.                    *
* Set up internal boundary condition object.                            *
* Share data to coordinate with FAC preconditioner and                  *
* Stokes FAC operator.                                                 *
*                                                                       *
*************************************************************************
*/

void Stokes::FACSolver::initializeSolverState
(const int p,
 const int cell_viscosity,
 const int edge_viscosity,
 const int dp,
 const int p_rhs,
 const int v,
 const int v_rhs,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
 const int coarse_level,
 const int fine_level)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

  if (!d_bc_object) {
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

  d_level_min = coarse_level;
  d_level_max = fine_level;
  if (d_level_min == -1) {
    d_level_min = 0;
  }
  if (d_level_max == -1) {
    d_level_max = d_hierarchy->getFinestLevelNumber();
  }

#ifdef DEBUG_CHECK_ASSERTIONS
  if (d_level_min < 0 || d_level_max < 0 || d_level_min > d_level_max) {
    TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
               << "inititialization.\n");
  }
#endif

  if (d_bc_object == &d_simple_bc) {
    d_simple_bc.setHierarchy(d_hierarchy,
                             d_level_min,
                             d_level_max);
  }

  d_fac_ops->set_viscosity_dp_id(cell_viscosity,edge_viscosity,dp);

  createVectorWrappers(p, p_rhs, v, v_rhs);

  d_fac_precond.initializeSolverState(*d_uv, *d_fv);

  d_solver_is_initialized = true;
}
