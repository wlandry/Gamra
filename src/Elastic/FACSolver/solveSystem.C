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
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an initialized solver state.                      *
* Before solving, set the final piece of the boundary condition,        *
* which is not known until now, and initialize some internal            *
* solver quantities.                                                    *
*                                                                       *
*************************************************************************
*/

bool Elastic::FACSolver::solveSystem(const int v_id, const int v_rhs_id)
{
  createVectorWrappers(v_id, v_rhs_id);
  return d_fac_precond.solveSystem(*d_uv, *d_fv);
}

/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an uninitialized solver state.                    *
* 1. Initialize the (currently uninitialized) solver state.             *
* 2. Solve.                                                             *
* 3. Deallocate the solver state.                                       *
*                                                                       *
*************************************************************************
*/

bool Elastic::FACSolver::solveSystem
(const int cell_moduli_id,
 const int edge_moduli_id,
 const int dv_diagonal_id,
 const int dv_mixed_id, 
 const int level_set_id,
 const int v_id,
 const int v_rhs_id,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
 hierarchy,
 int coarse_ln,
 int fine_ln)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

  if (d_enable_logging) {
    SAMRAI::tbox::plog << "Elastic::FACSolver::solveSystem (" << d_object_name
               << ")\n";
  }
#ifdef DEBUG_CHECK_ASSERTIONS
  if (d_solver_is_initialized) {
    TBOX_ERROR(
               d_object_name << ".solveSystem(int,int,...): initialized\n"
               <<
               "solver state.  This function can only used when the\n"
               <<
               "solver state is uninitialized.  You should deallocate\n"
               <<
               "the solver state or use solveSystem(int,int).\n");
  }
  if (!hierarchy) {
    TBOX_ERROR(d_object_name << ".solveSystem(): Null hierarchy\n"
               << "specified.\n");
  }
#endif
  initializeSolverState(cell_moduli_id, edge_moduli_id, dv_diagonal_id,
                        dv_mixed_id, level_set_id, v_id, v_rhs_id,
                        hierarchy, coarse_ln, fine_ln);

  bool solver_rval;
  solver_rval = solveSystem(v_id, v_rhs_id);

  deallocateSolverState();

  return solver_rval;
}
