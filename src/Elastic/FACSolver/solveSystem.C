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
#ifdef DEBUG_CHECK_ASSERTIONS
  if (!d_solver_is_initialized) {
    TBOX_ERROR(
               d_object_name << ".solveSystem(int,int): uninitialized\n"
               <<
               "solver state.  You must call initializeSolverState()\n"
               <<
               "before using this function.  Or you can use\n"
               <<
               "solveSystem(int,int,...) to initialize the solver,\n"
               << "solve and deallocate the solver.\n");
  }
  if (v_id < 0 || v_rhs_id < 0) {
    TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
  }
#endif

  createVectorWrappers(v_id, v_rhs_id);
  bool solver_rval;

  solver_rval = d_fac_precond.solveSystem(*d_uv, *d_fv);

  return solver_rval;
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
 const bool &have_faults,
 const int v_id,
 const int v_rhs_id,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
 hierarchy,
 int coarse_ln,
 int fine_ln)
{
  TBOX_ASSERT(hierarchy);
  TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

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
                        dv_mixed_id, have_faults, v_id, v_rhs_id,
                        hierarchy, coarse_ln, fine_ln);

  bool solver_rval;
  solver_rval = solveSystem(v_id, v_rhs_id);

  deallocateSolverState();

  return solver_rval;
}
