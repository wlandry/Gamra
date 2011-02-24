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
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an initialized solver state.                      *
* Before solving, set the final piece of the boundary condition,        *
* which is not known until now, and initialize some internal            *
* solver quantities.                                                    *
*                                                                       *
*************************************************************************
*/

bool SAMRAI::solv::StokesFACSolver::solveSystem(const int p,
                                                const int cell_viscosity,
                                                const int node_viscosity,
                                                const int dp, const int p_rhs,
                                                const int v, const int v_rhs)
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
  if (p < 0 || p_rhs < 0 || v < 0 || v_rhs < 0) {
    TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
  }
#endif
  if (d_bc_object == &d_simple_bc) {
    /*
     * Knowing that we are using the SimpelCellRobinBcCoefsX
     * implementation of RobinBcCoefStrategy, we must save
     * the ghost data in u before solving.
     * The solver overwrites it, but SimpleCellRobinBcCoefs
     * needs to get to access it repeatedly.
     */
    d_simple_bc.cacheDirichletData(p);
  }

  createVectorWrappers(p, p_rhs, v, v_rhs);
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

bool SAMRAI::solv::StokesFACSolver::solveSystem
(const int p,
 const int cell_viscosity,
 const int node_viscosity,
 const int dp,
 const int p_rhs,
 const int v,
 const int v_rhs,
 tbox::Pointer<hier::PatchHierarchy>
 hierarchy,
 int coarse_ln,
 int fine_ln)
{
  TBOX_ASSERT(!hierarchy.isNull());
  TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

  if (d_enable_logging) {
    tbox::plog << "StokesFACSolver::solveSystem (" << d_object_name
               << ")\n";
    d_stokes_spec.printClassData(tbox::plog);
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
  initializeSolverState(p, cell_viscosity, node_viscosity, dp, p_rhs, v, v_rhs,
                        hierarchy, coarse_ln, fine_ln);

  bool solver_rval;
  solver_rval = solveSystem(p, cell_viscosity, node_viscosity,
                            dp, p_rhs, v, v_rhs);

  deallocateSolverState();

  return solver_rval;
}
