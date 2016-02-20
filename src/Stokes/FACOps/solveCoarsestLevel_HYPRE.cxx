/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "Stokes/FACOps.hxx"

#include IOMANIP_HEADER_FILE

#include <SAMRAI/hier/BoundaryBoxUtils.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Index.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include <SAMRAI/pdat/CellDoubleConstantRefine.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/OutersideData.h>
#include <SAMRAI/pdat/OutersideVariable.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/solv/FACPreconditioner.h>
#include "Stokes/HypreSolver.hxx"
#include <SAMRAI/tbox/Array.h>
#include <SAMRAI/tbox/MathUtilities.h>
#include <SAMRAI/tbox/StartupShutdownManager.h>
#include <SAMRAI/tbox/Timer.h>
#include <SAMRAI/tbox/TimerManager.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/tbox/MathUtilities.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/xfer/PatchLevelFullFillPattern.h>

#ifdef HAVE_HYPRE
/*
********************************************************************
* Solve coarsest level using Hypre                                 *
* We only solve for the error, so we always use homogeneous bc.    *
********************************************************************
*/

int Stokes::FACOps::solveCoarsestLevel_HYPRE(
                                             SAMRAI::solv::SAMRAIVectorReal<double>& data,
                                             const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                                             int coarsest_ln) {

  NULL_USE(coarsest_ln);

#ifndef HAVE_HYPRE
  TBOX_ERROR(d_object_name << ": Coarse level solver choice '"
             << d_coarse_solver_choice
             << "' unavailable in "
             << "Stokes::FACOps::solveCoarsestLevel.");

  return 0;

#else

  d_hypre_solver.setStoppingCriteria(d_coarse_solver_max_iterations,
                                     d_coarse_solver_tolerance);
  const int solver_ret =
    d_hypre_solver.solveSystem(
                               data.getComponentDescriptorIndex(0),
                               residual.getComponentDescriptorIndex(0),
                               true);
  /*
   * Present data on the solve.
   * The Hypre solver returns 0 if converged.
   */
  if (d_enable_logging) SAMRAI::tbox::plog
                          << d_object_name << " Hypre solve " << (solver_ret ? "" : "NOT ")
                          << "converged\n"
                          << "\titerations: " << d_hypre_solver.getNumberOfIterations() << "\n"
                          << "\tresidual: " << d_hypre_solver.getRelativeResidualNorm() << "\n";

  return !solver_ret;

#endif

}
#endif

