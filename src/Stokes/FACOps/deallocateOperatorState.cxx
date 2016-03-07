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

/*
********************************************************************
* FACOperatorStrategy virtual deallocateOperatorState        *
* function.  Deallocate internal hierarchy-dependent data.         *
* State is allocated iff hierarchy is set.                         *
********************************************************************
*/

void Stokes::FACOps::deallocateOperatorState()
{
  if (d_hierarchy) {
    d_cf_boundary.clear();
#ifdef HAVE_HYPRE
    d_hypre_solver.deallocateSolverState();
#endif
    d_hierarchy.reset();
    d_level_min = -1;
    d_level_max = -1;

    p_prolongation_refine_schedules.clear();
    v_prolongation_refine_schedules.clear();
    p_urestriction_coarsen_schedules.clear();
    v_urestriction_coarsen_schedules.clear();
    p_rrestriction_coarsen_schedules.clear();
    v_rrestriction_coarsen_schedules.clear();
    p_ghostfill_refine_schedules.clear();
    v_ghostfill_refine_schedules.clear();
    p_nocoarse_refine_schedules.clear();
    v_nocoarse_refine_schedules.clear();
  }
}
