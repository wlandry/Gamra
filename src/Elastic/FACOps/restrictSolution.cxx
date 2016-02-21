/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.hxx"
#include "Elastic/V_Boundary_Refine.hxx"

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
* FACOperatorStrategy virtual restrictSolution function.     *
* After restricting solution, update ghost cells of the affected   *
* level.                                                           *
********************************************************************
*/

void Elastic::FACOps::restrictSolution
(const SAMRAI::solv::SAMRAIVectorReal<double>& s,
 SAMRAI::solv::SAMRAIVectorReal<double>& d,
 int dest_ln)
{
  t_restrict_solution->start();

  int v_src(s.getComponentDescriptorIndex(0)),
    v_dst(d.getComponentDescriptorIndex(0));

  /* Need to do a sync because the coarsening for v uses ghost zones. */
  v_coarsen_patch_strategy.data_id=v_src;
  v_coarsen_patch_strategy.is_residual=false;
  xeqScheduleGhostFillNoCoarse(v_src,dest_ln+1);

  xeqScheduleURestriction(v_dst,v_src,dest_ln);

  boost::shared_ptr<SAMRAI::hier::PatchLevel>
    level = s.getPatchHierarchy()->getPatchLevel(dest_ln);
  v_refine_patch_strategy.is_residual=false;
  v_refine_patch_strategy.data_id=d.getComponentDescriptorIndex(0);
  V_Boundary_Refine::is_residual=false;

  /// FIXME: Is this necessary?
  if (dest_ln == d_ln_min)
    xeqScheduleGhostFillNoCoarse(v_dst,dest_ln);
  else
    xeqScheduleGhostFill(v_dst,dest_ln);

  t_restrict_solution->stop();
}
