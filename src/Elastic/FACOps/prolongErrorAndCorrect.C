/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "Elastic/HypreSolver.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

/*
***********************************************************************
* FACOperatorStrategy virtual prolongErrorAndCorrect function.  *
* After the prolongation, we set the physical boundary condition      *
* for the correction, which is zero.  Other ghost cell values,        *
* which are preset to zero, need not be set.                          *
***********************************************************************
*/

void SAMRAI::solv::Elastic::FACOps::prolongErrorAndCorrect
(const SAMRAIVectorReal<double>& s,
 SAMRAIVectorReal<double>& d,
 int dest_ln)
{
  t_prolong->start();

#ifdef DEBUG_CHECK_ASSERTIONS
  if (s.getPatchHierarchy() != d_hierarchy
      || d.getPatchHierarchy() != d_hierarchy) {
    TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
               "internal state hierarchy.");
  }
#endif

  tbox::Pointer<hier::PatchLevel> coarse_level =
    d_hierarchy->getPatchLevel(dest_ln - 1);
  tbox::Pointer<hier::PatchLevel> fine_level =
    d_hierarchy->getPatchLevel(dest_ln);

  /*
   * Data is prolonged into the scratch space corresponding
   * to index d_cell_scratch_id and allocated here.
   */
  fine_level->allocatePatchData(d_cell_scratch_id);
  fine_level->allocatePatchData(d_side_scratch_id);

  // int p_src(s.getComponentDescriptorIndex(0)),
  //   v_src(s.getComponentDescriptorIndex(1)),
  //   p_dst(d.getComponentDescriptorIndex(0)),
  //   v_dst(d.getComponentDescriptorIndex(1));
  // xeqScheduleGhostFillNoCoarse(invalid_id,v_src,dest_ln+1);

  /*
   * Refine solution into scratch space to fill the fine level
   * interior in the scratch space, then use that refined data
   * to correct the fine level error.
   */
  p_refine_patch_strategy.setTargetDataId(d_cell_scratch_id);
  v_refine_patch_strategy.setTargetDataId(d_side_scratch_id);
  // v_refine_patch_strategy.setHomogeneousBc(true);
  xeqScheduleProlongation(d_cell_scratch_id,
                          s.getComponentDescriptorIndex(0),
                          d_cell_scratch_id,
                          d_side_scratch_id,
                          s.getComponentDescriptorIndex(1),
                          d_side_scratch_id,
                          dest_ln);

  set_boundaries(s.getComponentDescriptorIndex(0),
                 s.getComponentDescriptorIndex(1),fine_level,true);
  /*
   * Add the refined error in the scratch space to the error currently
   * residing in the destination level.
   */
  {
    math::HierarchyCellDataOpsReal<double>
      hierarchy_math_ops(d_hierarchy, dest_ln, dest_ln);
    const int p_dst = d.getComponentDescriptorIndex(0);
    hierarchy_math_ops.add(p_dst, p_dst, d_cell_scratch_id);
  }
  {
    math::HierarchySideDataOpsReal<double>
      hierarchy_math_ops(d_hierarchy, dest_ln, dest_ln);
    const int v_dst = d.getComponentDescriptorIndex(1);
    hierarchy_math_ops.add(v_dst, v_dst, d_side_scratch_id);
  }
  fine_level->deallocatePatchData(d_cell_scratch_id);
  fine_level->deallocatePatchData(d_side_scratch_id);

  t_prolong->stop();
}
