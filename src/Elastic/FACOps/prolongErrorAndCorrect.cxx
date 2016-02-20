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

void Elastic::FACOps::prolongErrorAndCorrect
(const SAMRAI::solv::SAMRAIVectorReal<double>& s,
 SAMRAI::solv::SAMRAIVectorReal<double>& d,
 int dest_ln)
{
  t_prolong->start();

  boost::shared_ptr<SAMRAI::hier::PatchLevel> coarse_level =
    d_hierarchy->getPatchLevel(dest_ln - 1);
  boost::shared_ptr<SAMRAI::hier::PatchLevel> fine_level =
    d_hierarchy->getPatchLevel(dest_ln);

  fine_level->allocatePatchData(d_side_scratch_id);

  v_refine_patch_strategy.data_id=d_side_scratch_id;
  v_refine_patch_strategy.is_residual=true;
  V_Boundary_Refine::is_residual=true;
  xeqScheduleProlongation(d_side_scratch_id,
                          s.getComponentDescriptorIndex(0),
                          d_side_scratch_id,
                          dest_ln);

  {
    SAMRAI::math::HierarchySideDataOpsReal<double>
      hierarchy_math_ops(d_hierarchy, dest_ln, dest_ln);
    const int v_dst = d.getComponentDescriptorIndex(0);
    hierarchy_math_ops.add(v_dst, v_dst, d_side_scratch_id);
  }
  fine_level->deallocatePatchData(d_side_scratch_id);

  t_prolong->stop();
}
