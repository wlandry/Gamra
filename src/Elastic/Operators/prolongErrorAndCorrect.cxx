/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"
#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"

void Elastic::Operators::prolongErrorAndCorrect
(const SAMRAI::solv::SAMRAIVectorReal<double>& s,
 SAMRAI::solv::SAMRAIVectorReal<double>& d,
 int dest_level)
{
  t_prolong->start();

  SAMRAI::hier::PatchLevel &fine_level =
    *(s.getPatchHierarchy()->getPatchLevel(dest_level));

  // FIXME: Is there a way to do this without having to allocate a
  // scratch space?
  fine_level.allocatePatchData(side_scratch_id);

  v_refine_patch_strategy.data_id=side_scratch_id;
  v_refine_patch_strategy.is_residual=true;
  /// FIXME: Get rid of this global variable
  Coarse_Fine_Boundary_Refine::is_residual=true;
  refine(side_scratch_id, s.getComponentDescriptorIndex(0),
         side_scratch_id, dest_level);
  {
    SAMRAI::math::HierarchySideDataOpsReal<double>
      hierarchy_math_ops(s.getPatchHierarchy(), dest_level, dest_level);
    const int v_dst = d.getComponentDescriptorIndex(0);
    hierarchy_math_ops.add(v_dst, v_dst, side_scratch_id);
  }
  fine_level.deallocatePatchData(side_scratch_id);

  t_prolong->stop();
}
