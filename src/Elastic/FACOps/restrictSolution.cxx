/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"
#include "Elastic/V_Boundary_Refine.hxx"

void Elastic::FACOps::restrictSolution
(const SAMRAI::solv::SAMRAIVectorReal<double>& s,
 SAMRAI::solv::SAMRAIVectorReal<double>& d,
 int dest_ln)
{
  t_restrict_solution->start();

  int v_src(s.getComponentDescriptorIndex(0)),
    v_dst(d.getComponentDescriptorIndex(0));

  /// Need to do a sync because the coarsening for v uses ghost zones.
  v_coarsen_patch_strategy.data_id=v_src;
  v_coarsen_patch_strategy.is_residual=false;
  ghostfill_nocoarse(v_src,dest_ln+1);
  coarsen_solution(v_dst,v_src,dest_ln);

  boost::shared_ptr<SAMRAI::hier::PatchLevel>
    level = s.getPatchHierarchy()->getPatchLevel(dest_ln);
  v_refine_patch_strategy.is_residual=false;
  v_refine_patch_strategy.data_id=d.getComponentDescriptorIndex(0);
  V_Boundary_Refine::is_residual=false;

  /// FIXME: Is this sync necessary?  The style is to sync before use,
  /// not after something is set.
  if (dest_ln == d_ln_min)
    ghostfill_nocoarse(v_dst,dest_ln);
  else
    ghostfill(v_dst,dest_ln);

  t_restrict_solution->stop();
}
