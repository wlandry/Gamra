/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::restrictResidual
(const SAMRAI::solv::SAMRAIVectorReal<double>& s,
 SAMRAI::solv::SAMRAIVectorReal<double>& d,
 int dest_ln)
{
  t_restrict_residual->start();

  int v_src(s.getComponentDescriptorIndex(0)),
    v_dst(d.getComponentDescriptorIndex(0));

  /// Need to do a sync because the coarsening for v uses ghost zones
  v_coarsen_patch_strategy.data_id=v_src;
  v_coarsen_patch_strategy.is_residual=true;
  xeqScheduleGhostFillNoCoarse(v_src,dest_ln+1);

  xeqScheduleRRestriction(v_dst,v_src,dest_ln);

  t_restrict_residual->stop();
}
