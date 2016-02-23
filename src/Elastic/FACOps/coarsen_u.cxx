/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include "Elastic/FACOps.hxx"

void Elastic::FACOps::coarsen_u(int v_dst, int v_src, int dest_ln)
{
  if (!coarsen_u_schedules[dest_ln])
    TBOX_ERROR("Expected schedule not found.");

  SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
  coarsener.registerCoarsen(v_dst, v_src, coarsen_u_operator);
  if(have_embedded_boundary())
    {
      coarsener.registerCoarsen
        (level_set_id,level_set_id,
         boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
    }
  coarsener.resetSchedule(coarsen_u_schedules[dest_ln]);
  coarsen_u_schedules[dest_ln]->coarsenData();
}
