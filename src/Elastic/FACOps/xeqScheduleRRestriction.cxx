/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/xfer/CoarsenAlgorithm.h>

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::xeqScheduleRRestriction(int v_dst, int v_src,
                                              int dest_ln)
{
  /* v */
  if (!v_rrestriction_coarsen_schedules[dest_ln]) {
    TBOX_ERROR("Expected schedule not found.");
  }

  SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
  coarsener.registerCoarsen(v_dst,v_src,v_rrestriction_coarsen_operator);
  if(have_embedded_boundary())
    {
      coarsener.registerCoarsen
        (level_set_id,level_set_id,
         boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
    }
  coarsener.resetSchedule(v_rrestriction_coarsen_schedules[dest_ln]);
  v_rrestriction_coarsen_schedules[dest_ln]->coarsenData();
}

