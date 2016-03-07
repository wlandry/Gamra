/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/xfer/CoarsenAlgorithm.h>

#include "Elastic/Operators.hxx"

void Elastic::Operators::coarsen_resid(int v_dst, int v_src, int dest_level)
{
  if (!coarsen_resid_schedules[dest_level])
    { TBOX_ERROR("Expected schedule not found."); }

  SAMRAI::xfer::CoarsenAlgorithm coarsener(dimension);
  coarsener.registerCoarsen(v_dst,v_src,coarsen_resid_operator);
  if(have_embedded_boundary())
    {
      coarsener.registerCoarsen
        (level_set_id,level_set_id,
         boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
    }
  coarsener.resetSchedule(coarsen_resid_schedules[dest_level]);
  coarsen_resid_schedules[dest_level]->coarsenData();
}

