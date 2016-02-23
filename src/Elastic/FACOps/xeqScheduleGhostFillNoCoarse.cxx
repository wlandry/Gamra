/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::xeqScheduleGhostFillNoCoarse(int v_id, int dest_ln)
{
  if (!ghostfill_nocoarse_schedules[dest_ln])
    TBOX_ERROR("Expected side schedule not found.");
  SAMRAI::xfer::RefineAlgorithm refiner;
  refiner.registerRefine(v_id,v_id,v_id,
                         boost::shared_ptr<SAMRAI::hier::RefineOperator>());
  refiner.resetSchedule(ghostfill_nocoarse_schedules[dest_ln]);
  ghostfill_nocoarse_schedules[dest_ln]->fillData(0.0,false);
}
