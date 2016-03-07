/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::ghostfill_nocoarse(int v_id, int dest_level)
{
  if (!ghostfill_nocoarse_schedules[dest_level])
    { TBOX_ERROR("Expected side schedule not found."); }
  SAMRAI::xfer::RefineAlgorithm refiner;
  refiner.registerRefine(v_id,v_id,v_id,
                         boost::shared_ptr<SAMRAI::hier::RefineOperator>());
  refiner.resetSchedule(ghostfill_nocoarse_schedules[dest_level]);
  ghostfill_nocoarse_schedules[dest_level]->fillData(0.0,false);
}
