/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::refine (int v_dst, int v_src, int v_scr, int dest_level)
{
  if (!refine_schedules[dest_level])
    { TBOX_ERROR("Expected schedule not found."); }
  SAMRAI::xfer::RefineAlgorithm refiner;
  refiner.registerRefine(v_dst, v_src, v_scr, refine_operator);
  if(have_embedded_boundary())
    {
      refiner.registerRefine(level_set_id,level_set_id,level_set_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());
    }
  refiner.resetSchedule(refine_schedules[dest_level]);
  refine_schedules[dest_level]->fillData(0.0,false);
}
