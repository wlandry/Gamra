/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::ghostfill(int v_id, int dest_level)
{
  if (!ghostfill_schedules[dest_level])
    TBOX_ERROR("Expected schedule not found.");
  SAMRAI::xfer::RefineAlgorithm refiner;

  refiner.registerRefine(v_id,v_id,v_id,ghostfill_operator);
  if(have_faults())
    {
      refiner.registerRefine(dv_diagonal_id,dv_diagonal_id,dv_diagonal_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      refiner.registerRefine(dv_mixed_id,dv_mixed_id,dv_mixed_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());
    }
  if(have_embedded_boundary())
    {
      refiner.registerRefine(level_set_id,level_set_id,level_set_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());
    }
  refiner.resetSchedule(ghostfill_schedules[dest_level]);
  ghostfill_schedules[dest_level]->fillData(0.0,false);
}

