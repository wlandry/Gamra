/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::deallocateOperatorState()
{
  coarse_fine_boundary.clear();
  level_min = -1;
  level_max = -1;

  refine_schedules.clear();
  coarsen_solution_schedules.clear();
  coarsen_resid_schedules.clear();
  ghostfill_schedules.clear();
  ghostfill_nocoarse_schedules.clear();
}
