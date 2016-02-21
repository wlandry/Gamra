/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

void Elastic::FACSolver::initializeStatics()
{
  for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d)
    {
      s_weight_id[d] = -1;
      s_instance_counter[d] = -1;
    }
  s_initialized = true;
}
