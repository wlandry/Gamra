/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

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
