/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::finalizeCallback()
{
  for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d)
    s_side_scratch_var[d].reset();
}
