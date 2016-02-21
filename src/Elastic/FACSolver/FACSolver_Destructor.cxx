/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

#include "Elastic/FACSolver.hxx"

Elastic::FACSolver::~FACSolver()
{
  s_instance_counter[d_dim.getValue() - 1]--;
  deallocateSolverState();
  if (s_instance_counter[d_dim.getValue() - 1] == 0)
    {
      SAMRAI::hier::VariableDatabase::getDatabase()->
        removeInternalSAMRAIVariablePatchDataIndex
        (s_weight_id[d_dim.getValue() - 1]);
      s_weight_id[d_dim.getValue() - 1] = -1;
    }
}
