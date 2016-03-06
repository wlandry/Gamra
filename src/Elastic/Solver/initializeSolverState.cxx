/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Solver.hxx"

void Elastic::Solver::initializeSolverState
(const int cell_moduli_id,
 const int edge_moduli_id,
 const int dv_diagonal_id,
 const int dv_mixed_id, 
 const int level_set_id,
 const int v_id,
 const int v_rhs_id,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &Hierarchy,
 const int coarse_level,
 const int fine_level)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dimension, *hierarchy);

  if (v_id < 0 || v_rhs_id < 0)
    { TBOX_ERROR(__FILE__ << ": Bad patch data id.\n"); }
  hierarchy = Hierarchy;

  level_min = (coarse_level == -1 ? 0 : coarse_level);
  level_max = (fine_level == -1 ? hierarchy->getFinestLevelNumber()
              : fine_level);

  if (level_min < 0 || level_max < 0 || level_min > level_max)
    { TBOX_ERROR(__FILE__ << ": Bad range of levels in\n"
                 << "inititialization.\n"); }

  operators->set_extra_ids(cell_moduli_id,edge_moduli_id,dv_diagonal_id,
                           dv_mixed_id,level_set_id);
  createVectorWrappers(v_id, v_rhs_id);
  preconditioner.initializeSolverState(*uv, *fv);
}
