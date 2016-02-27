/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

void Elastic::FACSolver::initializeSolverState
(const int cell_moduli_id,
 const int edge_moduli_id,
 const int dv_diagonal_id,
 const int dv_mixed_id, 
 const int level_set_id,
 const int v_id,
 const int v_rhs_id,
 boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
 const int coarse_level,
 const int fine_level)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

  if (v_id < 0 || v_rhs_id < 0)
    { TBOX_ERROR(d_object_name << ": Bad patch data id.\n"); }
  d_hierarchy = hierarchy;

  d_ln_min = (coarse_level == -1 ? 0 : coarse_level);
  d_ln_max = (fine_level == -1 ? d_hierarchy->getFinestLevelNumber()
              : fine_level);

  if (d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max)
    { TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
                 << "inititialization.\n"); }

  d_fac_ops->set_extra_ids(cell_moduli_id,edge_moduli_id,dv_diagonal_id,
                           dv_mixed_id,level_set_id);
  createVectorWrappers(v_id, v_rhs_id);
  d_fac_precond.initializeSolverState(*d_uv, *d_fv);
  d_solver_is_initialized = true;
}
