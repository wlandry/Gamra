/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
Elastic::FACOps::s_side_scratch_var[SAMRAI::MAX_DIM_VAL];

SAMRAI::tbox::StartupShutdownManager::Handler
Elastic::FACOps::s_finalize_handler
(0, 0, 0, Elastic::FACOps::finalizeCallback,
 SAMRAI::tbox::StartupShutdownManager::priorityVariables);

Elastic::FACOps::FACOps
(const SAMRAI::tbox::Dimension& dim,
 const boost::shared_ptr<SAMRAI::tbox::Database> &database,
 const Boundary_Conditions &bc):
  d_dim(dim),
  d_ln_min(-1),
  d_ln_max(-1),
  d_coarse_solver_tolerance(1.e-8),
  d_coarse_solver_max_iterations(10),
  d_residual_tolerance_during_smoothing(-1.0),
  cell_moduli_id(invalid_id),
  edge_moduli_id(invalid_id),
  dv_diagonal_id(invalid_id),
  dv_mixed_id(invalid_id),
  level_set_id(invalid_id),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext("Elastic::FACOps::PRIVATE_CONTEXT")),
  d_side_scratch_id(invalid_id),
  v_refine_patch_strategy("refine patch strategy",bc),
  v_coarsen_patch_strategy("coarsen patch strategy",bc),
  logging(false),
  d_boundary_conditions(bc)
{
  t_restrict_solution = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::restrictSolution()");
  t_restrict_residual = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::restrictResidual()");
  t_prolong = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::prolongErrorAndCorrect()");
  t_smooth_error = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::smoothError()");
  t_solve_coarsest = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::solveCoarsestLevel()");
  t_compute_composite_residual = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::computeCompositeResidualOnLevel()");
  t_compute_residual_norm = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::FACOps::computeResidualNorm()");

  if (d_dim == SAMRAI::tbox::Dimension(1)
      || d_dim > SAMRAI::tbox::Dimension(3))
    { TBOX_ERROR("Elastic::FACOps : DIM == 1 or > 3 not implemented yet.\n"); }

  if (!s_side_scratch_var[dim.getValue() - 1])
    {
      std::ostringstream ss;
      ss << "Elastic::FACOps::private_side_scratch" << dim.getValue();
      s_side_scratch_var[dim.getValue() - 1] =
        boost::make_shared<SAMRAI::pdat::SideVariable<double> >
        (dim, ss.str(),SAMRAI::hier::IntVector::getOne(d_dim));
    }

  d_side_scratch_id = SAMRAI::hier::VariableDatabase::getDatabase()->
    registerVariableAndContext(s_side_scratch_var[dim.getValue() - 1],
                               d_context,
                               SAMRAI::hier::IntVector::getOne(d_dim));

  if (database)
    {
      d_coarse_solver_tolerance =
        database->getDoubleWithDefault("coarse_solver_tolerance",
                                       d_coarse_solver_tolerance);
      d_coarse_solver_max_iterations =
        database->getIntegerWithDefault("coarse_solver_max_iterations",
                                        d_coarse_solver_max_iterations);
      logging = database->getBoolWithDefault("enable_logging", logging);
    }
}
