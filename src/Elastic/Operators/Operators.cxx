/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
Elastic::Operators::s_side_scratch_var[SAMRAI::MAX_DIM_VAL];

SAMRAI::tbox::StartupShutdownManager::Handler
Elastic::Operators::s_finalize_handler
(0, 0, 0, Elastic::Operators::finalizeCallback,
 SAMRAI::tbox::StartupShutdownManager::priorityVariables);

Elastic::Operators::Operators
(const SAMRAI::tbox::Dimension& dim,
 const boost::shared_ptr<SAMRAI::tbox::Database> &database,
 const Boundary_Conditions &bc):
  dimension(dim),
  level_min(-1),
  level_max(-1),
  coarse_solver_tolerance(1.e-8),
  coarse_solver_max_iterations(10),
  cell_moduli_id(invalid_id),
  edge_moduli_id(invalid_id),
  dv_diagonal_id(invalid_id),
  dv_mixed_id(invalid_id),
  level_set_id(invalid_id),
  context(SAMRAI::hier::VariableDatabase::getDatabase()
          ->getContext("Elastic::Operators::PRIVATE_CONTEXT")),
  side_scratch_id(invalid_id),
  v_refine_patch_strategy("refine patch strategy",bc),
  v_coarsen_patch_strategy("coarsen patch strategy",bc),
  logging(false),
  boundary_conditions(bc)
{
  t_restrict_solution = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::restrictSolution()");
  t_restrict_residual = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::restrictResidual()");
  t_prolong = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::prolongErrorAndCorrect()");
  t_smooth_error = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::smoothError()");
  t_solve_coarsest = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::solveCoarsestLevel()");
  t_compute_composite_residual = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::computeCompositeResidualOnLevel()");
  t_compute_residual_norm = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("solv::Elastic::Operators::computeResidualNorm()");

  if (dimension == SAMRAI::tbox::Dimension(1)
      || dimension > SAMRAI::tbox::Dimension(3))
    { TBOX_ERROR("Elastic::Operators : DIM == 1 or > 3 not implemented yet.\n"); }

  if (!s_side_scratch_var[dim.getValue() - 1])
    {
      std::ostringstream ss;
      ss << "Elastic::Operators::private_side_scratch" << dim.getValue();
      s_side_scratch_var[dim.getValue() - 1] =
        boost::make_shared<SAMRAI::pdat::SideVariable<double> >
        (dim, ss.str(),SAMRAI::hier::IntVector::getOne(dimension));
    }

  side_scratch_id = SAMRAI::hier::VariableDatabase::getDatabase()->
    registerVariableAndContext(s_side_scratch_var[dim.getValue() - 1],
                               context,
                               SAMRAI::hier::IntVector::getOne(dimension));

  if (database)
    {
      coarse_solver_tolerance =
        database->getDoubleWithDefault("coarse_solver_tolerance",
                                       coarse_solver_tolerance);
      coarse_solver_max_iterations =
        database->getIntegerWithDefault("coarse_solver_max_iterations",
                                        coarse_solver_max_iterations);
      logging = database->getBoolWithDefault("enable_logging", logging);
    }
}
