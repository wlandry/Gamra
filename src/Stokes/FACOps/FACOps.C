/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "Stokes/FACOps.h"

boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
Stokes::FACOps::s_cell_scratch_var[SAMRAI::MAX_DIM_VAL];

boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
Stokes::FACOps::s_side_scratch_var[SAMRAI::MAX_DIM_VAL];

SAMRAI::tbox::StartupShutdownManager::Handler
Stokes::FACOps::s_finalize_handler(
                                   0,
                                   0,
                                   0,
                                   Stokes::FACOps::finalizeCallback,
                                   SAMRAI::tbox::StartupShutdownManager::priorityVariables);

/*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/
Stokes::FACOps::FACOps(const SAMRAI::tbox::Dimension& dim,
                       const std::string& object_name,
                       boost::shared_ptr<SAMRAI::tbox::Database> database):
  d_dim(dim),
  d_object_name(object_name),
  d_hierarchy(),
  d_ln_min(-1),
  d_ln_max(-1),
  d_cf_boundary(),
  d_smoothing_choice("Tackley"),
  d_coarse_solver_choice(
#ifdef HAVE_HYPRE
                         "hypre"
#else
                         "Tackley"
#endif

                         ),
  d_cf_discretization("Ewing"),
  p_prolongation_method("P_REFINE"),
  v_prolongation_method("V_REFINE"),
  p_rrestriction_method("CONSERVATIVE_COARSEN"),
  d_coarse_solver_tolerance(1.e-8),
  d_coarse_solver_max_iterations(10),
  d_residual_tolerance_during_smoothing(-1.0),
  cell_viscosity_id(invalid_id),
  edge_viscosity_id(invalid_id),
  dp_id(invalid_id),
#ifdef HAVE_HYPRE
  d_hypre_solver(dim,
                 object_name + "::hypre_solver",
                 database && database->isDatabase("hypre_solver") ?
                 database->getDatabase("hypre_solver"):
                 boost::shared_ptr<SAMRAI::tbox::Database>()),
#endif
  // d_physical_bc_coef(),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::PRIVATE_CONTEXT")),
  d_cell_scratch_id(invalid_id),
  d_side_scratch_id(invalid_id),
  p_prolongation_refine_operator(),
  p_prolongation_refine_schedules(),
  v_prolongation_refine_operator(),
  v_prolongation_refine_schedules(),
  p_urestriction_coarsen_operator(),
  p_urestriction_coarsen_schedules(),
  v_urestriction_coarsen_operator(),
  v_urestriction_coarsen_schedules(),
  p_rrestriction_coarsen_operator(),
  p_rrestriction_coarsen_schedules(),
  v_rrestriction_coarsen_operator(),
  v_rrestriction_coarsen_schedules(),
  p_ghostfill_refine_operator(),
  p_ghostfill_refine_schedules(),
  v_ghostfill_refine_operator(),
  v_ghostfill_refine_schedules(),
  p_nocoarse_refine_schedules(),
  v_nocoarse_refine_schedules(),
  p_refine_patch_strategy(d_object_name + "::refine patch strategy"),
  v_refine_patch_strategy(d_object_name + "::refine patch strategy"),
  v_coarsen_patch_strategy(d_object_name + "::coarsen patch strategy"),
  d_enable_logging(false),
  d_preconditioner(),
  d_hopscell(),
  d_hopsside()
{

  t_restrict_solution = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::restrictSolution()");
  t_restrict_residual = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::restrictResidual()");
  t_prolong = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::prolongErrorAndCorrect()");
  t_smooth_error = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::smoothError()");
  t_solve_coarsest = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::solveCoarsestLevel()");
  t_compute_composite_residual = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::computeCompositeResidualOnLevel()");
  t_compute_residual_norm = SAMRAI::tbox::TimerManager::getManager()->
    getTimer("Stokes::FACOps::computeResidualNorm()");

  if (d_dim == SAMRAI::tbox::Dimension(1)
      || d_dim > SAMRAI::tbox::Dimension(3)) {
    TBOX_ERROR("Stokes::FACOps : DIM == 1 or > 3 not implemented yet.\n");
  }

  if (!s_cell_scratch_var[dim.getValue() - 1]) {
    std::ostringstream ss;
    ss << "Stokes::FACOps::private_cell_scratch" << dim.getValue();
    s_cell_scratch_var[dim.getValue() - 1] =
      boost::make_shared<SAMRAI::pdat::CellVariable<double> >(dim, ss.str());
    ss.str("");
    ss << "Stokes::FACOps::private_side_scratch" << dim.getValue();
    s_side_scratch_var[dim.getValue() - 1] =
      boost::make_shared<SAMRAI::pdat::SideVariable<double> >
      (dim, ss.str(),SAMRAI::hier::IntVector::getOne(d_dim));
  }

  SAMRAI::hier::VariableDatabase* vdb = SAMRAI::hier::VariableDatabase::getDatabase();
  d_cell_scratch_id = vdb->
    registerVariableAndContext(s_cell_scratch_var[dim.getValue() - 1],
                               d_context,
                               SAMRAI::hier::IntVector::getOne(dim));
  d_side_scratch_id = vdb->
    registerVariableAndContext(s_side_scratch_var[dim.getValue() - 1],
                               d_context,
                               SAMRAI::hier::IntVector::getOne(d_dim));

  /*
   * Some variables initialized by default are overriden by input.
   */

  if (database) {
    d_coarse_solver_choice =
      database->getStringWithDefault("coarse_solver_choice",
                                     d_coarse_solver_choice);
    d_coarse_solver_tolerance =
      database->getDoubleWithDefault("coarse_solver_tolerance",
                                     d_coarse_solver_tolerance);
    d_coarse_solver_max_iterations =
      database->getIntegerWithDefault("coarse_solver_max_iterations",
                                      d_coarse_solver_max_iterations);
    d_smoothing_choice =
      database->getStringWithDefault("smoothing_choice",
                                     d_smoothing_choice);

    d_cf_discretization =
      database->getStringWithDefault("cf_discretization",
                                     d_cf_discretization);

    p_prolongation_method =
      database->getStringWithDefault("p_prolongation_method",
                                     p_prolongation_method);

    v_prolongation_method =
      database->getStringWithDefault("v_prolongation_method",
                                     v_prolongation_method);

    p_rrestriction_method =
      database->getStringWithDefault("p_rrestriction_method",
                                     p_rrestriction_method);

    d_enable_logging =
      database->getBoolWithDefault("enable_logging",
                                   d_enable_logging);

  }
}
