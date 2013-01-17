/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.h"

boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
Elastic::FACOps::s_cell_scratch_var[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
Elastic::FACOps::s_side_scratch_var[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

SAMRAI::tbox::StartupShutdownManager::Handler Elastic::FACOps::s_finalize_handler
(0,
 0,
 0,
 Elastic::FACOps::finalizeCallback,
 SAMRAI::tbox::StartupShutdownManager::priorityVariables);

/*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/
Elastic::FACOps::FACOps(const SAMRAI::tbox::Dimension& dim,
                        const std::string& object_name,
                        boost::shared_ptr<SAMRAI::tbox::Database> database,
                        Boundary_Conditions &bc):
  d_dim(dim),
  d_object_name(object_name),
  d_hierarchy(),
  d_ln_min(-1),
  d_ln_max(-1),
  d_cf_boundary(),
  d_smoothing_choice("Tackley"),
  d_coarse_solver_choice("Tackley"),
  d_cf_discretization("Ewing"),
  v_prolongation_method("V_REFINE"),
  d_coarse_solver_tolerance(1.e-8),
  d_coarse_solver_max_iterations(10),
  d_residual_tolerance_during_smoothing(-1.0),
  cell_moduli_id(invalid_id),
  edge_moduli_id(invalid_id),
  dv_diagonal_id(invalid_id),
  dv_perpendicular_id(invalid_id),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::PRIVATE_CONTEXT")),
  d_side_scratch_id(invalid_id),
  v_prolongation_refine_operator(),
  v_prolongation_refine_schedules(),
  v_urestriction_coarsen_operator(),
  v_urestriction_coarsen_schedules(),
  v_rrestriction_coarsen_operator(),
  v_rrestriction_coarsen_schedules(),
  v_ghostfill_refine_operator(),
  v_ghostfill_refine_schedules(),
  v_nocoarse_refine_schedules(),
  v_refine_patch_strategy(dim,
                          d_object_name + "::refine patch strategy",bc),
  v_coarsen_patch_strategy(dim,
                           d_object_name + "::coarsen patch strategy",bc),
  d_enable_logging(false),
  d_preconditioner(),
  d_boundary_conditions(bc),
  d_hopsside()
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
      || d_dim > SAMRAI::tbox::Dimension(3)) {
    TBOX_ERROR("Elastic::FACOps : DIM == 1 or > 3 not implemented yet.\n");
  }

  if (!s_cell_scratch_var[dim.getValue() - 1]) {
    std::ostringstream ss;
    ss << "Elastic::FACOps::private_cell_scratch" << dim.getValue();
    s_cell_scratch_var[dim.getValue() - 1] =
      boost::make_shared<SAMRAI::pdat::CellVariable<double> >(dim, ss.str());
      
    ss.str("");
    ss << "Elastic::FACOps::private_side_scratch" << dim.getValue();
    s_side_scratch_var[dim.getValue() - 1] =
      boost::make_shared<SAMRAI::pdat::SideVariable<double> >(dim, ss.str());
  }

  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();
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

    v_prolongation_method =
      database->getStringWithDefault("v_prolongation_method",
                                     v_prolongation_method);

    d_enable_logging =
      database->getBoolWithDefault("enable_logging",
                                   d_enable_logging);

  }
}
