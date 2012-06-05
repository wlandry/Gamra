/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "ElasticFACOps.h"

#ifndef SAMRAI_INLINE
#include "ElasticFACOps.I"
#endif

namespace SAMRAI {
  namespace solv {

    tbox::Pointer<pdat::CellVariable<double> >
    ElasticFACOps::s_cell_scratch_var[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    tbox::Pointer<pdat::SideVariable<double> >
    ElasticFACOps::s_side_scratch_var[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    tbox::StartupShutdownManager::Handler
    ElasticFACOps::s_finalize_handler(
                                     0,
                                     0,
                                     0,
                                     ElasticFACOps::finalizeCallback,
                                     tbox::StartupShutdownManager::priorityVariables);

    /*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/
    ElasticFACOps::ElasticFACOps(const tbox::Dimension& dim,
                               const std::string& object_name,
                               tbox::Pointer<tbox::Database> database):
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
      cell_moduli_id(invalid_id),
      edge_moduli_id(invalid_id),
      dp_id(invalid_id),
#ifdef HAVE_HYPRE
      d_hypre_solver(dim,
                     object_name + "::hypre_solver",
                     database && database->isDatabase("hypre_solver") ?
                     database->getDatabase("hypre_solver"):
                     tbox::Pointer<tbox::Database>(NULL)),
#endif
      // d_physical_bc_coef(NULL),
      d_context(hier::VariableDatabase::getDatabase()
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
      p_refine_patch_strategy(dim,
                              d_object_name + "::refine patch strategy"),
      v_refine_patch_strategy(dim,
                              d_object_name + "::refine patch strategy"),
      v_coarsen_patch_strategy(dim,
                               d_object_name + "::coarsen patch strategy"),
      d_enable_logging(false),
      d_preconditioner(NULL),
      d_hopscell(),
      d_hopsside()
    {

      t_restrict_solution = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::restrictSolution()");
      t_restrict_residual = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::restrictResidual()");
      t_prolong = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::prolongErrorAndCorrect()");
      t_smooth_error = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::smoothError()");
      t_solve_coarsest = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::solveCoarsestLevel()");
      t_compute_composite_residual = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::computeCompositeResidualOnLevel()");
      t_compute_residual_norm = tbox::TimerManager::getManager()->
        getTimer("solv::ElasticFACOps::computeResidualNorm()");

      if (d_dim == tbox::Dimension(1) || d_dim > tbox::Dimension(3)) {
        TBOX_ERROR("ElasticFACOps : DIM == 1 or > 3 not implemented yet.\n");
      }

      if (s_cell_scratch_var[dim.getValue() - 1].isNull()) {
        TBOX_ASSERT(s_cell_scratch_var[dim.getValue() - 1].isNull());
        TBOX_ASSERT(s_cell_scratch_var[dim.getValue() - 1].isNull());

        std::ostringstream ss;
        ss << "ElasticFACOps::private_cell_scratch" << dim.getValue();
        s_cell_scratch_var[dim.getValue() - 1] = new pdat::CellVariable<double>
          (dim, ss.str());
        ss.str("");
        ss << "ElasticFACOps::private_side_scratch" << dim.getValue();
        s_side_scratch_var[dim.getValue() - 1] = new pdat::SideVariable<double>
          (dim, ss.str());
      }

      hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
      d_cell_scratch_id = vdb->
        registerVariableAndContext(s_cell_scratch_var[dim.getValue() - 1],
                                   d_context,
                                   hier::IntVector::getOne(dim));
      d_side_scratch_id = vdb->
        registerVariableAndContext(s_side_scratch_var[dim.getValue() - 1],
                                   d_context,
                                   hier::IntVector::getOne(d_dim));

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

  }
}
