/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "StokesFACOps.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "StokesHypreSolver.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

#ifndef SAMRAI_INLINE
#include "StokesFACOps.I"
#endif

namespace SAMRAI {
  namespace solv {

    tbox::Pointer<pdat::CellVariable<double> >
    StokesFACOps::s_cell_scratch_var[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    tbox::Pointer<pdat::SideVariable<double> >
    StokesFACOps::s_flux_scratch_var[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    tbox::Pointer<pdat::SideVariable<double> >
    StokesFACOps::s_side_scratch_var[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    tbox::Pointer<pdat::OutersideVariable<double> >
    StokesFACOps::s_oflux_scratch_var[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    tbox::StartupShutdownManager::Handler
    StokesFACOps::s_finalize_handler(
                                     0,
                                     0,
                                     0,
                                     StokesFACOps::finalizeCallback,
                                     tbox::StartupShutdownManager::priorityVariables);

    /*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/
    StokesFACOps::StokesFACOps(const tbox::Dimension& dim,
                               const std::string& object_name,
                               tbox::Pointer<tbox::Database> database):
      d_dim(dim),
      d_object_name(object_name),
      d_hierarchy(),
      d_ln_min(-1),
      d_ln_max(-1),
      d_cf_boundary(),
      d_stokes_spec(object_name + "::Stokes specs"),
      d_smoothing_choice("redblack"),
      d_coarse_solver_choice(
#ifdef HAVE_HYPRE
                             "hypre"
#else
                             "redblack"
#endif

                             ),
      d_cf_discretization("Ewing"),
      p_prolongation_method("P_REFINE"),
      v_prolongation_method("V_REFINE"),
      d_coarse_solver_tolerance(1.e-8),
      d_coarse_solver_max_iterations(10),
      d_residual_tolerance_during_smoothing(-1.0),
      cell_viscosity_id(-1),
      edge_viscosity_id(-1),
      dp_id(-1),
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
      d_cell_scratch_id(-1),
      d_side_scratch_id(-1),
      d_flux_scratch_id(-1),
      d_oflux_scratch_id(-1),
      invalid_id(-1),
      p_prolongation_refine_operator(),
      p_prolongation_refine_algorithm(),
      p_prolongation_refine_schedules(),
      v_prolongation_refine_operator(),
      v_prolongation_refine_algorithm(),
      v_prolongation_refine_schedules(),
      p_urestriction_coarsen_operator(),
      p_urestriction_coarsen_algorithm(),
      p_urestriction_coarsen_schedules(),
      v_urestriction_coarsen_operator(),
      v_urestriction_coarsen_algorithm(),
      v_urestriction_coarsen_schedules(),
      p_rrestriction_coarsen_operator(),
      p_rrestriction_coarsen_algorithm(),
      p_rrestriction_coarsen_schedules(),
      v_rrestriction_coarsen_operator(),
      v_rrestriction_coarsen_algorithm(),
      v_rrestriction_coarsen_schedules(),
      d_flux_coarsen_operator(),
      d_flux_coarsen_algorithm(),
      d_flux_coarsen_schedules(),
      p_ghostfill_refine_operator(),
      p_ghostfill_refine_algorithm(),
      p_ghostfill_refine_schedules(),
      v_ghostfill_refine_operator(),
      v_ghostfill_refine_algorithm(),
      v_ghostfill_refine_schedules(),
      p_nocoarse_refine_operator(),
      p_nocoarse_refine_algorithm(),
      p_nocoarse_refine_schedules(),
      v_nocoarse_refine_operator(),
      v_nocoarse_refine_algorithm(),
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
        getTimer("solv::StokesFACOps::restrictSolution()");
      t_restrict_residual = tbox::TimerManager::getManager()->
        getTimer("solv::StokesFACOps::restrictResidual()");
      t_prolong = tbox::TimerManager::getManager()->
        getTimer("solv::StokesFACOps::prolongErrorAndCorrect()");
      t_smooth_error = tbox::TimerManager::getManager()->
        getTimer("solv::StokesFACOps::smoothError()");
      t_solve_coarsest = tbox::TimerManager::getManager()->
        getTimer("solv::StokesFACOps::solveCoarsestLevel()");
      t_compute_composite_residual = tbox::TimerManager::getManager()->
        getTimer("solv::StokesFACOps::computeCompositeResidualOnLevel()");
      t_compute_residual_norm = tbox::TimerManager::getManager()->
        getTimer("solv::StokesFACOps::computeResidualNorm()");

      if (d_dim == tbox::Dimension(1) || d_dim > tbox::Dimension(3)) {
        TBOX_ERROR("StokesFACOps : DIM == 1 or > 3 not implemented yet.\n");
      }

      if (s_cell_scratch_var[dim.getValue() - 1].isNull()) {
        TBOX_ASSERT(s_cell_scratch_var[dim.getValue() - 1].isNull());
        TBOX_ASSERT(s_cell_scratch_var[dim.getValue() - 1].isNull());

        std::ostringstream ss;
        ss << "StokesFACOps::private_cell_scratch" << dim.getValue();
        s_cell_scratch_var[dim.getValue() - 1] = new pdat::CellVariable<double>
          (dim, ss.str());
        ss.str("");
        ss << "StokesFACOps::private_flux_scratch" << dim.getValue();
        s_flux_scratch_var[dim.getValue() - 1] = new pdat::SideVariable<double>
          (dim, ss.str());
        ss.str("");
        ss << "StokesFACOps::private_side_scratch" << dim.getValue();
        s_side_scratch_var[dim.getValue() - 1] = new pdat::SideVariable<double>
          (dim, ss.str());
        ss.str("");
        ss << "StokesFACOps::private_oflux_scratch" << dim.getValue();
        s_oflux_scratch_var[dim.getValue() - 1] = new pdat::OutersideVariable<double>
          (dim, ss.str());
      }

      hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
      d_cell_scratch_id = vdb->
        registerVariableAndContext(s_cell_scratch_var[dim.getValue() - 1],
                                   d_context,
                                   hier::IntVector::getOne(dim));
      d_flux_scratch_id = vdb->
        registerVariableAndContext(s_flux_scratch_var[dim.getValue() - 1],
                                   d_context,
                                   hier::IntVector::getZero(d_dim));
      d_side_scratch_id = vdb->
        registerVariableAndContext(s_side_scratch_var[dim.getValue() - 1],
                                   d_context,
                                   hier::IntVector::getOne(d_dim));
      d_oflux_scratch_id = vdb->
        registerVariableAndContext(s_oflux_scratch_var[dim.getValue() - 1],
                                   d_context,
                                   hier::IntVector::getZero(d_dim));

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

        d_enable_logging =
          database->getBoolWithDefault("enable_logging",
                                       d_enable_logging);

      }

      /*
       * Check input validity and correctness.
       */
      checkInputPatchDataIndices();

    }

  }
}
