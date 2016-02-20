/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "Stokes/FACOps.hxx"

#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/xfer/PatchLevelFullFillPattern.h>
/*
************************************************************************
* FACOperatorStrategy virtual initializeOperatorState function.  *
*                                                                      *
* Set internal variables to correspond to the solution passed in.      *
* Look up transfer operators.                                          *
************************************************************************
*/

void Stokes::FACOps::initializeOperatorState
(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& rhs)
{
  deallocateOperatorState();
  int ln;
  SAMRAI::hier::VariableDatabase* vdb = SAMRAI::hier::VariableDatabase::getDatabase();

  d_hierarchy = solution.getPatchHierarchy();
  d_ln_min = solution.getCoarsestLevelNumber();
  d_ln_max = solution.getFinestLevelNumber();
  d_hopscell = boost::make_shared<SAMRAI::math::HierarchyCellDataOpsReal<double> >
    (d_hierarchy,d_ln_min,d_ln_max);
  d_hopsside = boost::make_shared<SAMRAI::math::HierarchySideDataOpsReal<double> >
    (d_hierarchy,d_ln_min,d_ln_max);

#ifdef DEBUG_CHECK_ASSERTIONS

  // if (d_physical_bc_coef == NULL) {
  //   /*
  //    * It's an error not to have bc object set.
  //    * Note that the bc object cannot be passed in through
  //    * the argument because the interface is inherited.
  //    */
  //   TBOX_ERROR(
  //              d_object_name << ": No physical bc object in\n"
  //              <<
  //              "Stokes::FACOps::initializeOperatorState\n"
  //              << "You must use "
  //              <<
  //              "Stokes::FACOps::setPhysicalBcCoefObject\n"
  //              <<
  //              "to set one before calling initializeOperatorState\n");
  // }

  if (solution.getNumberOfComponents() != 1) {
    TBOX_WARNING(d_object_name
                 << ": Solution vector has multiple components.\n"
                 << "Solver is for component 0 only.\n");
  }
  if (rhs.getNumberOfComponents() != 1) {
    TBOX_WARNING(d_object_name
                 << ": RHS vector has multiple components.\n"
                 << "Solver is for component 0 only.\n");
  }

  /*
   * Make sure that solution and rhs data
   *   are of correct type
   *   are allocated
   *   has sufficient ghost width
   */
  boost::shared_ptr<SAMRAI::hier::Variable> var;
  {
    vdb->mapIndexToVariable(rhs.getComponentDescriptorIndex(0),
                            var);
    if (!var) {
      TBOX_ERROR(d_object_name << ": RHS component does not\n"
                 << "correspond to a variable.\n");
    }
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > cell_var =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
      (var);
    if (!cell_var) {
      TBOX_ERROR(d_object_name
                 << ": RHS variable is not cell-centered double\n");
    }
  }
  {
    vdb->mapIndexToVariable(solution.getComponentDescriptorIndex(0),
                            var);
    if (!var) {
      TBOX_ERROR(d_object_name << ": Solution component does not\n"
                 << "correspond to a variable.\n");
    }
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > cell_var =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
      (var);
    if (!cell_var) {
      TBOX_ERROR(d_object_name
                 << ": Solution variable is not cell-centered double\n");
    }
  }
  for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
    boost::shared_ptr<SAMRAI::hier::PatchLevel> level_ptr =
      d_hierarchy->getPatchLevel(ln);
    SAMRAI::hier::PatchLevel& level = *level_ptr;
    for (SAMRAI::hier::PatchLevel::Iterator pi(level.begin()); pi!=level.end(); ++pi) {
      SAMRAI::hier::Patch& patch = **pi;
      boost::shared_ptr<SAMRAI::hier::PatchData> fd =
        patch.getPatchData(rhs.getComponentDescriptorIndex(0));
      if (fd) {
        /*
         * Some data checks can only be done if the data already exists.
         */
        boost::shared_ptr<SAMRAI::pdat::CellData<double> > cd =
          boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
          (fd);
        if (!cd) {
          TBOX_ERROR(d_object_name
                     << ": RHS data is not cell-centered double\n");
        }
        if (cd->getDepth() > 1) {
          TBOX_WARNING(d_object_name
                       << ": RHS data has multiple depths.\n"
                       << "Solver is for depth 0 only.\n");
        }
      }
      boost::shared_ptr<SAMRAI::hier::PatchData> ud =
        patch.getPatchData(solution.getComponentDescriptorIndex(0));
      if (ud) {
        /*
         * Some data checks can only be done if the data already exists.
         */
        boost::shared_ptr<SAMRAI::pdat::CellData<double> > cd =
          boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
          (ud);
        if (!cd) {
          TBOX_ERROR(d_object_name
                     << ": Solution data is not cell-centered double\n");
        }
        if (cd->getDepth() > 1) {
          TBOX_WARNING(d_object_name
                       << ": Solution data has multiple depths.\n"
                       << "Solver is for depth 0 only.\n");
        }
        if (cd->getGhostCellWidth() < SAMRAI::hier::IntVector::getOne(d_dim)) {
          TBOX_ERROR(d_object_name
                     << ": Solution data has insufficient ghost width\n");
        }
      }
    }
  }

  /*
   * Solution and rhs must have some similar properties.
   */
  if (rhs.getPatchHierarchy() != d_hierarchy
      || rhs.getCoarsestLevelNumber() != d_ln_min
      || rhs.getFinestLevelNumber() != d_ln_max) {
    TBOX_ERROR(d_object_name << ": solution and rhs do not have\n"
               << "the same set of patch levels.\n");
  }

#endif

  /*
   * Initialize the coarse-fine boundary description for the
   * hierarchy.
   */
  d_cf_boundary.resize(d_hierarchy->getNumberOfLevels());

  SAMRAI::hier::IntVector max_gcw(d_dim, 1);
  for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
    d_cf_boundary[ln] = boost::make_shared<SAMRAI::hier::CoarseFineBoundary>
      (*d_hierarchy,ln,max_gcw);
  }

  v_coarsen_patch_strategy.coarse_fine=d_cf_boundary;
// #ifdef HAVE_HYPRE
//   if (d_coarse_solver_choice == "hypre") {
//     d_hypre_solver.initializeSolverState(d_hierarchy, d_ln_min);
//     /*
//      * Share the boundary condition object with the hypre solver
//      * to make sure that boundary condition settings are consistent
//      * between the two objects.
//      */
//     d_hypre_solver.setPhysicalBcCoefObject(d_physical_bc_coef);
//     d_hypre_solver.setMatrixCoefficients(d_stokes_spec);
//   }
// #endif

  /*
   * Get the transfer operators.
   * Cell (solution, error, etc) coarsening is conservative.
   */
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (d_hierarchy->getGridGeometry());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  p_prolongation_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   p_prolongation_method);

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_prolongation_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   v_prolongation_method);

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  p_urestriction_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "CONSERVATIVE_COARSEN");
  p_rrestriction_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                   p_rrestriction_method);

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_urestriction_coarsen_operator =
    v_rrestriction_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "V_COARSEN");

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  p_ghostfill_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   "P_BOUNDARY_REFINE");

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_ghostfill_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   "V_BOUNDARY_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
  if (!p_prolongation_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find p prolongation refine operator");
  }
  if (!v_prolongation_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find v prolongation refine operator");
  }
  if (!p_urestriction_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find p restriction coarsening operator");
  }
  if (!v_urestriction_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find v restriction coarsening operator");
  }
  if (!p_rrestriction_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find p restriction coarsening operator");
  }
  if (!v_rrestriction_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find v restriction coarsening operator");
  }
  if (!p_ghostfill_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find ghost filling refinement operator");
  }
  if (!v_ghostfill_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find ghost filling refinement operator");
  }
#endif

  /*
   * Make space for saving communication schedules.
   * There is no need to delete the old schedules first
   * because we have deallocated the solver state above.
   */
  p_prolongation_refine_schedules.resize(d_ln_max + 1);
  v_prolongation_refine_schedules.resize(d_ln_max + 1);
  p_ghostfill_refine_schedules.resize(d_ln_max + 1);
  v_ghostfill_refine_schedules.resize(d_ln_max + 1);
  p_nocoarse_refine_schedules.resize(d_ln_max + 1);
  v_nocoarse_refine_schedules.resize(d_ln_max + 1);
  p_urestriction_coarsen_schedules.resize(d_ln_max + 1);
  p_rrestriction_coarsen_schedules.resize(d_ln_max + 1);
  v_urestriction_coarsen_schedules.resize(d_ln_max + 1);
  v_rrestriction_coarsen_schedules.resize(d_ln_max + 1);

  SAMRAI::xfer::RefineAlgorithm p_prolongation_refine_algorithm,
    v_prolongation_refine_algorithm,
    p_ghostfill_refine_algorithm,
    v_ghostfill_refine_algorithm,
    p_nocoarse_refine_algorithm,
    v_nocoarse_refine_algorithm;
  SAMRAI::xfer::CoarsenAlgorithm p_urestriction_coarsen_algorithm(d_dim),
    p_rrestriction_coarsen_algorithm(d_dim),
    v_urestriction_coarsen_algorithm(d_dim),
    v_rrestriction_coarsen_algorithm(d_dim);

  /* This is a little confusing.  The only real purpose here is to
     create a communication schedule.  That communication schedule is
     then reused later when refining, though with a different source,
     scratch, and destination.  So the arguments to registerRefine are
     not all that important, because a different refineAlgorithm will
     be used then. */

  p_prolongation_refine_algorithm.
    registerRefine(d_cell_scratch_id,
                   solution.getComponentDescriptorIndex(0),
                   d_cell_scratch_id,
                   p_prolongation_refine_operator);
  v_prolongation_refine_algorithm.
    registerRefine(d_side_scratch_id,
                   solution.getComponentDescriptorIndex(1),
                   d_side_scratch_id,
                   v_prolongation_refine_operator);
  p_urestriction_coarsen_algorithm.
    registerCoarsen(solution.getComponentDescriptorIndex(0),
                    solution.getComponentDescriptorIndex(0),
                    p_urestriction_coarsen_operator);
  p_rrestriction_coarsen_algorithm.
    registerCoarsen(rhs.getComponentDescriptorIndex(0),
                    rhs.getComponentDescriptorIndex(0),
                    p_rrestriction_coarsen_operator);
  v_urestriction_coarsen_algorithm.
    registerCoarsen(solution.getComponentDescriptorIndex(1),
                    solution.getComponentDescriptorIndex(1),
                    v_urestriction_coarsen_operator);
  v_rrestriction_coarsen_algorithm.
    registerCoarsen(rhs.getComponentDescriptorIndex(1),
                    rhs.getComponentDescriptorIndex(1),
                    v_rrestriction_coarsen_operator);
  p_ghostfill_refine_algorithm.
    registerRefine(solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   p_ghostfill_refine_operator);
  v_ghostfill_refine_algorithm.
    registerRefine(solution.getComponentDescriptorIndex(1),
                   solution.getComponentDescriptorIndex(1),
                   solution.getComponentDescriptorIndex(1),
                   v_ghostfill_refine_operator);
  p_nocoarse_refine_algorithm.
    registerRefine(solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   boost::shared_ptr<SAMRAI::hier::RefineOperator>());
  v_nocoarse_refine_algorithm.
    registerRefine(solution.getComponentDescriptorIndex(1),
                   solution.getComponentDescriptorIndex(1),
                   solution.getComponentDescriptorIndex(1),
                   boost::shared_ptr<SAMRAI::hier::RefineOperator>());

  /* Refinement and ghost fill operators */
  for (int dest_ln = d_ln_min + 1; dest_ln <= d_ln_max; ++dest_ln) {

    boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern>
      fill_pattern(new SAMRAI::xfer::PatchLevelFullFillPattern());
    p_prolongation_refine_schedules[dest_ln] =
      p_prolongation_refine_algorithm.
      createSchedule(fill_pattern,
                     d_hierarchy->getPatchLevel(dest_ln),
                     boost::shared_ptr<SAMRAI::hier::PatchLevel>(),
                     dest_ln - 1,
                     d_hierarchy);
    if (!p_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for p prolongation!\n");
    }
    v_prolongation_refine_schedules[dest_ln] =
      v_prolongation_refine_algorithm.
      createSchedule(fill_pattern,
                     d_hierarchy->getPatchLevel(dest_ln),
                     boost::shared_ptr<SAMRAI::hier::PatchLevel>(),
                     dest_ln - 1,
                     d_hierarchy);
    if (!v_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for v prolongation!\n");
    }
    p_ghostfill_refine_schedules[dest_ln] =
      p_ghostfill_refine_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     dest_ln - 1,
                     d_hierarchy,
                     &p_refine_patch_strategy);
    if (!p_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling!\n");
    }
    v_ghostfill_refine_schedules[dest_ln] =
      v_ghostfill_refine_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     dest_ln - 1,
                     d_hierarchy,
                     &v_refine_patch_strategy);
    if (!v_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling!\n");
    }
    p_nocoarse_refine_schedules[dest_ln] =
      p_nocoarse_refine_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln));
    if (!p_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling on bottom level!\n");
    }
    v_nocoarse_refine_schedules[dest_ln] =
      v_nocoarse_refine_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln));
    if (!v_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling on bottom level!\n");
    }
  }

  /* Coarsening operators */
  for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln) {
    p_urestriction_coarsen_schedules[dest_ln] =
      p_urestriction_coarsen_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!p_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for U p restriction!\n");
    }
    p_rrestriction_coarsen_schedules[dest_ln] =
      p_rrestriction_coarsen_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!p_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for R p restriction!\n");
    }

    v_urestriction_coarsen_schedules[dest_ln] =
      v_urestriction_coarsen_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1),
                     &v_coarsen_patch_strategy);
    if (!v_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for U v restriction!\n");
    }
    v_rrestriction_coarsen_schedules[dest_ln] =
      v_rrestriction_coarsen_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1),
                     &v_coarsen_patch_strategy);
    if (!v_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for R v restriction!\n");
    }
  }

  /* Ordinary ghost fill operator on the coarsest level */
  p_nocoarse_refine_schedules[d_ln_min] =
    p_nocoarse_refine_algorithm.
    createSchedule(d_hierarchy->getPatchLevel(d_ln_min));
  if (!p_nocoarse_refine_schedules[d_ln_min]) {
    TBOX_ERROR(
               d_object_name
               <<
               ": Cannot create a refine schedule for p ghost filling on bottom level!\n");
  }
  v_nocoarse_refine_schedules[d_ln_min] =
    v_nocoarse_refine_algorithm.
    createSchedule(d_hierarchy->getPatchLevel(d_ln_min));
  if (!v_nocoarse_refine_schedules[d_ln_min]) {
    TBOX_ERROR(
               d_object_name
               <<
               ": Cannot create a refine schedule for v ghost filling on bottom level!\n");
  }
}
