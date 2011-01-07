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

/*
************************************************************************
* FACOperatorStrategy virtual initializeOperatorState function.  *
*                                                                      *
* Set internal variables to correspond to the solution passed in.      *
* Look up transfer operators.                                          *
************************************************************************
*/

void SAMRAI::solv::StokesFACOps::initializeOperatorState
(const SAMRAIVectorReal<double>& solution,
 const SAMRAIVectorReal<double>& rhs)
{
  deallocateOperatorState();
  int ln;
  hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

  d_hierarchy = solution.getPatchHierarchy();
  d_ln_min = solution.getCoarsestLevelNumber();
  d_ln_max = solution.getFinestLevelNumber();
  d_hopscell = new math::HierarchyCellDataOpsReal<double>(d_hierarchy,
                                                          d_ln_min,
                                                          d_ln_max);
  d_hopsside = new math::HierarchySideDataOpsReal<double>(d_hierarchy,
                                                          d_ln_min,
                                                          d_ln_max);

#ifdef DEBUG_CHECK_ASSERTIONS

  if (d_physical_bc_coef == NULL) {
    /*
     * It's an error not to have bc object set.
     * Note that the bc object cannot be passed in through
     * the argument because the interface is inherited.
     */
    TBOX_ERROR(
               d_object_name << ": No physical bc object in\n"
               <<
               "StokesFACOps::initializeOperatorState\n"
               << "You must use "
               <<
               "StokesFACOps::setPhysicalBcCoefObject\n"
               <<
               "to set one before calling initializeOperatorState\n");
  }

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
  tbox::Pointer<hier::Variable> var;
  {
    vdb->mapIndexToVariable(rhs.getComponentDescriptorIndex(0),
                            var);
    if (!var) {
      TBOX_ERROR(d_object_name << ": RHS component does not\n"
                 << "correspond to a variable.\n");
    }
    tbox::Pointer<pdat::CellVariable<double> > cell_var = var;
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
    tbox::Pointer<pdat::CellVariable<double> > cell_var = var;
    if (!cell_var) {
      TBOX_ERROR(d_object_name
                 << ": Solution variable is not cell-centered double\n");
    }
  }
  for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
    tbox::Pointer<hier::PatchLevel> level_ptr =
      d_hierarchy->getPatchLevel(ln);
    hier::PatchLevel& level = *level_ptr;
    for (hier::PatchLevel::Iterator pi(level); pi; pi++) {
      hier::Patch& patch = **pi;
      tbox::Pointer<hier::PatchData> fd =
        patch.getPatchData(rhs.getComponentDescriptorIndex(0));
      if (fd) {
        /*
         * Some data checks can only be done if the data already exists.
         */
        tbox::Pointer<pdat::CellData<double> > cd = fd;
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
      tbox::Pointer<hier::PatchData> ud =
        patch.getPatchData(solution.getComponentDescriptorIndex(0));
      if (ud) {
        /*
         * Some data checks can only be done if the data already exists.
         */
        tbox::Pointer<pdat::CellData<double> > cd = ud;
        if (!cd) {
          TBOX_ERROR(d_object_name
                     << ": Solution data is not cell-centered double\n");
        }
        if (cd->getDepth() > 1) {
          TBOX_WARNING(d_object_name
                       << ": Solution data has multiple depths.\n"
                       << "Solver is for depth 0 only.\n");
        }
        if (cd->getGhostCellWidth() < hier::IntVector::getOne(d_dim)) {
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
  d_cf_boundary.resizeArray(d_hierarchy->getNumberOfLevels());

  hier::IntVector max_gcw(d_dim, 1);
  for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
    d_cf_boundary[ln] = new hier::CoarseFineBoundary(*d_hierarchy,
                                                     ln,
                                                     max_gcw);
  }
#ifdef HAVE_HYPRE
  if (d_coarse_solver_choice == "hypre") {
    d_hypre_solver.initializeSolverState(d_hierarchy, d_ln_min);
    /*
     * Share the boundary condition object with the hypre solver
     * to make sure that boundary condition settings are consistent
     * between the two objects.
     */
    d_hypre_solver.setPhysicalBcCoefObject(d_physical_bc_coef);
    d_hypre_solver.setMatrixCoefficients(d_stokes_spec);
  }
#endif

  /*
   * Get the transfer operators.
   * Flux coarsening is conservative.
   * Cell (solution, error, etc) coarsening is conservative.
   * Cell refinement from same level is constant refinement.
   * Cell refinement from coarser level is chosen by the
   *   choice of coarse-fine discretization, d_cf_discretization,
   *   which should be set to either "Ewing" or one of the
   *   acceptable strings for looking up the refine operator.
   */
  tbox::Pointer<geom::CartesianGridGeometry> geometry =
    d_hierarchy->getGridGeometry();
  tbox::Pointer<hier::Variable> variable;

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  d_prolongation_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   d_prolongation_method);

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  d_urestriction_coarsen_operator =
    d_rrestriction_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "CONSERVATIVE_COARSEN");

  vdb->mapIndexToVariable(d_oflux_scratch_id, variable);
  d_flux_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "CONSERVATIVE_COARSEN");

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  d_ghostfill_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   d_cf_discretization == "Ewing" ?
                                   "CONSTANT_REFINE" : d_cf_discretization);

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  d_ghostfill_nocoarse_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   "CONSTANT_REFINE");

  vdb->mapIndexToVariable(d_cell_scratch_id, variable);
  p_nocoarse_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   "CONSTANT_REFINE");

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_nocoarse_refine_operator =
    geometry->lookupRefineOperator(variable,
                                   "CONSTANT_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
  if (!d_prolongation_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find prolongation refine operator");
  }
  if (!d_urestriction_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find restriction coarsening operator");
  }
  if (!d_rrestriction_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find restriction coarsening operator");
  }
  if (!d_flux_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find flux coarsening operator");
  }
  if (!d_ghostfill_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find ghost filling refinement operator");
  }
  if (!d_ghostfill_nocoarse_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find ghost filling refinement operator");
  }
  if (!p_nocoarse_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find p ghost filling refinement operator");
  }
  if (!v_nocoarse_refine_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find v ghost filling refinement operator");
  }
#endif

  for (ln = d_ln_min + 1; ln <= d_ln_max; ++ln) {
    d_hierarchy->getPatchLevel(ln)->
      allocatePatchData(d_oflux_scratch_id);
  }

  /*
   * Make space for saving communication schedules.
   * There is no need to delete the old schedules first
   * because we have deallocated the solver state above.
   */
  d_prolongation_refine_schedules.resizeArray(d_ln_max + 1);
  d_ghostfill_refine_schedules.resizeArray(d_ln_max + 1);
  d_ghostfill_nocoarse_refine_schedules.resizeArray(d_ln_max + 1);
  p_nocoarse_refine_schedules.resizeArray(d_ln_max + 1);
  v_nocoarse_refine_schedules.resizeArray(d_ln_max + 1);
  d_urestriction_coarsen_schedules.resizeArray(d_ln_max + 1);
  d_rrestriction_coarsen_schedules.resizeArray(d_ln_max + 1);
  d_flux_coarsen_schedules.resizeArray(d_ln_max + 1);

  d_prolongation_refine_algorithm = new xfer::RefineAlgorithm(d_dim);
  d_urestriction_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  d_rrestriction_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  d_flux_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  d_ghostfill_refine_algorithm = new xfer::RefineAlgorithm(d_dim);
  d_ghostfill_nocoarse_refine_algorithm = new xfer::RefineAlgorithm(d_dim);
  p_nocoarse_refine_algorithm = new xfer::RefineAlgorithm(d_dim);
  v_nocoarse_refine_algorithm = new xfer::RefineAlgorithm(d_dim);

  d_prolongation_refine_algorithm->
    registerRefine(d_cell_scratch_id,
                   solution.getComponentDescriptorIndex(0),
                   d_cell_scratch_id,
                   d_prolongation_refine_operator);
  d_urestriction_coarsen_algorithm->
    registerCoarsen(solution.getComponentDescriptorIndex(0),
                    solution.getComponentDescriptorIndex(0),
                    d_urestriction_coarsen_operator);
  d_rrestriction_coarsen_algorithm->
    registerCoarsen(rhs.getComponentDescriptorIndex(0),
                    rhs.getComponentDescriptorIndex(0),
                    d_rrestriction_coarsen_operator);
  d_ghostfill_refine_algorithm->
    registerRefine(solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   d_ghostfill_refine_operator);
  d_flux_coarsen_algorithm->
    registerCoarsen(((d_flux_id != -1) ? d_flux_id : d_flux_scratch_id),
                    d_oflux_scratch_id,
                    d_flux_coarsen_operator);
  d_ghostfill_nocoarse_refine_algorithm->
    registerRefine(solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   d_ghostfill_nocoarse_refine_operator);
  p_nocoarse_refine_algorithm->
    registerRefine(solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   solution.getComponentDescriptorIndex(0),
                   p_nocoarse_refine_operator);
  v_nocoarse_refine_algorithm->
    registerRefine(solution.getComponentDescriptorIndex(1),
                   solution.getComponentDescriptorIndex(1),
                   solution.getComponentDescriptorIndex(1),
                   v_nocoarse_refine_operator);

  for (int dest_ln = d_ln_min + 1; dest_ln <= d_ln_max; ++dest_ln) {

    tbox::Pointer<xfer::PatchLevelFullFillPattern> fill_pattern(
                                                                new xfer::PatchLevelFullFillPattern());
    d_prolongation_refine_schedules[dest_ln] =
      d_prolongation_refine_algorithm->
      createSchedule(fill_pattern,
                     d_hierarchy->getPatchLevel(dest_ln),
                     tbox::Pointer<hier::PatchLevel>(),
                     dest_ln - 1,
                     d_hierarchy,
                     &d_bc_helper);
    if (!d_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for prolongation!\n");
    }
    d_ghostfill_refine_schedules[dest_ln] =
      d_ghostfill_refine_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     dest_ln - 1,
                     d_hierarchy,
                     &d_bc_helper);
    if (!d_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling!\n");
    }
    d_ghostfill_nocoarse_refine_schedules[dest_ln] =
      d_ghostfill_nocoarse_refine_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     &d_bc_helper);
    if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR(
                 d_object_name
                 <<
                 ": Cannot create a refine schedule for ghost filling on bottom level!\n");
    }
    p_nocoarse_refine_schedules[dest_ln] =
      p_nocoarse_refine_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     &d_bc_helper);
    if (!p_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR(
                 d_object_name
                 <<
                 ": Cannot create a refine schedule for ghost filling on bottom level!\n");
    }
    v_nocoarse_refine_schedules[dest_ln] =
      v_nocoarse_refine_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     &d_bc_helper);
    if (!v_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR(
                 d_object_name
                 <<
                 ": Cannot create a refine schedule for ghost filling on bottom level!\n");
    }
  }
  for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln) {
    d_urestriction_coarsen_schedules[dest_ln] =
      d_urestriction_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!d_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for U restriction!\n");
    }
    d_rrestriction_coarsen_schedules[dest_ln] =
      d_rrestriction_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!d_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for R restriction!\n");
    }
    d_flux_coarsen_schedules[dest_ln] =
      d_flux_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!d_flux_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for flux transfer!\n");
    }
  }
  d_ghostfill_nocoarse_refine_schedules[d_ln_min] =
    d_ghostfill_nocoarse_refine_algorithm->
    createSchedule(d_hierarchy->getPatchLevel(d_ln_min),
                   &d_bc_helper);
  if (!d_ghostfill_nocoarse_refine_schedules[d_ln_min]) {
    TBOX_ERROR(
               d_object_name
               <<
               ": Cannot create a refine schedule for ghost filling on bottom level!\n");
  }
  p_nocoarse_refine_schedules[d_ln_min] =
    p_nocoarse_refine_algorithm->
    createSchedule(d_hierarchy->getPatchLevel(d_ln_min),
                   &d_bc_helper);
  if (!p_nocoarse_refine_schedules[d_ln_min]) {
    TBOX_ERROR(
               d_object_name
               <<
               ": Cannot create a refine schedule for p ghost filling on bottom level!\n");
  }
  v_nocoarse_refine_schedules[d_ln_min] =
    v_nocoarse_refine_algorithm->
    createSchedule(d_hierarchy->getPatchLevel(d_ln_min),
                   &d_bc_helper);
  if (!v_nocoarse_refine_schedules[d_ln_min]) {
    TBOX_ERROR(
               d_object_name
               <<
               ": Cannot create a refine schedule for v ghost filling on bottom level!\n");
  }
}
