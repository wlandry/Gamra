/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/xfer/PatchLevelFullFillPattern.h>

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::initializeOperatorState
(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& rhs)
{
  deallocateOperatorState();
  boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy
    = solution.getPatchHierarchy();
  int ln;
  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();

  d_ln_min = solution.getCoarsestLevelNumber();
  d_ln_max = solution.getFinestLevelNumber();
  const int v_id(solution.getComponentDescriptorIndex(0)),
    v_rhs_id(rhs.getComponentDescriptorIndex(0));

  /*
   * Initialize the coarse-fine boundary description for the
   * hierarchy.
   */
  d_cf_boundary.resize(hierarchy->getNumberOfLevels());

  SAMRAI::hier::IntVector max_gcw(d_dim, 1);
  for (ln = d_ln_min; ln <= d_ln_max; ++ln)
    d_cf_boundary[ln] = boost::make_shared<SAMRAI::hier::CoarseFineBoundary >
      (*hierarchy, ln, max_gcw);

  v_coarsen_patch_strategy.coarse_fine=d_cf_boundary;
  /*
   * Get the transfer operators.
   */
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (hierarchy->getGridGeometry());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_prolongation_refine_operator =
    geometry->lookupRefineOperator(variable,"V_REFINE");

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_ghostfill_refine_operator =
    geometry->lookupRefineOperator(variable,"V_BOUNDARY_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
  {
    if (!v_prolongation_refine_operator)
      TBOX_ERROR(d_object_name
                 << ": Cannot find v prolongation refine operator");

    if (!v_ghostfill_refine_operator)
      TBOX_ERROR(d_object_name
                 << ": Cannot find ghost filling refinement operator");
  }
#endif

  /*
   * Make space for saving communication schedules.
   * There is no need to delete the old schedules first
   * because we have deallocated the solver state above.
   */
  v_prolongation_refine_schedules.resizeArray(d_ln_max + 1);
  v_ghostfill_refine_schedules.resizeArray(d_ln_max + 1);
  v_nocoarse_refine_schedules.resizeArray(d_ln_max + 1);
  v_urestriction_coarsen_schedules.resizeArray(d_ln_max + 1);
  v_rrestriction_coarsen_schedules.resizeArray(d_ln_max + 1);

  SAMRAI::xfer::RefineAlgorithm v_prolongation_refine_algorithm,
    v_ghostfill_refine_algorithm, v_nocoarse_refine_algorithm;
    
  SAMRAI::xfer::CoarsenAlgorithm v_urestriction_coarsen_algorithm(d_dim),
    v_rrestriction_coarsen_algorithm(d_dim);

  /* This is a little confusing.  The only real purpose here is to
     create a communication schedule.  That communication schedule is
     then reused later when refining, though with a different source,
     scratch, and destination.  So the arguments to registerRefine are
     not all that important, because a different refineAlgorithm will
     be used then. */

  v_prolongation_refine_algorithm.
    registerRefine(d_side_scratch_id,v_id,d_side_scratch_id,
                   v_prolongation_refine_operator);
  v_urestriction_coarsen_algorithm.
    registerCoarsen(v_id,v_id,v_urestriction_coarsen_operator);
  v_rrestriction_coarsen_algorithm.
    registerCoarsen(v_rhs_id,v_rhs_id,v_rrestriction_coarsen_operator);
  v_ghostfill_refine_algorithm.
    registerRefine(v_id,v_id,v_id,
                   v_ghostfill_refine_operator);

  if(have_faults())
    {
      v_ghostfill_refine_algorithm.
        registerRefine(dv_diagonal_id,dv_diagonal_id,dv_diagonal_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      v_ghostfill_refine_algorithm.
        registerRefine(dv_mixed_id,dv_mixed_id,dv_mixed_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
    }

  if(have_embedded_boundary())
    {
      v_ghostfill_refine_algorithm.
        registerRefine(level_set_id,level_set_id,level_set_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      v_prolongation_refine_algorithm.
        registerRefine(level_set_id,level_set_id,level_set_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      v_urestriction_coarsen_algorithm.
        registerCoarsen(level_set_id,level_set_id,
                        boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
      v_rrestriction_coarsen_algorithm.
        registerCoarsen(level_set_id,level_set_id,
                        boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
    }

  v_nocoarse_refine_algorithm.
    registerRefine(v_id,v_id,v_id,
                   boost::shared_ptr<SAMRAI::hier::RefineOperator>());

  /* Refinement and ghost fill operators */
  for (int dest_ln = d_ln_min + 1; dest_ln <= d_ln_max; ++dest_ln)
    {
      boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern>
        fill_pattern(new SAMRAI::xfer::PatchLevelFullFillPattern());
      v_prolongation_refine_schedules[dest_ln] =
        v_prolongation_refine_algorithm.
        createSchedule(fill_pattern,hierarchy->getPatchLevel(dest_ln),
                       boost::shared_ptr<SAMRAI::hier::PatchLevel>(),
                       dest_ln - 1,hierarchy);
                       
    if (!v_prolongation_refine_schedules[dest_ln])
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for v prolongation!\n");
    v_ghostfill_refine_schedules[dest_ln] =
      v_ghostfill_refine_algorithm.
      createSchedule(hierarchy->getPatchLevel(dest_ln),
                     dest_ln - 1,hierarchy,
                     &v_refine_patch_strategy);
    if (!v_ghostfill_refine_schedules[dest_ln])
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling!\n");

    v_nocoarse_refine_schedules[dest_ln] =
      v_nocoarse_refine_algorithm.
      createSchedule(hierarchy->getPatchLevel(dest_ln));
    if (!v_nocoarse_refine_schedules[dest_ln])
      TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                 "ghost filling on bottom level!\n");
    }

  /* Coarsening operators */
  for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln)
    {
      v_urestriction_coarsen_schedules[dest_ln] =
        v_urestriction_coarsen_algorithm.
        createSchedule(hierarchy->getPatchLevel(dest_ln),
                       hierarchy->getPatchLevel(dest_ln + 1),
                       &v_coarsen_patch_strategy);
      if (!v_urestriction_coarsen_schedules[dest_ln])
        TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for "
                   "U v restriction!\n");

      v_rrestriction_coarsen_schedules[dest_ln] =
        v_rrestriction_coarsen_algorithm.
        createSchedule(hierarchy->getPatchLevel(dest_ln),
                       hierarchy->getPatchLevel(dest_ln + 1),
                       &v_coarsen_patch_strategy);
      if (!v_rrestriction_coarsen_schedules[dest_ln])
        TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for "
                   "R v restriction!\n");
    }

  /* Ordinary ghost fill operator on the coarsest level */
  v_nocoarse_refine_schedules[d_ln_min] =
    v_nocoarse_refine_algorithm.
    createSchedule(hierarchy->getPatchLevel(d_ln_min));
  if (!v_nocoarse_refine_schedules[d_ln_min])
    TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for v "
               "ghost filling on bottom level!\n");
  initialized=true;
}
