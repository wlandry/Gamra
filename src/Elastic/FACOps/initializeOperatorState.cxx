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
  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();

  level_min = solution.getCoarsestLevelNumber();
  level_max = solution.getFinestLevelNumber();
  const int v_id(solution.getComponentDescriptorIndex(0)),
    v_rhs_id(rhs.getComponentDescriptorIndex(0));

  /// Initialize the coarse-fine boundary description for the hierarchy.
  coarse_fine_boundary.resize(hierarchy->getNumberOfLevels());

  SAMRAI::hier::IntVector max_gcw(dimension, 1);
  for (int ln = level_min; ln <= level_max; ++ln)
    {
      coarse_fine_boundary[ln]
        = boost::make_shared<SAMRAI::hier::CoarseFineBoundary >
        (*hierarchy, ln, max_gcw);
    }

  v_coarsen_patch_strategy.coarse_fine=coarse_fine_boundary;
  /// Get the transfer operators.
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (hierarchy->getGridGeometry());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  vdb->mapIndexToVariable(side_scratch_id, variable);
  refine_operator =
    geometry->lookupRefineOperator(variable,"V_REFINE");
  if (!refine_operator)
    { TBOX_ERROR(__FILE__
                 << ": Cannot find v prolongation refine operator"); }

  vdb->mapIndexToVariable(side_scratch_id, variable);
  ghostfill_operator = 
    geometry->lookupRefineOperator(variable, "COARSE_FINE_BOUNDARY_REFINE");
  if (!ghostfill_operator)
    { TBOX_ERROR(__FILE__
                 << ": Cannot find ghost filling refinement operator"); }

  /// Make space for saving communication schedules. There is no need
  /// to delete the old schedules first because we have deallocated
  /// the solver state above.

  refine_schedules.resize(level_max + 1);
  ghostfill_schedules.resize(level_max + 1);
  ghostfill_nocoarse_schedules.resize(level_max + 1);
  coarsen_solution_schedules.resize(level_max + 1);
  coarsen_resid_schedules.resize(level_max + 1);

  SAMRAI::xfer::RefineAlgorithm refine_algorithm,
    ghostfill_algorithm, ghostfill_nocoarse_algorithm;
    
  SAMRAI::xfer::CoarsenAlgorithm coarsen_solution_algorithm(dimension),
    coarsen_resid_algorithm(dimension);

  /// This is a little confusing.  The only real purpose here is to
  /// create a communication schedule.  That communication schedule is
  /// then reused later when refining, though with a different source,
  /// scratch, and destination.  So the arguments to registerRefine
  /// are not all that important, because a different refineAlgorithm
  /// will be used then.

  refine_algorithm.registerRefine(side_scratch_id,v_id,side_scratch_id,
                                  refine_operator);
  coarsen_solution_algorithm.registerCoarsen(v_id,v_id,
                                             coarsen_solution_operator);
  coarsen_resid_algorithm.registerCoarsen(v_rhs_id,v_rhs_id,
                                          coarsen_resid_operator);
  ghostfill_algorithm.registerRefine(v_id,v_id,v_id,ghostfill_operator);
                                     

  if(have_faults())
    {
      ghostfill_algorithm.
        registerRefine(dv_diagonal_id,dv_diagonal_id,dv_diagonal_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      ghostfill_algorithm.
        registerRefine(dv_mixed_id,dv_mixed_id,dv_mixed_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
    }

  if(have_embedded_boundary())
    {
      ghostfill_algorithm.
        registerRefine(level_set_id,level_set_id,level_set_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      refine_algorithm.
        registerRefine(level_set_id,level_set_id,level_set_id,
                       boost::shared_ptr<SAMRAI::hier::RefineOperator>());
      coarsen_solution_algorithm.
        registerCoarsen(level_set_id,level_set_id,
                        boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
      coarsen_resid_algorithm.
        registerCoarsen(level_set_id,level_set_id,
                        boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());
    }

  ghostfill_nocoarse_algorithm.
    registerRefine(v_id,v_id,v_id,
                   boost::shared_ptr<SAMRAI::hier::RefineOperator>());

  /// Refinement and ghost fill operators
  for (int dest_ln = level_min + 1; dest_ln <= level_max; ++dest_ln)
    {
      boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern>
        fill_pattern
        (boost::make_shared<SAMRAI::xfer::PatchLevelFullFillPattern>());

      refine_schedules[dest_ln] =
        refine_algorithm.
        createSchedule(fill_pattern,hierarchy->getPatchLevel(dest_ln),
                       boost::shared_ptr<SAMRAI::hier::PatchLevel>(),
                       dest_ln - 1,hierarchy);
      if (!refine_schedules[dest_ln])
        { TBOX_ERROR(__FILE__
                     << ": Cannot create a refine schedule for refining!\n"); }

      ghostfill_schedules[dest_ln] = ghostfill_algorithm.
        createSchedule(hierarchy->getPatchLevel(dest_ln),
                       dest_ln - 1,hierarchy,
                       &v_refine_patch_strategy);
      if (!ghostfill_schedules[dest_ln])
        { TBOX_ERROR(__FILE__
                     << ": Cannot create a refine schedule for "
                     "ghost filling!\n"); }

      ghostfill_nocoarse_schedules[dest_ln] = ghostfill_nocoarse_algorithm.
        createSchedule(hierarchy->getPatchLevel(dest_ln));
      if (!ghostfill_nocoarse_schedules[dest_ln])
        { TBOX_ERROR(__FILE__ << ": Cannot create a refine schedule for "
                     "ghost filling on bottom level!\n"); }
    }

  /// Coarsening operators
  for (int dest_ln = level_min; dest_ln < level_max; ++dest_ln)
    {
      coarsen_solution_schedules[dest_ln] =
        coarsen_solution_algorithm.
        createSchedule(hierarchy->getPatchLevel(dest_ln),
                       hierarchy->getPatchLevel(dest_ln + 1),
                       &v_coarsen_patch_strategy);
      if (!coarsen_solution_schedules[dest_ln])
        TBOX_ERROR(__FILE__ << ": Cannot create a coarsen schedule for "
                   "U v restriction!\n");

      coarsen_resid_schedules[dest_ln] =
        coarsen_resid_algorithm.
        createSchedule(hierarchy->getPatchLevel(dest_ln),
                       hierarchy->getPatchLevel(dest_ln + 1),
                       &v_coarsen_patch_strategy);
      if (!coarsen_resid_schedules[dest_ln])
        TBOX_ERROR(__FILE__ << ": Cannot create a coarsen schedule for "
                   "R v restriction!\n");
    }

  /// Ordinary ghost fill operator on the coarsest level
  ghostfill_nocoarse_schedules[level_min] =
    ghostfill_nocoarse_algorithm.
    createSchedule(hierarchy->getPatchLevel(level_min));
  if (!ghostfill_nocoarse_schedules[level_min])
    { TBOX_ERROR(__FILE__ << ": Cannot create a refine schedule for v "
                 "ghost filling on bottom level!\n"); }
}
