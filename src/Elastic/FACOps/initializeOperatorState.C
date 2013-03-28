/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for staggered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
/*
************************************************************************
* FACOperatorStrategy virtual initializeOperatorState function.  *
*                                                                      *
* Set internal variables to correspond to the solution passed in.      *
* Look up transfer operators.                                          *
************************************************************************
*/

void Elastic::FACOps::initializeOperatorState
(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& rhs)
{
  deallocateOperatorState();
  int ln;
  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();

  d_hierarchy = solution.getPatchHierarchy();
  d_ln_min = solution.getCoarsestLevelNumber();
  d_ln_max = solution.getFinestLevelNumber();
  d_hopsside = boost::make_shared<SAMRAI::math::HierarchySideDataOpsReal<double> >
    (d_hierarchy, d_ln_min, d_ln_max);

  const int v_id(solution.getComponentDescriptorIndex(0)),
    v_rhs_id(rhs.getComponentDescriptorIndex(0));

#ifdef DEBUG_CHECK_ASSERTIONS

  {
    /*
     * Make sure that solution and rhs data
     *   are of correct type
     *   are allocated
     *   has sufficient ghost width
     */

    if (solution.getNumberOfComponents() != 1)
      TBOX_WARNING(d_object_name
                   << ": Solution vector has multiple components.\n"
                   << "Solver is for component 0 only.\n");
    if (rhs.getNumberOfComponents() != 1)
      TBOX_WARNING(d_object_name
                   << ": RHS vector has multiple components.\n"
                   << "Solver is for component 0 only.\n");

    boost::shared_ptr<SAMRAI::hier::Variable> var;
    {
      vdb->mapIndexToVariable(v_rhs_id,var);
      if (!var) {
        TBOX_ERROR(d_object_name << ": RHS component does not\n"
                   << "correspond to a variable.\n");
      }
      boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > side_var =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
        (var);
      if (!side_var) {
        TBOX_ERROR(d_object_name
                   << ": RHS variable is not side-centered double\n");
      }
    }
    {
      vdb->mapIndexToVariable(v_id,var);
      if (!var) {
        TBOX_ERROR(d_object_name << ": Solution component does not\n"
                   << "correspond to a variable.\n");
      }
      boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > side_var =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
        (var);
      if (!side_var) {
        TBOX_ERROR(d_object_name
                   << ": Solution variable is not side-centered double\n");
      }
    }
    for (ln = d_ln_min; ln <= d_ln_max; ++ln)
      {
        boost::shared_ptr<SAMRAI::hier::PatchLevel> level_ptr =
          d_hierarchy->getPatchLevel(ln);
        SAMRAI::hier::PatchLevel& level = *level_ptr;
        for (SAMRAI::hier::PatchLevel::Iterator pi(level.begin());
             pi!=level.end(); pi++)
          {
            SAMRAI::hier::Patch& patch = **pi;
            boost::shared_ptr<SAMRAI::hier::PatchData> fd=
              patch.getPatchData(v_rhs_id);
        
            if (fd)
              {
                /*
                 * Some data checks can only be done if the data already exists.
                 */
                boost::shared_ptr<SAMRAI::pdat::SideData<double> > cd =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                  (fd);
                if (!cd)
                  TBOX_ERROR(d_object_name
                             << ": RHS data is not side-centered double\n");
                if (cd->getDepth() > 1)
                  TBOX_WARNING(d_object_name
                               << ": RHS data has multiple depths.\n"
                               << "Solver is for depth 0 only.\n");
              }
            boost::shared_ptr<SAMRAI::hier::PatchData> ud =
              patch.getPatchData(v_id);
            if (ud)
              {
                /*
                 * Some data checks can only be done if the data already exists.
                 */
                boost::shared_ptr<SAMRAI::pdat::SideData<double> > cd =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                  (ud);
                if (!cd)
                  TBOX_ERROR(d_object_name
                             << ": Solution data is not side-centered double\n");
                if (cd->getDepth() > 1)
                  TBOX_WARNING(d_object_name
                               << ": Solution data has multiple depths.\n"
                               << "Solver is for depth 0 only.\n");
                if (cd->getGhostCellWidth()
                    < SAMRAI::hier::IntVector::getOne(d_dim))
                  TBOX_ERROR(d_object_name << ": Solution data has "
                             "insufficient ghost width\n");
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
  }
#endif

  /*
   * Initialize the coarse-fine boundary description for the
   * hierarchy.
   */
  d_cf_boundary.resizeArray(d_hierarchy->getNumberOfLevels());

  SAMRAI::hier::IntVector max_gcw(d_dim, 1);
  for (ln = d_ln_min; ln <= d_ln_max; ++ln)
    d_cf_boundary[ln] = boost::make_shared<SAMRAI::hier::CoarseFineBoundary >
      (*d_hierarchy, ln, max_gcw);

  v_coarsen_patch_strategy.coarse_fine=d_cf_boundary;
  /*
   * Get the transfer operators.
   */
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (d_hierarchy->getGridGeometry());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  vdb->mapIndexToVariable(d_side_scratch_id, variable);
  v_prolongation_refine_operator =
    geometry->lookupRefineOperator(variable,v_prolongation_method);

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

  SAMRAI::xfer::RefineAlgorithm v_prolongation_refine_algorithm(d_dim),
    v_ghostfill_refine_algorithm(d_dim), v_nocoarse_refine_algorithm(d_dim);
    
  SAMRAI::xfer::CoarsenAlgorithm v_urestriction_coarsen_algorithm(d_dim,true),
    v_rrestriction_coarsen_algorithm(d_dim,true);

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
        createSchedule(fill_pattern,d_hierarchy->getPatchLevel(dest_ln),
                       boost::shared_ptr<SAMRAI::hier::PatchLevel>(),
                       dest_ln - 1,d_hierarchy);
                       
    if (!v_prolongation_refine_schedules[dest_ln])
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for v prolongation!\n");
    v_ghostfill_refine_schedules[dest_ln] =
      v_ghostfill_refine_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     dest_ln - 1,d_hierarchy,
                     &v_refine_patch_strategy);
    if (!v_ghostfill_refine_schedules[dest_ln])
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling!\n");

    v_nocoarse_refine_schedules[dest_ln] =
      v_nocoarse_refine_algorithm.
      createSchedule(d_hierarchy->getPatchLevel(dest_ln));
    if (!v_nocoarse_refine_schedules[dest_ln])
      TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                 "ghost filling on bottom level!\n");
    }

  /* Coarsening operators */
  for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln)
    {
      v_urestriction_coarsen_schedules[dest_ln] =
        v_urestriction_coarsen_algorithm.
        createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                       d_hierarchy->getPatchLevel(dest_ln + 1),
                       &v_coarsen_patch_strategy);
      if (!v_urestriction_coarsen_schedules[dest_ln])
        TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for "
                   "U v restriction!\n");

      v_rrestriction_coarsen_schedules[dest_ln] =
        v_rrestriction_coarsen_algorithm.
        createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                       d_hierarchy->getPatchLevel(dest_ln + 1),
                       &v_coarsen_patch_strategy);
      if (!v_rrestriction_coarsen_schedules[dest_ln])
        TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for "
                   "R v restriction!\n");
    }

  /* Ordinary ghost fill operator on the coarsest level */
  v_nocoarse_refine_schedules[d_ln_min] =
    v_nocoarse_refine_algorithm.
    createSchedule(d_hierarchy->getPatchLevel(d_ln_min));
  if (!v_nocoarse_refine_schedules[d_ln_min])
    TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for v "
               "ghost filling on bottom level!\n");
}
