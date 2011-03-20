#include "StokesFACOps.h"

void SAMRAI::solv::StokesFACOps::computeCompositeResidualOnLevel
(SAMRAIVectorReal<double>& residual,
 const SAMRAIVectorReal<double>& solution,
 const SAMRAIVectorReal<double>& rhs,
 int ln,
 bool error_equation_indicator) {

  t_compute_composite_residual->start();

  checkInputPatchDataIndices();
#ifdef DEBUG_CHECK_ASSERTIONS
  if (residual.getPatchHierarchy() != d_hierarchy
      || solution.getPatchHierarchy() != d_hierarchy
      || rhs.getPatchHierarchy() != d_hierarchy) {
    TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
               "internal hierarchy.");
  }
#endif
  tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

  /*
   * Set up the bc helper so that when we use a refine schedule
   * to fill ghosts, the correct data is operated on.
   */
  const int p_id = solution.getComponentDescriptorIndex(0);
  const int v_id = solution.getComponentDescriptorIndex(1);
  p_refine_patch_strategy.setTargetDataId(p_id);
  v_refine_patch_strategy.setTargetDataId(v_id);
  // v_refine_patch_strategy.setHomogeneousBc(error_equation_indicator);

  /*
   * Assumptions:
   * 1. Data does not yet exist in ghost boundaries.
   * 2. Residual data on next finer grid (if any)
   *    has been computed already.
   * 3. Flux data from next finer grid (if any) has
   *    been computed but has not been coarsened to
   *    this level.
   *
   * Steps:
   * S1. Fill solution ghost data by refinement
   *     or setting physical boundary conditions.
   *     This also brings in information from coarser
   *     to form the composite grid flux.
   * S2. Compute flux on ln.
   * S3. If next finer is available,
   *     Coarsen flux data on next finer level,
   *     overwriting flux computed from coarse data.
   * S4. Compute residual data from flux.
   */

  /* S1. Fill solution ghost data. */

  set_boundaries(v_id,ln,error_equation_indicator);
  if (ln > d_ln_min) {
    /* Fill from current, next coarser level and physical boundary */
    xeqScheduleGhostFill(p_id, v_id, ln);
  } else {
    /* Fill from current and physical boundary */
    xeqScheduleGhostFillNoCoarse(p_id, v_id, ln);
  }

  /*
   * S4. Compute residual on patches in level.
   */

  for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
    tbox::Pointer<hier::Patch> patch = *pi;
    tbox::Pointer<pdat::CellData<double> >
      p_ptr = solution.getComponentPatchData(0, *patch);
    tbox::Pointer<pdat::SideData<double> >
      v_ptr = solution.getComponentPatchData(1, *patch);
    tbox::Pointer<pdat::CellData<double> >
      cell_viscosity_ptr = patch->getPatchData(cell_viscosity_id);
    tbox::Pointer<pdat::CellData<double> >
      p_rhs_ptr = rhs.getComponentPatchData(0, *patch);
    tbox::Pointer<pdat::SideData<double> >
      v_rhs_ptr = rhs.getComponentPatchData(1, *patch);
    tbox::Pointer<pdat::CellData<double> >
      p_resid_ptr = residual.getComponentPatchData(0, *patch);
    tbox::Pointer<pdat::SideData<double> >
      v_resid_ptr = residual.getComponentPatchData(1, *patch);

    hier::Box pbox=patch->getBox();
    pbox.growUpper(hier::IntVector::getOne(d_dim));
    tbox::Pointer<geom::CartesianPatchGeometry>
      geom = patch->getPatchGeometry();

    switch(d_dim.getValue())
      {
      case 2:
        residual_2D(*p_ptr,*v_ptr,*cell_viscosity_ptr,*p_rhs_ptr,*v_rhs_ptr,
                    *p_resid_ptr,*v_resid_ptr,*patch,pbox,*geom);
        break;
      case 3:
        residual_3D(*p_ptr,*v_ptr,*cell_viscosity_ptr,*p_rhs_ptr,*v_rhs_ptr,
                    *p_resid_ptr,*v_resid_ptr,*patch,pbox,*geom);
        break;
      default:
        abort();
      }
  }

  t_compute_composite_residual->stop();
}

