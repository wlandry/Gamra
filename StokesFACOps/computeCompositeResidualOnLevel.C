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
********************************************************************
* FACOperatorStrategy virtual                                *
* computeCompositeResidualOnLevel function                         *
********************************************************************
*/

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
  // d_bc_helper.setTargetDataId(soln_id);
  // d_bc_helper.setHomogeneousBc(error_equation_indicator);

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


  // /*
  //  * For the whole level, make sure the internal
  //  * side-centered data is allocated and note
  //  * whether that data should be deallocated when done.
  //  * We do this for the whole level because the data
  //  * undergoes transfer operations which require the
  //  * whole level data.
  //  */
  // bool deallocate_flux_data_when_done = false;
  // if (flux_id == d_flux_scratch_id) {
  //   if (!level->checkAllocated(flux_id)) {
  //     level->allocatePatchData(flux_id);
  //     deallocate_flux_data_when_done = true;
  //   }

  /* S1. Fill solution ghost data. */
  if (ln > d_ln_min) {
    /* Fill from current, next coarser level and physical boundary */
    xeqScheduleGhostFill(p_id, v_id, ln);
  } else {
    /* Fill from current and physical boundary */
    xeqScheduleGhostFillNoCoarse(p_id, v_id, ln);
  }

  // /*
  //  * S2. Compute flux on patches in level.
  //  */
  // for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
  //   tbox::Pointer<hier::Patch> patch = *pi;

  //   tbox::Pointer<pdat::CellData<double> >
  //     soln_data = solution.getComponentPatchData(0, *patch);
  //   tbox::Pointer<pdat::CellData<double> >
  //     rhs_data = rhs.getComponentPatchData(0, *patch);
  //   tbox::Pointer<pdat::CellData<double> >
  //     residual_data = residual.getComponentPatchData(0, *patch);
  //   tbox::Pointer<pdat::SideData<double> >
  //     flux_data = patch->getPatchData(flux_id);
  //   computeFluxOnPatch(
  //                      *patch,
  //                      level->getRatioToCoarserLevel(),
  //                      *soln_data,
  //                      *flux_data);

  // }

  // /*
  //  * S3. Coarsen oflux data from next finer level so that
  //  * the computed flux becomes the composite grid flux.
  //  */
  // if (ln < d_ln_max) {
  //   xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
  // }

  /*
   * S4. Compute residual on patches in level.
   */
  double viscosity=1;

  for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
    tbox::Pointer<hier::Patch> patch = *pi;
    tbox::Pointer<pdat::CellData<double> >
      p = solution.getComponentPatchData(0, *patch);
    tbox::Pointer<pdat::SideData<double> >
      v = solution.getComponentPatchData(1, *patch);
    tbox::Pointer<pdat::CellData<double> >
      p_rhs = rhs.getComponentPatchData(0, *patch);
    tbox::Pointer<pdat::SideData<double> >
      v_rhs = rhs.getComponentPatchData(1, *patch);
    tbox::Pointer<pdat::CellData<double> >
      p_resid = residual.getComponentPatchData(0, *patch);
    tbox::Pointer<pdat::SideData<double> >
      v_resid = residual.getComponentPatchData(1, *patch);

    hier::Box pbox=patch->getBox();
    tbox::Pointer<geom::CartesianPatchGeometry> geom = patch->getPatchGeometry();
    double dx = *(geom->getDx());
    double dy = *(geom->getDx());

    for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
      {
        for(int i=pbox.lower(0); i<=pbox.upper(0)+1; ++i)
          {
            pdat::CellIndex center(tbox::Dimension(2));
            center[0]=i;
            center[1]=j;

            pdat::CellIndex up(center), down(center), right(center),
              left(center);

            ++up[1];
            --down[1];
            ++right[0];
            --left[0];

            /* p */
            if(i!=pbox.upper(0)+1 && j!=pbox.upper(1)+1)
              {
                double dvx_dx=
                  ((*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                        pdat::SideIndex::Upper))
                   - (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                          pdat::SideIndex::Lower)))/dx;
                double dvy_dy=
                  ((*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                        pdat::SideIndex::Upper))
                   - (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                          pdat::SideIndex::Lower)))/dy;
                (*p_resid)(center)=(*p_rhs)(center) - dvx_dx - dvy_dy;
              }

            /* vx */
            if(j!=pbox.upper(1)+1)
              {
                /* If x==0 */
                if((center[0]==pbox.lower(0)
                    && geom->getTouchesRegularBoundary(0,0))
                   || (center[0]==pbox.upper(0)+1
                       && geom->getTouchesRegularBoundary(0,1)))
                  {
                    (*v_resid)(pdat::SideIndex(center,pdat::SideIndex::X,
                                               pdat::SideIndex::Lower))=0;
                  }
                else
                  {
                    double dp_dx, d2vx_dxx, d2vx_dyy, C_vx;
                    /* If y==0 */
                    if(center[1]==pbox.lower(1)
                       && geom->getTouchesRegularBoundary(1,0))
                      {
                        d2vx_dyy=
                          ((*v)(pdat::SideIndex(up,pdat::SideIndex::X,
                                                pdat::SideIndex::Lower))
                           - (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                  pdat::SideIndex::Lower)))
                          /(dy*dy);
                        C_vx=-viscosity*(2/(dx*dx) + 1/(dy*dy));
                      }
                    /* If y==max_y */
                    else if(center[1]==pbox.upper(1)
                            && geom->getTouchesRegularBoundary(1,1))
                      {
                        d2vx_dyy=
                          (-(*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                 pdat::SideIndex::Lower))
                           + (*v)(pdat::SideIndex(down,pdat::SideIndex::X,
                                                  pdat::SideIndex::Lower)))
                          /(dy*dy);
                        C_vx=-viscosity*(2/(dx*dx) + 1/(dy*dy));
                        // tbox::plog << "vx y1 boundary ";
                      }
                    else
                      {
                        d2vx_dyy=
                          ((*v)(pdat::SideIndex(up,pdat::SideIndex::X,
                                                pdat::SideIndex::Lower))
                           - 2*(*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                    pdat::SideIndex::Lower))
                           + (*v)(pdat::SideIndex(down,pdat::SideIndex::X,
                                                  pdat::SideIndex::Lower)))
                          /(dy*dy);

                        C_vx=-2*viscosity*(1/(dx*dx) + 1/(dy*dy));
                      }
                    d2vx_dxx=((*v)(pdat::SideIndex(left,pdat::SideIndex::X,
                                                   pdat::SideIndex::Lower))
                              - 2*(*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                       pdat::SideIndex::Lower))
                              + (*v)(pdat::SideIndex(right,pdat::SideIndex::X,
                                                     pdat::SideIndex::Lower)))
                      /(dx*dx);

                    dp_dx=((*p)(center)-(*p)(left))/dx;
                              
                    (*v_resid)(pdat::SideIndex(center,pdat::SideIndex::X,
                                               pdat::SideIndex::Lower))=
                      (*v_rhs)(pdat::SideIndex(center,pdat::SideIndex::X,
                                               pdat::SideIndex::Lower))
                      - viscosity*(d2vx_dxx + d2vx_dyy) + dp_dx;
                  }
              }

            /* vy */
            if(i!=pbox.upper(0)+1)
              {
                /* If y==0 */
                if((center[1]==pbox.lower(1)
                    && geom->getTouchesRegularBoundary(1,0))
                   || (center[1]==pbox.upper(1)+1
                       && geom->getTouchesRegularBoundary(1,1)))
                  {
                    (*v_resid)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                               pdat::SideIndex::Lower))=0;
                  }
                else
                  {
                    double dp_dy, d2vy_dxx, d2vy_dyy, C_vy;
                    /* If x==0 */
                    if(center[0]==pbox.lower(0)
                       && geom->getTouchesRegularBoundary(0,0))
                      {
                        d2vy_dxx=
                          ((*v)(pdat::SideIndex(right,pdat::SideIndex::Y,
                                                pdat::SideIndex::Lower))
                           - (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                  pdat::SideIndex::Lower)))
                          /(dx*dx);
                        C_vy=-viscosity*(1/(dx*dx) + 2/(dy*dy));
                      }
                    /* If x==max_x */
                    else if(center[0]==pbox.upper(0)
                            && geom->getTouchesRegularBoundary(0,1))
                      {
                        d2vy_dxx=
                          ((*v)(pdat::SideIndex(left,pdat::SideIndex::Y,
                                                pdat::SideIndex::Lower))
                           - (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                  pdat::SideIndex::Lower)))
                          /(dx*dx);
                        C_vy=-viscosity*(1/(dx*dx) + 2/(dy*dy));
                      }
                    else
                      {
                        d2vy_dxx=
                          ((*v)(pdat::SideIndex(left,pdat::SideIndex::Y,
                                                pdat::SideIndex::Lower))
                           - 2*(*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                    pdat::SideIndex::Lower))
                           + (*v)(pdat::SideIndex(right,pdat::SideIndex::Y,
                                                  pdat::SideIndex::Lower)))
                          /(dx*dx);

                        C_vy=-2*viscosity*(1/(dx*dx) + 1/(dy*dy));
                      }
                    d2vy_dyy=((*v)(pdat::SideIndex(up,pdat::SideIndex::Y,
                                                   pdat::SideIndex::Lower))
                              - 2*(*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                       pdat::SideIndex::Lower))
                              + (*v)(pdat::SideIndex(down,pdat::SideIndex::Y,
                                                     pdat::SideIndex::Lower)))
                      /(dy*dy);

                    dp_dy=((*p)(center)-(*p)(down))/dy;
                              
                    (*v_resid)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                               pdat::SideIndex::Lower))=
                      (*v_rhs)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                               pdat::SideIndex::Lower))
                      - viscosity*(d2vy_dxx + d2vy_dyy) + dp_dy;
                  }
              }
          }
      }

    // if (ln > d_ln_min) {
    //   /*
    //    * Save outerflux data so that next coarser level
    //    *  can compute its coarse-fine composite flux.
    //    *  This is not strictly needed in this "compute residual"
    //    *  loop through the patches, but we put it here to
    //    *  avoid writing another loop for it.
    //    */
    //   tbox::Pointer<pdat::OutersideData<double> >
    //     oflux_data = patch->getPatchData(d_oflux_scratch_id);

    //   TBOX_ASSERT(oflux_data);

    //   oflux_data->copy(*flux_data);
    // }
  }

  // if (deallocate_flux_data_when_done) {
  //   level->deallocatePatchData(flux_id);
  // }

  t_compute_composite_residual->stop();
}

