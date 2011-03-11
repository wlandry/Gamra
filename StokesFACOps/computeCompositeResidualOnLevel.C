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

#include "Boundary.h"
#include "dRc_dp.h"
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


  // /*
  //  * For the whole level, make sure the internal
  //  * side-centered data is allocated and note
  //  * whether that data should be deallocated when done.
  //  * We do this for the whole level because the data
  //  * undergoes transfer operations which require the
  //  * whole level data.
  //  */

  /* S1. Fill solution ghost data. */

  set_boundaries(v_id,ln,error_equation_indicator);
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

  /*
   * S4. Compute residual on patches in level.
   */

  const hier::Index ip(1,0), jp(0,1);
  for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
    tbox::Pointer<hier::Patch> patch = *pi;
    tbox::Pointer<pdat::CellData<double> >
      p_ptr = solution.getComponentPatchData(0, *patch);
    pdat::CellData<double> &p(*p_ptr);
    tbox::Pointer<pdat::SideData<double> >
      v_ptr = solution.getComponentPatchData(1, *patch);
    pdat::SideData<double> &v(*v_ptr);
    tbox::Pointer<pdat::CellData<double> >
      cell_viscosity_ptr = patch->getPatchData(cell_viscosity_id);
    pdat::CellData<double> &cell_viscosity(*cell_viscosity_ptr);
    tbox::Pointer<pdat::NodeData<double> >
      edge_viscosity_ptr = patch->getPatchData(edge_viscosity_id);
    pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);
    tbox::Pointer<pdat::CellData<double> >
      p_rhs_ptr = rhs.getComponentPatchData(0, *patch);
    pdat::CellData<double> &p_rhs(*p_rhs_ptr);
    tbox::Pointer<pdat::SideData<double> >
      v_rhs_ptr = rhs.getComponentPatchData(1, *patch);
    pdat::SideData<double> &v_rhs(*v_rhs_ptr);
    tbox::Pointer<pdat::CellData<double> >
      p_resid_ptr = residual.getComponentPatchData(0, *patch);
    pdat::CellData<double> &p_resid(*p_resid_ptr);
    tbox::Pointer<pdat::SideData<double> >
      v_resid_ptr = residual.getComponentPatchData(1, *patch);
    pdat::SideData<double> &v_resid(*v_resid_ptr);

    hier::Box pbox=patch->getBox();
    tbox::Pointer<geom::CartesianPatchGeometry> geom = patch->getPatchGeometry();
    double dx = geom->getDx()[0];
    double dy = geom->getDx()[1];

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

            const pdat::SideIndex
              center_x(center,0,pdat::SideIndex::Lower),
              left_x(left,0,pdat::SideIndex::Lower),
              right_x(right,0,pdat::SideIndex::Lower),
              down_x(down,0,pdat::SideIndex::Lower),
              up_x(up,0,pdat::SideIndex::Lower),
              center_y(center,1,pdat::SideIndex::Lower),
              left_y(left,1,pdat::SideIndex::Lower),
              right_y(right,1,pdat::SideIndex::Lower),
              down_y(down,1,pdat::SideIndex::Lower),
              up_y(up,1,pdat::SideIndex::Lower);
            const pdat::NodeIndex
              center_e(center,pdat::NodeIndex::LowerLeft),
              up_e(up,pdat::NodeIndex::LowerLeft),
              right_e(right,pdat::NodeIndex::LowerLeft);

            // tbox::plog << "resid "
            //            << ln << " "
            //            << i << " "
            //            << j << " ";
            /* p */
            if(i!=pbox.upper(0)+1 && j!=pbox.upper(1)+1)
              {
                double dvx_dx=(v(right_x) - v(center_x))/dx;
                double dvy_dy=(v(up_y) - v(center_y))/dy;
                p_resid(center)=p_rhs(center) - dvx_dx - dvy_dy;

                // tbox::plog << "p "
                //            << p_resid(center) << " ";
              }

            /* vx */
            if(j!=pbox.upper(1)+1)
              {
                /* If x==0 */
                if((center[0]==pbox.lower(0) && v(left_x)==boundary_value)
                   || (center[0]==pbox.upper(0)+1
                       && v(right_x)==boundary_value))
                       
                  {
                    v_resid(center_x)=0;
                  }
                else
                  {
                    v_resid(center_x)=v_rhs(center_x)
                      - v_operator(v,p,cell_viscosity,edge_viscosity,center,
                                   left,center_x,right_x,left_x,up_x,down_x,
                                   center_y,up_y,center_e,up_e,ip,dx,dy);
                // tbox::plog << "vx "
                //            << v_resid(center_x) << " "
                //            << v(center_x) << " "
                //            << v(right_x) << " "
                //            << v(left_x) << " ";
                  }
              }

            /* vy */
            if(i!=pbox.upper(0)+1)
              {
                /* If y==0 */
                if((center[1]==pbox.lower(1) && v(down_y)==boundary_value)
                   || (center[1]==pbox.upper(1)+1 && v(up_y)==boundary_value))
                  {
                    v_resid(center_y)=0;
                  }
                else
                  {
                    v_resid(center_y)=v_rhs(center_y)
                      - v_operator(v,p,cell_viscosity,edge_viscosity,center,
                                   down,center_y,up_y,down_y,right_y,left_y,
                                   center_x,right_x,center_e,right_e,jp,
                                   dy,dx);
                  }
                // tbox::plog << "vy "
                //            << v_resid(center_y) << " ";
              }
            // tbox::plog << "\n";
          }
      }
  }

  /* We also need to set the boundaries of the rhs so that coarsening
     works correctly. */
  const int v_rhs_id = rhs.getComponentDescriptorIndex(1);
  set_boundaries(v_rhs_id,ln,true);
  xeqScheduleGhostFillNoCoarse(invalid_id, v_rhs_id, ln);
  

  // if (deallocate_flux_data_when_done) {
  //   level->deallocatePatchData(flux_id);
  // }

  t_compute_composite_residual->stop();
}

