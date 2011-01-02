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

namespace SAMRAI {
  namespace solv {

    /*
********************************************************************
* FACOperatorStrategy virtual                                *
* computeCompositeResidualOnLevel function                         *
********************************************************************
*/

    void StokesFACOps::computeCompositeResidualOnLevel(
                                                       SAMRAIVectorReal<double>& residual,
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
      const int soln_id = solution.getComponentDescriptorIndex(0);
      d_bc_helper.setTargetDataId(soln_id);
      d_bc_helper.setHomogeneousBc(error_equation_indicator);

      const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

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
      {
        tbox::Pointer<xfer::RefineSchedule> ln_refine_schedule;
        if (ln > d_ln_min) {
          /* Fill from current, next coarser level and physical boundary */
          xeqScheduleGhostFill(soln_id, ln);
        } else {
          /* Fill from current and physical boundary */
          xeqScheduleGhostFillNoCoarse(soln_id, ln);
        }
      }

      /*
       * For the whole level, make sure the internal
       * side-centered data is allocated and note
       * whether that data should be deallocated when done.
       * We do this for the whole level because the data
       * undergoes transfer operations which require the
       * whole level data.
       */
      bool deallocate_flux_data_when_done = false;
      if (flux_id == d_flux_scratch_id) {
        if (!level->checkAllocated(flux_id)) {
          level->allocatePatchData(flux_id);
          deallocate_flux_data_when_done = true;
        }
      }

      /*
       * S2. Compute flux on patches in level.
       */
      for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
        tbox::Pointer<hier::Patch> patch = *pi;

        tbox::Pointer<pdat::CellData<double> >
          soln_data = solution.getComponentPatchData(0, *patch);
        tbox::Pointer<pdat::CellData<double> >
          rhs_data = rhs.getComponentPatchData(0, *patch);
        tbox::Pointer<pdat::CellData<double> >
          residual_data = residual.getComponentPatchData(0, *patch);
        tbox::Pointer<pdat::SideData<double> >
          flux_data = patch->getPatchData(flux_id);
        computeFluxOnPatch(
                           *patch,
                           level->getRatioToCoarserLevel(),
                           *soln_data,
                           *flux_data);

      }

      /*
       * S3. Coarsen oflux data from next finer level so that
       * the computed flux becomes the composite grid flux.
       */
      if (ln < d_ln_max) {
        xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
      }

      /*
       * S4. Compute residual on patches in level.
       */
      for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
        tbox::Pointer<hier::Patch> patch = *pi;
        tbox::Pointer<pdat::CellData<double> >
          soln_data = solution.getComponentPatchData(0, *patch);
        tbox::Pointer<pdat::CellData<double> >
          rhs_data = rhs.getComponentPatchData(0, *patch);
        tbox::Pointer<pdat::CellData<double> >
          residual_data = residual.getComponentPatchData(0, *patch);
        tbox::Pointer<pdat::SideData<double> >
          flux_data = patch->getPatchData(flux_id);
        computeResidualOnPatch(*patch,
                               *flux_data,
                               *soln_data,
                               *rhs_data,
                               *residual_data);

        if (ln > d_ln_min) {
          /*
           * Save outerflux data so that next coarser level
           *  can compute its coarse-fine composite flux.
           *  This is not strictly needed in this "compute residual"
           *  loop through the patches, but we put it here to
           *  avoid writing another loop for it.
           */
          tbox::Pointer<pdat::OutersideData<double> >
            oflux_data = patch->getPatchData(d_oflux_scratch_id);

          TBOX_ASSERT(oflux_data);

          oflux_data->copy(*flux_data);
        }
      }

      if (deallocate_flux_data_when_done) {
        level->deallocatePatchData(flux_id);
      }

      t_compute_composite_residual->stop();
    }

  }
}

