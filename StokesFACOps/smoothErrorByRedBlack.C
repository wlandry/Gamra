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
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void SAMRAI::solv::StokesFACOps::smoothErrorByRedBlack(SAMRAIVectorReal<double>& data,
                                                       const SAMRAIVectorReal<double>&
                                                       residual,
                                                       int ln,
                                                       int num_sweeps,
                                                       double residual_tolerance)
{

  checkInputPatchDataIndices();

#ifdef DEBUG_CHECK_ASSERTIONS
  if (data.getPatchHierarchy() != d_hierarchy
      || residual.getPatchHierarchy() != d_hierarchy) {
    TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
               "internal hierarchy.");
  }
#endif
  tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

  const int data_id = data.getComponentDescriptorIndex(0);

  const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

  d_bc_helper.setTargetDataId(data_id);
  d_bc_helper.setHomogeneousBc(true);
  xeqScheduleGhostFillNoCoarse(data_id, ln);

  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(data_id, ln);
  }

  /*
   * Smooth the number of sweeps specified or until
   * the convergence is satisfactory.
   */
  int isweep;
  double red_maxres, blk_maxres, maxres = 0;
  red_maxres = blk_maxres = residual_tolerance + 1;
  /*
   * Instead of checking residual convergence globally,
   * we check the not_converged flag.  This avoids possible
   * round-off errors affecting different processes differently,
   * leading to disagreement on whether to continue smoothing.
   */
  int not_converged = 1;
  for (isweep = 0; isweep < num_sweeps && not_converged; ++isweep) {
    red_maxres = blk_maxres = 0;

    // Red sweep.
    xeqScheduleGhostFillNoCoarse(data_id, ln);
    for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
      tbox::Pointer<hier::Patch> patch = *pi;

      bool deallocate_flux_data_when_done = false;
      if (flux_id == d_flux_scratch_id) {
        /*
         * Using internal temporary storage for flux.
         * For each patch, make sure the internal
         * side-centered data is allocated and note
         * whether that data should be deallocated when done.
         */
        if (!patch->checkAllocated(flux_id)) {
          patch->allocatePatchData(flux_id);
          deallocate_flux_data_when_done = true;
        }
      }

      tbox::Pointer<pdat::CellData<double> >
        scalar_field_data = d_stokes_spec.cIsVariable() ?
        patch->getPatchData(d_stokes_spec.getCPatchDataId()) :
        tbox::Pointer<hier::PatchData>(NULL);
      tbox::Pointer<pdat::CellData<double> >
        err_data = data.getComponentPatchData(0, *patch);
      tbox::Pointer<pdat::CellData<double> >
        residual_data = residual.getComponentPatchData(0, *patch);
      tbox::Pointer<pdat::SideData<double> >
        flux_data = patch->getPatchData(flux_id);

      computeFluxOnPatch(
                         *patch,
                         level->getRatioToCoarserLevel(),
                         *err_data,
                         *flux_data);

      redOrBlackSmoothingOnPatch(*patch,
                                 *flux_data,
                                 *residual_data,
                                 *err_data,
                                 'r',
                                 &red_maxres);

      if (deallocate_flux_data_when_done) {
        patch->deallocatePatchData(flux_id);
      }
    }        // End patch number *pi
    xeqScheduleGhostFillNoCoarse(data_id, ln);

    // Black sweep.
    for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
      tbox::Pointer<hier::Patch> patch = *pi;

      bool deallocate_flux_data_when_done = false;
      if (flux_id == d_flux_scratch_id) {
        /*
         * Using internal temporary storage for flux.
         * For each patch, make sure the internal
         * side-centered data is allocated and note
         * whether that data should be deallocated when done.
         */
        if (!patch->checkAllocated(flux_id)) {
          patch->allocatePatchData(flux_id);
          deallocate_flux_data_when_done = true;
        }
      }

      tbox::Pointer<pdat::CellData<double> >
        scalar_field_data = d_stokes_spec.cIsVariable() ?
        patch->getPatchData(d_stokes_spec.getCPatchDataId()) :
        tbox::Pointer<hier::PatchData>(NULL);
      tbox::Pointer<pdat::CellData<double> >
        err_data = data.getComponentPatchData(0, *patch);
      tbox::Pointer<pdat::CellData<double> >
        residual_data = residual.getComponentPatchData(0, *patch);
      tbox::Pointer<pdat::SideData<double> >
        flux_data = patch->getPatchData(flux_id);

      computeFluxOnPatch(
                         *patch,
                         level->getRatioToCoarserLevel(),
                         *err_data,
                         *flux_data);

      redOrBlackSmoothingOnPatch(*patch,
                                 *flux_data,
                                 *residual_data,
                                 *err_data,
                                 'b',
                                 &blk_maxres);

      if (deallocate_flux_data_when_done) {
        patch->deallocatePatchData(flux_id);
      }
    }        // End patch number *pi
    xeqScheduleGhostFillNoCoarse(data_id, ln);
    if (residual_tolerance >= 0.0) {
      /*
       * Check for early end of sweeps due to convergence
       * only if it is numerically possible (user gave a
       * non negative value for residual tolerance).
       */
      maxres = tbox::MathUtilities<double>::Max(red_maxres, blk_maxres);
      not_converged = maxres > residual_tolerance;
      const tbox::SAMRAI_MPI& mpi(d_hierarchy->getDomainMappedBoxLevel().getMPI());
      if (mpi.getSize() > 1) {
        mpi.AllReduce(&not_converged, 1, MPI_MAX);
      }
    }
  }        // End sweep number isweep
  if (d_enable_logging) tbox::plog
                          << d_object_name << " RBGS smoothing maxres = " << maxres << "\n"
                          << "  after " << isweep << " sweeps.\n";

}

