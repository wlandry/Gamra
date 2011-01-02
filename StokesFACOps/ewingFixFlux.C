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

extern "C" {
void F77_FUNC(ewingfixfluxvardc2d, EWINGFIXFLUXVARDC2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const int* dcgi,
   const int* dcgj,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);
void F77_FUNC(ewingfixfluxcondc2d, EWINGFIXFLUXCONDC2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double & diff_coef,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);

void F77_FUNC(ewingfixfluxvardc3d, EWINGFIXFLUXVARDC3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const double* zdiff_coef,
   const int* dcgi,
   const int* dcgj,
   const int* dcgk,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);
void F77_FUNC(ewingfixfluxcondc3d, EWINGFIXFLUXCONDC3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double & diff_coef,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);
}
    /*
********************************************************************
* Fix flux on coarse-fine boundaries computed from a               *
* constant-refine interpolation of coarse level data.              *
********************************************************************
*/

    void StokesFACOps::ewingFixFlux(
                                    const hier::Patch& patch,
                                    const pdat::CellData<double>& soln_data,
                                    pdat::SideData<double>& flux_data,
                                    const hier::IntVector& ratio_to_coarser) const
    {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(d_dim, patch, soln_data, flux_data,
                                      ratio_to_coarser);

      const int patch_ln = patch.getPatchLevelNumber();
      const hier::GlobalId id = patch.getGlobalId();
      tbox::Pointer<geom::CartesianPatchGeometry> patch_geom =
        patch.getPatchGeometry();
      const double* dx = patch_geom->getDx();
      const hier::Box& patch_box(patch.getBox());
      const hier::Index& plower = patch_box.lower();
      const hier::Index& pupper = patch_box.upper();

      const tbox::Array<hier::BoundaryBox>& bboxes =
        d_cf_boundary[patch_ln]->getBoundaries(id, 1);
      int bn, nboxes = bboxes.getSize();

      if (d_stokes_spec.dIsVariable()) {

        tbox::Pointer<pdat::SideData<double> > diffcoef_data;
        diffcoef_data = patch.getPatchData(d_stokes_spec.getDPatchDataId());

        for (bn = 0; bn < nboxes; ++bn) {
          const hier::BoundaryBox& boundary_box = bboxes[bn];

          TBOX_ASSERT(boundary_box.getBoundaryType() == 1);

          const hier::Box& bdry_box = boundary_box.getBox();
          const hier::Index& blower = bdry_box.lower();
          const hier::Index& bupper = bdry_box.upper();
          const int location_index = boundary_box.getLocationIndex();
          if (d_dim == tbox::Dimension(2)) {
            F77_FUNC(ewingfixfluxvardc2d, EWINGFIXFLUXVARDC2D) (
                                                                flux_data.getPointer(0), flux_data.getPointer(1),
                                                                &flux_data.getGhostCellWidth()[0],
                                                                &flux_data.getGhostCellWidth()[1],
                                                                diffcoef_data->getPointer(0), diffcoef_data->getPointer(1),
                                                                &diffcoef_data->getGhostCellWidth()[0],
                                                                &diffcoef_data->getGhostCellWidth()[1],
                                                                soln_data.getPointer(),
                                                                &soln_data.getGhostCellWidth()[0],
                                                                &soln_data.getGhostCellWidth()[1],
                                                                &plower[0], &pupper[0], &plower[1], &pupper[1],
                                                                &location_index,
                                                                &ratio_to_coarser[0],
                                                                &blower[0], &bupper[0],
                                                                dx);
          } else if (d_dim == tbox::Dimension(3)) {
            F77_FUNC(ewingfixfluxvardc3d, EWINGFIXFLUXVARDC3D) (
                                                                flux_data.getPointer(0),
                                                                flux_data.getPointer(1),
                                                                flux_data.getPointer(2),
                                                                &flux_data.getGhostCellWidth()[0],
                                                                &flux_data.getGhostCellWidth()[1],
                                                                &flux_data.getGhostCellWidth()[2],
                                                                diffcoef_data->getPointer(0),
                                                                diffcoef_data->getPointer(1),
                                                                diffcoef_data->getPointer(2),
                                                                &diffcoef_data->getGhostCellWidth()[0],
                                                                &diffcoef_data->getGhostCellWidth()[1],
                                                                &diffcoef_data->getGhostCellWidth()[2],
                                                                soln_data.getPointer(),
                                                                &soln_data.getGhostCellWidth()[0],
                                                                &soln_data.getGhostCellWidth()[1],
                                                                &soln_data.getGhostCellWidth()[2],
                                                                &plower[0], &pupper[0],
                                                                &plower[1], &pupper[1],
                                                                &plower[2], &pupper[2],
                                                                &location_index,
                                                                &ratio_to_coarser[0],
                                                                &blower[0], &bupper[0],
                                                                dx);
          } else {
            TBOX_ERROR("StokesFACOps : DIM > 3 not supported" << std::endl);
          }

        }
      } else {

        const double diffcoef_constant = d_stokes_spec.getDConstant();

        for (bn = 0; bn < nboxes; ++bn) {
          const hier::BoundaryBox& boundary_box = bboxes[bn];

          TBOX_ASSERT(boundary_box.getBoundaryType() == 1);

          const hier::Box& bdry_box = boundary_box.getBox();
          const hier::Index& blower = bdry_box.lower();
          const hier::Index& bupper = bdry_box.upper();
          const int location_index = boundary_box.getLocationIndex();
          if (d_dim == tbox::Dimension(2)) {
            F77_FUNC(ewingfixfluxcondc2d, EWINGFIXFLUXCONDC2D) (
                                                                flux_data.getPointer(0), flux_data.getPointer(1),
                                                                &flux_data.getGhostCellWidth()[0],
                                                                &flux_data.getGhostCellWidth()[1],
                                                                diffcoef_constant,
                                                                soln_data.getPointer(),
                                                                &soln_data.getGhostCellWidth()[0],
                                                                &soln_data.getGhostCellWidth()[1],
                                                                &plower[0], &pupper[0],
                                                                &plower[1], &pupper[1],
                                                                &location_index,
                                                                &ratio_to_coarser[0],
                                                                &blower[0], &bupper[0],
                                                                dx);
          } else if (d_dim == tbox::Dimension(3)) {
            F77_FUNC(ewingfixfluxcondc3d, EWINGFIXFLUXCONDC3D) (
                                                                flux_data.getPointer(0),
                                                                flux_data.getPointer(1),
                                                                flux_data.getPointer(2),
                                                                &flux_data.getGhostCellWidth()[0],
                                                                &flux_data.getGhostCellWidth()[1],
                                                                &flux_data.getGhostCellWidth()[2],
                                                                diffcoef_constant,
                                                                soln_data.getPointer(),
                                                                &soln_data.getGhostCellWidth()[0],
                                                                &soln_data.getGhostCellWidth()[1],
                                                                &soln_data.getGhostCellWidth()[2],
                                                                &plower[0], &pupper[0],
                                                                &plower[1], &pupper[1],
                                                                &plower[2], &pupper[2],
                                                                &location_index,
                                                                &ratio_to_coarser[0],
                                                                &blower[0], &bupper[0],
                                                                dx);
          }
        }
      }
    }

  }
}
