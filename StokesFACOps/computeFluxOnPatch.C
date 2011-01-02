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

void F77_FUNC(compfluxvardc2d, COMPFLUXVARDC2D) (
   double* xflux,
   double* yflux,
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
   const double* dx);
void F77_FUNC(compfluxcondc2d, COMPFLUXCONDC2D) (
   double* xflux,
   double* yflux,
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
   const double* dx);

void F77_FUNC(compfluxvardc3d, COMPFLUXVARDC3D) (
   double* xflux,
   double* yflux,
   double* zflux,
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
   const double* dx);
void F77_FUNC(compfluxcondc3d, COMPFLUXCONDC3D) (
   double* xflux,
   double* yflux,
   double* zflux,
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
   const double* dx);

}
    /*
*******************************************************************
*                                                                 *
* AMR-unaware patch-centered computational kernels.               *
*                                                                 *
*******************************************************************
*/

    void StokesFACOps::computeFluxOnPatch(
                                          const hier::Patch& patch,
                                          const hier::IntVector& ratio_to_coarser_level,
                                          const pdat::CellData<double>& w_data,
                                          pdat::SideData<double>& Dgradw_data) const
    {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(d_dim, patch, ratio_to_coarser_level, w_data,
                                      Dgradw_data);
      TBOX_ASSERT(patch.inHierarchy());
      TBOX_ASSERT(w_data.getGhostCellWidth() >=
                  hier::IntVector::getOne(ratio_to_coarser_level.getDim()));

      tbox::Pointer<geom::CartesianPatchGeometry> patch_geom =
        patch.getPatchGeometry();
      const hier::Box& box = patch.getBox();
      const int* lower = &box.lower()[0];
      const int* upper = &box.upper()[0];
      const double* dx = patch_geom->getDx();

      double D_value;
      tbox::Pointer<pdat::SideData<double> > D_data;
      if (d_stokes_spec.dIsConstant()) {
        D_value = d_stokes_spec.getDConstant();
      } else {
        D_data = patch.getPatchData(d_stokes_spec.getDPatchDataId());
      }

      if (d_stokes_spec.dIsConstant()) {
        if (d_dim == tbox::Dimension(2)) {
          F77_FUNC(compfluxcondc2d, COMPFLUXCONDC2D) (
                                                      Dgradw_data.getPointer(0),
                                                      Dgradw_data.getPointer(1),
                                                      &Dgradw_data.getGhostCellWidth()[0],
                                                      &Dgradw_data.getGhostCellWidth()[1],
                                                      D_value,
                                                      w_data.getPointer(),
                                                      &w_data.getGhostCellWidth()[0],
                                                      &w_data.getGhostCellWidth()[1],
                                                      &lower[0], &upper[0],
                                                      &lower[1], &upper[1],
                                                      dx);
        } else if (d_dim == tbox::Dimension(3)) {
          F77_FUNC(compfluxcondc3d, COMPFLUXCONDC3D) (
                                                      Dgradw_data.getPointer(0),
                                                      Dgradw_data.getPointer(1),
                                                      Dgradw_data.getPointer(2),
                                                      &Dgradw_data.getGhostCellWidth()[0],
                                                      &Dgradw_data.getGhostCellWidth()[1],
                                                      &Dgradw_data.getGhostCellWidth()[2],
                                                      D_value,
                                                      w_data.getPointer(),
                                                      &w_data.getGhostCellWidth()[0],
                                                      &w_data.getGhostCellWidth()[1],
                                                      &w_data.getGhostCellWidth()[2],
                                                      &lower[0], &upper[0],
                                                      &lower[1], &upper[1],
                                                      &lower[2], &upper[2],
                                                      dx);
        }
      } else {
        if (d_dim == tbox::Dimension(2)) {
          F77_FUNC(compfluxvardc2d, COMPFLUXVARDC2D) (
                                                      Dgradw_data.getPointer(0),
                                                      Dgradw_data.getPointer(1),
                                                      &Dgradw_data.getGhostCellWidth()[0],
                                                      &Dgradw_data.getGhostCellWidth()[1],
                                                      D_data->getPointer(0),
                                                      D_data->getPointer(1),
                                                      &D_data->getGhostCellWidth()[0],
                                                      &D_data->getGhostCellWidth()[1],
                                                      w_data.getPointer(),
                                                      &w_data.getGhostCellWidth()[0],
                                                      &w_data.getGhostCellWidth()[1],
                                                      &lower[0], &upper[0],
                                                      &lower[1], &upper[1],
                                                      dx);
        }
        if (d_dim == tbox::Dimension(3)) {
          F77_FUNC(compfluxvardc3d, COMPFLUXVARDC3D) (
                                                      Dgradw_data.getPointer(0),
                                                      Dgradw_data.getPointer(1),
                                                      Dgradw_data.getPointer(2),
                                                      &Dgradw_data.getGhostCellWidth()[0],
                                                      &Dgradw_data.getGhostCellWidth()[1],
                                                      &Dgradw_data.getGhostCellWidth()[2],
                                                      D_data->getPointer(0),
                                                      D_data->getPointer(1),
                                                      D_data->getPointer(2),
                                                      &D_data->getGhostCellWidth()[0],
                                                      &D_data->getGhostCellWidth()[1],
                                                      &D_data->getGhostCellWidth()[2],
                                                      w_data.getPointer(),
                                                      &w_data.getGhostCellWidth()[0],
                                                      &w_data.getGhostCellWidth()[1],
                                                      &w_data.getGhostCellWidth()[2],
                                                      &lower[0], &upper[0],
                                                      &lower[1], &upper[1],
                                                      &lower[2], &upper[2],
                                                      dx);
        }
      }

      const int patch_ln = patch.getPatchLevelNumber();

      if (d_cf_discretization == "Ewing" && patch_ln > d_ln_min) {
        ewingFixFlux(patch,
                     w_data,
                     Dgradw_data,
                     ratio_to_coarser_level);
      }

    }

  }
}
