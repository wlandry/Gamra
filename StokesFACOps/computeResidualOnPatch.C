/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#ifndef included_solv_StokesFACOps_C
#define included_solv_StokesFACOps_C

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
void F77_FUNC(compresvarsca2d, COMPRESVARSCA2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx);
void F77_FUNC(compresconsca2d, COMPRESCONSCA2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const double & scalar_field,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx);
void F77_FUNC(compresvarsca3d, COMPRESVARSCA3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const int* residualgk,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   const int* scalar_field_gk,
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
void F77_FUNC(compresconsca3d, COMPRESCONSCA3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const int* residualgk,
   const double & scalar_field,
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

    void StokesFACOps::computeResidualOnPatch(
                                              const hier::Patch& patch,
                                              const pdat::SideData<double>& flux_data,
                                              const pdat::CellData<double>& soln_data,
                                              const pdat::CellData<double>& rhs_data,
                                              pdat::CellData<double>& residual_data) const
    {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS5(d_dim, patch, flux_data, soln_data, rhs_data,
                                      residual_data);

      tbox::Pointer<geom::CartesianPatchGeometry> patch_geom =
        patch.getPatchGeometry();
      const hier::Box& box = patch.getBox();
      const int* lower = &box.lower()[0];
      const int* upper = &box.upper()[0];
      const double* dx = patch_geom->getDx();

      tbox::Pointer<pdat::CellData<double> > scalar_field_data;
      double scalar_field_constant;
      if (d_stokes_spec.cIsVariable()) {
        scalar_field_data =
          patch.getPatchData(d_stokes_spec.getCPatchDataId());
        if (d_dim == tbox::Dimension(2)) {
          F77_FUNC(compresvarsca2d, COMPRESVARSCA2D) (
                                                      flux_data.getPointer(0),
                                                      flux_data.getPointer(1),
                                                      &flux_data.getGhostCellWidth()[0],
                                                      &flux_data.getGhostCellWidth()[1],
                                                      rhs_data.getPointer(),
                                                      &rhs_data.getGhostCellWidth()[0],
                                                      &rhs_data.getGhostCellWidth()[1],
                                                      residual_data.getPointer(),
                                                      &residual_data.getGhostCellWidth()[0],
                                                      &residual_data.getGhostCellWidth()[1],
                                                      scalar_field_data->getPointer(),
                                                      &scalar_field_data->getGhostCellWidth()[0],
                                                      &scalar_field_data->getGhostCellWidth()[1],
                                                      soln_data.getPointer(),
                                                      &soln_data.getGhostCellWidth()[0],
                                                      &soln_data.getGhostCellWidth()[1],
                                                      &lower[0], &upper[0], &lower[1], &upper[1],
                                                      dx);
        } else if (d_dim == tbox::Dimension(3)) {
          F77_FUNC(compresvarsca3d, COMPRESVARSCA3D) (
                                                      flux_data.getPointer(0),
                                                      flux_data.getPointer(1),
                                                      flux_data.getPointer(2),
                                                      &flux_data.getGhostCellWidth()[0],
                                                      &flux_data.getGhostCellWidth()[1],
                                                      &flux_data.getGhostCellWidth()[2],
                                                      rhs_data.getPointer(),
                                                      &rhs_data.getGhostCellWidth()[0],
                                                      &rhs_data.getGhostCellWidth()[1],
                                                      &rhs_data.getGhostCellWidth()[2],
                                                      residual_data.getPointer(),
                                                      &residual_data.getGhostCellWidth()[0],
                                                      &residual_data.getGhostCellWidth()[1],
                                                      &residual_data.getGhostCellWidth()[2],
                                                      scalar_field_data->getPointer(),
                                                      &scalar_field_data->getGhostCellWidth()[0],
                                                      &scalar_field_data->getGhostCellWidth()[1],
                                                      &scalar_field_data->getGhostCellWidth()[2],
                                                      soln_data.getPointer(),
                                                      &soln_data.getGhostCellWidth()[0],
                                                      &soln_data.getGhostCellWidth()[1],
                                                      &soln_data.getGhostCellWidth()[2],
                                                      &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2],
                                                      dx);
        }
      } else if (d_stokes_spec.cIsConstant()) {
        scalar_field_constant = d_stokes_spec.getCConstant();
        if (d_dim == tbox::Dimension(2)) {
          F77_FUNC(compresconsca2d, COMPRESCONSCA2D) (
                                                      flux_data.getPointer(0),
                                                      flux_data.getPointer(1),
                                                      &flux_data.getGhostCellWidth()[0],
                                                      &flux_data.getGhostCellWidth()[1],
                                                      rhs_data.getPointer(),
                                                      &rhs_data.getGhostCellWidth()[0],
                                                      &rhs_data.getGhostCellWidth()[1],
                                                      residual_data.getPointer(),
                                                      &residual_data.getGhostCellWidth()[0],
                                                      &residual_data.getGhostCellWidth()[1],
                                                      scalar_field_constant,
                                                      soln_data.getPointer(),
                                                      &soln_data.getGhostCellWidth()[0],
                                                      &soln_data.getGhostCellWidth()[1],
                                                      &lower[0], &upper[0], &lower[1], &upper[1],
                                                      dx);
        } else if (d_dim == tbox::Dimension(3)) {
          F77_FUNC(compresconsca3d, COMPRESCONSCA3D) (
                                                      flux_data.getPointer(0),
                                                      flux_data.getPointer(1),
                                                      flux_data.getPointer(2),
                                                      &flux_data.getGhostCellWidth()[0],
                                                      &flux_data.getGhostCellWidth()[1],
                                                      &flux_data.getGhostCellWidth()[2],
                                                      rhs_data.getPointer(),
                                                      &rhs_data.getGhostCellWidth()[0],
                                                      &rhs_data.getGhostCellWidth()[1],
                                                      &rhs_data.getGhostCellWidth()[2],
                                                      residual_data.getPointer(),
                                                      &residual_data.getGhostCellWidth()[0],
                                                      &residual_data.getGhostCellWidth()[1],
                                                      &residual_data.getGhostCellWidth()[2],
                                                      scalar_field_constant,
                                                      soln_data.getPointer(),
                                                      &soln_data.getGhostCellWidth()[0],
                                                      &soln_data.getGhostCellWidth()[1],
                                                      &soln_data.getGhostCellWidth()[2],
                                                      &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2],
                                                      dx);
        }
      } else {
        scalar_field_constant = 0.0;
        if (d_dim == tbox::Dimension(2)) {
          F77_FUNC(compresconsca2d, COMPRESCONSCA2D) (
                                                      flux_data.getPointer(0),
                                                      flux_data.getPointer(1),
                                                      &flux_data.getGhostCellWidth()[0],
                                                      &flux_data.getGhostCellWidth()[1],
                                                      rhs_data.getPointer(),
                                                      &rhs_data.getGhostCellWidth()[0],
                                                      &rhs_data.getGhostCellWidth()[1],
                                                      residual_data.getPointer(),
                                                      &residual_data.getGhostCellWidth()[0],
                                                      &residual_data.getGhostCellWidth()[1],
                                                      0.0,
                                                      soln_data.getPointer(),
                                                      &soln_data.getGhostCellWidth()[0],
                                                      &soln_data.getGhostCellWidth()[1],
                                                      &lower[0], &upper[0], &lower[1], &upper[1],
                                                      dx);
        } else if (d_dim == tbox::Dimension(3)) {
          F77_FUNC(compresconsca3d, COMPRESCONSCA3D) (
                                                      flux_data.getPointer(0),
                                                      flux_data.getPointer(1),
                                                      flux_data.getPointer(2),
                                                      &flux_data.getGhostCellWidth()[0],
                                                      &flux_data.getGhostCellWidth()[1],
                                                      &flux_data.getGhostCellWidth()[2],
                                                      rhs_data.getPointer(),
                                                      &rhs_data.getGhostCellWidth()[0],
                                                      &rhs_data.getGhostCellWidth()[1],
                                                      &rhs_data.getGhostCellWidth()[2],
                                                      residual_data.getPointer(),
                                                      &residual_data.getGhostCellWidth()[0],
                                                      &residual_data.getGhostCellWidth()[1],
                                                      &residual_data.getGhostCellWidth()[2],
                                                      0.0,
                                                      soln_data.getPointer(),
                                                      &soln_data.getGhostCellWidth()[0],
                                                      &soln_data.getGhostCellWidth()[1],
                                                      &soln_data.getGhostCellWidth()[2],
                                                      &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2],
                                                      dx);
        }
      }
    }

  }
}
#endif
