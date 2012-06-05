/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "ElasticFACOps.h"

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
#include "ElasticHypreSolver.h"
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
* Compute the vector weight and put it at a specified patch data   *
* index.                                                           *
********************************************************************
*/

    void ElasticFACOps::computeVectorWeights(
                                            tbox::Pointer<hier::PatchHierarchy> hierarchy,
                                            int weight_id,
                                            int coarsest_ln,
                                            int finest_ln) const
    {
      TBOX_ASSERT(!hierarchy.isNull());
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

      if (coarsest_ln == -1) coarsest_ln = 0;
      if (finest_ln == -1) finest_ln = hierarchy->getFinestLevelNumber();
      if (finest_ln < coarsest_ln) {
        TBOX_ERROR(d_object_name
                   << ": Illegal level number range.  finest_ln < coarsest_ln.");
      }

      int ln;
      for (ln = finest_ln; ln >= coarsest_ln; --ln) {

        /*
         * On every level, first assign cell volume to vector weight.
         */

        tbox::Pointer<hier::PatchLevel> level =
          hierarchy->getPatchLevel(ln);
        for (hier::PatchLevel::Iterator p(level); p; p++) {
          tbox::Pointer<hier::Patch> patch = *p;
          tbox::Pointer<geom::CartesianPatchGeometry> patch_geometry =
            patch->getPatchGeometry();
          const double* dx = patch_geometry->getDx();
          double cell_vol = dx[0];
          if (d_dim > tbox::Dimension(1)) {
            cell_vol *= dx[1];
          }

          if (d_dim > tbox::Dimension(2)) {
            cell_vol *= dx[2];
          }

          tbox::Pointer<pdat::CellData<double> > w =
            patch->getPatchData(weight_id);
          if (!w) {
            TBOX_ERROR(d_object_name
                       << ": weight id must refer to a pdat::CellVariable");
          }
          w->fillAll(cell_vol);
        }

        /*
         * On all but the finest level, assign 0 to vector
         * weight to cells covered by finer cells.
         */

        if (ln < finest_ln) {

          /*
           * First get the boxes that describe index space of the next finer
           * level and coarsen them to describe corresponding index space
           * at this level.
           */

          tbox::Pointer<hier::PatchLevel> next_finer_level =
            hierarchy->getPatchLevel(ln + 1);
          hier::BoxArray coarsened_boxes = next_finer_level->getBoxes();
          hier::IntVector coarsen_ratio(next_finer_level->getRatioToLevelZero());
          coarsen_ratio /= level->getRatioToLevelZero();
          coarsened_boxes.coarsen(coarsen_ratio);

          /*
           * Then set vector weight to 0 wherever there is
           * a nonempty intersection with the next finer level.
           * Note that all assignments are local.
           */

          for (hier::PatchLevel::Iterator p(level); p; p++) {

            tbox::Pointer<hier::Patch> patch = *p;
            for (int i = 0; i < coarsened_boxes.getNumberOfBoxes(); i++) {

              hier::Box coarse_box = coarsened_boxes[i];
              hier::Box intersection = coarse_box * (patch->getBox());
              if (!intersection.empty()) {
                tbox::Pointer<pdat::CellData<double> > w =
                  patch->getPatchData(weight_id);
                w->fillAll(0.0, intersection);

              }  // assignment only in non-empty intersection
            }  // loop over coarsened boxes from finer level
          }  // loop over patches in level
        }  // all levels except finest
      }  // loop over levels
    }

  }
}
