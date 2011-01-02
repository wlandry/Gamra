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

    /*
********************************************************************
* Check the validity and correctness of input data for this class. *
********************************************************************
*/

    void StokesFACOps::checkInputPatchDataIndices() const {
      /*
       * Check input validity and correctness.
       */
      hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());

      if (!d_stokes_spec.dIsConstant()
          && d_stokes_spec.getDPatchDataId() != -1) {
        tbox::Pointer<hier::Variable> var;
        vdb.mapIndexToVariable(d_stokes_spec.getDPatchDataId(), var);
        tbox::Pointer<pdat::SideVariable<double> > diffcoef_var = var;
        if (!diffcoef_var) {
          TBOX_ERROR(d_object_name
                     << ": Bad diffusion coefficient patch data index.");
        }
      }

      if (!d_stokes_spec.cIsConstant() && !d_stokes_spec.cIsZero()) {
        tbox::Pointer<hier::Variable> var;
        vdb.mapIndexToVariable(d_stokes_spec.getCPatchDataId(), var);
        tbox::Pointer<pdat::CellVariable<double> > scalar_field_var = var;
        if (!scalar_field_var) {
          TBOX_ERROR(d_object_name << ": Bad linear term patch data index.");
        }
      }

      if (d_flux_id != -1) {
        tbox::Pointer<hier::Variable> var;
        vdb.mapIndexToVariable(d_flux_id, var);
        tbox::Pointer<pdat::SideVariable<double> > flux_var = var;

        TBOX_ASSERT(flux_var);
      }

    }

  }
}
#endif
