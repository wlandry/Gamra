/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "FACStokes.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "StokesSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

  /*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
  bool FACStokes::packDerivedDataIntoDoubleBuffer(
                                                  double* buffer,
                                                  const hier::Patch& patch,
                                                  const hier::Box& region,
                                                  const std::string& variable_name,
                                                  int depth_id) const
  {
    (void)depth_id;

    pdat::CellData<double>::Iterator icell(region);

    if (variable_name == "Error") {
      tbox::Pointer<pdat::CellData<double> > current_solution_ =
        patch.getPatchData(p_id);
      tbox::Pointer<pdat::CellData<double> > exact_solution_ =
        patch.getPatchData(p_exact_id);
      pdat::CellData<double>& current_solution = *current_solution_;
      pdat::CellData<double>& exact_solution = *exact_solution_;
      for ( ; icell; icell++) {
        double diff = (current_solution(*icell) - exact_solution(*icell));
        *buffer = diff;
        buffer = buffer + 1;
      }
    } else {
      // Did not register this name.
      TBOX_ERROR(
                 "Unregistered variable name '" << variable_name << "' in\n"
                 <<
                 "FACStokesX::packDerivedDataIntoDoubleBuffer");

    }
    // Return true if this patch has derived data on it.
    // False otherwise.
    return true;
  }

}
