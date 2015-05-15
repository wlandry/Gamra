/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "Stokes/FAC.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

/*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
bool Stokes::FAC::packDerivedDataIntoDoubleBuffer(
                                                  double* buffer,
                                                  const SAMRAI::hier::Patch& patch,
                                                  const SAMRAI::hier::Box& region,
                                                  const std::string&
                                                  variable_name,
                                                  int depth_id,
                                                  double) const
{
  (void)depth_id;

  SAMRAI::pdat::CellIterator icell(SAMRAI::pdat::CellGeometry::begin(region));
  SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));

  if (variable_name == "Error") {
    boost::shared_ptr<SAMRAI::pdat::CellData<double> > current_solution_ =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
      (patch.getPatchData(p_id));
    boost::shared_ptr<SAMRAI::pdat::CellData<double> > exact_solution_ =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
      (patch.getPatchData(p_exact_id));
    SAMRAI::pdat::CellData<double>& current_solution = *current_solution_;
    SAMRAI::pdat::CellData<double>& exact_solution = *exact_solution_;
    for ( ; icell!=iend; ++icell) {
      double diff = (current_solution(*icell) - exact_solution(*icell));
      *buffer = diff;
      buffer = buffer + 1;
    }
  } else {
    // Did not register this name.
    TBOX_ERROR(
               "Unregistered variable name '" << variable_name << "' in\n"
               <<
               "Stokes::FACX::packDerivedDataIntoDoubleBuffer");

  }
  // Return true if this patch has derived data on it.
  // False otherwise.
  return true;
}
