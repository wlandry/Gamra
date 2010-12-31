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

#ifdef HAVE_HDF5
  /*
*************************************************************************
* Set up external plotter to plot internal data from this class.        *
* Register variables appropriate for plotting.                          *
*************************************************************************
*/
  int FACStokes::setupPlotter(
                              appu::VisItDataWriter& plotter) const {
    if (d_hierarchy.isNull()) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                 << " FACStokes::setupPlotter\n"
                 << "The hierarchy must be set before calling\n"
                 << "this function.\n");
    }
    plotter.registerPlotQuantity("Computed solution",
                                 "SCALAR",
                                 p_id);
    plotter.registerDerivedPlotQuantity("Error",
                                        "SCALAR",
                                        (appu::VisDerivedDataStrategy *)this);
    plotter.registerPlotQuantity("Exact solution",
                                 "SCALAR",
                                 p_exact_id);
    plotter.registerPlotQuantity("Stokes source",
                                 "SCALAR",
                                 p_rhs_id);

    return 0;
  }
#endif

}
