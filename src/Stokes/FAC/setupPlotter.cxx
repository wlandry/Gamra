/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "Stokes/FAC.hxx"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/solv/SimpleCellRobinBcCoefs.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/math/HierarchyCellDataOpsReal.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/hier/VariableDatabase.h>

#ifdef HAVE_HDF5
/*
*************************************************************************
* Set up external plotter to plot internal data from this class.        *
* Register variables appropriate for plotting.                          *
*************************************************************************
*/
void Stokes::FAC::setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const {
  if (!d_hierarchy) {
    TBOX_ERROR("Stokes::FAC: No hierarchy in\n"
               << " Stokes::FAC::setupPlotter\n"
               << "The hierarchy must be set before calling\n"
               << "this function.\n");
  }
  plotter.registerPlotQuantity("Pressure",
                               "SCALAR",
                               p_id);
  plotter.registerDerivedPlotQuantity("Error",
                                      "SCALAR",
                                      (SAMRAI::appu::VisDerivedDataStrategy *)this);
  plotter.registerPlotQuantity("Exact solution",
                               "SCALAR",
                               p_exact_id);
  plotter.registerPlotQuantity("Cell Viscosity",
                               "SCALAR",
                               cell_viscosity_id);
  plotter.registerPlotQuantity("Stokes source",
                               "SCALAR",
                               p_rhs_id);
}
#endif
