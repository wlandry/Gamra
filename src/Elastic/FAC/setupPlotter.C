/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.h"

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

#ifdef HAVE_HDF5
/// Set up external plotter to plot internal data from this class.
/// Register variables appropriate for plotting.

void Elastic::FAC::setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const {
  if (!d_hierarchy) {
    TBOX_ERROR(d_object_name << ": No hierarchy in\n"
               << " Elastic::FAC::setupPlotter\n"
               << "The hierarchy must be set before calling\n"
               << "this function.\n");
  }
  plotter.registerDerivedPlotQuantity("Displacement",
                                      "VECTOR",
                                      (SAMRAI::appu::VisDerivedDataStrategy *)
                                      this);
  plotter.registerDerivedPlotQuantity("Fault Correction + RHS",
                                      "VECTOR",
                                      (SAMRAI::appu::VisDerivedDataStrategy *)
                                      this);
  plotter.registerPlotQuantity("Cell lambda",
                               "SCALAR",
                               cell_moduli_id,0);
  plotter.registerPlotQuantity("Cell mu",
                               "SCALAR",
                               cell_moduli_id,1);
  plotter.registerDerivedPlotQuantity("Strain","TENSOR",
                                      (SAMRAI::appu::VisDerivedDataStrategy *)
                                      this);
  if(v_initial[0].is_valid || v_initial[1].is_valid || v_initial[2].is_valid)
    plotter.registerDerivedPlotQuantity("Initial Displacement",
                                        "VECTOR",
                                        (SAMRAI::appu::VisDerivedDataStrategy *)
                                        this);

  if(have_embedded_boundary())
    {
      plotter.registerDerivedPlotQuantity("Level Set",
                                          "SCALAR",
                                          (SAMRAI::appu::VisDerivedDataStrategy *)
                                          this);
    }
}
#endif
