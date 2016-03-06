/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"

void Elastic::FAC::setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const
{
  if (!d_hierarchy)
    {
      TBOX_ERROR("Elastic::FAC: No hierarchy in\n"
                 << " Elastic::FAC::setupPlotter\n"
                 << "The hierarchy must be set before calling\n"
                 << "this function.\n");
    }
  plotter.registerDerivedPlotQuantity
    ("Displacement","VECTOR",(SAMRAI::appu::VisDerivedDataStrategy *) this);
  plotter.registerDerivedPlotQuantity
    ("Fault Correction + RHS","VECTOR",
     (SAMRAI::appu::VisDerivedDataStrategy *) this);
  plotter.registerPlotQuantity("Cell lambda","SCALAR",cell_moduli_id,0);
  plotter.registerPlotQuantity("Cell mu","SCALAR",cell_moduli_id,1);
  plotter.registerDerivedPlotQuantity
    ("Strain","TENSOR", (SAMRAI::appu::VisDerivedDataStrategy *) this);
  if(v_initial[0].is_valid || v_initial[1].is_valid || v_initial[2].is_valid)
    {
      plotter.registerDerivedPlotQuantity
        ("Initial Displacement", "VECTOR",
         (SAMRAI::appu::VisDerivedDataStrategy *) this);
    }

  if(have_embedded_boundary())
    {
      plotter.registerDerivedPlotQuantity
        ("Level Set","SCALAR",(SAMRAI::appu::VisDerivedDataStrategy *) this);
    }
}
