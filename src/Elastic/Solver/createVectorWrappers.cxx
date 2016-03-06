/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Solver.hxx"

void Elastic::Solver::createVectorWrappers(int v, int v_rhs)
{
  SAMRAI::hier::VariableDatabase&
    vdb(*SAMRAI::hier::VariableDatabase::getDatabase());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  if (!uv || uv->getComponentDescriptorIndex(0) != v)
    {
      uv.reset();
      uv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
        ("Elastic::Solver::uv", hierarchy, level_min, level_max);
      vdb.mapIndexToVariable(v, variable);
      if (!variable)
        { TBOX_ERROR(__FILE__ << ": No variable for patch data index "
                     << v << "\n"); }
      uv->addComponent(variable, v);
    }

  if (!fv || fv->getComponentDescriptorIndex(0) != v_rhs)
    {
      fv.reset();
      fv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
        ("Elastic::Solver::fv", hierarchy, level_min, level_max);
      vdb.mapIndexToVariable(v_rhs, variable);    
      if (!variable)
        { TBOX_ERROR(__FILE__ << ": No variable for patch data index "
                     << v_rhs << "\n"); }
      fv->addComponent(variable, v_rhs);
    }
}
