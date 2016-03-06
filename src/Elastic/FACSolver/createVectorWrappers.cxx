/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

void Elastic::FACSolver::createVectorWrappers(int v, int v_rhs)
{
  SAMRAI::hier::VariableDatabase&
    vdb(*SAMRAI::hier::VariableDatabase::getDatabase());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  if (!d_uv || d_uv->getComponentDescriptorIndex(0) != v)
    {
      d_uv.reset();
      d_uv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
        ("Elastic::FACSolver::uv", d_hierarchy, level_min, level_max);
      vdb.mapIndexToVariable(v, variable);
      if (!variable)
        { TBOX_ERROR(__FILE__ << ": No variable for patch data index "
                     << v << "\n"); }
      d_uv->addComponent(variable, v);
    }

  if (!d_fv || d_fv->getComponentDescriptorIndex(0) != v_rhs)
    {
      d_fv.reset();
      d_fv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
        ("Elastic::FACSolver::fv", d_hierarchy, level_min, level_max);
      vdb.mapIndexToVariable(v_rhs, variable);    
      if (!variable)
        { TBOX_ERROR(__FILE__ << ": No variable for patch data index "
                     << v_rhs << "\n"); }
      d_fv->addComponent(variable, v_rhs);
    }
}
