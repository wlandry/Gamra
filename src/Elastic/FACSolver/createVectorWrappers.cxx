/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

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
        (d_object_name + "::uv", d_hierarchy, d_ln_min, d_ln_max);
      /// Add v
      vdb.mapIndexToVariable(v, variable);
      if (!variable)
        TBOX_ERROR(d_object_name << ": No variable for patch data index "
                   << v << "\n");
      d_uv->addComponent(variable, v);
    }

  if (!d_fv || d_fv->getComponentDescriptorIndex(0) != v_rhs)
    {
      d_fv.reset();
      d_fv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
        (d_object_name + "::fv", d_hierarchy, d_ln_min, d_ln_max);
      /// Add v_rhs
      vdb.mapIndexToVariable(v_rhs, variable);    
      if (!variable)
        TBOX_ERROR(d_object_name << ": No variable for patch data index "
                   << v_rhs << "\n");
      d_fv->addComponent(variable, v_rhs);
    }
}