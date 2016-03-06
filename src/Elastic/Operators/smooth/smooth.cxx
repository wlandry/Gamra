/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::smooth
(SAMRAI::solv::SAMRAIVectorReal<double>& data,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 const int &ln,
 const int &num_sweeps,
 const double &tolerance)
{
  t_smooth_error->start();

  if(dimension.getValue()==2)
    { Gauss_Seidel_red_black_2D(data,residual,ln,num_sweeps,tolerance); }
  else if(dimension.getValue()==3)
    { Gauss_Seidel_red_black_3D(data,residual,ln,num_sweeps,tolerance); }
  else
    { TBOX_ERROR(__FILE__ << ": Invalid dimension in Elastic::Operators."); }
  t_smooth_error->stop();
}
