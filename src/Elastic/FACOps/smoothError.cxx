/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::smoothError
(SAMRAI::solv::SAMRAIVectorReal<double>& data,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int ln,
 int num_sweeps)
{
  t_smooth_error->start();

  if(d_dim.getValue()==2)
    { Gauss_Seidel_red_black_2D(data,residual,ln,num_sweeps,
                                d_residual_tolerance_during_smoothing); }
  else if(d_dim.getValue()==3)
    { Gauss_Seidel_red_black_3D(data,residual,ln,num_sweeps,
                                d_residual_tolerance_during_smoothing); }
  else
    { TBOX_ERROR(d_object_name << ": Invalid dimension in Elastic::FACOps."); }
  t_smooth_error->stop();
}
