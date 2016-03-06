/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <FTensor.hpp>

bool intersect_fault(const int &dim,
                     const FTensor::Tensor1<double,3> &c0,
                     const FTensor::Tensor1<double,3> &c1,
                     const double fault[])
{
  // FIXME: Need to check whether the intersection happens within
  // the boundary (level_set>=0)
  bool result(true);
  for(int d=1;d<dim;++d)
    {
      double y((c1(d)*c0(0) -  c1(0)*c0(d))/(c0(0) - c1(0)));
      result=result && ((y<=fault[d-1] && y>0) || (y>=fault[d-1] && y<0));
    }
  return result;
}

