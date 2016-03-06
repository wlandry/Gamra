/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

/// Returns whether the line between two grid points separated by dx
/// intersects the fault.

#include <FTensor.hpp>

bool intersect_fault(const int &dim,
                     const FTensor::Tensor1<double,3> &c0,
                     const FTensor::Tensor1<double,3> &c1,
                     const double fault[]);

int intersection(const FTensor::Tensor1<double,3> &ntt,
                 const FTensor::Tensor2<double,3,3> &rot,
                 const FTensor::Tensor1<double,3> &dx,
                 const double fault[],
                 const int &dim)
{
  FTensor::Tensor1<double,3> ntt_dp;
  FTensor::Index<'a',3> a;
  FTensor::Index<'b',3> b;

  ntt_dp(a)=ntt(a) + rot(a,b)*dx(b);

  int result=0;
  if(ntt(0)<=0 && ntt_dp(0)>0 && intersect_fault(dim,ntt,ntt_dp,fault))
    {
      result=1;
    }
  else if(ntt(0)>0 && ntt_dp(0)<=0 && intersect_fault(dim,ntt,ntt_dp,fault))
    {
      result=-1;
    }
  return result;
}
