#include "Elastic/FAC.h"

int Elastic::FAC::intersection(const FTensor::Tensor1<double,3> &ntt,
                               const FTensor::Tensor1<double,3> &xyz,
                               const FTensor::Tensor2<double,3,3> &rot,
                               const FTensor::Tensor1<double,3> &dx,
                               const double fault[],
                               const int &dim)
{
  FTensor::Tensor1<double,3> ntt_dp, ntt_dm;
  FTensor::Index<'a',3> a;
  FTensor::Index<'b',3> b;

  ntt_dp(a)=rot(a,b)*(xyz(b)+dx(b));

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
