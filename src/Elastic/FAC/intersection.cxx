/* Returns whether the line between two grid points separated by dx
   intersects the fault. */

#include "Elastic/FAC.hxx"

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

  /* ntt is already in the coordinate frame of the fault.  So we could
     write this as ntt_dp(a)=ntt(a)+rot(a,b)*dx(b). Hmm. Maybe I
     should. Then I would not have to pass xyz. */

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
