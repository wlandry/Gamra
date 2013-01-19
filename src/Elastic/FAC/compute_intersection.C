#include "Elastic/FAC.h"

void Elastic::FAC::compute_intersection(const FTensor::Tensor1<double,3> &ntt,
                                        const FTensor::Tensor1<double,3> &xyz,
                                        const FTensor::Tensor2<double,3,3> &rot,
                                        const FTensor::Tensor1<double,3> dx[],
                                        const double fault[],
                                        const int &dim,
                                        int intersect[][2])
{
  FTensor::Tensor1<double,3> ntt_dp, ntt_dm;
  FTensor::Index<'a',3> a;
  FTensor::Index<'b',3> b;
  for(int d=0;d<dim;++d)
    {
      ntt_dp(a)=rot(a,b)*(xyz(b)+dx[d](b));
      ntt_dm(a)=rot(a,b)*(xyz(b)-dx[d](b));

      intersect[d][0]=0;
      if(ntt(0)<=0 && ntt_dp(0)>0
         && intersect_fault(dim,ntt,ntt_dp,fault))
        {
          intersect[d][0]=1;
        }
      else if(ntt(0)>0 && ntt_dp(0)<=0
              && intersect_fault(dim,ntt,ntt_dp,fault))
        {
          intersect[d][0]=-1;
        }

      intersect[d][1]=0;
      if(ntt(0)<=0 && ntt_dm(0)>0
         && intersect_fault(dim,ntt,ntt_dm,fault))
        {
          intersect[d][1]=-1;
        }
      else if(ntt(0)>0 && ntt_dm(0)<=0
              && intersect_fault(dim,ntt,ntt_dm,fault))
        {
          intersect[d][1]=1;
        }
    }
}
