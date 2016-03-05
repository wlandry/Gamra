/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "intersection.hxx"

void compute_intersections_2D(const FTensor::Tensor1<double,3> &ntt,
                              const FTensor::Tensor1<double,3> &xyz,
                              const FTensor::Tensor2<double,3,3> &rot,
                              const FTensor::Tensor1<double,3> dx[],
                              const double fault[],
                              const int &dim,
                              const int &ix,
                              int &intersect_diagonal,
                              int intersect_mixed[2])
{
  intersect_diagonal=intersection(ntt,xyz,rot,dx[ix],fault,dim);

  FTensor::Tensor1<double,3> dx_2;
  FTensor::Index<'a',3> a;
  const int iy((ix+1)%dim);
  dx_2(a)=dx[iy](a)/2;
  intersect_mixed[0]=intersection(ntt,xyz,rot,dx_2,fault,dim);

  dx_2(a)=-dx_2(a);
  intersect_mixed[1]=-intersection(ntt,xyz,rot,dx_2,fault,dim);
}
