/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "intersection.hxx"

void compute_intersections_3D(const FTensor::Tensor1<double,3> &ntt,
                              const FTensor::Tensor1<double,3> &xyz,
                              const FTensor::Tensor2<double,3,3> &rot,
                              const FTensor::Tensor1<double,3> dx[],
                              const double fault[],
                              const int &dim,
                              const int &ix,
                              int &intersect_diagonal,
                              int intersect_mixed[4],
                              int intersect_corner[4])
{
  intersect_diagonal=intersection(ntt,xyz,rot,dx[ix],fault,dim);

  FTensor::Tensor1<double,3> dx_2_y, dx_2_z;
  FTensor::Index<'a',3> a;

  int iy((ix+1)%dim), iz((ix+2)%dim);
      
  dx_2_y(a)=dx[iy](a)/2;
  intersect_mixed[0]=intersection(ntt,xyz,rot,dx_2_y,fault,dim);
  dx_2_y(a)=-dx_2_y(a);
  intersect_mixed[1]=intersection(ntt,xyz,rot,dx_2_y,fault,dim);

  dx_2_z(a)=dx[iz](a)/2;
  intersect_mixed[2]=intersection(ntt,xyz,rot,dx_2_z,fault,dim);
  dx_2_z(a)=-dx_2_z(a);
  intersect_mixed[3]=intersection(ntt,xyz,rot,dx_2_z,fault,dim);

  FTensor::Tensor1<double,3> dx_corner;
  dx_corner(ix)=0;
  dx_corner(iy)=dx[iy](iy)/2;
  dx_corner(iz)=dx[iz](iz)/2;
  intersect_corner[0]=intersection(ntt,xyz,rot,dx_corner,fault,dim);
  dx_corner(iz)=-dx_corner(iz);
  intersect_corner[1]=intersection(ntt,xyz,rot,dx_corner,fault,dim);

  dx_corner(iy)=-dx_corner(iy);
  intersect_corner[2]=intersection(ntt,xyz,rot,dx_corner,fault,dim);
  dx_corner(iz)=-dx_corner(iz);
  intersect_corner[3]=intersection(ntt,xyz,rot,dx_corner,fault,dim);
}
