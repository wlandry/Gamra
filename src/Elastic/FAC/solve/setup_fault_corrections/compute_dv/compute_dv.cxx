/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <vector>

#include <FTensor.hpp>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/CellData.h>

#include "Dir.hxx"

void compute_intersections_2D(const FTensor::Tensor1<double,3> &ntt,
                              const FTensor::Tensor2<double,3,3> &rot,
                              const FTensor::Tensor1<double,3> dx[],
                              const double fault[],
                              const int &dim,
                              const int &ix,
                              int &intersect_diagonal,
                              int intersect_mixed[2]);

void compute_intersections_3D(const FTensor::Tensor1<double,3> &ntt,
                              const FTensor::Tensor2<double,3,3> &rot,
                              const FTensor::Tensor1<double,3> dx[],
                              const double fault[],
                              const int &dim,
                              const int &ix,
                              int &intersect_diagonal,
                              int intersect_mixed[4],
                              int intersect_corner[4]);

void compute_dv(const std::vector<double> &faults,
                const Gamra::Dir &dim, const double *dx,
                const double *geom_xlower,
                const SAMRAI::hier::Box &gbox,
                const SAMRAI::hier::Index &pbox_lower,
                SAMRAI::pdat::CellData<double> &dv_diagonal,
                SAMRAI::pdat::SideData<double> &dv_mixed)
{
  dv_diagonal.fillAll(0);
  dv_mixed.fillAll(0);

  const int params_per_fault(9);
  for(size_t fault_index=0;fault_index<faults.size();
      fault_index+=params_per_fault)
    {
      const double *fault_ptr=faults.data()+fault_index;
      const double pi=4*atan(1);
      double scale(-fault_ptr[0]), x(fault_ptr[1]),
        y(fault_ptr[2]), z(fault_ptr[3]),
        L(fault_ptr[4]), W(fault_ptr[5]),
        strike(pi/2-fault_ptr[6]*pi/180),
        dip(fault_ptr[7]*pi/180),
        rake(fault_ptr[8]*pi/180);

      const FTensor::Tensor1<double,3> center(x,y,z);
      const FTensor::Tensor2<double,3,3>
        rot_strike(std::cos(strike),-std::sin(strike),0,
                   std::sin(strike),std::cos(strike),0,0,0,1),
        rot_dip(std::sin(dip),0,std::cos(dip),
                0,1,0,
                -std::cos(dip),0,std::sin(dip));

      FTensor::Tensor2<double,3,3> rot;
      FTensor::Index<'a',3> a;
      FTensor::Index<'b',3> b;
      FTensor::Index<'c',3> c;
      rot(a,b)=rot_dip(a,c)*rot_strike(c,b);

      std::vector<FTensor::Tensor1<double,3> > Dx(dim);
      for(int d0=0;d0<dim;++d0)
        {
          Dx[d0](a)=0.0;
          Dx[d0](d0)=dx[d0];
        }
      FTensor::Tensor1<double,3> slip(0,scale,0);
      const FTensor::Tensor2<double,3,3>
        rot_rake(1,0,0,
                 0,std::cos(rake),std::sin(rake),
                 0,-std::sin(rake),std::cos(rake));
      FTensor::Tensor1<double,3> jump;
      jump(c)=rot(b,c)*rot_rake(b,a)*slip(a);

      double fault[]={L,W};

      for(Gamra::Dir ix=0;ix<dim;++ix)
        {
          double offset[]={0.5,0.5,0.5};
          offset[ix]=0;

          SAMRAI::pdat::SideIterator
            s_end(SAMRAI::pdat::SideGeometry::end(gbox,ix));
          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(gbox,ix));
              si!=s_end; ++si)
            {
              const SAMRAI::pdat::SideIndex &s(*si);

              FTensor::Tensor1<double,3> xyz(0,0,0);
              for(int d=0;d<dim;++d)
                {
                  xyz(d)=geom_xlower[d]
                    + dx[d]*(s[d]-pbox_lower[d]+offset[d]) - center(d);
                }

              /// Rotate the coordinates into the coordinates
              /// of the fault.  So in those coordinates, if
              /// x<0, you are on the left, and if x>0, you
              /// are on the right.
              FTensor::Tensor1<double,3> ntt;
              ntt(a)=rot(a,b)*xyz(b);
              if(dim==2)
                {
                  int intersect, intersect_mixed[2];
                  compute_intersections_2D(ntt,rot,Dx.data(),fault,dim,ix,
                                           intersect,intersect_mixed);

                  if(gbox.contains(s))
                    {
                      SAMRAI::pdat::CellIndex c(s);
                      dv_diagonal(c,ix)+=intersect*jump(ix);
                    }

                  dv_mixed(s,0)+=intersect_mixed[0]*jump(ix);
                  dv_mixed(s,1)-=intersect_mixed[1]*jump(ix);
                }
              else
                {
                  int intersect_diagonal, intersect_mixed[4],
                    intersect_corner[4];
                  compute_intersections_3D(ntt,rot,Dx.data(),fault,dim,ix,
                                           intersect_diagonal,intersect_mixed,
                                           intersect_corner);

                  if(gbox.contains(s))
                    {
                      SAMRAI::pdat::CellIndex c(s);
                      dv_diagonal(c,ix)+=intersect_diagonal*jump(ix);
                    }

                  for(int n=0;n<4;++n)
                    {
                      dv_mixed(s,n)+=intersect_mixed[n]*jump(ix);
                      dv_mixed(s,n+4)+=intersect_corner[n]*jump(ix);
                    }
                }
            }
        }
    }
}
