/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "FTensor.hpp"

bool intersect_fault(const int &dim,
                     const FTensor::Tensor1<double,3> &c0,
                     const FTensor::Tensor1<double,3> &c1,
                     const double fault[])
{
  bool result(true);
  for(int d=1;d<dim;++d)
    {
      double y((c1(d)*c0(0) -  c1(0)*c0(d))/(c0(0) - c1(0)));
      result=result && ((y<=fault[d-1] && y>0) || (y>=fault[d-1] && y<0));
    }
  return result;
}

/*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
void Elastic::FAC::initializeLevelData
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& patch_hierarchy,
 const int level_number,
 const double ,
 const bool ,
 const bool ,
 const boost::shared_ptr<SAMRAI::hier::PatchLevel>& ,
 const bool allocate_data)
{
  boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
    hierarchy = patch_hierarchy;
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (hierarchy->getGridGeometry());

  boost::shared_ptr<SAMRAI::hier::PatchLevel> level =
    hierarchy->getPatchLevel(level_number);
  const int dim=d_dim.getValue();

  if (allocate_data) {
    level->allocatePatchData(cell_moduli_id);
    level->allocatePatchData(edge_moduli_id);
    level->allocatePatchData(v_id);
    level->allocatePatchData(v_rhs_id);
  }

  /*
   * Initialize data in all patches in the level.
   */
  SAMRAI::hier::PatchLevel::Iterator p_i(level->begin());
  for (; p_i!=level->end(); p_i++) {

    boost::shared_ptr<SAMRAI::hier::Patch> patch = *p_i;
    if (!patch) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find patch.  Null patch pointer.");
    }
    boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
      boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
      (patch->getPatchGeometry());
    const double *dx=geom->getDx();

    /* Initialize cell moduli */
    boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
      (patch->getPatchData(cell_moduli_id));

    SAMRAI::hier::Box cell_moduli_box = cell_moduli->getBox();

    SAMRAI::pdat::CellIterator cend(cell_moduli->getGhostBox(),false);
    for(SAMRAI::pdat::CellIterator ci(cell_moduli->getGhostBox(),true);
        ci!=cend; ci++)
      {
        SAMRAI::pdat::CellIndex c=*ci;
        double xyz[3];
        for(int d=0;d<dim;++d)
          xyz[d]=geom->getXLower()[d]
            + dx[d]*(c[d]-cell_moduli_box.lower()[d] + 0.5);

        (*cell_moduli)(c,0)=lambda.eval(xyz);
        (*cell_moduli)(c,1)=mu.eval(xyz);
      }

    /* v_rhs */
    boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_data =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
      (patch->getPatchData(v_rhs_id));

    v_rhs_data->fill(0,0);
    /* FIXME: need to add in the v_rhs from the input file */

    /* Iterate over the faults */
    for(int fault_index=0;fault_index<faults.size();fault_index+=9)
      {
        const double pi=4*atan(1);
        /* The conventions for faults are different from the regular
         * xyz coordinates, so we have to convert.  Depth is a
         * coordinate, giving a left-handed coordinate system.  So we
         * invert the depth and width to convert to a right-handed
         * coordinate system.  Strike is opposite in direction, and
         * dip is measured from a plane lying flat. */
        double scale(faults[fault_index+0]), x(faults[fault_index+2]),
          y(faults[fault_index+1]), z(-faults[fault_index+3]),
          L(faults[fault_index+4]), W(-faults[fault_index+5]),
          strike(faults[fault_index+6]*pi/180),
          dip(faults[fault_index+7]*pi/180),
          rake(faults[fault_index+8]*pi/180);

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

        FTensor::Tensor1<double,3> Dx[dim];
        for(int d0=0;d0<dim;++d0)
          {
            for(int d1=0;d1<dim;++d1)
              Dx[d0](d1)=0;
            
            Dx[d0](d0)=dx[d0];
          }
        FTensor::Tensor1<double,3> slip(0,scale,0);
        const FTensor::Tensor2<double,3,3>
          rot_rake(1,0,0,
                   0,std::cos(rake),-std::sin(rake),
                   0,std::sin(rake),std::cos(rake));
        FTensor::Tensor1<double,3> jump;
        jump(c)=rot(b,c)*rot_rake(b,a)*slip(a);

        double fault[]={L,W};

        SAMRAI::hier::Box pbox = v_rhs_data->getBox();
        for(int ix=0;ix<dim;++ix)
          {
            double offset[]={0.5,0.5,0.5};
            offset[ix]=0;

            SAMRAI::pdat::SideIterator send(pbox,ix,false);
            for(SAMRAI::pdat::SideIterator si(pbox,ix,true); si!=send; si++)
              {
                SAMRAI::pdat::SideIndex s=*si;

                FTensor::Tensor1<double,3> xyz(0,0,0);
                for(int d=0;d<dim;++d)
                  xyz(d)=geom->getXLower()[d]
                    + dx[d]*(s[d]-pbox.lower()[d]+offset[d]);

                /* Rotate the coordinates into the coordinates of the
                   fault.  So in those coordinates, if x<0, you are on
                   the left, and if x>0, you are on the right. */
                FTensor::Tensor1<double,3> ntt;
                FTensor::Tensor1<double,3> ntt_dp[dim], ntt_dm[dim];
                ntt(a)=rot(a,b)*(xyz(b)-center(b));
                for(int d=0;d<dim;++d)
                  {
                    ntt_dp[d](a)=rot(a,b)*(xyz(b)+Dx[d](b)-center(b));
                    ntt_dm[d](a)=rot(a,b)*(xyz(b)-Dx[d](b)-center(b));
                  }

                /* d/dx^2, d/dy^2, d/dz^2 */
                for(int d=0;d<dim;++d)
                  {
                    int sign(0);
                    if(ntt(0)<=0 && ntt_dp[d](0)>0
                       && intersect_fault(dim,ntt,ntt_dp[d],fault))
                      sign=1;
                    else if(ntt(0)>0 && ntt_dp[d](0)<=0
                            && intersect_fault(dim,ntt,ntt_dp[d],fault))
                      sign=-1;
                    else if(ntt(0)<=0 && ntt_dm[d](0)>0
                            && intersect_fault(dim,ntt,ntt_dm[d],fault))
                      sign=1;
                    else if(ntt(0)>0 && ntt_dm[d](0)<=0
                            && intersect_fault(dim,ntt,ntt_dm[d],fault))
                      sign=-1;

                    if(sign!=0)
                      {
                        /* FIXME: lambda and mu should be calculated. */
                        const double lambda_here(1), mu_here(1);
                        double factor(mu_here);
                        if(ix==d)
                          factor=lambda_here+2*mu_here;
                        (*v_rhs_data)(s)+=sign*factor*jump(ix)/(dx[d]*dx[d]);
                      }
                  }

                /* d/dxy */
                for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                  {
                    FTensor::Tensor1<double,3> ntt_dxy[2][2];
                    ntt_dxy[0][0](a)=
                      rot(a,b)*(xyz(b)+Dx[ix](b)/2+Dx[iy](b)/2-center(b));
                    ntt_dxy[0][1](a)=
                      rot(a,b)*(xyz(b)+Dx[ix](b)/2-Dx[iy](b)/2-center(b));
                    ntt_dxy[1][0](a)=
                      rot(a,b)*(xyz(b)-Dx[ix](b)/2+Dx[iy](b)/2-center(b));
                    ntt_dxy[1][1](a)=
                      rot(a,b)*(xyz(b)-Dx[ix](b)/2-Dx[iy](b)/2-center(b));

                    /* FIXME: lambda and mu should be calculated. */
                    const double lambda_here(1), mu_here(1);
                    int l(0), m(0);
                    double j(0);

                    if(ntt_dxy[0][0](0)<=0 && ntt_dxy[1][0](0)>0
                       && intersect_fault(dim,ntt_dxy[0][0],ntt_dxy[1][0],fault))
                      m-=1;
                    if(ntt_dxy[0][0](0)>0 && ntt_dxy[1][0](0)<=0
                       && intersect_fault(dim,ntt_dxy[0][0],ntt_dxy[1][0],fault))
                      m+=1;
                    if(ntt_dxy[0][1](0)<=0 && ntt_dxy[1][1](0)>0
                       && intersect_fault(dim,ntt_dxy[0][1],ntt_dxy[1][1],fault))
                      m+=1;
                    if(ntt_dxy[0][1](0)>0 && ntt_dxy[1][1](0)<=0
                       && intersect_fault(dim,ntt_dxy[0][1],ntt_dxy[1][1],fault))
                      m-=1;

                    if(ntt_dxy[0][0](0)<=0 && ntt_dxy[0][1](0)>0
                       && intersect_fault(dim,ntt_dxy[0][0],ntt_dxy[0][1],fault))
                      l-=1;
                    if(ntt_dxy[0][0](0)>0 && ntt_dxy[0][1](0)<=0
                       && intersect_fault(dim,ntt_dxy[0][0],ntt_dxy[0][1],fault))
                      l+=1;
                    if(ntt_dxy[1][0](0)<=0 && ntt_dxy[1][1](0)>0
                       && intersect_fault(dim,ntt_dxy[1][0],ntt_dxy[1][1],fault))
                      l+=1;
                    if(ntt_dxy[1][0](0)>0 && ntt_dxy[1][1](0)<=0
                       && intersect_fault(dim,ntt_dxy[1][0],ntt_dxy[1][1],fault))
                      l-=1;

                    j=l*lambda_here + m*mu_here;

                    if(j!=0)
                      (*v_rhs_data)(s)+=j*jump(iy) / (dx[0]*dx[1]);
                  }
              }
          }
      }
  }    // End patch loop.
}
