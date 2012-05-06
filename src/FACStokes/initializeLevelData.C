/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "FACStokes.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "FTensor.hpp"

bool intersect_fault(const FTensor::Tensor1<double,3> &c0,
                     const FTensor::Tensor1<double,3> &c1, const double &L)
{
  double y((c1(1)*c0(0) -  c1(0)*c0(1))/(c1(0) - c0(0)));
  return y<=L/2 && y>-L/2;
}

void rotate(double cstrike, double sstrike, double cdip, double sdip, 
	      double x1, double x2, double x3, 
	      double * x1s, double * x1i, double * x2s, double * x3s, double * x3i)
{

	double x2r(cstrike*x1-sstrike*x2);

        (*x1s)= cdip*x2r-sdip*x3;
        (*x1i)= cdip*x2r+sdip*x3;
        (*x2s)= sstrike*x1+cstrike*x2;
        (*x3s)= sdip*x2r+cdip*x3;
        (*x3i)=-sdip*x2r+cdip*x3;

}

double gauss(double x, double sigma)
{
	const double pi2 = atan(1.0)*8;
	return exp(-0.5*(x/sigma)*(x/sigma))/sqrt(pi2)/sigma;
}

double omega(double x, double beta)
{
  const double pi=4*atan(1);

  if (fabs(x) <= (1-2*beta)/(1-beta)/2)
    {
      return 1;
    }
  else
    {
      if (fabs(x) < 1./(2.*(1.-beta)))
        {
          return pow(cos(pi*((1.-beta)*fabs(x)-0.5+beta)/(2.*beta)),2.);
        }
      else
        {
        return 0;
        }
    }
}

/*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
void SAMRAI::FACStokes::initializeLevelData
(const tbox::Pointer<hier::BasePatchHierarchy> patch_hierarchy,
 const int level_number,
 const double ,
 const bool ,
 const bool ,
 const tbox::Pointer<hier::BasePatchLevel> ,
 const bool allocate_data)
{
  tbox::Pointer<hier::PatchHierarchy> hierarchy = patch_hierarchy;
  tbox::Pointer<geom::CartesianGridGeometry> grid_geom =
    hierarchy->getGridGeometry();

  tbox::Pointer<hier::PatchLevel> level =
    hierarchy->getPatchLevel(level_number);
  const int dim=d_dim.getValue();

  if (allocate_data) {
    level->allocatePatchData(p_id);
    level->allocatePatchData(cell_moduli_id);
    level->allocatePatchData(edge_moduli_id);
    level->allocatePatchData(dp_id);
    level->allocatePatchData(p_rhs_id);
    level->allocatePatchData(p_exact_id);
    level->allocatePatchData(v_id);
    level->allocatePatchData(v_rhs_id);
  }

  /*
   * Initialize data in all patches in the level.
   */
  hier::PatchLevel::Iterator p_i(*level);
  for (p_i.initialize(*level); p_i; p_i++) {

    tbox::Pointer<hier::Patch> patch = *p_i;
    if (patch.isNull()) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find patch.  Null patch pointer.");
    }
    tbox::Pointer<geom::CartesianPatchGeometry>
      geom = patch->getPatchGeometry();
    const double *dx=geom->getDx();

    /* Initialize cell moduli */
    tbox::Pointer<pdat::CellData<double> > cell_moduli =
      patch->getPatchData(cell_moduli_id);

    hier::Box cell_moduli_box = cell_moduli->getBox();

    for(pdat::CellIterator ci(cell_moduli->getGhostBox()); ci; ci++)
      {
        pdat::CellIndex c=ci();
        double xyz[dim];
        for(int d=0;d<dim;++d)
          xyz[d]=geom->getXLower()[d]
            + dx[d]*(c[d]-cell_moduli_box.lower()[d] + 0.5);

        int ijk(0), factor(1);
        for(int d=0;d<dim;++d)
          {
            int i=static_cast<int>(xyz[d]*(lambda_ijk[d]-1)
                                   /(lambda_xyz_max[d]-lambda_xyz_min[d]));
            i=std::max(0,std::min(lambda_ijk[d]-1,i));
            ijk+=i*factor;
            factor*=lambda_ijk[d];
          }
        (*cell_moduli)(c,0)=lambda[ijk];
      }

    for(pdat::CellIterator ci(cell_moduli->getGhostBox()); ci; ci++)
      {
        pdat::CellIndex c=ci();
        double xyz[dim];
        for(int d=0;d<dim;++d)
          xyz[d]=geom->getXLower()[d]
            + dx[d]*(c[d]-cell_moduli_box.lower()[d] + 0.5);

        int ijk(0), factor(1);
        for(int d=0;d<dim;++d)
          {
            int i=static_cast<int>(xyz[d]*(mu_ijk[d]-1)
                                   /(mu_xyz_max[d]-mu_xyz_min[d]));
            i=std::max(0,std::min(mu_ijk[d]-1,i));
            ijk+=i*factor;
            factor*=mu_ijk[d];
          }
        (*cell_moduli)(c,1)=mu[ijk];
      }

    /* I do not think this is actually necessary. */
    tbox::Pointer<pdat::CellData<double> > dp_data =
      patch->getPatchData(dp_id);
    dp_data->fill(0.0);

    tbox::Pointer<pdat::CellData<double> > p_rhs_data =
      patch->getPatchData(p_rhs_id);
    p_rhs_data->fill(0.0);

    /* v_rhs */
    tbox::Pointer<pdat::SideData<double> > v_rhs_data =
      patch->getPatchData(v_rhs_id);

    if(v_rhs.empty())
      {
        v_rhs_data->fill(0,0);
      }
    else
      {
        double L(0.2),W(0.25);
        double cdip(0.707107),sdip(0.707107);
        //double cdip(0.),sdip(1.);
        double cstrike(0.),sstrike(1.);
        double cr(0.),sr(1.);
        double delta(0.003);
        double x(0.5),y(0.5),z(0.5);
        double beta(0.2);
	double scale(-10.);

        const double pi=4*atan(1);
        const double theta(pi/4);
        const FTensor::Tensor1<double,3> center(x,y,z);
        const FTensor::Tensor2<double,3,3> rot(std::cos(theta),std::sin(theta),0,
                                               -std::sin(theta),cos(theta),0,
                                               0,0,1);
        FTensor::Tensor1<double,3> Dx[dim];
        for(int d0=0;d0<dim;++d0)
          {
            for(int d1=0;d1<dim;++d1)
              Dx[d0](d1)=0;
            
            Dx[d0](d0)=dx[d0];
          }
        FTensor::Tensor1<double,3> slip(0,scale,0);
        FTensor::Index<'a',3> a;
        FTensor::Index<'b',3> b;
        FTensor::Tensor1<double,3> jump;
        jump(a)=slip(b)*rot(b,a);

        hier::Box pbox = v_rhs_data->getBox();
        for(int ix=0;ix<dim;++ix)
          {
            double offset[]={0.5,0.5,0.5};
            offset[ix]=0;

            for(pdat::SideIterator si(pbox,ix); si; si++)
              {
                pdat::SideIndex s=si();
                (*v_rhs_data)(s)=0;

                FTensor::Tensor1<double,3> xyz(0,0,0);
                for(int d=0;d<dim;++d)
                  xyz(d)=geom->getXLower()[d]
                    + dx[d]*(s[d]-pbox.lower()[d]+offset[d]);

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
                       && intersect_fault(ntt,ntt_dp[d],L))
                      sign=1;
                    else if(ntt(0)>0 && ntt_dp[d](0)<=0
                            && intersect_fault(ntt,ntt_dp[d],L))
                      sign=-1;
                    else if(ntt(0)<=0 && ntt_dm[d](0)>0
                            && intersect_fault(ntt,ntt_dm[d],L))
                      sign=1;
                    else if(ntt(0)>0 && ntt_dm[d](0)<=0
                            && intersect_fault(ntt,ntt_dm[d],L))
                      sign=-1;

                    if(sign!=0)
                      {
                        if(level_number==4)
                        tbox::pout << "ddv "
                                   << ix << " "
                                   << d << " "
                                   << sign << " "
                                   << jump(ix) << " "
                                   << ntt(0) << " "
                                   << ntt(1) << " "
                                   << ntt_dp[d](0) << " "
                                   << ntt_dp[d](1) << " "
                                   << ntt_dm[d](0) << " "
                                   << ntt_dm[d](1) << " "
                                   << "\n";

                        const double lambda_here(1), mu_here(1);
                        double factor(mu_here);
                        if(ix==d)
                          factor=lambda_here+2*mu_here;
                        (*v_rhs_data)(s)+=sign*factor*jump(ix)/(dx[d]*dx[d]);
                      }
                  }

                /* d/dxy */

                FTensor::Tensor1<double,3> ntt_dxy[2][2];
                ntt_dxy[0][0](a)=
                  rot(a,b)*(xyz(b)+Dx[0](b)/2+Dx[1](b)/2-center(b));
                ntt_dxy[0][1](a)=
                  rot(a,b)*(xyz(b)+Dx[0](b)/2-Dx[1](b)/2-center(b));
                ntt_dxy[1][0](a)=
                  rot(a,b)*(xyz(b)-Dx[0](b)/2+Dx[1](b)/2-center(b));
                ntt_dxy[1][1](a)=
                  rot(a,b)*(xyz(b)-Dx[0](b)/2-Dx[1](b)/2-center(b));

                const double lambda_here(1), mu_here(1);
                int l(0), m(0);
                double j(0);

                if(ntt_dxy[0][0](0)<=0 && ntt_dxy[1][0](0)>0
                   && intersect_fault(ntt_dxy[0][0],ntt_dxy[1][0],L))
                  m-=1;
                if(ntt_dxy[0][0](0)>0 && ntt_dxy[1][0](0)<=0
                   && intersect_fault(ntt_dxy[0][0],ntt_dxy[1][0],L))
                  m+=1;
                if(ntt_dxy[0][1](0)<=0 && ntt_dxy[1][1](0)>0
                   && intersect_fault(ntt_dxy[0][1],ntt_dxy[1][1],L))
                  m+=1;
                if(ntt_dxy[0][1](0)>0 && ntt_dxy[1][1](0)<=0
                   && intersect_fault(ntt_dxy[0][1],ntt_dxy[1][1],L))
                  m-=1;

                if(ntt_dxy[0][0](0)<=0 && ntt_dxy[0][1](0)>0
                   && intersect_fault(ntt_dxy[0][0],ntt_dxy[0][1],L))
                  l-=1;
                if(ntt_dxy[0][0](0)>0 && ntt_dxy[0][1](0)<=0
                   && intersect_fault(ntt_dxy[0][0],ntt_dxy[0][1],L))
                  l+=1;
                if(ntt_dxy[1][0](0)<=0 && ntt_dxy[1][1](0)>0
                   && intersect_fault(ntt_dxy[1][0],ntt_dxy[1][1],L))
                  l+=1;
                if(ntt_dxy[1][0](0)>0 && ntt_dxy[1][1](0)<=0
                   && intersect_fault(ntt_dxy[1][0],ntt_dxy[1][1],L))
                  l-=1;

                switch(ix)
                  {
                  case 0:
                    j=l*lambda_here + m*mu_here;
                    break;
                  case 1:
                    j=m*lambda_here + l*mu_here;
                    break;
                  default:
                    abort();
                  }

                if(j!=0)
                  {
                    (*v_rhs_data)(s)+=j*jump((ix+1)%dim) / (dx[0]*dx[1]);

                    tbox::pout << "intersect "
                               << level_number << " "
                               << ix << " "
                               << xyz(0) << " "
                               << xyz(1) << " "
                               << ntt(0) << " "
                               << ntt(1) << " "
                               << ntt_dxy[0][0](0) << " "
                               << ntt_dxy[0][0](1) << " "
                               << j << " "
                               << l << " "
                               << m << " "
                               << jump((ix+1)%dim) << " "
                               << (*v_rhs_data)(s) << " "
                               << "\n";
                  }                  

                // for(int dir0=0;dir0<2;++dir0)
                //   for(int dir1=0;dir1<2;++dir1)
                //     {

                //       if((ntt(0)<=0 && ntt_dxy[dir0][dir1](0)>0)
                //          || (ntt(0)>0 && ntt_dxy[dir0][dir1](0)<=0)
                //          && intersect_fault(ntt,ntt_dxy[dir0][dir1],L))
                //         {
                //           (*v_rhs_data)(s)+=(2*dir0-1)*(2*dir1-1)
                //             *jump((ix+1)%dim) / (dx[0]*dx[1]);

                //           tbox::pout << "intersect "
                //                      << level_number << " "
                //                      << ix << " "
                //                      << dir0 << " "
                //                      << dir1 << " "
                //                      << xyz(0) << " "
                //                      << xyz(1) << " "
                //                      << ntt(0) << " "
                //                      << ntt(1) << " "
                //                      << ntt_dxy[dir0][dir1](0) << " "
                //                      << ntt_dxy[dir0][dir1](1) << " "
                //                      << (*v_rhs_data)(s) << " "
                //                      << "\n";

                //         }
                //     }

                  //   switch(ix)
                  //     {
                  //     case 0:
                  //       if(ntt(1)<=L/2 && ntt_dp[1](1)>L/2)
                  //         sign=-1;
                  //       else if(ntt(1)<=-L/2 && ntt_dp[1](1)>-L/2)
                  //         sign=1;
                  //       else
                  //         sign=0;
                  //       (*v_rhs_data)(s)=scale*sign/(dx[0]*dx[1]);
                  //       break;
                  //     case 1:
                  //       if(sign!=0 && ntt(1)<=L/2 && ntt(1)>-L/2)
                  //         (*v_rhs_data)(s)=scale*sign/(dx[0]*dx[0]);
                  //       break;
                  //     case 2:
                  //       break;
                  //     default:
                  //       break;
                  //     }
                  // }
              }
          }
        if(0)
          {

	typedef struct {
		double x1s,x1i,x2s,x3s,x3i;
	} rpoint;

	rpoint xp00,xm00,x0p0,x0m0,x00p,x00m;
	double xr,yr,zr;
        double x2r;

	x2r= cstrike*x  -sstrike*y;
	xr = cdip   *x2r-sdip   *z;
	yr = sstrike*x  +cstrike*y;
	zr = sdip   *x2r+cdip   *z;

        hier::Box pbox = v_rhs_data->getBox();
        for(int ix=0;ix<dim;++ix)
          {
            double offset[]={0.5,0.5,0.5};
            offset[ix]=0;

            for(pdat::SideIterator si(pbox,ix); si; si++)
              {
                pdat::SideIndex s=si();
                double xyz[dim];
                for(int d=0;d<dim;++d)
                  xyz[d]=geom->getXLower()[d]
                    + dx[d]*(s[d]-pbox.lower()[d]+offset[d]);
            
		double x1,x2,x3;
		double dx1,dx2,dx3;

		if (2==d_dim.getValue())
		{
			x1=0;
			x2=xyz[0];
			x3=xyz[1];
			dx1=dx[0];
			dx2=dx[0];
			dx3=dx[1];
		} else if (3==d_dim.getValue())
		{
                        x1=xyz[0];
			x2=xyz[1];
			x3=xyz[2];
			dx1=dx[0];
			dx2=dx[1];
			dx3=dx[2];
		}

		double g0m0(1.),g0p0(1.),g00p(1.),g00m(1.),gm00(1.),gp00(1.);
		double x1s,x1i,x2s,x3s,x3i;

                x2r= cstrike*x1-sstrike*x2;
                x1s= cdip*x2r-sdip*x3;
                x1i= cdip*x2r+sdip*x3;
                x2s= sstrike*x1+cstrike*x2;
                x3s= sdip*x2r+cdip*x3;
                x3i=-sdip*x2r+cdip*x3;

		double n[]={cdip*cstrike,-cdip*sstrike,-sdip};
		double b[]={sstrike*cr+cstrike*sdip*sr,cstrike*cr-sstrike*sdip*sr,+cdip*sr};

		rotate(cstrike,sstrike,cdip,sdip,x1+dx1/2.,x2,x3,&(xp00.x1s),&(xp00.x1i),&(xp00.x2s),&(xp00.x3s),&(xp00.x3i));
		rotate(cstrike,sstrike,cdip,sdip,x1-dx1/2.,x2,x3,&(xm00.x1s),&(xm00.x1i),&(xm00.x2s),&(xm00.x3s),&(xm00.x3i));
		rotate(cstrike,sstrike,cdip,sdip,x1,x2+dx2/2.,x3,&(x0p0.x1s),&(x0p0.x1i),&(x0p0.x2s),&(x0p0.x3s),&(x0p0.x3i));
		rotate(cstrike,sstrike,cdip,sdip,x1,x2-dx2/2.,x3,&(x0m0.x1s),&(x0m0.x1i),&(x0m0.x2s),&(x0m0.x3s),&(x0m0.x3i));
		rotate(cstrike,sstrike,cdip,sdip,x1,x2,x3+dx3/2.,&(x00p.x1s),&(x00p.x1i),&(x00p.x2s),&(x00p.x3s),&(x00p.x3i));
		rotate(cstrike,sstrike,cdip,sdip,x1,x2,x3-dx3/2.,&(x00m.x1s),&(x00m.x1i),&(x00m.x2s),&(x00m.x3s),&(x00m.x3i));

		double temp1=gauss(x1s-xr,delta);
		double temp2(1.0);
		if (3==d_dim.getValue())
		{
			temp2=omega((x2s-yr)/L,beta);
		}
		double temp3=omega((x3s-zr)/W,beta);
		double sourc=(+(gp00*gauss(xp00.x1s-xr,delta)-gm00*gauss(xm00.x1s-xr,delta))*n[0]/dx1
                              +(g0p0*gauss(x0p0.x1s-xr,delta)-g0m0*gauss(x0m0.x1s-xr,delta))*n[1]/dx2
                              +(g00p*gauss(x00p.x1s-xr,delta)-g00m*gauss(x00m.x1s-xr,delta))*n[2]/dx3 )
			*temp2 
			*temp3;

                double dblcp=temp1 
                       *( (gp00*omega((xp00.x2s-yr)/L,beta)-gm00*omega((xm00.x2s-yr)/L,beta))*b[0]/dx1
                         +(g0p0*omega((x0p0.x2s-yr)/L,beta)-g0m0*omega((x0m0.x2s-yr)/L,beta))*b[1]/dx2
                         +(g00p*omega((x00p.x2s-yr)/L,beta)-g00m*omega((x00m.x2s-yr)/L,beta))*b[2]/dx3 ) 
                       *temp3;

		double dipcs=temp1 
			*temp2 
			*(+(gp00*omega((xp00.x3s-zr)/W,beta)-gm00*omega((xm00.x3s-zr)/W,beta))*b[0]/dx1
			  +(g0p0*omega((x0p0.x3s-zr)/W,beta)-g0m0*omega((x0m0.x3s-zr)/W,beta))*b[1]/dx2
			  +(g00p*omega((x00p.x3s-zr)/W,beta)-g00m*omega((x00m.x3s-zr)/W,beta))*b[2]/dx3 );

                temp1=gauss(x1i-xr,delta);
                temp3=omega((x3i+zr)/W,beta);
                // double image=( (gp00*gauss(xp00.x1i-xr,delta)-gm00*gauss(xm00.x1i-xr,delta))*n[0]/dx1
                //               +(g0p0*gauss(x0p0.x1i-xr,delta)-g0m0*gauss(x0m0.x1i-xr,delta))*n[1]/dx2
                //               +(g00p*gauss(x00p.x1i-xr,delta)-g00m*gauss(x00m.x1i-xr,delta))*n[2]/dx3 )
                //      *temp2 
                //      *temp3;
                // double cplei=temp1 
                //             *( (gp00*omega((xp00.x2s-yr)/L,beta)-gm00*omega((xp00.x2s-yr)/L,beta))*b[0]/dx1
                //               +(g0p0*omega((x0p0.x2s-yr)/L,beta)-g0m0*omega((x0p0.x2s-yr)/L,beta))*b[1]/dx2
                //               +(g00p*omega((x00p.x2s-yr)/L,beta)-g00m*omega((x00p.x2s-yr)/L,beta))*b[2]/dx3 ) 
                //             *temp3;
                // double dipci=temp1 
                //             *temp2 
                //             *( (gp00*omega((xp00.x3i+zr)/W,beta)-gm00*omega((xm00.x3i+zr)/W,beta))*b[0]/dx1
                //               +(g0p0*omega((x0p0.x3i+zr)/W,beta)-g0m0*omega((x0m0.x3i+zr)/W,beta))*b[1]/dx2
                //               +(g00p*omega((x00p.x3i+zr)/W,beta)-g00m*omega((x00m.x3i+zr)/W,beta))*b[2]/dx3 );
 
       		// force update
		switch (d_dim.getValue()-ix-1)
		{
			case 2:{
				// f_1
				double f1=+cr*sstrike*sourc
					  +cr*cdip*cstrike*dblcp
					  +sr*cdip*cstrike*dipcs
					  +sr*sdip*cstrike*sourc;
                		(*v_rhs_data)(s)=-f1*scale;
				break;
			       }
			case 1:{
			       // f_2
				double f2=+cr*cstrike*sourc
					  -cr*cdip*sstrike*dblcp
					  -sr*cdip*sstrike*dipcs
					  -sr*sdip*sstrike*sourc;
                		(*v_rhs_data)(s)=-f2*scale;
				break;
			       }
			case 0:{
				// f_3
				double f3=-cr*sdip*dblcp
					  +sr*cdip*sourc
					  -sr*sdip*dipcs;
				(*v_rhs_data)(s)=-f3*scale;
			       }
		}

              }
            int i=1;
            for(int d=0;d<dim;++d)
              i*=v_rhs_ijk[d];
          }
          }
      }
  }    // End patch loop.
}
