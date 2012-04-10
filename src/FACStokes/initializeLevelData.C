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

double gauss(double x, double sigma)
{
	const double pi2 = atan(1.0)*8;
	return exp(-0.5*(x/sigma)*(x/sigma))/sqrt(pi2)/sigma;
}

double omega(double x, double beta)
{
  double pi=4*atan(1);

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
  hier::PatchLevel::Iterator pi(*level);
  for (pi.initialize(*level); pi; pi++) {

    tbox::Pointer<hier::Patch> patch = *pi;
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
        hier::Box pbox = v_rhs_data->getBox();
        int ix_offset(0);
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
            
		double x(xyz[0]),y(xyz[1]);

		double xr(0.5),yr(0.5);
		double g0m0(1.),g0p0(1.),g00p(1.),g00m(1.);
		double L(0.3),cdip(1.),sdip(0.);

		double beta=0.2;
		double n[]={0,cdip,sdip};
		double b[]={0,-sdip,+cdip};

		double temp1=gauss(x-xr,dx[0]);
		double temp2=1.0;
		double temp3=omega((y-yr)/L,beta);
		double sourc=(+(g0p0*gauss(x+dx[0]-xr,0.005)-g0m0*gauss(x-xr,0.005))*n[1]/dx[0]
                              +(g00p*gauss(x+dx[0]-xr,0.005)-g00m*gauss(x-xr,0.005))*n[2]/dx[1] )
			*temp2 
			*temp3;

		double dipcs=temp1 
			*temp2 
			*(+(g0p0*omega((y+dx[1]-yr)/L,beta)-g0m0*omega((y-yr)/L,beta))*b[1]/dx[0]
			  +(g00p*omega((y+dx[1]-yr)/L,beta)-g00m*omega((y-yr)/L,beta))*b[2]/dx[1] );

       		// force update
		switch (ix)
		{
			case 0:{
				// f_2
				double f2=-cdip*dipcs-sdip*sourc;
				       //std::cout << x << " " << y << " "<< omega((y-yr)/L,beta)<< " " << temp2 << " " << dipcs << "\n";
                		(*v_rhs_data)(s)=+10*f2;
				break;
			       }
			case 1:{
			       // f_3
				double f3=+cdip*sourc-sdip*dipcs;

                		(*v_rhs_data)(s)=-10*f3;
				break;
			       }
		}

              }
            int i=1;
            for(int d=0;d<dim;++d)
              i*=v_rhs_ijk[d];
            ix_offset+=i;
          }
      }
  }    // End patch loop.
}
