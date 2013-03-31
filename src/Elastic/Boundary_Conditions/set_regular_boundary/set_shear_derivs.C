#include "Constants.h"
#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_shear_derivs
(SAMRAI::pdat::SideData<double> &v,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_ptr,
 const SAMRAI::hier::Index unit[],
 const int &dim,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::hier::Box &gbox,
 const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
 const double *dx,
 const bool &homogeneous)
{
  for(int ix=0; ix<dim; ++ix)
    {
      double offset[]={0.5,0.5,0.5};
      offset[ix]=0;

      int xmin((!is_dirichlet[ix][ix][0]
                || !geom->getTouchesRegularBoundary(ix,0)) ? pbox.lower(ix)
               : pbox.lower(ix)+1),
        xmax((!is_dirichlet[ix][ix][1] 
              || !geom->getTouchesRegularBoundary(ix,1)) ? pbox.upper(ix)
             : pbox.upper(ix)-1);
      SAMRAI::hier::Box x_box(gbox);
      x_box.lower(ix)=xmin;
      x_box.upper(ix)=xmax;
      for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
        {
          const int ix_iy(index_map(ix,iy,dim));
          const int iy_ix(index_map(iy,ix,dim));
          if(geom->getTouchesRegularBoundary(iy,0) && !is_dirichlet[ix][iy][0]
             && x_box.lower(iy)<pbox.lower(iy))
            {
              SAMRAI::hier::Box box(x_box);
              box.upper(iy)=pbox.lower(iy)-1;
              SAMRAI::pdat::SideIterator end(box,ix,false);
              for(SAMRAI::pdat::SideIterator si(box,ix,true); si!=end; ++si)
                {
                  const SAMRAI::pdat::SideIndex &x(*si);
                  for(int d=0;d<dim;++d)
                    coord[d]=geom->getXLower()[d]
                      + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                  coord[iy]=geom->getXLower()[iy];

                  SAMRAI::pdat::SideIndex
                    y(x+unit[iy],iy,SAMRAI::pdat::SideIndex::Lower);
                  const double duyx((v(y)-v(y-unit[ix]))/dx[ix]);
                  v(x)=v(x+unit[iy]) + duyx*dx[iy];
                  if(!homogeneous)
                    {
                      if(have_faults())
                        v(x)+=((*dv_mixed_ptr)(y,iy_ix+1)
                               - (*dv_mixed_ptr)(y-unit[ix],iy_ix))
                          *dx[iy]/dx[ix]
                          + (*dv_mixed_ptr)(x+unit[iy],ix_iy+1);
                      v(x)-=expression[ix][iy][0].eval(coord)*dx[iy];
                    }
                }
            }
          if(geom->getTouchesRegularBoundary(iy,1) && !is_dirichlet[ix][iy][1]
             && x_box.upper(iy)>pbox.upper(iy))
            {
              SAMRAI::hier::Box box(x_box);
              box.lower(iy)=pbox.upper(iy)+1;
              SAMRAI::pdat::SideIterator end(box,ix,false);
              for(SAMRAI::pdat::SideIterator si(box,ix,true); si!=end; ++si)
                {
                  const SAMRAI::pdat::SideIndex &x(*si);
                  for(int d=0;d<dim;++d)
                    coord[d]=geom->getXLower()[d]
                      + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                  coord[iy]=geom->getXUpper()[iy];

                  SAMRAI::pdat::SideIndex
                    y(x,iy,SAMRAI::pdat::SideIndex::Lower);
                  const double duyx((v(y)- v(y-unit[ix]))/dx[ix]);
                  v(x)=v(x-unit[iy]) - duyx*dx[iy];

                  if(!homogeneous)
                    {
                      if(have_faults())
                        v(x)-=((*dv_mixed_ptr)(y,iy_ix+1)
                               - (*dv_mixed_ptr)(y-unit[ix],iy_ix))
                          *dx[iy]/dx[ix]
                          - (*dv_mixed_ptr)(x-unit[iy],ix_iy);
                      v(x)+=expression[ix][iy][1].eval(coord)*dx[iy];
                    }
                }
            }
        }
    }
}

