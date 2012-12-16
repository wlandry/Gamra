#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_shear_derivs
(SAMRAI::pdat::SideData<double> &v,
 SAMRAI::hier::Index pp[], const int &dim,
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
      SAMRAI::pdat::SideIterator send(gbox,ix,false);
      for(SAMRAI::pdat::SideIterator si(gbox,ix,true); si!=send; si++)
        {
          SAMRAI::pdat::SideIndex x(*si);
              
          /* Set the derivative.  Do not set a shear derivative at the
           * corner where a dirichlet condition is set. */
          if((x[ix]>pbox.lower(ix)
              || (x[ix]>=pbox.lower(ix)
                  && (!is_dirichlet[ix][ix][0]
                      || !geom->getTouchesRegularBoundary(ix,0))))
             && (x[ix]<pbox.upper(ix)+1
                 || (x[ix]<=pbox.upper(ix)+1
                     && (!is_dirichlet[ix][ix][1] 
                         || !geom->getTouchesRegularBoundary(ix,1)))))
            {
              for(int d=0;d<dim;++d)
                coord[d]=geom->getXLower()[d]
                  + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  if(x[iy]<pbox.lower(iy)
                     && geom->getTouchesRegularBoundary(iy,0))
                    {
                      if(!is_dirichlet[ix][iy][0])
                        {
                          double coord_save(geom->getXLower()[iy]);
                          std::swap(coord[iy],coord_save);
                          SAMRAI::pdat::SideIndex
                            y(x+pp[iy],iy,SAMRAI::pdat::SideIndex::Lower);
                          const double duyx((v(y)-v(y-pp[ix]))/dx[ix]);
                          v(x)=v(x+pp[iy]) + duyx*dx[iy];
                          if(!homogeneous)
                            {
                              v(x)-=expression[ix][iy][0].eval(coord)*dx[iy];
                            }
                          std::swap(coord[iy],coord_save);
                        }
                    }
                  else if(x[iy]>pbox.upper(iy)
                          && geom->getTouchesRegularBoundary(iy,1))
                    {
                      if(!is_dirichlet[ix][iy][1])
                        {
                          double coord_save(geom->getXUpper()[iy]);
                          std::swap(coord[iy],coord_save);
                          SAMRAI::pdat::SideIndex
                            y(x,iy,SAMRAI::pdat::SideIndex::Lower);
                          const double duyx((v(y)-v(y-pp[ix]))/dx[ix]);
                          v(x)=v(x-pp[iy]) - duyx*dx[iy];
                          if(!homogeneous)
                            {
                              v(x)+=expression[ix][iy][1].eval(coord)*dx[iy];
                            }
                          std::swap(coord[iy],coord_save);
                        }
                    }
                }
            }
        }
    }
}

