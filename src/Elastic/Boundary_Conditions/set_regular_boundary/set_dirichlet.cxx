#include "Constants.hxx"
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/EdgeData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include "Elastic/Boundary_Conditions.hxx"

void Elastic::Boundary_Conditions::set_dirichlet
(SAMRAI::pdat::SideData<double> &v,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_ptr,
 const SAMRAI::hier::Index unit[],
 const Gamra::Dir &dim,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::hier::Box &gbox,
 const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
 const double *dx,
 const bool &homogeneous) const
{
  for(Gamra::Dir ix=0; ix<dim; ++ix)
    {
      double offset[]={0.5,0.5,0.5};
      offset[ix]=0;

      for(Gamra::Dir iy=ix.next(dim); iy!=ix; iy=iy.next(dim))
        {
          const Gamra::Dir ix_iy(index_map(ix,iy,dim));

          SAMRAI::hier::Box x_box(gbox);
          x_box.setLower(ix,(geom->getTouchesRegularBoundary(ix,0)
                             && is_dirichlet[ix][ix][0]) ?
                         pbox.lower(ix)+1 : pbox.lower(ix));
          x_box.setUpper(ix,(geom->getTouchesRegularBoundary(ix,1)
                             && is_dirichlet[ix][ix][1]) ?
                         pbox.upper(ix)-1 : pbox.upper(ix));

          if(geom->getTouchesRegularBoundary(iy,0) && is_dirichlet[ix][iy][0])
            {
              SAMRAI::hier::Box y_box(x_box);
              y_box.setUpper(iy,y_box.lower(iy));

              SAMRAI::pdat::SideIterator
                end(SAMRAI::pdat::SideGeometry::end(y_box,ix));
              for(SAMRAI::pdat::SideIterator
                    si(SAMRAI::pdat::SideGeometry::begin(y_box,ix));
                  si!=end; ++si)
                {
                  const SAMRAI::pdat::SideIndex &x(*si);
                  v(x)=-v(x+unit[iy]);
                  if(!homogeneous)
                    {
                      std::vector<double> coord(dim);
                      coord[iy]=geom->getXLower()[iy];
                      for(int d=(iy+1)%dim;d!=iy;d=(d+1)%dim)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                      if(have_faults())
                        v(x)+= -(*dv_mixed_ptr)(x+unit[iy],ix_iy+1);
                      v(x)+=2*expression[ix][iy][0].eval(coord.data());
                    }
                }
            }
          if(geom->getTouchesRegularBoundary(iy,1) && is_dirichlet[ix][iy][1])
            {
              SAMRAI::hier::Box y_box(x_box);
              y_box.setLower(iy,y_box.upper(iy));

              SAMRAI::pdat::SideIterator
                end(SAMRAI::pdat::SideGeometry::end(y_box,ix));
              for(SAMRAI::pdat::SideIterator
                    si(SAMRAI::pdat::SideGeometry::begin(y_box,ix));
                  si!=end; ++si)
                {
                  const SAMRAI::pdat::SideIndex &x(*si);
                  v(x)=-v(x-unit[iy]);
                  if(!homogeneous)
                    {
                      std::vector<double> coord(dim);
                      coord[iy]=geom->getXUpper()[iy];
                      for(int d=(iy+1)%dim;d!=iy;d=(d+1)%dim)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                      if(have_faults())
                        v(x)+= -(*dv_mixed_ptr)(x-unit[iy],ix_iy);
                      v(x)+=2*expression[ix][iy][1].eval(coord.data());
                    }
                }
            }
        }
      
      SAMRAI::hier::Box box(gbox);

      if(geom->getTouchesRegularBoundary(ix,0) && is_dirichlet[ix][ix][0])
        {
          SAMRAI::hier::Box x_box(box);
          x_box.setUpper(ix,x_box.lower(ix));

          SAMRAI::pdat::SideIterator
            end(SAMRAI::pdat::SideGeometry::end(x_box,ix));
          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(x_box,ix)); si!=end; ++si)
            {
              const SAMRAI::pdat::SideIndex &x(*si);
              if(x[ix]<pbox.lower(ix))
                v(x)=boundary_value;
              else
                {
                  if(homogeneous)
                    {
                      v(x)=0;
                    }
                  else
                    {
                      std::vector<double> coord(dim);
                      for(int d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                      v(x)=expression[ix][ix][0].eval(coord.data());
                    }
                }
            }
        }
      if(geom->getTouchesRegularBoundary(ix,1) && is_dirichlet[ix][ix][1])
        {
          SAMRAI::hier::Box x_box(box);
          x_box.setLower(ix,x_box.upper(ix));

          SAMRAI::pdat::SideIterator
            end(SAMRAI::pdat::SideGeometry::end(x_box,ix));
          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(x_box,ix)); si!=end; ++si)
            {
              const SAMRAI::pdat::SideIndex &x(*si);
              if(x[ix]>pbox.upper(ix)+1)
                v(x)=boundary_value;
              else
                {
                  if(homogeneous)
                    {
                      v(x)=0;
                    }
                  else
                    {
                      std::vector<double> coord(dim);
                      for(int d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                      v(x)=expression[ix][ix][1].eval(coord.data());
                    }
                }
            }
        }
    }    
}
