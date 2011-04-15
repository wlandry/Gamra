#include "set_boundary.h"
#include "Boundary.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

void set_boundary(const SAMRAI::hier::Patch& patch, const int &p_id,
                    const int &v_id, const bool &rhs)
{
  hier::Box pbox=patch.getBox();

  tbox::Pointer<geom::CartesianPatchGeometry> geom = patch.getPatchGeometry();
  const tbox::Dimension Dim(patch.getDim());
  const int dim(patch.getDim().getValue());

  const hier::Index zero(hier::Index::getZeroIndex(Dim));
  hier::Index pp[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    pp[i][i]=1;
  /* This should really get read from the input file. */
  double lower_boundary[]={0,0,0};
  bool lower_dirichlet[]={true,true,true};

  // double upper_boundary[]={-6.94444444444e4,0,0};
  // bool upper_dirichlet[]={true,false,true};

  // bool upper_dirichlet[]={true,true,true};
  // double upper_boundary[]={-1,1,0};

  bool upper_dirichlet[]={true,true,true};
  double upper_boundary[]={0,0,0};

  if(p_id!=-1)
    {
      tbox::Pointer<pdat::CellData<double> > p_ptr = patch.getPatchData(p_id);
      pdat::CellData<double> &p(*p_ptr);

      hier::Box gbox=p.getGhostBox();
      
      for(pdat::CellIterator ci(gbox); ci; ci++)
        {
          pdat::CellIndex center(*ci);
          for(int ix=0; ix<dim; ++ix)
            {
              if(center[ix]<pbox.lower(ix)
                 && geom->getTouchesRegularBoundary(ix,0))
                {
                  p(center)=2*p(center+pp[ix])-p(center+pp[ix]*2);
                  if(!lower_dirichlet[ix])
                    p(center)=p(center+pp[ix]);
                }
              else if(center[ix]>pbox.upper(ix)
                      && geom->getTouchesRegularBoundary(ix,1))
                {
                  p(center)=2*p(center-pp[ix])-p(center-pp[ix]*2);
                  if(!upper_dirichlet[ix])
                    {
                      // p(center)=p(center-pp[ix]);
                      p(center)=0;
                    }
                }
            }
        }
    }


  tbox::Pointer<pdat::SideData<double> > v_ptr = patch.getPatchData(v_id);
  pdat::SideData<double> &v(*v_ptr);

  hier::Box gbox=v.getGhostBox();
  for(int ix=0; ix<dim; ++ix)
    {
      for(pdat::SideIterator si(gbox,ix); si; si++)
        {
          pdat::SideIndex x(*si);

          // double pos_x=geom->getXLower()[0]
          //   + geom->getDx()[0]*(x[0]-pbox.lower()[0]);

          /* Set a sentinel value for normal components */
          if(x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0))
            {
              if(lower_dirichlet[ix])
                v(x)=boundary_value;
              else
                v(x)=v(x+pp[ix]*2);
            }
          else if(x[ix]>pbox.upper(ix)+1 && geom->getTouchesRegularBoundary(ix,1))
            {
              if(upper_dirichlet[ix])
                v(x)=boundary_value;
              else
                v(x)=v(x-pp[ix]*2);
            }
          /* Set values for normal components */
          else if(x[ix]==pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0)
                  && !rhs && lower_dirichlet[ix])
            {
              v(x)=lower_boundary[ix];
            }
          else if(x[ix]==pbox.upper(ix)+1 && geom->getTouchesRegularBoundary(ix,1)
                  && !rhs && upper_dirichlet[ix])
            {
              v(x)=upper_boundary[ix];
            }
          /* Set derivatives for tangential component */
          else
            {
              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  if(x[iy]<pbox.lower(iy)
                     && geom->getTouchesRegularBoundary(iy,0))
                    {
                      v(x)=v(x+pp[iy]);

                      // if(ix==0 && iy==1)
                      //   {
                      //     if(pos_x<0.1 || rhs)
                      //       {
                      //         v(x)=-v(x+pp[iy]);
                      //       }
                      //     else
                      //       {
                      //         v(x)=-v(x+pp[iy]) + 2*upper_boundary[0];
                      //       }
                      //   }
                    }
                  else if(x[iy]>pbox.upper(iy)
                          && geom->getTouchesRegularBoundary(iy,1))
                    {
                      v(x)=v(x-pp[iy]);

                      // if(ix==1 && iy==0)
                      //   {
                      //     v(x)=-v(x-pp[iy]);
                      //   }
                    }
                }
            }
        }
    }
}
