#include "set_V_boundary.h"
#include "Boundary.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

void set_V_boundary(const SAMRAI::hier::Patch& patch, const int &v_id,
                    const bool &rhs)
{
  tbox::Pointer<pdat::SideData<double> > v_ptr = patch.getPatchData(v_id);
  pdat::SideData<double> &v(*v_ptr);

  hier::Box pbox=patch.getBox();
  hier::Box gbox=v.getGhostBox();

  tbox::Pointer<geom::CartesianPatchGeometry> geom = patch.getPatchGeometry();
  const tbox::Dimension Dim(patch.getDim());
  const int dim(patch.getDim().getValue());

  const hier::Index zero(hier::Index::getZeroIndex(Dim));
  hier::Index pp[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    pp[i][i]=1;
  /* This should really get read from the input file. */
  double lower_boundary[]={0,0,0};
  double upper_boundary[]={-1,1,0};

  for(int ix=0; ix<dim; ++ix)
    {
      for(pdat::SideIterator si(gbox,ix); si; si++)
        {
          pdat::SideIndex x(*si);

          /* Set a sentinel value for normal components */
          if((x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0))
             || (x[ix]>pbox.upper(ix)+1 && geom->getTouchesRegularBoundary(ix,1)))
            {
              v(x)=boundary_value;
            }
          /* Set values for normal components */
          else if(x[ix]==pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0)
                  && !rhs)
            {
              v(x)=lower_boundary[ix];
            }
          else if(x[ix]==pbox.upper(ix)+1 && geom->getTouchesRegularBoundary(ix,1)
                  && !rhs)
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
                    }
                  else if(x[iy]>pbox.upper(iy)
                          && geom->getTouchesRegularBoundary(iy,1))
                    {
                      v(x)=v(x-pp[iy]);
                    }
                }
            }
        }
    }
}
