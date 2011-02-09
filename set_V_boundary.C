#include "set_V_boundary.h"
#include "Boundary.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

void set_V_boundary(const SAMRAI::hier::Patch& patch, const int &v_id)
{
  tbox::Pointer<pdat::SideData<double> >
    v = patch.getPatchData(v_id);

  hier::Box pbox=patch.getBox();
  hier::Box gbox=v->getGhostBox();

  tbox::Pointer<geom::CartesianPatchGeometry>
    geom = patch.getPatchGeometry();

  for(int j=gbox.lower(1); j<=gbox.upper(1)+1; ++j)
    for(int i=gbox.lower(0); i<=gbox.upper(0)+1; ++i)
      {
        pdat::CellIndex center(tbox::Dimension(2));
        center[0]=i;
        center[1]=j;
        hier::Index ip(1,0), jp(0,1);

        /* vx */
        if(j<=gbox.upper(1))
          {
            /* Set a sentinel value */
            if((i<pbox.lower(0) && geom->getTouchesRegularBoundary(0,0))
               || (i>pbox.upper(0)+1 && geom->getTouchesRegularBoundary(0,1)))
              {
                (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                     pdat::SideIndex::Lower))=
                  boundary_value;
              }
            /* Set the value so the derivative=0 */
            else if(j<pbox.lower(1) && geom->getTouchesRegularBoundary(1,0))
              {
                (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                     pdat::SideIndex::Lower))=
                  (*v)(pdat::SideIndex(center+jp,pdat::SideIndex::X,
                                       pdat::SideIndex::Lower));
              }
            else if(j>pbox.upper(1) && geom->getTouchesRegularBoundary(1,1))
              {
                (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                     pdat::SideIndex::Lower))=
                  (*v)(pdat::SideIndex(center-jp,pdat::SideIndex::X,
                                       pdat::SideIndex::Lower));
              }
          }
        /* vy */
        if(i<=gbox.upper(0))
          {
            if((j<pbox.lower(1) && geom->getTouchesRegularBoundary(1,0))
               || (j>pbox.upper(1)+1 && geom->getTouchesRegularBoundary(1,1)))
              {
                (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                     pdat::SideIndex::Lower))=
                  boundary_value;
              }
            else if(i<pbox.lower(0) && geom->getTouchesRegularBoundary(0,0))
              {
                (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                     pdat::SideIndex::Lower))=
                  (*v)(pdat::SideIndex(center+ip,pdat::SideIndex::Y,
                                       pdat::SideIndex::Lower));
              }
            else if(i>pbox.upper(0) && geom->getTouchesRegularBoundary(0,1))
              {
                (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                     pdat::SideIndex::Lower))=
                  (*v)(pdat::SideIndex(center-ip,pdat::SideIndex::Y,
                                       pdat::SideIndex::Lower));
              }
          }
      }
}
