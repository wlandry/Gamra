#include "V_Coarsen_Patch_Strategy.h"

void
SAMRAI::solv::V_Coarsen_Patch_Strategy::postprocessCoarsen_3D
(hier::Patch& coarse,
 const hier::Patch& fine,
 const hier::Box& ,
 const hier::IntVector& )
{
  /* Fix up the boundary elements by iterating through the boundary
     boxes */

  /* We only care about faces, not edges or corners, so we only
     iterate over face boundary boxes. */
  const tbox::Array<hier::BoundaryBox> &boundaries
    (coarse_fine[fine.getPatchLevelNumber()]
     ->getFaceBoundaries(coarse.getGlobalId()));
     
  tbox::Pointer<pdat::SideData<double> >
    v_fine = fine.getPatchData(v_id);
  tbox::Pointer<pdat::SideData<double> >
    v = coarse.getPatchData(v_id);

  TBOX_ASSERT(!v.isNull());
  TBOX_ASSERT(!v_fine.isNull());
  TBOX_ASSERT(v_fine->getDepth() == v->getDepth());
  TBOX_ASSERT(v->getDepth() == 1);

  hier::Box gbox(v_fine->getGhostBox());
  hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  for(int mm=0; mm<boundaries.size(); ++mm)
    {
      hier::Box bbox=boundaries[mm].getBox();
      /* location_index tells where, in relation to the box, the boundary is.
         0: x lower
         1: x upper
         2: y lower
         3: y upper
         4: z lower
         5: z upper

         Therefore, if location_index==3, then we need to set vy on
         the __lower__ side of that boundary box. */
         
      int location_index=boundaries[mm].getLocationIndex();
      int direction(location_index/2);
      int side(location_index%2==0 ? pdat::SideIndex::Upper
               : pdat::SideIndex::Lower);
      int dir2((direction+1)%3), dir3((direction+2)%3);
      hier::Index yp(ip), zp(ip);
      switch(direction)
        {
        case 0:
          yp=jp;
          zp=kp;
          break;
        case 1:
          yp=kp;
          zp=ip;
          break;
        case 2:
          yp=ip;
          zp=jp;
          break;
        }      

      hier::Index lower=hier::Index::coarsen(bbox.lower(),hier::Index(2,2,2)),
        upper=hier::Index::coarsen(bbox.upper(),hier::Index(2,2,2));

      for(int k=lower(2); k<=upper(2); ++k)
        for(int j=lower(1); j<=upper(1); ++j)
          for(int i=lower(0); i<=upper(0); ++i)
            {
              pdat::SideIndex coarse(hier::Index(i,j,k),direction,side);
              pdat::SideIndex center(coarse*2);
              if(center[dir2]>=gbox.lower(dir2)
                 && center[dir2]<gbox.upper(dir2)
                 && center[dir3]>=gbox.lower(dir3)
                 && center[dir3]<gbox.upper(dir3))
                {
                  (*v)(coarse)=
                    ((*v_fine)(center) + (*v_fine)(center+yp)
                     + (*v_fine)(center+zp) + (*v_fine)(center+yp+zp))/4;
                }
            }
    }
}
