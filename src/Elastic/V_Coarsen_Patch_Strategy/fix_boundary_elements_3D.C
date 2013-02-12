#include "Elastic/V_Coarsen_Patch_Strategy.h"

void
Elastic::V_Coarsen_Patch_Strategy::fix_boundary_elements_3D
(SAMRAI::pdat::SideData<double>& v,
 const SAMRAI::pdat::SideData<double>& v_fine,
 const SAMRAI::pdat::SideData<double>& dv_mixed,
 const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox> &boundaries) const
{
  /* Fix up the boundary elements by iterating through the boundary
     boxes */

  /* We only care about faces, not edges or corners, so we only
     iterate over face boundary boxes. */

  SAMRAI::hier::Box gbox(v_fine.getGhostBox());
  SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  for(int mm=0; mm<boundaries.size(); ++mm)
    {
      SAMRAI::hier::Box bbox=boundaries[mm].getBox();
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
      int side(location_index%2==0 ? SAMRAI::pdat::SideIndex::Upper
               : SAMRAI::pdat::SideIndex::Lower);
      int dir2((direction+1)%3), dir3((direction+2)%3);
      SAMRAI::hier::Index yp(ip), zp(ip);
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

      SAMRAI::hier::Index
        lower=SAMRAI::hier::Index::coarsen(bbox.lower(),
                                           SAMRAI::hier::Index(2,2,2)),
        upper=SAMRAI::hier::Index::coarsen(bbox.upper(),
                                           SAMRAI::hier::Index(2,2,2));

      for(int k=lower(2); k<=upper(2); ++k)
        for(int j=lower(1); j<=upper(1); ++j)
          for(int i=lower(0); i<=upper(0); ++i)
            {
              SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j,k),
                                             direction,side);
              SAMRAI::pdat::SideIndex center(coarse*2);
              if(center[dir2]>=gbox.lower(dir2)
                 && center[dir2]<gbox.upper(dir2)
                 && center[dir3]>=gbox.lower(dir3)
                 && center[dir3]<gbox.upper(dir3))
                {
                  v(coarse)=
                    (v_fine(center) + v_fine(center+yp)
                     + v_fine(center+zp) + v_fine(center+yp+zp))/4;
                }
            }
    }
}
