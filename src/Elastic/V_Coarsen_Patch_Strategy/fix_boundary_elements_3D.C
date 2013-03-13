#include "Elastic/V_Coarsen_Patch_Strategy.h"

void
Elastic::V_Coarsen_Patch_Strategy::fix_boundary_elements_3D
(SAMRAI::pdat::SideData<double>& v,
 const SAMRAI::pdat::SideData<double>& v_fine,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
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

      SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
      const SAMRAI::hier::Index unit[]={ip,jp,kp};
      const int dim(3);

      const int ix(location_index/2);
      const int iy((ix+1)%dim), iz((ix+2)%dim);
      int side(location_index%2==0 ? SAMRAI::pdat::SideIndex::Upper
               : SAMRAI::pdat::SideIndex::Lower);

      SAMRAI::hier::Index
        lower=SAMRAI::hier::Index::coarsen(bbox.lower(),
                                           SAMRAI::hier::Index(2,2,2)),
        upper=SAMRAI::hier::Index::coarsen(bbox.upper(),
                                           SAMRAI::hier::Index(2,2,2));

      int ijk[3];
      for(ijk[2]=lower(2); ijk[2]<=upper(2); ++ijk[2])
        for(ijk[1]=lower(1); ijk[1]<=upper(1); ++ijk[1])
          for(ijk[0]=lower(0); ijk[0]<=upper(0); ++ijk[0])
            {
              SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(ijk[0],ijk[1],
                                                                 ijk[2]),
                                             ix,side);
              SAMRAI::pdat::SideIndex fine(coarse*2);
              if(fine[iy]>=gbox.lower(iy) && fine[iy]<gbox.upper(iy)
                 && fine[iz]>=gbox.lower(iz) && fine[iz]<gbox.upper(iz))
                {
                  v(coarse)=coarsen_plane(v_fine,fine,unit[iy],unit[iz]);
                  if(have_faults)
                    v(coarse)+=coarsen_plane_correction(*dv_mixed,fine,unit[iy],
                                                        unit[iz]);
                }
            }
    }
}
