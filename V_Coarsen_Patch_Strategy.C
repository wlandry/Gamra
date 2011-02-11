#include "V_Coarsen_Patch_Strategy.h"

void
SAMRAI::solv::V_Coarsen_Patch_Strategy::postprocessCoarsen
(hier::Patch& coarse,
 const hier::Patch& fine,
 const hier::Box& coarse_box,
 const hier::IntVector& ratio)
{
  /* Fix up the boundary elements by iterating through the boundary
     boxes */

  /* We only care about edges, not corners, so we only iterate over
     edge boundary boxes. */
  const tbox::Array<hier::BoundaryBox>
    &boundaries=coarse_fine[fine.getPatchLevelNumber()]->getEdgeBoundaries(coarse.getGlobalId());
     
  tbox::Pointer<pdat::SideData<double> >
    v_fine = fine.getPatchData(v_id);
  tbox::Pointer<pdat::SideData<double> >
    v = coarse.getPatchData(v_id);

  TBOX_ASSERT(!v.isNull());
  TBOX_ASSERT(!v_fine.isNull());
  TBOX_ASSERT(v_fine->getDepth() == v->getDepth());
  TBOX_ASSERT(v->getDepth() == 1);

  hier::Box gbox(v_fine->getGhostBox());
  hier::Index ip(1,0), jp(0,1);
  for(int mm=0; mm<boundaries.size(); ++mm)
    {
      hier::Box bbox=boundaries[mm].getBox();
      int location_index=boundaries[mm].getLocationIndex();

      // hier::Index unit(1,1);
      hier::Index lower=hier::Index::coarsen(bbox.lower(),hier::Index(2,2)),
        upper=hier::Index::coarsen(bbox.upper(),hier::Index(2,2));

      tbox::plog << "VCPS "
                 << mm << " "
                 << bbox.lower(0) << " "
                 << bbox.upper(0) << " "
                 << bbox.lower(1) << " "
                 << bbox.upper(1) << " "
                 << lower(0) << " "
                 << upper(0) << " "
                 << lower(1) << " "
                 << upper(1) << " "
                 << v_fine->getGhostBox().lower(0) << " "
                 << v_fine->getGhostBox().upper(0) << " "
                 << v_fine->getGhostBox().lower(1) << " "
                 << v_fine->getGhostBox().upper(1) << " "
                 << "\n";

      for(int j=lower(1); j<=upper(1); ++j)
        for(int i=lower(0); i<=upper(0); ++i)
          {
            tbox::plog << "V Patch "
                       << i << " "
                       << j << " ";
            /* Fix vx */
            if(location_index==0)
              {
                pdat::SideIndex coarse(hier::Index(i,j),0,
                                       pdat::SideIndex::Upper);
                pdat::SideIndex center(coarse*2);
                if(center[1]>=gbox.lower(1) && center[1]<gbox.upper(1))
                  {
                    (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/2;
                    tbox::plog << (*v)(coarse) << " ";
                  }
              }
            else if(location_index==1)
              {
                pdat::SideIndex coarse(hier::Index(i,j),0,
                                       pdat::SideIndex::Lower);
                pdat::SideIndex center(coarse*2);
                if(center[1]>=gbox.lower(1) && center[1]<gbox.upper(1))
                  {
                    (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/2;
                    tbox::plog << (*v)(coarse) << " ";
                  }
              }
            /* Fix vy */
            else if(location_index==2)
              {
                pdat::SideIndex coarse(hier::Index(i,j),1,
                                       pdat::SideIndex::Upper);
                pdat::SideIndex center(coarse*2);
                if(center[0]>=gbox.lower(0) && center[0]<gbox.upper(0))
                  {
                    (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+ip))/2;
                    tbox::plog << (*v)(coarse) << " ";
                  }
              }
            else if(location_index==3)
              {
                pdat::SideIndex coarse(hier::Index(i,j),1,
                                       pdat::SideIndex::Lower);
                pdat::SideIndex center(coarse*2);
                if(center[0]>=gbox.lower(0) && center[0]<gbox.upper(0))
                  {
                    (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+ip))/2;
                    tbox::plog << (*v)(coarse) << " ";
                  }
              }
            else
              {
                abort();
              }
            tbox::plog << "\n";
          }
    }
}
