#include "V_Coarsen_Patch_Strategy.h"
#include "set_V_boundary.h"

void
SAMRAI::solv::V_Coarsen_Patch_Strategy::preprocessCoarsen
(hier::Patch& coarse,
 const hier::Patch& fine,
 const hier::Box& coarse_box,
 const hier::IntVector& ratio)
{
  set_V_boundary(fine,v_id);

  // tbox::Pointer<pdat::SideData<double> >
  //   v = fine.getPatchData(v_id);

  // hier::Box pbox=fine.getBox();

  // hier::Box gbox=v->getGhostBox();

  // // tbox::plog << "boundary "
  // //            << gbox.lower(0) << " "
  // //            << gbox.upper(0) << " "
  // //            << gbox.lower(1) << " "
  // //            << gbox.upper(1) << " "
  // //            << pbox.lower(0) << " "
  // //            << pbox.upper(0) << " "
  // //            << pbox.lower(1) << " "
  // //            << pbox.upper(1) << " "
  // //            << "\n";
  // for(int j=gbox.lower(1); j<=gbox.upper(1)+1; ++j)
  //   for(int i=gbox.lower(0); i<=gbox.upper(0)+1; ++i)
  //     {
  //       pdat::CellIndex center(tbox::Dimension(2));
  //       center[0]=i;
  //       center[1]=j;
  //       hier::Index ip(1,0), jp(0,1);

  //       /* vx */
  //       if(j<=gbox.upper(1))
  //         {
  //           /* Set a sentinel value */
  //           if(i<pbox.lower(0) || i>pbox.upper(0)+1)
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                    pdat::SideIndex::Lower))=
  //                 boundary_value;
  //             }
  //           /* Set the value so the derivative=0 */
  //           else if(j<pbox.lower(0))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center+jp,pdat::SideIndex::X,
  //                                      pdat::SideIndex::Lower));
  //           tbox::plog << "V Coarsen Patch vx lower "
  //                      << fine.getPatchLevelNumber() << " "
  //                      << i << " "
  //                      << j << " "
  //                      << (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                              pdat::SideIndex::Lower))
  //                      << " "
  //                      << (*v)(pdat::SideIndex(center+jp,pdat::SideIndex::X,
  //                                              pdat::SideIndex::Lower))
  //                      << " "
  //                      << gbox.lower(0) << " "
  //                      << pbox.lower(0) << " "
  //                      << "\n";
  //             }
  //           else if(j>pbox.upper(0))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center-jp,pdat::SideIndex::X,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           tbox::plog << "V Coarsen Patch vx "
  //                      << fine.getPatchLevelNumber() << " "
  //                      << i << " "
  //                      << j << " "
  //                      << (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                              pdat::SideIndex::Lower))
  //                      << " "
  //                      << (&(*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                                pdat::SideIndex::Lower)))
  //                      << " "
  //                      << "\n";
  //         }
  //       /* vy */
  //       if(i<=gbox.upper(0))
  //         {
  //           if(j<pbox.lower(1) || j>pbox.upper(1)+1)
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                    pdat::SideIndex::Lower))=
  //                 boundary_value;
  //             }
  //           else if(i<pbox.lower(0))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center+ip,pdat::SideIndex::Y,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           else if(i>pbox.upper(0))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center-ip,pdat::SideIndex::Y,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           tbox::plog << "V Coarsen Patch vy "
  //                      << fine.getPatchLevelNumber() << " "
  //                      << i << " "
  //                      << j << " "
  //                      << (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                              pdat::SideIndex::Lower))
  //                      << " "
  //                      << (&(*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                                pdat::SideIndex::Lower)))
  //                      << " "
  //                      << "\n";
  //         }
  //     }
}


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

  tbox::plog << "V Coarsen Patch Strategy "
             << boundaries.size() << " "
             << "\n";

  hier::Index ip(1,0), jp(0,1);
  for(int mm=0; mm<boundaries.size(); ++mm)
    {
      hier::Box bbox=boundaries[mm].getBox();
      int location_index=boundaries[mm].getLocationIndex();

      hier::Index lower=hier::Index::coarsen(bbox.lower(),hier::Index(2,2)),
        upper=hier::Index::coarsen(bbox.upper(),hier::Index(2,2));

      tbox::plog << "BBox "
                 << mm << " "
                 << bbox.lower(0) << " "
                 << bbox.upper(0) << " "
                 << bbox.lower(1) << " "
                 << bbox.upper(1) << " "
                 << location_index << " "
                 << lower(0) << " "
                 << upper(0) << " "
                 << lower(1) << " "
                 << upper(1) << " "
                 << "\n";

      for(int j=lower(1); j<=upper(1); ++j)
        for(int i=lower(0); i<=upper(0); ++i)
          {
            tbox::plog << "VCPS "
                       << i << " "
                       << j << " ";

            /* Fix vx */
            if(location_index==0)
              {
                // if(!cgeom->getTouchesRegularBoundary(0,0))
                {
                  pdat::SideIndex coarse(hier::Index(i,j),0,
                                         pdat::SideIndex::Upper);
                  pdat::SideIndex center(coarse*2);
                  (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/2;
                  tbox::plog << (*v)(coarse) << " ";
                }
              }
            else if(location_index==1)
              {
                // if(!cgeom->getTouchesRegularBoundary(0,1))
                {
                  pdat::SideIndex coarse(hier::Index(i,j),0,
                                         pdat::SideIndex::Lower);
                  pdat::SideIndex center(coarse*2);
                  (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/2;
                  tbox::plog << (*v)(coarse) << " ";
                }
              }
            /* Fix vy */
            else if(location_index==2)
              {
                // if(!cgeom->getTouchesRegularBoundary(1,0))
                {
                  pdat::SideIndex coarse(hier::Index(i,j),1,
                                         pdat::SideIndex::Upper);
                  pdat::SideIndex center(coarse*2);
                  (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+ip))/2;
                  tbox::plog << (*v)(coarse) << " ";
                }
              }
            else if(location_index==3)
              {
                // if(!cgeom->getTouchesRegularBoundary(1,1))
                {
                  pdat::SideIndex coarse(hier::Index(i,j),1,
                                         pdat::SideIndex::Lower);
                  pdat::SideIndex center(coarse*2);
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
