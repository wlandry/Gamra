#include "V_Refine_Patch_Strategy.h"
#include "set_V_boundary.h"

void
SAMRAI::solv::V_Refine_Patch_Strategy::preprocessRefine
(hier::Patch& fine,
 const hier::Patch& coarse,
 const hier::Box& fine_box,
 const hier::IntVector& ratio)
{
  // tbox::Pointer<pdat::SideData<double> >
  //   v = coarse.getPatchData(v_id);

  // hier::Box pbox=coarse.getBox();
  // hier::Box gbox=v->getGhostBox();

  // tbox::Pointer<geom::CartesianPatchGeometry>
  //   geom = coarse->getPatchGeometry();
  // tbox::plog << "V Refine Patch preprocess "
  //            << gbox.lower(0) << " "
  //            << gbox.upper(0) << " "
  //            << gbox.lower(1) << " "
  //            << gbox.upper(1) << " "
  //            << pbox.lower(0) << " "
  //            << pbox.upper(0) << " "
  //            << pbox.lower(1) << " "
  //            << pbox.upper(1) << " "
  //            << "\n";

  set_V_boundary(coarse,v_id,true);


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
  //           if((i<pbox.lower(0) && geom->getTouchesPhysicalBoundary(0,0))
  //              || (i>pbox.upper(0)+1 && geom->getTouchesPhysicalBoundary(0,1)))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                    pdat::SideIndex::Lower))=
  //                 boundary_value;
  //             }
  //           /* Set the value so the derivative=0 */
  //           else if(j<pbox.lower(1) && geom->getTouchesPhysicalBoundary(1,0))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center+jp,pdat::SideIndex::X,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           else if(j>pbox.upper(1) && geom->getTouchesPhysicalBoundary(1,1))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center-jp,pdat::SideIndex::X,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           // tbox::plog << "V Refine Patch vx "
  //           //            << coarse.getPatchLevelNumber() << " "
  //           //            << i << " "
  //           //            << j << " "
  //           //            << (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //           //                                    pdat::SideIndex::Lower))
  //           //            << " "
  //           //            << (&(*v)(pdat::SideIndex(center,pdat::SideIndex::X,
  //           //                                      pdat::SideIndex::Lower)))
  //           //            << " "
  //           //            << "\n";
  //         }
  //       /* vy */
  //       if(i<=gbox.upper(0))
  //         {
  //           if((j<pbox.lower(1) && geom->getTouchesPhysicalBoundary(1,0))
  //              || (j>pbox.upper(1)+1 && geom->getTouchesPhysicalBoundary(1,1)))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                    pdat::SideIndex::Lower))=
  //                 boundary_value;
  //             }
  //           else if(i<pbox.lower(0) && geom->getTouchesPhysicalBoundary(0,0))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center+ip,pdat::SideIndex::Y,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           else if(i>pbox.upper(0) && geom->getTouchesPhysicalBoundary(0,1))
  //             {
  //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //                                    pdat::SideIndex::Lower))=
  //                 (*v)(pdat::SideIndex(center-ip,pdat::SideIndex::Y,
  //                                      pdat::SideIndex::Lower));
  //             }
  //           // tbox::plog << "V Refine Patch vy "
  //           //            << coarse.getPatchLevelNumber() << " "
  //           //            << i << " "
  //           //            << j << " "
  //           //            << (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //           //                                    pdat::SideIndex::Lower))
  //           //            << " "
  //           //            << (&(*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
  //           //                                      pdat::SideIndex::Lower)))
  //           //            << " "
  //           //            << "\n";
  //         }
  //     }
}
