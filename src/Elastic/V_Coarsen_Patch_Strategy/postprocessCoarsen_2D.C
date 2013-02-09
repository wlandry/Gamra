#include "Elastic/V_Coarsen_Patch_Strategy.h"

void
Elastic::V_Coarsen_Patch_Strategy::postprocessCoarsen_2D
(SAMRAI::hier::Patch& coarse,
 const SAMRAI::hier::Patch& fine,
 const SAMRAI::hier::Box& ,
 const SAMRAI::hier::IntVector& )
{
  /* FIXME: Why is this required?  Shouldn't the boundary points be
     ok?  They are partially interpolated from the coarse level, but
     that should not be a problem.  Unfortunately, it is a problem,
     because removing this routine generates nan's. */

  /* Fix up the boundary elements by iterating through the boundary
     boxes */

  /* We only care about edges, not corners, so we only iterate over
     edge boundary boxes. */
  const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>
    &boundaries=coarse_fine[fine.getPatchLevelNumber()]->getEdgeBoundaries(coarse.getGlobalId());
     
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(v_id));
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse.getPatchData(v_id));

  const bool is_residual(v_id!=7);
  const int dv_mixed_id(5);
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(dv_mixed_id));

  TBOX_ASSERT(v);
  TBOX_ASSERT(v_fine);
  TBOX_ASSERT(v_fine->getDepth() == v->getDepth());
  TBOX_ASSERT(v->getDepth() == 1);

  SAMRAI::hier::Box gbox(v_fine->getGhostBox());
  SAMRAI::hier::Index ip(1,0), jp(0,1);
  for(int mm=0; mm<boundaries.size(); ++mm)
    {
      SAMRAI::hier::Box bbox=boundaries[mm].getBox();
      int location_index=boundaries[mm].getLocationIndex();

      SAMRAI::hier::Index
        lower=SAMRAI::hier::Index::coarsen(bbox.lower(),
                                           SAMRAI::hier::Index(2,2)),
        upper=SAMRAI::hier::Index::coarsen(bbox.upper(),
                                           SAMRAI::hier::Index(2,2));

      for(int j=lower(1); j<=upper(1); ++j)
        for(int i=lower(0); i<=upper(0); ++i)
          {
            /* Fix vx */
            if(location_index==0)
              {
                SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j),0,
                                               SAMRAI::pdat::SideIndex::Upper);
                SAMRAI::pdat::SideIndex x(coarse*2);
                if(x[1]>=gbox.lower(1) && x[1]<gbox.upper(1))
                  {
                    (*v)(coarse)=((*v_fine)(x) + (*v_fine)(x+jp))/2;
                    if(!is_residual)
                      (*v)(coarse)+=((*dv_mixed)(x,0) + (*dv_mixed)(x+jp,1))/2;
                  }
              }
            else if(location_index==1)
              {
                SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j),0,
                                               SAMRAI::pdat::SideIndex::Lower);
                SAMRAI::pdat::SideIndex x(coarse*2);
                if(x[1]>=gbox.lower(1) && x[1]<gbox.upper(1))
                  {
                    (*v)(coarse)=((*v_fine)(x) + (*v_fine)(x+jp))/2;
                    if(!is_residual)
                      (*v)(coarse)+=((*dv_mixed)(x,0) + (*dv_mixed)(x+jp,1))/2;
                  }
              }
            /* Fix vy */
            else if(location_index==2)
              {
                SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j),1,
                                               SAMRAI::pdat::SideIndex::Upper);
                SAMRAI::pdat::SideIndex y(coarse*2);
                if(y[0]>=gbox.lower(0) && y[0]<gbox.upper(0))
                  {
                    (*v)(coarse)=((*v_fine)(y) + (*v_fine)(y+ip))/2;
                    if(!is_residual)
                      (*v)(coarse)+=((*dv_mixed)(y,0) + (*dv_mixed)(y+ip,1))/2;
                  }
              }
            else if(location_index==3)
              {
                SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j),1,
                                               SAMRAI::pdat::SideIndex::Lower);
                SAMRAI::pdat::SideIndex y(coarse*2);
                if(y[0]>=gbox.lower(0) && y[0]<gbox.upper(0))
                  {
                    (*v)(coarse)=((*v_fine)(y) + (*v_fine)(y+ip))/2;
                    if(!is_residual)
                      (*v)(coarse)+=((*dv_mixed)(y,0) + (*dv_mixed)(y+ip,1))/2;
                  }
              }
            else
              {
                abort();
              }
          }
    }
}
