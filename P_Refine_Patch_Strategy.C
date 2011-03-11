/* Extrapolate pressure onto ghost cells at physical boundaries.
   Everywhere else the pressure should already be set correctly. */ 

#include "P_Refine_Patch_Strategy.h"

void
SAMRAI::solv::P_Refine_Patch_Strategy::preprocessRefine
(hier::Patch& ,
 const hier::Patch& coarse,
 const hier::Box& ,
 const hier::IntVector& )
{
  tbox::Pointer<pdat::CellData<double> >
    p = coarse.getPatchData(p_id);

  hier::Box pbox=coarse.getBox();
  hier::Box gbox=p->getGhostBox();

  tbox::Pointer<geom::CartesianPatchGeometry>
    geom = coarse.getPatchGeometry();
  for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
    for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
      {
        pdat::CellIndex center(tbox::Dimension(2));
        center[0]=i;
        center[1]=j;
        hier::Index ip(1,0), jp(0,1);

        /* This is a bit complicated because we have to deal with
           corners.  I have a feeling that corners are impossible on
           these boundary boxes, but this is just to be safe. */

        /* x=0 */
        if(i<pbox.lower(0) && geom->getTouchesRegularBoundary(0,0))
          {
            /* y=0 */
            if(j<pbox.lower(1) && geom->getTouchesRegularBoundary(1,0))
              {
                (*p)(center)=4*(*p)(center+ip+jp) + (*p)(center+ip*2+jp*2)
                  - 2*((*p)(center+ip*2+jp) + (*p)(center+ip+jp*2));
              }
            /* y=y_max */
            else if(j>pbox.upper(1) && geom->getTouchesRegularBoundary(1,1))
              {
                (*p)(center)=4*(*p)(center+ip-jp) + (*p)(center+ip*2-jp*2)
                  - 2*((*p)(center+ip*2-jp) + (*p)(center+ip-jp*2));
              }
            else
              {
                (*p)(center)=2*(*p)(center+ip) - (*p)(center+ip*2);
              }
          }
        /* x=x_max */
        else if(i>pbox.upper(0) && geom->getTouchesRegularBoundary(0,1))
          {
            /* y=0 */
            if(j<pbox.lower(1) && geom->getTouchesRegularBoundary(1,0))
              {
                (*p)(center)=4*(*p)(center-ip+jp) + (*p)(center-ip*2+jp*2)
                  - 2*((*p)(center-ip*2+jp) + (*p)(center-ip+jp*2));
              }
            /* y=y_max */
            else if(j>pbox.upper(1) && geom->getTouchesRegularBoundary(1,1))
              {
                (*p)(center)=4*(*p)(center-ip-jp) + (*p)(center-ip*2-jp*2)
                  - 2*((*p)(center-ip*2-jp) + (*p)(center-ip-jp*2));
              }
            else
              {
                (*p)(center)=2*(*p)(center-ip) - (*p)(center-ip*2);
              }
          }
        /* These boundaries are simpler because we have already
           handled the corners */
        /* y=0 */
        else if(j<pbox.lower(1) && geom->getTouchesRegularBoundary(1,0))
          {
            (*p)(center)=2*(*p)(center+jp) - (*p)(center+jp*2);
          }
        /* y=y_max */
        else if(j>pbox.upper(1) && geom->getTouchesRegularBoundary(1,1))
          {
            (*p)(center)=2*(*p)(center-jp) - (*p)(center-jp*2);
          }
      }
}
