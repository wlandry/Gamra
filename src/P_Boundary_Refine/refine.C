#include "P_Boundary_Refine.h"
#include "set_boundary.h"
#include "Constants.h"

void SAMRAI::geom::P_Boundary_Refine::refine(hier::Patch& fine,
                                             const hier::Patch& coarse,
                                             const int dst_component,
                                             const int src_component,
                                             const hier::BoxOverlap& overlap,
                                             const hier::IntVector& ratio)
  const
{
  const pdat::CellOverlap* t_overlap =
    dynamic_cast<const pdat::CellOverlap *>(&overlap);
  const hier::BoxList& boxes = t_overlap->getDestinationBoxList();
  const tbox::Dimension& dimension(getDim());
  const int dim(dimension.getValue());

  set_boundary(coarse,src_component,invalid_id,true);

  for (hier::BoxList::Iterator b(boxes); b; b++)
    {
      hier::Box &overlap_box=b();
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension,fine,coarse,overlap_box,ratio);

      tbox::Pointer<pdat::CellData<double> >
        p = coarse.getPatchData(src_component);
      tbox::Pointer<pdat::CellData<double> >
        p_fine = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!p.isNull());
      TBOX_ASSERT(!p_fine.isNull());
      TBOX_ASSERT(p->getDepth() == p_fine->getDepth());
      TBOX_ASSERT(p->getDepth() == 1);
#endif

      hier::Box fine_box=fine.getBox();
      hier::Box gbox=p_fine->getGhostBox();

      /* We have to infer where the boundary is from the boxes */
      int boundary_direction;
      bool boundary_positive;

      for(int d=0;d<dim;++d)
        {
          if(overlap_box.lower(d)-overlap_box.upper(d)==0)
            {
              boundary_direction=d;
              if(fine_box.upper(d)<overlap_box.lower(d))
                boundary_positive=true;
              else if(fine_box.lower(d)>overlap_box.upper(d))
                boundary_positive=false;
              else
                abort();
              break;
            }
        }

      hier::Index p_min(dimension), p_max(dimension);
      for(int d=0;d<dim;++d)
        {
          p_min[d]=std::max(overlap_box.lower(d),gbox.lower(d));
          p_max[d]=std::min(overlap_box.upper(d),gbox.upper(d));
        }

      hier::Index ip(hier::Index::getZeroIndex(dimension)), jp(ip), kp(ip);
      ip[0]=1;
      jp[1]=1;
      if(dim>2)
        kp[2]=1;

      if(dim==2)
        {
          /* This odd stride is because we handle all of the fine
           * boundary cells in a coarse cell at once.  However,
           * sometimes there is only one fine cell in a coarse cell,
           * so the starting point does not align with the coarse
           * cell.  The stride ensures that, if we start not aligned,
           * the next step will be aligned. */

          for(int j=p_min[1]; j<=p_max[1]; j=(j/2)*2+2)
            for(int i=p_min[0]; i<=p_max[0]; i=(i/2)*2+2)
              {
                pdat::CellIndex fine(hier::Index(i,j));
        
                switch(boundary_direction)
                  {
                  case 0:
                    Update_P_2D(fine,boundary_positive ? ip : -ip,jp,j,p_max[1],
                                *p,*p_fine);
                    break;
                  case 1:
                    Update_P_2D(fine,boundary_positive ? jp : -jp,ip,i,p_max[0],
                                *p,*p_fine);
                    break;
                  }
              }
        }
      else
        {
          for(int k=p_min[2]; k<=p_max[2]; ++k)
            for(int j=p_min[1]; j<=p_max[1]; ++j)
              for(int i=p_min[0]; i<=p_max[0]; ++i)
                {
                  pdat::CellIndex fine(hier::Index(i,j,k));
        
                  switch(boundary_direction)
                    {
                    case 0:
                      Update_P_3D(fine,boundary_positive ? ip : -ip,jp,kp,
                                  j,k,p_max[1],p_max[2],*p,*p_fine);
                      break;
                    case 1:
                      Update_P_3D(fine,boundary_positive ? jp : -jp,kp,ip,
                                  k,i,p_max[2],p_max[0],*p,*p_fine);
                      break;
                    case 2:
                      Update_P_3D(fine,boundary_positive ? kp : -kp,ip,jp,
                                  i,j,p_max[0],p_max[1],*p,*p_fine);
                      break;
                    }
                }
        }
    }
}
