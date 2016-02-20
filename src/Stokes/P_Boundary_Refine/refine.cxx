#include "Stokes/P_Boundary_Refine.hxx"
#include "Stokes/set_boundary.hxx"
#include "Constants.hxx"

void Stokes::P_Boundary_Refine::refine
(SAMRAI::hier::Patch& fine,
 const SAMRAI::hier::Patch& coarse,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::BoxOverlap& overlap,
 const SAMRAI::hier::IntVector&) const
{
  const SAMRAI::pdat::CellOverlap* t_overlap =
    dynamic_cast<const SAMRAI::pdat::CellOverlap *>(&overlap);
  const SAMRAI::hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
  const SAMRAI::tbox::Dimension& dimension(fine.getDim());
  const int dim(dimension.getValue());

  Stokes_set_boundary(coarse,src_component,invalid_id,true);

  for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin()); b!=boxes.end(); ++b)
    {
      const SAMRAI::hier::Box &overlap_box=*b;

      boost::shared_ptr<SAMRAI::pdat::CellData<double> > p =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (coarse.getPatchData(src_component));
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_fine =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (fine.getPatchData(dst_component));
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(p);
      TBOX_ASSERT(p_fine);
      TBOX_ASSERT(p->getDepth() == p_fine->getDepth());
      TBOX_ASSERT(p->getDepth() == 1);
#endif

      SAMRAI::hier::Box fine_box=fine.getBox();
      SAMRAI::hier::Box gbox=p_fine->getGhostBox();

      /* We have to infer where the boundary is from the boxes */
      int boundary_direction(-1);
      bool boundary_positive(true);

      for(Gamra::Dir d=0;d<dim;++d)
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

      SAMRAI::hier::Index p_min(dimension), p_max(dimension);
      for(Gamra::Dir d=0;d<dim;++d)
        {
          p_min[d]=std::max(overlap_box.lower(d),gbox.lower(d));
          p_max[d]=std::min(overlap_box.upper(d),gbox.upper(d));
        }

      SAMRAI::hier::IntVector ip(SAMRAI::hier::IntVector::getZero(dimension)), jp(ip), kp(ip);
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
                SAMRAI::pdat::CellIndex fine(SAMRAI::hier::Index(i,j));
        
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
                  SAMRAI::pdat::CellIndex fine(SAMRAI::hier::Index(i,j,k));
        
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
