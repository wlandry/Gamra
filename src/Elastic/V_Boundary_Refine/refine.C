/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#include "Elastic/V_Boundary_Refine.h"
#include "Constants.h"

void Elastic::V_Boundary_Refine::refine
(SAMRAI::hier::Patch& fine,
 const SAMRAI::hier::Patch& coarse,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::BoxOverlap& fine_overlap,
 const SAMRAI::hier::IntVector& ratio) const
{
   const SAMRAI::pdat::SideOverlap* t_overlap =
      dynamic_cast<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   /* Commenting this out since it should be set by the refine patch
      strategy */
   // boundary_conditions.set_boundary(coarse,src_component,true);

   for(int axis=0; axis<getDim().getValue(); ++axis)
     {
       const SAMRAI::hier::BoxContainer&
         boxes = t_overlap->getDestinationBoxContainer(axis);
       for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin());
            b!=boxes.end(); b++)
         {
           refine(fine,coarse,dst_component,src_component,*b,ratio,axis);
         }
     }
}

void Elastic::V_Boundary_Refine::refine
(SAMRAI::hier::Patch &fine,
 const SAMRAI::hier::Patch &coarse,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::Box &overlap_box,
 const SAMRAI::hier::IntVector &,
 const int &axis) const
{
  const SAMRAI::tbox::Dimension &dimension(getDim());
  const int dim(dimension.getValue());

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse.getPatchData(src_component));
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(dst_component));

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal;
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_fine;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_fine;

  if(!is_residual)
    {
      dv_mixed=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (coarse.getPatchData(dv_mixed_id));
      dv_diagonal=boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (coarse.getPatchData(dv_diagonal_id));
      dv_mixed_fine=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (fine.getPatchData(dv_mixed_id));
      dv_diagonal_fine=
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (fine.getPatchData(dv_diagonal_id));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(v);
   TBOX_ASSERT(v_fine);
   TBOX_ASSERT(v->getDepth() == v_fine->getDepth());
   TBOX_ASSERT(v->getDepth() == 1);
#endif

   SAMRAI::hier::Box fine_box=fine.getBox();

   /* We have to infer where the boundary is from the boxes */
   int boundary_direction;
   bool boundary_positive(false);

   for(int d=0;d<dim;++d)
     {
       if(std::abs(overlap_box.lower(d)-overlap_box.upper(d))==(axis==d ? 1 : 0))
         {
           boundary_direction=d;
           if(fine_box.upper(d)<=overlap_box.lower(d))
             boundary_positive=true;
           else if(fine_box.lower(d)>=overlap_box.upper(d))
             boundary_positive=false;
           else
             abort();
           break;
         }
     }

   SAMRAI::hier::Index p_min(overlap_box.lower()), p_max(overlap_box.upper());

   if(boundary_direction==axis)
     {
       if(boundary_positive)
         {
           p_min[axis]=p_max[axis];
         }
       else
         {
           p_max[axis]=p_min[axis];
         }
     }

   SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dimension)),
     jp(ip), kp(ip);
   ip[0]=1;
   jp[1]=1;
   if(dim>2)
     kp[2]=1;

   if(dim==2)
     {
       for(int j=p_min[1]; j<=p_max[1]; ++j)
         for(int i=p_min[0]; i<=p_max[0]; ++i)
           {
             SAMRAI::pdat::SideIndex fine_index(SAMRAI::hier::Index(i,j),axis,
                                                SAMRAI::pdat::SideIndex::Lower);
             switch(axis)
               {
               case 0:
                 Update_V_2D(axis,boundary_direction,boundary_positive,
                             fine_index,ip,jp,i,j,*v,*v_fine);
                 if(!is_residual)
                   Correction_2D(axis,boundary_direction,boundary_positive,
                                 fine_index,ip,jp,i,j,p_min[0],p_max[0],*dv_diagonal,
                                 *dv_diagonal_fine,*dv_mixed,*dv_mixed_fine,
                                 *v_fine);
                 break;
               case 1:
                 Update_V_2D(axis,boundary_direction,boundary_positive,
                             fine_index,jp,ip,j,i,*v,*v_fine);

                 if(!is_residual)
                   Correction_2D(axis,boundary_direction,boundary_positive,
                                 fine_index,jp,ip,j,i,p_min[1],p_max[1],*dv_diagonal,
                                 *dv_diagonal_fine,*dv_mixed,*dv_mixed_fine,
                                 *v_fine);
                 break;
               default:
                 abort();
                 break;
               }
         }
     }
   else
     {
       SAMRAI::hier::Index unit[]={ip,jp,kp};
       SAMRAI::hier::Index ijk(dimension);
       SAMRAI::hier::Box coarse_box(coarse.getBox());
       boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
         boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
         (coarse.getPatchGeometry());

       for(ijk[2]=p_min[2]; ijk[2]<=p_max[2]; ijk[2]+=1)
         for(ijk[1]=p_min[1]; ijk[1]<=p_max[1]; ijk[1]+=1)
           for(ijk[0]=p_min[0]; ijk[0]<=p_max[0]; ijk[0]+=1)
             {
               SAMRAI::pdat::SideIndex
                 fine(ijk,axis,SAMRAI::pdat::SideIndex::Lower);
               Update_V_3D(axis,boundary_direction,boundary_positive,fine,
                           unit,ijk,coarse_box,*geom,*v,*v_fine);
             }
     }
}
