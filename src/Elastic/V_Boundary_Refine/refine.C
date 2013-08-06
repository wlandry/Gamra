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
(SAMRAI::hier::Patch& fine_patch,
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

   for(int ix=0; ix<fine_patch.getDim().getValue(); ++ix)
     {
       const SAMRAI::hier::BoxContainer&
         boxes = t_overlap->getDestinationBoxContainer(ix);
       for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin());
            b!=boxes.end(); ++b)
         {
           refine(fine_patch,coarse,dst_component,src_component,*b,ratio,ix);
         }
     }
}

void Elastic::V_Boundary_Refine::refine
(SAMRAI::hier::Patch &fine_patch,
 const SAMRAI::hier::Patch &coarse_patch,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::Box &overlap_box,
 const SAMRAI::hier::IntVector &,
 const int &ix) const
{
  const SAMRAI::tbox::Dimension &dimension(fine_patch.getDim());
  const int dim(dimension.getValue());

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse_patch.getPatchData(src_component));
  SAMRAI::pdat::SideData<double> &v(*v_ptr);
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine_ptr =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine_patch.getPatchData(dst_component));
  SAMRAI::pdat::SideData<double> &v_fine(*v_fine_ptr);

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr;
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_fine_ptr;

  if(have_embedded_boundary())
    {
      level_set_ptr=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (coarse_patch.getPatchData(level_set_id));
      level_set_fine_ptr=
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (fine_patch.getPatchData(level_set_id));
    }

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal;
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_fine;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_fine;

  if(have_faults() && !is_residual)
    {
      dv_mixed=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (coarse_patch.getPatchData(dv_mixed_id));
      dv_diagonal=boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (coarse_patch.getPatchData(dv_diagonal_id));
      dv_mixed_fine=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (fine_patch.getPatchData(dv_mixed_id));
      dv_diagonal_fine=
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (fine_patch.getPatchData(dv_diagonal_id));
    }

   SAMRAI::hier::Box fine_box=fine_patch.getBox();

   /* We have to infer where the boundary is from the boxes */
   int boundary_direction;
   bool boundary_positive(false);

   for(int d=0;d<dim;++d)
     {
       if(std::abs(overlap_box.lower(d)-overlap_box.upper(d))==(ix==d ? 1 : 0))
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

   SAMRAI::hier::Index fine_min(overlap_box.lower()),
     fine_max(overlap_box.upper());

   if(boundary_direction==ix)
     {
       if(boundary_positive)
         {
           fine_min[ix]=fine_max[ix];
         }
       else
         {
           fine_max[ix]=fine_min[ix];
         }
     }

   SAMRAI::hier::Index unit[]={SAMRAI::hier::Index::getZeroIndex(dimension),
                               SAMRAI::hier::Index::getZeroIndex(dimension),
                               SAMRAI::hier::Index::getZeroIndex(dimension)};
   for(int d=0;d<dim;++d)
     unit[d][d]=1;
   SAMRAI::hier::Index ijk(dimension);

   if(dim==2)
     {
       const int iy((ix+1)%dim);
       if(!have_embedded_boundary())
         {
           for(ijk[1]=fine_min[1]; ijk[1]<=fine_max[1]; ijk[1]+=1)
             for(ijk[0]=fine_min[0]; ijk[0]<=fine_max[0]; ijk[0]+=1)
               {
                 SAMRAI::pdat::SideIndex
                   fine_index(ijk,ix,SAMRAI::pdat::SideIndex::Lower);

                 Update_V_2D(ix,boundary_direction,boundary_positive,fine_index,
                             unit[ix],unit[iy],ijk[ix],ijk[iy],v,v_fine);
                 if(have_faults() && !is_residual)
                   Correction_2D(ix,boundary_direction,boundary_positive,
                                 fine_index,unit[ix],unit[iy],ijk[ix],ijk[iy],
                                 fine_min[ix],fine_max[ix],
                                 *dv_diagonal,*dv_diagonal_fine,*dv_mixed,
                                 *dv_mixed_fine,v_fine);
               }
         }
       else
         {
           SAMRAI::pdat::SideData<double> &level_set(*level_set_ptr);
           SAMRAI::pdat::SideData<double> &level_set_fine(*level_set_fine_ptr);

           for(ijk[1]=fine_min[1]; ijk[1]<=fine_max[1]; ijk[1]+=1)
             for(ijk[0]=fine_min[0]; ijk[0]<=fine_max[0]; ijk[0]+=1)
               {
                 SAMRAI::pdat::SideIndex
                   fine_index(ijk,ix,SAMRAI::pdat::SideIndex::Lower);

                 if(level_set_fine(fine_index)>=0)
                   {
                     Update_V_embedded_2D(ix,boundary_direction,
                                          boundary_positive,
                                          fine_index,unit[ix],unit[iy],ijk[ix],
                                          ijk[iy],level_set,level_set_fine,v,v_fine);
                     if(have_faults() && !is_residual)
                       Correction_2D(ix,boundary_direction,boundary_positive,
                                     fine_index,unit[ix],unit[iy],ijk[ix],ijk[iy],
                                     fine_min[ix],fine_max[ix],
                                     *dv_diagonal,*dv_diagonal_fine,*dv_mixed,
                                     *dv_mixed_fine,v_fine);
                   }
               }
         }
     }
   else
     {
       SAMRAI::hier::Box coarse_box(coarse_patch.getBox());
       boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
         boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
         (coarse_patch.getPatchGeometry());

       for(ijk[2]=fine_min[2]; ijk[2]<=fine_max[2]; ijk[2]+=1)
         for(ijk[1]=fine_min[1]; ijk[1]<=fine_max[1]; ijk[1]+=1)
           for(ijk[0]=fine_min[0]; ijk[0]<=fine_max[0]; ijk[0]+=1)
             {
               SAMRAI::pdat::SideIndex
                 fine(ijk,ix,SAMRAI::pdat::SideIndex::Lower);
               Update_V_3D(ix,boundary_direction,boundary_positive,fine,
                           unit,ijk,coarse_box,fine_min,fine_max,*coarse_geom,
                           dv_diagonal,dv_diagonal_fine,dv_mixed,
                           dv_mixed_fine,v,v_fine);
             }
     }
}
