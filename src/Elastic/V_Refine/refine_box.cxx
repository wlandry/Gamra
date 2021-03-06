/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/pdat/CellIterator.h>

#include "Elastic/V_Refine.hxx"

void Elastic::V_Refine::refine_box
(SAMRAI::hier::Patch &fine_patch,
 const SAMRAI::hier::Patch &coarse_patch,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::Box &fine_box,
 const SAMRAI::hier::IntVector &,
 const Gamra::Dir &ix) const
{
   const SAMRAI::tbox::Dimension &dimension(fine_patch.getDim());
   const Gamra::Dir dim(dimension.getValue());

   boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
     boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
     (coarse_patch.getPatchData(src_component));
   SAMRAI::pdat::SideData<double> &v(*v_ptr);
   boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine_ptr =
     boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
     (fine_patch.getPatchData(dst_component));
   SAMRAI::pdat::SideData<double> &v_fine(*v_fine_ptr);

   SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dimension)),
     jp(ip), kp(ip);
   ip[0]=1;
   jp[1]=1;
   if(dim>2)
     kp[2]=1;
   SAMRAI::hier::Index unit[]={ip,jp,kp};

   SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(fine_box));

   if(have_embedded_boundary())
     {
       boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_coarse_ptr =
         boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
         (coarse_patch.getPatchData(level_set_id));
       SAMRAI::pdat::SideData<double> &level_set_coarse(*level_set_coarse_ptr);

       boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_fine_ptr =
         boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
         (fine_patch.getPatchData(level_set_id));
       SAMRAI::pdat::SideData<double> &level_set_fine(*level_set_fine_ptr);
       
       boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
         boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
         (coarse_patch.getPatchGeometry());

       const SAMRAI::hier::Box coarse_box=level_set_coarse.getBox();
       const SAMRAI::hier::Box coarse_ghost_box=level_set_coarse.getGhostBox();

       /// The copy of the coarse values does not include the edges,
       /// so we have so manually set the level set on the edges
       SAMRAI::pdat::CellIterator
         coarse_end(SAMRAI::pdat::CellGeometry::end(coarse_ghost_box));
       for(SAMRAI::pdat::CellIterator
             ci(SAMRAI::pdat::CellGeometry::begin(coarse_ghost_box));
           ci!=coarse_end; ++ci)
         {
           SAMRAI::pdat::SideIndex
             coarse(*ci,ix,SAMRAI::pdat::SideIndex::Lower);
           for(Gamra::Dir d=ix.next(dim);d!=ix;d=d.next(dim))
             {
               if((coarse[d]<coarse_box.lower(d)
                   && geom->getTouchesRegularBoundary(d,0))
                 || (coarse[d]>coarse_box.upper(d)
                     && geom->getTouchesRegularBoundary(d,1)))
                 {
                   level_set_coarse(coarse)=-boundary_value;
                 }
             }
         }

       for(SAMRAI::pdat::CellIterator
             ci(SAMRAI::pdat::CellGeometry::begin(fine_box)); ci!=cend; ++ci)
         {
           SAMRAI::pdat::SideIndex
             fine(*ci,ix,SAMRAI::pdat::SideIndex::Lower);

           if(level_set_fine(fine)<0)
             {
               v_fine(fine)=invalid_value;
               continue;
             }

           SAMRAI::pdat::SideIndex coarse(fine);
           coarse.coarsen(SAMRAI::hier::Index::getOneIndex(dimension)*2);

           double temp(refine_along_line(v,ix,dim,unit,fine,coarse,
                                         level_set_coarse));

           if(fine[ix]%2==0)
             {
               v_fine(fine)=temp;
               if(temp==std::numeric_limits<double>::max())
                 TBOX_ERROR("Can not find a valid coarse point to refine from.\n"
                            << "  direction: " << ix 
                            << "\n  coarse: " << coarse
                            << "\n  fine: " << fine << "\n");
             }
           else
             {
               double temp2(refine_along_line(v,ix,dim,unit,fine,
                                              coarse+unit[ix],
                                              level_set_coarse));

               v_fine(fine)=(temp+temp2)/2;
               if(temp==std::numeric_limits<double>::max())
                 {
                   if(temp2==std::numeric_limits<double>::max())
                     {
                       TBOX_ERROR("Can not find a valid coarse point to "
                                  "refine from.\n"
                                  << "  direction: " << ix 
                                  << "\n  coarse: " << coarse
                                  << "\n  fine: " << fine << "\n");
                     }
                   else
                     {
                       v_fine(fine)=temp2;
                     }
                 }
               else if(temp2==std::numeric_limits<double>::max())
                 {
                   v_fine(fine)=temp;
                 }
             }
         }
     }
   else
     {
       SAMRAI::hier::Box coarse_box=coarse_patch.getBox();
       boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
         boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
         (coarse_patch.getPatchGeometry());

       for(SAMRAI::pdat::CellIterator
             ci(SAMRAI::pdat::CellGeometry::begin(fine_box)); ci!=cend; ++ci)
         {
           SAMRAI::pdat::SideIndex fine(*ci,ix,SAMRAI::pdat::SideIndex::Lower);

           SAMRAI::pdat::SideIndex coarse(fine);
           coarse.coarsen(SAMRAI::hier::Index::getOneIndex(dimension)*2);

           if(fine[ix]%2==0)
             {
               v_fine(fine)=refine_along_line(v,ix,dim,unit,fine,coarse,
                                              coarse_box,*coarse_geom);
             }
           else
             {
               v_fine(fine)=(refine_along_line(v,ix,dim,unit,fine,coarse,
                                               coarse_box,*coarse_geom)
                             + refine_along_line(v,ix,dim,unit,fine,
                                                 coarse+unit[ix],
                                                 coarse_box,*coarse_geom))/2;
             }
         }
     }
}
