#include "Elastic/V_Coarsen_Patch_Strategy.hxx"

void
Elastic::V_Coarsen_Patch_Strategy::postprocessCoarsen
(SAMRAI::hier::Patch &coarse_patch,
 const SAMRAI::hier::Patch &fine_patch,
 const SAMRAI::hier::Box &coarse_box,
 const SAMRAI::hier::IntVector& )
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine_patch.getPatchData(data_id));
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse_patch.getPatchData(data_id));

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr;

  if(have_faults())
    {
      dv_mixed_ptr=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (fine_patch.getPatchData(dv_mixed_id));
      dv_diagonal_ptr=boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (fine_patch.getPatchData(dv_diagonal_id));
    }

  if(fine_patch.getDim().getValue()==2)
    {
      coarsen_2D(*v,*v_fine,dv_mixed_ptr,dv_diagonal_ptr,
                 coarse_patch,fine_patch,coarse_box);

      const std::vector<SAMRAI::hier::BoundaryBox>
        &boundaries=coarse_fine[fine_patch.getPatchLevelNumber()]
        ->getEdgeBoundaries(coarse_patch.getGlobalId());

      fix_boundary_elements_2D(*v,*v_fine,dv_mixed_ptr,boundaries);
    }
  else
    {
      const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (coarse_patch.getPatchGeometry());

      coarsen_3D(*v,*v_fine,dv_mixed_ptr,dv_diagonal_ptr,*coarse_geom,coarse_box);

      const std::vector<SAMRAI::hier::BoundaryBox>
        &boundaries=coarse_fine[fine_patch.getPatchLevelNumber()]
        ->getFaceBoundaries(coarse_patch.getGlobalId());

      fix_boundary_elements_3D(*v,*v_fine,dv_mixed_ptr,boundaries);
    }
}
