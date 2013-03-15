#include "Elastic/V_Coarsen_Patch_Strategy.h"

void
Elastic::V_Coarsen_Patch_Strategy::postprocessCoarsen
(SAMRAI::hier::Patch &coarse,
 const SAMRAI::hier::Patch &fine,
 const SAMRAI::hier::Box &coarse_box,
 const SAMRAI::hier::IntVector& )
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(data_id));
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse.getPatchData(data_id));

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr;

  if(have_faults())
    {
      dv_mixed_ptr=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (fine.getPatchData(dv_mixed_id));
      dv_diagonal_ptr=boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (fine.getPatchData(dv_diagonal_id));
    }

  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (coarse.getPatchGeometry());

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr;
  if(have_embedded_boundary())
    level_set_ptr=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
      (fine.getPatchData(level_set_id));

  if(getDim().getValue()==2)
    {
      coarsen_2D(*v,*v_fine,dv_mixed_ptr,dv_diagonal_ptr,*coarse_geom,coarse_box);

      const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>
        &boundaries=coarse_fine[fine.getPatchLevelNumber()]
        ->getEdgeBoundaries(coarse.getGlobalId());

      fix_boundary_elements_2D(*v,*v_fine,dv_mixed_ptr,boundaries);
    }
  else
    {
      coarsen_3D(*v,*v_fine,dv_mixed_ptr,dv_diagonal_ptr,*coarse_geom,coarse_box);

      const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>
        &boundaries=coarse_fine[fine.getPatchLevelNumber()]
        ->getFaceBoundaries(coarse.getGlobalId());

      fix_boundary_elements_3D(*v,*v_fine,dv_mixed_ptr,boundaries);
    }
}
