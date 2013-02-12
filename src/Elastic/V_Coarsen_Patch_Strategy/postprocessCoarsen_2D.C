#include "Elastic/V_Coarsen_Patch_Strategy.h"

void
Elastic::V_Coarsen_Patch_Strategy::postprocessCoarsen_2D
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

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(dv_mixed_id));

  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal =
    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
    (fine.getPatchData(dv_diagonal_id));

  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (coarse.getPatchGeometry());

  coarsen_2D(*v,*v_fine,*dv_mixed,*dv_diagonal,*coarse_geom,coarse_box);

  const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>
    &boundaries=coarse_fine[fine.getPatchLevelNumber()]
    ->getEdgeBoundaries(coarse.getGlobalId());

  fix_boundary_elements_2D(*v,*v_fine,*dv_mixed,boundaries);
}
