/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"
#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"

namespace Elastic
{
  void residual_2D
  (SAMRAI::pdat::SideData<double> &v,
   SAMRAI::pdat::CellData<double> &cell_moduli,
   SAMRAI::pdat::NodeData<double> &edge_moduli,
   SAMRAI::pdat::SideData<double> &v_rhs,
   SAMRAI::pdat::SideData<double> &v_resid,
   const SAMRAI::hier::Box &pbox,
   const double dxy[]);

  void residual_embedded_2D
  (SAMRAI::pdat::SideData<double> &v,
   SAMRAI::pdat::CellData<double> &cell_moduli,
   SAMRAI::pdat::NodeData<double> &edge_moduli,
   SAMRAI::pdat::SideData<double> &v_rhs,
   SAMRAI::pdat::SideData<double> &v_resid,
   SAMRAI::pdat::SideData<double> &level_set,
   const SAMRAI::hier::Box &pbox,
   const double dxy[]);

  void residual_3D
  (SAMRAI::pdat::SideData<double> &v,
   SAMRAI::pdat::CellData<double> &cell_moduli,
   SAMRAI::pdat::EdgeData<double> &edge_moduli,
   SAMRAI::pdat::SideData<double> &v_rhs,
   SAMRAI::pdat::SideData<double> &v_resid,
   const SAMRAI::hier::Box &pbox,
   const double dxyz[]);
}

void Elastic::Operators::computeCompositeResidualOnLevel
(SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& rhs,
 int level,
 bool error_equation_indicator)
{
  t_compute_composite_residual->start();

  if (residual.getPatchHierarchy() != solution.getPatchHierarchy()
      || rhs.getPatchHierarchy() != solution.getPatchHierarchy())
    { TBOX_ERROR(__FILE__ << ": residual, solution, and rhs hierarchies "
                 << "are not consistent."); }
  const SAMRAI::hier::PatchHierarchy &hierarchy=*residual.getPatchHierarchy();

  boost::shared_ptr<SAMRAI::hier::PatchLevel>
    patch_level = hierarchy.getPatchLevel(level);

  const int v_id = solution.getComponentDescriptorIndex(0);
  v_refine_patch_strategy.data_id=v_id;
  v_refine_patch_strategy.is_residual=error_equation_indicator;
  Coarse_Fine_Boundary_Refine::is_residual=error_equation_indicator;

  if (level > level_min)
    { ghostfill(v_id, level); }
  else
    { ghostfill_nocoarse(v_id, level); }

  set_physical_boundaries(v_id,hierarchy,level,error_equation_indicator);

  for (SAMRAI::hier::PatchLevel::Iterator pi(patch_level->begin());
       pi!=patch_level->end(); ++pi)
    {
      boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (solution.getComponentPatchData(0,*patch));
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch->getPatchData(cell_moduli_id));
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (rhs.getComponentPatchData(0,*patch));
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_resid_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (residual.getComponentPatchData(0,*patch));

      SAMRAI::hier::Box pbox=patch->getBox();
      pbox.growUpper(SAMRAI::hier::IntVector::getOne(dimension));
      boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (patch->getPatchGeometry());

      switch(dimension.getValue())
        {
        case 2:
          {
            boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_moduli_ptr =
              boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
              (patch->getPatchData(edge_moduli_id));
            
            if(!have_embedded_boundary())
              {
                residual_2D(*v_ptr,*cell_moduli_ptr,*edge_moduli_ptr,*v_rhs_ptr,
                            *v_resid_ptr,pbox,geom->getDx());
              }
            else
              {
                boost::shared_ptr<SAMRAI::pdat::SideData<double> >
                  level_set_ptr =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                  (patch->getPatchData(level_set_id));
                residual_embedded_2D(*v_ptr,*cell_moduli_ptr,*edge_moduli_ptr,
                                     *v_rhs_ptr,*v_resid_ptr,*level_set_ptr,
                                     pbox,geom->getDx());
              }
          }
          break;
        case 3:
          {
            boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_moduli_ptr =
              boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
              (patch->getPatchData(edge_moduli_id));
            residual_3D(*v_ptr,*cell_moduli_ptr,*edge_moduli_ptr,*v_rhs_ptr,
                        *v_resid_ptr,pbox,geom->getDx());
          }
          break;
        default:
          abort();
        }
    }

  t_compute_composite_residual->stop();
}

