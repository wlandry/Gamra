/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"
#include "Elastic/V_Boundary_Refine.hxx"

void Elastic::FACOps::computeCompositeResidualOnLevel
(SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& rhs,
 int ln,
 bool error_equation_indicator)
{
  t_compute_composite_residual->start();

  if (residual.getPatchHierarchy() != solution.getPatchHierarchy()
      || rhs.getPatchHierarchy() != solution.getPatchHierarchy())
    TBOX_ERROR(d_object_name << ": residual, solution, and rhs hierarchies "
               << "are not consistent.");
  const SAMRAI::hier::PatchHierarchy &hierarchy=*residual.getPatchHierarchy();

  boost::shared_ptr<SAMRAI::hier::PatchLevel>
    level = hierarchy.getPatchLevel(ln);

  const int v_id = solution.getComponentDescriptorIndex(0);
  v_refine_patch_strategy.data_id=v_id;
  v_refine_patch_strategy.is_residual=error_equation_indicator;
  V_Boundary_Refine::is_residual=error_equation_indicator;

  if (ln > d_ln_min)
    xeqScheduleGhostFill(v_id, ln);
  else
    xeqScheduleGhostFillNoCoarse(v_id, ln);

  set_boundaries(v_id,hierarchy,ln,error_equation_indicator);

  for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
       pi!=level->end(); ++pi)
    {
      boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (solution.getComponentPatchData(0,(*patch)));
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
      pbox.growUpper(SAMRAI::hier::IntVector::getOne(d_dim));
      boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (patch->getPatchGeometry());

      switch(d_dim.getValue())
        {
        case 2:
          residual_2D(*v_ptr,*cell_moduli_ptr,*v_rhs_ptr,
                      *v_resid_ptr,*patch,pbox,*geom);
          break;
        case 3:
          residual_3D(*v_ptr,*cell_moduli_ptr,*v_rhs_ptr,
                      *v_resid_ptr,*patch,pbox,*geom);
          break;
        default:
          abort();
        }
    }

  t_compute_composite_residual->stop();
}

