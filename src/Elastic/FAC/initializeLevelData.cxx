/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"
#include "Constants.hxx"

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>

/// Initialize data on a level, including allocating memory.  It does
/// not set the terms due to the faults, since the moduli have to be
/// fixed first.

void Elastic::FAC::initializeLevelData
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &patch_hierarchy,
 const int level_number,
 const double ,
 const bool ,
 const bool ,
 const boost::shared_ptr<SAMRAI::hier::PatchLevel>& ,
 const bool allocate_data)
{
  boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
    hierarchy = patch_hierarchy;
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (hierarchy->getGridGeometry());

  SAMRAI::hier::PatchLevel &level = *hierarchy->getPatchLevel(level_number);
  const int dim=dimension.getValue();

  if(allocate_data)
    {
      level.allocatePatchData(cell_moduli_id);
      level.allocatePatchData(edge_moduli_id);
      level.allocatePatchData(v_id);
      level.allocatePatchData(v_rhs_id);
      if(!faults.empty())
        {
          level.allocatePatchData(dv_diagonal_id);
          level.allocatePatchData(dv_mixed_id);
        }
      if(have_embedded_boundary())
        {
          level.allocatePatchData(level_set_id);
        }
    }
  /// Initialize data in all patches in the level.
  for (SAMRAI::hier::PatchLevel::Iterator p(level.begin()); p!=level.end();
       ++p)
    {
      SAMRAI::hier::Patch &patch = **p;
      boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (patch.getPatchGeometry());
      const double *dx=geom->getDx();

      /// cell moduli
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch.getPatchData(cell_moduli_id));

      SAMRAI::hier::Box cell_moduli_box = cell_moduli->getBox();

      SAMRAI::pdat::CellIterator
        cend(SAMRAI::pdat::CellGeometry::end(cell_moduli->getGhostBox()));
      for(SAMRAI::pdat::CellIterator
            ci(SAMRAI::pdat::CellGeometry::begin(cell_moduli->getGhostBox()));
          ci!=cend; ++ci)
        {
          const SAMRAI::pdat::CellIndex &c(*ci);
          double coord[3];
          for(int d=0;d<dim;++d)
            coord[d]=geom->getXLower()[d]
              + dx[d]*(c[d]-cell_moduli_box.lower()[d] + 0.5);

          (*cell_moduli)(c,0)=lambda.eval(coord);
          (*cell_moduli)(c,1)=mu.eval(coord);
        }

      /// v_rhs
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_data =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(v_rhs_id));

      for(Gamra::Dir ix=0; ix<dim; ++ix)
        {
          double offset[]={0.5,0.5,0.5};
          offset[ix]=0;

          SAMRAI::hier::Box v_rhs_box = v_rhs_data->getBox();
          SAMRAI::pdat::SideIterator
            end(SAMRAI::pdat::SideGeometry::end(v_rhs_box,ix));
          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(v_rhs_box,ix));
              si!=end; ++si)
            {
              const SAMRAI::pdat::SideIndex &x(*si);

              std::vector<double> coord(dim);
              for(int d=0;d<dim;++d)
                { coord[d]=geom->getXLower()[d]
                    + dx[d]*(x[d]-v_rhs_box.lower()[d]+offset[d]); }
              (*v_rhs_data)(x)=v_rhs[ix].eval(coord.data());
            }
        }

      /// Embedded boundaries
      if(have_embedded_boundary())
        {
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch.getPatchData(level_set_id));
          SAMRAI::hier::Box level_set_box = level_set_ptr->getBox();

          double dx_max(0);
          for(int d=0;d<dim;++d)
            { dx_max=std::max(dx_max,dx[d]); }
          double epsilon((grid_geom->getXUpper()[0]-grid_geom->getXLower()[0])
                         /1e10);
          for(Gamra::Dir ix=0;ix<dim;++ix)
            {
              double offset[]={0.5,0.5,0.5};
              offset[ix]=0;
              SAMRAI::pdat::SideIterator
                end(SAMRAI::pdat::SideGeometry::end
                    (level_set_ptr->getGhostBox(),ix));
              for(SAMRAI::pdat::SideIterator
                    si(SAMRAI::pdat::SideGeometry::begin
                       (level_set_ptr->getGhostBox(),ix));
                  si!=end; ++si)
                {
                  const SAMRAI::pdat::SideIndex &x(*si);
                  double coord[3];
                  for(int d=0;d<dim;++d)
                    { coord[d]=geom->getXLower()[d]
                        + dx[d]*(x[d]-level_set_box.lower()[d]+offset[d]); }

                  double distance(level_set.eval(coord));
                  for(int d=0;d<dim;++d)
                    distance=
                      std::min
                      (distance,
                       std::min(coord[d]-grid_geom->getXLower()[d]-epsilon,
                                grid_geom->getXUpper()[d]-coord[d]-epsilon));
                  (*level_set_ptr)(x)=distance/dx_max;
                }
            }
        }
    }
}
