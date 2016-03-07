#include "Stokes/FAC.hxx"
#include <SAMRAI/geom/CartesianGridGeometry.h>

void Stokes::FAC::applyGradientDetector
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy_,
 const int level,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy__ = hierarchy_;
  SAMRAI::hier::PatchHierarchy& hierarchy = *hierarchy__;
  SAMRAI::hier::PatchLevel& patch_level =
    (SAMRAI::hier::PatchLevel &) * hierarchy.getPatchLevel(level);
  
  size_t ntag = 0, ntotal = 0;
  double max_curvature = 0;
  for(SAMRAI::hier::PatchLevel::Iterator p(patch_level.begin());
      p!=patch_level.end(); ++p)
    {
      SAMRAI::hier::Patch& patch = **p;
      boost::shared_ptr<SAMRAI::hier::PatchData>
        tag_data = patch.getPatchData(tag_index);
      ntotal += patch.getBox().numberCells().getProduct();
      if (!tag_data)
        {
          TBOX_ERROR("Data index "
                     << tag_index << " does not exist for patch.\n");
        }
      boost::shared_ptr<SAMRAI::pdat::CellData<int> > tag_cell_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<int> >
        (patch.getPatchData(tag_index));
      SAMRAI::pdat::CellData<int>& tag_cell(*tag_cell_ptr);

      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(v_id));
      SAMRAI::pdat::SideData<double>& v(*v_ptr);

      boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (patch.getPatchGeometry());

      tag_cell.fill(0);
      const SAMRAI::hier::Box &box(patch.getBox());
      SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(box));
      const int dim(d_dim.getValue());

      for(SAMRAI::pdat::CellIterator ci(SAMRAI::pdat::CellGeometry::begin(box));
          ci!=cend; ++ci)
        {
          const SAMRAI::pdat::CellIndex &cell(*ci);

          double curvature(0);
	  for (Gamra::Dir ix=0; ix<dim; ++ix)
            {
              const SAMRAI::pdat::SideIndex x(cell,ix,
                                              SAMRAI::pdat::SideIndex::Lower);
              for (int d=0; d<dim; ++d)
                {
                  SAMRAI::hier::Index ip(d_dim,0),
                    jp(d_dim,0),
                    kp(d_dim,0);
                  ip(0)=1;
                  jp(1)=1;
                  if (3==dim)
                    {
                      kp(2)=1;
                    }
                  const SAMRAI::hier::Index unit[]={ip,jp,kp};

                  /* Special treatment near the boundary.  For Dirichlet
                     boundaries, the ghost point may not be valid. */

                  double curve(0);

                  if(cell[ix]==box.lower(ix)
                     && geom->getTouchesRegularBoundary(ix,0))
                    {
                      curve=v(x+unit[ix]*2) - 2*v(x+unit[ix]) + v(x);
                    }
                  else if(cell[ix]==box.upper(ix)
                          && geom->getTouchesRegularBoundary(ix,1))
                    {
                      curve=v(x+unit[ix]) - 2*v(x) + v(x-unit[ix]);
                    }
                  else
                    {
                      curve=v(x+unit[ix]*2) - v(x+unit[ix])
                        - v(x) + v(x-unit[ix]);
                    }
                  curvature=std::max(curvature,std::abs(curve));
                }
            }

          if(max_curvature < curvature)
            max_curvature=curvature;

          if (curvature > d_adaption_threshold
              || level<min_full_refinement_level)
            {
              tag_cell(cell) = 1;
              ++ntag;
            }
        }
    }
  SAMRAI::tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  SAMRAI::tbox::plog << "Number of cells tagged on level " << level << " is "
             << ntag << "/" << ntotal << "\n";
  SAMRAI::tbox::plog << "Max estimate is " << max_curvature << "\n";
}
