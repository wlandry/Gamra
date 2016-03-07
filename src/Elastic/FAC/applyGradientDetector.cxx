/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"
#include "Constants.hxx"

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

void Elastic::FAC::applyGradientDetector
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy_,
 const int level,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
    hierarchy__ = hierarchy_;
  SAMRAI::hier::PatchHierarchy& hierarchy = *hierarchy__;
  SAMRAI::hier::PatchLevel& patch_level =
    (SAMRAI::hier::PatchLevel &) * hierarchy.getPatchLevel(level);
  
  size_t ntag = 0, ntotal = 0;
  double max_curvature(0);
  for(SAMRAI::hier::PatchLevel::Iterator p(patch_level.begin());
      p!=patch_level.end(); ++p)
    {
      SAMRAI::hier::Patch& patch = **p;
      ntotal += patch.getBox().numberCells().getProduct();
      boost::shared_ptr<SAMRAI::pdat::CellData<int> > tag_cell_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<int> >
        (patch.getPatchData(tag_index));
      SAMRAI::pdat::CellData<int>& tag_cell(*tag_cell_ptr);

      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(v_id));
      SAMRAI::pdat::SideData<double>& v(*v_ptr);
                              
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr;
      if(!faults.empty())
        {
          dv_diagonal_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(dv_diagonal_id));
        }

      boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr;
      if(have_embedded_boundary())
        {
          level_set_ptr=
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch.getPatchData(level_set_id));
        }

      boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (patch.getPatchGeometry());

      tag_cell.fill(0);
      const SAMRAI::hier::Box &box(patch.getBox());
      SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(box));
      const int dim(dimension.getValue());
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
                  SAMRAI::hier::Index ip(dimension,0), jp(dimension,0),
                    kp(dimension,0);
                  ip(0)=1;
                  jp(1)=1;
                  if (3==dim)
                    {
                      kp(2)=1;
                    }
                  const SAMRAI::hier::Index unit[]={ip,jp,kp};

                  /// Special treatment near the boundary.  For
                  /// Dirichlet boundaries, the ghost point may not be
                  /// valid.

                  double curve(0);
                  if(have_embedded_boundary())
                    {
                      SAMRAI::pdat::SideData<double> &level_set(*level_set_ptr);
                      // FIXME: Need to do the correct interpolation
                      // at the boundary.  For now assume v=0.
                      const double boundary(0);
                      double dv_plus, dv_minus;
                      if(level_set(x)<0)
                        {
                          if(!(level_set(x+unit[ix])<0))
                            {
                              dv_minus=v(x+unit[ix]) - boundary;
                              dv_plus=v(x+unit[ix]*2) - v(x+unit[ix]);
                              if(!faults.empty())
                                {
                                  dv_minus-=(*dv_diagonal_ptr)(cell,ix);
                                  dv_plus-=(*dv_diagonal_ptr)(cell+unit[ix],ix);
                                }
                            }
                          else
                            {
                              dv_minus=dv_plus=0;
                            }
                        }
                      else
                        {
                          dv_minus=v(x) - (level_set(x-unit[ix])<0 ? boundary
                                           : v(x-unit[ix]));
                          if(!faults.empty())
                            { dv_minus-=(*dv_diagonal_ptr)(cell-unit[ix],ix); }
                          if(level_set(x+unit[ix])<0)
                            {
                              dv_plus=boundary-v(x);
                              if(!faults.empty())
                                { dv_plus-=(*dv_diagonal_ptr)(cell,ix); }
                            }
                          else
                            {
                              dv_plus=(level_set(x+unit[ix]*2)<0 ? boundary
                                       : v(x+unit[ix]*2))
                                - v(x+unit[ix]);
                              if(!faults.empty())
                                {
                                  dv_plus-=(*dv_diagonal_ptr)(cell+unit[ix],ix);
                                }
                            }
                        }
                      curve=dv_plus-dv_minus;
                    }
                  else
                    {
                      if(cell[ix]==box.lower(ix)
                         && geom->getTouchesRegularBoundary(ix,0))
                        {
                          curve=v(x+unit[ix]*2) - 2*v(x+unit[ix]) + v(x);
                          if(!faults.empty())
                            curve+=-(*dv_diagonal_ptr)(cell+unit[ix],ix)
                              + (*dv_diagonal_ptr)(cell,ix);
                        }
                      else if(cell[ix]==box.upper(ix)
                              && geom->getTouchesRegularBoundary(ix,1))
                        {
                          curve=v(x+unit[ix]) - 2*v(x) + v(x-unit[ix]);
                          if(!faults.empty())
                            curve+=-(*dv_diagonal_ptr)(cell,ix)
                              + (*dv_diagonal_ptr)(cell-unit[ix],ix);
                        }
                      else
                        {
                          curve=v(x+unit[ix]*2) - v(x+unit[ix])
                            - v(x) + v(x-unit[ix]);
                          if(!faults.empty())
                            curve+=-(*dv_diagonal_ptr)(cell+unit[ix],ix)
                              + (*dv_diagonal_ptr)(cell-unit[ix],ix);
                        }
                    }
                  curvature=std::max(curvature,std::abs(curve));
                }
            }

          if(max_curvature < curvature)
            { max_curvature=curvature; }

          if (curvature > d_adaption_threshold
              || level<min_full_refinement_level)
            {
              tag_cell(cell) = 1;
              ++ntag;
            }
        }
      for(size_t i=0;i<refinement_points.size();i+=3)
        {
          SAMRAI::pdat::CellIndex cell(dimension);
          for(int d=0;d<dim;++d)
            cell[d]=static_cast<int>((refinement_points[i+d]
                                      - geom->getXLower()[d])/geom->getDx()[d]
                                     + 0.5);
          cell+=box.lower();
          if(box.contains(cell) && tag_cell(cell)!=1)
            {
              tag_cell(cell)=1;
              ++ntag;
            }
        }
    }
  SAMRAI::tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  SAMRAI::tbox::plog << "Number of cells tagged on level " << level << " is "
             << ntag << "/" << ntotal << "\n";
  SAMRAI::tbox::plog << "Max estimate is " << max_curvature << "\n";
}
