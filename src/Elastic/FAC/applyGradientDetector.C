#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"

void Elastic::FAC::applyGradientDetector
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy_,
 const int ln,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
    hierarchy__ = hierarchy_;
  SAMRAI::hier::PatchHierarchy& hierarchy = *hierarchy__;
  SAMRAI::hier::PatchLevel& level =
    (SAMRAI::hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
  
  int ntag = 0, ntotal = 0;
  double maxestimate = 0;
  for(SAMRAI::hier::PatchLevel::Iterator pi(level.begin());
      pi!=level.end(); pi++)
    {
      SAMRAI::hier::Patch& patch = **pi;
      boost::shared_ptr<SAMRAI::hier::PatchData>
        tag_data = patch.getPatchData(tag_index);
      ntotal += patch.getBox().numberCells().getProduct();
      if (!tag_data)
        {
          TBOX_ERROR("Data index "
                     << tag_index << " does not exist for patch.\n");
        }
      boost::shared_ptr<SAMRAI::pdat::CellData<int> > tag_cell_data_ =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<int> >
        (tag_data);
      if (!tag_cell_data_)
        {
          TBOX_ERROR("Data index " << tag_index << " is not cell int data.\n");
        }
      boost::shared_ptr<SAMRAI::hier::PatchData>
        soln_data = patch.getPatchData(v_id);
      if (!soln_data)
        {
          TBOX_ERROR("Data index " << v_id << " does not exist for patch.\n");
        }
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > soln_side_data_ =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (soln_data);
      if (!soln_side_data_)
        {
          TBOX_ERROR("Data index " << v_id << " is not side data.\n");
        }
      SAMRAI::pdat::SideData<double>& v = *soln_side_data_;
      SAMRAI::pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
                              
      boost::shared_ptr<SAMRAI::hier::PatchData>
        dv_diagonal_data = patch.getPatchData(dv_diagonal_id);
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (dv_diagonal_data);
      if (!dv_diagonal_ptr)
        {
          TBOX_ERROR("Can not find dv_diagonal in applyGradientDetector.\n");
        }
      SAMRAI::pdat::CellData<double>& dv_diagonal = *dv_diagonal_ptr;

      boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (patch.getPatchGeometry());

      tag_cell_data.fill(0);
      SAMRAI::pdat::CellIterator cend(patch.getBox(),false);
      for (SAMRAI::pdat::CellIterator ci(patch.getBox(),true); ci!=cend; ci++)
        {
          const SAMRAI::pdat::CellIndex cell_index(*ci);

	  double curve(0.);
	  for (int ix=0; ix<d_dim.getValue(); ++ix)
	  {
            const SAMRAI::pdat::SideIndex x(cell_index,ix,
                                            SAMRAI::pdat::SideIndex::Lower);
	    for (int d=0; d<d_dim.getValue(); ++d)
              {
                SAMRAI::hier::Index ip(d_dim,0),
                  jp(d_dim,0),
                  kp(d_dim,0);
                ip(0)=1;
                jp(1)=1;
                if (3==d_dim.getValue())
                  {
                    kp(2)=1;
                  }
              const SAMRAI::hier::Index pp[]={ip,jp,kp};

              /* Special treatment near the boundary.  For Dirichlet
                 boundaries, the ghost point may not be valid. */
              if(cell_index[ix]==patch.getBox().lower(ix)
                 && geom->getTouchesRegularBoundary(ix,0))
                {
                  curve=std::max(curve,
                                 std::abs(v(x+pp[ix]+pp[ix])
                                          - 2*v(x+pp[ix]) + v(x)
                                          - dv_diagonal(cell_index+pp[ix],ix)
                                          + dv_diagonal(cell_index,ix)));
                    
                }
              else if(cell_index[ix]==patch.getBox().upper(ix)
                      && geom->getTouchesRegularBoundary(ix,1))
                {
                  curve=std::max(curve,
                                 std::abs(v(x+pp[ix])
                                          - 2*v(x) + v(x-pp[ix])
                                          - dv_diagonal(cell_index,ix)
                                          + dv_diagonal(cell_index-pp[ix],ix)));
                }
	      else
                {
                  curve=std::max(curve,
                                 std::abs(v(x+pp[ix]+pp[ix]) - v(x+pp[ix])
                                          - v(x) + v(x-pp[ix])
                                          - dv_diagonal(cell_index+pp[ix],ix)
                                          + dv_diagonal(cell_index-pp[ix],ix)));
                }
	    }
	  }

          if (maxestimate < curve)
               maxestimate=curve;
          if (curve > d_adaption_threshold || ln<min_full_refinement_level)
            {
              tag_cell_data(cell_index) = 1;
              ++ntag;
            }
        }
    }
  SAMRAI::tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  SAMRAI::tbox::plog << "Number of cells tagged on level " << ln << " is "
             << ntag << "/" << ntotal << "\n";
  SAMRAI::tbox::plog << "Max estimate is " << maxestimate << "\n";
}

