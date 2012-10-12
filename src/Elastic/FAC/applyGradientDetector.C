#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"

void Elastic::FAC::applyGradientDetector
(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy> hierarchy_,
 const int ln,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy>
    hierarchy__ = hierarchy_;
  SAMRAI::hier::PatchHierarchy& hierarchy = *hierarchy__;
  SAMRAI::hier::PatchLevel& level =
    (SAMRAI::hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
  
  int ntag = 0, ntotal = 0;
  double maxestimate = 0;
  for(SAMRAI::hier::PatchLevel::Iterator pi(level); pi; pi++)
    {
      SAMRAI::hier::Patch& patch = **pi;
      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData>
        tag_data = patch.getPatchData(tag_index);
      ntotal += patch.getBox().numberCells().getProduct();
      if (tag_data.isNull())
        {
          TBOX_ERROR("Data index "
                     << tag_index << " does not exist for patch.\n");
        }
      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<int> >
        tag_cell_data_ = tag_data;
      if (tag_cell_data_.isNull())
        {
          TBOX_ERROR("Data index " << tag_index << " is not cell int data.\n");
        }
      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData>
        soln_data = patch.getPatchData(v_id);
      if (soln_data.isNull())
        {
          TBOX_ERROR("Data index " << v_id << " does not exist for patch.\n");
        }
      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
        soln_side_data_ = soln_data;
      if (soln_side_data_.isNull())
        {
          TBOX_ERROR("Data index " << v_id << " is not side data.\n");
        }
      SAMRAI::pdat::SideData<double>& v = *soln_side_data_;
      SAMRAI::pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
                              
      SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
        geom = patch.getPatchGeometry();

      tag_cell_data.fill(0);
      for (SAMRAI::pdat::CellIterator ci(patch.getBox()); ci; ci++)
        {
          const SAMRAI::pdat::CellIndex cell_index(*ci);

	  double curve(0.);
	  for (int ix=0; ix<d_dim.getValue(); ++ix)
	  {
            const SAMRAI::pdat::SideIndex x(cell_index,ix,
                                            SAMRAI::pdat::SideIndex::Lower);
	    for (int d=0; d<d_dim.getValue(); ++d){
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

              if(cell_index[ix]==patch.getBox().lower(ix)
                 && geom->getTouchesRegularBoundary(ix,0))
	      {
	      curve=
                std::max(curve,std::abs(v(x+pp[ix]+pp[ix])-2*v(x+pp[ix])+v(x)));
	      } else if(cell_index[ix]==patch.getBox().upper(ix)
                 && geom->getTouchesRegularBoundary(ix,1))
	      {
	      curve=std::max(curve,std::abs(v(x+pp[ix])-2*v(x)+v(x-pp[ix])));
	      }
	      else
	      {
	      curve=std::max(curve,std::abs(v(x+pp[ix]+pp[ix]) - v(x+pp[ix])
                                            - v(x) + v(x-pp[ix])));
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

