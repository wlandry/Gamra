#include "FACStokes.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"

void SAMRAI::FACStokes::applyGradientDetector
(const tbox::Pointer<hier::BasePatchHierarchy> hierarchy_,
 const int ln,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const tbox::Pointer<hier::PatchHierarchy> hierarchy__ = hierarchy_;
  hier::PatchHierarchy& hierarchy = *hierarchy__;
  hier::PatchLevel& level =
    (hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
  
  int ntag = 0, ntotal = 0;
  double maxestimate = 0;
  for(hier::PatchLevel::Iterator pi(level); pi; pi++)
    {
      hier::Patch& patch = **pi;
      tbox::Pointer<hier::PatchData>
        tag_data = patch.getPatchData(tag_index);
      ntotal += patch.getBox().numberCells().getProduct();
      if (tag_data.isNull())
        {
          TBOX_ERROR("Data index "
                     << tag_index << " does not exist for patch.\n");
        }
      tbox::Pointer<pdat::CellData<int> > tag_cell_data_ = tag_data;
      if (tag_cell_data_.isNull())
        {
          TBOX_ERROR("Data index " << tag_index << " is not cell int data.\n");
        }
      tbox::Pointer<hier::PatchData> soln_data = patch.getPatchData(v_id);
      if (soln_data.isNull())
        {
          TBOX_ERROR("Data index " << v_id << " does not exist for patch.\n");
        }
      tbox::Pointer<pdat::SideData<double> > soln_side_data_ = soln_data;
      if (soln_side_data_.isNull())
        {
          TBOX_ERROR("Data index " << v_id << " is not side data.\n");
        }
      pdat::SideData<double>& v = *soln_side_data_;
      pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
                              
      tbox::Pointer<geom::CartesianPatchGeometry> geom = patch.getPatchGeometry();

      tag_cell_data.fill(0);
      for (pdat::CellIterator ci(patch.getBox()); ci; ci++)
        {
          const pdat::CellIndex cell_index(*ci);

	  double curve(0.);
	  for (int ix=0; ix<d_dim.getValue(); ++ix)
	  {
            const pdat::SideIndex x(cell_index,ix,pdat::SideIndex::Lower);
	    for (int d=0; d<d_dim.getValue(); ++d){
		    hier::Index ip(d_dim,0),
			        jp(d_dim,0),
			        kp(d_dim,0);
		    ip(0)=1;
		    jp(1)=1;
		    if (3==d_dim.getValue())
		    {
			kp(2)=1;
		    }
              const hier::Index pp[]={ip,jp,kp};

              if(cell_index[ix]==patch.getBox().lower(ix)
                 && geom->getTouchesRegularBoundary(ix,0))
	      {
	      curve=std::max(curve,std::abs(v(x+pp[ix]+pp[ix])-2*v(x+pp[ix])+v(x)));
	      } else if(cell_index[ix]==patch.getBox().upper(ix)
                 && geom->getTouchesRegularBoundary(ix,1))
	      {
	      curve=std::max(curve,std::abs(v(x+pp[ix])-2*v(x)+v(x-pp[ix])));
	      }
	      else
	      {
	      curve=std::max(curve,std::abs(v(x+pp[ix]+pp[ix])-v(x+pp[ix])-v(x)+v(x-pp[ix])));
	      }
	    }
	  }
          /*
	   * tbox::plog << "estimate "
                     << cell_index << " "
                     << d_adaption_threshold << " "
                     << curve << " "
                     << std::boolalpha
                     << (curve > d_adaption_threshold)
                     << " "
                     << "\n";
		     */

          if (maxestimate < curve)
               maxestimate=curve;
          if (curve > d_adaption_threshold || ln<min_full_refinement_level)
            {
              tag_cell_data(cell_index) = 1;
              ++ntag;
            }
        }
    }
  tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  tbox::plog << "Number of cells tagged on level " << ln << " is "
             << ntag << "/" << ntotal << "\n";
  tbox::plog << "Max estimate is " << maxestimate << "\n";
}

