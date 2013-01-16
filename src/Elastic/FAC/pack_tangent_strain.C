#include "Elastic/FAC.h"

/* Strain correction is stored on the cell centers and edges */

bool
Elastic::FAC::pack_tangent_strain(double* buffer,
                                  const SAMRAI::hier::Patch& patch,
                                  const SAMRAI::hier::Box& region,
                                  const int &depth) const
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(v_id));
  SAMRAI::pdat::SideData<double>& v = *v_ptr;

  const int dim=d_dim.getValue();
  int ix(depth/dim), d(depth%dim);
  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(d_dim));
  SAMRAI::hier::Index ip(zero), jp(zero);
  ip[ix]=1;
  jp[d]=1;

  if(d_dim.getValue()==2)
    {
      boost::shared_ptr<SAMRAI::pdat::NodeData<double> > dv_perpendicular_ptr=
        boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
        (patch.getPatchData(dv_perpendicular_id));
      SAMRAI::pdat::NodeData<double> &dv_perpendicular(*dv_perpendicular_ptr);
      SAMRAI::pdat::NodeData<double>::iterator iend(region,false);
      for(SAMRAI::pdat::NodeData<double>::iterator inode(region,true);
          inode!=iend; inode++)
        {
          if(ix==d)
            {
              *buffer=0;
            }
          else
            {
              const SAMRAI::pdat::SideIndex
                s(*inode,ix,SAMRAI::pdat::SideIndex::Lower);
              *buffer=v(s)-v(s-jp)-dv_perpendicular(*inode,index_map[ix][d]);
            }
          ++buffer;
        }
    }
  else
    {
      int other_dim((ix+1)%dim == d ? (d+1)%dim : (ix+1)%dim);
      boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > dv_perpendicular_ptr=
        boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
        (patch.getPatchData(dv_perpendicular_id));
      SAMRAI::pdat::EdgeData<double> &dv_perpendicular(*dv_perpendicular_ptr);
      SAMRAI::pdat::EdgeData<double>::iterator iend(region,other_dim,false);
      for(SAMRAI::pdat::EdgeData<double>::iterator iedge(region,other_dim,true);
          iedge!=iend; iedge++)
        {
          if(ix==d)
            {
              *buffer=0;
            }
          else
            {
              const SAMRAI::pdat::SideIndex
                s(*iedge,ix,SAMRAI::pdat::SideIndex::Lower);
              *buffer=v(s)-v(s-jp)-dv_perpendicular(*iedge,index_map[ix][d]);
            }
          ++buffer;
        }
    }
  return true;
}
