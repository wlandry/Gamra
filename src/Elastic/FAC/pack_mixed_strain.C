#include "Elastic/FAC.h"

/* Store the mixed derivatives of the displacement.  The derivative
   correction terms for mixed derivatives are stored on the faces,
   with one component being the correction 'up' and one component
   being the correction 'down.  */

bool
Elastic::FAC::pack_mixed_strain(double* buffer,
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

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(dv_mixed_id));
  SAMRAI::pdat::SideData<double> &dv_mixed(*dv_mixed_ptr);

  if(d_dim.getValue()==2)
    {
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
              *buffer=v(s) - v(s-jp) + dv_mixed(s,1) - dv_mixed(s-jp,0);
            }
          ++buffer;
        }
    }
  else
    {
      /* This is not going to work because Visit can not handle edge
         data.  So we have to cell centered for it anyway. */
      int other_dim((ix+1)%dim == d ? (d+1)%dim : (ix+1)%dim);
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
              *buffer=v(s) - v(s-jp) + dv_mixed(s,2*((d-ix)%(dim-1))+1)
                - dv_mixed(s,2*((d-ix)%(dim-1)));
            }
          ++buffer;
        }
    }
  return true;
}
