#include "Elastic/FAC.h"

/* Strain correction is stored on the cell centers and edges */

bool
Elastic::FAC::pack_diagonal_strain(double* buffer,
                                   const SAMRAI::hier::Patch& patch,
                                   const SAMRAI::hier::Box& region,
                                   const int &ix) const
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(v_id));
  SAMRAI::pdat::SideData<double>& v = *v_ptr;

  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_aligned_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
    (patch.getPatchData(dv_aligned_id));
  SAMRAI::pdat::CellData<double> &dv_aligned(*dv_aligned_ptr);

  const int dim=d_dim.getValue();
  SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(d_dim));
  ip[ix]=1;

  if(d_dim.getValue()==2)
    {
      SAMRAI::pdat::CellData<double>::iterator iend(region,false);
      for(SAMRAI::pdat::CellData<double>::iterator icell(region,true);
          icell!=iend; icell++)
        {
          const SAMRAI::pdat::SideIndex
            s(*icell,ix,SAMRAI::pdat::SideIndex::Lower);
          *buffer=v(s+ip)-v(s)-dv_aligned(*icell,ix);
          ++buffer;
        }
    }
  else
    {
      SAMRAI::pdat::CellData<double>::iterator iend(region,false);
      for(SAMRAI::pdat::CellData<double>::iterator icell(region,true);
          icell!=iend; icell++)
        {
          const SAMRAI::pdat::SideIndex
            s(*icell,ix,SAMRAI::pdat::SideIndex::Lower);
          *buffer=v(s+ip)-v(s)-dv_aligned(*icell,ix);
          ++buffer;
        }
    }
  return true;
}
