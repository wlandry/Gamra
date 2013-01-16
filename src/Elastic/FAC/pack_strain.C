#include "Elastic/FAC.h"

/* Strain correction is stored on the cell centers and edges */

bool
Elastic::FAC::pack_strain(double* buffer,
                          const SAMRAI::hier::Patch& patch,
                          const SAMRAI::hier::Box& region,
                          int depth) const
{
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_aligned_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
    (patch.getPatchData(dv_aligned_id));
  SAMRAI::pdat::CellData<double> &dv_aligned(*dv_aligned_ptr);

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
      SAMRAI::pdat::CellData<double>::iterator iend(region,false);
      for(SAMRAI::pdat::CellData<double>::iterator icell(region,true);
          icell!=iend; icell++)
        {
          if(ix==d)
            {
              *buffer=dv_aligned(*icell,ix);
            }
          else
            {
              const SAMRAI::pdat::NodeIndex
                corner(*icell,SAMRAI::pdat::NodeIndex::LowerLeft);
              *buffer=(dv_perpendicular(corner,index_map[ix][d])
                       + dv_perpendicular(corner+ip,index_map[ix][d])
                       + dv_perpendicular(corner+jp,index_map[ix][d])
                       + dv_perpendicular(corner+ip+jp,index_map[ix][d]))/4;

            }
          buffer = buffer + 1;
        }
    }
  else
    {
      int other_dim((ix+1)%dim == d ? (d+1)%dim : (ix+1)%dim);
      boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > dv_perpendicular_ptr=
        boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
        (patch.getPatchData(dv_perpendicular_id));
      SAMRAI::pdat::EdgeData<double> &dv_perpendicular(*dv_perpendicular_ptr);
      SAMRAI::pdat::CellData<double>::iterator iend(region,false);
      for(SAMRAI::pdat::CellData<double>::iterator icell(region,true);
          icell!=iend; icell++)
        {
          if(ix==d)
            {
              *buffer=dv_aligned(*icell,ix);
            }
          else
            {
              const SAMRAI::pdat::EdgeIndex
                corner(*icell,other_dim,SAMRAI::pdat::EdgeIndex::LowerLeft);
              *buffer=(dv_perpendicular(corner,index_map[ix][d])
                       + dv_perpendicular(corner+ip,index_map[ix][d])
                       + dv_perpendicular(corner+jp,index_map[ix][d])
                       + dv_perpendicular(corner+ip+jp,index_map[ix][d]))/4;
            }
          buffer = buffer + 1;
        }
    }
  return true;
}
