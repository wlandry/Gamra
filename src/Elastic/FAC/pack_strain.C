#include "Elastic/FAC.h"

/* Strain correction is stored on the cell centers and edges */

bool
Elastic::FAC::pack_strain(double* buffer,
                          const SAMRAI::hier::Patch& patch,
                          const SAMRAI::hier::Box& region,
                          const int &depth) const
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(v_id));
  SAMRAI::pdat::SideData<double>& v = *v_ptr;

  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal;
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed;
  if(!faults.empty())
    {
      dv_diagonal=boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch.getPatchData(dv_diagonal_id));
      dv_mixed=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(dv_mixed_id));
    }

  const int dim=d_dim.getValue();
  int ix(depth/dim), iy(depth%dim);
  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(d_dim));
  SAMRAI::hier::Index ip(zero), jp(zero);
  ip[ix]=1;
  jp[iy]=1;

  const double *dx(boost::dynamic_pointer_cast
                   <SAMRAI::geom::CartesianPatchGeometry>
                   (patch.getPatchGeometry())->getDx());

  SAMRAI::pdat::CellData<double>::iterator iend(region,false);
  for(SAMRAI::pdat::CellData<double>::iterator icell(region,true);
      icell!=iend; icell++)
    {
      const SAMRAI::pdat::SideIndex
        s(*icell,ix,SAMRAI::pdat::SideIndex::Lower);
      if(ix==iy)
        {
          double diff(v(s+ip)-v(s));
          if(!faults.empty())
            diff-=(*dv_diagonal)(*icell,ix);
          *buffer=diff/dx[ix];
        }
      else
        {
          const int ix_iy(index_map(ix,iy,dim));
          double diff(v(s) - v(s-jp) + v(s+ip) - v(s+ip-jp)
                      + v(s+jp) - v(s) + v(s+ip+jp) - v(s+ip));
          if(!faults.empty())
            diff+=(*dv_mixed)(s,ix_iy+1) - (*dv_mixed)(s-jp,ix_iy)
              + (*dv_mixed)(s+ip,ix_iy+1) - (*dv_mixed)(s+ip-jp,ix_iy)
              + (*dv_mixed)(s+jp,ix_iy+1) - (*dv_mixed)(s,ix_iy)
              + (*dv_mixed)(s+ip+jp,ix_iy+1) - (*dv_mixed)(s+ip,ix_iy);
          *buffer=diff/(4*dx[iy]);
        }
      ++buffer;
    }
  return true;
}
