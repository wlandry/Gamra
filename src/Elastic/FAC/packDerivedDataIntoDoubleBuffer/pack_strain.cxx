/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"
#include "Constants.hxx"

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

/// Strain correction is stored on the cell centers and edges

void pack_strain(double* buffer,
                 const SAMRAI::hier::Patch& patch,
                 const SAMRAI::hier::Box& region,
                 const int &depth,
                 const SAMRAI::tbox::Dimension &dimension,
                 const bool &have_faults,
                 const bool &have_embedded_boundary,
                 const int &v_id,
                 const int &dv_diagonal_id,
                 const int &dv_mixed_id,
                 const int &level_set_id)
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(v_id));
  SAMRAI::pdat::SideData<double>& v = *v_ptr;

  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal;
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed;
  if(have_faults)
    {
      dv_diagonal=boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch.getPatchData(dv_diagonal_id));
      dv_mixed=boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(dv_mixed_id));
    }

  const Gamra::Dir dim=dimension.getValue();
  if(have_embedded_boundary)
    {
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr;
      if(have_embedded_boundary)
        {
          level_set_ptr=
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch.getPatchData(level_set_id));
        }
    }

  Gamra::Dir ix(Gamra::Dir::from_int(depth/dim)),
    iy(Gamra::Dir::from_int(depth%dim));
  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(dimension));
  SAMRAI::hier::Index ip(zero), jp(zero);
  ip[ix]=1;
  jp[iy]=1;

  const double *dx(boost::dynamic_pointer_cast
                   <SAMRAI::geom::CartesianPatchGeometry>
                   (patch.getPatchGeometry())->getDx());

  SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
  for(SAMRAI::pdat::CellIterator
        icell(SAMRAI::pdat::CellGeometry::begin(region));
      icell!=iend; ++icell)
    {
      const SAMRAI::pdat::SideIndex
        s(*icell,ix,SAMRAI::pdat::SideIndex::Lower);
      if(ix==iy)
        {
          double diff(v(s+ip)-v(s));
          if(have_faults)
            { diff-=(*dv_diagonal)(*icell,ix); }
          *buffer=diff/dx[ix];
        }
      else
        {
          const int ix_iy(index_map(ix,iy,dim));
          double diff(- v(s-jp) - v(s+ip-jp) + v(s+jp)  + v(s+ip+jp));

          if(have_faults)
            { diff+=(*dv_mixed)(s,ix_iy+1) - (*dv_mixed)(s-jp,ix_iy)
                + (*dv_mixed)(s+ip,ix_iy+1) - (*dv_mixed)(s+ip-jp,ix_iy)
                + (*dv_mixed)(s+jp,ix_iy+1) - (*dv_mixed)(s,ix_iy)
                + (*dv_mixed)(s+ip+jp,ix_iy+1) - (*dv_mixed)(s+ip,ix_iy); }
          *buffer=diff/(4*dx[iy]);
        }
      ++buffer;
    }
}
