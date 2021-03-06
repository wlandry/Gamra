/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "../v_operator_3D.hxx"
#include "Constants.hxx"

namespace Elastic
{
  void residual_3D
  (SAMRAI::pdat::SideData<double> &v,
   SAMRAI::pdat::CellData<double> &cell_moduli,
   SAMRAI::pdat::EdgeData<double> &edge_moduli,
   SAMRAI::pdat::SideData<double> &v_rhs,
   SAMRAI::pdat::SideData<double> &v_resid,
   const SAMRAI::hier::Box &pbox,
   const double dxyz[])
  {
    const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
    const SAMRAI::hier::Index unit[]={ip,jp,kp};

    SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(pbox));
    for(SAMRAI::pdat::CellIterator ci(SAMRAI::pdat::CellGeometry::begin(pbox));
        ci!=cend; ++ci)
      {
        const SAMRAI::pdat::CellIndex &center(*ci);
        SAMRAI::pdat::CellIndex up(center), down(center), right(center),
          left(center), front(center), back(center);

        ++right[0];
        --left[0];
        ++up[1];
        --down[1];
        ++front[2];
        --back[2];

        for(Gamra::Dir ix=0;ix<3;++ix)
          {
            const Gamra::Dir iy(ix.next(3));
            const Gamra::Dir iz(iy.next(3));
            const SAMRAI::pdat::SideIndex
              x(center,ix,SAMRAI::pdat::SideIndex::Lower),
              y(center,iy,SAMRAI::pdat::SideIndex::Lower),
              z(center,iz,SAMRAI::pdat::SideIndex::Lower);
            const SAMRAI::pdat::EdgeIndex
              edge_y(center,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
              edge_z(center,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);

            if(center[iy]!=pbox.upper(iy) && center[iz]!=pbox.upper(iz))
              {
                if((center[ix]==pbox.lower(ix) && v(x-unit[ix])==boundary_value)
                   || (center[ix]==pbox.upper(ix)
                       && v(x+unit[ix])==boundary_value))
                  {
                    v_resid(x)=0;
                  }
                else
                  {
                    v_resid(x)=v_rhs(x)
                      - v_operator_3D(v,cell_moduli,edge_moduli,
                                      center,edge_y,edge_z,x,y,z,
                                      unit[ix],unit[iy],unit[iz],
                                      dxyz[ix],dxyz[iy],dxyz[iz]);
                  }
              }
          }          
      }
  }
}
