/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "../v_operator_2D.hxx"
#include "../v_level_set_operator_2D.hxx"
#include "Constants.hxx"

namespace Elastic
{
  void residual_embedded_2D
  (SAMRAI::pdat::SideData<double> &v,
   SAMRAI::pdat::CellData<double> &cell_moduli,
   SAMRAI::pdat::NodeData<double> &edge_moduli,
   SAMRAI::pdat::SideData<double> &v_rhs,
   SAMRAI::pdat::SideData<double> &v_resid,
   SAMRAI::pdat::SideData<double> &level_set,
   const SAMRAI::hier::Box &pbox,
   const double dxy[])
  {
    const SAMRAI::hier::Index unit[]={SAMRAI::hier::Index(1,0),
                                      SAMRAI::hier::Index(0,1)};
    double dx = dxy[0];
    double dy = dxy[1];
    SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(pbox));
    for(SAMRAI::pdat::CellIterator
          ci(SAMRAI::pdat::CellGeometry::begin(pbox)); ci!=cend; ++ci)
      {
        const SAMRAI::pdat::CellIndex &cell(*ci);

        const SAMRAI::pdat::NodeIndex
          edge(cell,SAMRAI::pdat::NodeIndex::LowerLeft);

        const int dim(2);
        for(Gamra::Dir ix=0;ix<dim;++ix)
          {
            const Gamra::Dir iy(ix.next(dim));
            const SAMRAI::pdat::SideIndex
              x(cell,ix,SAMRAI::pdat::SideIndex::Lower),
              y(cell,iy,SAMRAI::pdat::SideIndex::Lower);
            const SAMRAI::hier::Index ip(unit[ix]), jp(unit[iy]);

            if(cell[iy]!=pbox.upper(iy))
              {
                if(level_set(x)<0)
                  {
                    v_resid(x)=0;
                  }
                else if(level_set(x)>1)
                  {
                    v_resid(x)=v_rhs(x)
                      - v_operator_2D(v,cell_moduli,edge_moduli,cell,
                                      edge,x,y,ip,jp,dx,dy);
                  }
                else
                  {
                    v_resid(x)=v_rhs(x)
                      - v_level_set_operator_2D(level_set,v,cell_moduli,
                                                edge_moduli,cell,
                                                edge,x,y,ip,jp,dx,dy);
                  }
              }
          }
      }
  }
}
