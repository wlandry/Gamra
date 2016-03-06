#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"

template<class T>
void correct_rhs(const SAMRAI::tbox::Dimension &dimension,
                 const Gamra::Dir &dim, const double *dx,
                 const SAMRAI::hier::Box &pbox,
                 const SAMRAI::pdat::CellData<double> &cell_moduli,
                 const T &edge_moduli,
                 const SAMRAI::pdat::CellData<double> &dv_diagonal,
                 const SAMRAI::pdat::SideData<double> &dv_mixed,
                 SAMRAI::pdat::SideData<double> &v_rhs)
{
  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(dimension));
  SAMRAI::hier::Index unit[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    unit[i][i]=1;

  for(Gamra::Dir ix=0;ix<dim;++ix)
    {
      SAMRAI::pdat::SideIterator
        s_end(SAMRAI::pdat::SideGeometry::end(pbox,ix));
      for(SAMRAI::pdat::SideIterator
            si(SAMRAI::pdat::SideGeometry::begin(pbox,ix));
          si!=s_end; ++si)
        {
          const SAMRAI::pdat::SideIndex &s(*si);
          SAMRAI::pdat::CellIndex c(s);

          /// d/dx^2, d/dy^2, d/dz^2

          double lambda_plus(cell_moduli(c,0)),
            lambda_minus(cell_moduli(c-unit[ix],0)),
            mu_plus(cell_moduli(c,1)),
            mu_minus(cell_moduli(c-unit[ix],1));
          v_rhs(s)+=
            (dv_diagonal(c,ix)*(lambda_plus + 2*mu_plus)
             - dv_diagonal(c-unit[ix],ix)*(lambda_minus + 2*mu_minus))
            /(dx[ix]*dx[ix]);

          for(Gamra::Dir iy=ix.next(dim);iy!=ix;
              iy=iy.next(dim))
            {
              const Gamra::Dir
                iz(ix.next(dim)!=iy ? ix.next(dim)
                   : ix.next(dim).next(dim));
              mu_plus=edge_node_eval(edge_moduli,s+unit[iy],iz,1);
              mu_minus=edge_node_eval(edge_moduli,s,iz,1);

              const Gamra::Dir
                ix_iy(index_map(ix,iy,dim));
              v_rhs(s)+=
                (mu_plus*(dv_mixed(s,ix_iy)
                          - dv_mixed(s+unit[iy],ix_iy+1))
                 + mu_minus*(dv_mixed(s,ix_iy+1)
                             - dv_mixed(s-unit[iy],ix_iy)))
                /(dx[iy]*dx[iy]);

              /// d/dxy
              lambda_plus=cell_moduli(c,0);
              lambda_minus=cell_moduli(c-unit[ix],0);
              const SAMRAI::pdat::SideIndex
                s_y(c,iy,SAMRAI::pdat::SideIndex::Lower);
              const Gamra::Dir iy_ix(index_map(iy,ix,dim));

              v_rhs(s)+=
                (lambda_plus*dv_diagonal(c,iy)
                 - lambda_minus*dv_diagonal(c-unit[ix],iy)
                 - mu_plus*(dv_mixed(s_y+unit[iy],iy_ix+1)
                            - dv_mixed(s_y+unit[iy]-unit[ix],iy_ix))
                 + mu_minus*(dv_mixed(s_y,iy_ix+1)
                             - dv_mixed(s_y-unit[ix],iy_ix)))
                /(dx[ix]*dx[iy]);
            }
        }
    }
}
