#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

/// The action of the velocity operator. It is written from the
/// perspective of vx, but pass in different values for cell_x
/// etc. to get vy.

#include "lame_mixed.hxx"
#include "aligned_terms.hxx"
#include "shear_noncell.hxx"

namespace Elastic
{
  inline double v_operator_2D(const SAMRAI::pdat::SideData<double> &v,
                              const SAMRAI::pdat::CellData<double> &cell_moduli,
                              const SAMRAI::pdat::NodeData<double> &edge_moduli,
                              const SAMRAI::pdat::CellIndex &cell,
                              const SAMRAI::pdat::NodeIndex &edge,
                              const SAMRAI::pdat::SideIndex &x,
                              const SAMRAI::pdat::SideIndex &y,
                              const SAMRAI::hier::Index &ip,
                              const SAMRAI::hier::Index &jp,
                              const double &dx,
                              const double &dy)
  {
    return aligned_terms(v,cell_moduli,cell,x,ip,dx)
      +lame_mixed(v,cell_moduli,cell,y,ip,jp,dx,dy)
      +shear_noncell(v,edge_moduli,x,y,edge,ip,jp,dx,dy);
  }
}
