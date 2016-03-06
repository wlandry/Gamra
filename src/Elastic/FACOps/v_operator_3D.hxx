#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "lame_mixed.hxx"
#include "aligned_terms.hxx"
#include "shear_noncell.hxx"

#include <SAMRAI/pdat/EdgeData.h>

namespace Elastic
{
  inline double v_operator_3D(const SAMRAI::pdat::SideData<double> &v,
                              const SAMRAI::pdat::CellData<double> &cell_moduli,
                              const SAMRAI::pdat::EdgeData<double> &edge_moduli,
                              const SAMRAI::pdat::CellIndex &cell,
                              const SAMRAI::pdat::EdgeIndex &edge_y,
                              const SAMRAI::pdat::EdgeIndex &edge_z,
                              const SAMRAI::pdat::SideIndex &x,
                              const SAMRAI::pdat::SideIndex &y,
                              const SAMRAI::pdat::SideIndex &z,
                              const SAMRAI::hier::Index &ip,
                              const SAMRAI::hier::Index &jp,
                              const SAMRAI::hier::Index &kp,
                              const double &dx,
                              const double &dy,
                              const double &dz)
  {
    return aligned_terms(v,cell_moduli,cell,x,ip,dx)
      +lame_mixed(v,cell_moduli,cell,y,ip,jp,dx,dy)
      +lame_mixed(v,cell_moduli,cell,z,ip,kp,dx,dz)
      +shear_noncell(v,edge_moduli,x,y,edge_z,ip,jp,dx,dy)
      +shear_noncell(v,edge_moduli,x,z,edge_y,ip,kp,dx,dz);
  }
}
