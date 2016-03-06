#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/CellData.h>

namespace Elastic
{
  inline double lame_mixed(const SAMRAI::pdat::SideData<double> &v,
                           const SAMRAI::pdat::CellData<double> &cell_moduli,
                           const SAMRAI::pdat::CellIndex &cell,
                           const SAMRAI::pdat::SideIndex &y,
                           const SAMRAI::hier::Index &ip,
                           const SAMRAI::hier::Index &jp,
                           const double &dx,
                           const double &dy)
  {
    return (  cell_moduli(cell   ,0)*(v(y+jp   )-v(y   ))/dy
              - cell_moduli(cell-ip,0)*(v(y+jp-ip)-v(y-ip))/dy)/dx;
  }
}
