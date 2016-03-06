#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

namespace Elastic
{
  inline double aligned_terms(const SAMRAI::pdat::SideData<double> &v,
                              const SAMRAI::pdat::CellData<double> &cell_moduli,
                              const SAMRAI::pdat::CellIndex &cell,
                              const SAMRAI::pdat::SideIndex &x,
                              const SAMRAI::hier::Index &ip,
                              const double &dx)
  {
    return (( v(x+ip)-v(x   )) *(cell_moduli(cell   ,0)+2*cell_moduli(cell   ,1))
            -(v(x   )-v(x-ip)) *(cell_moduli(cell-ip,0)+2*cell_moduli(cell-ip,1)))
      /(dx*dx);
  }
}
