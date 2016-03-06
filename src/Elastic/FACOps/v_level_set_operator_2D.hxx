#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/NodeData.h>

namespace Elastic
{
  double v_level_set_operator_2D
  (const SAMRAI::pdat::SideData<double> &level_set,
   const SAMRAI::pdat::SideData<double> &v,
   const SAMRAI::pdat::CellData<double> &cell_moduli,
   const SAMRAI::pdat::NodeData<double> &edge_moduli,
   const SAMRAI::pdat::CellIndex &cell,
   const SAMRAI::pdat::NodeIndex &edge,
   const SAMRAI::pdat::SideIndex &x,
   const SAMRAI::pdat::SideIndex &y,
   const SAMRAI::hier::Index &ip,
   const SAMRAI::hier::Index &jp,
   const double &dx,
   const double &dy);
}

