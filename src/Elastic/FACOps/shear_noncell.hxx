#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

/// The mixed derivative of the stress.  We have to use a template
/// because 2D uses Node's for the edge moduli, while 3D uses
/// Edge's.  Written as if it is dtau_xy_dy.

namespace Elastic
{
  template<class E_data, class E_index>
  double shear_noncell(const SAMRAI::pdat::SideData<double> &v,
                       const E_data &edge_moduli,
                       const SAMRAI::pdat::SideIndex &x,
                       const SAMRAI::pdat::SideIndex &y,
                       const E_index &edge,
                       const SAMRAI::hier::Index &ip,
                       const SAMRAI::hier::Index &jp,
                       const double &dx,
                       const double &dy)
  {
    return 
      edge_moduli(edge+jp,1)*(v(x+jp)-v(x   ))/(dy*dy)
      -edge_moduli(edge   ,1)*(v(x   )-v(x-jp))/(dy*dy)
      +edge_moduli(edge+jp,1)*(v(y+jp)-v(y+jp-ip))/(dx*dy) 
      -edge_moduli(edge   ,1)*(v(y   )-v(y-ip   ))/(dx*dy);
  }
}
