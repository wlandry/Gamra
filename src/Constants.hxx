/* Special values such as whether we are at a boundary, plus the
   function index_map */

#pragma once

#include <limits>
#include "Dir.hxx"

const double boundary_value=1.1111e100;
const double invalid_value=3.1415926*boundary_value;
const int invalid_id=-1;
inline Gamra::Dir index_map(const Gamra::Dir &ix,
                            const Gamra::Dir &iy,
                            const int &dim)
{
  if(dim==2)
    return 0;

  const Gamra::Dir max=std::numeric_limits<Gamra::Dir>::max();
  const Gamra::Dir
    index[3][3]={{max,0,2},
                 {2,max,0},
                 {0,2,max}};
  return index[ix][iy];
}

