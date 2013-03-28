#ifndef GAMRA_CONSTANTS_H
#define GAMRA_CONSTANTS_H

#include <limits>
const double boundary_value=1.1111e100;
const double invalid_value=3.1415926*boundary_value;
const int invalid_id=-1;
inline int index_map(const int &ix, const int &iy, const int &dim)
{
  if(dim==2)
    return 0;

  const int index[3][3]={{std::numeric_limits<int>::max(),0,2},
                         {2,std::numeric_limits<int>::max(),0},
                         {0,2,std::numeric_limits<int>::max()}};
  return index[ix][iy];
}
#endif
