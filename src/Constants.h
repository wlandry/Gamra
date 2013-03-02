#ifndef GAMRA_CONSTANTS_H
#define GAMRA_CONSTANTS_H

#include <limits>
const double boundary_value=1e100;
const int invalid_id=-1;
const int index_map[3][3]={{std::numeric_limits<int>::max(),0,2},
                           {2,std::numeric_limits<int>::max(),0},
                           {0,2,std::numeric_limits<int>::max()}};
#endif
