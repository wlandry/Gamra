#ifndef SET_V_BOUNDARY_H
#define SET_V_BOUNDARY_H

#include "SAMRAI/hier/Patch.h"

void set_boundary(const SAMRAI::hier::Patch& patch,
                  const int &p_id, const int &v_id, const bool &rhs);

#endif
