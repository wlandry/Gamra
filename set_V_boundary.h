#ifndef SET_V_BOUNDARY_H
#define SET_V_BOUNDARY_H

#include "SAMRAI/hier/Patch.h"

void set_V_boundary(const SAMRAI::hier::Patch& patch, const int &id,
                    const bool &rhs);

#endif
