#ifndef GAMRA_ELASTIC_SET_BOUNDARY_H
#define GAMRA_ELASTIC_SET_BOUNDARY_H

#include "SAMRAI/hier/Patch.h"

void set_boundary(const SAMRAI::hier::Patch& patch,
                  const int &v_id, const bool &rhs);

#endif
