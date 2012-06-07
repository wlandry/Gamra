#ifndef GAMRA_STOKES_SET_BOUNDARY_H
#define GAMRA_STOKES_SET_BOUNDARY_H

#include "SAMRAI/hier/Patch.h"

void Stokes_set_boundary(const SAMRAI::hier::Patch& patch,
                         const int &p_id, const int &v_id, const bool &rhs);

#endif
