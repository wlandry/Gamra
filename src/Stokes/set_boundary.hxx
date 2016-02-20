#pragma once

#include <SAMRAI/hier/Patch.h>

void Stokes_set_boundary(const SAMRAI::hier::Patch& patch,
                         const int &p_id, const int &v_id, const bool &rhs);


