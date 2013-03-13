#include "Elastic/V_Boundary_Refine.h"

bool Elastic::V_Boundary_Refine::is_residual;
int Elastic::V_Boundary_Refine::dv_diagonal_id(invalid_id),
  Elastic::V_Boundary_Refine::dv_mixed_id(invalid_id),
  Elastic::V_Boundary_Refine::level_set_id(invalid_id);

