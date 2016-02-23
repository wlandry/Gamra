#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"

bool Elastic::Coarse_Fine_Boundary_Refine::is_residual;
int Elastic::Coarse_Fine_Boundary_Refine::dv_diagonal_id(invalid_id),
  Elastic::Coarse_Fine_Boundary_Refine::dv_mixed_id(invalid_id),
  Elastic::Coarse_Fine_Boundary_Refine::level_set_id(invalid_id);

