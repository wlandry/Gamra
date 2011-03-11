#include "V_Refine_Patch_Strategy.h"
#include "set_V_boundary.h"

void
SAMRAI::solv::V_Refine_Patch_Strategy::preprocessRefine
(hier::Patch& ,
 const hier::Patch& coarse,
 const hier::Box& ,
 const hier::IntVector& )
{
  set_V_boundary(coarse,v_id,true);
}
