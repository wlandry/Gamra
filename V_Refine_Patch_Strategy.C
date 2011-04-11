#include "V_Refine_Patch_Strategy.h"
#include "set_boundary.h"

void
SAMRAI::solv::V_Refine_Patch_Strategy::preprocessRefine
(hier::Patch& ,
 const hier::Patch& coarse,
 const hier::Box& ,
 const hier::IntVector& )
{
  set_boundary(coarse,-1,v_id,true);
}
