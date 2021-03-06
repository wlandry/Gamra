/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Boundary_Conditions.hxx"

void Elastic::Boundary_Conditions::set_physical_boundary
(const SAMRAI::hier::Patch& patch, const int &v_id, const bool &homogeneous,
 const bool &apply_normal_stress) const
{
  try
    {
      if(have_embedded_boundary())
        {
          set_embedded_boundary(patch,level_set_id,dv_mixed_id);
        }
      else
        {
          set_regular_boundary(patch,v_id,homogeneous,apply_normal_stress);
        }
    }
  catch(mu::Parser::exception_type &e)
    {
      TBOX_ERROR("Error in input formula\n"
                 << "Message:  " << e.GetMsg() << "\n"
                 << "Formula:  " << e.GetExpr() << "\n"
                 << "Token:    " << e.GetToken() << "\n"
                 << "Position: " << e.GetPos() << "\n"
                 << "Errc:     " << e.GetCode() << "\n");
    }
}
