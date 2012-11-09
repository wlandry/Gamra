#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Database.h"
#include "Elastic/Boundary_Conditions.h"
#include "Constants.h"
#include <string>

#include "SAMRAI/tbox/PIO.h"

Elastic::Boundary_Conditions::Boundary_Conditions
(const SAMRAI::tbox::Dimension& dimension,
 const std::string& object_name,
 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database):
  d_object_name(object_name), edge_moduli_id(invalid_id)
{
  const int dim(dimension.getValue());
  std::string v("v");
  std::string xyz("xyz");
  std::string upper_string[]={"lower","upper"};

  SAMRAI::tbox::Array<std::string> keys(database->getAllKeys());

  for(int vxyz=0;vxyz<dim;++vxyz)
    for(int direction=0;direction<dim;++direction)
      for(int upper_lower=0;upper_lower<2;++upper_lower)
        {
          std::string bc_name(v+xyz[vxyz]+"_"+xyz[direction]+"_"
                              +upper_string[upper_lower]);
          if(database->keyExists("dirichlet_"+bc_name))
            {
              is_dirichlet[vxyz][direction][upper_lower]=true;
              dirichlet[vxyz][direction][upper_lower].DefineVar("x",&coord[0]);
              dirichlet[vxyz][direction][upper_lower].DefineVar("y",&coord[1]);
              dirichlet[vxyz][direction][upper_lower].DefineVar("z",&coord[2]);
              dirichlet[vxyz][direction][upper_lower].SetVarFactory(variable_factory, NULL);
              dirichlet[vxyz][direction][upper_lower].SetExpr(database->getString("dirichlet_"+bc_name));
            }
          else if(database->keyExists("normal_stress_"+bc_name))
            {
              if(vxyz!=direction)
                TBOX_ERROR(d_object_name
                           << ": normal_stress boundaries must be normal (e.g. vx_x, not vx_y) '");
              is_dirichlet[vxyz][direction][upper_lower]=false;
              normal_stress[vxyz][upper_lower].DefineVar("x",&coord[0]);
              normal_stress[vxyz][upper_lower].DefineVar("y",&coord[1]);
              normal_stress[vxyz][upper_lower].DefineVar("z",&coord[2]);
              normal_stress[vxyz][upper_lower].SetVarFactory(variable_factory, NULL);
              normal_stress[vxyz][upper_lower].SetExpr(database->getString("normal_stress_"+bc_name));
            }
          else if(database->keyExists("shear_deriv_"+bc_name))
            {
              if(vxyz==direction)
                TBOX_ERROR(d_object_name
                           << ": shear_deriv boundaries must be mixed (e.g. vx_y, not vx_x) '");
              is_dirichlet[vxyz][direction][upper_lower]=false;
              shear_derivs[vxyz][direction][upper_lower].DefineVar("x",&coord[0]);
              shear_derivs[vxyz][direction][upper_lower].DefineVar("y",&coord[1]);
              shear_derivs[vxyz][direction][upper_lower].DefineVar("z",&coord[2]);
              shear_derivs[vxyz][direction][upper_lower].SetVarFactory(variable_factory, NULL);
              shear_derivs[vxyz][direction][upper_lower].SetExpr(database->getString("shear_deriv_"+bc_name));
            }
          else
            {
              TBOX_ERROR(d_object_name
                         << ": Can not find a boundary condition for '"
                         << bc_name
                         << "' in Elastic::FACOps::setSmoothingChoice.");
            }
        }
}
