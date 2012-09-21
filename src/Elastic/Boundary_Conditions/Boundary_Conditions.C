#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Database.h"
#include "Elastic/Boundary_Conditions.h"
#include <string>

#include "SAMRAI/tbox/PIO.h"

Elastic::Boundary_Conditions::Boundary_Conditions
(const SAMRAI::tbox::Dimension& dimension,
 const std::string& object_name,
 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database):
  d_object_name(object_name)
{
  const int dim(dimension.getValue());
  std::string v("v");
  std::string xyz("xyz");
  std::string upper_string[]={"lower","upper"};

  SAMRAI::tbox::Array<std::string> keys(database->getAllKeys());
  SAMRAI::tbox::plog << "keys\n";
  for(int i=0;i<keys.size();++i)
    SAMRAI::tbox::plog << keys[i] << "\n";

  for(int vxyz=0;vxyz<dim;++vxyz)
    for(int direction=0;direction<dim;++direction)
      for(int upper_lower=0;upper_lower<2;++upper_lower)
        {
          std::string bc_name(v+xyz[vxyz]+"_"+xyz[direction]+"_"
                              +upper_string[upper_lower]);
          if(database->keyExists("dirichlet_"+bc_name))
            {
              is_dirichlet[direction][upper_lower][vxyz]=true;
              dirichlet[direction][upper_lower][vxyz].DefineVar("x",&coord[0]);
              dirichlet[direction][upper_lower][vxyz].DefineVar("y",&coord[1]);
              dirichlet[direction][upper_lower][vxyz].DefineVar("z",&coord[2]);
              dirichlet[direction][upper_lower][vxyz].SetVarFactory(variable_factory, NULL);
              dirichlet[direction][upper_lower][vxyz].SetExpr(database->getString("dirichlet_"+bc_name));
            }
          else if(database->keyExists("neumann_"+bc_name))
            {
              is_dirichlet[direction][upper_lower][vxyz]=false;
              neumann[direction][upper_lower][vxyz].DefineVar("x",&coord[0]);
              neumann[direction][upper_lower][vxyz].DefineVar("y",&coord[1]);
              neumann[direction][upper_lower][vxyz].DefineVar("z",&coord[2]);
              neumann[direction][upper_lower][vxyz].SetVarFactory(variable_factory, NULL);
              neumann[direction][upper_lower][vxyz].SetExpr(database->getString("neumann_"+bc_name));
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
