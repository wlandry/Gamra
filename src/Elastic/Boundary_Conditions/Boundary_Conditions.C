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
          expression[vxyz][direction][upper_lower]=
            Input_Expression(bc_name,database,dimension);
          is_dirichlet[vxyz][direction][upper_lower]=
            database->getBool(bc_name + "_is_dirichlet");
        }
}
