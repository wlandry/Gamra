/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Boundary_Conditions.hxx"

Elastic::Boundary_Conditions::Boundary_Conditions
(const SAMRAI::tbox::Dimension& dimension,
 const std::string& object_name,
 SAMRAI::tbox::Database &database):
  d_object_name(object_name), edge_moduli_id(invalid_id),
  dv_diagonal_id(invalid_id), dv_mixed_id(invalid_id),
  level_set_id(invalid_id)
{
  const int dim(dimension.getValue());
  std::string v("v");
  std::string xyz("xyz");
  std::string upper_string[]={"lower","upper"};

  for(int vxyz=0;vxyz<dim;++vxyz)
    for(int direction=0;direction<dim;++direction)
      for(int upper_lower=0;upper_lower<2;++upper_lower)
        {
          std::string bc_name(v+xyz[vxyz]+"_"+xyz[direction]+"_"
                              +upper_string[upper_lower]);
          expression[vxyz][direction][upper_lower]=
            Input_Expression(bc_name,database,dimension,1,direction);
          is_dirichlet[vxyz][direction][upper_lower]=
            database.getBool(bc_name + "_is_dirichlet");
        }
}
