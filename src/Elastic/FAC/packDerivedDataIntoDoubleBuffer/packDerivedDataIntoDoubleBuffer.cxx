/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"

void pack_strain(double* buffer,
                 const SAMRAI::hier::Patch& patch,
                 const SAMRAI::hier::Box& region,
                 const int &depth,
                 const SAMRAI::tbox::Dimension &d_dim,
                 const bool &have_faults,
                 const bool &have_embedded_boundary,
                 const int &v_id,
                 const int &dv_diagonal_id,
                 const int &dv_mixed_id,
                 const int &level_set_id);

void pack_level_set(double* buffer,
                    const SAMRAI::hier::Patch& patch,
                    const SAMRAI::hier::Box& region,
                    const SAMRAI::tbox::Dimension &d_dim,
                    const int &level_set_id);

void pack_v_v_rhs(double* buffer,
                  const SAMRAI::hier::Patch& patch,
                  const SAMRAI::hier::Box& region,
                  const std::string& variable_name,
                  const int &depth,
                  const SAMRAI::tbox::Dimension &d_dim,
                  const bool &offset_vector_on_output,
                  const int &v_id,
                  const int &v_rhs_id);

void pack_v_initial(double* buffer,
                    const SAMRAI::hier::Patch& patch,
                    const SAMRAI::hier::Box& region,
                    const int &depth,
                    const SAMRAI::tbox::Dimension &d_dim,
                    const bool &offset_vector_on_output,
                    const Input_Expression v_initial[]);

bool Elastic::FAC::packDerivedDataIntoDoubleBuffer
(double* buffer,
 const SAMRAI::hier::Patch& patch,
 const SAMRAI::hier::Box& region,
 const std::string& variable_name,
 int depth, double) const
{
  if(variable_name=="Strain")
    {
      pack_strain(buffer,patch,region,depth,d_dim,!faults.empty(),
                  have_embedded_boundary(),v_id,dv_diagonal_id,dv_mixed_id,
                  level_set_id);
    }
  else if(variable_name=="Level Set")
    {
      pack_level_set(buffer,patch,region,d_dim,level_set_id);
    }
  else if(variable_name=="Initial Displacement")
    {
      pack_v_initial(buffer,patch,region,depth,d_dim,offset_vector_on_output,
                     v_initial);
    }
  else
    {
      pack_v_v_rhs(buffer,patch,region,variable_name,depth,d_dim,
                   offset_vector_on_output,v_id,v_rhs_id);
    }
  // Always return true, since every patch has derived data.
  return true;
}
