#include "Elastic/FAC.h"
#include <muParser.h>

double Elastic::FAC::evaluate_moduli(mu::Parser moduli_equation[2],
                                     const double xyz[3], const int &m,
                                     const int &dim)
{
  double result;
  if(!moduli_expression[m].empty())
    {
      result=moduli_equation[m].Eval();
    }
  else
    {
      int ijk(0), factor(1);
      for(int d=0;d<dim;++d)
        {
          int i=static_cast<int>(xyz[d]*(moduli_ijk[m][d]-1)
                                 /(moduli_xyz_max[m][d]-moduli_xyz_min[m][d]));
          i=std::max(0,std::min(moduli_ijk[m][d]-1,i));
          ijk+=i*factor;
          factor*=moduli_ijk[m][d];
        }
      result=moduli[m][ijk];
    }
  return result;
}
