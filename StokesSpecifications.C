/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Stokes equation 
 *
 ************************************************************************/
#include "StokesSpecifications.h"

#ifndef SAMRAI_INLINE
#include "StokesSpecifications.I"
#endif

namespace SAMRAI {
namespace solv {

void StokesSpecifications::printClassData(
   std::ostream& stream) const
{
   stream << "StokesSpecifications " << d_object_name << "\n"
          << "   D is ";
   if (d_D_id != -1) {
      stream << "variable with patch id " << d_D_id << "\n";
   } else {
      stream << "constant with value " << d_D_constant << "\n";
   }
   stream << "   C is ";
   if (d_C_zero) {
      stream << "zero\n";
   } else if (d_C_id != -1) {
      stream << "variable with patch id " << d_C_id << "\n";
   } else {
      stream << "constant with value " << d_C_constant << "\n";
   }
}

} // namespace solv
} // namespace SAMRAI
