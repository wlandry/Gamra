#ifndef GAMRA_ELASTIC_BOUNDARY_CONDITIONS_H
#define GAMRA_ELASTIC_BOUNDARY_CONDITIONS_H

#include <muParser.h>
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/hier/Patch.h"
#include <string>
#include <vector>

namespace Elastic {
  class Boundary_Conditions
  {
  public:
    Boundary_Conditions(const SAMRAI::tbox::Dimension& dim,
                        const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);
    void set_boundary(const SAMRAI::hier::Patch& patch,
                      const int &v_id, const bool &rhs);

    static double* variable_factory(const char *, void *)
    {
      static std::list<double> variables;
      variables.push_back(0);
      return &variables.back();
    }
    mu::Parser dirichlet[3][2][3], neumann[3][2][3];
    bool is_dirichlet[3][2][3];
    double coord[3];
    std::string d_object_name;
  };
}

#endif
