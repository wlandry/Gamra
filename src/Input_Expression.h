#ifndef GAMRA_INPUT_EXPRESSION_H
#define GAMRA_INPUT_EXPRESSION_H

#include <muParser.h>
#include <string>
#include <list>
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Pointer.h"

class Input_Expression
{
public:
  mu::Parser equation;
  SAMRAI::tbox::Array<double> data, xyz_min, xyz_max;
  SAMRAI::tbox::Array<int> ijk;

  int dim, slice;
  double coord[3];
  bool use_equation;

  static double* variable_factory(const char *, void *)
  {
    static std::list<double> variables;
    variables.push_back(0);
    return &variables.back();
  }

  Input_Expression() {}

  Input_Expression(const std::string &name,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database,
                   const SAMRAI::tbox::Dimension& dimension,
                   const int num_components=1,
                   const int Slice=-1):
    dim(dimension.getValue()), slice(Slice)
  {
    if(database->keyExists(name))
      {
        equation.DefineVar("x",&coord[0]);
        equation.DefineVar("y",&coord[1]);
        equation.DefineVar("z",&coord[2]);
        equation.SetVarFactory(variable_factory, NULL);
        equation.SetExpr(database->getString(name));
        use_equation=true;
      }
    else if(database->keyExists(name+"_data"))
      {
        ijk=database->getIntegerArray(name+"_ijk");
        xyz_min=database->getDoubleArray(name+"_coord_min");
        xyz_max=database->getDoubleArray(name+"_coord_max");
        data=database->getDoubleArray(name+"_data");
        check_array_sizes(name,num_components);
        use_equation=false;
      }
    else
      {
        TBOX_ERROR("Could not find an entry for " + name);
      }
  }

  Input_Expression(const Input_Expression &e): dim(e.dim)
  {
    *this=e;
  }

  Input_Expression& operator=(const Input_Expression &e)
  {
    equation=e.equation;
    data=e.data;
    xyz_min=e.xyz_min;
    xyz_max=e.xyz_max;
    ijk=e.ijk;
    use_equation=e.use_equation;
    dim=e.dim;
    slice=e.slice;

    if(use_equation)
      {
        equation.DefineVar("x",&coord[0]);
        equation.DefineVar("y",&coord[1]);
        equation.DefineVar("z",&coord[2]);
      }
    return *this;
  }

  /* A little utility routine to validate the sizes of input arrays */
  void check_array_sizes(const std::string &name, const int &num_components)
  {
    const int array_dim(slice==-1 ? dim : dim-1);
    if(ijk.size()!=array_dim)
      TBOX_ERROR("Bad number of elements in " << name << "_ijk.  Expected "
                 << array_dim << " but got " << ijk.size());
    if(xyz_min.size()!=array_dim)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_coord_min.  Expected "
                 << array_dim << " but got " << xyz_min.size());
    if(xyz_max.size()!=array_dim)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_coord_max.  Expected "
                 << array_dim << " but got " << xyz_max.size());
    int data_size(1);
    for(int d=0; d<dim; ++d)
      if(d!=slice)
        data_size*=ijk[d];
    
    if(data.size()!=data_size*num_components)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_data.  Expected "
                 << data_size*num_components << " but got " << data.size());
  }

  double eval(const double Coord[3])
  {
    double result;

    if(use_equation)
      {
        coord[0]=Coord[0];
        coord[1]=Coord[1];
        coord[2]=Coord[2];
        result=equation.Eval();
      }
    else
      {
        int index(0), factor(1), d(0);
        for(int dd=0; dd<dim; ++dd)
          {
            if(dd==slice)
              continue;
            int i=static_cast<int>(Coord[dd]*(ijk[d]-1)
                                   /(xyz_max[d]-xyz_min[d])+0.5);
            i=std::max(0,std::min(ijk[d]-1,i));
            index+=i*factor;
            factor*=ijk[d];
            ++d;
          }
        result=data[index];
      }
    return result;
  }

};

#endif
