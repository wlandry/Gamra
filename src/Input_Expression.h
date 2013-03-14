#ifndef GAMRA_INPUT_EXPRESSION_H
#define GAMRA_INPUT_EXPRESSION_H

#include <muParser.h>
#include <string>
#include <list>
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

class Input_Expression
{
public:
  mu::Parser equation;
  SAMRAI::tbox::Array<double> data, xyz_min, xyz_max;
  SAMRAI::tbox::Array<int> ijk;

  int dim, slice;
  double coord[3];
  bool use_equation, is_valid;

  static double* variable_factory(const char *, void *)
  {
    static std::list<double> variables;
    variables.push_back(0);
    return &variables.back();
  }

  Input_Expression(): is_valid(false) {}

  Input_Expression(const std::string &name,
                   boost::shared_ptr<SAMRAI::tbox::Database> database,
                   const SAMRAI::tbox::Dimension& dimension,
                   const bool &error_if_missing,
                   const int num_components=1,
                   const int Slice=-1):
    dim(dimension.getValue()), slice(Slice), is_valid(true)
  {
    Init(name,database,num_components,error_if_missing);
  }

  Input_Expression(const std::string &name,
                   boost::shared_ptr<SAMRAI::tbox::Database> database,
                   const SAMRAI::tbox::Dimension& dimension,
                   const int num_components=1,
                   const int Slice=-1):
    dim(dimension.getValue()), slice(Slice), is_valid(true)
  {
    Init(name,database,num_components,false);
  }

  void Init(const std::string &name,
            boost::shared_ptr<SAMRAI::tbox::Database> database,
            const int num_components,
            const bool &error_if_missing)
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
        is_valid=false;
        if(error_if_missing)
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
        int d(0);
        int ix[3], ixp[3];
        double dx[3];
        for(int dd=0; dd<dim; ++dd)
          {
            if(dd==slice)
              continue;

            /* Use max(1,ijk-1) rather than (ijk-1) to handle the case
               when the array is only one element wide */
            double delta=(xyz_max[d]-xyz_min[d])/std::max(1,ijk[d]-1);
            ix[d]=static_cast<int>((Coord[dd]-xyz_min[d])/delta);
            ix[d]=std::max(0,std::min(ijk[d]-2,ix[d]));
            ixp[d]=std::max(0,std::min(ijk[d]-1,ix[d]+1));

            if(ijk[d]==1)
              {
                dx[d]=0;
              }
            else
              {
                dx[d]=(Coord[dd]-xyz_min[d]-ix[d]*delta)/delta;
              }
            ++d;
          }

        switch (d)
          {
          case 1:
            result=data[ix[0]]*(1-dx[0]) + data[ixp[0]]*dx[0];
            break;
          case 2:
            result=data[ix[0] + ix[1]*ijk[0]]*(1-dx[0]-dx[1]+dx[0]*dx[1])
              + data[ixp[0] + ix[1] *ijk[0]]*dx[0]*(1-dx[1])
              + data[ix[0]  + ixp[1]*ijk[0]]*dx[1]*(1-dx[0])
              + data[ixp[0] + ixp[1]*ijk[0]]*dx[0]*dx[1];
            break;
          case 3:
            result=data[ix[0] + ijk[0]*(ix[1]+ ijk[1]*ix[2])] *(1-dx[0])*(1-dx[1])*(1-dx[2])
              + data[ixp[0] + ijk[0]*(ix[1]  + ijk[1]*ix[2])] *dx[0]    *(1-dx[1])*(1-dx[2])
              + data[ix[0]  + ijk[0]*(ixp[1] + ijk[1]*ix[2])] *(1-dx[0])*dx[1]    *(1-dx[2])
              + data[ixp[0] + ijk[0]*(ixp[1] + ijk[1]*ix[2])] *dx[0]    *dx[1]    *(1-dx[2])
              + data[ix[0]  + ijk[0]*(ix[1]  + ijk[1]*ixp[2])]*(1-dx[0])*(1-dx[1])*dx[2]
              + data[ixp[0] + ijk[0]*(ix[1]  + ijk[1]*ixp[2])]*dx[0]    *(1-dx[1])*dx[2]
              + data[ix[0]  + ijk[0]*(ixp[1] + ijk[1]*ixp[2])]*(1-dx[0])*dx[1]    *dx[2]
              + data[ixp[0] + ijk[0]*(ixp[1] + ijk[1]*ixp[2])]*dx[0]    *dx[1]    *dx[2];
            break;
          default:
            abort();
          }
      }
    return result;
  }

};

#endif
