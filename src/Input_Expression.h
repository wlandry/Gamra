/* Class to parse equations and interpolate arrays in input files. */

#ifndef GAMRA_INPUT_EXPRESSION_H
#define GAMRA_INPUT_EXPRESSION_H

#include <muParser.h>
#include <string>
#include <list>
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "Okada.hxx"

class Input_Expression
{
public:
  mu::Parser equation;
  std::vector<double> data, xyz_min, xyz_max;
  std::vector<int> ijk;

  int dim, slice;
  mutable double coord[3];
  bool use_equation, is_valid;

  static double* variable_factory(const char *, void *)
  {
    static std::list<double> variables;
    variables.push_back(0);
    return &variables.back();
  }

  static double okada(const double* args, const int n)
  {
    if(n!=15)
      TBOX_ERROR("Wrong number of arguments for okada.  Expected 15, got "
                 << n << "\n");
      
    int index(static_cast<int>(args[14]));
    if(index<0 || index>2)
      TBOX_ERROR("Bad index for okada.  Expected 0, 1, or 2, got "
                 << args[14] << "\n");

    return okada_internal(args).first(index);
  }

  static double d_okada(const double* args, const int n)
  {
    if(n!=16)
      TBOX_ERROR("Wrong number of arguments for d_okada.  Expected 16, got "
                 << n << "\n");
      
    int index0(static_cast<int>(args[14])), index1(static_cast<int>(args[15]));
    if(index0<0 || index0>2 || index1<0 || index1>2)
      TBOX_ERROR("Bad index for d_okada.  Expected 0, 1, or 2, got "
                 << args[14] << " and " << args[15] << "\n");

    return okada_internal(args).second(index0,index1);
  }

  static Okada::Displacement okada_internal(const double* args)
  {
    FTensor::Tensor1<double,3> origin(args[3],args[4],args[5]),
      coord(args[12],args[13],args[14]);
    /* Opening=0 */
    return Okada(args[0],args[1],args[2],0.0,args[7],args[6],args[8],args[9],
                 args[10],origin).displacement(coord);
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
        equation.DefineFun("okada",okada);
        equation.SetVarFactory(variable_factory, NULL);
        equation.SetExpr(database->getString(name));
        use_equation=true;
      }
    else if(database->keyExists(name+"_data"))
      {
        ijk=database->getIntegerVector(name+"_ijk");
        xyz_min=database->getDoubleVector(name+"_coord_min");
        xyz_max=database->getDoubleVector(name+"_coord_max");
        data=database->getDoubleVector(name+"_data");
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
    const size_t array_dim(slice==-1 ? dim : dim-1);
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
    size_t data_size(1);
    for(int d=0; d<dim; ++d)
      if(d!=slice)
        data_size*=ijk[d];
    
    if(data.size()!=data_size*num_components)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_data.  Expected "
                 << data_size*num_components << " but got " << data.size());
  }

  double eval(const double Coord[3]) const
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
                dx[d]=std::min(1.,std::max(0.,(Coord[dd]-xyz_min[d]-ix[d]*delta)
                                           /delta));
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
