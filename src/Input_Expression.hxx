/* Class to parse equations and interpolate arrays in input files. */

#pragma once

#include <muParser.h>
#include <string>
#include <list>
#include <SAMRAI/tbox/Array.h>
#include <SAMRAI/tbox/Database.h>
#include <Okada.hxx>
#include "Patch.hxx"

class Input_Expression
{
public:
  mu::Parser equation;
  std::vector<Patch> patches;

  int dim, slice;
  mutable double coord[3];
  bool use_equation, is_valid;

  static double* variable_factory(const char *, void *)
  {
    static std::list<double> variables;
    variables.push_back(0);
    return &variables.back();
  }

  static double okada(const double* arg_array, const int n)
  {
    std::vector<double> args(arg_array,arg_array+n);
    if(args.size()!=15)
      TBOX_ERROR("Wrong number of arguments for okada.  Expected 15, got "
                 << args.size() << "\n");
      
    int index(static_cast<int>(args[14]));
    if(index<0 || index>2)
      TBOX_ERROR("Bad index for okada.  Expected 0, 1, or 2, got "
                 << args[14] << "\n");

    /* If trying to get displacement above the surface, then compute
       the derivative and use that to extend the solution.  This works
       since we only need the displacement above the surface in order
       to compute a numerical derivative at the surface. */
    double z(args[13]), result;
    if(z<0)
      {
        std::vector<double> diff_args(args);
        args[13]=-z;
        diff_args[13]=0;
        result=okada_internal(args).first(index)
          + 2*z*okada_internal(diff_args).second(index,2);
      }
    else
      {
        result=okada_internal(args).first(index);
      }
    return result;
  }

  static double d_okada(const double* arg_array, const int n)
  {
    std::vector<double> args(arg_array,arg_array+n);
    if(args.size()!=16)
      TBOX_ERROR("Wrong number of arguments for d_okada.  Expected 16, got "
                 << args.size() << "\n");
      
    int index0(static_cast<int>(args[14])), index1(static_cast<int>(args[15]));
    if(index0<0 || index0>2 || index1<0 || index1>2)
      TBOX_ERROR("Bad index for d_okada.  Expected 0, 1, or 2, got "
                 << args[14] << " and " << args[15] << "\n");

    return okada_internal(args).second(index0,index1);
  }

  static Okada::Displacement okada_internal(const std::vector<double> &args)
  {
    FTensor::Tensor1<double,3> origin(args[3],args[4],args[5]),
      coord(args[11],args[12],args[13]);
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
    dim(dimension.getValue()), slice(Slice)
  {
    Init(name,database,num_components,error_if_missing);
  }

  Input_Expression(const std::string &name,
                   boost::shared_ptr<SAMRAI::tbox::Database> database,
                   const SAMRAI::tbox::Dimension& dimension,
                   const int num_components=1,
                   const int Slice=-1):
    dim(dimension.getValue()), slice(Slice)
  {
    Init(name,database,num_components,false);
  }

  /// Make this private so that I do not accidentally use it and not
  /// set dim or slice.
private:
  void Init(const std::string &name,
            boost::shared_ptr<SAMRAI::tbox::Database> database,
            const int num_components,
            const bool &error_if_missing)
  {
    is_valid=true;
    if(database->keyExists(name))
      {
        equation.DefineVar("x",&coord[0]);
        equation.DefineVar("y",&coord[1]);
        equation.DefineVar("z",&coord[2]);
        equation.DefineFun("okada",okada);
        equation.DefineFun("d_okada",d_okada);
        equation.SetVarFactory(variable_factory, NULL);
        equation.SetExpr(database->getString(name));
        use_equation=true;
      }
    else if(database->keyExists(name+"_data"))
      {
        patches.push_back(Patch(database,name,num_components,dim,slice));
        use_equation=false;
      }
    else if(database->keyExists(name+"_patches"))
      {
        if(!database->isDatabase(name+"_patches"))
          TBOX_ERROR("The entry for " + name + "_patches must be a struct.");

        boost::shared_ptr<SAMRAI::tbox::Database>
          patches_database(database->getDatabase(name+"_patches"));

        /// There is no way to get an iterator.  You have to get all
        /// of the keys and look up each one.
        std::vector<std::string>
          patch_names(database->getDatabase(name+"_patches")->getAllKeys());
        for(std::vector<std::string>::iterator patch=patch_names.begin();
            patch!=patch_names.end(); ++patch)
          {
            if(!patches_database->isDatabase(*patch))
              TBOX_ERROR("The entry for " + *patch + " in " + name
                         + "_patches must be a struct.");
            patches.push_back(Patch(patches_database->getDatabase(*patch),
                                    num_components,dim,slice));
          }
        use_equation=false;
      }
    else
      {
        is_valid=false;
        if(error_if_missing)
          TBOX_ERROR("Could not find an entry for " + name);
      }
  }

public:
  Input_Expression(const Input_Expression &e): dim(e.dim)
  {
    *this=e;
  }

  Input_Expression& operator=(const Input_Expression &e)
  {
    is_valid=e.is_valid;
    if(is_valid)
      {
        equation=e.equation;
        patches=e.patches;
        use_equation=e.use_equation;
        dim=e.dim;
        slice=e.slice;

        if(use_equation)
          {
            equation.DefineVar("x",&coord[0]);
            equation.DefineVar("y",&coord[1]);
            equation.DefineVar("z",&coord[2]);
          }
      }
    return *this;
  }

  double eval(const double Coord[3]) const
  {
    if(!is_valid)
      TBOX_ERROR("INTERNAL ERROR: Tried to use an invalid Input_Expression");

    double result;
    if(use_equation)
      {
        for(int d=0; d<dim; ++d)
          coord[d]=Coord[d];
        result=equation.Eval();
      }
    else
      {
        /// If there is only a single patch, then allow extrapolation
        /// from the edge of the patch.
        if (patches.size()==1)
          {
            result=patches[0].eval(Coord);
          }
        else
          {
            std::vector<Patch>::const_iterator patch=patches.begin();
            for (; patch!=patches.end(); ++patch)
              {
                if (patch->contains(Coord))
                  {
                    result=patch->eval(Coord);
                    break;
                  }
              }
            if (patch==patches.end())
              {
                std::stringstream ss;
                ss << "None of the patches contain the point ("
                   << Coord[0];
                for (int d=1; d<dim; ++d)
                  ss << ", " << Coord[d];
                ss << ")";
                TBOX_ERROR(ss.str());
                /// Add an abort() to remove the compiler warning
                /// about using uninitialized variables.
                abort();
              }
          }
      }
    return result;
  }

};


