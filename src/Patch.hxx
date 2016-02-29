#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <vector>
#include <SAMRAI/tbox/Database.h>


class Patch
{
public:
  std::vector<double> data, xyz_min, xyz_max;
  std::vector<int> ijk;
  int dim, slice;
  
  Patch(SAMRAI::tbox::Database &database,
        const int &num_components, const int &Dim, const int &Slice) :
    data(database.getDoubleVector("data")),
    xyz_min(database.getDoubleVector("x_lo")),
    xyz_max(database.getDoubleVector("x_up")),
    ijk(database.getIntegerVector("ijk")),
    dim(Dim), slice(Slice)
  {
    check_array_sizes(database.getName() + "::",num_components);
  }

  Patch(SAMRAI::tbox::Database &database,
        const std::string &name, const int &num_components,
        const int &Dim, const int &Slice) :
    data(database.getDoubleVector(name+"_data")),
    ijk(database.getIntegerVector(name+"_ijk")),
    dim(Dim), slice(Slice)
  {
    if(database.keyExists(name+"_coord_min"))
      {
        xyz_min=database.getDoubleVector(name+"_coord_min");
        xyz_max=database.getDoubleVector(name+"_coord_max");
      }
    else
      {
        xyz_min=database.getDoubleVector(name+"_x_lo");
        xyz_max=database.getDoubleVector(name+"_x_up");
      }
    check_array_sizes(name+"_",num_components);
  }

  
  /// A little utility routine to validate the sizes of input arrays
  void check_array_sizes(const std::string &name, const int &num_components)
  {
    const size_t array_dim(slice==-1 ? dim : dim-1);
    if(ijk.size()!=array_dim)
      { TBOX_ERROR("Bad number of elements in " << name << "ijk.  Expected "
                   << array_dim << " but got " << ijk.size()); }
    if(xyz_min.size()!=array_dim)
      { TBOX_ERROR("Bad number of elements in "
                   << name << "coord_min or "
                   << name << "x_lo.  Expected "
                   << array_dim << " but got " << xyz_min.size()); }
    if(xyz_max.size()!=array_dim)
      { TBOX_ERROR("Bad number of elements in "
                   << name << "coord_max or "
                   << name << "x_up.  Expected "
                   << array_dim << " but got " << xyz_max.size()); }
    size_t data_size(1);
    for(int d=0; d<dim; ++d)
      { if(d!=slice)
          { data_size*=ijk[d]; } }
    
    if(data.size()!=data_size*num_components)
      { TBOX_ERROR("Bad number of elements in "
                   << name << "data.  Expected "
                   << data_size*num_components << " but got " << data.size()); }
  }

  bool contains(const double Coord[3]) const
  {
    bool result=true;
    int d(0);
    for (int dd=0; dd<dim; ++dd)
      {
        if (dd!=slice)
          {
            result=result && (Coord[dd]>=xyz_min[d]) && (Coord[dd]<=xyz_max[d]);
            ++d;
          }
      }
    return result;
  }
  
  double eval(const double Coord[3]) const
  {
    int d(0);
    int ix[3], ixp[3];
    double dx[3];
    for(int dd=0; dd<dim; ++dd)
      {
        if(dd==slice)
          { continue; }

        /* Use max(1,ijk-1) rather than (ijk-1) to handle the case
           when the array is only one element wide */
        double delta=(xyz_max[d]-xyz_min[d])/std::max(1,ijk[d]-1);
        ix[d]=static_cast<int>((Coord[dd]-xyz_min[d])/delta);
        ix[d]=std::max(0,std::min(ijk[d]-2,ix[d]));
        ixp[d]=std::max(0,std::min(ijk[d]-1,ix[d]+1));

        if(ijk[d]==1)
          { dx[d]=0; }
        else
          { dx[d]=std::min(1.,std::max(0.,(Coord[dd]-xyz_min[d]-ix[d]*delta)
                                       /delta)); }
        ++d;
      }
    double result;
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
        result=data[ix[0] + ijk[0]*(ix[1]+ ijk[1]*ix[2])]
          *(1-dx[0])*(1-dx[1])*(1-dx[2])
          + data[ixp[0] + ijk[0]*(ix[1]  + ijk[1]*ix[2])]
          *dx[0]    *(1-dx[1])*(1-dx[2])
          + data[ix[0]  + ijk[0]*(ixp[1] + ijk[1]*ix[2])]
          *(1-dx[0])*dx[1]    *(1-dx[2])
          + data[ixp[0] + ijk[0]*(ixp[1] + ijk[1]*ix[2])]
          *dx[0]    *dx[1]    *(1-dx[2])
          + data[ix[0]  + ijk[0]*(ix[1]  + ijk[1]*ixp[2])]
          *(1-dx[0])*(1-dx[1])*dx[2]
          + data[ixp[0] + ijk[0]*(ix[1]  + ijk[1]*ixp[2])]
          *dx[0]    *(1-dx[1])*dx[2]
          + data[ix[0]  + ijk[0]*(ixp[1] + ijk[1]*ixp[2])]
          *(1-dx[0])*dx[1]    *dx[2]
          + data[ixp[0] + ijk[0]*(ixp[1] + ijk[1]*ixp[2])]
          *dx[0]    *dx[1]    *dx[2];
        break;
      default:
        abort();
      }
    return result;
  }
};
