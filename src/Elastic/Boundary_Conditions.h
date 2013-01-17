#ifndef GAMRA_ELASTIC_BOUNDARY_CONDITIONS_H
#define GAMRA_ELASTIC_BOUNDARY_CONDITIONS_H

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include <string>
#include <vector>
#include "Input_Expression.h"
#include "edge_node_eval.h"

namespace Elastic {
  class Boundary_Conditions
  {
  public:
    Boundary_Conditions(const SAMRAI::tbox::Dimension& dim,
                        const std::string& object_name,
                        boost::shared_ptr<SAMRAI::tbox::Database> database);
    void set_boundary(const SAMRAI::hier::Patch& patch,
                      const int &v_id, const bool &rhs);
    void set_dirichlet(const SAMRAI::hier::Patch& patch,
                       const int &v_id, const bool &rhs);
    void set_extra_ids(const int &e_id, const int &dva, const int &dvp)
    {
      edge_moduli_id=e_id;
      dv_diagonal_id=dva;
      dv_perpendicular_id=dvp;
    }

    Input_Expression expression[3][3][2];
    bool is_dirichlet[3][3][2];
    double coord[3];
    std::string d_object_name;
    int edge_moduli_id,dv_diagonal_id,dv_perpendicular_id;

    void set_dirichlet
    (SAMRAI::pdat::SideData<double> &v,
     SAMRAI::hier::Index pp[], const int &dim,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::hier::Box &gbox,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
     const double *dx,
     const bool &homogeneous);

    void set_shear_derivs
    (SAMRAI::pdat::SideData<double> &v,
     SAMRAI::hier::Index pp[], const int &dim,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::hier::Box &gbox,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
     const double *dx,
     const bool &homogeneous);

    template<class T>
    void set_normal_stress
    (SAMRAI::pdat::SideData<double> &v,
     const T &edge_moduli,
     SAMRAI::hier::Index pp[], const int &dim,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::hier::Box &gbox,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
     const double *dx,
     const bool &homogeneous)
    {
      for(int ix=0; ix<dim; ++ix)
        {
          double offset[]={0.5,0.5,0.5};
          offset[ix]=0;
          SAMRAI::pdat::SideIterator send(gbox,ix,false);
          for(SAMRAI::pdat::SideIterator si(gbox,ix,true); si!=send; si++)
            {
              SAMRAI::pdat::SideIndex x(*si);
              
              bool on_corner(false);
              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  on_corner=on_corner
                    || !(x[iy]>=pbox.lower(iy) && x[iy]<=pbox.upper(iy));
                }
              
              if(!on_corner)
                {
                  for(int d=0;d<dim;++d)
                    coord[d]=geom->getXLower()[d]
                      + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

                  /* For normal BC's, for the point just outside the
                     boundary, set a sentinel value for normal dirichlet
                     BC or the derivative for normal stress BC. */
                  if(x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0))
                    {
                      if(!is_dirichlet[ix][ix][0])
                        {
                          double duyy=0;
                          for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                            {
                              SAMRAI::pdat::SideIndex
                                y(x,iy,SAMRAI::pdat::SideIndex::Lower);
                              duyy+=(v(y+pp[iy]) + v(y+pp[ix]+pp[iy])
                                     - v(y+pp[ix]) - v(y))
                                /(2*dx[ix]);
                            }
                          double lambda=
                            edge_node_average(edge_moduli,x+pp[ix],ix,pp,0);
                          double mu=
                            edge_node_average(edge_moduli,x+pp[ix],ix,pp,1);
                          v(x)=v(x+pp[ix]*2)
                            + lambda*duyy*2*dx[ix]/(lambda+2*mu);
                        
                          if(!homogeneous)
                            {
                              double coord_save(geom->getXLower()[ix]);
                              std::swap(coord[ix],coord_save);
                              v(x)-=expression[ix][ix][0].eval(coord)*2*dx[ix]
                                /(lambda+2*mu);
                              std::swap(coord[ix],coord_save);
                            }
                        }
                    }
                  else if(x[ix]>pbox.upper(ix)+1
                          && geom->getTouchesRegularBoundary(ix,1))
                    {
                      if(!is_dirichlet[ix][ix][1])
                        {
                          double duyy=0;
                          for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                            {
                              SAMRAI::pdat::SideIndex
                                y(x-pp[ix],iy,SAMRAI::pdat::SideIndex::Lower);
                              duyy+=(v(y+pp[iy]) + v(y-pp[ix]+pp[iy])
                                     -v(y) - v(y-pp[ix]))
                                /(2*dx[iy]);
                            }
                          double lambda=
                            edge_node_average(edge_moduli,x-pp[ix],ix,pp,0);
                          double mu=
                            edge_node_average(edge_moduli,x-pp[ix],ix,pp,1);

                          v(x)=v(x-pp[ix]*2) - lambda*duyy*2*dx[ix]/(lambda+2*mu);
                        
                          if(!homogeneous)
                            {
                              double coord_save(geom->getXUpper()[ix]);
                              std::swap(coord[ix],coord_save);
                              v(x)+=expression[ix][ix][1].eval(coord)*2*dx[ix]
                                /(lambda+2*mu);
                              std::swap(coord[ix],coord_save);
                            }
                        }
                    }
                }
            }
        }
    }
  };
}

#endif
