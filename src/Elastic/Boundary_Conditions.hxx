#ifndef GAMRA_ELASTIC_BOUNDARY_CONDITIONS_H
#define GAMRA_ELASTIC_BOUNDARY_CONDITIONS_H

#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/pdat/NodeData.h>
#include <SAMRAI/pdat/EdgeData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <string>
#include <vector>
#include "Input_Expression.hxx"
#include "edge_node_eval.hxx"
#include "Constants.hxx"

namespace Elastic {
  class Boundary_Conditions
  {
  public:
    Boundary_Conditions(const SAMRAI::tbox::Dimension& dim,
                        const std::string& object_name,
                        boost::shared_ptr<SAMRAI::tbox::Database> database);
    void set_boundary(const SAMRAI::hier::Patch& patch,
                      const int &v_id, const bool &homogeneous) const
    {
      set_boundary(patch,v_id,homogeneous,true);
    }
    void set_boundary(const SAMRAI::hier::Patch& patch,
                      const int &v_id, const bool &homogeneous,
                      const bool &apply_normal_stress)  const;
    void set_regular_boundary(const SAMRAI::hier::Patch& patch,
                              const int &v_id, const bool &homogeneous,
                              const bool &apply_normal_stress)  const;
    static void set_embedded_boundary(const SAMRAI::hier::Patch& patch,
                                      const int &level_id,
                                      const int &dv_mixed_id);

    void set_dirichlet(const SAMRAI::hier::Patch& patch,
                       const int &v_id, const bool &rhs) const;
    void set_extra_ids(const int &Edge_moduli_id, const int &Dv_diagonal_id,
                       const int &Dv_mixed_id, const int &Level_set_id)
    {
      edge_moduli_id=Edge_moduli_id;
      dv_diagonal_id=Dv_diagonal_id;
      dv_mixed_id=Dv_mixed_id;
      level_set_id=Level_set_id;
    }

    Input_Expression expression[3][3][2];
    bool is_dirichlet[3][3][2];
    std::string d_object_name;
    int edge_moduli_id, dv_diagonal_id, dv_mixed_id, level_set_id;

    bool have_faults() const
    {
      return dv_diagonal_id!=invalid_id;
    }

    bool have_embedded_boundary() const
    {
      return level_set_id!=invalid_id;
    }

    void set_dirichlet
    (SAMRAI::pdat::SideData<double> &v,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_ptr,
     const SAMRAI::hier::Index unit[],
     const Gamra::Dir &dim,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::hier::Box &gbox,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
     const double *dx,
     const bool &homogeneous) const;

    void set_shear_stress
    (SAMRAI::pdat::SideData<double> &v,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_ptr,
     const SAMRAI::hier::Index unit[],
     const Gamra::Dir &dim,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::hier::Box &gbox,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
     const double *dx,
     const bool &homogeneous) const;

    bool at_corner
    (const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> &geom,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::pdat::SideIndex &x,
     const Gamra::Dir &ix,
     const Gamra::Dir &iz) const
    {
      return (geom->getTouchesRegularBoundary(ix,0)
              && (x[ix]<pbox.lower(ix)
                  || (x[ix]<=pbox.lower(ix) && is_dirichlet[ix][ix][0])))
        || (geom->getTouchesRegularBoundary(ix,1)
            && (x[ix]>pbox.upper(ix)+1
                || (x[ix]>pbox.upper(ix) && is_dirichlet[ix][ix][1])))
        || (geom->getTouchesRegularBoundary(iz,0) && x[iz]<pbox.lower(iz))
        || (geom->getTouchesRegularBoundary(iz,1) && x[iz]>pbox.upper(iz));
    }

    template<class T>
    void set_normal_stress
    (SAMRAI::pdat::SideData<double> &v,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> > &dv_diagonal_ptr,
     const T &edge_moduli,
     const SAMRAI::hier::Index unit[],
     const Gamra::Dir &dim,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::hier::Box &gbox,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
     const double *dx,
     const bool &homogeneous) const
    {
      for(Gamra::Dir ix=0; ix<dim; ++ix)
        {
          double offset[]={0.5,0.5,0.5};
          offset[ix]=0;

          if(geom->getTouchesRegularBoundary(ix,0) && !is_dirichlet[ix][ix][0])
            {
              SAMRAI::hier::Box box(gbox);
              box.setUpper(ix,pbox.lower(ix)-1);

              SAMRAI::pdat::CellIterator
                end(SAMRAI::pdat::CellGeometry::end(box));
              for(SAMRAI::pdat::CellIterator
                    ci(SAMRAI::pdat::CellGeometry::begin(box)); ci!=end; ++ci)
                {
                  const SAMRAI::pdat::SideIndex
                    x(*ci,ix,SAMRAI::pdat::SideIndex::Lower);

                  bool on_corner(false);
                  for(Gamra::Dir iy=ix.next(dim);
                      iy!=ix; iy=iy.next(dim))
                    {
                      on_corner=on_corner
                        || !(x[iy]>=pbox.lower(iy) && x[iy]<=pbox.upper(iy));
                    }

                  if(!on_corner)
                    {
                      double coord[3];
                      for(int d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                      double duyy=0;
                      for(Gamra::Dir iy=ix.next(dim);
                          iy!=ix; iy=iy.next(dim))
                        {
                          SAMRAI::pdat::SideIndex
                            y(x,iy,SAMRAI::pdat::SideIndex::Lower);
                          duyy+=(v(y+unit[iy]) + v(y+unit[ix]+unit[iy])
                                 - v(y+unit[ix]) - v(y))
                            /(2*dx[iy]);
                          if(have_faults() && !homogeneous)
                            {
                              /* We only have to correct for one
                                 of the derivatives because the
                                 other half is outside the domain
                                 and defined to be regular. */
                              SAMRAI::pdat::CellIndex c(y+unit[ix]);
                              duyy-=(*dv_diagonal_ptr)(c,iy)/(2*dx[iy]);
                            }
                        }
                      double lambda=
                        edge_node_average(edge_moduli,x+unit[ix],ix,unit,0);
                      double mu=
                        edge_node_average(edge_moduli,x+unit[ix],ix,unit,1);
                      v(x)=v(x+unit[ix]*2)
                        + lambda*duyy*2*dx[ix]/(lambda+2*mu);
                        
                      if(!homogeneous)
                        {
                          SAMRAI::pdat::CellIndex c(x+unit[ix]);
                          double coord_save(geom->getXLower()[ix]);
                          std::swap(coord[ix],coord_save);
                          if(have_faults())
                            v(x)-=(*dv_diagonal_ptr)(c,ix);
                          v(x)-=expression[ix][ix][0].eval(coord)*2*dx[ix]
                            /(lambda+2*mu);
                          std::swap(coord[ix],coord_save);
                        }
                    }
                }
            }
          if(geom->getTouchesRegularBoundary(ix,1) && !is_dirichlet[ix][ix][1])
            {
              SAMRAI::hier::Box box(gbox);
              box.setLower(ix,pbox.upper(ix)+2);
              box.setUpper(ix,std::max(box.upper(ix),box.lower(ix)));

              SAMRAI::pdat::CellIterator
                end(SAMRAI::pdat::CellGeometry::end(box));
              for(SAMRAI::pdat::CellIterator
                    ci(SAMRAI::pdat::CellGeometry::begin(box)); ci!=end; ++ci)
                {
                  const SAMRAI::pdat::SideIndex
                    x(*ci,ix,SAMRAI::pdat::SideIndex::Lower);

                  bool on_corner(false);
                  for(Gamra::Dir iy=ix.next(dim); iy!=ix; iy=iy.next(dim))
                    {
                      on_corner=on_corner
                        || !(x[iy]>=pbox.lower(iy) && x[iy]<=pbox.upper(iy));
                    }

                  if(!on_corner)
                    {
                      double coord[3];
                      for(Gamra::Dir d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

                      double duyy=0;
                      for(Gamra::Dir iy=ix.next(dim); iy!=ix; iy=iy.next(dim))
                        {
                          SAMRAI::pdat::SideIndex
                            y(x-unit[ix],iy,SAMRAI::pdat::SideIndex::Lower);
                          duyy+=(v(y+unit[iy]) + v(y-unit[ix]+unit[iy])
                                 -v(y) - v(y-unit[ix]))
                            /(2*dx[iy]);
                          if(have_faults() && !homogeneous)
                            {
                              /* We only have to correct for one
                                 of the derivatives because the
                                 other half is outside the domain
                                 and defined to regular. */
                              SAMRAI::pdat::CellIndex c(y-unit[ix]);
                              duyy-=(*dv_diagonal_ptr)(c,iy)/(2*dx[iy]);
                            }
                        }
                      double lambda=
                        edge_node_average(edge_moduli,x-unit[ix],ix,unit,0);
                      double mu=
                        edge_node_average(edge_moduli,x-unit[ix],ix,unit,1);

                      v(x)=v(x-unit[ix]*2)
                        - lambda*duyy*2*dx[ix]/(lambda+2*mu);
                        
                      if(!homogeneous)
                        {
                          SAMRAI::pdat::CellIndex c(x-unit[ix]*2);
                          double coord_save(geom->getXUpper()[ix]);
                          std::swap(coord[ix],coord_save);
                          if(have_faults())
                            v(x)+=(*dv_diagonal_ptr)(c,ix);
                          v(x)+=expression[ix][ix][1].eval(coord)*2*dx[ix]
                            /(lambda+2*mu);
                          std::swap(coord[ix],coord_save);
                        }
                    }
                }
            }
        }
    }
  };
}

#endif
