/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#ifndef GAMRA_ELASTIC_FAC_H
#define GAMRA_ELASTIC_FAC_H

#include "Elastic/FACSolver.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "Elastic/Boundary_Conditions.h"
#include "Input_Expression.h"
#include "FTensor.hpp"

namespace Elastic {
  /*!
   * @brief Class to solve a sample Elastic equation on a SAMR grid.
   */
  class FAC:
    public SAMRAI::mesh::StandardTagAndInitStrategy,
    public SAMRAI::appu::VisDerivedDataStrategy
  {

  public:
    /*!
     * @brief Constructor.
     *
     * If you want standard output and logging,
     * pass in valid pointers for those streams.
     *
     * @param object_name Ojbect name
     * @param database Input database (may be NULL)
     */
    FAC(const std::string& object_name,
        const SAMRAI::tbox::Dimension& dim,
        boost::shared_ptr<SAMRAI::tbox::Database> database =
        boost::shared_ptr<SAMRAI::tbox::Database>());

    virtual ~FAC() {}

    //@{ @name mesh::StandardTagAndInitStrategy virtuals

    /*!
     * @brief Allocate and initialize data for a new level
     * in the patch hierarchy.
     *
     * This is where you implement the code for initialize data on
     * the grid.  All the information needed to initialize the grid
     * are in the arguments.
     *
     * @see mesh::StandardTagAndInitStrategy::initializeLevelData()
     */
    virtual void
    initializeLevelData(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>&
                        hierarchy,
                        const int level_number,
                        const double init_data_time,
                        const bool can_be_refined,
                        const bool initial_time,
                        const boost::shared_ptr<SAMRAI::hier::PatchLevel>&
                        old_level,
                        const bool allocate_data);

    /*!
     * @brief Reset any internal hierarchy-dependent information.
     */
    virtual void
    resetHierarchyConfiguration(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>&
                                new_hierarchy,
                                int coarsest_level,
                                int finest_level);

    //@}

    virtual void
    applyGradientDetector(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                          const int level_number,
                          const double error_data_time,
                          const int tag_index,
                          const bool initial_time,
                          const bool uses_richardson_extrapolation);

    void
    computeAdaptionEstimate(SAMRAI::pdat::CellData<double>& estimate_data,
                            const SAMRAI::pdat::CellData<double>& soln_cell_data)
      const;

    //@{ @name appu::VisDerivedDataStrategy virtuals

    virtual bool
    packDerivedDataIntoDoubleBuffer(double* buffer,
                                    const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& region,
                                    const std::string& variable_name,
                                    int depth) const
    {
      if(variable_name=="Strain")
        return pack_strain(buffer,patch,region,depth);
      else if(variable_name=="Level Set")
        return pack_level_set(buffer,patch,region);
      else
        return pack_v_v_rhs(buffer,patch,region,variable_name,depth);
    }

    bool
    pack_strain(double* buffer,
                const SAMRAI::hier::Patch& patch,
                const SAMRAI::hier::Box& region,
                const int &depth) const;

    bool
    pack_level_set(double* buffer,
                   const SAMRAI::hier::Patch& patch,
                   const SAMRAI::hier::Box& region) const;

    bool
    pack_v_v_rhs(double* buffer,
                 const SAMRAI::hier::Patch& patch,
                 const SAMRAI::hier::Box& region,
                 const std::string& variable_name,
                 const int &depth) const;
    //@}

    /*!
     * @brief Solve using HYPRE Elastic solver
     *
     * Set up the linear algebra problem and use a
     * solv::Elastic::FACSolver object to solve it.
     * -# Set initial guess
     * -# Set boundary conditions
     * -# Specify Elastic equation parameters
     * -# Call solver
     */
    int
    solve();

#ifdef HAVE_HDF5
    /*!
     * @brief Set up external plotter to plot internal
     * data from this class.
     *
     * After calling this function, the external
     * data writer may be used to write the
     * viz file for this object.
     *
     * The internal hierarchy is used and must be
     * established before calling this function.
     * (This is commonly done by building a hierarchy
     * with the mesh::StandardTagAndInitStrategy virtual
     * functions implemented by this class.)
     *
     * @param viz_writer VisIt writer
     */
    int
    setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const;
#endif

  private:
    void fix_moduli();
    std::string d_object_name;

    const SAMRAI::tbox::Dimension d_dim;

    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;

    //@{
    /*!
     * @name Major algorithm objects.
     */
  public:
    Boundary_Conditions d_boundary_conditions;
  private:
    /*!
     * @brief FAC Elastic solver.
     */
    Elastic::FACSolver d_elastic_fac_solver;

    //@}

    //@{

    /*!
     * @name Private state variables for solution.
     */

    /*!
     * @brief Context owned by this object.
     */
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
  
    /*!
     * @brief Descriptor indices of internal data.
     *
     * These are initialized in the constructor and never change.
     */

    double d_adaption_threshold;
    int min_full_refinement_level;
  public:
    int cell_moduli_id, edge_moduli_id, v_id, v_rhs_id, dv_diagonal_id,
      dv_mixed_id, level_set_id;

    Input_Expression lambda, mu, v_rhs, level_set;

    SAMRAI::tbox::Array<double> faults;
    //@}

    bool have_embedded_boundary() const
    {
      return level_set.is_valid;
    }

    template<class T> void add_faults();

    bool intersect_fault(const int &dim,
                         const FTensor::Tensor1<double,3> &c0,
                         const FTensor::Tensor1<double,3> &c1,
                         const double fault[])
    {
      bool result(true);
      for(int d=1;d<dim;++d)
        {
          double y((c1(d)*c0(0) -  c1(0)*c0(d))/(c0(0) - c1(0)));
          result=result && ((y<=fault[d-1] && y>0) || (y>=fault[d-1] && y<0));
        }
      return result;
    }

    int intersection(const FTensor::Tensor1<double,3> &ntt,
                     const FTensor::Tensor1<double,3> &xyz,
                     const FTensor::Tensor2<double,3,3> &rot,
                     const FTensor::Tensor1<double,3> &dx,
                     const double fault[],
                     const int &dim);

    void compute_intersections_2D(const FTensor::Tensor1<double,3> &ntt,
                                  const FTensor::Tensor1<double,3> &xyz,
                                  const FTensor::Tensor2<double,3,3> &rot,
                                  const FTensor::Tensor1<double,3> dx[],
                                  const double fault[],
                                  const int &dim,
                                  const int &ix,
                                  int &intersect_diagonal,
                                  int intersect_mixed[2])
    {
      intersect_diagonal=intersection(ntt,xyz,rot,dx[ix],fault,dim);

      FTensor::Tensor1<double,3> dx_2;
      FTensor::Index<'a',3> a;
      const int iy((ix+1)%dim);
      dx_2(a)=dx[iy](a)/2;
      intersect_mixed[0]=intersection(ntt,xyz,rot,dx_2,fault,dim);

      dx_2(a)=-dx_2(a);
      intersect_mixed[1]=-intersection(ntt,xyz,rot,dx_2,fault,dim);
    }

    void compute_intersections_3D(const FTensor::Tensor1<double,3> &ntt,
                                  const FTensor::Tensor1<double,3> &xyz,
                                  const FTensor::Tensor2<double,3,3> &rot,
                                  const FTensor::Tensor1<double,3> dx[],
                                  const double fault[],
                                  const int &dim,
                                  const int &ix,
                                  int &intersect_diagonal,
                                  int intersect_mixed[4],
                                  int intersect_corner[4])
    {
      intersect_diagonal=intersection(ntt,xyz,rot,dx[ix],fault,dim);

      FTensor::Tensor1<double,3> dx_2_y, dx_2_z;
      FTensor::Index<'a',3> a;

      int iy((ix+1)%dim), iz((ix+2)%dim);
      
      dx_2_y(a)=dx[iy](a)/2;
      intersect_mixed[0]=intersection(ntt,xyz,rot,dx_2_y,fault,dim);
      dx_2_y(a)=-dx_2_y(a);
      intersect_mixed[1]=intersection(ntt,xyz,rot,dx_2_y,fault,dim);

      dx_2_z(a)=dx[iz](a)/2;
      intersect_mixed[2]=intersection(ntt,xyz,rot,dx_2_z,fault,dim);
      dx_2_z(a)=-dx_2_z(a);
      intersect_mixed[3]=intersection(ntt,xyz,rot,dx_2_z,fault,dim);

      FTensor::Tensor1<double,3> dx_corner;
      dx_corner(ix)=0;
      dx_corner(iy)=dx[iy](iy)/2;
      dx_corner(iz)=dx[iz](iz)/2;
      intersect_corner[0]=intersection(ntt,xyz,rot,dx_corner,fault,dim);
      dx_corner(iz)=-dx_corner(iz);
      intersect_corner[1]=intersection(ntt,xyz,rot,dx_corner,fault,dim);

      dx_corner(iy)=-dx_corner(iy);
      intersect_corner[2]=intersection(ntt,xyz,rot,dx_corner,fault,dim);
      dx_corner(iz)=-dx_corner(iz);
      intersect_corner[3]=intersection(ntt,xyz,rot,dx_corner,fault,dim);
    }
  };
}

/* Add corrections due to faults */

template<class T>
void Elastic::FAC::add_faults()
{
  const int max_level(d_hierarchy->getFinestLevelNumber());
  const int dim=d_dim.getValue();

  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(d_dim));
  SAMRAI::hier::Index unit[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    unit[i][i]=1;

  for(int l=0; l<=max_level; ++l)
    {
      boost::shared_ptr<SAMRAI::hier::PatchLevel>
        level(d_hierarchy->getPatchLevel(l));
      
      for(SAMRAI::hier::PatchLevel::Iterator p(level->begin());
          p!=level->end(); p++)
        {
          boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
            boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
            ((*p)->getPatchGeometry());
          const double *dx=geom->getDx();

          /* v_rhs */
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            ((*p)->getPatchData(v_rhs_id));
          SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);


          /* dv */
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            ((*p)->getPatchData(dv_diagonal_id));
          SAMRAI::pdat::CellData<double> &dv_diagonal(*dv_diagonal_ptr);
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            ((*p)->getPatchData(dv_mixed_id));
          SAMRAI::pdat::SideData<double> &dv_mixed(*dv_mixed_ptr);
          dv_diagonal.fillAll(0);
          dv_mixed.fillAll(0);

          /* moduli */
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            ((*p)->getPatchData(cell_moduli_id));
          SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);

          boost::shared_ptr<T> edge_moduli_ptr =
            boost::dynamic_pointer_cast<T>((*p)->getPatchData(edge_moduli_id));
          T &edge_moduli(*edge_moduli_ptr);

          SAMRAI::hier::Box pbox = v_rhs.getBox();
          SAMRAI::hier::Box gbox = v_rhs.getGhostBox();
          for(int ix=0;ix<dim;++ix)
            {
              if(geom->getTouchesRegularBoundary(ix,0))
                gbox.shorten(ix,-1);
              if(geom->getTouchesRegularBoundary(ix,1))
                gbox.shorten(ix,1);
            }

          /* Iterate over the faults */
          const int params_per_fault(9);
          for(int fault_index=0;fault_index<faults.size();
              fault_index+=params_per_fault)
            {
              const double pi=4*atan(1);
              double scale(-faults[fault_index+0]), x(faults[fault_index+1]),
                y(faults[fault_index+2]), z(faults[fault_index+3]),
                L(faults[fault_index+4]), W(faults[fault_index+5]),
                strike(pi/2-faults[fault_index+6]*pi/180),
                dip(faults[fault_index+7]*pi/180),
                rake(faults[fault_index+8]*pi/180);

              const FTensor::Tensor1<double,3> center(x,y,z);
              const FTensor::Tensor2<double,3,3>
                rot_strike(std::cos(strike),-std::sin(strike),0,
                           std::sin(strike),std::cos(strike),0,0,0,1),
                rot_dip(std::sin(dip),0,std::cos(dip),
                        0,1,0,
                        -std::cos(dip),0,std::sin(dip));

              FTensor::Tensor2<double,3,3> rot;
              FTensor::Index<'a',3> a;
              FTensor::Index<'b',3> b;
              FTensor::Index<'c',3> c;
              rot(a,b)=rot_dip(a,c)*rot_strike(c,b);

              FTensor::Tensor1<double,3> Dx[dim];
              for(int d0=0;d0<dim;++d0)
                {
                  Dx[d0](a)=0.0;
                  Dx[d0](d0)=dx[d0];
                }
              FTensor::Tensor1<double,3> slip(0,scale,0);
              const FTensor::Tensor2<double,3,3>
                rot_rake(1,0,0,
                         0,std::cos(rake),-std::sin(rake),
                         0,std::sin(rake),std::cos(rake));
              FTensor::Tensor1<double,3> jump;
              jump(c)=rot(b,c)*rot_rake(b,a)*slip(a);

              double fault[]={L,W};

              for(int ix=0;ix<dim;++ix)
                {
                  double offset[]={0.5,0.5,0.5};
                  offset[ix]=0;

                  SAMRAI::pdat::SideIterator s_end(gbox,ix,false);
                  for(SAMRAI::pdat::SideIterator si(gbox,ix,true); si!=s_end;
                      si++)
                    {
                      const SAMRAI::pdat::SideIndex s(*si);

                      FTensor::Tensor1<double,3> xyz(0,0,0);
                      for(int d=0;d<dim;++d)
                        xyz(d)=geom->getXLower()[d]
                          + dx[d]*(s[d]-pbox.lower()[d]+offset[d]) - center(d);

                      /* Rotate the coordinates into the coordinates of the
                         fault.  So in those coordinates, if x<0, you are on
                         the left, and if x>0, you are on the right. */
                      FTensor::Tensor1<double,3> ntt;
                      ntt(a)=rot(a,b)*xyz(b);
                      if(dim==2)
                        {
                          int intersect, intersect_mixed[2];
                          compute_intersections_2D(ntt,xyz,rot,Dx,fault,dim,ix,
                                                   intersect,intersect_mixed);

                          /* d/dx, d/dy, d/dz */
                          if(gbox.contains(s))
                            {
                              SAMRAI::pdat::CellIndex c(s);
                              dv_diagonal(c,ix)+=intersect*jump(ix);
                            }

                          dv_mixed(s,0)+=intersect_mixed[0]*jump(ix);
                          dv_mixed(s,1)-=intersect_mixed[1]*jump(ix);
                        }
                      else
                        {
                          int intersect_diagonal, intersect_mixed[4],
                            intersect_corner[4];
                          compute_intersections_3D(ntt,xyz,rot,Dx,fault,dim,ix,
                                                   intersect_diagonal,
                                                   intersect_mixed,
                                                   intersect_corner);

                          /* d/dx, d/dy, d/dz */
                          if(gbox.contains(s))
                            {
                              SAMRAI::pdat::CellIndex c(s);
                              dv_diagonal(c,ix)+=intersect_diagonal*jump(ix);
                            }

                          for(int n=0;n<4;++n)
                            {
                              dv_mixed(s,n)+=intersect_mixed[n]*jump(ix);
                              dv_mixed(s,n+4)+=intersect_corner[n]*jump(ix);
                            }
                        }
                    }
                }
            }

          /* Corrections to the rhs */
          for(int ix=0;ix<dim;++ix)
            {
              SAMRAI::pdat::SideIterator s_end(pbox,ix,false);
              for(SAMRAI::pdat::SideIterator si(pbox,ix,true); si!=s_end;
                  si++)
                {
                  const SAMRAI::pdat::SideIndex s(*si);
                  SAMRAI::pdat::CellIndex c(s);

                  /* d/dx^2, d/dy^2, d/dz^2 */

                  double lambda_plus(cell_moduli(c,0)),
                    lambda_minus(cell_moduli(c-unit[ix],0)),
                    mu_plus(cell_moduli(c,1)),
                    mu_minus(cell_moduli(c-unit[ix],1));
                  v_rhs(s)+=
                    (dv_diagonal(c,ix)*(lambda_plus + 2*mu_plus)
                     - dv_diagonal(c-unit[ix],ix)*(lambda_minus + 2*mu_minus))
                    /(dx[ix]*dx[ix]);

                  for(int iy=(ix+1)%dim;iy!=ix;iy=(iy+1)%dim)
                    {
                      const int iz((ix+1)%dim!=iy ? (ix+1)%dim :
                                   (ix+2)%dim);
                      mu_plus=edge_node_eval(edge_moduli,s+unit[iy],iz,1);
                      mu_minus=edge_node_eval(edge_moduli,s,iz,1);

                      const int ix_iy(index_map(ix,iy,dim));
                      v_rhs(s)+=
                        (mu_plus*(dv_mixed(s,ix_iy)
                                  - dv_mixed(s+unit[iy],ix_iy+1))
                         + mu_minus*(dv_mixed(s,ix_iy+1)
                                     - dv_mixed(s-unit[iy],ix_iy)))
                        /(dx[iy]*dx[iy]);

                      /* d/dxy */
                      lambda_plus=cell_moduli(c,0);
                      lambda_minus=cell_moduli(c-unit[ix],0);
                      const SAMRAI::pdat::SideIndex
                        s_y(c,iy,SAMRAI::pdat::SideIndex::Lower);
                      const int iy_ix(index_map(iy,ix,dim));

                      v_rhs(s)+=
                        (lambda_plus*dv_diagonal(c,iy)
                         - lambda_minus*dv_diagonal(c-unit[ix],iy)
                         - mu_plus*(dv_mixed(s_y+unit[iy],iy_ix+1)
                                    - dv_mixed(s_y+unit[iy]-unit[ix],iy_ix))
                         + mu_minus*(dv_mixed(s_y,iy_ix+1)
                                     - dv_mixed(s_y-unit[ix],iy_ix)))
                        /(dx[ix]*dx[iy]);
                    }
                }
            }
        }
    }
}

#endif  // included_FACElastic
