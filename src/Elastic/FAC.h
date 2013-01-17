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
                                    int depth_id) const
    {
      if(variable_name=="Diagonal_Strain")
        return pack_diagonal_strain(buffer,patch,region,depth_id);
      else if(variable_name=="Tangent_Strain")
        return pack_tangent_strain(buffer,patch,region,depth_id);
      else
        return pack_v_v_rhs(buffer,patch,region,variable_name,depth_id);
    }

    bool
    pack_tangent_strain(double* buffer,
                        const SAMRAI::hier::Patch& patch,
                        const SAMRAI::hier::Box& region,
                        const int &depth_id) const;

    bool
    pack_diagonal_strain(double* buffer,
                         const SAMRAI::hier::Patch& patch,
                         const SAMRAI::hier::Box& region,
                         const int &depth_id) const;

    bool
    pack_v_v_rhs(double* buffer,
                 const SAMRAI::hier::Patch& patch,
                 const SAMRAI::hier::Box& region,
                 const std::string& variable_name,
                 const int &depth_id) const;
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
      dv_perpendicular_id;

    Input_Expression lambda, mu, v_rhs;

    SAMRAI::tbox::Array<double> faults;
    //@}

    static const int index_map[3][3];

    template<class T> void add_faults();

    template<class T>
    void compute_dv_correction
    (const double fault[],
     const FTensor::Tensor1<double,3> &jump,
     const SAMRAI::pdat::SideIndex &s,
     const SAMRAI::hier::Index unit[],
     const int &ix,
     const int &d,
     const int &dim,
     const FTensor::Tensor1<double,3> &ntt,
     const FTensor::Tensor1<double,3> &ntt_p,
     SAMRAI::pdat::CellData<double> &dv_diagonal,
     T &dv_perpendicular);

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
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_data =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            ((*p)->getPatchData(v_rhs_id));

          /* dv */
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            ((*p)->getPatchData(dv_diagonal_id));
          boost::shared_ptr<T> dv_perpendicular =
            boost::dynamic_pointer_cast<T>
            ((*p)->getPatchData(dv_perpendicular_id));
          dv_diagonal->fillAll(0);
          dv_perpendicular->fillAll(0);

          /* moduli */
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            ((*p)->getPatchData(cell_moduli_id));
          boost::shared_ptr<T> edge_moduli =
            boost::dynamic_pointer_cast<T>((*p)->getPatchData(edge_moduli_id));


          /* Iterate over the faults */
          const int params_per_fault(9);
          for(int fault_index=0;fault_index<faults.size();
              fault_index+=params_per_fault)
            {
              const double pi=4*atan(1);
              /* The conventions for faults are different from the regular
               * xyz coordinates, so we have to convert.  Depth is a
               * coordinate, giving a left-handed coordinate system.  So we
               * invert the depth and width to convert to a right-handed
               * coordinate system.  Strike is opposite in direction, and
               * dip is measured from a plane lying flat. */
              double scale(faults[fault_index+0]), x(faults[fault_index+2]),
                y(faults[fault_index+1]), z(-faults[fault_index+3]),
                L(faults[fault_index+4]), W(-faults[fault_index+5]),
                strike(faults[fault_index+6]*pi/180),
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

              SAMRAI::hier::Box pbox = v_rhs_data->getBox();
              for(int ix=0;ix<dim;++ix)
                {
                  double offset[]={0.5,0.5,0.5};
                  offset[ix]=0;

                  SAMRAI::pdat::SideIterator send(pbox,ix,false);
                  for(SAMRAI::pdat::SideIterator si(pbox,ix,true); si!=send;
                      si++)
                    {
                      const SAMRAI::pdat::SideIndex s(*si);

                      FTensor::Tensor1<double,3> xyz(0,0,0);
                      for(int d=0;d<dim;++d)
                        xyz(d)=geom->getXLower()[d]
                          + dx[d]*(s[d]-pbox.lower()[d]+offset[d]);

                      /* Rotate the coordinates into the coordinates of the
                         fault.  So in those coordinates, if x<0, you are on
                         the left, and if x>0, you are on the right. */
                      FTensor::Tensor1<double,3> ntt;
                      FTensor::Tensor1<double,3> ntt_dp[dim], ntt_dm[dim];
                      ntt(a)=rot(a,b)*(xyz(b)-center(b));
                      for(int d=0;d<dim;++d)
                        {
                          ntt_dp[d](a)=rot(a,b)*(xyz(b)+Dx[d](b)-center(b));
                          ntt_dm[d](a)=rot(a,b)*(xyz(b)-Dx[d](b)-center(b));
                        }

                      /* d/dx, d/dy, d/dz */
                      for(int d=0;d<dim;++d)
                        {
                          compute_dv_correction(fault,jump,s,unit,ix,d,dim,
                                                ntt,ntt_dp[d],
                                                *dv_diagonal,*dv_perpendicular);
                          if(s[d]==pbox.lower()[d])
                            compute_dv_correction(fault,jump,s-unit[d],unit,ix,d,
                                                  dim,ntt_dm[d],ntt,
                                                  *dv_diagonal,
                                                  *dv_perpendicular);

                        }

                      /* d/dx^2, d/dy^2, d/dz^2 */
                      for(int d=0;d<dim;++d)
                        {
                          int sign(0);
                          SAMRAI::hier::Index cell_offset(zero),
                            edge_offset(zero);
                          if(((ntt(0)<=0 && ntt_dp[d](0)>0)
                              || (ntt(0)>0 && ntt_dp[d](0)<=0))
                             && intersect_fault(dim,ntt,ntt_dp[d],fault))
                            {
                              sign=ntt(0)<=0 ? 1 : -1;
                              edge_offset=unit[d];
                            }
                          else if(((ntt(0)<=0 && ntt_dm[d](0)>0)
                                   || (ntt(0)>0 && ntt_dm[d](0)<=0))
                                  && intersect_fault(dim,ntt,ntt_dm[d],
                                                     fault))
                            {
                              sign=ntt(0)<=0 ? 1 : -1;
                              cell_offset=-unit[d];
                            }

                          if(sign!=0)
                            {
                              /* This is the lambda and mu that multiply the
                                 jump term.  So we need the lambda and mu
                                 that are on the side where the jump
                                 occurs. */
                              double lambda_here, mu_here;
                              if(ix==d)
                                {
                                  SAMRAI::pdat::CellIndex c(s+cell_offset);
                                  lambda_here=(*cell_moduli)(c,0);
                                  mu_here=(*cell_moduli)(c,1);
                                }
                              else
                                {
                                  const int iz((ix+1)%dim!=d ? (ix+1)%dim :
                                               (ix+2)%dim);
                                  lambda_here=edge_node_eval(*edge_moduli,
                                                             s+edge_offset,
                                                             iz,0);
                                  mu_here=edge_node_eval(*edge_moduli,
                                                         s+edge_offset,
                                                         iz,1);
                                }
                              double factor(mu_here);
                              if(ix==d)
                                factor=lambda_here+2*mu_here;
                              (*v_rhs_data)(s)+=sign*factor*jump(ix)
                                /(dx[d]*dx[d]);
                            }
                        }

                      /* d/dxy */
                      for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                        {
                          const int iz((iy+1)%dim==ix ? (iy+2)%dim : (iy+1)%dim);
                          FTensor::Tensor1<double,3> ntt_dxy[2][2];
                          ntt_dxy[0][0](a)=rot(a,b)*(xyz(b)+Dx[ix](b)/2
                                                     +Dx[iy](b)/2-center(b));
                          ntt_dxy[0][1](a)=rot(a,b)*(xyz(b)+Dx[ix](b)/2
                                                     -Dx[iy](b)/2-center(b));
                          ntt_dxy[1][0](a)=rot(a,b)*(xyz(b)-Dx[ix](b)/2
                                                     +Dx[iy](b)/2-center(b));
                          ntt_dxy[1][1](a)=rot(a,b)*(xyz(b)-Dx[ix](b)/2
                                                     -Dx[iy](b)/2-center(b));

                          double lambda_mu(0);

                          /* mu terms */
                          if(((ntt_dxy[0][0](0)<=0 && ntt_dxy[1][0](0)>0)
                              || (ntt_dxy[0][0](0)>0 && ntt_dxy[1][0](0)<=0))
                             && intersect_fault(dim,ntt_dxy[0][0],
                                                ntt_dxy[1][0],fault))
                            {
                              lambda_mu+=(ntt_dxy[0][0](0)<=0 ? -1 : 1)*
                                edge_node_eval(*edge_moduli,s+unit[iy],iz,1);
                            }
                          if(((ntt_dxy[0][1](0)<=0 && ntt_dxy[1][1](0)>0)
                              || (ntt_dxy[0][1](0)>0 && ntt_dxy[1][1](0)<=0))
                             && intersect_fault(dim,ntt_dxy[0][1],
                                                ntt_dxy[1][1],fault))
                            {
                              lambda_mu+=(ntt_dxy[0][1](0)<=0 ? 1 : -1)
                                *edge_node_eval(*edge_moduli,s,iz,1);
                            }

                          /* lambda terms */
                          if(((ntt_dxy[0][0](0)<=0 && ntt_dxy[0][1](0)>0)
                              || (ntt_dxy[0][0](0)>0 && ntt_dxy[0][1](0)<=0))
                             && intersect_fault(dim,ntt_dxy[0][0],
                                                ntt_dxy[0][1],fault))
                            {
                              SAMRAI::pdat::CellIndex c(s);
                              lambda_mu+=(ntt_dxy[0][0](0)<=0 ? -1 : 1)
                                *(*cell_moduli)(c,0);
                            }
                          if(((ntt_dxy[1][0](0)<=0 && ntt_dxy[1][1](0)>0)
                              || (ntt_dxy[1][0](0)>0 && ntt_dxy[1][1](0)<=0))
                             && intersect_fault(dim,ntt_dxy[1][0],
                                                ntt_dxy[1][1],fault))
                            {
                              SAMRAI::pdat::CellIndex c(s-unit[ix]);
                              lambda_mu+=(ntt_dxy[1][0](0)<=0 ? 1 : -1)
                                *(*cell_moduli)(c,0);
                            }
                          (*v_rhs_data)(s)+=lambda_mu*jump(iy)
                            /(dx[0]*dx[1]);
                        }
                    }
                }
            }
        }
    }
}

template<class T>
void Elastic::FAC::compute_dv_correction
(const double fault[],
 const FTensor::Tensor1<double,3> &jump,
 const SAMRAI::pdat::SideIndex &s,
 const SAMRAI::hier::Index unit[],
 const int &ix,
 const int &d,
 const int &dim,
 const FTensor::Tensor1<double,3> &ntt,
 const FTensor::Tensor1<double,3> &ntt_p,
 SAMRAI::pdat::CellData<double> &dv_diagonal,
 T &dv_perpendicular)
{
  int sign(0);
  if(ntt(0)<=0 && ntt_p(0)>0 && intersect_fault(dim,ntt,ntt_p,fault))
    sign=1;
  else if(ntt(0)>0 && ntt_p(0)<=0 && intersect_fault(dim,ntt,ntt_p,fault))
    sign=-1;

  if(sign!=0)
    {
      if(ix==d)
        {
          SAMRAI::pdat::CellIndex c(s);
          dv_diagonal(c,ix)+=sign*jump(ix);
        }
      else
        {
          edge_node_eval(dv_perpendicular,s+unit[d],
                         (ix+1)%dim!=d ? (ix+1)%dim : (ix+2)%dim,
                         index_map[ix][d])
            +=sign*jump(ix);
        }
    }
}


#endif  // included_FACElastic
