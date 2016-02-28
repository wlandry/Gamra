/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "Stokes/FAC.hxx"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/solv/SimpleCellRobinBcCoefs.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/math/HierarchyCellDataOpsReal.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/hier/VariableDatabase.h>

namespace Stokes {
  /* A little utility routine to validate the sizes of input arrays */
  void check_array_sizes(const std::vector<int> ijk,
                         const std::vector<double> min,
                         const std::vector<double> max,
                         const std::vector<double> data,
                         const unsigned dim, const std::string &name,
                         const int num_components=1)
  {
    if(ijk.size()!=dim)
      TBOX_ERROR("Bad number of elements in " << name << "_ijk.  Expected "
                 << dim << " but got " << ijk.size());
    if(min.size()!=dim)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_coord_min.  Expected "
                 << dim << " but got " << min.size());
    if(max.size()!=dim)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_coord_max.  Expected "
                 << dim << " but got " << max.size());
    size_t data_size(1);
    for(unsigned d=0; d<dim; ++d)
      data_size*=ijk[d];
    if(data.size()!=data_size*num_components)
      TBOX_ERROR("Bad number of elements in "
                 << name << "_data.  Expected "
                 << data_size << " but got " << data.size());
  }
}
/*
*************************************************************************
* Constructor creates a unique context for the object and register      *
* all its internal variables with the variable database.                *
*************************************************************************
*/
Stokes::FAC::FAC(const SAMRAI::tbox::Dimension& dimension,
                 boost::shared_ptr<SAMRAI::tbox::Database> database):
  d_dim(dimension),
  d_stokes_fac_solver((d_dim),
                      "Stokes::FAC::stokes_hypre",
                      (database &&
                       database->isDatabase("fac_solver")) ?
                      database->getDatabase("fac_solver"):
                      boost::shared_ptr<SAMRAI::tbox::Database>()),
  d_bc_coefs(d_dim,
             "Stokes::FAC::bc_coefs",
             (database &&
              database->isDatabase("bc_coefs")) ?
             database->getDatabase("bc_coefs"):
             boost::shared_ptr<SAMRAI::tbox::Database>()),
  d_context()
{
  const unsigned dim(d_dim.getValue());
  SAMRAI::hier::VariableDatabase* vdb =
    SAMRAI::hier::VariableDatabase::getDatabase();

  /*
   * Get a unique context for variables owned by this object.
   */
  d_context = vdb->getContext("Stokes::FAC:Context");

  /*
   * Register variables with SAMRAI::hier::VariableDatabase
   * and get the descriptor indices for those variables.
   */

  boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    p_ptr(new SAMRAI::pdat::CellVariable<double>(d_dim, "Stokes::FAC:p", 1));
  p_id = vdb->registerVariableAndContext(p_ptr, d_context,
                                         SAMRAI::hier::IntVector::getOne(d_dim)
                                         /* ghost cell width is 1 for
                                            stencil widths */);

  boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    cell_viscosity_ptr(new SAMRAI::pdat::CellVariable<double>
                       (d_dim, "Stokes::FAC:cell_viscosity"));
  cell_viscosity_id = vdb->registerVariableAndContext
    (cell_viscosity_ptr, d_context, SAMRAI::hier::IntVector::getOne(d_dim)
                                                      /* ghost cell width is
                                                         1 in case needed */);

  if(dim==2)
    {
      boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> >
        edge_viscosity_ptr(new SAMRAI::pdat::NodeVariable<double>
                           (d_dim, "Stokes::FAC:edge_viscosity"));
      edge_viscosity_id =
        vdb->registerVariableAndContext(edge_viscosity_ptr,d_context,
                                        SAMRAI::hier::IntVector(d_dim,1)
                                        /* ghost cell width is 1 in
                                           case needed */);
    }
  else if(dim==3)
    {
      boost::shared_ptr<SAMRAI::pdat::EdgeVariable<double> >
        edge_viscosity_ptr(new SAMRAI::pdat::EdgeVariable<double>
                           (d_dim, "Stokes::FAC:edge_viscosity"));
      edge_viscosity_id =
        vdb->registerVariableAndContext(edge_viscosity_ptr,d_context,
                                        SAMRAI::hier::IntVector(d_dim,1)
                                        /* ghost cell width is 1 in
                                           case needed */);
    }

  boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    dp_ptr(new SAMRAI::pdat::CellVariable<double>(d_dim, "Stokes::FAC:dp"));
  dp_id = vdb->registerVariableAndContext(dp_ptr,d_context,
                                          SAMRAI::hier::IntVector::getOne(d_dim)
                                          /* ghost cell width is
                                             1 in case needed */);

  boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    p_exact_ptr(new SAMRAI::pdat::CellVariable<double>
                (d_dim, "Stokes::FAC:p exact"));
  p_exact_id = vdb->registerVariableAndContext(p_exact_ptr,d_context,
                                               SAMRAI::hier::IntVector::getOne(d_dim)
                                               /* ghost cell width is
                                                  1 in case needed */);

  boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    p_rhs_ptr(new SAMRAI::pdat::CellVariable<double>
              (d_dim,"Stokes::FAC:p right hand side"));
  p_rhs_id = vdb->registerVariableAndContext(p_rhs_ptr,d_context,
                                             SAMRAI::hier::IntVector::getOne(d_dim));

  boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    v_ptr(new SAMRAI::pdat::SideVariable<double>
          (d_dim, "Stokes::FAC:v",SAMRAI::hier::IntVector::getOne(d_dim),1));
  v_id = vdb->registerVariableAndContext(v_ptr, d_context,
                                         SAMRAI::hier::IntVector::getOne(d_dim)
                                         /* ghost cell width is 1 for
                                            stencil widths */);

  boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    v_rhs_ptr(new SAMRAI::pdat::SideVariable<double>
              (d_dim,"Stokes::FAC:v right hand side",
               SAMRAI::hier::IntVector::getOne(d_dim)));
  v_rhs_id = vdb->registerVariableAndContext(v_rhs_ptr,d_context,
                                             SAMRAI::hier::IntVector::getOne(d_dim)
                                             /* ghost cell width is
                                                1 for coarsening
                                                operator */);

  d_adaption_threshold=database->getDoubleWithDefault("adaption_threshold",
                                                      1.0e-15);
  min_full_refinement_level
    =database->getIntegerWithDefault("min_full_refinement_level",0);

  if(database->keyExists("viscosity_data"))
    {
      viscosity_ijk=database->getIntegerVector("viscosity_ijk");
      viscosity_xyz_min=database->getDoubleVector("viscosity_coord_min");
      viscosity_xyz_max=database->getDoubleVector("viscosity_coord_max");
      viscosity=database->getDoubleVector("viscosity_data");
      check_array_sizes(viscosity_ijk,viscosity_xyz_min,viscosity_xyz_max,
                        viscosity,dim,"viscosity");
    }

  if(database->keyExists("v_rhs_data"))
    {
      v_rhs_ijk=database->getIntegerVector("v_rhs_ijk");
      v_rhs_xyz_min=database->getDoubleVector("v_rhs_coord_min");
      v_rhs_xyz_max=database->getDoubleVector("v_rhs_coord_max");
      v_rhs=database->getDoubleVector("v_rhs_data");
      check_array_sizes(v_rhs_ijk,v_rhs_xyz_min,v_rhs_xyz_max,
                        v_rhs,dim,"v_rhs",dim);
    }

  if(database->keyExists("p_initial_data"))
    {
      p_initial_ijk=database->getIntegerVector("p_initial_ijk");
      p_initial_xyz_min=database->getDoubleVector("p_initial_coord_min");
      p_initial_xyz_max=database->getDoubleVector("p_initial_coord_max");
      p_initial=database->getDoubleVector("p_initial_data");
      check_array_sizes(p_initial_ijk,p_initial_xyz_min,p_initial_xyz_max,
                        p_initial,dim,"p_initial");
    }

  /*
   * Specify an implementation of solv::RobinBcCoefStrategy for the
   * solver to use.  We use the implementation
   * solv::LocationIndexRobinBcCoefs, but other implementations are
   * possible, including user-implemented.
   */
  d_stokes_fac_solver.setBcObject(&d_bc_coefs);
}
