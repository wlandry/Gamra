/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

  namespace Elastic {
    /* A little utility routine to validate the sizes of input arrays */
    void check_array_sizes(const tbox::Array<int> ijk,
                           const tbox::Array<double> min,
                           const tbox::Array<double> max,
                           const tbox::Array<double> data,
                           const int dim, const std::string &name,
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
      int data_size(1);
      for(int d=0; d<dim; ++d)
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
  Elastic::FAC::FAC(const std::string& object_name,
                    const tbox::Dimension& dimension,
                    tbox::Pointer<tbox::Database> database):
    d_object_name(object_name),
    d_dim(dimension),
    d_hierarchy(NULL),
    d_elastic_fac_solver((d_dim),
                         object_name + "::elastic_hypre",
                         (!database.isNull() &&
                          database->isDatabase("fac_solver")) ?
                         database->getDatabase("fac_solver"):
                         tbox::Pointer<tbox::Database>(NULL)),
    d_bc_coefs(d_dim,
               object_name + "::bc_coefs",
               (!database.isNull() &&
                database->isDatabase("bc_coefs")) ?
               database->getDatabase("bc_coefs"):
               tbox::Pointer<tbox::Database>(NULL)),
    d_context()
  {
    const int dim(d_dim.getValue());
    hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();

    /*
     * Get a unique context for variables owned by this object.
     */
    d_context = vdb->getContext(d_object_name + ":Context");

    /*
     * Register variables with hier::VariableDatabase
     * and get the descriptor indices for those variables.
     */

    tbox::Pointer<pdat::CellVariable<double> >
      p_ptr(new pdat::CellVariable<double>(d_dim, object_name + ":p", 1));
    p_id = vdb->registerVariableAndContext(p_ptr, d_context,
                                           hier::IntVector(d_dim, 1)
                                           /* ghost cell width is 1 for
                                              stencil widths */);

	int depth=2;
    tbox::Pointer<pdat::CellVariable<double> >
      cell_moduli_ptr(new pdat::CellVariable<double>(d_dim,
                                                        object_name
                                                        + ":cell_moduli",depth));
    cell_moduli_id = vdb->registerVariableAndContext(cell_moduli_ptr,
                                                        d_context,
                                                        hier::IntVector(d_dim, 1)
                                                        /* ghost cell width is
                                                           1 in case needed */);

    if(dim==2)
      {
        tbox::Pointer<pdat::NodeVariable<double> >
          edge_moduli_ptr(new pdat::NodeVariable<double>(d_dim,
                                                            object_name
                                                            + ":edge_moduli",depth));
        edge_moduli_id =
          vdb->registerVariableAndContext(edge_moduli_ptr,d_context,
                                          hier::IntVector(d_dim,1)
                                          /* ghost cell width is 1 in
                                             case needed */);
      }
    else if(dim==3)
      {
        tbox::Pointer<pdat::EdgeVariable<double> >
          edge_moduli_ptr(new pdat::EdgeVariable<double>(d_dim,
                                                            object_name
                                                            + ":edge_moduli",depth));
        edge_moduli_id =
          vdb->registerVariableAndContext(edge_moduli_ptr,d_context,
                                          hier::IntVector(d_dim,1)
                                          /* ghost cell width is 1 in
                                             case needed */);
      }

    tbox::Pointer<pdat::CellVariable<double> >
      dp_ptr(new pdat::CellVariable<double>(d_dim, object_name + ":dp"));
    dp_id = vdb->registerVariableAndContext(dp_ptr,d_context,
                                            hier::IntVector(d_dim, 1)
                                            /* ghost cell width is
                                                    1 in case needed */);

    tbox::Pointer<pdat::CellVariable<double> >
      p_exact_ptr(new pdat::CellVariable<double>(d_dim, object_name + ":p exact"));
    p_exact_id = vdb->registerVariableAndContext(p_exact_ptr,d_context,
                                                 hier::IntVector(d_dim, 1)
                                                 /* ghost cell width is
                                                    1 in case needed */);

    tbox::Pointer<pdat::CellVariable<double> >
      p_rhs_ptr(new pdat::CellVariable<double>(d_dim,object_name
                                               + ":p right hand side"));
    p_rhs_id = vdb->registerVariableAndContext(p_rhs_ptr,d_context,
                                               hier::IntVector(d_dim, 1));

    tbox::Pointer<pdat::SideVariable<double> >
      v_ptr(new pdat::SideVariable<double>(d_dim, object_name + ":v", 1));
    v_id = vdb->registerVariableAndContext(v_ptr, d_context,
                                           hier::IntVector(d_dim, 1)
                                           /* ghost cell width is 1 for
                                              stencil widths */);

    tbox::Pointer<pdat::SideVariable<double> >
      v_rhs_ptr(new pdat::SideVariable<double>(d_dim,object_name
                                               + ":v right hand side"));
    v_rhs_id = vdb->registerVariableAndContext(v_rhs_ptr,d_context,
                                               hier::IntVector(d_dim, 1)
                                               /* ghost cell width is
                                                  1 for coarsening
                                                  operator */);

    d_adaption_threshold=database->getDoubleWithDefault("adaption_threshold",
                                                        1.0e-15);
    min_full_refinement_level
      =database->getIntegerWithDefault("min_full_refinement_level",0);

    if(database->keyExists("lambda_data"))
      {
        lambda_ijk=database->getIntegerArray("lambda_ijk");
        lambda_xyz_min=database->getDoubleArray("lambda_coord_min");
        lambda_xyz_max=database->getDoubleArray("lambda_coord_max");
        lambda=database->getDoubleArray("lambda_data");
        check_array_sizes(lambda_ijk,lambda_xyz_min,lambda_xyz_max,
                          lambda,dim,"lambda");
      }

    if(database->keyExists("mu_data"))
      {
        mu_ijk=database->getIntegerArray("mu_ijk");
        mu_xyz_min=database->getDoubleArray("mu_coord_min");
        mu_xyz_max=database->getDoubleArray("mu_coord_max");
        mu=database->getDoubleArray("mu_data");
        check_array_sizes(mu_ijk,mu_xyz_min,mu_xyz_max,
                          mu,dim,"mu");
      }

    if(database->keyExists("v_rhs_data"))
      {
        v_rhs_ijk=database->getIntegerArray("v_rhs_ijk");
        v_rhs_xyz_min=database->getDoubleArray("v_rhs_coord_min");
        v_rhs_xyz_max=database->getDoubleArray("v_rhs_coord_max");
        v_rhs=database->getDoubleArray("v_rhs_data");
        check_array_sizes(v_rhs_ijk,v_rhs_xyz_min,v_rhs_xyz_max,
                          v_rhs,dim,"v_rhs",dim);
      }

    if(database->keyExists("p_initial_data"))
      {
        p_initial_ijk=database->getIntegerArray("p_initial_ijk");
        p_initial_xyz_min=database->getDoubleArray("p_initial_coord_min");
        p_initial_xyz_max=database->getDoubleArray("p_initial_coord_max");
        p_initial=database->getDoubleArray("p_initial_data");
        check_array_sizes(p_initial_ijk,p_initial_xyz_min,p_initial_xyz_max,
                          p_initial,dim,"p_initial");
      }

    /*
     * Specify an implementation of solv::RobinBcCoefStrategy for the
     * solver to use.  We use the implementation
     * solv::LocationIndexRobinBcCoefs, but other implementations are
     * possible, including user-implemented.
     */
    d_elastic_fac_solver.setBcObject(&d_bc_coefs);
  }
}
