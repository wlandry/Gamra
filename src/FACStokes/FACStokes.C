/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "FACStokes.h"

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

  /*
*************************************************************************
* Constructor creates a unique context for the object and register      *
* all its internal variables with the variable database.                *
*************************************************************************
*/
  FACStokes::FACStokes(const std::string& object_name,
                       const tbox::Dimension& dim,
                       tbox::Pointer<tbox::Database> database):
    d_object_name(object_name),
    d_dim(dim),
    d_hierarchy(NULL),
    d_stokes_fac_solver((d_dim),
                        object_name + "::stokes_hypre",
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
      p(new pdat::CellVariable<double>(dim, object_name + ":p", 1));
    p_id = vdb->registerVariableAndContext(p, d_context, hier::IntVector(dim, 1)
                                           /* ghost cell width is 1 for
                                              stencil widths */);

    tbox::Pointer<pdat::CellVariable<double> >
      cell_viscosity(new pdat::CellVariable<double>(dim,
                                                    object_name
                                                    + ":cell_viscosity"));
    cell_viscosity_id = vdb->registerVariableAndContext(cell_viscosity,
                                                        d_context,
                                                        hier::IntVector(dim, 1)
                                                        /* ghost cell width is
                                                           1 in case needed */);

    if(dim.getValue()==2)
      {
        tbox::Pointer<pdat::NodeVariable<double> >
          edge_viscosity(new pdat::NodeVariable<double>(dim,
                                                        object_name
                                                        + ":edge_viscosity"));
        edge_viscosity_id =
          vdb->registerVariableAndContext(edge_viscosity,d_context,
                                          hier::IntVector(dim,1)
                                          /* ghost cell width is 1 in
                                             case needed */);
      }
    else if(dim.getValue()==3)
      {
        tbox::Pointer<pdat::EdgeVariable<double> >
          edge_viscosity(new pdat::EdgeVariable<double>(dim,
                                                        object_name
                                                        + ":edge_viscosity"));
        edge_viscosity_id =
          vdb->registerVariableAndContext(edge_viscosity,d_context,
                                          hier::IntVector(dim,1)
                                          /* ghost cell width is 1 in
                                             case needed */);
      }

    tbox::Pointer<pdat::CellVariable<double> >
      dp(new pdat::CellVariable<double>(dim, object_name + ":dp"));
    dp_id = vdb->registerVariableAndContext(dp,d_context,
                                            hier::IntVector(dim, 1)
                                            /* ghost cell width is
                                                    1 in case needed */);

    tbox::Pointer<pdat::CellVariable<double> >
      p_exact(new pdat::CellVariable<double>(dim, object_name + ":p exact"));
    p_exact_id = vdb->registerVariableAndContext(p_exact,d_context,
                                                 hier::IntVector(dim, 1)
                                                 /* ghost cell width is
                                                    1 in case needed */);

    tbox::Pointer<pdat::CellVariable<double> >
      p_rhs(new pdat::CellVariable<double>(dim,object_name
                                           + ":p right hand side"));
    p_rhs_id = vdb->registerVariableAndContext(p_rhs,d_context,
                                               hier::IntVector(dim, 1));

    tbox::Pointer<pdat::SideVariable<double> >
      v(new pdat::SideVariable<double>(dim, object_name + ":v", 1));
    v_id = vdb->registerVariableAndContext(v, d_context, hier::IntVector(dim, 1)
                                           /* ghost cell width is 1 for
                                              stencil widths */);

    tbox::Pointer<pdat::SideVariable<double> >
      v_rhs(new pdat::SideVariable<double>(dim,object_name
                                           + ":v right hand side"));
    v_rhs_id = vdb->registerVariableAndContext(v_rhs,d_context,
                                               hier::IntVector(dim, 1)
                                               /* ghost cell width is
                                                  1 for coarsening
                                                  operator */);

    d_adaptation_threshold=database->getDoubleWithDefault("adaption_threshold",
                                                          // 1.0e-15);
                                                          2e-3);

    /*
     * Specify an implementation of solv::RobinBcCoefStrategy for the
     * solver to use.  We use the implementation
     * solv::LocationIndexRobinBcCoefs, but other implementations are
     * possible, including user-implemented.
     */
    d_stokes_fac_solver.setBcObject(&d_bc_coefs);
  }
}
