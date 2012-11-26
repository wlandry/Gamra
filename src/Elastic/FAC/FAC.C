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

Elastic::FAC::FAC(const std::string& object_name,
                  const SAMRAI::tbox::Dimension& dimension,
                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database):
  d_object_name(object_name),
  d_dim(dimension),
  d_hierarchy(NULL),
  d_boundary_conditions(dimension,d_object_name + "::boundary conditions",
                        database->getDatabase("boundary_conditions")),
  d_elastic_fac_solver((d_dim),
                       object_name + "::fac_solver",
                       (!database.isNull() &&
                        database->isDatabase("fac_solver")) ?
                       database->getDatabase("fac_solver"):
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL),
                       d_boundary_conditions),
  d_context(),
  lambda("lambda",database,dimension),
  mu("mu",database,dimension),
  v_rhs("v_rhs",database,dimension,dimension.getValue())
{
  const int dim(d_dim.getValue());
  SAMRAI::hier::VariableDatabase* vdb =
    SAMRAI::hier::VariableDatabase::getDatabase();

  /*
   * Get a unique context for variables owned by this object.
   */
  d_context = vdb->getContext(d_object_name + ":Context");

  /*
   * Register variables with SAMRAI::hier::VariableDatabase
   * and get the descriptor indices for those variables.
   */

  int depth=2;
  SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<double> >
    cell_moduli_ptr(new SAMRAI::pdat::CellVariable<double>(d_dim,
                                                   object_name
                                                   + ":cell_moduli",depth));
  cell_moduli_id = vdb->registerVariableAndContext(cell_moduli_ptr,
                                                   d_context,
                                                   SAMRAI::hier::IntVector(d_dim, 1)
                                                   /* ghost cell width is
                                                      1 in case needed */);

  if(dim==2)
    {
      SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<double> >
        edge_moduli_ptr(new SAMRAI::pdat::NodeVariable<double>(d_dim,
                                                       object_name
                                                       + ":edge_moduli",depth));
      edge_moduli_id =
        vdb->registerVariableAndContext(edge_moduli_ptr,d_context,
                                        SAMRAI::hier::IntVector(d_dim,1)
                                        /* ghost cell width is 1 in
                                           case needed */);
    }
  else if(dim==3)
    {
      SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<double> >
        edge_moduli_ptr(new SAMRAI::pdat::EdgeVariable<double>(d_dim,
                                                       object_name
                                                       + ":edge_moduli",depth));
      edge_moduli_id =
        vdb->registerVariableAndContext(edge_moduli_ptr,d_context,
                                        SAMRAI::hier::IntVector(d_dim,1)
                                        /* ghost cell width is 1 in
                                           case needed */);
    }

  SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<double> >
    v_ptr(new SAMRAI::pdat::SideVariable<double>(d_dim, object_name + ":v", 1));
  v_id = vdb->registerVariableAndContext(v_ptr, d_context,
                                         SAMRAI::hier::IntVector(d_dim, 1)
                                         /* ghost cell width is 1 for
                                            stencil widths */);

  SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<double> >
    v_rhs_ptr(new SAMRAI::pdat::SideVariable<double>(d_dim,object_name
                                             + ":v right hand side"));
  v_rhs_id = vdb->registerVariableAndContext(v_rhs_ptr,d_context,
                                             SAMRAI::hier::IntVector(d_dim, 1)
                                             /* ghost cell width is
                                                1 for coarsening
                                                operator */);

  d_adaption_threshold=database->getDoubleWithDefault("adaption_threshold",
                                                      1.0e-15);
  min_full_refinement_level
    =database->getIntegerWithDefault("min_full_refinement_level",0);

  if(database->keyExists("faults"))
    {
      faults=database->getDoubleArray("faults");
    }
}
