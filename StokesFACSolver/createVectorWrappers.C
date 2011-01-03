/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellVariable.h"
#include "StokesFACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

void SAMRAI::solv::StokesFACSolver::createVectorWrappers(int p, int p_rhs,
                                                         int v, int v_rhs) {

  hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());
  tbox::Pointer<hier::Variable> variable;

  if (!d_uv || d_uv->getComponentDescriptorIndex(0) != p) {
    d_uv.setNull();
    d_uv = new SAMRAIVectorReal<double>(d_object_name + "::uv",
                                        d_hierarchy,
                                        d_ln_min,
                                        d_ln_max);
    /* Add p */
    vdb.mapIndexToVariable(p, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << p << "\n");
    }
    tbox::Pointer<pdat::CellVariable<double> > cell_variable = variable;
    if (!cell_variable) {
      TBOX_ERROR(d_object_name << ": hier::Patch data index " << p
                 << " is not a cell-double variable.\n");
    }
#endif
    d_uv->addComponent(variable, p, s_weight_id[d_dim.getValue() - 1]);

    /* Add v */
    vdb.mapIndexToVariable(v, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << v << "\n");
    }
    tbox::Pointer<pdat::SideVariable<double> > side_variable = variable;
    if (!side_variable) {
      TBOX_ERROR(d_object_name << ": hier::Patch data index " << v
                 << " is not a side-double variable.\n");
    }
#endif
    d_uv->addComponent(variable, v, s_weight_id[d_dim.getValue() - 2]);
  }

  if (!d_fv || d_fv->getComponentDescriptorIndex(0) != p_rhs) {
    d_fv.setNull();
    d_fv = new SAMRAIVectorReal<double>(d_object_name + "::fv",
                                        d_hierarchy,
                                        d_ln_min,
                                        d_ln_max);
    /* Add p_rhs */
    vdb.mapIndexToVariable(p_rhs, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << p_rhs << "\n");
    }
    tbox::Pointer<pdat::CellVariable<double> > cell_variable = variable;
    if (!cell_variable) {
      TBOX_ERROR(d_object_name << ": hier::Patch data index " << p_rhs
                 << " is not a cell-double variable.\n");
    }
#endif
    d_fv->addComponent(variable, p_rhs, s_weight_id[d_dim.getValue() - 1]);

    /* Add v_rhs */
    vdb.mapIndexToVariable(v_rhs, variable);    
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << v_rhs << "\n");
    }
    tbox::Pointer<pdat::SideVariable<double> > side_variable = variable;
    if (!side_variable) {
      TBOX_ERROR(d_object_name << ": hier::Patch data index " << v_rhs
                 << " is not a cell-double variable.\n");
    }
#endif
    d_fv->addComponent(variable, v_rhs, s_weight_id[d_dim.getValue() - 2]);
  }
}
