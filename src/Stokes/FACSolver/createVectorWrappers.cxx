/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#include <SAMRAI/pdat/CellVariable.h>
#include "Stokes/FACSolver.hxx"
#include <SAMRAI/tbox/PIO.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/tbox/StartupShutdownManager.h>

#include IOMANIP_HEADER_FILE

void Stokes::FACSolver::createVectorWrappers(int p, int p_rhs,
                                                           int v, int v_rhs) {

  SAMRAI::hier::VariableDatabase& vdb(*SAMRAI::hier::VariableDatabase::getDatabase());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;

  if (!d_uv || d_uv->getComponentDescriptorIndex(0) != p) {
    d_uv.reset();
    d_uv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
      (d_object_name + "::uv", d_hierarchy, d_ln_min, d_ln_max);
    /* Add p */
    vdb.mapIndexToVariable(p, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << p << "\n");
    }
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > cell_variable =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
      (variable);
    if (!cell_variable) {
      TBOX_ERROR(d_object_name << ": SAMRAI::hier::Patch data index " << p
                 << " is not a cell-double variable.\n");
    }
#endif

    /* Add v */
    vdb.mapIndexToVariable(v, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << v << "\n");
    }
    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > side_variable =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
      (variable);
    if (!side_variable) {
      TBOX_ERROR(d_object_name << ": SAMRAI::hier::Patch data index " << v
                 << " is not a side-double variable.\n");
    }
#endif
    d_uv->addComponent(variable, v);
  }

  if (!d_fv || d_fv->getComponentDescriptorIndex(0) != p_rhs) {
    d_fv.reset();
    d_fv = boost::make_shared<SAMRAI::solv::SAMRAIVectorReal<double> >
      (d_object_name + "::fv", d_hierarchy, d_ln_min, d_ln_max);
    /* Add p_rhs */
    vdb.mapIndexToVariable(p_rhs, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << p_rhs << "\n");
    }
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > cell_variable =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
      (variable);
    if (!cell_variable) {
      TBOX_ERROR(d_object_name << ": SAMRAI::hier::Patch data index " << p_rhs
                 << " is not a cell-double variable.\n");
    }
#endif

    /* Add v_rhs */
    vdb.mapIndexToVariable(v_rhs, variable);    
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!variable) {
      TBOX_ERROR(d_object_name << ": No variable for patch data index "
                 << v_rhs << "\n");
    }
    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > side_variable =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
      (variable);
    if (!side_variable) {
      TBOX_ERROR(d_object_name << ": SAMRAI::hier::Patch data index " << v_rhs
                 << " is not a cell-double variable.\n");
    }
#endif
    d_fv->addComponent(variable, v_rhs);
  }
}
