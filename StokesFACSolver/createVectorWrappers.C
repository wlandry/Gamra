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

namespace SAMRAI {
  namespace solv {

    void StokesFACSolver::createVectorWrappers(
                                               int u,
                                               int f) {

      hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());
      tbox::Pointer<hier::Variable> variable;

      if (!d_uv || d_uv->getComponentDescriptorIndex(0) != u) {
        d_uv.setNull();
        d_uv = new SAMRAIVectorReal<double>(d_object_name + "::uv",
                                            d_hierarchy,
                                            d_ln_min,
                                            d_ln_max);
        vdb.mapIndexToVariable(u, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
        if (!variable) {
          TBOX_ERROR(d_object_name << ": No variable for patch data index "
                     << u << "\n");
        }
        tbox::Pointer<pdat::CellVariable<double> > cell_variable = variable;
        if (!cell_variable) {
          TBOX_ERROR(d_object_name << ": hier::Patch data index " << u
                     << " is not a cell-double variable.\n");
        }
#endif
        d_uv->addComponent(variable, u, s_weight_id[d_dim.getValue() - 1]);
      }

      if (!d_fv || d_fv->getComponentDescriptorIndex(0) != f) {
        d_fv.setNull();
        d_fv = new SAMRAIVectorReal<double>(d_object_name + "::fv",
                                            d_hierarchy,
                                            d_ln_min,
                                            d_ln_max);
        vdb.mapIndexToVariable(f, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
        if (!variable) {
          TBOX_ERROR(d_object_name << ": No variable for patch data index "
                     << f << "\n");
        }
        tbox::Pointer<pdat::CellVariable<double> > cell_variable = variable;
        if (!cell_variable) {
          TBOX_ERROR(d_object_name << ": hier::Patch data index " << f
                     << " is not a cell-double variable.\n");
        }
#endif
        d_fv->addComponent(variable, f, s_weight_id[d_dim.getValue() - 1]);
      }
    }
  }
}
