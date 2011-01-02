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
#include "StokesSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

  /*
*************************************************************************
* Set up the initial guess and problem parameters                       *
* and solve the Stokes problem.  We explicitly initialize and          *
* deallocate the solver state in this example.                          *
*************************************************************************
*/
  int FACStokes::solveStokes()
  {

    if (d_hierarchy.isNull()) {
      TBOX_ERROR(d_object_name
                 << "Cannot solve using an uninitialized object.\n");
    }

    int ln;
    /*
     * Fill in the initial guess.
     */
    for (ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::Iterator ip(*level);
      for ( ; ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;
        tbox::Pointer<pdat::CellData<double> > data = patch->getPatchData(
                                                                          p_id);
        data->fill(0.0);
        tbox::Pointer<pdat::SideData<double> > vdata = patch->getPatchData(
                                                                           v_id);
        vdata->fill(0.0);
      }
    }

    /*
     * Set the parameters for the Stokes equation.
     * See classes solv::StokesFACSolver or
     * solv::StokesSpecifications.
     * (D is the diffusion coefficient.
     * C is the source term which is not used in this example.)
     */
    d_stokes_fac_solver.setDConstant(1.0);
    d_stokes_fac_solver.setCConstant(0.0);

    d_stokes_fac_solver.initializeSolverState
      (p_id,p_rhs_id,d_hierarchy,0,
       d_hierarchy->getFinestLevelNumber());

    tbox::plog << "solving..." << std::endl;
    int solver_ret;
    solver_ret = d_stokes_fac_solver.solveSystem(p_id,p_rhs_id,v_id,v_rhs_id);
    /*
     * Present data on the solve.
     */
    double avg_factor, final_factor;
    d_stokes_fac_solver.getConvergenceFactors(avg_factor, final_factor);
    tbox::plog << "\t" << (solver_ret ? "" : "NOT ") << "converged " << "\n"
               << "	iterations: "
               << d_stokes_fac_solver.getNumberOfIterations() << "\n"
               << "	residual: "<< d_stokes_fac_solver.getResidualNorm()
               << "\n"
               << "	average convergence: "<< avg_factor << "\n"
               << "	final convergence: "<< final_factor << "\n"
               << std::flush;

    d_stokes_fac_solver.deallocateSolverState();

    return 0;
  }

}
