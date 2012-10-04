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
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

/*
*************************************************************************
* Set up the initial guess and problem parameters                       *
* and solve the Elastic problem.  We explicitly initialize and          *
* deallocate the solver state in this example.                          *
*************************************************************************
*/
int Elastic::FAC::solve()
{

  if (d_hierarchy.isNull()) {
    TBOX_ERROR(d_object_name
               << "Cannot solve using an uninitialized object.\n");
  }

  /*
   * Fill in the initial guess.
   */
  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel>
      level = d_hierarchy->getPatchLevel(ln);
    SAMRAI::hier::PatchLevel::Iterator ip(*level);
    for ( ; ip; ip++) {
      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch> patch = *ip;

      SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
        geom = patch->getPatchGeometry();

      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
        v = patch->getPatchData(v_id);
      v->fill(0.0);
      // v->fill(13.0);
    }
    d_elastic_fac_solver.set_boundaries(v_id,level,false);
  }

  fix_moduli();

  d_elastic_fac_solver.initializeSolverState
    (cell_moduli_id,edge_moduli_id,v_id,v_rhs_id,
     d_hierarchy,0,d_hierarchy->getFinestLevelNumber());

  SAMRAI::tbox::plog << "solving..." << std::endl;
  int solver_ret;
  solver_ret = d_elastic_fac_solver.solveSystem(v_id,v_rhs_id);
  /*
   * Present data on the solve.
   */
  double avg_factor, final_factor;
  d_elastic_fac_solver.getConvergenceFactors(avg_factor, final_factor);
  SAMRAI::tbox::plog << "\t" << (solver_ret ? "" : "NOT ") << "converged " << "\n"
             << "	iterations: "
             << d_elastic_fac_solver.getNumberOfIterations() << "\n"
             << "	residual: "<< d_elastic_fac_solver.getResidualNorm()
             << "\n"
             << "	average convergence: "<< avg_factor << "\n"
             << "	final convergence: "<< final_factor << "\n"
             << std::flush;

  d_elastic_fac_solver.deallocateSolverState();

  return 0;
}
