/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.hxx"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/hier/VariableDatabase.h>

/*
*************************************************************************
* Set up the initial guess and problem parameters                       *
* and solve the Elastic problem.  We explicitly initialize and          *
* deallocate the solver state in this example.                          *
*************************************************************************
*/
bool Elastic::FAC::solve()
{

  if (!d_hierarchy)
    { TBOX_ERROR("Elastic::FAC: Cannot solve using an uninitialized "
                 "object.\n"); }

  d_boundary_conditions.set_extra_ids(edge_moduli_id,dv_diagonal_id,
                                      dv_mixed_id,level_set_id);

  fix_moduli();
  const int dim(d_dim.getValue());
  if(!faults.empty())
    {
      if(dim==2)
        add_faults<SAMRAI::pdat::NodeData<double> >();
      else
        add_faults<SAMRAI::pdat::EdgeData<double> >();
    }

  /* Fill in the initial guess. */
  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
      boost::shared_ptr<SAMRAI::hier::PatchLevel>
        level = d_hierarchy->getPatchLevel(ln);
    
      for (SAMRAI::hier::PatchLevel::Iterator ip(level->begin());
           ip!=level->end(); ++ip)
        {
          const boost::shared_ptr<SAMRAI::hier::Patch> patch(*ip);

          const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v_ptr
            (boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
             (patch->getPatchData(v_id)));
          SAMRAI::pdat::SideData<double> &v(*v_ptr);
          v.fill(0.0);

          const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry>
            &geom(boost::dynamic_pointer_cast
                  <SAMRAI::geom::CartesianPatchGeometry>
                  (patch->getPatchGeometry()));
          const SAMRAI::hier::Box &pbox(v.getBox());
          const SAMRAI::hier::Box &gbox(v.getGhostBox());
          const double *dx=geom->getDx();
          for(Gamra::Dir ix=0;ix<dim;++ix)
            {
              if(v_initial[ix].is_valid)
                {
                  double offset[]={0.5,0.5,0.5};
                  offset[ix]=0;
                  
                  SAMRAI::pdat::SideIterator
                    s_end(SAMRAI::pdat::SideGeometry::end(gbox,ix));
                  for(SAMRAI::pdat::SideIterator
                        si(SAMRAI::pdat::SideGeometry::begin(gbox,ix));
                      si!=s_end; ++si)
                    {
                      const SAMRAI::pdat::SideIndex &s(*si);

                      double coord[3];
                      for(Gamra::Dir d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(s[d]-pbox.lower()[d]+offset[d]);
                      v(s)=v_initial[ix].eval(coord);
                    }
                }
            }
        }
    d_elastic_fac_solver.set_physical_boundaries(v_id,level,false);
  }

  d_elastic_fac_solver.initializeSolverState
    (cell_moduli_id,edge_moduli_id,dv_diagonal_id,dv_mixed_id,level_set_id,
     v_id,v_rhs_id,d_hierarchy,0,d_hierarchy->getFinestLevelNumber());

  SAMRAI::tbox::plog << "solving..." << std::endl;
  bool solver_ret(d_elastic_fac_solver.solveSystem(v_id,v_rhs_id));

  /// Write out convergence data
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

  return solver_ret;
}
