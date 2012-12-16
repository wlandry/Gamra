/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "Stokes/FAC.h"

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

/*
*************************************************************************
* Set up the initial guess and problem parameters                       *
* and solve the Stokes problem.  We explicitly initialize and          *
* deallocate the solver state in this example.                          *
*************************************************************************
*/
int SAMRAI::Stokes::FAC::solve()
{

  if (!d_hierarchy) {
    TBOX_ERROR(d_object_name
               << "Cannot solve using an uninitialized object.\n");
  }

  int ln;
  /*
   * Fill in the initial guess.
   */
  for (ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
    boost::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
    hier::PatchLevel::Iterator ip(level->begin());
    hier::PatchLevel::Iterator iend(level->end());
    for ( ; ip!=iend; ip++) {
      boost::shared_ptr<hier::Patch> patch = *ip;
      boost::shared_ptr<pdat::CellData<double> > p =
        boost::dynamic_pointer_cast<pdat::CellData<double> >
        (patch->getPatchData(p_id));

      boost::shared_ptr<geom::CartesianPatchGeometry> geom =
        boost::dynamic_pointer_cast<geom::CartesianPatchGeometry>
        (patch->getPatchGeometry());

      if(p_initial.empty())
        {
          p->fill(0.0);
        }
      else
        {
          const int dim=d_dim.getValue();
          const double *dx=geom->getDx();
          double dx_p[dim];
          for(int d=0;d<dim;++d)
            dx_p[d]=(p_initial_xyz_max[d]
                     - p_initial_xyz_min[d])/(p_initial_ijk[d]-1);
          int di[dim];
          di[0]=1;
          for(int d=1;d<dim;++d)
            di[d]=di[d-1]*p_initial_ijk[d-1];

          hier::Box pbox = p->getBox();
          pdat::CellIterator cend(p->getGhostBox(),false);
          for(pdat::CellIterator ci(p->getGhostBox(),true); ci!=cend; ci++)
            {
              pdat::CellIndex c=*ci;
              double xyz[dim], weight[dim][2];
              for(int d=0;d<dim;++d)
                xyz[d]=geom->getXLower()[d]
                  + dx[d]*(c[d]-pbox.lower()[d] + 0.5);

              int ijk(0);
              int ddi[dim];
              for(int d=0;d<dim;++d)
                {
                  int i=static_cast<int>(xyz[d]*(p_initial_ijk[d]-1)
                                         /(p_initial_xyz_max[d]
                                           - p_initial_xyz_min[d]));
                  i=std::max(0,std::min(p_initial_ijk[d]-1,i));
                  ijk+=i*di[d];

                  if(i==p_initial_ijk[d]-1)
                    {
                      weight[d][0]=1;
                      weight[d][1]=0;
                      ddi[d]=0;
                    }
                  else
                    {
                      weight[d][1]=
                        (xyz[d]-(i*dx_p[d] + p_initial_xyz_min[d]))/dx_p[d];
                      weight[d][0]=1-weight[d][1];
                      ddi[d]=di[d];
                    }
                }

              if(dim==2)
                {
                  (*p)(c)=p_initial[ijk]*weight[0][0]*weight[1][0]
                    + p_initial[ijk+ddi[0]]*weight[0][1]*weight[1][0]
                    + p_initial[ijk+ddi[1]]*weight[0][0]*weight[1][1]
                    + p_initial[ijk+ddi[0]+ddi[1]]*weight[0][1]*weight[1][1];
                }
              else
                {
                  (*p)(c)=p_initial[ijk]*weight[0][0]*weight[1][0]*weight[2][0]
                    + p_initial[ijk+ddi[0]]*weight[0][1]*weight[1][0]*weight[2][0]
                    + p_initial[ijk+ddi[1]]*weight[0][0]*weight[1][1]*weight[2][0]
                    + p_initial[ijk+ddi[0]+ddi[1]]*weight[0][1]*weight[1][1]*weight[2][0]
                    
                    + p_initial[ijk+ddi[2]]*weight[0][0]*weight[1][0]*weight[2][1]
                    + p_initial[ijk+ddi[0]+ddi[2]]*weight[0][1]*weight[1][0]*weight[2][1]
                    + p_initial[ijk+ddi[1]+ddi[2]]*weight[0][0]*weight[1][1]*weight[2][1]
                    + p_initial[ijk+ddi[0]+ddi[1]+ddi[2]]*weight[0][1]*weight[1][1]*weight[2][1];
                }
            }
        }

      boost::shared_ptr<pdat::SideData<double> > v =
        boost::dynamic_pointer_cast<pdat::SideData<double> >
        (patch->getPatchData(v_id));
      v->fill(0.0);
    }
    d_stokes_fac_solver.set_boundaries(p_id,v_id,level,false);
  }

  fix_viscosity();

  d_stokes_fac_solver.initializeSolverState
    (p_id,cell_viscosity_id,edge_viscosity_id,dp_id,p_rhs_id,v_id,v_rhs_id,
     d_hierarchy,0,d_hierarchy->getFinestLevelNumber());

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
