#include "Elastic/FACOps.h"
#include "Constants.h"
/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void Elastic::FACOps::smooth_Tackley_2D
(SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int ln,
 int num_sweeps,
 double residual_tolerance)
{
  const int v_id(solution.getComponentDescriptorIndex(0)),
    v_rhs_id(residual.getComponentDescriptorIndex(0));

#ifdef DEBUG_CHECK_ASSERTIONS
  if (solution.getPatchHierarchy() != d_hierarchy
      || residual.getPatchHierarchy() != d_hierarchy)
    {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal hierarchy.");
    }
#endif
  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel>
    level = d_hierarchy->getPatchLevel(ln);

  /* Only need to sync the rhs once. This sync is needed because
     calculating a new pressure update requires computing in the ghost
     region so that the update for the velocity inside the box will be
     correct. */
  v_refine_patch_strategy.setTargetDataId(v_id);
  set_boundaries(v_id,level,true);
  xeqScheduleGhostFillNoCoarse(v_rhs_id,ln);

  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(v_id, ln);
  }

  double theta_momentum=1.0;

  /*
   * Smooth the number of sweeps specified or until
   * the convergence is satisfactory.
   */
  double maxres;
  /*
   * Instead of checking residual convergence globally, we check the
   * converged flag.  This avoids possible round-off errors affecting
   * different processes differently, leading to disagreement on
   * whether to continue smoothing.
   */
  const SAMRAI::hier::Index ip(1,0), jp(0,1);
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps*(1<<(d_ln_max-ln)) && !converged;
       ++sweep)
    {
      maxres=0;

      /* vx sweep */
      for(int rb=0;rb<2;++rb)
        {
          xeqScheduleGhostFillNoCoarse(v_id,ln);
          for (SAMRAI::hier::PatchLevel::Iterator pi(*level); pi; pi++)
            {
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch> patch = *pi;

              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> > v_ptr =
                patch->getPatchData(v_id);
              SAMRAI::pdat::SideData<double> &v(*v_ptr);
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
                patch->getPatchData(v_rhs_id);
              SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<double> >
                cell_visc_ptr= patch->getPatchData(cell_moduli_id);
              SAMRAI::pdat::CellData<double> &cell_moduli(*cell_visc_ptr);
              SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<double> >
                edge_visc_ptr= patch->getPatchData(edge_moduli_id);
              SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_visc_ptr);

              SAMRAI::hier::Box pbox=patch->getBox();
              SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
                geom = patch->getPatchGeometry();
              double dx = geom->getDx()[0];
              double dy = geom->getDx()[1];

              for(int j=pbox.lower(1); j<=pbox.upper(1); ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      SAMRAI::pdat::CellIndex
                        center(SAMRAI::tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      /* Update v */
                      smooth_V_2D(0,pbox,geom,center,ip,jp,
                                  v,v_rhs,maxres,dx,dy,cell_moduli,
                                  edge_moduli,theta_momentum);
                    }
                }
            }
          set_boundaries(v_id,level,true);
  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(v_id, ln);
  }
        }


      /* vy sweep */

      for(int rb=0;rb<2;++rb)
        {
          xeqScheduleGhostFillNoCoarse(v_id,ln);
          for (SAMRAI::hier::PatchLevel::Iterator pi(*level); pi; pi++)
            {
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch> patch = *pi;

              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> > v_ptr =
                patch->getPatchData(v_id);
              SAMRAI::pdat::SideData<double> &v(*v_ptr);
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
                patch->getPatchData(v_rhs_id);
              SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<double> >
                cell_visc_ptr= patch->getPatchData(cell_moduli_id);
              SAMRAI::pdat::CellData<double> &cell_moduli(*cell_visc_ptr);
              SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<double> >
                edge_visc_ptr= patch->getPatchData(edge_moduli_id);
              SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_visc_ptr);

              SAMRAI::hier::Box pbox=patch->getBox();
              SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
                geom = patch->getPatchGeometry();
              double dx = geom->getDx()[0];
              double dy = geom->getDx()[1];

              for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0); i+=2)
                    {
                      SAMRAI::pdat::CellIndex center(SAMRAI::tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      /* Update v */
                      smooth_V_2D(1,pbox,geom,center,jp,ip,
                                  v,v_rhs,maxres,dy,dx,cell_moduli,
                                  edge_moduli,theta_momentum);
                    }
                }
            }
          set_boundaries(v_id,level,true);
  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(v_id, ln);
  }
        }

      // if (residual_tolerance >= 0.0) {
        /*
         * Check for early end of sweeps due to convergence
         * only if it is numerically possible (user gave a
         * non negative value for residual tolerance).
         */
        converged = maxres < residual_tolerance;
        const SAMRAI::tbox::SAMRAI_MPI&
          mpi(d_hierarchy->getDomainMappedBoxLevel().getMPI());
        int tmp= converged ? 1 : 0;
        if (mpi.getSize() > 1)
          {
            mpi.AllReduce(&tmp, 1, MPI_MIN);
          }
        converged=(tmp==1);
        // if (d_enable_logging)
        //   SAMRAI::tbox::plog
        //     // << d_object_name << "\n"
        //     << "Tackley  " << ln << " " << sweep << " : " << maxres << "\n";
      // }
    }
}

