#include "Stokes/FACOps.hxx"
#include "Constants.hxx"
/* Smooth the error using a red-black Gauss-Seidel smoother inspired
   by Introduction to Numerical Geodynamic Modelling, Taras Gerya,
   2010

   This does not give the same answers in serial and parallel, because
   it smooths both vx and vy at the same time.  The stencil for the vx
   uses diagonal elements from vy and vice-versa.  So the red and
   black updates are not entirely separate.  */

void Stokes::FACOps::smooth_Gerya
(SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int ln,
 int num_sweeps,
 double residual_tolerance)
{
  const int p_id(solution.getComponentDescriptorIndex(0)),
    p_rhs_id(residual.getComponentDescriptorIndex(0)),
    v_id(solution.getComponentDescriptorIndex(1)),
    v_rhs_id(residual.getComponentDescriptorIndex(1));

#ifdef DEBUG_CHECK_ASSERTIONS
  if (solution.getPatchHierarchy() != d_hierarchy
      || residual.getPatchHierarchy() != d_hierarchy)
    {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal hierarchy.");
    }
#endif
  boost::shared_ptr<SAMRAI::hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

  /* Only need to sync the rhs once. This sync is needed because
     calculating a new pressure update requires computing in the ghost
     region so that the update for the velocity inside the box will be
     correct. */
  p_refine_patch_strategy.setTargetDataId(p_id);
  v_refine_patch_strategy.setTargetDataId(v_id);
  set_physical_boundaries(p_id,v_id,level,true);
  xeqScheduleGhostFillNoCoarse(p_rhs_id,v_rhs_id,ln);

  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(p_id, v_id, ln);
  }

  double theta_momentum=1.2;
  double theta_continuity=0.3;

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
      for(int rb=0;rb<2;++rb)
        {
          xeqScheduleGhostFillNoCoarse(p_id,v_id,ln);
          for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
               pi!=level->end(); ++pi)
            {
              boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;

              boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch->getPatchData(p_id));
              SAMRAI::pdat::CellData<double> &p(*p_ptr);
              boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_rhs_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch->getPatchData(p_rhs_id));
              SAMRAI::pdat::CellData<double> &p_rhs(*p_rhs_ptr);
                
              boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                (patch->getPatchData(v_id));
              SAMRAI::pdat::SideData<double> &v(*v_ptr);
              boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                (patch->getPatchData(v_rhs_id));
              SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch->getPatchData(cell_viscosity_id));
              SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_visc_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                (patch->getPatchData(edge_viscosity_id));
              SAMRAI::pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              SAMRAI::hier::Box pbox=patch->getBox();
              boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
                (patch->getPatchGeometry());
              double dx = geom->getDx()[0];
              double dy = geom->getDx()[1];

              for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      SAMRAI::pdat::CellIndex center(SAMRAI::tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      SAMRAI::pdat::CellIndex up(center), down(center), right(center),
                        left(center);

                      ++up[1];
                      --down[1];
                      ++right[0];
                      --left[0];

                      const SAMRAI::pdat::SideIndex
                        x(center,0,SAMRAI::pdat::SideIndex::Lower),
                        y(center,1,SAMRAI::pdat::SideIndex::Lower);
                      const SAMRAI::pdat::NodeIndex
                        edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);

                      /* Update p */
                      if(j<pbox.upper(1)+1 && i<pbox.upper(0)+1)
                        {
                          double dvx_dx=(v(x+ip)-v(x))/dx;
                          double dvy_dy=(v(y+jp)-v(y))/dy;

                          double delta_R_continuity=p_rhs(center)-dvx_dx-dvy_dy;

                          /* No scaling here, though there should be. */
                          maxres=std::max(maxres,std::fabs(delta_R_continuity));

                          p(center)+=cell_viscosity(center)
                            *delta_R_continuity*theta_continuity;
                        }
                      /* Update v */
                      if(j<pbox.upper(1)+1)
                        {
                          smooth_V_2D(0,pbox,geom,center,ip,jp,
                                      p,v,v_rhs,maxres,dx,dy,cell_viscosity,
                                      edge_viscosity,theta_momentum);
                        }
                      if(i<pbox.upper(0)+1)
                        {
                          smooth_V_2D(1,pbox,geom,center,jp,ip,
                                      p,v,v_rhs,maxres,dy,dx,cell_viscosity,
                                      edge_viscosity,theta_momentum);
                        }
                    }
                }
            }
          set_physical_boundaries(p_id,v_id,level,true);
        }
      // if (residual_tolerance >= 0.0) {
        /*
         * Check for early end of sweeps due to convergence
         * only if it is numerically possible (user gave a
         * non negative value for residual tolerance).
         */
        converged = maxres < residual_tolerance;
        const SAMRAI::tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
        int tmp= converged ? 1 : 0;
        if (mpi.getSize() > 1)
          {
            mpi.AllReduce(&tmp, 1, MPI_MIN);
          }
        converged=(tmp==1);
        if (d_enable_logging)
          SAMRAI::tbox::plog
            // << d_object_name << "\n"
            << "Gerya " << ln << " " << sweep << " : " << maxres << "\n";
      // }
    }
}

