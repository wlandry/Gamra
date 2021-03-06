#include "Stokes/FACOps.hxx"
#include "Constants.hxx"
#include "Stokes/dRc_dp.hxx"
/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void Stokes::FACOps::smooth_Tackley_2D
(SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int level,
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
  boost::shared_ptr<SAMRAI::hier::PatchLevel> patch_level = d_hierarchy->getPatchLevel(level);

  /* Only need to sync the rhs once. This sync is needed because
     calculating a new pressure update requires computing in the ghost
     region so that the update for the velocity inside the box will be
     correct. */
  p_refine_patch_strategy.setTargetDataId(p_id);
  v_refine_patch_strategy.setTargetDataId(v_id);
  set_physical_boundaries(p_id,v_id,patch_level,true);
  xeqScheduleGhostFillNoCoarse(p_rhs_id,v_rhs_id,level);

  if (level > d_level_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(p_id, v_id, level);
  }

  double theta_momentum=0.7;
  double theta_continuity=1.0;

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
  for (int sweep=0; sweep < num_sweeps*(1<<(d_level_max-level)) && !converged;
       ++sweep)
    {
      maxres=0;

      /* vx sweep */
      xeqScheduleGhostFillNoCoarse(p_id,invalid_id,level);
      for(int rb=0;rb<2;++rb)
        {
          xeqScheduleGhostFillNoCoarse(invalid_id,v_id,level);
          for (SAMRAI::hier::PatchLevel::Iterator patch_iter(patch_level->begin());
               patch_iter!=patch_level->end(); ++patch_iter)
            {
              SAMRAI::hier::Patch &patch = **patch_iter;

              boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch.getPatchData(p_id));
              SAMRAI::pdat::CellData<double> &p(*p_ptr);
                
              boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                (patch.getPatchData(v_id));
              SAMRAI::pdat::SideData<double> &v(*v_ptr);
              boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                (patch.getPatchData(v_rhs_id));
              SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch.getPatchData(cell_viscosity_id));
              SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_visc_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                (patch.getPatchData(edge_viscosity_id));
              SAMRAI::pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              SAMRAI::hier::Box pbox=patch.getBox();
              boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
                (patch.getPatchGeometry());
              double dx = geom->getDx()[0];
              double dy = geom->getDx()[1];

              for(int j=pbox.lower(1); j<=pbox.upper(1); ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      SAMRAI::pdat::CellIndex center(SAMRAI::tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      /* Update v */
                      smooth_V_2D(0,pbox,geom,center,ip,jp,
                                  p,v,v_rhs,maxres,dx,dy,cell_viscosity,
                                  edge_viscosity,theta_momentum);
                    }
                }
            }
          set_physical_boundaries(invalid_id,v_id,patch_level,true);
        }


      /* vy sweep */

      for(int rb=0;rb<2;++rb)
        {
          xeqScheduleGhostFillNoCoarse(invalid_id,v_id,level);
          for (SAMRAI::hier::PatchLevel::Iterator patch_iter(patch_level->begin());
               patch_iter!=patch_level->end(); ++patch_iter)
            {
              SAMRAI::hier::Patch &patch = **patch_iter;

              boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch.getPatchData(p_id));
              SAMRAI::pdat::CellData<double> &p(*p_ptr);
                
              boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                (patch.getPatchData(v_id));
              SAMRAI::pdat::SideData<double> &v(*v_ptr);
              boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                (patch.getPatchData(v_rhs_id));
              SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                (patch.getPatchData(cell_viscosity_id));
              SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_visc_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                (patch.getPatchData(edge_viscosity_id));
              SAMRAI::pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              SAMRAI::hier::Box pbox=patch.getBox();
              boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
                (patch.getPatchGeometry());
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
                                  p,v,v_rhs,maxres,dy,dx,cell_viscosity,
                                  edge_viscosity,theta_momentum);
                    }
                }
            }
          set_physical_boundaries(invalid_id,v_id,patch_level,true);
        }



      /* p sweep
         No need for red-black, because dp does not depend on
         the pressure. */
      xeqScheduleGhostFillNoCoarse(invalid_id,v_id,level);

      for (SAMRAI::hier::PatchLevel::Iterator patch_iter(patch_level->begin());
           patch_iter!=patch_level->end(); ++patch_iter)
        {
          SAMRAI::hier::Patch &patch = **patch_iter;

          boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(p_id));
          SAMRAI::pdat::CellData<double> &p(*p_ptr);
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dp_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(dp_id));
          SAMRAI::pdat::CellData<double> &dp(*dp_ptr);
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_rhs_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(p_rhs_id));
          SAMRAI::pdat::CellData<double> &p_rhs(*p_rhs_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch.getPatchData(v_id));
          SAMRAI::pdat::SideData<double> &v(*v_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(cell_viscosity_id));
          SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
          boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
            (patch.getPatchData(edge_viscosity_id));
          SAMRAI::pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

          SAMRAI::hier::Box pbox=patch.getBox();
          boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
            boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
            (patch.getPatchGeometry());
          double dx = geom->getDx()[0];
          double dy = geom->getDx()[1];

          SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(pbox));
          for(SAMRAI::pdat::CellIterator
                ci(SAMRAI::pdat::CellGeometry::begin(pbox)); ci!=cend; ++ci)
            {
              const SAMRAI::pdat::CellIndex &center(*ci);
              const SAMRAI::pdat::SideIndex
                x(center,0,SAMRAI::pdat::SideIndex::Lower),
                y(center,1,SAMRAI::pdat::SideIndex::Lower);

              /* Update p */
              double dvx_dx=(v(x+ip) - v(x))/dx;
              double dvy_dy=(v(y+jp) - v(y))/dy;

              double delta_R_continuity=
                p_rhs(center) - dvx_dx - dvy_dy;

              /* No scaling here, though there should be. */
              maxres=std::max(maxres,std::fabs(delta_R_continuity));

              dp(center)=delta_R_continuity*theta_continuity
                /Stokes_dRc_dp_2D(pbox,center,x,y,cell_viscosity,edge_viscosity,v,dx,dy);
              p(center)+=dp(center);
            }
        }
      set_physical_boundaries(p_id,invalid_id,patch_level,true);


      /* fix v sweep */
      xeqScheduleGhostFillNoCoarse(dp_id,invalid_id,level);

      for (SAMRAI::hier::PatchLevel::Iterator p(patch_level->begin());
           p!=patch_level->end(); ++p)
        {
          SAMRAI::hier::Patch &patch = **p;

          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dp_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(dp_id));
          SAMRAI::pdat::CellData<double> &dp(*dp_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch.getPatchData(v_id));
          SAMRAI::pdat::SideData<double> &v(*v_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(cell_viscosity_id));
          SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
          boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
            (patch.getPatchData(edge_viscosity_id));
          SAMRAI::pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

          SAMRAI::hier::Box pbox=patch.getBox();
          boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
            boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
            (patch.getPatchGeometry());
          double dx = geom->getDx()[0];
          double dy = geom->getDx()[1];

          pbox.growUpper(SAMRAI::hier::IntVector::getOne(d_dim));

          SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(pbox));
          for(SAMRAI::pdat::CellIterator
                ci(SAMRAI::pdat::CellGeometry::begin(pbox)); ci!=cend; ++ci)
            {
              const SAMRAI::pdat::CellIndex &center(*ci);

              const SAMRAI::pdat::SideIndex
                x(center,0,SAMRAI::pdat::SideIndex::Lower),
                y(center,1,SAMRAI::pdat::SideIndex::Lower);
              const SAMRAI::pdat::NodeIndex
                edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);

              /* Update v */
              if(center[1]<pbox.upper(1))
                {
                  if(!((center[0]==pbox.lower(0) && v(x-ip)==boundary_value)
                       || (center[0]==pbox.upper(0)
                           && v(x+ip)==boundary_value)))
                    v(x)+=(dp(center) - dp(center-ip))
                      /(dx*Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,center,
                                             center-ip,edge+jp,edge,dx,dy));
                }
              if(center[0]<pbox.upper(0))
                {
                  if(!((center[1]==pbox.lower(1) && v(y-jp)==boundary_value)
                       || (center[1]==pbox.upper(1)
                           && v(y+jp)==boundary_value)))
                    v(y)+=(dp(center) - dp(center-jp))
                      /(dy*Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,center,
                                             center-jp,edge+ip,edge,dy,dx));
                }
            }
        }
      set_physical_boundaries(invalid_id,v_id,patch_level,true);

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
      // if (d_enable_logging)
      //   SAMRAI::tbox::plog
      //     // << d_object_name << "\n"
      //     << "Tackley  " << level << " " << sweep << " : " << maxres << "\n";
      // }
    }
}

