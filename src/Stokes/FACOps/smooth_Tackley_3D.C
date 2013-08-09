#include "Stokes/FACOps.h"
#include "Constants.h"
#include "Stokes/dRc_dp.h"
/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void Stokes::FACOps::smooth_Tackley_3D
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
  set_boundaries(p_id,v_id,level,true);
  xeqScheduleGhostFillNoCoarse(p_rhs_id,v_rhs_id,ln);

  p_refine_patch_strategy.setTargetDataId(p_id);
  v_refine_patch_strategy.setTargetDataId(v_id);
  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(p_id, v_id, ln);
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
  const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const SAMRAI::hier::Index pp[]={ip,jp,kp};
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps*(1<<(2*(d_ln_max-ln))) && !converged;
       ++sweep)
    {
      maxres=0;

      /* v sweeps */
      xeqScheduleGhostFillNoCoarse(p_id,invalid_id,ln);

      for(int ix=0;ix<3;++ix)
        for(int rb=0;rb<2;++rb)
          {
            xeqScheduleGhostFillNoCoarse(invalid_id,v_id,ln);
            for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
                 pi!=level->end(); ++pi)
              {
                boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;

                boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                  (patch->getPatchData(p_id));
                SAMRAI::pdat::CellData<double> &p(*p_ptr);
                
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
                boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_visc_ptr =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
                  (patch->getPatchData(edge_viscosity_id));
                SAMRAI::pdat::EdgeData<double> &edge_viscosity(*edge_visc_ptr);

                SAMRAI::hier::Box pbox=patch->getBox();
                boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                  boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
                  (patch->getPatchGeometry());
                const double *Dx = geom->getDx();

                for(int k=pbox.lower(2); k<=pbox.upper(2)+pp[ix][2]; ++k)
                  for(int j=pbox.lower(1); j<=pbox.upper(1)+pp[ix][1]; ++j)
                    {
                      /* Do the red-black skip */
                      int i_min=pbox.lower(0)
                        + (abs(pbox.lower(0) + j + k + rb))%2;
                      for(int i=i_min; i<=pbox.upper(0)+pp[ix][0]; i+=2)
                        {
                          SAMRAI::pdat::CellIndex center(SAMRAI::hier::Index(i,j,k));

                          /* Update v */
                          smooth_V_3D(ix,pbox,geom,p,v,v_rhs,cell_viscosity,
                                      edge_viscosity,center,
                                      Dx,theta_momentum,pp,maxres);
                        }
                    }
              }
            set_boundaries(invalid_id,v_id,level,true);
          }

      /* p sweep
         No need for red-black, because dp does not depend on
         the pressure. */
      xeqScheduleGhostFillNoCoarse(invalid_id,v_id,ln);

      for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
           pi!=level->end(); ++pi)
        {
          boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;

          boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(p_id));
          SAMRAI::pdat::CellData<double> &p(*p_ptr);
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dp_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(dp_id));
          SAMRAI::pdat::CellData<double> &dp(*dp_ptr);
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_rhs_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(p_rhs_id));
          SAMRAI::pdat::CellData<double> &p_rhs(*p_rhs_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch->getPatchData(v_id));
          SAMRAI::pdat::SideData<double> &v(*v_ptr);

          // boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
          //   patch->getPatchData(v_rhs_id);
          // SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(cell_viscosity_id));
          SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
          boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
            (patch->getPatchData(edge_viscosity_id));
          SAMRAI::pdat::EdgeData<double> &edge_viscosity(*edge_visc_ptr);

          SAMRAI::hier::Box pbox=patch->getBox();
          boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
            boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
            (patch->getPatchGeometry());
          const double *Dx = geom->getDx();

          SAMRAI::pdat::CellIterator
            cend(SAMRAI::pdat::CellGeometry::end(pbox));
          for(SAMRAI::pdat::CellIterator
                ci(SAMRAI::pdat::CellGeometry::begin(pbox)); ci!=cend; ++ci)
            {
              const SAMRAI::pdat::CellIndex &center(*ci);

              double delta_R_continuity=p_rhs(center);
              for(int ix=0;ix<3;++ix)
                {
                  const SAMRAI::pdat::SideIndex
                    x(center,ix,SAMRAI::pdat::SideIndex::Lower);
                  delta_R_continuity-=(v(x+pp[ix]) - v(x))/Dx[ix];;
                }

              /* No scaling here, though there should be. */
              maxres=std::max(maxres,std::fabs(delta_R_continuity));

              dp(center)=delta_R_continuity*theta_continuity
                /Stokes_dRc_dp_3D(pbox,center,cell_viscosity,edge_viscosity,v,Dx,pp);
              p(center)+=dp(center);
            }
        }
      set_boundaries(p_id,invalid_id,level,true);

      /* fix v sweep */
      xeqScheduleGhostFillNoCoarse(dp_id,invalid_id,ln);

      for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
           pi!=level->end(); ++pi)
        {
          boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;

          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dp_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(dp_id));
          SAMRAI::pdat::CellData<double> &dp(*dp_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            (patch->getPatchData(v_id));
          SAMRAI::pdat::SideData<double> &v(*v_ptr);
                
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(cell_viscosity_id));
          SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
          boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_visc_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
            (patch->getPatchData(edge_viscosity_id));
          SAMRAI::pdat::EdgeData<double> &edge_viscosity(*edge_visc_ptr);

          SAMRAI::hier::Box pbox=patch->getBox();
          boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
            boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
            (patch->getPatchGeometry());
          const double *Dx=geom->getDx();

          pbox.growUpper(SAMRAI::hier::IntVector::getOne(d_dim));

          SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(pbox));
          for(SAMRAI::pdat::CellIterator
                ci(SAMRAI::pdat::CellGeometry::begin(pbox)); ci!=cend; ++ci)
            {
              const SAMRAI::pdat::CellIndex &center(*ci);

              /* Update v */
              for(int ix=0;ix<3;++ix)
                {
                  const int iy((ix+1)%3), iz((ix+2)%3);
                  if(center[iy]<pbox.upper(iy) && center[iz]<pbox.upper(iz))
                    {
                      const SAMRAI::pdat::SideIndex
                        x(center,ix,SAMRAI::pdat::SideIndex::Lower);
                      const SAMRAI::pdat::EdgeIndex
                        edge_y(center,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
                        edge_z(center,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);

                      if(!((center[ix]==pbox.lower(ix)
                            && v(x-pp[ix])==boundary_value)
                           || (center[ix]==pbox.upper(ix)
                               && v(x+pp[ix])==boundary_value)))
                        v(x)+=(dp(center) - dp(center-pp[ix]))
                          /(Dx[ix]*Stokes_dRm_dv_3D(cell_viscosity,edge_viscosity,center,
                                                     center-pp[ix],edge_y+pp[iz],edge_y,
                                                     edge_z+pp[iy],edge_z,
                                                     Dx[ix],Dx[iy],Dx[iz]));
                    }
                }
            }
        }
      /* This is probably not necessary, since everyone always makes
         sure that everything is set before use. */
      set_boundaries(invalid_id,v_id,level,true);

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
      //     << "Tackley  " << ln << " " << sweep << " : " << maxres << "\n";
      // }
    }
}

