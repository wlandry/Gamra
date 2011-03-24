#include "StokesFACOps.h"
#include "Boundary.h"
#include "dRc_dp.h"
/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void SAMRAI::solv::StokesFACOps::smooth_Tackley_3D
(SAMRAIVectorReal<double>& solution,
 const SAMRAIVectorReal<double>& residual,
 int ln,
 int num_sweeps,
 double residual_tolerance)
{
  const int p_id(solution.getComponentDescriptorIndex(0)),
    p_rhs_id(residual.getComponentDescriptorIndex(0)),
    v_id(solution.getComponentDescriptorIndex(1)),
    v_rhs_id(residual.getComponentDescriptorIndex(1));

  checkInputPatchDataIndices();

#ifdef DEBUG_CHECK_ASSERTIONS
  if (solution.getPatchHierarchy() != d_hierarchy
      || residual.getPatchHierarchy() != d_hierarchy)
    {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal hierarchy.");
    }
#endif
  tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

  /* Only need to sync the rhs once. This sync is needed because
     calculating a new pressure update requires computing in the ghost
     region so that the update for the velocity inside the box will be
     correct. */
  set_boundaries(v_id,level,true);
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
  const hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const hier::Index pp[]={ip,jp,kp};
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps*(1<<(d_ln_max-ln)) && !converged;
       ++sweep)
    {
      maxres=0;

      /* v sweeps */
      xeqScheduleGhostFillNoCoarse(p_id,invalid_id,ln);

      for(int ix=0;ix<3;++ix)
        for(int rb=0;rb<2;++rb)
          {
            xeqScheduleGhostFillNoCoarse(invalid_id,v_id,ln);
            for (hier::PatchLevel::Iterator pi(*level); pi; pi++)
              {
                tbox::Pointer<hier::Patch> patch = *pi;

                tbox::Pointer<pdat::CellData<double> > p_ptr =
                  patch->getPatchData(p_id);
                pdat::CellData<double> &p(*p_ptr);
                
                tbox::Pointer<pdat::SideData<double> > v_ptr =
                  patch->getPatchData(v_id);
                pdat::SideData<double> &v(*v_ptr);
                tbox::Pointer<pdat::SideData<double> > v_rhs_ptr =
                  patch->getPatchData(v_rhs_id);
                pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
                tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
                  = patch->getPatchData(cell_viscosity_id);
                pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
                tbox::Pointer<pdat::EdgeData<double> > edge_visc_ptr
                  = patch->getPatchData(edge_viscosity_id);
                pdat::EdgeData<double> &edge_viscosity(*edge_visc_ptr);

                hier::Box pbox=patch->getBox();
                tbox::Pointer<geom::CartesianPatchGeometry>
                  geom = patch->getPatchGeometry();
                const double *Dx = geom->getDx();

                for(int k=pbox.lower(2); k<=pbox.upper(2)+pp[ix][2]; ++k)
                  for(int j=pbox.lower(1); j<=pbox.upper(1)+pp[ix][1]; ++j)
                    {
                      /* Do the red-black skip */
                      int i_min=pbox.lower(0)
                        + (abs(pbox.lower(0) + j + k + rb))%2;
                      for(int i=i_min; i<=pbox.upper(0)+pp[ix][0]; i+=2)
                        {
                          pdat::CellIndex center(hier::Index(i,j,k));

                          /* Update v */
                          smooth_V_3D(ix,pbox,geom,p,v,v_rhs,cell_viscosity,
                                      edge_viscosity,center,
                                      Dx,theta_momentum,pp,maxres);
                        }
                    }
              }
            set_boundaries(v_id,level,true);
          }

      /* p sweep
         No need for red-black, because dp does not depend on
         the pressure. */
      xeqScheduleGhostFillNoCoarse(invalid_id,v_id,ln);

      for (hier::PatchLevel::Iterator pi(*level); pi; pi++)
        {
          tbox::Pointer<hier::Patch> patch = *pi;

          tbox::Pointer<pdat::CellData<double> > p_ptr =
            patch->getPatchData(p_id);
          pdat::CellData<double> &p(*p_ptr);
          tbox::Pointer<pdat::CellData<double> > dp_ptr =
            patch->getPatchData(dp_id);
          pdat::CellData<double> &dp(*dp_ptr);
          tbox::Pointer<pdat::CellData<double> > p_rhs_ptr =
            patch->getPatchData(p_rhs_id);
          pdat::CellData<double> &p_rhs(*p_rhs_ptr);
                
          tbox::Pointer<pdat::SideData<double> > v_ptr =
            patch->getPatchData(v_id);
          pdat::SideData<double> &v(*v_ptr);

          // tbox::Pointer<pdat::SideData<double> > v_rhs_ptr =
          //   patch->getPatchData(v_rhs_id);
          // pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
          tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
            = patch->getPatchData(cell_viscosity_id);
          pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
          tbox::Pointer<pdat::EdgeData<double> > edge_visc_ptr
            = patch->getPatchData(edge_viscosity_id);
          pdat::EdgeData<double> &edge_viscosity(*edge_visc_ptr);

          hier::Box pbox=patch->getBox();
          tbox::Pointer<geom::CartesianPatchGeometry>
            geom = patch->getPatchGeometry();
          const double *Dx = geom->getDx();

          for(pdat::CellIterator ci(pbox); ci; ci++)
            {
              pdat::CellIndex center(*ci);

              // tbox::plog << "smooth "
              //            << sweep << " "
              //            << center[0] << " "
              //            << center[1] << " "
              //            << center[2] << " ";

              double delta_R_continuity=p_rhs(center);
              for(int ix=0;ix<3;++ix)
                {
                  const pdat::SideIndex x(center,ix,pdat::SideIndex::Lower);
                  delta_R_continuity-=(v(x+pp[ix]) - v(x))/Dx[ix];;

                  // tbox::plog << v(x) << " "
                  //            << v_rhs(x) << " ";
                }

              /* No scaling here, though there should be. */
              maxres=std::max(maxres,std::fabs(delta_R_continuity));

              dp(center)=delta_R_continuity*theta_continuity
                /dRc_dp_3D(pbox,center,cell_viscosity,edge_viscosity,v,Dx,pp);
              p(center)+=dp(center);
              // tbox::plog << p(center) << " "
              //            << dp(center) << " "
              //            << p_rhs(center) << " "
              //            // << maxres << " "
              //            << "\n";
            }
        }


      /* fix v sweep */
      xeqScheduleGhostFillNoCoarse(dp_id,invalid_id,ln);

      for (hier::PatchLevel::Iterator pi(*level); pi; pi++)
        {
          tbox::Pointer<hier::Patch> patch = *pi;

          tbox::Pointer<pdat::CellData<double> > dp_ptr =
            patch->getPatchData(dp_id);
          pdat::CellData<double> &dp(*dp_ptr);
                
          tbox::Pointer<pdat::SideData<double> > v_ptr =
            patch->getPatchData(v_id);
          pdat::SideData<double> &v(*v_ptr);
                
          tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
            = patch->getPatchData(cell_viscosity_id);
          pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
          tbox::Pointer<pdat::EdgeData<double> > edge_visc_ptr
            = patch->getPatchData(edge_viscosity_id);
          pdat::EdgeData<double> &edge_viscosity(*edge_visc_ptr);

          hier::Box pbox=patch->getBox();
          tbox::Pointer<geom::CartesianPatchGeometry>
            geom = patch->getPatchGeometry();
          const double *Dx=geom->getDx();

          pbox.growUpper(hier::IntVector::getOne(d_dim));

          for(pdat::CellIterator ci(pbox); ci; ci++)
            {
              pdat::CellIndex center(*ci);

              /* Update v */
              for(int ix=0;ix<3;++ix)
                {
                  const int iy((ix+1)%3), iz((ix+2)%3);
                  if(center[iy]<pbox.upper(iy) && center[iz]<pbox.upper(iz))
                    {
                      const pdat::SideIndex x(center,ix,pdat::SideIndex::Lower);
                      const pdat::EdgeIndex
                        edge_y(center,iy,pdat::EdgeIndex::LowerLeft),
                        edge_z(center,iz,pdat::EdgeIndex::LowerLeft);

                      if(!((center[ix]==pbox.lower(ix)
                            && v(x-pp[ix])==boundary_value)
                           || (center[ix]==pbox.upper(ix)
                               && v(x+pp[ix])==boundary_value)))
                        v(x)+=(dp(center) - dp(center-pp[ix]))
                          /(Dx[ix]*dRm_dv_3D(cell_viscosity,edge_viscosity,center,
                                             center-pp[ix],edge_y+pp[ix],edge_y,
                                             edge_z+pp[ix],edge_z,
                                             Dx[ix],Dx[iy],Dx[iz]));
                    }
                }
            }
          set_boundaries(v_id,level,true);
        }
      // if (residual_tolerance >= 0.0) {
        /*
         * Check for early end of sweeps due to convergence
         * only if it is numerically possible (user gave a
         * non negative value for residual tolerance).
         */
        converged = maxres < residual_tolerance;
        const tbox::SAMRAI_MPI&
          mpi(d_hierarchy->getDomainMappedBoxLevel().getMPI());
        int tmp= converged ? 1 : 0;
        if (mpi.getSize() > 1)
          {
            mpi.AllReduce(&tmp, 1, MPI_MIN);
          }
        converged=(tmp==1);
        if (d_enable_logging)
          tbox::plog
            // << d_object_name << "\n"
            << "Tackley  " << ln << " " << sweep << " : " << maxres << "\n";
      // }
    }
}

