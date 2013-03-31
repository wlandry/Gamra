#include "Elastic/FACOps.h"
#include "Elastic/V_Boundary_Refine.h"
#include "Constants.h"
#include "Elastic/dRm_dv.h"

/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void Elastic::FACOps::smooth_2D
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
  boost::shared_ptr<SAMRAI::hier::PatchLevel>
    level = d_hierarchy->getPatchLevel(ln);

  v_refine_patch_strategy.data_id=v_id;
  v_refine_patch_strategy.is_residual=true;
  V_Boundary_Refine::is_residual=true;
  xeqScheduleGhostFillNoCoarse(v_rhs_id,ln);

  if (ln > d_ln_min) {
    xeqScheduleGhostFill(v_id, ln);
  }

  double theta_momentum=1.0;

  /*
   * Smooth the number of sweeps specified or until the convergence is
   * satisfactory.
   */

  const SAMRAI::hier::Index unit[]={SAMRAI::hier::Index(1,0),SAMRAI::hier::Index(0,1)};
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps*(1<<(d_ln_max-ln)) && !converged;
       ++sweep)
    {
      double max_residual(0);

      const int dim(2);
      for(int ix=0; ix<dim; ++ix)
        {
          const int iy((ix+1)%dim);
          const SAMRAI::hier::Index ip(unit[ix]), jp(unit[iy]);
          for(int rb=0;rb<2;++rb)
            {
              xeqScheduleGhostFillNoCoarse(v_id,ln);
              if (ln > d_ln_min)
                {
                  xeqScheduleGhostFill(v_id, ln);
                }
              set_boundaries(v_id,level,true);
              for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
                   pi!=level->end(); ++pi)
                {
                  boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;

                  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
                    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                    (patch->getPatchData(v_id));
                  SAMRAI::pdat::SideData<double> &v(*v_ptr);
                  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_ptr =
                    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
                    (patch->getPatchData(v_rhs_id));
                  SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
                  boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli_ptr=
                    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                    (patch->getPatchData(cell_moduli_id));
                  SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
                  boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_moduli_ptr=
                    boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                    (patch->getPatchData(edge_moduli_id));
                  SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

                  SAMRAI::hier::Box pbox=patch->getBox();
                  boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
                    (patch->getPatchGeometry());
                  double dx = geom->getDx()[ix];
                  double dy = geom->getDx()[iy];

                  if(have_embedded_boundary())
                    {
                      boost::shared_ptr<SAMRAI::pdat::SideData<double> >
                        level_set_ptr =boost::dynamic_pointer_cast
                        <SAMRAI::pdat::SideData<double> >
                        (patch->getPatchData(level_set_id));
                      SAMRAI::pdat::SideData<double> &level_set(*level_set_ptr);
                      for(int j=pbox.lower(1); j<=pbox.upper(1)+unit[ix][1]; ++j)
                        {
                          /* Do the red-black skip */
                          int i_min=pbox.lower(0)
                            + (abs(pbox.lower(0) + j + rb))%2;
                          for(int i=i_min; i<=pbox.upper(0)+unit[ix][0]; i+=2)
                            {
                              SAMRAI::pdat::CellIndex
                                center(SAMRAI::hier::Index(i,j));
                              const SAMRAI::pdat::SideIndex
                                x(center,ix,SAMRAI::pdat::SideIndex::Lower),
                                y(center,iy,SAMRAI::pdat::SideIndex::Lower);
                              if(level_set(x)>1)
                                {
                                  smooth_V_2D(ix,pbox,center,unit[ix],unit[iy],
                                              v,v_rhs,max_residual,dx,dy,
                                              cell_moduli,edge_moduli,
                                              theta_momentum);
                                }
                              else if(level_set(x)>0)
                                {
                                  const SAMRAI::pdat::CellIndex cell(x);
                                  const SAMRAI::pdat::NodeIndex
                                    edge(x,SAMRAI::pdat::NodeIndex::LowerLeft);

                                  double delta_Rx=v_rhs(x)
                                    - v_level_set_operator_2D(level_set,v,
                                                              cell_moduli,
                                                              edge_moduli,cell,
                                                              edge,x,y,ip,jp,
                                                              dx,dy);

                                  max_residual=std::max(max_residual,
                                                        std::fabs(delta_Rx));
                                  double C_vx(dRm_dv_2D(cell_moduli,edge_moduli,
                                                        cell,cell-ip,
                                                        edge+jp,edge,dx,dy));

                                  v(x)+=delta_Rx*theta_momentum/C_vx;
                                }
                            }
                        }
                    }
                  else
                    {
                      for(int j=pbox.lower(1); j<=pbox.upper(1)+unit[ix][1]; ++j)
                        {
                          /* Do the red-black skip */
                          int i_min=pbox.lower(0)
                            + (abs(pbox.lower(0) + j + rb))%2;
                          for(int i=i_min; i<=pbox.upper(0)+unit[ix][0]; i+=2)
                            {
                              SAMRAI::pdat::CellIndex
                                center(SAMRAI::hier::Index(i,j));
                              smooth_V_2D(ix,pbox,center,unit[ix],unit[iy],
                                          v,v_rhs,max_residual,dx,dy,cell_moduli,
                                          edge_moduli,theta_momentum);
                            }
                        }
                    }
                }
            }
        }

      if (residual_tolerance >= 0.0)
        {
          /*
           * Check for early end of sweeps due to convergence
           * only if it is numerically possible (user gave a
           * non negative value for residual tolerance).
           */
          /*
           * Instead of checking residual convergence globally, we check the
           * converged flag.  This avoids possible round-off errors affecting
           * different processes differently, leading to disagreement on
           * whether to continue smoothing.
           */
          converged = max_residual < residual_tolerance;
          const SAMRAI::tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
          int tmp= converged ? 1 : 0;
          if (mpi.getSize() > 1)
            {
              mpi.AllReduce(&tmp, 1, MPI_MIN);
            }
          converged=(tmp==1);
        }

      // SAMRAI::tbox::plog
      //   << "smooth_2D  " << ln << " " << sweep << " : " << max_residual << "\n";
    }

  xeqScheduleGhostFillNoCoarse(v_id,ln);
  if (ln > d_ln_min)
    {
      xeqScheduleGhostFill(v_id, ln);
    }
  set_boundaries(v_id,level,true);
}

