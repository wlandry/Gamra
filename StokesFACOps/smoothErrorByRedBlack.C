/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "StokesFACOps.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "StokesHypreSolver.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

#include "Boundary.h"
/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void SAMRAI::solv::StokesFACOps::smoothErrorByRedBlack
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
  p_refine_patch_strategy.setTargetDataId(p_id);
  v_refine_patch_strategy.setTargetDataId(v_id);
  set_boundaries(v_id,level,true);
  xeqScheduleGhostFillNoCoarse(p_rhs_id,v_rhs_id,ln);

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

  const hier::Index ip(1,0), jp(0,1);

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
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps*(1<<(d_ln_max-ln)) && !converged; ++sweep)
    {

      maxres=0;

      /* vx sweep */
      for(int rb=0;rb<2;++rb)
        {
          // Need to sync
          xeqScheduleGhostFillNoCoarse(p_id,v_id,ln);
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
              tbox::Pointer<pdat::SideData<double> > v_rhs_ptr =
                patch->getPatchData(v_rhs_id);
              pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
                = patch->getPatchData(cell_viscosity_id);
              pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              tbox::Pointer<pdat::NodeData<double> > edge_visc_ptr
                = patch->getPatchData(edge_viscosity_id);
              pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              hier::Box pbox=patch->getBox();
              tbox::Pointer<geom::CartesianPatchGeometry>
                geom = patch->getPatchGeometry();
              double dx = *(geom->getDx());
              double dy = *(geom->getDx());

              /* Set an array of bools that tells me whether a point
                 should set the pressure or just let it be.  This is
                 needed at coarse/fine boundaries where the pressure
                 is fixed. */
              hier::Box gbox=p.getGhostBox();
              std::vector<bool> set_p(gbox.size(),true);

              const tbox::Array<hier::BoundaryBox >&edges
                =d_cf_boundary[ln]->getEdgeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<edges.size(); ++mm)
                for(int j=edges[mm].getBox().lower(1);
                    j<=edges[mm].getBox().upper(1); ++j)
                  for(int i=edges[mm].getBox().lower(0);
                      i<=edges[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              const tbox::Array<hier::BoundaryBox >&nodes
                =d_cf_boundary[ln]->getNodeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<nodes.size(); ++mm)
                for(int j=nodes[mm].getBox().lower(1);
                    j<=nodes[mm].getBox().upper(1); ++j)
                  for(int i=nodes[mm].getBox().lower(0);
                      i<=nodes[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(0,0))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;
                  
              if(geom->getTouchesRegularBoundary(0,1))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(1,0))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[i-gbox.lower(0)]=false;

              if(geom->getTouchesRegularBoundary(1,1))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[(i-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(gbox.upper(1)-gbox.lower(1))]=
                    false;

              for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      pdat::CellIndex center(tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      pdat::CellIndex up(center), down(center), right(center),
                        left(center);

                      ++up[1];
                      --down[1];
                      ++right[0];
                      --left[0];

                      /* Update v */
                      if(set_p[(i-gbox.lower(0))
                               + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]
                         || (i==pbox.upper(0)+1 && j<pbox.upper(1)+1))
                        {
                          Update_V(0,j,pbox,geom,center,left,right,down,up,p,
                                   v,v_rhs,maxres,dx,dy,cell_viscosity,
                                   edge_viscosity,theta_momentum);
                        }
                    }
                }
            }
          set_boundaries(v_id,level,true);
        }


      /* vy sweep */
      for(int rb=0;rb<2;++rb)
        {
          // Need to sync
          xeqScheduleGhostFillNoCoarse(p_id,v_id,ln);
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
              tbox::Pointer<pdat::SideData<double> > v_rhs_ptr =
                patch->getPatchData(v_rhs_id);
              pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
                = patch->getPatchData(cell_viscosity_id);
              pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              tbox::Pointer<pdat::NodeData<double> > edge_visc_ptr
                = patch->getPatchData(edge_viscosity_id);
              pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              hier::Box pbox=patch->getBox();
              tbox::Pointer<geom::CartesianPatchGeometry>
                geom = patch->getPatchGeometry();
              double dx = *(geom->getDx());
              double dy = *(geom->getDx());

              /* Set an array of bools that tells me whether a point
                 should set the pressure or just let it be.  This is
                 needed at coarse/fine boundaries where the pressure
                 is fixed. */
              hier::Box gbox=p.getGhostBox();
              std::vector<bool> set_p(gbox.size(),true);

              const tbox::Array<hier::BoundaryBox >&edges
                =d_cf_boundary[ln]->getEdgeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<edges.size(); ++mm)
                for(int j=edges[mm].getBox().lower(1);
                    j<=edges[mm].getBox().upper(1); ++j)
                  for(int i=edges[mm].getBox().lower(0);
                      i<=edges[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              const tbox::Array<hier::BoundaryBox >&nodes
                =d_cf_boundary[ln]->getNodeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<nodes.size(); ++mm)
                for(int j=nodes[mm].getBox().lower(1);
                    j<=nodes[mm].getBox().upper(1); ++j)
                  for(int i=nodes[mm].getBox().lower(0);
                      i<=nodes[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(0,0))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;
                  
              if(geom->getTouchesRegularBoundary(0,1))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(1,0))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[i-gbox.lower(0)]=false;

              if(geom->getTouchesRegularBoundary(1,1))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[(i-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(gbox.upper(1)-gbox.lower(1))]=
                    false;

              for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      pdat::CellIndex center(tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      pdat::CellIndex up(center), down(center), right(center),
                        left(center);

                      ++up[1];
                      --down[1];
                      ++right[0];
                      --left[0];

                      /* Update v */
                      if(set_p[(i-gbox.lower(0))
                               + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]
                         || (i<pbox.upper(0)+1 && j==pbox.upper(1)+1))
                        {
                          Update_V(1,i,pbox,geom,center,down,up,left,right,p,
                                   v,v_rhs,maxres,dy,dx,cell_viscosity,
                                   edge_viscosity,theta_momentum);
                        }
                    }
                }
            }
          set_boundaries(v_id,level,true);
        }



      /* p sweep */
      for(int rb=0;rb<2;++rb)
        {
          // Need to sync
          xeqScheduleGhostFillNoCoarse(-1,v_id,ln);
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
              tbox::Pointer<pdat::SideData<double> > v_rhs_ptr =
                patch->getPatchData(v_rhs_id);
              pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
                = patch->getPatchData(cell_viscosity_id);
              pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              tbox::Pointer<pdat::NodeData<double> > edge_visc_ptr
                = patch->getPatchData(edge_viscosity_id);
              pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              hier::Box pbox=patch->getBox();
              tbox::Pointer<geom::CartesianPatchGeometry>
                geom = patch->getPatchGeometry();
              double dx = *(geom->getDx());
              double dy = *(geom->getDx());

              /* Set an array of bools that tells me whether a point
                 should set the pressure or just let it be.  This is
                 needed at coarse/fine boundaries where the pressure
                 is fixed. */
              hier::Box gbox=p.getGhostBox();
              std::vector<bool> set_p(gbox.size(),true);

              const tbox::Array<hier::BoundaryBox >&edges
                =d_cf_boundary[ln]->getEdgeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<edges.size(); ++mm)
                for(int j=edges[mm].getBox().lower(1);
                    j<=edges[mm].getBox().upper(1); ++j)
                  for(int i=edges[mm].getBox().lower(0);
                      i<=edges[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              const tbox::Array<hier::BoundaryBox >&nodes
                =d_cf_boundary[ln]->getNodeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<nodes.size(); ++mm)
                for(int j=nodes[mm].getBox().lower(1);
                    j<=nodes[mm].getBox().upper(1); ++j)
                  for(int i=nodes[mm].getBox().lower(0);
                      i<=nodes[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(0,0))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;
                  
              if(geom->getTouchesRegularBoundary(0,1))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(1,0))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[i-gbox.lower(0)]=false;

              if(geom->getTouchesRegularBoundary(1,1))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[(i-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(gbox.upper(1)-gbox.lower(1))]=
                    false;

              for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      pdat::CellIndex center(tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      pdat::CellIndex up(center), down(center), right(center),
                        left(center);

                      ++up[1];
                      --down[1];
                      ++right[0];
                      --left[0];

                      const pdat::NodeIndex
                        center_e(center,pdat::NodeIndex::LowerLeft),
                        up_e(up,pdat::NodeIndex::LowerLeft),
                        right_e(right,pdat::NodeIndex::LowerLeft);

                      /* Update p */
                      if(set_p[(i-gbox.lower(0))
                               + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))])
                        {
                          double dvx_dx=
                            (v(pdat::SideIndex(center,pdat::SideIndex::X,
                                                  pdat::SideIndex::Upper))
                             - v(pdat::SideIndex(center,pdat::SideIndex::X,
                                                    pdat::SideIndex::Lower)))/dx;
                          double dvy_dy=
                            (v(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                  pdat::SideIndex::Upper))
                             - v(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                    pdat::SideIndex::Lower)))/dy;

                          double delta_R_continuity=
                            p_rhs(center) - dvx_dx - dvy_dy;

                          /* No scaling here, though there should be. */
                          maxres=std::max(maxres,std::fabs(delta_R_continuity));

                          const double dRm_dp_xp(1/dx), dRm_dp_xm(-1/dx),
                            dRm_dp_yp(1/dy), dRm_dp_ym(-1/dy),
                            dRc_dvx_p(-1/dx), dRc_dvx_m(1/dx),
                            dRc_dvy_p(-1/dy), dRc_dvy_m(1/dy);

                          const double dRm_dvx_p =
                            dRm_dv(cell_viscosity,edge_viscosity,
                                   right,center,up_e+ip,center_e+ip,dx,dy);

                          const double dRm_dvx_m =
                            dRm_dv(cell_viscosity,edge_viscosity,
                                   center,left,up_e,center_e,dx,dy);

                          const double dRm_dvy_p =
                            dRm_dv(cell_viscosity,edge_viscosity,
                                   up,center,right_e+jp,center_e+jp,dy,dx);

                          const double dRm_dvy_m =
                            dRm_dv(cell_viscosity,edge_viscosity,
                                   center,down,right_e,center_e,dy,dx);

                          const double dRc_dp=dRc_dvx_p * dRm_dp_xp/dRm_dvx_p
                            + dRc_dvx_m * dRm_dp_xm/dRm_dvx_m
                            + dRc_dvy_p * dRm_dp_yp/dRm_dvy_p
                            + dRc_dvy_m * dRm_dp_ym/dRm_dvy_m;
                          

                          dp(center)=
                            delta_R_continuity*theta_continuity/dRc_dp;
                          // dp(center)=
                          //   delta_R_continuity*theta_continuity;



                          // if(ln==2)
                          //   tbox::plog << "smooth p "
                          //              << i << " "
                          //              << j << " "
                          //              << dRc_dp << " "
                          //              // << dp(center) << " "
                          //              // << delta_R_continuity << " "
                          //              // << dvx_dx << " "
                          //              // << dvy_dy << " "
                          //              // << p_rhs(center) << " "
                          //              // << v(pdat::SideIndex(center,pdat::SideIndex::X,
                          //              //                         pdat::SideIndex::Upper)) << " "
                          //              // << v(pdat::SideIndex(center,pdat::SideIndex::X,
                          //              //                         pdat::SideIndex::Lower)) << " "
                          //              // << v(pdat::SideIndex(center,pdat::SideIndex::Y,
                          //              //                         pdat::SideIndex::Upper)) << " "
                          //              // <<  v(pdat::SideIndex(center,pdat::SideIndex::Y,
                          //              //                          pdat::SideIndex::Lower)) << " "

                          //              << "\n";

                          p(center)+=dp(center);
                        }
                    }
                }
            }
        }


      /* fix v sweep */
      for(int rb=0;rb<2;++rb)
        {
          // Need to sync
          xeqScheduleGhostFillNoCoarse(dp_id,invalid_id,ln);
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
              tbox::Pointer<pdat::SideData<double> > v_rhs_ptr =
                patch->getPatchData(v_rhs_id);
              pdat::SideData<double> &v_rhs(*v_rhs_ptr);
                
              tbox::Pointer<pdat::CellData<double> > cell_visc_ptr
                = patch->getPatchData(cell_viscosity_id);
              pdat::CellData<double> &cell_viscosity(*cell_visc_ptr);
              tbox::Pointer<pdat::NodeData<double> > edge_visc_ptr
                = patch->getPatchData(edge_viscosity_id);
              pdat::NodeData<double> &edge_viscosity(*edge_visc_ptr);

              hier::Box pbox=patch->getBox();
              tbox::Pointer<geom::CartesianPatchGeometry>
                geom = patch->getPatchGeometry();
              double dx = *(geom->getDx());
              double dy = *(geom->getDx());

              /* Set an array of bools that tells me whether a point
                 should set the pressure or just let it be.  This is
                 needed at coarse/fine boundaries where the pressure
                 is fixed. */
              hier::Box gbox=p.getGhostBox();
              std::vector<bool> set_p(gbox.size(),true);

              const tbox::Array<hier::BoundaryBox >&edges
                =d_cf_boundary[ln]->getEdgeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<edges.size(); ++mm)
                for(int j=edges[mm].getBox().lower(1);
                    j<=edges[mm].getBox().upper(1); ++j)
                  for(int i=edges[mm].getBox().lower(0);
                      i<=edges[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              const tbox::Array<hier::BoundaryBox >&nodes
                =d_cf_boundary[ln]->getNodeBoundaries(patch->getGlobalId());
              for(int mm=0; mm<nodes.size(); ++mm)
                for(int j=nodes[mm].getBox().lower(1);
                    j<=nodes[mm].getBox().upper(1); ++j)
                  for(int i=nodes[mm].getBox().lower(0);
                      i<=nodes[mm].getBox().upper(0); ++i)
                    set_p[(i-gbox.lower(0))
                          + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(0,0))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;
                  
              if(geom->getTouchesRegularBoundary(0,1))
                for(int j=gbox.lower(1); j<=gbox.upper(1); ++j)
                  set_p[(gbox.upper(0)-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]=false;

              if(geom->getTouchesRegularBoundary(1,0))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[i-gbox.lower(0)]=false;

              if(geom->getTouchesRegularBoundary(1,1))
                for(int i=gbox.lower(0); i<=gbox.upper(0); ++i)
                  set_p[(i-gbox.lower(0))
                        + (gbox.upper(0)-gbox.lower(0)+1)*(gbox.upper(1)-gbox.lower(1))]=
                    false;

              for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
                {
                  /* Do the red-black skip */
                  int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
                  for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                    {
                      pdat::CellIndex center(tbox::Dimension(2));
                      center[0]=i;
                      center[1]=j;

                      pdat::CellIndex up(center), down(center), right(center),
                        left(center);

                      ++up[1];
                      --down[1];
                      ++right[0];
                      --left[0];

                      const pdat::SideIndex
                        center_x(center,0,pdat::SideIndex::Lower),
                        left_x(left,0,pdat::SideIndex::Lower),
                        right_x(right,0,pdat::SideIndex::Lower),
                        center_y(center,1,pdat::SideIndex::Lower),
                        up_y(up,1,pdat::SideIndex::Lower),
                        down_y(down,1,pdat::SideIndex::Lower);
                      const pdat::NodeIndex
                        center_e(center,pdat::NodeIndex::LowerLeft),
                        up_e(up,pdat::NodeIndex::LowerLeft),
                        right_e(right,pdat::NodeIndex::LowerLeft);

                      /* Update v */
                      if(set_p[(i-gbox.lower(0))
                               + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]
                         || (i==pbox.upper(0)+1 && j<pbox.upper(1)+1))
                        {
                          if(!((center[0]==pbox.lower(0)
                                && v(left_x)==boundary_value)
                               || (center[0]==pbox.upper(0)+1
                                   && v(right_x)==boundary_value)))

                            // v(pdat::SideIndex(center,0,pdat::SideIndex::Lower))
                            //   -=(dp(center) - dp(left))
                            //   /(dx*3*(1/(dx*dx) + 1/(dy*dy)));

                          // tbox::plog << "dRm_dv "
                          //            << i << " "
                          //            << j << " "
                          //            << -3*(1/(dx*dx) + 1/(dy*dy)) << " "
                          //            << dRm_dv(cell_viscosity,edge_viscosity,center,
                          //                      left,up_e,center_e,dx,dy) << " "
                          //            << (dp(center) - dp(left))
                          //   /(dx*2*(1/(dx*dx) + 1/(dy*dy))) << " "
                          //            << (dp(center) - dp(left))
                          //   /(dx*dRm_dv(cell_viscosity,edge_viscosity,center,
                          //               left,up_e,center_e,dx,dy)) << " "
                          //            << "\n";

                            v(center_x)+=(dp(center) - dp(left))
                              /(dx*dRm_dv(cell_viscosity,edge_viscosity,center,
                                          left,up_e,center_e,dx,dy));
                        }
                      if(set_p[(i-gbox.lower(0))
                               + (gbox.upper(0)-gbox.lower(0)+1)*(j-gbox.lower(1))]
                         || (i<pbox.upper(0)+1 && j==pbox.upper(1)+1))
                        {
                          if(!((center[1]==pbox.lower(1)
                                && v(down_y)==boundary_value)
                               || (center[1]==pbox.upper(1)+1
                                   && v(up_y)==boundary_value)))

                            // v(pdat::SideIndex(center,1,pdat::SideIndex::Lower))
                            //   -=(dp(center) - dp(down))
                            //   /(dy*3*(1/(dx*dx) + 1/(dy*dy)));
                            v(center_y)+=(dp(center) - dp(down))
                              /(dy*dRm_dv(cell_viscosity,edge_viscosity,center,
                                          down,right_e,center_e,dy,dx));
                        }
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
            << "Smooth " << ln << " " << sweep << " : " << maxres << "\n";
      // }
    }
}

