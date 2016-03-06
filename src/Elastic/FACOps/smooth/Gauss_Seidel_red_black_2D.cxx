/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "../v_level_set_operator_2D.hxx"
#include "../v_operator_2D.hxx"

#include "Elastic/FACOps.hxx"
#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"
#include "Elastic/dRm_dv.hxx"

void Elastic::FACOps::Gauss_Seidel_red_black_2D
(SAMRAI::solv::SAMRAIVectorReal<double>& solution,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int ln,
 int num_sweeps,
 double residual_tolerance)
{
  const int v_id(solution.getComponentDescriptorIndex(0)),
    v_rhs_id(residual.getComponentDescriptorIndex(0));

  const SAMRAI::hier::PatchHierarchy &hierarchy=*residual.getPatchHierarchy();
  boost::shared_ptr<SAMRAI::hier::PatchLevel>
    level = hierarchy.getPatchLevel(ln);

  v_refine_patch_strategy.data_id=v_id;
  v_refine_patch_strategy.is_residual=true;
  Coarse_Fine_Boundary_Refine::is_residual=true;
  ghostfill_nocoarse(v_rhs_id,ln);

  if (ln > level_min)
    { ghostfill(v_id, ln); }

  double theta_momentum=1.0;

  const SAMRAI::hier::Index
    unit[]={SAMRAI::hier::Index(1,0),SAMRAI::hier::Index(0,1)};
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps && !converged; ++sweep)
    {
      double max_residual(0);

      const Gamra::Dir dim(2);
      for(Gamra::Dir ix=0; ix<dim; ++ix)
        {
          const Gamra::Dir iy(ix.next(dim));
          const SAMRAI::hier::Index ip(unit[ix]), jp(unit[iy]);
          for(int rb=0;rb<2;++rb)
            {
              ghostfill_nocoarse(v_id,ln);
              if (ln > level_min)
                { ghostfill(v_id, ln); }
              set_physical_boundaries(v_id,level,true);
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
                
                  boost::shared_ptr<SAMRAI::pdat::CellData<double> >
                    cell_moduli_ptr=
                    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                    (patch->getPatchData(cell_moduli_id));
                  SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
                  boost::shared_ptr<SAMRAI::pdat::NodeData<double> >
                    edge_moduli_ptr=
                    boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                    (patch->getPatchData(edge_moduli_id));
                  SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

                  SAMRAI::hier::Box pbox=patch->getBox();
                  boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                    boost::dynamic_pointer_cast
                    <SAMRAI::geom::CartesianPatchGeometry>
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
                          /// Do the red-black skip
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
                                  update_V_2D(ix,pbox,center,unit[ix],unit[iy],
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
                          /// Do the red-black skip
                          int i_min=pbox.lower(0)
                            + (abs(pbox.lower(0) + j + rb))%2;
                          for(int i=i_min; i<=pbox.upper(0)+unit[ix][0]; i+=2)
                            {
                              SAMRAI::pdat::CellIndex
                                center(SAMRAI::hier::Index(i,j));
                              update_V_2D(ix,pbox,center,unit[ix],unit[iy],
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
          converged = max_residual < residual_tolerance;
          const SAMRAI::tbox::SAMRAI_MPI& mpi(hierarchy.getMPI());
          int tmp= converged ? 1 : 0;
          if (mpi.getSize() > 1)
            { mpi.AllReduce(&tmp, 1, MPI_MIN); }
          converged=(tmp==1);
        }
    }

  ghostfill_nocoarse(v_id,ln);
  if (ln > level_min)
    { ghostfill(v_id, ln); }
  set_physical_boundaries(v_id,level,true);
}

