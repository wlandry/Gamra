/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"
#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"
#include "Constants.hxx"

void Elastic::Operators::Gauss_Seidel_red_black_3D
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

  double theta_momentum=1.0;
  double maxres;
  const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const SAMRAI::hier::Index pp[]={ip,jp,kp};
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps && !converged; ++sweep)
    {
      maxres=0;
      for(Gamra::Dir ix=0;ix<3;++ix)
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
                  cell_moduli_ptr =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
                  (patch->getPatchData(cell_moduli_id));
                SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
                boost::shared_ptr<SAMRAI::pdat::EdgeData<double> >
                  edge_moduli_ptr =
                  boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
                  (patch->getPatchData(edge_moduli_id));
                SAMRAI::pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);

                SAMRAI::hier::Box pbox=patch->getBox();
                boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
                  boost::dynamic_pointer_cast
                  <SAMRAI::geom::CartesianPatchGeometry>
                  (patch->getPatchGeometry());
                const double *Dx = geom->getDx();

                for(int k=pbox.lower(2); k<=pbox.upper(2)+pp[ix][2]; ++k)
                  for(int j=pbox.lower(1); j<=pbox.upper(1)+pp[ix][1]; ++j)
                    {
                      /// Do the red-black skip
                      int i_min=pbox.lower(0)
                        + (abs(pbox.lower(0) + j + k + rb))%2;
                      for(int i=i_min; i<=pbox.upper(0)+pp[ix][0]; i+=2)
                        {
                          SAMRAI::pdat::CellIndex
                            center(SAMRAI::hier::Index(i,j,k));
                          update_V_3D(ix,pbox,v,v_rhs,cell_moduli,
                                      edge_moduli,center,
                                      Dx,theta_momentum,pp,maxres);
                        }
                    }
              }
          }

      if (residual_tolerance >= 0.0)
        {
          converged = maxres < residual_tolerance;
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

