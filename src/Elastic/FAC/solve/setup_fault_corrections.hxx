#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "setup_fault_corrections/correct_rhs.hxx"

void compute_dv(const std::vector<double> &faults,
                const Gamra::Dir &dim, const double *dx,
                const double *geom_xlower,
                const SAMRAI::hier::Box &gbox,
                const SAMRAI::hier::Index &pbox_lower,
                SAMRAI::pdat::CellData<double> &dv_diagonal,
                SAMRAI::pdat::SideData<double> &dv_mixed);

template<class T>
void Elastic::FAC::setup_fault_corrections()
{
  const int max_level(d_hierarchy->getFinestLevelNumber());
  const Gamra::Dir dim=d_dim.getValue();

  for(int l=0; l<=max_level; ++l)
    {
      boost::shared_ptr<SAMRAI::hier::PatchLevel>
        level(d_hierarchy->getPatchLevel(l));
      
      for(SAMRAI::hier::PatchLevel::Iterator p(level->begin());
          p!=level->end(); ++p)
        {
          const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v_rhs_ptr
            (boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
             ((*p)->getPatchData(v_rhs_id)));
          SAMRAI::pdat::SideData<double> &v_rhs(*v_rhs_ptr);

          SAMRAI::hier::Box pbox = v_rhs.getBox();
          SAMRAI::hier::Box gbox = v_rhs.getGhostBox();
          const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> &geom
            (boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
             ((*p)->getPatchGeometry()));
          for(Gamra::Dir ix=0;ix<dim;++ix)
            {
              if(geom->getTouchesRegularBoundary(ix,0))
                { gbox.shorten(ix,-1); }
              if(geom->getTouchesRegularBoundary(ix,1))
                { gbox.shorten(ix,1); }
            }
          const double *dx=geom->getDx();

          boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            ((*p)->getPatchData(dv_diagonal_id));
          SAMRAI::pdat::CellData<double> &dv_diagonal(*dv_diagonal_ptr);
          boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
            ((*p)->getPatchData(dv_mixed_id));
          SAMRAI::pdat::SideData<double> &dv_mixed(*dv_mixed_ptr);

          compute_dv(faults,dim,dx,geom->getXLower(),gbox,
                     pbox.lower(),*dv_diagonal_ptr,*dv_mixed_ptr);

          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            ((*p)->getPatchData(cell_moduli_id));

          boost::shared_ptr<T> edge_moduli_ptr =
            boost::dynamic_pointer_cast<T>((*p)->getPatchData(edge_moduli_id));

          correct_rhs(d_dim,dim,dx,pbox,*cell_moduli_ptr,*edge_moduli_ptr,
                      *dv_diagonal_ptr,*dv_mixed_ptr,v_rhs);
        }
    }
}


