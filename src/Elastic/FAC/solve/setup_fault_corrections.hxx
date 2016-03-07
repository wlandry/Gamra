#pragma once

/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "setup_fault_corrections/correct_rhs.hxx"

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

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
  const int max_level(hierarchy->getFinestLevelNumber());
  const Gamra::Dir dim=dimension.getValue();

  for(int l=0; l<=max_level; ++l)
    {
      SAMRAI::hier::PatchLevel &patch_level(*hierarchy->getPatchLevel(l));
      for(SAMRAI::hier::PatchLevel::Iterator p(patch_level.begin());
          p!=patch_level.end(); ++p)
        {
          const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v_rhs
            (boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
             ((*p)->getPatchData(v_rhs_id)));

          SAMRAI::hier::Box pbox = v_rhs->getBox();
          SAMRAI::hier::Box gbox = v_rhs->getGhostBox();
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

          const boost::shared_ptr<SAMRAI::pdat::CellData<double> >
            &dv_diagonal_ptr
            (boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
             ((*p)->getPatchData(dv_diagonal_id)));
          const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_ptr
            (boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
             ((*p)->getPatchData(dv_mixed_id)));

          compute_dv(faults,dim,dx,geom->getXLower(),gbox,pbox.lower(),
                     *dv_diagonal_ptr,*dv_mixed_ptr);

          const boost::shared_ptr<SAMRAI::pdat::CellData<double> >
            &cell_moduli_ptr
            (boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
             ((*p)->getPatchData(cell_moduli_id)));

          const boost::shared_ptr<T> &edge_moduli_ptr
            (boost::dynamic_pointer_cast<T>((*p)->getPatchData(edge_moduli_id)));

          correct_rhs(dimension,dim,dx,pbox,*cell_moduli_ptr,*edge_moduli_ptr,
                      *dv_diagonal_ptr,*dv_mixed_ptr,*v_rhs);
        }
    }
}


