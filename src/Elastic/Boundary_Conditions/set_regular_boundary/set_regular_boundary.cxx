/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Boundary_Conditions.hxx"

void Elastic::Boundary_Conditions::set_regular_boundary
(const SAMRAI::hier::Patch& patch, const int &v_id, const bool &homogeneous,
 const bool &apply_normal_stress) const
{
  /// We need to set dirichlet boundaries before traction boundaries,
  /// and tangential traction boudaries before normal traction
  /// boundaries.

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(v_id));
  SAMRAI::pdat::SideData<double> &v(*v_ptr);

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr;
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr;
  if(have_faults() && !homogeneous)
    {
      dv_mixed_ptr=boost::dynamic_pointer_cast
        <SAMRAI::pdat::SideData<double> >(patch.getPatchData(dv_mixed_id));
      dv_diagonal_ptr=
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch.getPatchData(dv_diagonal_id));
    }

  const SAMRAI::tbox::Dimension Dim(patch.getDim());
  const Gamra::Dir dim(Dim.getValue());
  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(Dim));

  SAMRAI::hier::Index unit[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    unit[i][i]=1;

  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (patch.getPatchGeometry());
  const double *dx=geom->getDx();

  const SAMRAI::hier::Box pbox=patch.getBox();
  const SAMRAI::hier::Box gbox=v.getGhostBox();

  /// FIXME: Why is this even required?  Shouldn't the boundaries be
  /// set once and then forgotten?  On coarse levels, the boundaries
  /// may not be copied over. Taking this part out certainly breaks
  /// the code.

  /// FIXME: This looping seems really excessive.  It seems like there
  /// should be a better way using getBoundaryBoxes.
  if(have_faults() && !homogeneous)
    {
      SAMRAI::pdat::SideData<double> &dv_mixed(*dv_mixed_ptr);
      for(Gamra::Dir ix=0; ix<dim; ++ix)
        {
          SAMRAI::pdat::SideIterator
            s_end(SAMRAI::pdat::SideGeometry::end(gbox,ix));
          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(gbox,ix)); si!=s_end; ++si)
            {
              const SAMRAI::pdat::SideIndex &s(*si);
              if((s[ix]<pbox.lower(ix)
                  && geom->getTouchesRegularBoundary(ix,0))
                 || (s[ix]>pbox.upper(ix)+1
                     && geom->getTouchesRegularBoundary(ix,1)))
                {
                  for(int d=0;d<(dim==2 ? 2 : 8);++d)
                    { dv_mixed(s,d)=0; }
                }
              else
                {
                  for(Gamra::Dir iy=ix.next(dim); iy!=ix;
                      iy=iy.next(dim))
                    {
                      if((s[iy]<pbox.lower(iy)
                          && geom->getTouchesRegularBoundary(iy,0))
                         || (s[iy]>pbox.upper(iy)
                             && geom->getTouchesRegularBoundary(iy,1)))
                        {
                          for(int d=0;d<(dim==2 ? 2 : 8);++d)
                            { dv_mixed(s,d)=0; }
                        }
                    }
                }
            }
        }

      boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr=
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch.getPatchData(dv_diagonal_id));
      SAMRAI::pdat::CellData<double>& dv_diagonal(*dv_diagonal_ptr);

      SAMRAI::pdat::CellIterator c_end(SAMRAI::pdat::CellGeometry::end(gbox));
      for(SAMRAI::pdat::CellIterator
            ci(SAMRAI::pdat::CellGeometry::begin(gbox)); ci!=c_end; ++ci)
        {
          const SAMRAI::pdat::CellIndex &c(*ci);
          bool is_boundary(false);
          for(Gamra::Dir d=0;d<dim;++d)
            { is_boundary= is_boundary
                || ((c[d]<pbox.lower(d)
                     && geom->getTouchesRegularBoundary(d,0))
                    || (c[d]>pbox.upper(d)
                        && geom->getTouchesRegularBoundary(d,1))); }

          if(is_boundary)
            {
              for(Gamra::Dir ix=0;ix<dim;++ix)
                { dv_diagonal(c,ix)=0; }
            }
        }
    }

  set_dirichlet(v,dv_mixed_ptr,unit,dim,pbox,gbox,geom,dx,homogeneous);
  set_shear_stress(v,dv_mixed_ptr,unit,dim,pbox,gbox,geom,dx,homogeneous);

  if(apply_normal_stress)
    {
      if(dim==2)
        {
          boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_moduli_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
            (patch.getPatchData(edge_moduli_id));
          set_normal_stress(v,dv_diagonal_ptr,*edge_moduli_ptr,unit,dim,
                            pbox,gbox,geom,dx,homogeneous);
        }
      else
        {
          boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_moduli_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
            (patch.getPatchData(edge_moduli_id));
          set_normal_stress(v,dv_diagonal_ptr,*edge_moduli_ptr,unit,dim,pbox,
                            gbox,geom,dx,homogeneous);
        }
    }
}
