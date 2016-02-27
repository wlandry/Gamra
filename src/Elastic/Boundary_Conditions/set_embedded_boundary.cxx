/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Boundary_Conditions.hxx"

namespace {
  void set_embedded_values
  (const SAMRAI::pdat::SideIndex &x,
   const int &dim,
   SAMRAI::pdat::SideData<double> &level_set,
   boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr)
  {
    level_set(x)=-boundary_value;
    if(dv_mixed_ptr)
      for(int d=0;d<(dim==2 ? 2 : 8);++d)
        (*dv_mixed_ptr)(x,d)=0;
  }
}

void Elastic::Boundary_Conditions::set_embedded_boundary
(const SAMRAI::hier::Patch& patch, const int &level_set_id,
 const int &dv_mixed_id)
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(level_set_id));
  SAMRAI::pdat::SideData<double> &level_set(*level_set_ptr);

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr;
  if(dv_mixed_id!=invalid_id)
    dv_mixed_ptr=boost::dynamic_pointer_cast
      <SAMRAI::pdat::SideData<double> >(patch.getPatchData(dv_mixed_id));

  const SAMRAI::tbox::Dimension dimension(patch.getDim());
  const Gamra::Dir dim(dimension.getValue());

  const SAMRAI::hier::Box gbox=level_set.getGhostBox();
  const SAMRAI::hier::Box box=level_set.getBox();

  boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (patch.getPatchGeometry());

  /// Corners are not synced, so we set the level set to a generic
  /// negative value and dv_mixed to zero.  Setting dv_mixed to zero
  /// implies that we do not care about faults that intersect that
  /// corner (nor should we).
  for(Gamra::Dir ix=0; ix<dim; ++ix)
    {
      const Gamra::Dir iy(ix.next(dim));
      const Gamra::Dir iz(iy.next(dim));

      std::vector<std::vector<int> > corners(dim);
      corners[ix].push_back(box.lower(ix));
      corners[ix].push_back(box.upper(ix)+1);

      if(geom->getTouchesRegularBoundary(iy,0))
        corners[iy].push_back(gbox.lower(iy));
      if(geom->getTouchesRegularBoundary(iy,1))
        corners[iy].push_back(gbox.upper(iy));

      if(dim==3)
        {
          if(geom->getTouchesRegularBoundary(iz,0))
            corners[iz].push_back(gbox.lower(iz));
          if(geom->getTouchesRegularBoundary(iz,1))
            corners[iz].push_back(gbox.upper(iz));
        }

      SAMRAI::pdat::SideIndex x(SAMRAI::hier::Index(dimension),ix,
                                SAMRAI::pdat::SideIndex::Lower);
      for(std::vector<int>::iterator i=corners[ix].begin();
          i!=corners[ix].end(); ++i)
        {
          x[ix]=*i;
          if(dim==2)
            {
              for(std::vector<int>::iterator j=corners[iy].begin();
                  j!=corners[iy].end();++j)
                {
                  x[iy]=*j;
                  set_embedded_values(x,dim,level_set,dv_mixed_ptr);
                }
            }
          else
            {
              for(std::vector<int>::iterator j=corners[iy].begin();
                  j!=corners[iy].end();++j)
                {
                  for(int kk=gbox.lower(iz); kk<=gbox.upper(iz); ++kk)
                    {
                      x[iy]=*j;
                      x[iz]=kk;
                      set_embedded_values(x,dim,level_set,dv_mixed_ptr);
                    }
                }
              for(std::vector<int>::iterator k=corners[iz].begin();
                  k!=corners[iz].end();++k)
                {
                  for(int jj=gbox.lower(iy); jj<=gbox.upper(iy); ++jj)
                    {
                      x[iy]=jj;
                      x[iz]=*k;
                      set_embedded_values(x,dim,level_set,dv_mixed_ptr);
                    }
                }
            }
        }
    }
}
