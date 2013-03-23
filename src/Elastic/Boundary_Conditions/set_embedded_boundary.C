#include "Constants.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_embedded_boundary
(const SAMRAI::hier::Patch& patch, const bool &homogeneous)
{
  if(!have_faults() || homogeneous)
    return;

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(level_set_id));
  SAMRAI::pdat::SideData<double> &level_set(*level_set_ptr);

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed_ptr=
    boost::dynamic_pointer_cast
    <SAMRAI::pdat::SideData<double> >(patch.getPatchData(dv_mixed_id));
  SAMRAI::pdat::SideData<double> &dv_mixed(*dv_mixed_ptr);
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal_ptr=
    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
    (patch.getPatchData(dv_diagonal_id));
  SAMRAI::pdat::CellData<double> &dv_diagonal(*dv_diagonal_ptr);

  const SAMRAI::tbox::Dimension Dim(patch.getDim());
  const int dim(Dim.getValue());

  const SAMRAI::hier::Box gbox=level_set.getGhostBox();

  /* FIXME: Why is this even required?  Shouldn't the boundaries
     be set once and then forgotten?  On coarse levels, the
     boundaries may not be copied over. Taking this part out
     certainly breaks the code. */

  /* FIXME: This looping seems really excessive.  It seems like
     there should be a better way using getBoundaryBoxes. */

  for(int ix=0; ix<dim; ++ix)
    {
      SAMRAI::pdat::SideIterator s_end(gbox,ix,false);
      for(SAMRAI::pdat::SideIterator si(gbox,ix,true); si!=s_end; si++)
        {
          SAMRAI::pdat::SideIndex s(*si);

          if(level_set(s)<0)
            for(int d=0;d<(dim==2 ? 2 : 8);++d)
              dv_mixed(s,d)=0;
        }
    }

  SAMRAI::pdat::CellIterator c_end(gbox,false);
  for(SAMRAI::pdat::CellIterator ci(gbox,true); ci!=c_end; ci++)
    {
      SAMRAI::pdat::CellIndex c(*ci);
      for(int ix=0;ix<dim;++ix)
        {
          SAMRAI::pdat::SideIndex x_m(c,ix,SAMRAI::pdat::SideIndex::Lower),
            x_p(c,ix,SAMRAI::pdat::SideIndex::Upper);
          if(level_set(x_m) + level_set(x_p)<0)
            dv_diagonal(c,ix)=0;
        }
    }
}
