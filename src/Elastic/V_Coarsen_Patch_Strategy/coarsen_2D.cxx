/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/V_Coarsen_Patch_Strategy.hxx"

namespace
{
  void
  coarsen_point_2D(const SAMRAI::pdat::SideIndex& coarse,
                   const SAMRAI::pdat::SideIndex& fine,
                   const SAMRAI::hier::Index& ip, const SAMRAI::hier::Index& jp,
                   SAMRAI::pdat::SideData<double>& v,
                   const SAMRAI::pdat::SideData<double>& v_fine )
  {
    v(coarse)=(v_fine(fine) + v_fine(fine+jp))/4
      + (v_fine(fine-ip) + v_fine(fine-ip+jp)
         + v_fine(fine+ip) + v_fine(fine+jp+ip))/8;
  }

  double
  coarsen_correction_2D(const SAMRAI::pdat::SideIndex& fine,
                        const int& axis,
                        const SAMRAI::hier::Index& ip,
                        const SAMRAI::hier::Index& jp,
                        const SAMRAI::pdat::CellData<double>& dv_diagonal,
                        const SAMRAI::pdat::SideData<double>& dv_mixed)
  {
    SAMRAI::pdat::CellIndex cell(fine);
    return (dv_diagonal(cell-ip,axis) - dv_diagonal(cell,axis)
            + dv_diagonal(cell-ip+jp,axis) - dv_diagonal(cell+jp,axis))/8
      + (dv_mixed(fine,0) + dv_mixed(fine+jp,1))/2;
  }
}

void Elastic::V_Coarsen_Patch_Strategy::coarsen_2D
(SAMRAI::pdat::SideData<double>& v,
 const SAMRAI::pdat::SideData<double>& v_fine,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
 const boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal,
 SAMRAI::hier::Patch &coarse_patch,
 const SAMRAI::hier::Patch &fine_patch,
 const SAMRAI::hier::Box& coarse_box) const
{
  /// Numbering of v nodes is

  ///    x--x--x--x--x  Fine
  ///    0  1  2  3  4

  ///    x-----x-----x Coarse
  ///    0     1     2

  ///    So the i'th coarse point is affected by the i*2-1,
  ///    i*2, and i*2+1 fine points.  So, for example, i_fine=3
  ///    affects both i_coarse=1 and i_coarse=2.

  ///    |---------------------------------------------------------------|
  ///    |               |               |               |               |
  ///    f       f       f       f       f       f       f       f       f
  ///    |               |               |               |               |
  ///    c               c               c               c               c
  ///    |               |               |               |               |
  ///    f       f       f       f       f       f       f       f       f
  ///    |               |               |               |               |
  ///    |---------------------------------------------------------------|
  ///    |               |               |               |               |
  ///    f       f       f       f       f       f       f       f       f
  ///    |               |               |               |               |
  ///    c               c               c               c               c
  ///    |               |               |               |               |
  ///    f       f       f       f       f       f       f       f       f
  ///    |               |               |               |               |
  ///    |---------------------------------------------------------------|

  ///    In 2D, a coarse point depends on six points.  In this
  ///    case, (i*2,j*2), (i*2,j*2+1), (i*2-1,j*2),
  ///    (i*2-1,j*2+1), (i*2+1,j*2), (i*2+1,j*2+1).

  ///    The coarse/fine boundaries get fixed up later in
  ///    V_Coarsen_Patch_Strategy::postprocessCoarsen.

  SAMRAI::hier::Index ip(1,0), jp(0,1);
  const SAMRAI::hier::Index unit[]={ip,jp};

  /// From reading CoarsenSchedule::coarsenSourceData in
  /// SAMRAI/source/SAMRAI/xfer/CoarsenSchedule.C, it seems that the
  /// coarse box is created from the fine box.  So the coarse box is
  /// always covered by the fine box, meaning we do not have to do an
  /// intersection.

  SAMRAI::hier::Box big_box(coarse_box);
  big_box.growUpper(SAMRAI::hier::IntVector::getOne(fine_patch.getDim()));
  const SAMRAI::pdat::CellIterator end(SAMRAI::pdat::CellGeometry::end(big_box));
  const int dim(2);

  if(have_embedded_boundary())
    {
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_coarse_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (coarse_patch.getPatchData(level_set_id));
      SAMRAI::pdat::SideData<double> &level_set_coarse(*level_set_coarse_ptr);

      boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_fine_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (fine_patch.getPatchData(level_set_id));
      SAMRAI::pdat::SideData<double> &level_set_fine(*level_set_fine_ptr);

      for(SAMRAI::pdat::CellIterator
            ci(SAMRAI::pdat::CellGeometry::begin(big_box)); ci!=end; ++ci)
        {
          const SAMRAI::pdat::CellIndex &coarse_cell(*ci);
          for(Gamra::Dir ix=0;ix<dim;++ix)
            {
              const Gamra::Dir iy(ix.next(dim));
              if(coarse_cell[iy]!=big_box.upper(iy))
                {
                  SAMRAI::pdat::SideIndex
                    coarse(coarse_cell,ix,SAMRAI::pdat::SideIndex::Lower);
                  SAMRAI::pdat::SideIndex fine(coarse*2);
                  SAMRAI::pdat::CellIndex fine_cell(fine);

                  if(level_set_coarse(coarse)<0)
                    {
                      // When restricting the residual, we can not set
                      //    points outside to invalid numbers.  The norm
                      //    is computed with SAMRAIVectorReal::maxNorm,
                      //    which does not know about level sets and
                      //    looks at all points on all levels.
                      if(!is_residual)
                        { v(coarse)=invalid_value; }
                      else
                        { v(coarse)=0; }
                    }
                  else
                    {
                      double plus(0), minus(0);
                      if(level_set_fine(fine)<0)
                        {
                          // FIXME: need to correctly interpolate
                          minus=0;
                          if(have_faults() && !is_residual)
                            { minus+=(*dv_mixed)(fine,0); }
                        }
                      else
                        {
                          if(level_set_fine(fine+unit[ix])<0
                             || level_set_fine(fine-unit[ix])<0)
                            {
                              minus=v_fine(fine);
                              if(have_faults() && !is_residual)
                                { minus+=(*dv_mixed)(fine,0); }
                            }
                          else
                            {
                              minus=v_fine(fine)/2
                                + (v_fine(fine+unit[ix])
                                   + v_fine(fine-unit[ix]))/4;
                              if(have_faults() && !is_residual)
                                { minus+=(*dv_mixed)(fine,0)
                                    + ((*dv_diagonal)(fine_cell-unit[ix],ix)
                                       - (*dv_diagonal)(fine_cell,ix))/4; }
                            }
                        }
                      if(level_set_fine(fine+unit[iy])<0)
                        {
                          // FIXME: need to correctly interpolate
                          plus=0;
                          if(have_faults() && !is_residual)
                            { plus+=(*dv_mixed)(fine+unit[iy],1); }
                        }
                      else
                        {
                          if(level_set_fine(fine+unit[iy]+unit[ix])<0
                             || level_set_fine(fine+unit[iy]-unit[ix])<0)
                            {
                              plus=v_fine(fine+unit[iy]);
                              if(have_faults() && !is_residual)
                                { plus+=(*dv_mixed)(fine+unit[iy],1); }
                            }
                          else
                            {
                              plus=v_fine(fine+unit[iy])/2
                                + (v_fine(fine+unit[iy]+unit[ix])
                                   + v_fine(fine+unit[iy]-unit[ix]))/4;
                              if(have_faults() && !is_residual)
                                { plus+=(*dv_mixed)(fine+unit[iy],1)
                                    + ((*dv_diagonal)
                                       (fine_cell+unit[iy]-unit[ix],ix)
                                       - (*dv_diagonal)
                                       (fine_cell+unit[iy],ix))/4; }
                            }
                        }
                      v(coarse)=(plus+minus)/2;
                    }
                }
            }
        }
    }
  else
    {
      const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
        boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
        (coarse_patch.getPatchGeometry());

      for(SAMRAI::pdat::CellIterator
            ci(SAMRAI::pdat::CellGeometry::begin(big_box)); ci!=end; ++ci)
        {
          const SAMRAI::pdat::CellIndex &cell(*ci);
          for(Gamra::Dir ix=0;ix<2;++ix)
            {
              const Gamra::Dir iy(ix.next(dim));
              if(cell[iy]!=coarse_box.upper(iy)+1)
                {
                  SAMRAI::pdat::SideIndex
                    coarse(cell,ix,SAMRAI::pdat::SideIndex::Lower);
                  SAMRAI::pdat::SideIndex fine(coarse*2);
                  if((cell[ix]==coarse_box.lower(ix)
                      && coarse_geom->getTouchesRegularBoundary(ix,0))
                     || (cell[ix]==coarse_box.upper(ix)+1
                         && coarse_geom->getTouchesRegularBoundary(ix,1)))
                    {
                      v(coarse)=
                        (v_fine(fine) + v_fine(fine+unit[iy]))/2;
                      if(have_faults() && !is_residual)
                        { v(coarse)+=((*dv_mixed)(fine,0)
                                      + (*dv_mixed)(fine+unit[iy],1))/2; }
                    }
                  else
                    {
                      coarsen_point_2D(coarse,fine,unit[ix],unit[iy],
                                       v,v_fine);
                      if(have_faults() && !is_residual)
                        {
                          v(coarse)+=
                            coarsen_correction_2D(fine,ix,unit[ix],
                                                  unit[iy],*dv_diagonal,
                                                  *dv_mixed);
                        }
                    }
                }
            }
        }
    }
}
