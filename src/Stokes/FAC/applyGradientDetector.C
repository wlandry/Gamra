#include "Stokes/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

void Stokes::FAC::applyGradientDetector
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_,
 const int ln,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy__ = hierarchy_;
  SAMRAI::hier::PatchHierarchy& hierarchy = *hierarchy__;
  SAMRAI::hier::PatchLevel& level =
    (SAMRAI::hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
  
  int ntag = 0, ntotal = 0;
  double maxestimate = 0;
  for(SAMRAI::hier::PatchLevel::Iterator pi(level.begin()); pi!=level.end(); ++pi)
    {
      SAMRAI::hier::Patch& patch = **pi;
      boost::shared_ptr<SAMRAI::hier::PatchData>
        tag_data = patch.getPatchData(tag_index);
      ntotal += patch.getBox().numberCells().getProduct();
      if (!tag_data)
        {
          TBOX_ERROR("Data index "
                     << tag_index << " does not exist for patch.\n");
        }
      boost::shared_ptr<SAMRAI::pdat::CellData<int> > tag_cell_data_ = 
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<int> >(tag_data);
      if (!tag_cell_data_)
        {
          TBOX_ERROR("Data index " << tag_index << " is not cell int data.\n");
        }
      boost::shared_ptr<SAMRAI::hier::PatchData> soln_data = patch.getPatchData(p_id);
      if (!soln_data)
        {
          TBOX_ERROR("Data index " << p_id << " does not exist for patch.\n");
        }
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > soln_cell_data_ = 
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >(soln_data);
      if (!soln_cell_data_)
        {
          TBOX_ERROR("Data index " << p_id << " is not cell int data.\n");
        }
      SAMRAI::pdat::CellData<double>& soln_cell_data = *soln_cell_data_;
      SAMRAI::pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
      SAMRAI::pdat::CellData<double> estimate_data
        (patch.getBox(),1,SAMRAI::hier::IntVector::getZero(d_dim));
      computeAdaptionEstimate(estimate_data,soln_cell_data);
                              
      tag_cell_data.fill(0);
      SAMRAI::pdat::CellIterator
        cend(SAMRAI::pdat::CellGeometry::end(patch.getBox()));
      for (SAMRAI::pdat::CellIterator
             ci(SAMRAI::pdat::CellGeometry::begin(patch.getBox()));
           ci!=cend; ++ci)
        {
          const SAMRAI::pdat::CellIndex cell_index(*ci);
          if (maxestimate < estimate_data(cell_index))
            maxestimate=estimate_data(cell_index);

          // SAMRAI::tbox::plog << "estimate "
          //            << cell_index << " "
          //            << d_adaption_threshold << " "
          //            << estimate_data(cell_index) << " "
          //            << std::boolalpha
          //            << (estimate_data(cell_index) > d_adaption_threshold)
          //            << " "
          //            << "\n";
          if (estimate_data(cell_index) > d_adaption_threshold
              || ln<min_full_refinement_level)
            {
              tag_cell_data(cell_index) = 1;
              ++ntag;
            }
        }
    }
  SAMRAI::tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  SAMRAI::tbox::plog << "Number of cells tagged on level " << ln << " is "
             << ntag << "/" << ntotal << "\n";
  SAMRAI::tbox::plog << "Max estimate is " << maxestimate << "\n";
}


void Stokes::FAC::computeAdaptionEstimate
(SAMRAI::pdat::CellData<double>& estimate_data,
 const SAMRAI::pdat::CellData<double>& soln_cell_data) const
{
  TBOX_ERROR("Need to implement computeAdaptionEstimate without MDA_AccessConst");

  // const int* lower = &estimate_data.getBox().lower()[0];
  // const int* upper = &estimate_data.getBox().upper()[0];
  // if (d_dim == SAMRAI::tbox::Dimension(2))
  //   {
  //     MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > co =
  //       SAMRAI::pdat::ArrayDataAccess::access<2, double>(soln_cell_data.getArrayData());
  //     MDA_Access<double, 2, MDA_OrderColMajor<2> > es =
  //       SAMRAI::pdat::ArrayDataAccess::access<2, double>(estimate_data.getArrayData());
  //     int i, j;
  //     double estimate, est0, est1, est2, est3, est4, est5;
  //     for (j = lower[1]; j <= upper[1]; ++j)
  //       {
  //         for (i = lower[0]; i <= upper[0]; ++i)
  //           {
  //             est0=SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j) + co(i-1,j)
  //                                                   - 2*co(i,j));
                                                 
  //             est1=SAMRAI::tbox::MathUtilities<double>::Abs(co(i,j+1) + co(i,j-1)
  //                                                   - 2*co(i,j));
  //             est2=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j+1)
  //                                                         + co(i-1,j-1)
  //                                                         - 2*co(i,j));
  //             est3=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j-1)
  //                                                         + co(i-1,j+1)
  //                                                         - 2*co(i,j));
  //             est4=SAMRAI::tbox::MathUtilities<double>::Max(est0,est1);
  //             est5=SAMRAI::tbox::MathUtilities<double>::Max(est2,est3);
  //             estimate=SAMRAI::tbox::MathUtilities<double>::Max(est4,est5);
  //             es(i,j)=estimate;
  //           }
  //       }
  //   }
  // else if (d_dim == SAMRAI::tbox::Dimension(3))
  //   {
  //     MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > co =
  //       SAMRAI::pdat::ArrayDataAccess::access<3, double>(soln_cell_data.getArrayData());
  //     MDA_Access<double, 3, MDA_OrderColMajor<3> > es =
  //       SAMRAI::pdat::ArrayDataAccess::access<3, double>(estimate_data.getArrayData());
  //     int i, j, k;
  //     double estimate, est0, est1, est2, est3, est4, est5, est6, est7, est8,
  //       esta, estb, estc, estd, este, estf, estg;
  //     for (k = lower[2]; k <= upper[2]; ++k)
  //       for (j = lower[1]; j <= upper[1]; ++j)
  //         for (i = lower[0]; i <= upper[0]; ++i)
  //           {
  //             est0=SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j,k) + co(i-1,j,k)
  //                                                   - 2*co(i,j,k));
  //             est1=SAMRAI::tbox::MathUtilities<double>::Abs(co(i,j+1,k) + co(i,j-1,k)
  //                                                   - 2*co(i,j,k));
  //             est2=SAMRAI::tbox::MathUtilities<double>::Abs(co(i,j,k+1) + co(i,j,k-1)
  //                                                   - 2*co(i,j,k));
  //             est3=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i,j+1,k+1)
  //                                                         + co(i,j-1,k-1)
  //                                                         - 2*co(i,j,k));
  //             est4=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i,j+1,k-1)
  //                                                         + co(i,j-1,k+1)
  //                                                         - 2*co(i,j,k));
  //             est5=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j,k+1)
  //                                                         + co(i-1,j,k-1)
  //                                                         - 2*co(i,j,k));
  //             est6=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j,k-1)
  //                                                         + co(i-1,j,k+1)
  //                                                         - 2*co(i,j,k));
  //             est7=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j+1,k)
  //                                                         + co(i-1,j-1,k)
  //                                                         - 2*co(i,j,k));
  //             est8=0.5 * SAMRAI::tbox::MathUtilities<double>::Abs(co(i+1,j-1,k)
  //                                                         + co(i-1,j+1,k)
  //                                                         - 2*co(i,j,k));
  //             esta=SAMRAI::tbox::MathUtilities<double>::Max(est0,est1);
  //             estb=SAMRAI::tbox::MathUtilities<double>::Max(est2,est3);
  //             estc=SAMRAI::tbox::MathUtilities<double>::Max(est4,est5);
  //             estd=SAMRAI::tbox::MathUtilities<double>::Max(est6,est7);
  //             este=SAMRAI::tbox::MathUtilities<double>::Max(esta,estb);
  //             estf=SAMRAI::tbox::MathUtilities<double>::Max(estc,estd);
  //             estg=SAMRAI::tbox::MathUtilities<double>::Max(este,estf);
  //             estimate=SAMRAI::tbox::MathUtilities<double>::Max(estg,est8);
  //             es(i,j,k)=estimate;
  //           }
  //   }
}
