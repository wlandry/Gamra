#include "Stokes/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"

void SAMRAI::Stokes::FAC::applyGradientDetector
(const tbox::Pointer<hier::BasePatchHierarchy> hierarchy_,
 const int ln,
 const double ,
 const int tag_index,
 const bool ,
 const bool )
{
  const tbox::Pointer<hier::PatchHierarchy> hierarchy__ = hierarchy_;
  hier::PatchHierarchy& hierarchy = *hierarchy__;
  hier::PatchLevel& level =
    (hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
  
  int ntag = 0, ntotal = 0;
  double maxestimate = 0;
  for(hier::PatchLevel::Iterator pi(level); pi; pi++)
    {
      hier::Patch& patch = **pi;
      tbox::Pointer<hier::PatchData>
        tag_data = patch.getPatchData(tag_index);
      ntotal += patch.getBox().numberCells().getProduct();
      if (tag_data.isNull())
        {
          TBOX_ERROR("Data index "
                     << tag_index << " does not exist for patch.\n");
        }
      tbox::Pointer<pdat::CellData<int> > tag_cell_data_ = tag_data;
      if (tag_cell_data_.isNull())
        {
          TBOX_ERROR("Data index " << tag_index << " is not cell int data.\n");
        }
      tbox::Pointer<hier::PatchData> soln_data = patch.getPatchData(p_id);
      if (soln_data.isNull())
        {
          TBOX_ERROR("Data index " << p_id << " does not exist for patch.\n");
        }
      tbox::Pointer<pdat::CellData<double> > soln_cell_data_ = soln_data;
      if (soln_cell_data_.isNull())
        {
          TBOX_ERROR("Data index " << p_id << " is not cell int data.\n");
        }
      pdat::CellData<double>& soln_cell_data = *soln_cell_data_;
      pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
      pdat::CellData<double> estimate_data(patch.getBox(),1,
                                           hier::IntVector(d_dim, 0));
      computeAdaptionEstimate(estimate_data,soln_cell_data);
                              
      tag_cell_data.fill(0);
      for (pdat::CellIterator ci(patch.getBox()); ci; ci++)
        {
          const pdat::CellIndex cell_index(*ci);
          if (maxestimate < estimate_data(cell_index))
            maxestimate=estimate_data(cell_index);

          // tbox::plog << "estimate "
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
  tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  tbox::plog << "Number of cells tagged on level " << ln << " is "
             << ntag << "/" << ntotal << "\n";
  tbox::plog << "Max estimate is " << maxestimate << "\n";
}


void SAMRAI::Stokes::FAC::computeAdaptionEstimate
(pdat::CellData<double>& estimate_data,
 const pdat::CellData<double>& soln_cell_data) const
{
  const int* lower = &estimate_data.getBox().lower()[0];
  const int* upper = &estimate_data.getBox().upper()[0];
  if (d_dim == tbox::Dimension(2))
    {
      MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > co =
        pdat::ArrayDataAccess::access<2, double>(soln_cell_data.getArrayData());
      MDA_Access<double, 2, MDA_OrderColMajor<2> > es =
        pdat::ArrayDataAccess::access<2, double>(estimate_data.getArrayData());
      int i, j;
      double estimate, est0, est1, est2, est3, est4, est5;
      for (j = lower[1]; j <= upper[1]; ++j)
        {
          for (i = lower[0]; i <= upper[0]; ++i)
            {
              est0=tbox::MathUtilities<double>::Abs(co(i+1,j) + co(i-1,j)
                                                    - 2*co(i,j));
                                                 
              est1=tbox::MathUtilities<double>::Abs(co(i,j+1) + co(i,j-1)
                                                    - 2*co(i,j));
              est2=0.5 * tbox::MathUtilities<double>::Abs(co(i+1,j+1)
                                                          + co(i-1,j-1)
                                                          - 2*co(i,j));
              est3=0.5 * tbox::MathUtilities<double>::Abs(co(i+1,j-1)
                                                          + co(i-1,j+1)
                                                          - 2*co(i,j));
              est4=tbox::MathUtilities<double>::Max(est0,est1);
              est5=tbox::MathUtilities<double>::Max(est2,est3);
              estimate=tbox::MathUtilities<double>::Max(est4,est5);
              es(i,j)=estimate;
            }
        }
    }
  else if (d_dim == tbox::Dimension(3))
    {
      MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > co =
        pdat::ArrayDataAccess::access<3, double>(soln_cell_data.getArrayData());
      MDA_Access<double, 3, MDA_OrderColMajor<3> > es =
        pdat::ArrayDataAccess::access<3, double>(estimate_data.getArrayData());
      int i, j, k;
      double estimate, est0, est1, est2, est3, est4, est5, est6, est7, est8,
        esta, estb, estc, estd, este, estf, estg;
      for (k = lower[2]; k <= upper[2]; ++k)
        for (j = lower[1]; j <= upper[1]; ++j)
          for (i = lower[0]; i <= upper[0]; ++i)
            {
              est0=tbox::MathUtilities<double>::Abs(co(i+1,j,k) + co(i-1,j,k)
                                                    - 2*co(i,j,k));
              est1=tbox::MathUtilities<double>::Abs(co(i,j+1,k) + co(i,j-1,k)
                                                    - 2*co(i,j,k));
              est2=tbox::MathUtilities<double>::Abs(co(i,j,k+1) + co(i,j,k-1)
                                                    - 2*co(i,j,k));
              est3=0.5 * tbox::MathUtilities<double>::Abs(co(i,j+1,k+1)
                                                          + co(i,j-1,k-1)
                                                          - 2*co(i,j,k));
              est4=0.5 * tbox::MathUtilities<double>::Abs(co(i,j+1,k-1)
                                                          + co(i,j-1,k+1)
                                                          - 2*co(i,j,k));
              est5=0.5 * tbox::MathUtilities<double>::Abs(co(i+1,j,k+1)
                                                          + co(i-1,j,k-1)
                                                          - 2*co(i,j,k));
              est6=0.5 * tbox::MathUtilities<double>::Abs(co(i+1,j,k-1)
                                                          + co(i-1,j,k+1)
                                                          - 2*co(i,j,k));
              est7=0.5 * tbox::MathUtilities<double>::Abs(co(i+1,j+1,k)
                                                          + co(i-1,j-1,k)
                                                          - 2*co(i,j,k));
              est8=0.5 * tbox::MathUtilities<double>::Abs(co(i+1,j-1,k)
                                                          + co(i-1,j+1,k)
                                                          - 2*co(i,j,k));
              esta=tbox::MathUtilities<double>::Max(est0,est1);
              estb=tbox::MathUtilities<double>::Max(est2,est3);
              estc=tbox::MathUtilities<double>::Max(est4,est5);
              estd=tbox::MathUtilities<double>::Max(est6,est7);
              este=tbox::MathUtilities<double>::Max(esta,estb);
              estf=tbox::MathUtilities<double>::Max(estc,estd);
              estg=tbox::MathUtilities<double>::Max(este,estf);
              estimate=tbox::MathUtilities<double>::Max(estg,est8);
              es(i,j,k)=estimate;
            }
    }
}
