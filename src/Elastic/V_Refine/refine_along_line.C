#include "Elastic/V_Refine.h"

/* This assumes that the levels are always properly nested, so that we
   always have an extra grid space for interpolation.  So we only have
   to have a special case for physical boundaries, where we do not
   have an extra grid space. */

/* Maybe this has to be fixed when dvx/dy != 0 on the outer boundary
   because the approximation to the derivative is not accurate
   enough? */

/* Maybe in 3D we should include cross derivatives? */

double SAMRAI::geom::V_Refine::refine_along_line
(pdat::SideData<double> &v,
 const int &axis,
 const int &dim,
 const hier::Index pp[],
 const pdat::SideIndex &fine,
 const pdat::SideIndex &coarse,
 const hier::Box &coarse_box,
 const CartesianPatchGeometry &coarse_geom,
 const FTensor::Tensor1<double,3> &xyz,
 const double *dx) const
{
  double result=v(coarse);

  for(int d=(axis+1)%dim;d!=axis;d=(d+1)%dim)
    {
      const int sgn(fine[d]%2==0 ? -1 : 1);

      double dvx_dy;
      if(coarse[d]==coarse_box.lower(d)
         && coarse_geom.getTouchesRegularBoundary(d,0))
        {
          dvx_dy=sgn*(v(coarse+pp[d])-v(coarse))/4;
        }
      else if(coarse[d]==coarse_box.upper(d)
              && coarse_geom.getTouchesRegularBoundary(d,1))
        {
          dvx_dy=sgn*(v(coarse)-v(coarse-pp[d]))/4;
        }
      else
        {
          dvx_dy=sgn*(v(coarse+pp[d])-v(coarse-pp[d]))/8;

          if(axis==0 && xyz(0)-dx[0]<0.5+1e-6 && xyz(0)+dx[0]>0.5+1e-6)
            {
              const double y_min(0.5-sgn*0.1), y_max(0.5+sgn*0.1);
              /* Top tip */
              /* Singularity in coarse+1 */
              if(sgn*(xyz(1)+sgn*dx[1]*0.5-y_max)<0 && sgn*(xyz(1)+sgn*dx[1]*2.5-y_max)>0)
                {
                  dvx_dy=(v(coarse)-v(coarse-pp[d]*sgn))/4;
                  // tbox::pout << "refine p cp1 "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + dvx_dy << " "
                  //            << "\n";

                }
              /* Singularity in coarse */
              else if(sgn*(xyz(1)-sgn*dx[1]*1.5-y_max)<0 && sgn*(xyz(1)+sgn*dx[1]*0.5-y_max)>0)
                {
                  dvx_dy=0;
                  // tbox::pout << "refine p c "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + dvx_dy << " "
                  //            << "\n";
                }
              /* Singularity in coarse-1 */
              else if(sgn*(xyz(1)-sgn*dx[1]*3.5-y_max)<0 && sgn*(xyz(1)-sgn*dx[1]*1.5-y_max)>0)
                {
                  dvx_dy=(v(coarse+pp[d]*sgn)-v(coarse))/4;
                  // tbox::pout << "refine p cm1 "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + dvx_dy << " "
                  //            << "\n";
                }

              /* Bottom tip */
              /* Singularity in coarse+1 */
              else if(sgn*(xyz(1)+sgn*dx[1]*0.5-y_min)<0 && sgn*(xyz(1)+sgn*dx[1]*2.5-y_min)>0)
                {
                  dvx_dy=(v(coarse)-v(coarse-pp[d]*sgn))/4;
                  // tbox::pout << "refine m cp1 "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + dvx_dy << " "
                  //            << "\n";
                }
              /* Singularity in coarse */
              else if(sgn*(xyz(1)-sgn*dx[1]*1.5-y_min)<0 && sgn*(xyz(1)+sgn*dx[1]*0.5-y_min)>0)
                {
                  dvx_dy=0;
                  // tbox::pout << "refine m c "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + dvx_dy << " "
                  //            << "\n";
                }
              /* Singularity in coarse-1 */
              else if(sgn*(xyz(1)-sgn*dx[1]*3.5-y_min)<0 && sgn*(xyz(1)-sgn*dx[1]*1.5-y_min)>0)
                {
                  dvx_dy=(v(coarse+pp[d]*sgn)-v(coarse))/4;
                  // tbox::pout << "refine m cm1 "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + dvx_dy << " "
                  //            << "\n";
                }

              // /* Bottom tip */
              // /* Singularity in coarse+1 */
              // else if(sgn*(xyz(1)-sgn*dx[1]*1.5-y_min)<0 && sgn*(xyz(1)-sgn*dx[1]*3.5-y_min)>0)
              //   {
              //     dvx_dy=(v(coarse)-v(coarse-pp[d]*sgn))/4;
              //     tbox::pout << "refine m cp1 "
              //                << 1/dx[0] << " "
              //                << xyz(0) << " "
              //                << xyz(1) << " "
              //                << fine << " "
              //                << coarse << " "
              //                << sgn << " "
              //                << v(coarse+pp[d]*sgn) << " "
              //                << v(coarse) << " "
              //                << v(coarse-pp[d]*sgn) << " "
              //                << result + dvx_dy << " "
              //                << "\n";
              //   }
              // /* Jump in coarse */
              // else if(sgn*(xyz(1)-sgn*dx[1]*1.5-y_min)>0 && sgn*(xyz(1)+sgn*dx[1]*0.5-y_min)<0)
              //   {
              //     dvx_dy=(v(coarse+pp[d]*sgn)-v(coarse))/4;
              //     tbox::pout << "refine m c "
              //                << 1/dx[0] << " "
              //                << xyz(0) << " "
              //                << xyz(1) << " "
              //                << fine << " "
              //                << coarse << " "
              //                << sgn << " "
              //                << v(coarse+pp[d]*sgn) << " "
              //                << v(coarse) << " "
              //                << v(coarse-pp[d]*sgn) << " "
              //                << result + dvx_dy << " "
              //                << "\n";
              //   }

            }

          // if(coarse[1]==6)
          //       {
          //         tbox::pout << "refine "
          //                    << 1/dx[0] << " "
          //                    << xyz(0) << " "
          //                    << xyz(1) << " "
          //                    << fine << " "
          //                    << coarse << " "
          //                    << sgn << " "
          //                    << v(coarse+pp[d]*sgn) << " "
          //                    << v(coarse) << " "
          //                    << v(coarse-pp[d]*sgn) << " "
          //                    << result + dvx_dy << " "
          //                    << "\n";

          //       }



          if(axis==1 && xyz(1)>=0.4 && xyz(1)<=0.6)
            {
              /* Interface between fine and coarse+1, do a one-sided
                 interpolation. */
              if(sgn*(xyz(0)-0.5)<0 && sgn*(xyz(0)+sgn*1.5*dx[0]-0.5)>0)
                {
                  dvx_dy=(v(coarse)-v(coarse-pp[d]*sgn))/4;

                  // tbox::pout << "refine fcp1 "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + sgn*dvx_dy << " "
                  //            << "\n";
                }
              /* Interface between fine and coarse, use derivative
                 from the other side. */
              else if(sgn*(xyz(0)-0.5)>0 && sgn*(xyz(0)-sgn*dx[0]/2-0.5)<0)
                {
                  dvx_dy=(v(coarse+pp[d]*sgn)-v(coarse)
                          - (v(coarse)-v(coarse-pp[d]*sgn))*0.75);

                  // tbox::pout << "refine fc "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + sgn*dvx_dy << " "
                  //            << "\n";
                }
              /* Interface between coarse and coarse-1, do a one-sided
                 interpolation. */
              else if(sgn*(xyz(0)-0.5)>0 && sgn*(xyz(0)-sgn*2.5*dx[0]-0.5)<0)
                {
                  dvx_dy=(v(coarse+pp[d]*sgn)-v(coarse))/4;

                  // tbox::pout << "refine ccm1 "
                  //            << 1/dx[0] << " "
                  //            << xyz(0) << " "
                  //            << xyz(1) << " "
                  //            << fine << " "
                  //            << coarse << " "
                  //            << sgn << " "
                  //            << v(coarse+pp[d]*sgn) << " "
                  //            << v(coarse) << " "
                  //            << v(coarse-pp[d]*sgn) << " "
                  //            << result + sgn*dvx_dy << " "
                  //            << "\n";
                }
            }
        }
      result+=dvx_dy;
    }
  return result;
}

