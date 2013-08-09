#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

namespace {
  double complete_correction
  (const bool &boundary_positive,
   const int &axis,
   const SAMRAI::pdat::SideIndex &fine,
   const SAMRAI::hier::IntVector &ip_s,
   const double &coarse_correction,
   const SAMRAI::pdat::CellData<double> &dv_diagonal_fine)
  {
    double result;
    const SAMRAI::pdat::CellIndex cell(fine);
    if(boundary_positive)
      {
        result=dv_diagonal_fine(cell-ip_s,axis)
          + (coarse_correction - dv_diagonal_fine(cell-ip_s-ip_s,axis))/3;
      }
    else
      {
        result=-dv_diagonal_fine(cell,axis)
          + (coarse_correction + dv_diagonal_fine(cell-ip_s,axis))/3;
      }
    return result;
  }

  double linear_coarse_correction
  (const bool &boundary_positive,
   const SAMRAI::pdat::SideIndex &fine,
   const SAMRAI::pdat::SideIndex &coarse,
   const SAMRAI::pdat::SideData<double> &dv_mixed,
   const SAMRAI::pdat::SideData<double> &dv_mixed_fine)
  {
    return boundary_positive ? 
      (8*(dv_mixed(coarse,1) - dv_mixed_fine(fine,1))/15) :
      (8*(dv_mixed(coarse,0) - dv_mixed_fine(fine,0))/15) ;
  }

  double in_between_coarse_correction
  (const bool &boundary_positive,
   const int &axis,
   const SAMRAI::pdat::SideIndex &fine,
   const SAMRAI::pdat::SideIndex &coarse,
   const SAMRAI::hier::IntVector &ip,
   const int &i_min,
   const int &i_max,
   const SAMRAI::pdat::SideData<double> &dv_mixed,
   const SAMRAI::pdat::SideData<double> &dv_mixed_fine,
   const SAMRAI::pdat::CellData<double> &dv_diagonal,
   const SAMRAI::pdat::CellData<double> &dv_diagonal_fine)
  {
    const SAMRAI::pdat::CellIndex cell(fine), cell_coarse(coarse);
    
    if(fine[axis]==i_max)
      {
        return linear_coarse_correction(boundary_positive,fine-ip,
                                        coarse,dv_mixed,dv_mixed_fine)
          + 4*(-dv_diagonal(cell_coarse,axis)
               + 2*dv_diagonal_fine(cell-ip,axis))/15;
      }
    else if(fine[axis]==i_min)
      {
        return linear_coarse_correction(boundary_positive,fine+ip,
                                        coarse+ip,dv_mixed,dv_mixed_fine)
          + 4*(dv_diagonal(cell_coarse,axis)
               - 2*dv_diagonal_fine(cell,axis))/15;
      }

    return (linear_coarse_correction(boundary_positive,fine-ip,coarse,
                                     dv_mixed,dv_mixed_fine)
            + linear_coarse_correction(boundary_positive,fine+ip,coarse+ip,
                                       dv_mixed,dv_mixed_fine))/2
      + 4*(dv_diagonal_fine(cell-ip,axis) - dv_diagonal_fine(cell,axis))/15;
  }

  double linear_fine_correction
  (const bool &boundary_positive,
   const SAMRAI::pdat::SideIndex &fine,
   const SAMRAI::hier::IntVector &jp_s,
   const SAMRAI::pdat::SideData<double> &dv_mixed_fine)
  {
    double result;
    if(boundary_positive)
      {
        result=(7*(dv_mixed_fine(fine-jp_s,0) - dv_mixed_fine(fine,1))
                -3*(dv_mixed_fine(fine-jp_s-jp_s,0)
                    - dv_mixed_fine(fine-jp_s,1)))/15;
          
      }
    else
      {
        result=(7*(dv_mixed_fine(fine-jp_s,1) - dv_mixed_fine(fine,0))
                -3*(dv_mixed_fine(fine-jp_s-jp_s,1)
                    - dv_mixed_fine(fine-jp_s,0)))/15;
      }
    return result;
  }
}

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void Elastic::V_Boundary_Refine::Correction_2D
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::hier::IntVector &ip,
 const SAMRAI::hier::IntVector &jp,
 const int &i,
 const int &j,
 const int &i_min,
 const int &i_max, 
 const SAMRAI::pdat::CellData<double> &dv_diagonal,
 const SAMRAI::pdat::CellData<double> &dv_diagonal_fine,
 const SAMRAI::pdat::SideData<double> &dv_mixed,
 const SAMRAI::pdat::SideData<double> &dv_mixed_fine,
 SAMRAI::pdat::SideData<double> &v_fine) const
{
  /* Quadratic interpolation involving both coarse and fine grids for
     the normal direction

      i-1      i      i+1

        ------- -------
       |   f   f   F   |
   j-1 C       C       C
       |   f   f   F   |
        ------- -------
       |   f   f   F   |
   j   C       C       C
       |   f   f   F   |
        ------- -------
       |   f   f   F   |
   j+1 C       C       C
       |   f   f   F   |
        ------- -------
               |
               |
               |
        Coarse-Fine Boundary

      Interpolate to F.
  */

  if(boundary_direction==axis)
    {
      SAMRAI::hier::IntVector ip_s(boundary_positive ? ip : -ip);
      SAMRAI::pdat::SideIndex coarse(fine-ip_s);
      coarse.coarsen(SAMRAI::hier::Index(2,2));

      double correction[2],correction_diagonal;
      quad_offset_correction(dv_mixed(coarse+ip_s+jp,1),
                             dv_mixed(coarse+ip_s,0),
                             dv_mixed(coarse+ip_s,1),
                             dv_mixed(coarse+ip_s-jp,0),
                             correction[1], correction[0]);

      SAMRAI::pdat::CellIndex cell(coarse);
      correction_diagonal=(boundary_positive ? -dv_diagonal(cell,axis) :
                           dv_diagonal(cell+ip_s,axis));

      double coarse_correction;
      coarse_correction=correction_diagonal + correction[j%2]
        - dv_mixed_fine(fine-ip_s,j%2);

      v_fine(fine)+=complete_correction(boundary_positive,axis,fine,ip_s,
                                        coarse_correction,dv_diagonal_fine);
    }
  /* Quadratic interpolation involving both coarse and fine grids for
     the tangential direction

      i-1      i      i+1

        ------- -------
       |       |       |
   j-1 C       C       C
       F   F   F   F   F
        ------- -------    --- Coarse-Fine Boundary
       f   f   f   f   f
   j   C       C       C
       f   f   f   f   f
        ------- -------
       f   f   f   f   f
   j+1 C       C       C
       f   f   f   f   f
        ------- -------

    C are the coarse velocities, f are the interior fine velocities,
    and F are the boundary fine velocities that we need to set.
 */


  else
    {
      SAMRAI::pdat::SideIndex coarse(fine);
      coarse.coarsen(SAMRAI::hier::Index(2,2));

      SAMRAI::hier::Index jp_s(boundary_positive ? jp : -jp);

      double fine_correction(linear_fine_correction(boundary_positive,fine,jp_s,
                                                    dv_mixed_fine));
      if(i%2==0)
        {
          v_fine(fine)+=(fine_correction
            + linear_coarse_correction(boundary_positive,fine,coarse,
                                       dv_mixed,dv_mixed_fine));
        }
      else
        {
          v_fine(fine)+=(fine_correction
            + in_between_coarse_correction(boundary_positive,axis,fine,coarse,ip,
                                           i_min,i_max,dv_mixed,dv_mixed_fine,
                                           dv_diagonal,dv_diagonal_fine));
        }
    }
}
