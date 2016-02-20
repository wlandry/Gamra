#include "Stokes/P_Boundary_Refine.hxx"

/* Interpolate the pressure from the coarse (C) to the fine (f+/-)
   coordinates using the intermediate fine points (F+/-, F_).

         i-1      i       i+1

        ------- -------
       |       |       |
   j-1 |   C   |   C   |   C
       |       |       |
        ------- -------
       |       |f- F-  |
   j   |   C   |F_ C   |   C
       |       |f+ F+  |
        ------- -------
       |       |       |
   j+1 |   C   |   C   |   C
       |       |       |
        ------- -------

   C= a + b*x + c*x^2 + d*y + e*y^2 + f*x*y

   C(0,0)=a
   C(+,0)=a+b+c
   C(-,0)=a-b+c
   C(0,+)=a+d+e
   C(0,-)=a-d+e
   C(-,-)=a-b+c-d+e+f

   f(-,-) = a - b/4 + c/16 - d/4 + e/16 + f/16
          = C(-,-)/16 + (15/16)*C(0,0)
            + (3/32)*(-C(+,0) - C(0,+) + C(-,0) + C(0,-))
   
   This example show a boundary in the positive x direction.  To reverse
   the direction, pass in ip -> -ip.  To do the y direction, switch ip
   and jp, and replace j with i.

*/
     

void Stokes::P_Boundary_Refine::Update_P_2D
(const SAMRAI::pdat::CellIndex &fine,
 const SAMRAI::hier::IntVector &ip, const SAMRAI::hier::IntVector &jp,
 const int &j, const int &j_max,
 SAMRAI::pdat::CellData<double> &p,
 SAMRAI::pdat::CellData<double> &p_fine) const
{
  SAMRAI::pdat::CellIndex center(fine);
  center.coarsen(SAMRAI::hier::Index(2,2));

    
  /* If we are at an even index, update both of the elements in the cell */
  if(j%2==0)
    {
      const double p_low=p(center-ip-jp)/16 + (15.0/16)*p(center)
        + (3.0/32)*(-p(center+ip) - p(center+jp) + p(center-ip) + p(center-jp));
      p_fine(fine)=p_low;
      if(j<j_max)
        {
          const double p_high=p(center-ip+jp)/16 + (15.0/16)*p(center)
            + (3.0/32)*(-p(center+ip) - p(center-jp)
                        + p(center-ip) + p(center+jp));
          p_fine(fine+jp)=p_high;
        }
    }
  else
    {
      const double p_high=p(center-ip+jp)/16 + (15.0/16)*p(center)
        + (3.0/32)*(-p(center+ip) - p(center-jp)
                    + p(center-ip) + p(center+jp));
      p_fine(fine)=p_high;
    }
}
