#include "P_Boundary_Refine.h"
#include "quad_offset_interpolate.h"

/* Interpolate the pressure from the coarse (C) to the fine (f+/-)
   coordinates using the intermediate fine points (F+/-, F_).

   i-1      i       i+1

   ------- -------
   |       |       |
   j-1 |  C    |   C   |   C
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

   F+/- = interpolate(C(i,j-1),C(i,j),C(i,j+1))
   F_   = interpolate(C(i-1,j),C(i,j),C(i+1,j))

   then

   f+/- = F+/- + F_ - C(i,j) 
   + (C(i-1,j-1) + C(i,j) + C(i-1,j) + C(i,j-1))/16

   This example show a boundary in the positive x direction.  To reverse
   the direction, pass in ip -> -ip.  To do the y direction, switch ip
   and jp, and replace j with i.

*/
     

void SAMRAI::geom::P_Boundary_Refine::Update_P_2D
(const pdat::CellIndex &fine,
 const hier::Index &ip, const hier::Index &jp,
 int &j,
 tbox::Pointer<SAMRAI::pdat::CellData<double> > &p,
 tbox::Pointer<SAMRAI::pdat::CellData<double> > &p_fine) const
{
  double p_plus, p_minus, p_offset;
  pdat::CellIndex center(fine);
  center.coarsen(hier::Index(2,2));

  quad_offset_interpolate((*p)(center+jp),(*p)(center),(*p)(center-jp),
                          p_plus,p_minus);
  p_offset=
    quad_offset_interpolate((*p)(center-ip),(*p)(center),(*p)(center+ip));

  const double p_low=p_minus + p_offset - (*p)(center)
    + ((*p)(center-ip-jp) + (*p)(center)
       - (*p)(center-ip) - (*p)(center-jp))/16;

  const double p_high=p_plus + p_offset - (*p)(center)
    + ((*p)(center-ip+jp) + (*p)(center)
       - (*p)(center-ip) - (*p)(center+jp))/16;


  /* If we are at an even index, update both of the elements in the cell */
  if(j%2==0)
    {
      (*p_fine)(fine)=p_low;

      (*p_fine)(fine+jp)=p_high;

      /* Since we update two points on j at once, we increment j again.
         This is ok, since the box in the 'i' direction is defined to be
         only one cell wide */
      ++j;
    }
  else
    {
      (*p_fine)(fine)=p_high;
    }

  // tbox::plog << "p bc "
  //            << fine[0] << " "
  //            << fine[1] << " "
  //            << center[0] << " "
  //            << center[1] << " "
  //            << jp[0] << " "
  //            << jp[1] << " "
  //            << (*p_fine)(fine) << " "
  //            << (*p_fine)(fine+jp) << " "
  //            << p_minus << " "
  //            << p_plus << " "
  //            << p_offset << " "
  //            << (*p)(center+jp) << " "
  //            << (*p)(center) << " "
  //            << (*p)(center-jp) << " "
  //            << (*p)(center+ip) << " "
  //            << (*p)(center-ip) << " "
  //            << (*p)(center-ip+jp) << " "
  //            << (*p)(center-ip-jp) << " "
  //            << "\n";
}
