#include "Stokes/P_Boundary_Refine.hxx"

/* Interpolate the pressure from the coarse (C) to the fine (f+/-)
   coordinates using the intermediate fine points (F+/-, F_).

          i-1      i       i+1

        ------- -------
       |       |       |
   j-1 |  C    |   C   |   C
       |       |       |
       ------- -------
       |       |f- Y-  |
   j   |   C   |X_ C   |   C
       |       |f+ Y+  |
       ------- -------
       |       |       |
   j+1 |   C   |   C   |   C
       |       |       |
       ------- -------


   C= a + b*x + c*x^2 + d*y + e*y^2 + f*z + g*z^2
        + h*x*y + i*x*z + j*y*z + k*x*y*z

   C(0,0,0)=a
   C(+,0,0)=a+b+c
   C(-,0,0)=a-b+c
   C(0,+,0)=a+d+e
   C(0,-,0)=a-d+e
   C(-,-,0)=a-b+c-d+e+h
   C(-,-,-)=a-b+c-d+e-f+g+h+i+j-k

   f(-,-) = a - b/4 + c/16 - d/4 + e/16 - f/4 + g/16 + h/16 + i/16 + j/16 - k/64
          = C(-,-)/64 + (63/64)*C(0,0)
            - (6/64)*(C(+,0,0) + C(0,+,0) + C(0,0,+))
            + (3/64)*(C(-,0,0) + C(0,-,0) + C(0,0,-)
                      + C(-,-,0) + C(-,0,-) + C(0,-,-))

   Note that if, for example, C is constant in the Z direction, then
   this formula reduces to the one used in Update_P_2D.

   This example show a boundary in the positive x direction.  To
   reverse the direction, pass in ip -> -ip.  To do the y direction,
   rotate [ip,jp,kp], and switch j with k.  To do z, rotate again and
   switch k with i.
*/
     

void Stokes::P_Boundary_Refine::Update_P_3D
(const SAMRAI::pdat::CellIndex &fine,
 const SAMRAI::hier::IntVector &ip, const SAMRAI::hier::IntVector &jp,
 const SAMRAI::hier::IntVector &kp,
 const int &j, const int &k, const int &j_max, const int &k_max,
 SAMRAI::pdat::CellData<double> &p,
 SAMRAI::pdat::CellData<double> &p_fine) const
{
  SAMRAI::pdat::CellIndex center(fine);
  center.coarsen(SAMRAI::hier::Index(2,2,2));

  const double p_mmm=p(center-ip-jp-kp)/64 + (63.0/64)*p(center)
    - (3.0/32) * (p(center+ip) + p(center+jp) + p(center+kp))
    + (3.0/64) * (p(center-ip) + p(center-jp) + p(center-kp)
                  + p(center-ip-jp) + p(center-jp-kp) + p(center-ip-kp));

  const double p_mmp=p(center-ip-jp+kp)/64 + (63.0/64)*p(center)
    - (3.0/32) * (p(center+ip) + p(center+jp) + p(center-kp))
    + (3.0/64) * (p(center-ip) + p(center-jp) + p(center+kp)
                  + p(center-ip-jp) + p(center-jp+kp) + p(center-ip+kp));

  const double p_mpm=p(center-ip+jp-kp)/64 + (63.0/64)*p(center)
    - (3.0/32) * (p(center+ip) + p(center-jp) + p(center+kp))
    + (3.0/64) * (p(center-ip) + p(center+jp) + p(center-kp)
                  + p(center-ip+jp) + p(center+jp-kp) + p(center-ip-kp));

  const double p_mpp=p(center-ip+jp+kp)/64 + (63.0/64)*p(center)
    - (3.0/32) * (p(center+ip) + p(center-jp) + p(center-kp))
    + (3.0/64) * (p(center-ip) + p(center+jp) + p(center+kp)
                  + p(center-ip+jp) + p(center+jp+kp) + p(center-ip+kp));

  /* If we are at an even index, update both of the elements in the cell */
  if(j%2==0)
    {
      if(k%2==0)
        {
          p_fine(fine)=p_mmm;
          if(j<j_max)
            p_fine(fine+jp)=p_mpm;
          if(k<k_max)
            p_fine(fine+kp)=p_mmp;
          if(j<j_max && k<k_max)
            p_fine(fine+jp+kp)=p_mpp;
        }
      else
        {
          p_fine(fine)=p_mmp;
          if(j<j_max)
            p_fine(fine+jp)=p_mpp;
        }
    }
  else
    {
      if(k%2==0)
        {
          p_fine(fine)=p_mpm;
          if(k<k_max)
            p_fine(fine+kp)=p_mpp;
        }
      else
        {
          p_fine(fine)=p_mpp;
        }
    }
}
