#include <cmath>

double okada_tangent_deriv(const double xyz[],
                           const int &ix, const int &)
{
  const double x(xyz[0]-0.000001), yp(xyz[1]+0.1), ym(xyz[1]-0.1), b(10),
    pi(4*std::atan(1.0)), lam(1), mu(1);

  if(ix==0)
    {
      const double duxp_y=(b/(2*pi))*yp*(mu*yp*yp + (2*lam + 3*mu)*x*x)
        /((lam + 2*mu)*(x*x + yp*yp)*(x*x + yp*yp));
      const double duxm_y=(b/(2*pi))*ym*(mu*ym*ym + (2*lam + 3*mu)*x*x)
        /((lam + 2*mu)*(x*x + ym*ym)*(x*x + ym*ym));
      return -duxp_y+duxm_y;
    }
  else
    {
      const double duyp_x=(-b/(2*pi))*yp*(mu*x*x + (2*lam + 3*mu)*yp*yp)
        /((lam + 2*mu)*(x*x + yp*yp)*(x*x + yp*yp));
      const double duym_x=(-b/(2*pi))*ym*(mu*x*x + (2*lam + 3*mu)*ym*ym)
        /((lam + 2*mu)*(x*x + ym*ym)*(x*x + ym*ym));
      return -duyp_x+duym_x;
    }
}
