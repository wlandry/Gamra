#include <cmath>

double okada(const double xyz[], const int &ix)
{
  const double x(xyz[0]-0.000001), yp(xyz[1]+0.1), ym(xyz[1]-0.1), b(10),
    pi(4*std::atan(1.0)), lam(1), mu(1);

  if(ix==0)
    {
      const double uxp=
        -(b/(2*pi))*(-mu*std::log(x*x+yp*yp)/(2*(lam+2*mu))
                     + ((lam + mu)/(lam+2*mu))*x*x/(x*x + yp*yp));
      const double uxm=
        -(b/(2*pi))*(-mu*std::log(x*x+ym*ym)/(2*(lam+2*mu))
                     + ((lam + mu)/(lam+2*mu))*x*x/(x*x + ym*ym));

      return -uxp+uxm;
    }
  else
    {
      const double uyp=
        -(b/(2*pi))*(-std::atan(yp/x)
                     + ((lam + mu)/(lam+2*mu))*(x*yp)/(x*x + yp*yp));
      const double uym=
        -(b/(2*pi))*(-std::atan(ym/x)
                     + ((lam + mu)/(lam+2*mu))*(x*ym)/(x*x + ym*ym));
      return -uyp+uym;
    }
}
